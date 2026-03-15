"""
Blender render endpoints — async job queue for high-quality geometry renders.

POST /render                      submit a render job → {job_id, status}
GET  /render/{job_id}/status      poll job status + progress
GET  /render/{job_id}/result      download PNG when done
"""
from __future__ import annotations

import os
import queue
import shutil
import struct
import subprocess
import tempfile
import textwrap
import threading
import time
import uuid
from pathlib import Path
from typing import Any, Dict, List, Optional

from fastapi import APIRouter, HTTPException
from fastapi.responses import FileResponse, JSONResponse
from pydantic import BaseModel

router = APIRouter(prefix="/render", tags=["render"])

# ── Blender binary discovery ──────────────────────────────────────────────────

def _find_blender() -> Optional[str]:
    candidates = [
        shutil.which("blender"),
        "/Applications/Blender.app/Contents/MacOS/Blender",
        "/opt/homebrew/bin/blender",
    ]
    for c in candidates:
        if c and Path(c).exists():
            return c
    return None

BLENDER_BIN = _find_blender()

# ── In-memory job store ───────────────────────────────────────────────────────

_jobs: Dict[str, Dict[str, Any]] = {}   # job_id → job dict
_job_queue: "queue.Queue[str]" = queue.Queue()


def _new_job(config: Dict[str, Any]) -> str:
    job_id = str(uuid.uuid4())
    _jobs[job_id] = {
        "job_id": job_id,
        "status": "pending",
        "progress": 0,
        "error": None,
        "result_path": None,
        "config": config,
        "created_at": time.time(),
    }
    _job_queue.put(job_id)
    return job_id


# ── STL writer (from geometry JSON triangles) ─────────────────────────────────

def _write_stl(geometry: Dict[str, Any], path: str) -> None:
    """Write a binary STL from the geometry JSON returned by /dsl/eval."""
    entities = geometry.get("entities", [])
    all_triangles: List[List] = []
    for ent in entities:
        tris = ent.get("triangles", [])
        verts = ent.get("vertices", [])
        if tris and verts:
            for tri in tris:
                pts = [verts[i] for i in tri]
                all_triangles.append(pts)

    with open(path, "wb") as f:
        f.write(b"\x00" * 80)   # header
        f.write(struct.pack("<I", len(all_triangles)))
        for tri in all_triangles:
            # Normal (0,0,0 → Blender recalculates)
            f.write(struct.pack("<fff", 0.0, 0.0, 0.0))
            for v in tri:
                f.write(struct.pack("<fff", float(v[0]), float(v[1]), float(v[2])))
            f.write(struct.pack("<H", 0))   # attribute byte count


# ── Blender render script template ───────────────────────────────────────────

_BLENDER_SCRIPT = textwrap.dedent("""
import bpy, sys, math

stl_path   = {stl_path!r}
output_png = {output_png!r}
camera_preset  = {camera!r}
lighting_preset = {lighting!r}
quality    = {quality!r}
background = {background!r}

# Clear scene
bpy.ops.wm.read_factory_settings(use_empty=True)
bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete()

# Import STL
bpy.ops.import_mesh.stl(filepath=stl_path)
obj = bpy.context.selected_objects[0]
bpy.context.view_layer.objects.active = obj

# Centre geometry
bpy.ops.object.origin_set(type='GEOMETRY')
obj.location = (0, 0, 0)

# Calculate bounding box for camera placement
bbox = [obj.matrix_world @ v.co for v in obj.data.vertices]
if not bbox:
    bbox_min = (-1,-1,-1)
    bbox_max = (1,1,1)
else:
    bbox_min = [min(v[i] for v in bbox) for i in range(3)]
    bbox_max = [max(v[i] for v in bbox) for i in range(3)]

cx = (bbox_min[0] + bbox_max[0]) / 2
cy = (bbox_min[1] + bbox_max[1]) / 2
cz = (bbox_min[2] + bbox_max[2]) / 2
size = max(bbox_max[i] - bbox_min[i] for i in range(3))
dist = size * 2.0

# Material
mat = bpy.data.materials.new("DML_Material")
mat.use_nodes = True
nodes = mat.node_tree.nodes
bsdf = nodes.get("Principled BSDF")
if bsdf:
    if background == "dark":
        bsdf.inputs["Base Color"].default_value = (0.4, 0.45, 0.55, 1.0)
    else:
        bsdf.inputs["Base Color"].default_value = (0.7, 0.72, 0.8, 1.0)
    bsdf.inputs["Metallic"].default_value = 0.3
    bsdf.inputs["Roughness"].default_value = 0.4
obj.data.materials.append(mat)
obj.active_material = mat

# Camera
bpy.ops.object.camera_add()
cam = bpy.context.object
bpy.context.scene.camera = cam

cam_presets = {{
    "auto":  (cx + dist*0.7, cy - dist*0.7, cz + dist*0.5),
    "front": (cx, cy - dist, cz),
    "side":  (cx + dist, cy, cz),
    "top":   (cx, cy, cz + dist),
    "iso":   (cx + dist*0.6, cy - dist*0.6, cz + dist*0.6),
}}
cam.location = cam_presets.get(camera_preset, cam_presets["auto"])
# Point at centre
dx = cx - cam.location.x
dy = cy - cam.location.y
dz = cz - cam.location.z
import mathutils
direction = mathutils.Vector((dx, dy, dz))
rot_quat = direction.to_track_quat('-Z', 'Y')
cam.rotation_euler = rot_quat.to_euler()

# Lighting
world = bpy.data.worlds.new("DML_World")
bpy.context.scene.world = world
world.use_nodes = True
bg_node = world.node_tree.nodes.get("Background")
if lighting_preset == "studio":
    bg_node.inputs["Color"].default_value = (0.05, 0.05, 0.06, 1.0)
    bg_node.inputs["Strength"].default_value = 1.5
    # Add soft area lights
    for loc, energy in [((dist,0,dist*1.5), 200), ((-dist*0.5,dist,dist), 80)]:
        bpy.ops.object.light_add(type='AREA', location=(loc[0]+cx, loc[1]+cy, loc[2]+cz))
        light = bpy.context.object
        light.data.energy = energy * (size / 100) ** 2
        light.data.size = size * 0.8
elif lighting_preset == "outdoor":
    bg_node.inputs["Color"].default_value = (0.4, 0.6, 0.9, 1.0)
    bg_node.inputs["Strength"].default_value = 0.8
    bpy.ops.object.light_add(type='SUN', location=(cx + dist, cy, cz + dist))
    sun = bpy.context.object
    sun.data.energy = 4.0
    sun.rotation_euler = (math.radians(-45), 0, math.radians(45))
elif lighting_preset == "dramatic":
    bg_node.inputs["Color"].default_value = (0.01, 0.01, 0.01, 1.0)
    bg_node.inputs["Strength"].default_value = 0.0
    bpy.ops.object.light_add(type='SPOT', location=(cx + dist*0.3, cy - dist*0.3, cz + dist))
    spot = bpy.context.object
    spot.data.energy = 1000 * (size / 100) ** 2
    spot.data.spot_size = math.radians(40)

# Render settings
scene = bpy.context.scene
scene.render.engine = 'CYCLES'
scene.render.filepath = output_png
scene.render.image_settings.file_format = 'PNG'

quality_settings = {{
    "preview":     {{"samples": 64,  "x": 800,  "y": 600,  "denoise": False}},
    "review":      {{"samples": 256, "x": 1920, "y": 1080, "denoise": True}},
    "publication": {{"samples": 512, "x": 2560, "y": 1440, "denoise": True}},
}}
qs = quality_settings.get(quality, quality_settings["preview"])
scene.cycles.samples = qs["samples"]
scene.render.resolution_x = qs["x"]
scene.render.resolution_y = qs["y"]
if qs["denoise"]:
    scene.cycles.use_denoising = True

# Try Metal GPU on Apple Silicon
prefs = bpy.context.preferences.addons.get("cycles")
if prefs:
    try:
        cprefs = prefs.preferences
        cprefs.compute_device_type = 'METAL'
        for dev in cprefs.get_devices_for_type('METAL'):
            dev.use = True
        scene.cycles.device = 'GPU'
    except Exception:
        scene.cycles.device = 'CPU'

# Background colour
if background == "white":
    scene.render.film_transparent = False
    if bg_node:
        bg_node.inputs["Color"].default_value = (1.0, 1.0, 1.0, 1.0)
        bg_node.inputs["Strength"].default_value = 1.0
elif background == "transparent":
    scene.render.film_transparent = True

# Render
bpy.ops.render.render(write_still=True)
print("RENDER_COMPLETE:" + output_png)
""")


# ── Worker thread ─────────────────────────────────────────────────────────────

def _worker():
    while True:
        job_id = _job_queue.get()
        job = _jobs.get(job_id)
        if not job:
            continue

        job["status"] = "running"
        job["progress"] = 5

        tmpdir = tempfile.mkdtemp(prefix="yapcad_render_")
        try:
            if not BLENDER_BIN:
                raise RuntimeError("Blender not found on this machine")

            stl_path = os.path.join(tmpdir, "geometry.stl")
            script_path = os.path.join(tmpdir, "render.py")
            output_png = os.path.join(tmpdir, "render.png")

            cfg = job["config"]
            _write_stl(cfg["geometry"], stl_path)
            job["progress"] = 20

            script = _BLENDER_SCRIPT.format(
                stl_path=stl_path,
                output_png=output_png,
                camera=cfg.get("camera", "auto"),
                lighting=cfg.get("lighting", "studio"),
                quality=cfg.get("quality", "preview"),
                background=cfg.get("background", "dark"),
            )
            with open(script_path, "w") as f:
                f.write(script)

            job["progress"] = 30

            proc = subprocess.run(
                [BLENDER_BIN, "--background", "--python", script_path],
                capture_output=True, text=True, timeout=600,
            )

            job["progress"] = 90

            if proc.returncode != 0 or not Path(output_png).exists():
                raise RuntimeError(
                    f"Blender exited {proc.returncode}:\n{proc.stderr[-1000:]}"
                )

            # Move result to a persistent temp location (cleaned up after download)
            result_dir = tempfile.mkdtemp(prefix="yapcad_result_")
            result_png = os.path.join(result_dir, "render.png")
            shutil.move(output_png, result_png)

            job["result_path"] = result_png
            job["status"] = "done"
            job["progress"] = 100

        except Exception as e:
            job["status"] = "error"
            job["error"] = str(e)
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)


_worker_thread = threading.Thread(target=_worker, daemon=True)
_worker_thread.start()


# ── Models ────────────────────────────────────────────────────────────────────

class RenderRequest(BaseModel):
    geometry: Dict[str, Any]
    camera: str = "auto"       # auto | front | side | top | iso
    lighting: str = "studio"   # studio | outdoor | dramatic
    quality: str = "preview"   # preview | review | publication
    background: str = "dark"   # dark | white | transparent


class RenderStatus(BaseModel):
    job_id: str
    status: str
    progress: int
    error: Optional[str]


# ── Routes ────────────────────────────────────────────────────────────────────

@router.post("", status_code=202)
async def submit_render(body: RenderRequest):
    if not BLENDER_BIN:
        raise HTTPException(503, "Blender is not installed on this server")
    job_id = _new_job(body.dict())
    return {"job_id": job_id, "status": "pending"}


@router.get("/{job_id}/status", response_model=RenderStatus)
async def get_render_status(job_id: str):
    job = _jobs.get(job_id)
    if not job:
        raise HTTPException(404, "Job not found")
    return RenderStatus(
        job_id=job_id,
        status=job["status"],
        progress=job["progress"],
        error=job.get("error"),
    )


@router.get("/{job_id}/result")
async def get_render_result(job_id: str):
    job = _jobs.get(job_id)
    if not job:
        raise HTTPException(404, "Job not found")
    if job["status"] != "done":
        return JSONResponse(
            {"job_id": job_id, "status": job["status"], "progress": job["progress"]},
            status_code=202,
        )
    result_path = job.get("result_path")
    if not result_path or not Path(result_path).exists():
        raise HTTPException(500, "Result file missing")
    return FileResponse(result_path, media_type="image/png", filename="render.png")
