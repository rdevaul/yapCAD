"""Interactive viewer for yapCAD `.ycpkg` packages."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import pyglet
from pyglet import graphics
from pyglet.window import key
from pyglet.gl import (
    GL_COLOR_BUFFER_BIT,
    GL_DEPTH_BUFFER_BIT,
    GL_DEPTH_TEST,
    GL_COLOR_MATERIAL,
    GL_LIGHT0,
    GL_LIGHTING,
    GL_MODELVIEW,
    GL_ONE_MINUS_SRC_ALPHA,
    GL_PROJECTION,
    GL_SCISSOR_TEST,
    GL_LINES,
    GL_LINE_LOOP,
    GL_QUADS,
    GL_ENABLE_BIT,
    GL_DEPTH_BUFFER_BIT,
    GL_SRC_ALPHA,
    GL_TRIANGLES,
    GL_LINE,
    GL_FILL,
    GL_FRONT_AND_BACK,
    GL_POLYGON_OFFSET_LINE,
    glPolygonOffset,
    glBegin,
    glBlendFunc,
    glClear,
    glClearColor,
    glColorMaterial,
    glDisable,
    glEnable,
    glEnd,
    glLightfv,
    glLoadIdentity,
    glMatrixMode,
    glNormal3f,
    glOrtho,
    glPopMatrix,
    glPopAttrib,
    glPushMatrix,
    glPushAttrib,
    glScalef,
    glScissor,
    glTranslatef,
    glVertex2f,
    glVertex3f,
    glViewport,
    gluLookAt,
    gluPerspective,
    glColor4f,
    glPolygonMode,
    glLineWidth,
    GLfloat,
)

from yapcad.io.geometry_json import geometry_from_json
from .core import PackageManifest
from .validator import validate_package


def _tessellate_brep(brep_base64: str) -> List[Tuple[Tuple[float, float, float], Tuple[float, float, float]]]:
    """Tessellate BREP data to triangles for display.

    Args:
        brep_base64: Base64-encoded BREP ASCII string

    Returns:
        List of (vertex, normal) tuples for OpenGL rendering
    """
    import base64
    try:
        from OCC.Core.BRepMesh import BRepMesh_IncrementalMesh
        from OCC.Core.TopExp import TopExp_Explorer
        from OCC.Core.TopAbs import TopAbs_FACE
        from OCC.Core.TopLoc import TopLoc_Location
        from OCC.Core.BRep import BRep_Tool
        from OCC.Core.TopoDS import topods
        from OCC.Core.BRepTools import breptools
    except ImportError:
        print("Warning: pythonocc-core not available for BREP tessellation")
        return []

    # Decode BREP
    brep_ascii = base64.b64decode(brep_base64).decode('utf-8')

    # Restore shape from BREP string
    shape = breptools.ReadFromString(brep_ascii)
    if shape is None or shape.IsNull():
        return []

    # Mesh the shape
    mesh = BRepMesh_IncrementalMesh(shape, 0.5)  # 0.5mm linear deflection
    mesh.Perform()

    # Extract triangles from all faces
    triangles = []
    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    while explorer.More():
        face = topods.Face(explorer.Current())
        location = TopLoc_Location()
        triangulation = BRep_Tool.Triangulation(face, location)

        if triangulation is not None:
            transform = location.Transformation()

            for i in range(1, triangulation.NbTriangles() + 1):
                tri = triangulation.Triangle(i)
                n1, n2, n3 = tri.Get()

                # Get vertices and transform
                pts = []
                for n in [n1, n2, n3]:
                    p = triangulation.Node(n)
                    p.Transform(transform)
                    pts.append((p.X(), p.Y(), p.Z()))

                # Compute face normal from triangle
                v1 = (pts[1][0] - pts[0][0], pts[1][1] - pts[0][1], pts[1][2] - pts[0][2])
                v2 = (pts[2][0] - pts[0][0], pts[2][1] - pts[0][1], pts[2][2] - pts[0][2])
                nx = v1[1]*v2[2] - v1[2]*v2[1]
                ny = v1[2]*v2[0] - v1[0]*v2[2]
                nz = v1[0]*v2[1] - v1[1]*v2[0]
                mag = (nx*nx + ny*ny + nz*nz) ** 0.5
                if mag > 1e-10:
                    nx, ny, nz = nx/mag, ny/mag, nz/mag
                else:
                    nx, ny, nz = 0, 0, 1

                normal = (nx, ny, nz)
                for pt in pts:
                    triangles.append((pt, normal))

        explorer.Next()

    return triangles


def _load_geometry_doc(package_root: Path, manifest_data: Dict[str, any]) -> Dict[str, any]:
    primary_path = package_root / manifest_data["geometry"]["primary"]["path"]
    with primary_path.open("r", encoding="utf-8") as fp:
        return json.load(fp)


def _collect_triangles(doc: Dict[str, any], materials: Dict[str, any] = None) -> Tuple[Dict[str, List[Tuple[Tuple[float, float, float], Tuple[float, float, float]]]], Dict[str, Tuple[float, float, float]], Tuple[float, float, float, float, float, float]]:
    """Collect triangles organized by material (or layer if no material).

    Returns:
        - Dictionary mapping material/layer ID to list of (vertex, normal) tuples
        - Dictionary mapping material/layer ID to RGB color tuple
        - Bounding box tuple (xmin, ymin, zmin, xmax, ymax, zmax)
    """
    materials = materials or {}
    # Use material ID as bucket key when available, otherwise layer
    bucket_tris: Dict[str, List[Tuple[Tuple[float, float, float], Tuple[float, float, float]]]] = {}
    bucket_colors: Dict[str, Tuple[float, float, float]] = {}

    # Default color (yapCAD blue)
    default_color = (0.6, 0.85, 1.0)

    surfaces: Dict[str, Dict[str, any]] = {}
    for entry in doc.get("entities", []):
        if entry.get("type") == "surface":
            surfaces[entry["id"]] = entry

    for entry in doc.get("entities", []):
        if entry.get("type") not in {"solid", "surface"}:
            continue
        meta = entry.get("metadata", {})
        material_id = meta.get("material")
        layer = meta.get("layer", "default")

        # Determine bucket key and color
        if material_id and material_id in materials:
            bucket_key = f"mat:{material_id}"
            mat_def = materials[material_id]
            visual = mat_def.get("visual", {})
            color = visual.get("color", default_color)
            if isinstance(color, list) and len(color) >= 3:
                bucket_colors[bucket_key] = (float(color[0]), float(color[1]), float(color[2]))
            else:
                bucket_colors[bucket_key] = default_color
        else:
            bucket_key = f"layer:{layer}"
            if bucket_key not in bucket_colors:
                bucket_colors[bucket_key] = default_color

        if entry.get("type") == "surface":
            surf_entry = entry
            verts = surf_entry.get("vertices", [])
            norms = surf_entry.get("normals", [])
            faces = surf_entry.get("faces", [])
            bucket = bucket_tris.setdefault(bucket_key, [])
            for tri in faces:
                for idx in tri:
                    pt = verts[idx]
                    normal = norms[idx] if idx < len(norms) else [0.0, 0.0, 1.0, 0.0]
                    bucket.append(((float(pt[0]), float(pt[1]), float(pt[2])), (float(normal[0]), float(normal[1]), float(normal[2]))))
            continue

        shell_ids = entry.get("shell", [])

        # If shell is empty but BREP data exists, tessellate the BREP
        if not shell_ids and meta.get("brep"):
            brep_data = meta.get("brep", {})
            if brep_data.get("encoding") == "brep-ascii-base64" and brep_data.get("data"):
                try:
                    tris = _tessellate_brep(brep_data["data"])
                    if tris:
                        bucket = bucket_tris.setdefault(bucket_key, [])
                        bucket.extend(tris)
                except Exception as e:
                    print(f"Warning: BREP tessellation failed: {e}")
            continue

        for sid in shell_ids:
            surf_entry = surfaces.get(sid)
            if not surf_entry:
                continue
            verts = surf_entry.get("vertices", [])
            norms = surf_entry.get("normals", [])
            faces = surf_entry.get("faces", [])
            bucket = bucket_tris.setdefault(bucket_key, [])
            for tri in faces:
                for idx in tri:
                    pt = verts[idx]
                    normal = norms[idx] if idx < len(norms) else [0.0, 0.0, 1.0, 0.0]
                    bucket.append(((float(pt[0]), float(pt[1]), float(pt[2])), (float(normal[0]), float(normal[1]), float(normal[2]))))

    if not bucket_tris:
        return bucket_tris, bucket_colors, (0, 0, 0, 0, 0, 0)
    # List comprehensions + min/max are faster than explicit loops in Python
    # due to C-level optimization of built-ins
    xs = [v[0][0] for tris in bucket_tris.values() for v in tris]
    ys = [v[0][1] for tris in bucket_tris.values() for v in tris]
    zs = [v[0][2] for tris in bucket_tris.values() for v in tris]
    return bucket_tris, bucket_colors, (min(xs), min(ys), min(zs), max(xs), max(ys), max(zs))


def _collect_polylines(doc: Dict[str, any]) -> Tuple[Dict[str, List[List[List[float]]]], Tuple[float, float, float, float]]:
    sketches = [entry for entry in doc.get("entities", []) if entry.get("type") == "sketch"]
    layer_polys: Dict[str, List[List[List[float]]]] = {}
    xs: List[float] = []
    ys: List[float] = []
    for entry in sketches:
        layer = entry.get("metadata", {}).get("layer", "default")
        bucket = layer_polys.setdefault(layer, [])
        for poly in entry.get("polylines", []):
            points = []
            for pt in poly:
                x, y = float(pt[0]), float(pt[1])
                points.append([x, y])
                xs.append(x)
                ys.append(y)
            if points:
                bucket.append(points)
    if not layer_polys:
        return layer_polys, (0, 0, 0, 0)
    return layer_polys, (min(xs), min(ys), max(xs), max(ys))


def _compute_grid_step(span: float) -> float:
    if span <= 0 or math.isnan(span):
        return 1.0
    raw = span / 10.0
    if raw <= 0:
        raw = span or 1.0
    exp = math.floor(math.log10(raw)) if raw > 0 else 0
    base = 10 ** exp
    for factor in (1, 2, 5, 10):
        step = factor * base
        if raw <= step:
            return step
    return 10 * base


class FourViewWindow(pyglet.window.Window):
    """Four-view 3D window with perspective + orthographic panels."""

    def __init__(self, bucket_triangles: Dict[str, List[Tuple[Tuple[float, float, float], Tuple[float, float, float]]]], bucket_colors: Dict[str, Tuple[float, float, float]], bbox: Tuple[float, float, float, float, float, float], units: str = ""):
        super().__init__(width=1200, height=800, caption="yapCAD Package Viewer")
        self.bucket_triangles = bucket_triangles
        self.bucket_colors = bucket_colors
        self.bucket_names = sorted(bucket_triangles.keys()) or ["default"]
        self.visible_buckets = {bucket: True for bucket in self.bucket_names}
        # For display, extract human-readable names
        self.layer_names = [name.split(":", 1)[1] if ":" in name else name for name in self.bucket_names]
        self.visible_layers = {self.layer_names[i]: self.visible_buckets[self.bucket_names[i]] for i in range(len(self.bucket_names))}
        self.bbox = bbox
        self.units = units
        self.show_help = False
        self.azimuth = 35.0
        self.elevation = 25.0
        self.distance = max(bbox[3] - bbox[0], bbox[4] - bbox[1], bbox[5] - bbox[2]) * 1.5 or 10.0
        self.pan_x = 0.0
        self.pan_y = 0.0
        self._dragging = False
        self._last = (0, 0)
        self.render_mode = 0  # 0 = solid, 1 = solid+mesh, 2 = mesh only
        self.view_mode = 0  # 0 = quad, 1 = perspective, 2 = front, 3 = top, 4 = side
        self._view_mode_names = ["Quad", "Perspective", "Front", "Top", "Side"]
        self._bucket_vertex_lists: Dict[str, graphics.vertexdomain.VertexList] = {}
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_LIGHTING)
        glEnable(GL_LIGHT0)
        glEnable(GL_COLOR_MATERIAL)
        glColorMaterial(pyglet.gl.GL_FRONT_AND_BACK, pyglet.gl.GL_AMBIENT_AND_DIFFUSE)
        light_position = (GLfloat * 4)(0.6, 0.8, 1.2, 0.0)
        glLightfv(GL_LIGHT0, pyglet.gl.GL_POSITION, light_position)
        light_diffuse = (GLfloat * 4)(0.8, 0.8, 0.8, 1.0)
        glLightfv(GL_LIGHT0, pyglet.gl.GL_DIFFUSE, light_diffuse)
        light_ambient = (GLfloat * 4)(0.15, 0.15, 0.15, 1.0)
        glLightfv(GL_LIGHT0, pyglet.gl.GL_AMBIENT, light_ambient)
        glClearColor(0.05, 0.05, 0.07, 1.0)
        self._build_triangle_cache()

    def _build_triangle_cache(self) -> None:
        for vlist in self._bucket_vertex_lists.values():
            vlist.delete()
        self._bucket_vertex_lists = {}
        for bucket in self.bucket_names:
            tris = self.bucket_triangles.get(bucket, [])
            if not tris:
                continue
            coords: List[float] = []
            normals: List[float] = []
            for vertex, normal in tris:
                coords.extend((vertex[0], vertex[1], vertex[2]))
                normals.extend((normal[0], normal[1], normal[2]))
            if coords:
                self._bucket_vertex_lists[bucket] = graphics.vertex_list(
                    len(tris),
                    ('v3f/static', coords),
                    ('n3f/static', normals),
                )

    def on_draw(self):
        self.clear()
        fb_width, fb_height = self.get_framebuffer_size()
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

        if self.view_mode == 0:  # Quad view
            w2 = fb_width // 2
            h2 = fb_height // 2
            self._draw_viewport(0, h2, w2, h2, "Perspective", perspective=True, fb_dims=(fb_width, fb_height))
            self._draw_viewport(w2, h2, w2, h2, "Front", orientation="front", fb_dims=(fb_width, fb_height))
            self._draw_viewport(0, 0, w2, h2, "Top", orientation="top", fb_dims=(fb_width, fb_height))
            self._draw_viewport(w2, 0, w2, h2, "Side", orientation="side", fb_dims=(fb_width, fb_height))
        elif self.view_mode == 1:  # Perspective only
            self._draw_viewport(0, 0, fb_width, fb_height, "Perspective", perspective=True, fb_dims=(fb_width, fb_height))
        elif self.view_mode == 2:  # Front only
            self._draw_viewport(0, 0, fb_width, fb_height, "Front", orientation="front", fb_dims=(fb_width, fb_height))
        elif self.view_mode == 3:  # Top only
            self._draw_viewport(0, 0, fb_width, fb_height, "Top", orientation="top", fb_dims=(fb_width, fb_height))
        elif self.view_mode == 4:  # Side only
            self._draw_viewport(0, 0, fb_width, fb_height, "Side", orientation="side", fb_dims=(fb_width, fb_height))

        self._draw_help_overlay(fb_width, fb_height)

    def _apply_camera(self, orientation: str | None, width: int, height: int, perspective: bool):
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        aspect = max(width / height, 0.1)
        if perspective:
            gluPerspective(45.0, aspect, 0.1, 1000.0)
        else:
            span = self.distance
            glOrtho(-span * aspect, span * aspect, -span, span, -1000.0, 1000.0)

        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        cx = (self.bbox[0] + self.bbox[3]) / 2.0
        cy = (self.bbox[1] + self.bbox[4]) / 2.0
        cz = (self.bbox[2] + self.bbox[5]) / 2.0

        if perspective:
            theta = math.radians(self.azimuth)
            phi = math.radians(self.elevation)
            eye_x = cx + self.distance * math.cos(phi) * math.cos(theta)
            eye_y = cy + self.distance * math.sin(phi)
            eye_z = cz + self.distance * math.cos(phi) * math.sin(theta)
            gluLookAt(eye_x, eye_y, eye_z, cx + self.pan_x, cy + self.pan_y, cz, 0, 1, 0)
        elif orientation == "front":
            gluLookAt(cx, cy, cz + self.distance, cx + self.pan_x, cy + self.pan_y, cz, 0, 1, 0)
        elif orientation == "top":
            gluLookAt(cx, cy + self.distance, cz, cx + self.pan_x, cy, cz + self.pan_y, 0, 0, -1)
        elif orientation == "side":
            gluLookAt(cx + self.distance, cy, cz, cx, cy + self.pan_y, cz + self.pan_x, 0, 1, 0)
        else:
            gluLookAt(cx, cy, cz + self.distance, cx, cy, cz, 0, 1, 0)

    def _draw_triangles(self):
        mode = self.render_mode

        if mode in (0, 1):
            glEnable(GL_LIGHTING)
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
            for bucket in self.bucket_names:
                if not self.visible_buckets.get(bucket, True):
                    continue
                vlist = self._bucket_vertex_lists.get(bucket)
                if vlist:
                    # Set color based on material
                    color = self.bucket_colors.get(bucket, (0.6, 0.85, 1.0))
                    glColor4f(color[0], color[1], color[2], 1.0)
                    vlist.draw(GL_TRIANGLES)

        if mode in (1, 2):
            glDisable(GL_LIGHTING)
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)
            glEnable(GL_POLYGON_OFFSET_LINE)
            glPolygonOffset(-1.0, 1.0)
            glColor4f(1.0, 1.0, 1.0, 1.0)
            glLineWidth(1.0)
            for bucket in self.bucket_names:
                if not self.visible_buckets.get(bucket, True):
                    continue
                vlist = self._bucket_vertex_lists.get(bucket)
                if vlist:
                    vlist.draw(GL_TRIANGLES)
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
            glDisable(GL_POLYGON_OFFSET_LINE)
            glEnable(GL_LIGHTING)

    def _draw_grid(self, orientation: str | None, perspective: bool):
        xmin, ymin, zmin, xmax, ymax, zmax = self.bbox
        span_x = max(xmax - xmin, 1e-6)
        span_y = max(ymax - ymin, 1e-6)
        span_z = max(zmax - zmin, 1e-6)
        cx = (xmin + xmax) / 2.0
        cy = (ymin + ymax) / 2.0
        cz = (zmin + zmax) / 2.0

        if perspective:
            plane_z = zmin
            step = _compute_grid_step(max(span_x, span_y))
            min_x = (math.floor(xmin / step) - 2) * step
            max_x = (math.ceil(xmax / step) + 2) * step
            min_y = (math.floor(ymin / step) - 2) * step
            max_y = (math.ceil(ymax / step) + 2) * step
            glDisable(GL_LIGHTING)
            glColor4f(0.25, 0.25, 0.3, 0.7)
            glBegin(GL_LINES)
            value = min_x
            while value <= max_x:
                glVertex3f(value, min_y, plane_z)
                glVertex3f(value, max_y, plane_z)
                value += step
            value = min_y
            while value <= max_y:
                glVertex3f(min_x, value, plane_z)
                glVertex3f(max_x, value, plane_z)
                value += step
            glEnd()
            glEnable(GL_LIGHTING)
            return

        glDisable(GL_LIGHTING)
        glColor4f(0.2, 0.2, 0.28, 0.8)
        glBegin(GL_LINES)
        if orientation == "front":
            step_x = step_y = _compute_grid_step(max(span_x, span_y))
            min_x = (math.floor(xmin / step_x) - 2) * step_x
            max_x = (math.ceil(xmax / step_x) + 2) * step_x
            min_y = (math.floor(ymin / step_y) - 2) * step_y
            max_y = (math.ceil(ymax / step_y) + 2) * step_y
            value = min_x
            while value <= max_x:
                glVertex3f(value, min_y, cz)
                glVertex3f(value, max_y, cz)
                value += step_x
            value = min_y
            while value <= max_y:
                glVertex3f(min_x, value, cz)
                glVertex3f(max_x, value, cz)
                value += step_y
        elif orientation == "top":
            step_x = step_z = _compute_grid_step(max(span_x, span_z))
            min_x = (math.floor(xmin / step_x) - 2) * step_x
            max_x = (math.ceil(xmax / step_x) + 2) * step_x
            min_z = (math.floor(zmin / step_z) - 2) * step_z
            max_z = (math.ceil(zmax / step_z) + 2) * step_z
            value = min_x
            while value <= max_x:
                glVertex3f(value, cy, min_z)
                glVertex3f(value, cy, max_z)
                value += step_x
            value = min_z
            while value <= max_z:
                glVertex3f(min_x, cy, value)
                glVertex3f(max_x, cy, value)
                value += step_z
        elif orientation == "side":
            step_y = step_z = _compute_grid_step(max(span_y, span_z))
            min_y = (math.floor(ymin / step_y) - 2) * step_y
            max_y = (math.ceil(ymax / step_y) + 2) * step_y
            min_z = (math.floor(zmin / step_z) - 2) * step_z
            max_z = (math.ceil(zmax / step_z) + 2) * step_z
            value = min_y
            while value <= max_y:
                glVertex3f(cx, value, min_z)
                glVertex3f(cx, value, max_z)
                value += step_y
            value = min_z
            while value <= max_z:
                glVertex3f(cx, min_y, value)
                glVertex3f(cx, max_y, value)
                value += step_z
        glEnd()
        glEnable(GL_LIGHTING)

    def _draw_viewport(self, x: int, y: int, width: int, height: int, label: str, perspective: bool = False, orientation: str | None = None, fb_dims: Tuple[int, int] | None = None):
        fb_width, fb_height = fb_dims if fb_dims else (self.width, self.height)
        glViewport(x, y, width, height)
        glEnable(GL_SCISSOR_TEST)
        glScissor(x, y, width, height)
        self._apply_camera(orientation, width, height, perspective)
        self._draw_grid(orientation, perspective)
        self._draw_triangles()
        glDisable(GL_LIGHTING)
        label_text = label
        if self.units:
            label_text = f"{label} ({self.units})"
        self._draw_axes_overlay(x, y, width, height, label_text, fb_width, fb_height, orientation, perspective)
        glEnable(GL_LIGHTING)
        glDisable(GL_SCISSOR_TEST)

    def _draw_axes_overlay(self, x: int, y: int, width: int, height: int, text: str, fb_width: int, fb_height: int, orientation: str | None, perspective: bool):
        screen_x = x + max(12, width // 30)
        screen_y = y + max(12, height // 30)
        glMatrixMode(GL_PROJECTION)
        glPushMatrix()
        glLoadIdentity()
        glOrtho(0, fb_width, 0, fb_height, -1, 1)
        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()
        glLoadIdentity()
        font_px = max(16, min(width, height) // 10)
        label = pyglet.text.Label(
            text,
            font_size=font_px,
            x=screen_x,
            y=screen_y,
            anchor_x="left",
            anchor_y="bottom",
            color=(255, 255, 255, 220),
        )
        label.draw()
        axis_lines: List[str] = []
        if perspective:
            pass
            #axis_lines.append("+X →, +Y ↑, +Z out")
        elif orientation == "front":
            axis_lines.append("+X →, +Y ↑")
        elif orientation == "top":
            axis_lines.append("+X →, +Z ↑")
        elif orientation == "side":
            axis_lines.append("+Y →, +Z ↑")

        if len(self.layer_names) > 1:
            layer_display = []
            for idx, layer in enumerate(self.layer_names, start=1):
                state = "ON" if self.visible_layers.get(layer, True) else "off"
                layer_display.append(f"{idx}:{layer}({state})")
            axis_lines.append("Layers " + ", ".join(layer_display))
        elif self.layer_names:
            axis_lines.append(f"Layer {self.layer_names[0]}")

        y_offset = font_px + 6
        for line in axis_lines:
            axis_label = pyglet.text.Label(
                line,
                font_size=max(12, font_px - 2),
                x=screen_x,
                y=screen_y + y_offset,
                anchor_x="left",
                anchor_y="bottom",
                color=(200, 200, 220, 200),
            )
            axis_label.draw()
            y_offset += axis_label.content_height + 2
        glMatrixMode(GL_MODELVIEW)
        glPopMatrix()
        glMatrixMode(GL_PROJECTION)
        glPopMatrix()
        glMatrixMode(GL_MODELVIEW)

    def _draw_help_overlay(self, fb_width: int, fb_height: int) -> None:
        if not self.show_help:
            return
        glPushAttrib(GL_ENABLE_BIT | GL_DEPTH_BUFFER_BIT)
        try:
            glDisable(GL_LIGHTING)
            glDisable(GL_DEPTH_TEST)
            glMatrixMode(GL_PROJECTION)
            glPushMatrix()
            glLoadIdentity()
            glOrtho(0, fb_width, 0, fb_height, -1, 1)
            glMatrixMode(GL_MODELVIEW)
            glPushMatrix()
            glLoadIdentity()

            # margin = 40
            margin = fb_width // 20
            h_margin = fb_height // 30
            panel_width = fb_width - 2 * margin
            panel_height = min(fb_height // 1.5 , fb_height - 2 * h_margin)
            left = margin
            bottom = fb_height - h_margin - panel_height

            glColor4f(0.05, 0.05, 0.08, 0.9)
            pyglet.graphics.draw(
                4,
                GL_QUADS,
                ("v2f", [
                    left, bottom,
                    left + panel_width, bottom,
                    left + panel_width, bottom + panel_height,
                    left, bottom + panel_height,
                ]),
            )

            active_layers = ", ".join(layer for layer, vis in self.visible_layers.items() if vis) or "none"
            view_name = self._view_mode_names[self.view_mode]
            help_lines = [
                "Viewer Controls",
                f"Current view: {view_name}",
                "Navigation:",
                "  Left drag – rotate (perspective views)",
                "  Right drag – pan",
                "  Scroll/Swipe up – zoom out, down – zoom in",
                "View Modes:",
                "  V – cycle views (Quad → Perspective → Front → Top → Side)",
                "  Tab – cycle single views only",
                "Layers:",
                "  Number keys 1-9 toggle layers, 0 resets",
                f"  Active: {active_layers}",
                "General:",
                "  H or F1 – help, ESC – close, L – cycle shading",
            ]

#            title_font = max(28, min(fb_width, fb_height) // 13)
#            body_font = max(20, title_font // 1.5)
            title_font = max(28, min(fb_width, fb_height) // 23)
            body_font = max(20, title_font // 1.5)
            y = bottom + panel_height - (h_margin + title_font)
            for idx, line in enumerate(help_lines):
                size = title_font if idx == 0 else body_font
                label = pyglet.text.Label(
                    line,
                    font_size=int(size),
                    x=left + 26,
                    y=y,
                    anchor_x="left",
                    anchor_y="baseline",
                    color=(255, 255, 255, 235),
                )
                label.draw()
                y -= label.content_height + 10

            glMatrixMode(GL_MODELVIEW)
            glPopMatrix()
            glMatrixMode(GL_PROJECTION)
            glPopMatrix()
            glMatrixMode(GL_MODELVIEW)
        finally:
            glPopAttrib()

    def on_mouse_drag(self, x, y, dx, dy, buttons, modifiers):
        if buttons & pyglet.window.mouse.LEFT:
            # In quad mode, only rotate when in perspective quadrant (top-left)
            # In single perspective mode, rotate anywhere
            if self.view_mode == 1 or (self.view_mode == 0 and x < self.width // 2 and y > self.height // 2):
                self.azimuth += dx * 0.5
                self.elevation = max(-89.0, min(89.0, self.elevation - dy * 0.5))
        elif buttons & pyglet.window.mouse.RIGHT:
            self.pan_x += dx * 0.01
            self.pan_y += dy * 0.01

    def on_mouse_scroll(self, x, y, scroll_x, scroll_y):
        if scroll_y > 0:
            self.distance = max(0.5, self.distance * 1.1)
        elif scroll_y < 0:
            self.distance = max(0.5, self.distance * 0.9)

    def on_key_press(self, symbol, modifiers):
        if symbol == key.ESCAPE:
            self.close()
        elif key._1 <= symbol <= key._9:
            idx = symbol - key._1
            if idx < len(self.bucket_names):
                bucket = self.bucket_names[idx]
                self.visible_buckets[bucket] = not self.visible_buckets.get(bucket, True)
                # Keep display names in sync
                self.visible_layers[self.layer_names[idx]] = self.visible_buckets[bucket]
        elif symbol == key._0:
            for bucket in self.bucket_names:
                self.visible_buckets[bucket] = True
            for layer in self.layer_names:
                self.visible_layers[layer] = True
        elif symbol in (key.H, key.F1):
            self.show_help = not self.show_help
        elif symbol == key.L:
            self.render_mode = (self.render_mode + 1) % 3
        elif symbol == key.V:
            self.view_mode = (self.view_mode + 1) % 5
        elif symbol == key.TAB:
            # Tab cycles through single views only (1-4), skipping quad
            if self.view_mode == 0:
                self.view_mode = 1
            else:
                self.view_mode = (self.view_mode % 4) + 1

    def on_close(self):
        for vlist in self._bucket_vertex_lists.values():
            vlist.delete()
        self._bucket_vertex_lists.clear()
        super().on_close()
#        else:
#            print(f"Key pressed: {key.symbol_string(symbol)}")


class SketchWindow(pyglet.window.Window):
    """2D viewer for sketch entities."""

    def __init__(self, layer_polylines: Dict[str, List[List[List[float]]]], bounds: Tuple[float, float, float, float], units: str = ""):
        super().__init__(width=900, height=700, caption="yapCAD Package Viewer (2D)")
        self.layer_polylines = layer_polylines
        self.layer_names = sorted(layer_polylines.keys()) or ["default"]
        self.active_layers = {layer: True for layer in self.layer_names}
        self.bounds = bounds
        self.units = units
        self.offset_x = -(bounds[0] + bounds[2]) / 2.0
        self.offset_y = -(bounds[1] + bounds[3]) / 2.0
        span_x = bounds[2] - bounds[0]
        span_y = bounds[3] - bounds[1]
        self.scale = 1.8 / max(span_x, span_y, 1.0)
        self.pan = [0.0, 0.0]
        self.dragging = False
        self.drag_start = (0, 0)
        self.show_help = False
        glClearColor(0.05, 0.05, 0.07, 1.0)
        glEnable(pyglet.gl.GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)

    def _draw_grid(self):
        xmin, ymin, xmax, ymax = self.bounds
        span_x = max(xmax - xmin, 1e-6)
        span_y = max(ymax - ymin, 1e-6)
        step = _compute_grid_step(max(span_x, span_y))
        min_x = (math.floor(xmin / step) - 5) * step
        max_x = (math.ceil(xmax / step) + 5) * step
        min_y = (math.floor(ymin / step) - 5) * step
        max_y = (math.ceil(ymax / step) + 5) * step
        pyglet.gl.glColor4f(0.2, 0.2, 0.28, 0.5)
        pyglet.gl.glBegin(pyglet.gl.GL_LINES)
        value = min_x
        while value <= max_x:
            pyglet.gl.glVertex2f(value, min_y)
            pyglet.gl.glVertex2f(value, max_y)
            value += step
        value = min_y
        while value <= max_y:
            pyglet.gl.glVertex2f(min_x, value)
            pyglet.gl.glVertex2f(max_x, value)
            value += step
        pyglet.gl.glEnd()

    def on_draw(self):
        self.clear()
        fb_width, fb_height = self.get_framebuffer_size()
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        glOrtho(-fb_width / 2, fb_width / 2, -fb_height / 2, fb_height / 2, -1, 1)
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        glScalef(self.scale * fb_width / 2, self.scale * fb_height / 2, 1)
        glTranslatef(self.pan[0] + self.offset_x, self.pan[1] + self.offset_y, 0)

        self._draw_grid()

        pyglet.gl.glColor4f(0.6, 0.85, 1.0, 1.0)
        for layer in self.layer_names:
            if not self.active_layers.get(layer, True):
                continue
            for poly in self.layer_polylines.get(layer, []):
                if len(poly) < 2:
                    continue
                closed = poly[0] == poly[-1]
                if closed:
                    glBegin(GL_LINE_LOOP)
                    for pt in poly[:-1]:
                        glVertex2f(pt[0], pt[1])
                else:
                    glBegin(GL_LINES)
                    for i in range(len(poly) - 1):
                        x0, y0 = poly[i]
                        x1, y1 = poly[i + 1]
                        glVertex2f(x0, y0)
                        glVertex2f(x1, y1)
                glEnd()

        self._draw_overlay()
        self._draw_help_overlay(fb_width, fb_height)

    def _draw_overlay(self):
        glMatrixMode(GL_PROJECTION)
        glPushMatrix()
        glLoadIdentity()
        fb_width, fb_height = self.get_framebuffer_size()
        glOrtho(0, fb_width, 0, fb_height, -1, 1)
        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()
        glLoadIdentity()

        overlay_text = "Press H for help"
        if self.units:
            overlay_text += f" • Units: {self.units}"
        if len(self.layer_names) > 1:
            layer_display = []
            for idx, layer in enumerate(self.layer_names, start=1):
                state = "ON" if self.active_layers.get(layer, True) else "off"
                layer_display.append(f"{idx}:{layer}({state})")
            overlay_text += " • Layers " + ", ".join(layer_display)
        elif self.layer_names:
            overlay_text += f" • Layer {self.layer_names[0]}"
        info = pyglet.text.Label(
            overlay_text,
            font_size=12,
            x=16,
            y=16,
            anchor_x="left",
            anchor_y="bottom",
            color=(255, 255, 255, 200),
        )
        info.draw()

        glMatrixMode(GL_MODELVIEW)
        glPopMatrix()
        glMatrixMode(GL_PROJECTION)
        glPopMatrix()
        glMatrixMode(GL_MODELVIEW)

    def _draw_help_overlay(self, fb_width: int, fb_height: int) -> None:
        if not self.show_help:
            return
        glMatrixMode(GL_PROJECTION)
        glPushMatrix()
        glLoadIdentity()
        glOrtho(0, fb_width, 0, fb_height, -1, 1)
        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()
        glLoadIdentity()

        margin = fb_width // 20
        h_margin = fb_height // 30
        panel_width = fb_width/2 - 2 * margin
        panel_height = min(fb_height // 1.5 , fb_height // 3 - h_margin)
        left = margin
        bottom = fb_height - margin - panel_height

        pyglet.gl.glColor4f(0.08, 0.08, 0.1, 0.9)
        pyglet.graphics.draw(
            4,
            GL_QUADS,
            ("v2f", [
                left, bottom,
                left + panel_width, bottom,
                left + panel_width, bottom + panel_height,
                left, bottom + panel_height,
            ]),
        )

        help_lines = [
            "Sketch Viewer Controls",
            "  Scroll/Swipe up – zoom out",
            "  Scroll/Swipe down – zoom in",
            "  Right drag – pan",
            "  1-9 – toggle layers (0 resets)",
            "  H or F1 – toggle help, ESC – close",
        ]
        title_font = max(36, min(fb_width, fb_height) // 30)
        body_font = max(24, title_font // 1.5)
        y = bottom + panel_height - (h_margin + title_font)
        for idx, line in enumerate(help_lines):
            size = title_font if idx == 0 else body_font
            label = pyglet.text.Label(
                line,
                font_size=int(size),
                x=left + 22,
                y=y,
                anchor_x="left",
                anchor_y="baseline",
                color=(255, 255, 255, 230),
            )
            label.draw()
            y -= label.content_height + 10

        glMatrixMode(GL_MODELVIEW)
        glPopMatrix()
        glMatrixMode(GL_PROJECTION)
        glPopMatrix()
        glMatrixMode(GL_MODELVIEW)

    def on_mouse_scroll(self, x, y, scroll_x, scroll_y):
        if scroll_y > 0:
            self.scale *= 0.9
        elif scroll_y < 0:
            self.scale *= 1.1

    def on_mouse_drag(self, x, y, dx, dy, buttons, modifiers):
        if buttons & pyglet.window.mouse.RIGHT:
            self.pan[0] += dx / (self.width * self.scale)
            self.pan[1] += dy / (self.height * self.scale)

    def on_key_press(self, symbol, modifiers):
        if symbol == key.ESCAPE:
            self.close()
        elif pyglet.window.key._1 <= symbol <= pyglet.window.key._9:
            idx = symbol - pyglet.window.key._1
            if idx < len(self.layer_names):
                layer = self.layer_names[idx]
                self.active_layers[layer] = not self.active_layers.get(layer, True)
        elif symbol == pyglet.window.key._0:
            for layer in self.layer_names:
                self.active_layers[layer] = True
        elif symbol in (key.H, key.F1):
            self.show_help = not self.show_help
#        else:
#            print(f"Key pressed: {key.symbol_string(symbol)}")


def view_package(package_path: Path | str, *, strict: bool = False) -> bool:
    """Validate a package and launch appropriate viewer."""
    pkg_path = Path(package_path)
    ok, messages = validate_package(pkg_path, strict=strict)
    for msg in messages:
        print(msg)
    if not ok:
        return False

    manifest = PackageManifest.load(pkg_path)
    doc = _load_geometry_doc(manifest.root, manifest.data)
    try:
        geometry_from_json(doc)
    except Exception as exc:
        print(f"Failed to load geometry: {exc}")
        return False

    units = manifest.data.get("units", "")
    materials = manifest.get_materials()
    bucket_tris, bucket_colors, bbox = _collect_triangles(doc, materials)
    print(f"bounding box: {bbox}")
    if bucket_tris:
        window = FourViewWindow(bucket_tris, bucket_colors, bbox, units=units)
    else:
        layer_polys, bounds = _collect_polylines(doc)
        if not layer_polys:
            print("No geometry to display.")
            return False
        window = SketchWindow(layer_polys, bounds, units=units)

    pyglet.app.run()
    return True


__all__ = ["view_package"]
