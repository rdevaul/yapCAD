"""
Profile evaluation endpoints — parametric 2D curve profiles for the workbench.

GET  /profiles                    list available profiles
GET  /profiles/{id}/schema        parameter schema (derived from Python signature)
POST /profiles/{id}/eval          evaluate profile, return points + metadata
"""
from __future__ import annotations

import inspect
import sys
from typing import Any, Dict, List, Optional

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

router = APIRouter(prefix="/profiles", tags=["profiles"])

# ── Workspace profile path ────────────────────────────────────────────────────
_WORKSPACE_PROJECTS = "/Users/jarvis/.openclaw/workspace-jarvis/projects"


def _load_tank_dome():
    """Import tank_dome_profile from the workspace, caching the result."""
    if _WORKSPACE_PROJECTS not in sys.path:
        sys.path.insert(0, _WORKSPACE_PROJECTS)
    from tank_dome_profile import tank_dome_profile  # type: ignore
    return tank_dome_profile


# ── Profile registry ──────────────────────────────────────────────────────────

_PROFILES = {
    "tank_dome": {
        "id": "tank_dome",
        "label": "Tank Dome Profile",
        "description": "Parametric sigmoid cubic-Bézier dome for rocket tank cross-section",
        "loader": _load_tank_dome,
        # Hardcoded range hints (supplement inspect.signature defaults)
        "ranges": {
            "neck_radius":     {"min": 5.0,  "max": 200.0, "step": 0.5,  "unit": "mm"},
            "tank_radius":     {"min": 20.0, "max": 500.0, "step": 0.5,  "unit": "mm"},
            "total_height":    {"min": 20.0, "max": 500.0, "step": 1.0,  "unit": "mm"},
            "lead_in_length":  {"min": 0.0,  "max": 100.0, "step": 0.5,  "unit": "mm"},
            "lead_out_length": {"min": 0.0,  "max": 100.0, "step": 0.5,  "unit": "mm"},
            "tension":         {"min": 0.0,  "max": 1.0,   "step": 0.05, "unit": ""},
            "max_overhang":    {"min": 10.0, "max": 90.0,  "step": 1.0,  "unit": "deg",
                                "optional": True},
        },
    }
}


# ── Models ────────────────────────────────────────────────────────────────────

class ProfileEntry(BaseModel):
    id: str
    label: str
    description: str


class ParamSchema(BaseModel):
    name: str
    type: str
    default: Optional[Any]
    min: Optional[float]
    max: Optional[float]
    step: Optional[float]
    unit: str
    optional: bool


class ProfileSchema(BaseModel):
    id: str
    label: str
    params: List[ParamSchema]


class EvalRequest(BaseModel):
    params: Dict[str, Any]


class EvalResponse(BaseModel):
    points: List[List[float]]   # [[r, z], ...]
    metadata: Dict[str, Any]


# ── Helpers ───────────────────────────────────────────────────────────────────

def _build_schema(profile_id: str) -> ProfileSchema:
    entry = _PROFILES[profile_id]
    fn = entry["loader"]()
    sig = inspect.signature(fn)
    ranges = entry.get("ranges", {})

    params: List[ParamSchema] = []
    for name, param in sig.parameters.items():
        if name in ("args", "kwargs"):
            continue
        default = param.default if param.default is not inspect.Parameter.empty else None
        type_name = "float"
        if param.annotation is not inspect.Parameter.empty:
            ann = param.annotation
            if ann is int:
                type_name = "int"
            elif ann is bool:
                type_name = "bool"
            elif ann is str:
                type_name = "str"
            else:
                type_name = "float"

        r = ranges.get(name, {})
        params.append(ParamSchema(
            name=name,
            type=type_name,
            default=default,
            min=r.get("min"),
            max=r.get("max"),
            step=r.get("step"),
            unit=r.get("unit", ""),
            optional=r.get("optional", False),
        ))

    return ProfileSchema(
        id=profile_id,
        label=entry["label"],
        params=params,
    )


def _eval_profile(profile_id: str, params: Dict[str, Any]) -> EvalResponse:
    entry = _PROFILES[profile_id]
    fn = entry["loader"]()
    sig = inspect.signature(fn)

    # Fill defaults for omitted optional params
    call_kwargs: Dict[str, Any] = {}
    for name, param in sig.parameters.items():
        if name in params:
            call_kwargs[name] = params[name]
        elif param.default is not inspect.Parameter.empty:
            call_kwargs[name] = param.default
        # else: required param missing — let Python raise TypeError

    result = fn(**call_kwargs)

    # tank_dome_profile returns (points, metadata) or just points
    if isinstance(result, tuple):
        points_raw, meta = result
    else:
        points_raw = result
        meta = {}

    # Normalise points to [[r, z], ...]
    points: List[List[float]] = []
    for pt in points_raw:
        if hasattr(pt, "__iter__"):
            coords = list(pt)
            points.append([float(coords[0]), float(coords[1])])
        else:
            points.append([float(pt), 0.0])

    # Build metadata from call_kwargs + any returned meta
    metadata: Dict[str, Any] = {
        "neck_radius":  float(call_kwargs.get("neck_radius", 0)),
        "tank_radius":  float(call_kwargs.get("tank_radius", 0)),
        "total_height": float(call_kwargs.get("total_height", 0)),
    }
    metadata.update({k: v for k, v in meta.items() if isinstance(v, (int, float, bool, str))})

    # Compute max overhang angle from profile points if not in meta
    if "max_overhang_deg" not in metadata and len(points) >= 2:
        import math
        max_angle = 0.0
        for i in range(1, len(points)):
            dr = points[i][0] - points[i - 1][0]
            dz = points[i][1] - points[i - 1][1]
            if abs(dz) > 1e-9:
                angle = abs(math.degrees(math.atan2(abs(dr), abs(dz))))
                max_angle = max(max_angle, angle)
        metadata["max_overhang_deg"] = round(max_angle, 1)

    max_overhang_param = call_kwargs.get("max_overhang", 45.0)
    metadata["overhang_ok"] = metadata.get("max_overhang_deg", 0) <= (max_overhang_param or 45.0)

    return EvalResponse(points=points, metadata=metadata)


# ── Routes ────────────────────────────────────────────────────────────────────

@router.get("", response_model=List[ProfileEntry])
async def list_profiles():
    return [
        ProfileEntry(id=v["id"], label=v["label"], description=v["description"])
        for v in _PROFILES.values()
    ]


@router.get("/{profile_id}/schema", response_model=ProfileSchema)
async def get_profile_schema(profile_id: str):
    if profile_id not in _PROFILES:
        raise HTTPException(404, f"Profile '{profile_id}' not found")
    try:
        return _build_schema(profile_id)
    except Exception as e:
        raise HTTPException(500, f"Schema build failed: {e}")


@router.post("/{profile_id}/eval", response_model=EvalResponse)
async def eval_profile(profile_id: str, body: EvalRequest):
    if profile_id not in _PROFILES:
        raise HTTPException(404, f"Profile '{profile_id}' not found")
    try:
        return _eval_profile(profile_id, body.params)
    except TypeError as e:
        raise HTTPException(422, f"Parameter error: {e}")
    except Exception as e:
        raise HTTPException(500, f"Evaluation failed: {e}")
