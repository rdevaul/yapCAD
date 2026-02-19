from __future__ import annotations

import os
from typing import Any, Dict, List

from fastapi import APIRouter, HTTPException

from yapcad.boolean import get_engine
from yapcad.io.geometry_json import geometry_from_json, geometry_to_json
from yapcad.geom3d import issolid

from ..core.tessellator import solidish_to_mesh
from ..models.schemas import (
    BooleanRequest,
    BooleanResponse,
    MeshResponse,
    TessellateRequest,
)

router = APIRouter(tags=["geometry"])


@router.post("/tessellate", response_model=MeshResponse)
def tessellate(req: TessellateRequest) -> MeshResponse:
    try:
        entities = geometry_from_json(req.geometry)
        mesh = solidish_to_mesh(entities)
        return MeshResponse(**mesh)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


def _first_solid(entities: List[Any]) -> Any:
    for ent in entities:
        if issolid(ent):
            return ent
    raise ValueError("No solid entity found in geometry JSON")


@router.post("/boolean", response_model=BooleanResponse)
def boolean_op(req: BooleanRequest) -> BooleanResponse:
    engine_name = req.engine or os.environ.get("YAPCAD_BOOLEAN_ENGINE") or "native"
    engine = get_engine(engine_name)
    if engine is None:
        raise HTTPException(status_code=400, detail=f"Unknown boolean engine '{engine_name}'")

    try:
        solid_a = _first_solid(geometry_from_json(req.a))
        solid_b = _first_solid(geometry_from_json(req.b))

        kwargs: Dict[str, Any] = {}
        if engine_name == "trimesh":
            backend = req.backend or os.environ.get("YAPCAD_TRIMESH_BACKEND")
            kwargs["backend"] = backend
            # stitch isn't supported in trimesh_engine signature (it accepts stitch but ignores); pass anyway
            kwargs["stitch"] = req.stitch
        elif engine_name == "native":
            kwargs["stitch"] = req.stitch

        result = engine.solid_boolean(solid_a, solid_b, req.operation, **kwargs)
        doc = geometry_to_json(
            [result],
            generator={
                "name": "yapcad-service",
                "boolean": {"engine": engine_name, "operation": req.operation},
            },
        )
        return BooleanResponse(geometry=doc)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
