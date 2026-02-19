from __future__ import annotations

import os
from typing import Any, Dict

from fastapi import APIRouter

import yapcad
from yapcad.boolean import ENGINE_REGISTRY
from yapcad.boolean import trimesh_engine
from yapcad.boolean import occ_engine
from yapcad.dsl.introspection import get_api_reference

from ..models.schemas import BuiltinsResponse, HealthResponse, InfoResponse

router = APIRouter(tags=["system"])


@router.get("/health", response_model=HealthResponse)
def health() -> HealthResponse:
    return HealthResponse()


@router.get("/info", response_model=InfoResponse)
def info() -> InfoResponse:
    engines: Dict[str, Any] = {"registered": sorted(list(ENGINE_REGISTRY.keys()))}

    # trimesh engine availability
    try:
        available = [b for b in trimesh_engine.engines_available() if b]
        engines["trimesh"] = {
            "installed": trimesh_engine.trimesh is not None,
            "available_backends": sorted([str(b) for b in available]),
        }
    except Exception as exc:
        engines["trimesh"] = {"error": str(exc)}

    # occ engine availability
    try:
        engines["occ"] = {"available": bool(occ_engine.is_available())}
    except Exception as exc:
        engines["occ"] = {"error": str(exc)}

    env = {
        "YAPCAD_BOOLEAN_ENGINE": os.environ.get("YAPCAD_BOOLEAN_ENGINE"),
        "YAPCAD_TRIMESH_BACKEND": os.environ.get("YAPCAD_TRIMESH_BACKEND"),
        "YAPCAD_DSL_RECURSION_LIMIT": os.environ.get("YAPCAD_DSL_RECURSION_LIMIT"),
    }

    return InfoResponse(
        yapcad_version=getattr(yapcad, "__version__", "unknown"),
        boolean_engines=engines,
        env=env,
    )


@router.get("/builtins", response_model=BuiltinsResponse)
def builtins() -> BuiltinsResponse:
    return BuiltinsResponse(api=get_api_reference())
