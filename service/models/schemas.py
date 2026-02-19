from __future__ import annotations

from typing import Any, Dict, List, Literal, Optional

from pydantic import BaseModel, Field


class ErrorResponse(BaseModel):
    detail: str
    error_type: Optional[str] = None


class DslEvalRequest(BaseModel):
    source: str = Field(..., description="DSL source code")
    command: str = Field(..., description="Command/function to run")
    parameters: Dict[str, Any] = Field(default_factory=dict)
    format: Literal["json"] = Field("json", description="Output format")


class DslEvalResponse(BaseModel):
    success: bool
    geometry: Optional[Dict[str, Any]] = None
    metadata: Dict[str, Any] = Field(default_factory=dict)
    provenance: Optional[Dict[str, Any]] = None
    require_failures: List[Dict[str, Any]] = Field(default_factory=list)
    error_message: Optional[str] = None


class DslParseRequest(BaseModel):
    source: str


class DslParseResponse(BaseModel):
    ast: Dict[str, Any]


class TessellateRequest(BaseModel):
    geometry: Dict[str, Any] = Field(..., description="Geometry JSON document")


class MeshResponse(BaseModel):
    vertices: List[float]
    normals: List[float]
    indices: List[int]


class BooleanRequest(BaseModel):
    operation: Literal["union", "difference", "intersection"]
    a: Dict[str, Any] = Field(..., description="Geometry JSON for first solid")
    b: Dict[str, Any] = Field(..., description="Geometry JSON for second solid")
    engine: Optional[str] = Field(
        None,
        description="Boolean engine override (native|trimesh|occ). Defaults to YAPCAD_BOOLEAN_ENGINE env var or 'native'",
    )
    backend: Optional[str] = Field(
        None,
        description="Backend override for trimesh engine (e.g. 'blender'). Defaults to YAPCAD_TRIMESH_BACKEND",
    )
    stitch: bool = Field(False, description="Whether to attempt to stitch open edges (native engine only)")


class BooleanResponse(BaseModel):
    geometry: Dict[str, Any]


class BuiltinsResponse(BaseModel):
    api: Dict[str, Any]


class HealthResponse(BaseModel):
    status: Literal["ok"] = "ok"


class InfoResponse(BaseModel):
    yapcad_version: str
    boolean_engines: Dict[str, Any]
    env: Dict[str, Optional[str]]
