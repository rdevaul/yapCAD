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
    assembly_snapshot: Optional[str] = Field(
        None,
        description="JSON snapshot string emitted by assembly commands (-> string return type)",
    )
    metadata: Dict[str, Any] = Field(default_factory=dict)
    provenance: Optional[Dict[str, Any]] = None
    require_failures: List[Dict[str, Any]] = Field(default_factory=list)
    error_message: Optional[str] = None


class DslParseRequest(BaseModel):
    source: str


class DslParseResponse(BaseModel):
    ast: Dict[str, Any]


class DslUiEvalRequest(BaseModel):
    """Request a DSL command evaluation that returns a scalar or list of scalars.

    Used by the workbench to evaluate ``@ui``-driven parameter commands —
    e.g. a command that computes a derived value from widget state and returns
    a float or list of floats for display.
    """
    source: str = Field(..., description="DSL source code")
    command: str = Field(..., description="Command/function to run")
    parameters: Dict[str, Any] = Field(default_factory=dict)


class DslUiEvalResponse(BaseModel):
    """Response from a ``/dsl/ui_eval`` request."""
    success: bool
    values: Optional[List[Any]] = None
    type: Optional[str] = None  # "int", "float", "bool", "string"
    error_message: Optional[str] = None


class DslAssemblyEvalRequest(BaseModel):
    """Request evaluation of a ``-> string`` DSL command that emits an assembly snapshot.

    The service executes the command, expects ``emit_result`` to be a JSON
    string produced by ``emit_assembly()``, and returns it in the response
    without attempting geometry serialisation.
    """
    source: str = Field(..., description="DSL source code")
    command: str = Field(..., description="Command to run (must return -> string)")
    parameters: Dict[str, Any] = Field(default_factory=dict)


class DslAssemblyEvalResponse(BaseModel):
    """Response from a ``/dsl/assembly_eval`` request."""
    success: bool
    assembly_snapshot: Optional[str] = Field(
        None,
        description="Raw JSON string produced by emit_assembly()",
    )
    metadata: Dict[str, Any] = Field(default_factory=dict)
    provenance: Optional[Dict[str, Any]] = None
    require_failures: List[Dict[str, Any]] = Field(default_factory=list)
    error_message: Optional[str] = None


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
