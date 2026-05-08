from __future__ import annotations

from typing import Any, Dict

from fastapi import APIRouter, HTTPException

from yapcad.io.geometry_json import geometry_to_json
from yapcad.geom import isgeomlist
from yapcad.geom3d import issolid, issurface

from ..core.dsl_runner import DslTimeoutError, eval_dsl, parse_dsl
from ..models.schemas import (
    DslEvalRequest,
    DslEvalResponse,
    DslParseRequest,
    DslParseResponse,
    DslUiEvalRequest,
    DslUiEvalResponse,
)

router = APIRouter(prefix="/dsl", tags=["dsl"])


def _normalize_emitted_geometry(geom: Any):
    if geom is None:
        return []

    # Single entity
    if issolid(geom) or issurface(geom) or isgeomlist(geom):
        return [geom]

    # List of entities
    if isinstance(geom, list) and geom and all(
        (issolid(g) or issurface(g) or isgeomlist(g)) for g in geom
    ):
        return geom

    raise ValueError(f"Unsupported emitted geometry type: {type(geom).__name__}")


@router.post("/eval", response_model=DslEvalResponse)
async def dsl_eval(req: DslEvalRequest) -> DslEvalResponse:
    if req.format != "json":
        raise HTTPException(status_code=400, detail="Only format='json' is supported")

    try:
        result = await eval_dsl(
            source=req.source,
            command=req.command,
            parameters=req.parameters,
            timeout_s=30.0,
        )
    except DslTimeoutError as exc:
        raise HTTPException(status_code=408, detail=str(exc)) from exc
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc

    resp = DslEvalResponse(
        success=bool(result.success),
        metadata=getattr(result, "metadata", {}) or {},
        provenance=result.provenance.to_dict() if getattr(result, "provenance", None) else None,
        require_failures=[
            {"message": rf.message, "expression_text": rf.expression_text}
            for rf in (result.require_failures or [])
        ],
        error_message=result.error_message,
    )

    if result.emit_result is not None:
        try:
            entities = _normalize_emitted_geometry(result.geometry)
            resp.geometry = geometry_to_json(
                entities,
                units=None,
                generator={"name": "yapcad-service", "dslCommand": req.command},
            )
        except Exception as exc:
            # If serialization fails, still return execution info.
            resp.success = False
            resp.error_message = f"Geometry serialization failed: {exc}"

    return resp


@router.post("/commands")
def dsl_commands(req: DslParseRequest):
    """Extract command names, parameters, types, and defaults from DSL source."""
    try:
        ast_dict, _module = parse_dsl(source=req.source)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc

    commands = []
    for fn in ast_dict.get("functions", []):
        params = []
        for p in fn.get("parameters", []):
            ta = p.get("type_annotation")
            ptype = ta.get("name", "any") if isinstance(ta, dict) else "any"
            dv = p.get("default_value")
            default = dv.get("value") if isinstance(dv, dict) and "value" in dv else None
            param_info: Dict[str, Any] = {"name": p["name"], "type": ptype, "default": default}
            ui_hint = p.get("ui_hint")
            if ui_hint:
                param_info["ui_hint"] = ui_hint
            params.append(param_info)
        commands.append({"name": fn["name"], "params": params})
    return {"commands": commands}


@router.post("/ui_eval", response_model=DslUiEvalResponse)
async def dsl_ui_eval(req: DslUiEvalRequest) -> DslUiEvalResponse:
    """Evaluate a DSL command and return its scalar/list result.

    Intended for ``@ui``-driven commands that compute derived values (e.g. a
    command that returns a float from widget state).  Unlike ``/eval``, this
    endpoint does not attempt geometry serialisation — it expects the command
    to ``emit`` a scalar or list of scalars.
    """
    try:
        result = await eval_dsl(
            source=req.source,
            command=req.command,
            parameters=req.parameters,
            timeout_s=10.0,
        )
    except DslTimeoutError as exc:
        raise HTTPException(status_code=408, detail=str(exc)) from exc
    except Exception as exc:
        return DslUiEvalResponse(success=False, error_message=str(exc))

    if not result.success:
        return DslUiEvalResponse(success=False, error_message=result.error_message)

    raw = result.emit_result
    if raw is None:
        return DslUiEvalResponse(success=False, error_message="Command did not emit a value")

    # Coerce to list of scalars
    values = raw if isinstance(raw, list) else [raw]
    type_name = type(values[0]).__name__ if values else "unknown"
    return DslUiEvalResponse(success=True, values=values, type=type_name)


@router.post("/parse", response_model=DslParseResponse)
def dsl_parse(req: DslParseRequest) -> DslParseResponse:
    try:
        ast_dict, _module = parse_dsl(source=req.source)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    return DslParseResponse(ast=ast_dict)
