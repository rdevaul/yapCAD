from __future__ import annotations

import time
from typing import Any, Dict

from fastapi import APIRouter, HTTPException, WebSocket, WebSocketDisconnect

from yapcad.io.geometry_json import geometry_to_json
from yapcad.geom import isgeomlist
from yapcad.geom3d import issolid, issurface

from ..core.dsl_runner import DslTimeoutError, eval_dsl, parse_dsl
from ..models.schemas import (
    DslAssemblyEvalRequest,
    DslAssemblyEvalResponse,
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
        # If the emitted value is a string, treat it as an assembly snapshot
        # rather than trying geometry serialization (which would always fail).
        # result.emit_result is an EmitResult wrapper; .data is the raw Python value.
        if isinstance(result.emit_result.data, str):
            resp.assembly_snapshot = result.emit_result.data
        else:
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
    """Extract command names, parameters, types, defaults, and return types from DSL source."""
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
        cmd_info: Dict[str, Any] = {"name": fn["name"], "params": params}
        # Surface the return type so clients can route accordingly
        # (e.g. '-> string' commands should use /dsl/assembly_eval, not /dsl/eval)
        rt = fn.get("return_type")
        if isinstance(rt, dict):
            cmd_info["return_type"] = rt.get("name", "unknown")
        elif isinstance(rt, str):
            cmd_info["return_type"] = rt
        else:
            cmd_info["return_type"] = None
        meta_hint = fn.get("meta_hint")
        if meta_hint:
            cmd_info["meta_hint"] = meta_hint
        commands.append(cmd_info)
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


@router.post("/assembly_eval", response_model=DslAssemblyEvalResponse)
async def dsl_assembly_eval(req: DslAssemblyEvalRequest) -> DslAssemblyEvalResponse:
    """Evaluate a DSL command that emits an assembly snapshot (``-> string``).

    Runs the command with full provenance tracking and a generous timeout
    (assembly evaluation involves multiple solid builds), then returns the
    raw JSON string produced by ``emit_assembly()`` without any attempt at
    geometry serialisation.

    Use this endpoint for commands declared ``-> string``.  For geometry
    commands use ``/dsl/eval``; for scalar/list commands use ``/dsl/ui_eval``.
    """
    try:
        result = await eval_dsl(
            source=req.source,
            command=req.command,
            parameters=req.parameters,
            timeout_s=120.0,  # assembly eval can be slow — multiple solids
        )
    except DslTimeoutError as exc:
        raise HTTPException(status_code=408, detail=str(exc)) from exc
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc

    resp = DslAssemblyEvalResponse(
        success=bool(result.success),
        metadata=getattr(result, "metadata", {}) or {},
        provenance=result.provenance.to_dict() if getattr(result, "provenance", None) else None,
        require_failures=[
            {"message": rf.message, "expression_text": rf.expression_text}
            for rf in (result.require_failures or [])
        ],
        error_message=result.error_message,
    )

    if not result.success:
        return resp

    if result.emit_result is None:
        resp.success = False
        resp.error_message = "Command did not emit a value"
        return resp

    # emit_result is an EmitResult wrapper; .data is the raw Python value
    raw = result.emit_result.data

    if not isinstance(raw, str):
        resp.success = False
        resp.error_message = (
            f"Expected command to emit a string (assembly snapshot), "
            f"got {type(raw).__name__}. "
            f"Use /dsl/eval for geometry commands."
        )
        return resp

    resp.assembly_snapshot = raw
    return resp


@router.post("/parse", response_model=DslParseResponse)
def dsl_parse(req: DslParseRequest) -> DslParseResponse:
    try:
        ast_dict, _module = parse_dsl(source=req.source)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    return DslParseResponse(ast=ast_dict)


# ── WebSocket eval endpoint ───────────────────────────────────────────────────
# Mirrors the REST /dsl/eval endpoint over a persistent WebSocket connection.
# Protocol:
#   client → {"type": "eval", "request_id": str, "source": str,
#              "command": str, "parameters": {}, "format": "json"}
#   client → {"type": "cancel", "request_id": str}  (abort in-flight eval)
#   server → {"type": "result", "request_id": str, "success": bool,
#              "geometry": ...|null, "assembly_snapshot": ...|null,
#              "error_message": ...|null, "elapsed_ms": float}

_WS_ROUTER = APIRouter(tags=["dsl"])


@_WS_ROUTER.websocket("/ws/dsl/eval")
async def ws_dsl_eval(websocket: WebSocket) -> None:
    """WebSocket endpoint for DSL evaluation.

    Accepts the same eval logic as POST /dsl/eval, with request-id
    multiplexing and in-flight cancellation via 'cancel' messages.
    """
    import asyncio
    import json

    await websocket.accept()
    current_task: asyncio.Task | None = None
    current_request_id: str | None = None

    try:
        while True:
            raw = await websocket.receive_text()
            try:
                msg = json.loads(raw)
            except json.JSONDecodeError:
                await websocket.send_json({"type": "error", "error_message": "Invalid JSON"})
                continue

            msg_type = msg.get("type")

            if msg_type == "cancel":
                rid = msg.get("request_id")
                if current_task and not current_task.done() and rid == current_request_id:
                    current_task.cancel()
                    current_request_id = None
                continue

            if msg_type != "eval":
                await websocket.send_json({"type": "error", "error_message": f"Unknown message type: {msg_type!r}"})
                continue

            request_id: str = msg.get("request_id", "")
            source: str = msg.get("source", "")
            command: str = msg.get("command", "")
            parameters: dict = msg.get("parameters") or {}

            # Cancel any in-flight eval before starting a new one
            if current_task and not current_task.done():
                current_task.cancel()

            async def _run_eval(rid: str, src: str, cmd: str, params: dict) -> None:
                t0 = time.perf_counter()
                response: dict = {"type": "result", "request_id": rid}
                try:
                    result = await eval_dsl(
                        source=src,
                        command=cmd,
                        parameters=params,
                        timeout_s=300.0,
                    )
                    response["success"] = bool(result.success)
                    response["error_message"] = result.error_message
                    response["assembly_snapshot"] = None
                    response["geometry"] = None

                    if result.emit_result is not None:
                        if isinstance(result.emit_result.data, str):
                            response["assembly_snapshot"] = result.emit_result.data
                        else:
                            try:
                                from yapcad.io.geometry_json import geometry_to_json
                                entities = _normalize_emitted_geometry(result.geometry)
                                response["geometry"] = geometry_to_json(
                                    entities,
                                    units=None,
                                    generator={"name": "yapcad-service", "dslCommand": cmd},
                                )
                            except Exception as geo_exc:
                                response["success"] = False
                                response["error_message"] = f"Geometry serialization failed: {geo_exc}"

                except DslTimeoutError as exc:
                    response["success"] = False
                    response["error_message"] = str(exc)
                except asyncio.CancelledError:
                    return  # silently drop cancelled eval
                except Exception as exc:
                    response["success"] = False
                    response["error_message"] = str(exc)

                response["elapsed_ms"] = (time.perf_counter() - t0) * 1000
                try:
                    await websocket.send_json(response)
                except Exception:
                    pass  # connection may have closed

            current_request_id = request_id
            current_task = asyncio.ensure_future(_run_eval(request_id, source, command, parameters))

    except WebSocketDisconnect:
        if current_task and not current_task.done():
            current_task.cancel()


# Register the WS router on the module so main.py can include it
ws_router = _WS_ROUTER
