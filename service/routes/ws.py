"""WebSocket DSL evaluation endpoint with progress + cancellation."""

from __future__ import annotations

import asyncio
import time
from typing import Any, Dict, Optional

from fastapi import APIRouter, WebSocket, WebSocketDisconnect

from yapcad.io.geometry_json import geometry_to_json
from yapcad.geom import isgeomlist
from yapcad.geom3d import issolid, issurface

from ..core.dsl_runner import eval_dsl, DslTimeoutError

router = APIRouter(tags=["ws"])

_SCALAR_TYPES = (int, float, bool, str)


def _is_scalar_or_scalar_list(value: Any) -> bool:
    if isinstance(value, _SCALAR_TYPES):
        return True
    if isinstance(value, list) and all(isinstance(v, _SCALAR_TYPES) for v in value):
        return True
    return False


def _scalar_type_name(value: Any) -> str:
    if isinstance(value, bool):
        return "bool"
    if isinstance(value, int):
        return "int"
    if isinstance(value, float):
        return "float"
    if isinstance(value, str):
        return "str"
    if isinstance(value, list):
        if not value:
            return "list"
        return f"list[{_scalar_type_name(value[0])}]"
    return "unknown"


def _is_curve_value(geom: Any) -> bool:
    """Check if value is a yapCAD curve structure (catmullrom, bezier, nurbs, etc.)"""
    return (
        isinstance(geom, list)
        and len(geom) >= 2
        and isinstance(geom[0], str)
        and geom[0] in ("catmullrom", "bezier", "nurbs", "line_segment", "arc",
                        "circle", "ellipse", "parabola", "hyperbola", "path2d",
                        "path3d", "profile2d", "region2d", "loop3d")
    )


def _normalize_geometry(geom: Any):
    if geom is None:
        return []
    if issolid(geom) or issurface(geom):
        return [geom]
    # Curve/spline — wrap in a geomlist so _serialize_sketch sees the full curve structure
    if _is_curve_value(geom):
        return [[geom]]  # [[curve]] = geomlist containing the curve
    if isgeomlist(geom):
        return [geom]
    if isinstance(geom, list) and geom and all(
        (issolid(g) or issurface(g) or isgeomlist(g) or _is_curve_value(g)) for g in geom
    ):
        result = []
        for g in geom:
            result.extend(_normalize_geometry(g))
        return result
    raise ValueError(f"Unsupported geometry type: {type(geom).__name__}")


@router.websocket("/ws/dsl/eval")
async def ws_dsl_eval(ws: WebSocket):
    await ws.accept()

    # Track running evaluations for cancellation
    running: Dict[str, asyncio.Task] = {}

    try:
        while True:
            msg = await ws.receive_json()
            msg_type = msg.get("type")

            if msg_type == "eval":
                request_id = msg.get("request_id", "")
                command = msg.get("command", "")
                params = msg.get("parameters", msg.get("params", {}))
                source = msg.get("source", "")

                # Cancel any existing task with same request_id
                if request_id in running:
                    running[request_id].cancel()
                    del running[request_id]

                task = asyncio.create_task(
                    _run_eval(ws, request_id, source, command, params)
                )
                running[request_id] = task

                # Clean up when done
                def _on_done(t: asyncio.Task, rid=request_id):
                    running.pop(rid, None)
                task.add_done_callback(_on_done)

            elif msg_type == "cancel":
                request_id = msg.get("request_id", "")
                task = running.pop(request_id, None)
                if task and not task.done():
                    task.cancel()
                    try:
                        await ws.send_json({
                            "type": "result",
                            "request_id": request_id,
                            "success": False,
                            "error_message": "Cancelled by client",
                            "geometry": None,
                            "volume": None,
                        })
                    except Exception:
                        pass

    except WebSocketDisconnect:
        # Cancel all running tasks on disconnect
        for task in running.values():
            task.cancel()
    except Exception:
        for task in running.values():
            task.cancel()


async def _run_eval(
    ws: WebSocket,
    request_id: str,
    source: str,
    command: str,
    params: Dict[str, Any],
):
    """Run a DSL eval with progress reporting."""
    start_time = time.monotonic()

    # Start progress reporter
    progress_task = asyncio.create_task(
        _send_progress(ws, request_id, start_time)
    )

    try:
        result = await eval_dsl(
            source=source,
            command=command,
            parameters=params,
            timeout_s=30.0,
        )
    except asyncio.CancelledError:
        progress_task.cancel()
        raise
    except DslTimeoutError as exc:
        progress_task.cancel()
        await _send_result(ws, request_id, success=False, error_message=str(exc))
        return
    except Exception as exc:
        progress_task.cancel()
        await _send_result(ws, request_id, success=False, error_message=str(exc))
        return

    # Stop progress
    progress_task.cancel()

    if not result.success:
        await _send_result(
            ws, request_id,
            success=False,
            error_message=result.error_message,
        )
        return

    if result.emit_result is None:
        await _send_result(
            ws, request_id,
            success=False,
            error_message="Command did not emit any value",
        )
        return

    # Build response
    raw = result.geometry
    geometry_json = None
    scalar_result = None
    volume = None

    if _is_scalar_or_scalar_list(raw):
        values = raw if isinstance(raw, list) else [raw]
        scalar_result = {"values": values, "type": _scalar_type_name(raw)}
    else:
        try:
            entities = _normalize_geometry(raw)
            geometry_json = geometry_to_json(
                entities,
                units=None,
                generator={"name": "yapcad-service", "dslCommand": command},
            )
            # Extract volume if available
            if geometry_json and "entities" in geometry_json:
                for ent in geometry_json["entities"]:
                    v = ent.get("properties", {}).get("volume")
                    if v is not None:
                        volume = v
                        break
        except Exception as exc:
            await _send_result(
                ws, request_id,
                success=False,
                error_message=f"Geometry serialization failed: {exc}",
            )
            return

    await _send_result(
        ws, request_id,
        success=True,
        geometry=geometry_json,
        scalar_result=scalar_result,
        volume=volume,
    )


async def _send_progress(ws: WebSocket, request_id: str, start_time: float):
    """Send periodic progress updates."""
    try:
        while True:
            await asyncio.sleep(0.5)
            elapsed_ms = int((time.monotonic() - start_time) * 1000)
            await ws.send_json({
                "type": "progress",
                "request_id": request_id,
                "elapsed_ms": elapsed_ms,
            })
    except asyncio.CancelledError:
        pass
    except Exception:
        pass  # Connection may have closed


async def _send_result(
    ws: WebSocket,
    request_id: str,
    *,
    success: bool,
    geometry: Optional[Dict] = None,
    scalar_result: Optional[Dict] = None,
    volume: Optional[float] = None,
    error_message: Optional[str] = None,
):
    """Send an eval result."""
    msg: Dict[str, Any] = {
        "type": "result",
        "request_id": request_id,
        "success": success,
    }
    if geometry is not None:
        msg["geometry"] = geometry
    if scalar_result is not None:
        msg["scalar_result"] = scalar_result
    if volume is not None:
        msg["volume"] = volume
    if error_message is not None:
        msg["error_message"] = error_message
    try:
        await ws.send_json(msg)
    except Exception:
        pass
