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


_SCALAR_TYPES = (int, float, bool, str)


def _is_scalar_or_scalar_list(value: Any) -> bool:
    """Check if value is a scalar or a list of scalars."""
    if isinstance(value, _SCALAR_TYPES):
        return True
    if isinstance(value, list) and all(isinstance(v, _SCALAR_TYPES) for v in value):
        return True
    return False


def _scalar_type_name(value: Any) -> str:
    """Get a type name string for a scalar or list-of-scalar value."""
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
        elem = value[0]
        return f"list[{_scalar_type_name(elem)}]"
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


def _normalize_emitted_geometry(geom: Any):
    if geom is None:
        return []

    # Single solid/surface
    if issolid(geom) or issurface(geom):
        return [geom]

    # Curve/spline — wrap in a geomlist so _serialize_sketch sees the full curve structure
    if _is_curve_value(geom):
        return [[geom]]  # [[curve]] = geomlist containing the curve

    # geomlist (list of points/vectors/segments)
    if isgeomlist(geom):
        return [geom]

    # List of entities
    if isinstance(geom, list) and geom and all(
        (issolid(g) or issurface(g) or isgeomlist(g) or _is_curve_value(g)) for g in geom
    ):
        result = []
        for g in geom:
            result.extend(_normalize_emitted_geometry(g))
        return result

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
        raw = result.geometry
        if _is_scalar_or_scalar_list(raw):
            # Scalar/list emit — populate scalar_result instead of geometry
            values = raw if isinstance(raw, list) else [raw]
            resp.scalar_result = {
                "values": values,
                "type": _scalar_type_name(raw),
            }
        else:
            try:
                entities = _normalize_emitted_geometry(raw)
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
            param_info = {"name": p["name"], "type": ptype, "default": default}
            ui_hint = p.get("ui_hint")
            if ui_hint:
                param_info["ui_hint"] = ui_hint
            params.append(param_info)
        commands.append({"name": fn["name"], "params": params})
    return {"commands": commands}


@router.post("/ui_eval", response_model=DslUiEvalResponse)
async def dsl_ui_eval(req: DslUiEvalRequest) -> DslUiEvalResponse:
    """Evaluate a DSL command expecting scalar/list emit (for UI parameter population)."""
    try:
        result = await eval_dsl(
            source=req.source,
            command=req.command,
            parameters=req.parameters,
            timeout_s=10.0,
        )
    except DslTimeoutError as exc:
        return DslUiEvalResponse(success=False, error_message=str(exc))
    except Exception as exc:
        return DslUiEvalResponse(success=False, error_message=str(exc))

    if not result.success:
        return DslUiEvalResponse(success=False, error_message=result.error_message)

    if result.emit_result is None:
        return DslUiEvalResponse(success=False, error_message="Command did not emit a value")

    raw = result.geometry
    if not _is_scalar_or_scalar_list(raw):
        return DslUiEvalResponse(
            success=False,
            error_message=f"ui_eval expects scalar/list result, got {type(raw).__name__}",
        )

    values = raw if isinstance(raw, list) else [raw]
    return DslUiEvalResponse(success=True, values=values, type=_scalar_type_name(raw))


@router.post("/callgraph")
def dsl_callgraph(req: DslParseRequest):
    """Extract call graph from DSL source: nodes, edges, and cycle detection."""
    try:
        _ast_dict, module = parse_dsl(source=req.source)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc

    from yapcad.dsl.ast import (
        FunctionDef, FunctionCall, Identifier, MethodCall,
        Block, Statement, Expression,
        VarDecl, AssignmentStatement, EmitStatement, ForStatement,
        IfStatement, ExpressionStatement, ReturnStatement,
        BinaryOp, UnaryOp, MemberAccess, IndexAccess,
        ListLiteral, ListComprehension, ComprehensionClause,
        RangeExpr, ConditionalExpr, IfExpr, MatchExpr,
        LambdaExpr, DictLiteral, ElifBranch, AssertStatement,
    )

    # Collect all command names defined in this module
    command_names = {fn.name for fn in module.functions}

    def _collect_calls_expr(expr: Expression, calls: set):
        """Walk an expression tree collecting function call references."""
        if expr is None:
            return
        if isinstance(expr, FunctionCall):
            if isinstance(expr.callee, Identifier) and expr.callee.name in command_names:
                calls.add(expr.callee.name)
            _collect_calls_expr(expr.callee, calls)
            for arg in expr.arguments:
                _collect_calls_expr(arg, calls)
            for v in expr.named_arguments.values():
                _collect_calls_expr(v, calls)
        elif isinstance(expr, MethodCall):
            _collect_calls_expr(expr.object, calls)
            for arg in expr.arguments:
                _collect_calls_expr(arg, calls)
            for v in expr.named_arguments.values():
                _collect_calls_expr(v, calls)
        elif isinstance(expr, BinaryOp):
            _collect_calls_expr(expr.left, calls)
            _collect_calls_expr(expr.right, calls)
        elif isinstance(expr, UnaryOp):
            _collect_calls_expr(expr.operand, calls)
        elif isinstance(expr, MemberAccess):
            _collect_calls_expr(expr.object, calls)
        elif isinstance(expr, IndexAccess):
            _collect_calls_expr(expr.object, calls)
            _collect_calls_expr(expr.index, calls)
        elif isinstance(expr, ListLiteral):
            for e in expr.elements:
                _collect_calls_expr(e, calls)
        elif isinstance(expr, ListComprehension):
            _collect_calls_expr(expr.element_expr, calls)
            for clause in expr.clauses:
                _collect_calls_expr(clause.iterable, calls)
                for cond in clause.conditions:
                    _collect_calls_expr(cond, calls)
        elif isinstance(expr, RangeExpr):
            _collect_calls_expr(expr.start, calls)
            _collect_calls_expr(expr.end, calls)
            if expr.step:
                _collect_calls_expr(expr.step, calls)
        elif isinstance(expr, ConditionalExpr):
            _collect_calls_expr(expr.condition, calls)
            _collect_calls_expr(expr.true_branch, calls)
            _collect_calls_expr(expr.false_branch, calls)
        elif isinstance(expr, IfExpr):
            _collect_calls_expr(expr.condition, calls)
            _collect_calls_block(expr.then_branch, calls)
            for eb in expr.elif_branches:
                _collect_calls_expr(eb.condition, calls)
                _collect_calls_block(eb.body, calls)
            if expr.else_branch:
                _collect_calls_block(expr.else_branch, calls)
        elif isinstance(expr, MatchExpr):
            _collect_calls_expr(expr.subject, calls)
            for arm in expr.arms:
                _collect_calls_expr(arm.body, calls)
        elif isinstance(expr, LambdaExpr):
            _collect_calls_expr(expr.body, calls)
        elif isinstance(expr, DictLiteral):
            for v in expr.entries.values():
                _collect_calls_expr(v, calls)

    def _collect_calls_stmt(stmt: Statement, calls: set):
        """Walk a statement collecting function call references."""
        if isinstance(stmt, VarDecl):
            if stmt.initializer:
                _collect_calls_expr(stmt.initializer, calls)
        elif isinstance(stmt, AssignmentStatement):
            _collect_calls_expr(stmt.value, calls)
        elif isinstance(stmt, EmitStatement):
            _collect_calls_expr(stmt.value, calls)
            if isinstance(stmt.metadata, dict):
                for v in stmt.metadata.values():
                    _collect_calls_expr(v, calls)
        elif isinstance(stmt, ForStatement):
            _collect_calls_expr(stmt.iterable, calls)
            _collect_calls_block(stmt.body, calls)
        elif isinstance(stmt, IfStatement):
            _collect_calls_expr(stmt.condition, calls)
            _collect_calls_block(stmt.then_branch, calls)
            for eb in stmt.elif_branches:
                _collect_calls_expr(eb.condition, calls)
                _collect_calls_block(eb.body, calls)
            if stmt.else_branch:
                _collect_calls_block(stmt.else_branch, calls)
        elif isinstance(stmt, ExpressionStatement):
            _collect_calls_expr(stmt.expression, calls)
        elif isinstance(stmt, ReturnStatement):
            if stmt.value:
                _collect_calls_expr(stmt.value, calls)
        elif isinstance(stmt, AssertStatement):
            _collect_calls_expr(stmt.condition, calls)
        elif isinstance(stmt, Block):
            _collect_calls_block(stmt, calls)

    def _collect_calls_block(block: Block, calls: set):
        """Walk a block collecting function call references."""
        for stmt in block.statements:
            _collect_calls_stmt(stmt, calls)

    # Build adjacency: command_name -> set of commands it calls
    adjacency: dict[str, set] = {}
    for fn in module.functions:
        calls: set = set()
        _collect_calls_block(fn.body, calls)
        adjacency[fn.name] = calls

    # Build edges
    edges = []
    for caller, callees in adjacency.items():
        for callee in sorted(callees):
            edges.append({"from": caller, "to": callee})

    # Detect cycles using DFS with coloring
    WHITE, GRAY, BLACK = 0, 1, 2
    color = {name: WHITE for name in command_names}
    cycles = []
    path: list[str] = []

    def _dfs(node: str):
        color[node] = GRAY
        path.append(node)
        for neighbor in sorted(adjacency.get(node, set())):
            if neighbor not in color:
                continue
            if color[neighbor] == GRAY:
                # Back edge — extract cycle
                idx = path.index(neighbor)
                cycles.append(list(path[idx:]))
            elif color[neighbor] == WHITE:
                _dfs(neighbor)
        path.pop()
        color[node] = BLACK

    for node in sorted(command_names):
        if color[node] == WHITE:
            _dfs(node)

    return {
        "nodes": sorted(command_names),
        "edges": edges,
        "cycles": cycles,
    }


@router.post("/parse", response_model=DslParseResponse)
def dsl_parse(req: DslParseRequest) -> DslParseResponse:
    try:
        ast_dict, _module = parse_dsl(source=req.source)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    return DslParseResponse(ast=ast_dict)
