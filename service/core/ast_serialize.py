from __future__ import annotations

from dataclasses import is_dataclass, fields
from enum import Enum
from typing import Any, Dict


def _span_to_dict(span: Any) -> Dict[str, Any]:
    # yapcad.dsl.tokens.SourceSpan
    try:
        start = span.start
        end = span.end
        return {
            "start": {
                "line": start.line,
                "column": start.column,
                "offset": start.offset,
                "filename": getattr(start, "filename", None),
            },
            "end": {
                "line": end.line,
                "column": end.column,
                "offset": end.offset,
                "filename": getattr(end, "filename", None),
            },
        }
    except Exception:
        return {"repr": str(span)}


def ast_to_dict(obj: Any) -> Any:
    """Best-effort JSON-serializable representation of yapCAD DSL AST.

    The AST is a graph of dataclasses with Enums and SourceSpan objects.
    """
    if obj is None:
        return None

    if isinstance(obj, (str, int, float, bool)):
        return obj

    if isinstance(obj, Enum):
        return obj.name

    if isinstance(obj, (list, tuple)):
        return [ast_to_dict(v) for v in obj]

    if isinstance(obj, dict):
        return {str(k): ast_to_dict(v) for k, v in obj.items()}

    if is_dataclass(obj):
        out: Dict[str, Any] = {"_type": obj.__class__.__name__}
        for f in fields(obj):
            val = getattr(obj, f.name)
            if f.name == "span":
                out["span"] = _span_to_dict(val)
            else:
                out[f.name] = ast_to_dict(val)
        return out

    # fallback
    if hasattr(obj, "__dict__"):
        return {"_type": obj.__class__.__name__, **{k: ast_to_dict(v) for k, v in obj.__dict__.items()}}

    return str(obj)
