"""Apply @meta decorator hints from the DSL to yapCAD geometry metadata.

This module bridges the gap between the DSL's ``@meta(...)`` command decorator
(which produces a flat ``meta_hint`` dict on a ``FunctionDef`` node) and the
v1.1 metadata namespace helpers in ``yapcad.metadata``.

Typical usage — called by the DSL evaluator or service layer after a command
emits a solid::

    from yapcad.dsl.meta_apply import apply_meta_hint

    if fn.meta_hint:
        apply_meta_hint(solid, fn.meta_hint)

The function is **idempotent and additive**: calling it multiple times with
the same hints merges cleanly because the underlying namespace helpers use
``setdefault`` / field-level writes, not wholesale replacement.

Namespace routing
-----------------
Dotted keys in the hint dict are split on the first dot to determine the
namespace, then the remainder becomes the field name:

    ``assembly.joint_kind``  → assembly namespace, field ``joint_kind``
    ``operation.kind``       → operation namespace, field ``kind``
    ``layer``                → root metadata field ``layer``
    ``tags``                 → root metadata ``tags`` list (value appended)

Unknown namespaces or unrecognised field names within a known namespace are
stored verbatim in the namespace dict rather than raising an error.  This
keeps the DSL forward-compatible: a field defined in a future v1.2 spec won't
break existing parsers.

Limitations
-----------
The ``@meta`` decorator carries only **scalar** values (str, int, float, bool).
Structured sub-objects like ``bolt_patterns``, ``datums``, ``surfaces``, and
``keepouts`` require the full ``add_bolt_pattern()`` / ``add_datum()`` /
``add_surface()`` / ``add_keepout()`` helpers and cannot be expressed as flat
``@meta`` kwargs.  They remain available for direct Python use.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from yapcad.metadata import (
    get_solid_metadata,
    _ensure_root,
    set_assembly,
    set_operation,
    get_assembly_metadata,
    get_operation_metadata,
)

# ---------------------------------------------------------------------------
# Assembly namespace field routing
# ---------------------------------------------------------------------------

# Fields accepted as kwargs by set_assembly() — passed through directly.
_ASSEMBLY_SCALAR_FIELDS = {"joint_kind", "no_cut"}

# Fields that are list-of-dicts — can only be set with the dedicated helpers.
_ASSEMBLY_LIST_FIELDS = {"bolt_patterns", "datums", "surfaces", "keepouts"}


def _apply_assembly(meta: Dict[str, Any], field: str, value: Any) -> None:
    """Route a single ``assembly.<field>`` hint into the metadata dict."""
    if field in _ASSEMBLY_SCALAR_FIELDS:
        # Convert bool-like strings for convenience ("true" → True)
        if field == "no_cut" and isinstance(value, str):
            value = value.lower() in ("true", "1", "yes")
        set_assembly(meta, **{field: value})
    elif field in _ASSEMBLY_LIST_FIELDS:
        raise ValueError(
            f"assembly.{field} requires structured data; "
            "use add_bolt_pattern() / add_datum() / add_surface() / add_keepout() directly."
        )
    else:
        # Unknown field — store verbatim for forward-compatibility
        section = get_assembly_metadata(meta, create=True)
        section[field] = value


# ---------------------------------------------------------------------------
# Operation namespace field routing
# ---------------------------------------------------------------------------

# Fields accepted as kwargs by set_operation() — passed through directly.
# (target_filter is a list, handled separately)
_OPERATION_SCALAR_FIELDS = {
    "kind", "priority", "through", "consume", "policy",
    "feature_id", "feature_kind", "stage",
}


def _apply_operation(meta: Dict[str, Any], field: str, value: Any) -> None:
    """Route a single ``operation.<field>`` hint into the metadata dict.

    ``set_operation()`` requires ``kind`` to be present on every call.  To
    support per-field application from a flat dict (where ``kind`` may have
    been set in a previous iteration), we use ``set_operation`` only when
    ``kind`` is the field being set; for all other scalar fields we write
    directly to the namespace section after enum validation.
    """
    from yapcad.metadata import _validate_enum, _OPERATION_KINDS, _OPERATION_POLICIES, _FEATURE_KINDS

    if field == "target_filter":
        # Accept a list value or a comma-separated string
        if isinstance(value, str):
            value = [v.strip() for v in value.split(",") if v.strip()]
        elif not isinstance(value, list):
            value = [str(value)]
        section = get_operation_metadata(meta, create=True)
        section["target_filter"] = value
    elif field == "kind":
        # kind is required by set_operation — call it directly
        _validate_enum(value, _OPERATION_KINDS, field="operation.kind")
        section = get_operation_metadata(meta, create=True)
        section["kind"] = value
    elif field == "policy":
        _validate_enum(value, _OPERATION_POLICIES, field="operation.policy")
        section = get_operation_metadata(meta, create=True)
        section["policy"] = value
    elif field == "feature_kind":
        _validate_enum(value, _FEATURE_KINDS, field="operation.feature_kind")
        section = get_operation_metadata(meta, create=True)
        section["feature_kind"] = value
    elif field in ("through", "consume"):
        if isinstance(value, str):
            value = value.lower() in ("true", "1", "yes")
        section = get_operation_metadata(meta, create=True)
        section[field] = bool(value)
    elif field == "priority":
        section = get_operation_metadata(meta, create=True)
        section["priority"] = float(value)
    elif field in ("feature_id", "stage"):
        section = get_operation_metadata(meta, create=True)
        section[field] = str(value)
    else:
        # Unknown field — store verbatim for forward-compat
        section = get_operation_metadata(meta, create=True)
        section[field] = value


# ---------------------------------------------------------------------------
# Root-level field routing
# ---------------------------------------------------------------------------

def _apply_root(meta: Dict[str, Any], field: str, value: Any) -> None:
    """Route a plain (non-namespaced) hint key into the metadata root."""
    root = _ensure_root(meta)
    if field == "layer":
        root["layer"] = str(value)
    elif field == "tags":
        # Append a single tag string or extend with a list
        tags: List[str] = root.setdefault("tags", [])
        if isinstance(value, list):
            tags.extend(str(t) for t in value)
        else:
            tags.append(str(value))
    else:
        # Free-form root field — store verbatim
        root[field] = value


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def apply_meta_hint(solid: list, meta_hint: Dict[str, Any]) -> None:
    """Apply a ``FunctionDef.meta_hint`` dict to a yapCAD solid's metadata.

    Args:
        solid:     A yapCAD solid (list representation from ``geom3d``).
        meta_hint: The flat dict produced by one or more ``@meta(...)``
                   decorators.  Keys are plain strings or dotted namespace
                   paths (``"assembly.joint_kind"``, ``"operation.kind"``,
                   ``"layer"``, etc.).

    Raises:
        TypeError:  If ``solid`` is not a solid or ``meta_hint`` is not a dict.
        ValueError: If a structured list field (bolt_patterns, datums, …) is
                    attempted via the flat interface.

    Example::

        @meta(assembly.joint_kind="revolute", layer="kinematics")
        @meta(operation.kind="subtract", operation.feature_kind="pocket")
        command MAKE_POCKET(depth: float = 10.0) -> solid:
            ...

        # After evaluation:
        apply_meta_hint(result_solid, fn.meta_hint)
        # → assembly.joint_kind = "revolute"
        # → operation.kind = "subtract", operation.feature_kind = "pocket"
        # → root layer = "kinematics"
    """
    if not isinstance(meta_hint, dict):
        raise TypeError(f"meta_hint must be a dict, got {type(meta_hint).__name__}")

    # Get or create the metadata dict for this solid
    meta = get_solid_metadata(solid, create=True)

    for key, value in meta_hint.items():
        if not isinstance(key, str) or not key:
            continue  # skip malformed keys silently

        if "." in key:
            namespace, _, field = key.partition(".")
            if namespace == "assembly":
                _apply_assembly(meta, field, value)
            elif namespace == "operation":
                _apply_operation(meta, field, value)
            else:
                # Unknown namespace — store verbatim under that namespace dict
                root = _ensure_root(meta)
                root.setdefault(namespace, {})[field] = value
        else:
            _apply_root(meta, key, value)


def apply_meta_hint_to_raw(meta: Dict[str, Any], meta_hint: Dict[str, Any]) -> None:
    """Apply a ``meta_hint`` dict directly to an already-extracted metadata dict.

    Convenience variant for callers that have already retrieved the metadata
    dict (e.g. via ``get_solid_metadata``) rather than the raw solid.

    Args:
        meta:      Metadata dict (mutable, modified in place).
        meta_hint: Flat hint dict from ``@meta`` decorators.
    """
    if not isinstance(meta_hint, dict):
        raise TypeError(f"meta_hint must be a dict, got {type(meta_hint).__name__}")

    for key, value in meta_hint.items():
        if not isinstance(key, str) or not key:
            continue

        if "." in key:
            namespace, _, field = key.partition(".")
            if namespace == "assembly":
                _apply_assembly(meta, field, value)
            elif namespace == "operation":
                _apply_operation(meta, field, value)
            else:
                root = _ensure_root(meta)
                root.setdefault(namespace, {})[field] = value
        else:
            _apply_root(meta, key, value)


__all__ = [
    "apply_meta_hint",
    "apply_meta_hint_to_raw",
]
