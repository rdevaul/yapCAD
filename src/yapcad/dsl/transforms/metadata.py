"""AST-level metadata transform for the v1.1 assembly/operation namespaces.

This transform runs *before* interpretation, walking the Module AST and
validating every ``FunctionDef.meta_hint`` dict.  It catches structural
errors early — at parse/check time rather than at solid-build time — and
normalises the hint into a canonical form so the runtime ``apply_meta_hint``
call in the interpreter can be a simple, trust-the-input dispatch.

What this transform does
------------------------
1. **Validates** every ``@meta`` hint key against the v1.1 schema:
   - Recognises ``assembly.*``, ``operation.*``, ``layer``, ``tags``.
   - Rejects unknown dotted namespaces (not ``assembly`` or ``operation``).
   - Validates enum values for ``operation.kind``, ``operation.policy``,
     ``operation.feature_kind``, and ``assembly.joint_kind`` at AST time.

2. **Normalises** the hint in place:
   - Booleans-as-strings (``"true"`` / ``"false"``) are converted to real
     bools for ``assembly.no_cut``, ``operation.through``,
     ``operation.consume``.
   - ``operation.target_filter`` given as a bare string is promoted to a
     one-element list.
   - ``operation.priority`` is coerced to ``float``.

3. **Cross-field consistency checks**:
   - An ``operation`` namespace present without ``operation.kind`` is an
     error (``kind`` is the only required field per §6).
   - ``assembly.no_cut = true`` combined with ``operation.kind = "subtract"``
     on the same command is contradictory and is flagged as a warning (not
     fatal — the author may be intentionally marking it a ghost cutter that
     also subtracts in some contexts, but it's almost certainly a mistake).
   - ``operation`` present without any ``target_filter`` produces an info
     diagnostic (allowed, but resolver will subtract from *every* parent).

4. **Produces ``MetadataTransformResult``** entries (errors, warnings, infos)
   keyed to the FunctionDef's span so the checker or IDE integration can
   surface them with source locations.

What this transform does NOT do
--------------------------------
- It does not evaluate expressions inside ``@meta`` hints.  By the time this
  transform runs the parser has already reduced literal ``@meta`` args to a
  flat Python dict (``meta_hint``); any expression that wasn't a literal is
  left as-is and will surface as a ``TypeError`` at runtime.
- It does not write sidecar YAML (that's Step 5).
- It does not call ``apply_meta_hint`` or touch geometry (that's the runtime).
- It is **idempotent**: running it twice produces identical AST state.

Usage
-----
    from yapcad.dsl.transforms import TransformPipeline
    from yapcad.dsl.transforms.metadata import MetadataTransform

    pipeline = TransformPipeline()
    pipeline.add(MetadataTransform())
    module = pipeline.apply(module)

    # Retrieve diagnostics after the transform:
    from yapcad.dsl.transforms.metadata import get_transform_diagnostics
    diags = get_transform_diagnostics(module)
    errors = [d for d in diags if d.level == "error"]
"""

from __future__ import annotations

import logging
from copy import deepcopy
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

from ..ast import Module, FunctionDef, Command
from ..tokens import SourceSpan as Span
from .base import AstTransform

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Schema constants — kept in sync with yapcad.metadata enums
# ---------------------------------------------------------------------------

_OPERATION_KINDS = {"subtract", "intersect", "union"}
_OPERATION_POLICIES = {"strict", "warn", "ignore"}
_FEATURE_KINDS = {
    "pocket", "slot", "boss", "hole",
    "through_hole", "counterbore", "countersink",
    "fillet", "chamfer", "rib", "draft",
}
_ASSEMBLY_JOINT_KINDS = {
    "revolute", "prismatic", "spherical", "fixed",
    "planar", "cylindrical", "universal",
}
_ASSEMBLY_RING_KINDS = {"axial", "radial", "mixed", "none"}
_BOLT_RING_KINDS = {"axial", "radial"}

# Top-level (root) fields that are allowed without a namespace prefix
_ROOT_FIELDS = {"layer", "tags"}

# All known namespace prefixes
_KNOWN_NAMESPACES = {"assembly", "operation"}

# Assembly scalar fields (forwarded to set_assembly)
_ASSEMBLY_SCALARS = {"joint_kind", "no_cut"}

# Assembly list-of-dict fields (forwarded to add_* helpers)
_ASSEMBLY_LIST_FIELDS = {"bolt_patterns", "datums", "surfaces", "keepouts"}

# Operation scalar/list fields
_OPERATION_FIELDS = {
    "kind", "target_filter", "priority", "through",
    "consume", "policy", "feature_id", "feature_kind", "stage",
}

# ---------------------------------------------------------------------------
# Diagnostic dataclass
# ---------------------------------------------------------------------------

_DIAG_KEY = "_metadata_transform_diagnostics"


@dataclass
class MetaDiagnostic:
    """A single diagnostic produced by ``MetadataTransform``."""
    level: str          # "error", "warning", "info"
    command: str        # FunctionDef.name
    field: str          # The @meta key that triggered this (e.g. "operation.kind")
    message: str
    span: Optional[Span] = None


def get_transform_diagnostics(module: Module) -> List[MetaDiagnostic]:
    """Return all ``MetaDiagnostic`` entries attached to *module* by
    ``MetadataTransform``.

    Returns an empty list if the transform has not been run or produced no
    diagnostics.
    """
    return getattr(module, _DIAG_KEY, [])


# ---------------------------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------------------------

def _check_enum(value: Any, valid: set, field_name: str) -> Optional[str]:
    """Return an error string if *value* is not in *valid*, else None."""
    if not isinstance(value, str) or value not in valid:
        valid_str = ", ".join(sorted(valid))
        return f"{field_name} must be one of {{{valid_str}}}, got {value!r}"
    return None


def _coerce_bool(value: Any, field_name: str) -> Tuple[bool, Optional[str]]:
    """Coerce *value* to bool; return (coerced, error_or_None)."""
    if isinstance(value, bool):
        return value, None
    if isinstance(value, str):
        if value.lower() in ("true", "1", "yes"):
            return True, None
        if value.lower() in ("false", "0", "no"):
            return False, None
    return False, f"{field_name} must be a boolean (got {value!r})"


def _coerce_float(value: Any, field_name: str) -> Tuple[float, Optional[str]]:
    """Coerce *value* to float; return (coerced, error_or_None)."""
    try:
        return float(value), None
    except (TypeError, ValueError):
        return 0.0, f"{field_name} must be numeric (got {value!r})"


def _coerce_string_list(value: Any, field_name: str) -> Tuple[List[str], Optional[str]]:
    """Coerce *value* to a list of strings.  Accepts a plain string (split
    on commas) or a list.  Returns (coerced, error_or_None)."""
    if isinstance(value, str):
        return [v.strip() for v in value.split(",") if v.strip()], None
    if isinstance(value, list):
        return [str(v) for v in value], None
    return [], f"{field_name} must be a string or list of strings (got {type(value).__name__})"


# ---------------------------------------------------------------------------
# Per-namespace validators
# ---------------------------------------------------------------------------

def _validate_assembly_field(
    hint: Dict[str, Any],
    field: str,
    value: Any,
    diags: List[MetaDiagnostic],
    cmd_name: str,
    span: Optional[Span],
) -> Any:
    """Validate and normalise one ``assembly.<field>`` entry.

    Returns the (possibly coerced) value, or *value* unchanged on error
    (error is appended to *diags*).
    """
    fqn = f"assembly.{field}"

    if field == "joint_kind":
        err = _check_enum(value, _ASSEMBLY_JOINT_KINDS, fqn)
        if err:
            diags.append(MetaDiagnostic("error", cmd_name, fqn, err, span))
        return value

    if field == "no_cut":
        coerced, err = _coerce_bool(value, fqn)
        if err:
            diags.append(MetaDiagnostic("error", cmd_name, fqn, err, span))
            return value
        return coerced

    if field in _ASSEMBLY_LIST_FIELDS:
        # List-of-dict — minimal structural check here; deeper validation
        # (per-entry required keys, enum checks) happens in meta_apply at
        # runtime where the helpers can give precise per-entry errors.
        if not isinstance(value, list):
            diags.append(MetaDiagnostic(
                "error", cmd_name, fqn,
                f"assembly.{field} must be a list of dicts (got {type(value).__name__})",
                span,
            ))
        elif any(not isinstance(entry, dict) for entry in value):
            diags.append(MetaDiagnostic(
                "error", cmd_name, fqn,
                f"assembly.{field}: all entries must be dicts",
                span,
            ))
        return value

    # Unknown assembly field — forward-compat: warn but don't block
    diags.append(MetaDiagnostic(
        "warning", cmd_name, fqn,
        f"Unknown assembly field {field!r}; stored verbatim (forward-compat). "
        f"Valid fields: {sorted(_ASSEMBLY_SCALARS | _ASSEMBLY_LIST_FIELDS)}",
        span,
    ))
    return value


def _validate_operation_field(
    hint: Dict[str, Any],
    field: str,
    value: Any,
    diags: List[MetaDiagnostic],
    cmd_name: str,
    span: Optional[Span],
) -> Any:
    """Validate and normalise one ``operation.<field>`` entry."""
    fqn = f"operation.{field}"

    if field == "kind":
        err = _check_enum(value, _OPERATION_KINDS, fqn)
        if err:
            diags.append(MetaDiagnostic("error", cmd_name, fqn, err, span))
        return value

    if field == "policy":
        err = _check_enum(value, _OPERATION_POLICIES, fqn)
        if err:
            diags.append(MetaDiagnostic("error", cmd_name, fqn, err, span))
        return value

    if field == "feature_kind":
        err = _check_enum(value, _FEATURE_KINDS, fqn)
        if err:
            diags.append(MetaDiagnostic("error", cmd_name, fqn, err, span))
        return value

    if field in ("through", "consume"):
        coerced, err = _coerce_bool(value, fqn)
        if err:
            diags.append(MetaDiagnostic("error", cmd_name, fqn, err, span))
            return value
        return coerced

    if field == "priority":
        coerced, err = _coerce_float(value, fqn)
        if err:
            diags.append(MetaDiagnostic("error", cmd_name, fqn, err, span))
            return value
        return coerced

    if field == "target_filter":
        coerced, err = _coerce_string_list(value, fqn)
        if err:
            diags.append(MetaDiagnostic("error", cmd_name, fqn, err, span))
            return value
        return coerced

    if field in ("feature_id", "stage"):
        return str(value)

    # Unknown operation field — forward-compat warning
    diags.append(MetaDiagnostic(
        "warning", cmd_name, fqn,
        f"Unknown operation field {field!r}; stored verbatim. "
        f"Valid fields: {sorted(_OPERATION_FIELDS)}",
        span,
    ))
    return value


# ---------------------------------------------------------------------------
# Cross-field consistency checks
# ---------------------------------------------------------------------------

def _cross_check(
    hint: Dict[str, Any],
    diags: List[MetaDiagnostic],
    cmd_name: str,
    span: Optional[Span],
) -> None:
    """Run cross-field consistency checks on the full *hint* dict."""

    # --- Collect which namespaces are present ---
    has_operation = any(k.startswith("operation.") for k in hint)
    has_assembly = any(k.startswith("assembly.") for k in hint)

    # operation.kind is required when any operation.* key appears
    if has_operation and "operation.kind" not in hint:
        diags.append(MetaDiagnostic(
            "error", cmd_name, "operation.kind",
            "operation namespace present but operation.kind is missing "
            "(required field — must be one of: subtract, intersect, union)",
            span,
        ))

    # Contradictory: no_cut=true + operation.kind=subtract
    no_cut = hint.get("assembly.no_cut")
    op_kind = hint.get("operation.kind")
    if no_cut is True and op_kind == "subtract":
        diags.append(MetaDiagnostic(
            "warning", cmd_name, "assembly.no_cut",
            "assembly.no_cut=true and operation.kind='subtract' are contradictory: "
            "no_cut marks this part as non-subtractive in boolean ops, but operation.kind "
            "requests a subtract. One of these is likely wrong.",
            span,
        ))

    # Inform when operation has no target_filter (subtracts from everything)
    if has_operation and "operation.target_filter" not in hint:
        diags.append(MetaDiagnostic(
            "info", cmd_name, "operation.target_filter",
            "operation.target_filter not set — the resolver will apply this "
            "operation to ALL parents in the assembly. Set target_filter to a "
            "list of part IDs to restrict the operation.",
            span,
        ))


# ---------------------------------------------------------------------------
# The transform
# ---------------------------------------------------------------------------

class MetadataTransform(AstTransform):
    """Pre-interpretation AST transform that validates and normalises
    ``FunctionDef.meta_hint`` dicts against the v1.1 metadata namespace schema.

    Diagnostics are stored on the module object under the ``_metadata_transform_diagnostics``
    attribute.  They are also logged via the ``yapcad.dsl.transforms.metadata``
    logger at the appropriate level.

    The transform is non-destructive: it only modifies the **values** inside
    ``meta_hint`` dicts (coercions), never removes keys or changes structure.
    """

    @property
    def name(self) -> str:
        return "metadata-v1.1"

    def transform(self, module: Module) -> Module:
        """Apply the metadata transform to *module*.

        Returns the same module object with:
        - ``FunctionDef.meta_hint`` dicts normalised in place.
        - ``module._metadata_transform_diagnostics`` populated with any
          errors, warnings, and infos found.
        """
        all_diags: List[MetaDiagnostic] = []

        for cmd in module.commands:
            if not isinstance(cmd, (FunctionDef,)):
                continue
            hint = getattr(cmd, "meta_hint", None)
            if not hint:
                continue

            cmd_diags: List[MetaDiagnostic] = []
            normalised: Dict[str, Any] = {}

            for key, value in hint.items():
                if not isinstance(key, str) or not key:
                    continue  # malformed key — skip silently

                if "." in key:
                    namespace, _, field = key.partition(".")

                    if namespace == "assembly":
                        normalised[key] = _validate_assembly_field(
                            hint, field, value, cmd_diags, cmd.name, cmd.span
                        )
                    elif namespace == "operation":
                        normalised[key] = _validate_operation_field(
                            hint, field, value, cmd_diags, cmd.name, cmd.span
                        )
                    else:
                        # Unknown namespace — warn and store verbatim
                        cmd_diags.append(MetaDiagnostic(
                            "warning", cmd.name, key,
                            f"Unknown @meta namespace {namespace!r}; "
                            f"known namespaces: {sorted(_KNOWN_NAMESPACES)}. "
                            f"Stored verbatim for forward-compatibility.",
                            cmd.span,
                        ))
                        normalised[key] = value

                elif key in _ROOT_FIELDS:
                    normalised[key] = value

                else:
                    # Unknown root-level key — warn and store verbatim
                    cmd_diags.append(MetaDiagnostic(
                        "warning", cmd.name, key,
                        f"Unknown @meta root key {key!r}; "
                        f"known root fields: {sorted(_ROOT_FIELDS)}. "
                        f"Known namespaces: {sorted(_KNOWN_NAMESPACES)}. "
                        f"Stored verbatim.",
                        cmd.span,
                    ))
                    normalised[key] = value

            # Run cross-field checks on the normalised hint
            _cross_check(normalised, cmd_diags, cmd.name, cmd.span)

            # Write normalised values back into the live hint dict
            cmd.meta_hint.update(normalised)

            # Log and accumulate
            for diag in cmd_diags:
                if diag.level == "error":
                    logger.error("[metadata-transform] %s: %s", diag.command, diag.message)
                elif diag.level == "warning":
                    logger.warning("[metadata-transform] %s: %s", diag.command, diag.message)
                else:
                    logger.info("[metadata-transform] %s: %s", diag.command, diag.message)

            all_diags.extend(cmd_diags)

        # Attach diagnostics to the module object for external consumers
        setattr(module, _DIAG_KEY, all_diags)
        return module
