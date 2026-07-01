"""Tests for src/yapcad/dsl/transforms/metadata.py — Step 4 of the v1.1 migration.

The MetadataTransform runs *before* interpretation, walking the Module AST
and validating / normalising every FunctionDef.meta_hint dict.

Test strategy:
- Build minimal real FunctionDef/Module AST objects (no DSL parse needed).
- Verify: valid hints pass with no errors; invalid hints produce errors;
  coercions happen; cross-field checks fire; unknown fields get warnings.
"""

from __future__ import annotations

import pytest
from typing import Any, Dict, List, Optional

# ---------------------------------------------------------------------------
# Helpers to build minimal real AST objects
# ---------------------------------------------------------------------------

from yapcad.dsl.tokens import SourceSpan, SourceLocation
from yapcad.dsl.ast import FunctionDef, Block, Module

_LOC = SourceLocation(1, 0, 0)
_SPAN = SourceSpan(_LOC, _LOC)
_BLOCK = Block(statements=[], span=_SPAN)


def _make_fn(name: str, meta_hint: Optional[Dict[str, Any]] = None) -> FunctionDef:
    """Build a minimal FunctionDef with the given meta_hint."""
    return FunctionDef(
        span=_SPAN,
        name=name,
        parameters=[],
        return_type=None,
        body=_BLOCK,
        meta_hint=meta_hint,
    )


def _run(fns: List[FunctionDef]):
    """Apply MetadataTransform to a module built from *fns*; return (module, diags)."""
    module = Module(span=_SPAN, name="test_module", functions=fns)
    from yapcad.dsl.transforms.metadata import MetadataTransform, get_transform_diagnostics
    t = MetadataTransform()
    result = t.transform(module)
    return result, get_transform_diagnostics(result)


# ---------------------------------------------------------------------------
# Import the public API under test
# ---------------------------------------------------------------------------

from yapcad.dsl.transforms.metadata import (
    MetadataTransform,
    MetaDiagnostic,
    get_transform_diagnostics,
)


# ---------------------------------------------------------------------------
# 1. No meta_hint → no diagnostics
# ---------------------------------------------------------------------------

def test_no_meta_hint_is_clean():
    _, diags = _run([_make_fn("PLAIN_CMD", meta_hint=None)])
    assert diags == []


def test_empty_meta_hint_is_clean():
    _, diags = _run([_make_fn("EMPTY_CMD", meta_hint={})])
    assert diags == []


# ---------------------------------------------------------------------------
# 2. Valid hints — no errors
# ---------------------------------------------------------------------------

def test_valid_operation_subtract():
    fn = _make_fn("CUTTER", meta_hint={
        "operation.kind": "subtract",
        "operation.policy": "strict",
        "operation.feature_kind": "pocket",
        "operation.through": True,
        "operation.consume": False,
        "operation.target_filter": ["forward_bulkhead", "intertank"],
        "operation.stage": "clearance_cut",
        "operation.priority": 1.0,
    })
    _, diags = _run([fn])
    errors = [d for d in diags if d.level == "error"]
    assert errors == [], errors


def test_valid_operation_intersect():
    fn = _make_fn("INTER", meta_hint={"operation.kind": "intersect", "operation.target_filter": ["x"]})
    _, diags = _run([fn])
    errors = [d for d in diags if d.level == "error"]
    assert errors == [], errors


def test_valid_operation_union():
    fn = _make_fn("UNION", meta_hint={"operation.kind": "union", "operation.target_filter": ["x"]})
    _, diags = _run([fn])
    errors = [d for d in diags if d.level == "error"]
    assert errors == [], errors


def test_valid_assembly_scalars():
    fn = _make_fn("JOINT_PART", meta_hint={
        "assembly.joint_kind": "revolute",
        "assembly.no_cut": True,
    })
    _, diags = _run([fn])
    errors = [d for d in diags if d.level == "error"]
    assert errors == [], errors


def test_all_valid_joint_kinds():
    for jk in ("revolute", "prismatic", "spherical", "fixed", "planar", "cylindrical", "universal"):
        fn = _make_fn(f"PART_{jk}", meta_hint={"assembly.joint_kind": jk})
        _, diags = _run([fn])
        errors = [d for d in diags if d.level == "error"]
        assert errors == [], f"joint_kind={jk!r} should be valid, got: {errors}"


def test_valid_assembly_bolt_patterns():
    fn = _make_fn("BULKHEAD", meta_hint={
        "assembly.bolt_patterns": [
            {"id": "primary", "ring": "axial", "R_mm": 149.4, "z_mm": 317.3, "count": 8},
        ],
    })
    _, diags = _run([fn])
    errors = [d for d in diags if d.level == "error"]
    assert errors == [], errors


def test_valid_assembly_datums():
    fn = _make_fn("TANK", meta_hint={
        "assembly.datums": [
            {"id": "top_face", "kind": "axis", "ring": "upper", "R_mm": 150.0, "z_mm": 500.0},
        ],
    })
    _, diags = _run([fn])
    errors = [d for d in diags if d.level == "error"]
    assert errors == [], errors


def test_valid_assembly_surfaces_keepouts():
    fn = _make_fn("PART", meta_hint={
        "assembly.surfaces": [{"id": "mate_top", "kind": "mating", "mate_to": "nosecone.base"}],
        "assembly.keepouts": [{"id": "zone", "kind": "volume", "reason": "approach"}],
    })
    _, diags = _run([fn])
    errors = [d for d in diags if d.level == "error"]
    assert errors == [], errors


def test_valid_root_fields():
    fn = _make_fn("PART", meta_hint={
        "layer": "kinematics",
        "tags": ["structural", "aluminium"],
    })
    _, diags = _run([fn])
    errors = [d for d in diags if d.level == "error"]
    assert errors == [], errors


# ---------------------------------------------------------------------------
# 3. Enum validation → errors on bad values
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("kind", ["add", "merge", "remove", "boolean", ""])
def test_bad_operation_kind(kind):
    fn = _make_fn("BAD_CUTTER", meta_hint={"operation.kind": kind})
    _, diags = _run([fn])
    errors = [d for d in diags if d.level == "error"]
    assert any("operation.kind" in d.field for d in errors), \
        f"Expected operation.kind error for kind={kind!r}, got: {diags}"


@pytest.mark.parametrize("policy", ["hard", "soft", "none", "error"])
def test_bad_operation_policy(policy):
    fn = _make_fn("BAD_POLICY", meta_hint={
        "operation.kind": "subtract",
        "operation.policy": policy,
    })
    _, diags = _run([fn])
    errors = [d for d in diags if d.level == "error"]
    assert any("operation.policy" in d.field for d in errors), \
        f"Expected operation.policy error for policy={policy!r}, got: {diags}"


@pytest.mark.parametrize("jk", ["hinge", "slider", "free", "weld"])
def test_bad_assembly_joint_kind(jk):
    fn = _make_fn("BAD_JOINT", meta_hint={"assembly.joint_kind": jk})
    _, diags = _run([fn])
    errors = [d for d in diags if d.level == "error"]
    assert any("assembly.joint_kind" in d.field for d in errors), \
        f"Expected assembly.joint_kind error for {jk!r}, got: {diags}"


@pytest.mark.parametrize("fk", ["groove", "tapped", "knurl", "emboss"])
def test_bad_feature_kind(fk):
    fn = _make_fn("BAD_FEATURE", meta_hint={
        "operation.kind": "subtract",
        "operation.feature_kind": fk,
    })
    _, diags = _run([fn])
    errors = [d for d in diags if d.level == "error"]
    assert any("operation.feature_kind" in d.field for d in errors), \
        f"Expected operation.feature_kind error for {fk!r}, got: {diags}"


# ---------------------------------------------------------------------------
# 4. Boolean coercions
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("raw,expected", [
    (True, True), (False, False),
    ("true", True), ("True", True), ("1", True), ("yes", True),
    ("false", False), ("False", False), ("0", False), ("no", False),
])
def test_bool_coercion_no_cut(raw, expected):
    fn = _make_fn("PART", meta_hint={"assembly.no_cut": raw})
    _, diags = _run([fn])
    errors = [d for d in diags if d.level == "error"]
    assert errors == [], errors
    assert fn.meta_hint["assembly.no_cut"] == expected


@pytest.mark.parametrize("raw,expected", [
    ("true", True), ("false", False), (True, True), (False, False),
])
def test_bool_coercion_operation_through(raw, expected):
    fn = _make_fn("CUTTER", meta_hint={
        "operation.kind": "subtract",
        "operation.through": raw,
    })
    _, diags = _run([fn])
    errors = [d for d in diags if d.level == "error"]
    assert errors == [], errors
    assert fn.meta_hint["operation.through"] == expected


def test_bad_bool_coercion():
    fn = _make_fn("BAD", meta_hint={
        "operation.kind": "subtract",
        "operation.consume": "maybe",  # not a valid bool
    })
    _, diags = _run([fn])
    errors = [d for d in diags if d.level == "error"]
    assert any("operation.consume" in d.field for d in errors)


# ---------------------------------------------------------------------------
# 5. Float coercion for operation.priority
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("raw,expected", [
    (1, 1.0), (2.5, 2.5), ("3.0", 3.0), ("10", 10.0),
])
def test_priority_coercion(raw, expected):
    fn = _make_fn("CUTTER", meta_hint={
        "operation.kind": "subtract",
        "operation.priority": raw,
    })
    _, diags = _run([fn])
    errors = [d for d in diags if d.level == "error"]
    assert errors == [], errors
    assert fn.meta_hint["operation.priority"] == pytest.approx(expected)


def test_bad_priority_coercion():
    fn = _make_fn("BAD", meta_hint={
        "operation.kind": "subtract",
        "operation.priority": "high",  # not numeric
    })
    _, diags = _run([fn])
    errors = [d for d in diags if d.level == "error"]
    assert any("operation.priority" in d.field for d in errors)


# ---------------------------------------------------------------------------
# 6. target_filter coercion (string → list)
# ---------------------------------------------------------------------------

def test_target_filter_string_promoted_to_list():
    fn = _make_fn("CUTTER", meta_hint={
        "operation.kind": "subtract",
        "operation.target_filter": "forward_bulkhead",
    })
    _, diags = _run([fn])
    errors = [d for d in diags if d.level == "error"]
    assert errors == [], errors
    assert fn.meta_hint["operation.target_filter"] == ["forward_bulkhead"]


def test_target_filter_csv_string():
    fn = _make_fn("CUTTER", meta_hint={
        "operation.kind": "subtract",
        "operation.target_filter": "forward_bulkhead, intertank, nosecone",
    })
    _, diags = _run([fn])
    assert fn.meta_hint["operation.target_filter"] == [
        "forward_bulkhead", "intertank", "nosecone"
    ]


def test_target_filter_list_passthrough():
    ids = ["forward_bulkhead", "intertank"]
    fn = _make_fn("CUTTER", meta_hint={
        "operation.kind": "subtract",
        "operation.target_filter": list(ids),
    })
    _, diags = _run([fn])
    assert fn.meta_hint["operation.target_filter"] == ids


# ---------------------------------------------------------------------------
# 7. Cross-field checks
# ---------------------------------------------------------------------------

def test_operation_without_kind_is_error():
    fn = _make_fn("NO_KIND", meta_hint={
        "operation.policy": "strict",  # operation.* present but no kind
    })
    _, diags = _run([fn])
    errors = [d for d in diags if d.level == "error"]
    assert any("operation.kind" in d.field for d in errors), \
        f"Expected operation.kind required-field error, got: {diags}"


def test_no_cut_plus_subtract_is_warning():
    fn = _make_fn("CONTRADICTORY", meta_hint={
        "assembly.no_cut": True,
        "operation.kind": "subtract",
    })
    _, diags = _run([fn])
    errors = [d for d in diags if d.level == "error"]
    warnings = [d for d in diags if d.level == "warning"]
    assert errors == [], errors
    assert any("no_cut" in d.field or "contradict" in d.message.lower() for d in warnings), \
        f"Expected no_cut contradiction warning, got: {diags}"


def test_operation_without_target_filter_produces_info():
    fn = _make_fn("BROAD_CUTTER", meta_hint={
        "operation.kind": "subtract",
        # No target_filter → info diagnostic
    })
    _, diags = _run([fn])
    errors = [d for d in diags if d.level == "error"]
    infos = [d for d in diags if d.level == "info"]
    assert errors == [], errors
    assert any("target_filter" in d.field for d in infos), \
        f"Expected target_filter info, got: {diags}"


def test_target_filter_present_suppresses_info():
    fn = _make_fn("PRECISE_CUTTER", meta_hint={
        "operation.kind": "subtract",
        "operation.target_filter": ["forward_bulkhead"],
    })
    _, diags = _run([fn])
    infos = [d for d in diags if d.level == "info" and "target_filter" in d.field]
    assert infos == [], f"Unexpected target_filter info: {infos}"


# ---------------------------------------------------------------------------
# 8. Unknown fields → warnings, not errors
# ---------------------------------------------------------------------------

def test_unknown_namespace_is_warning():
    fn = _make_fn("PART", meta_hint={"thermal.expansion": 12.5})
    _, diags = _run([fn])
    errors = [d for d in diags if d.level == "error"]
    warnings = [d for d in diags if d.level == "warning"]
    assert errors == [], errors
    assert any("thermal" in d.message or "thermal" in d.field for d in warnings)


def test_unknown_root_key_is_warning():
    fn = _make_fn("PART", meta_hint={"colour": "red"})
    _, diags = _run([fn])
    errors = [d for d in diags if d.level == "error"]
    warnings = [d for d in diags if d.level == "warning"]
    assert errors == [], errors
    assert any("colour" in d.message or "colour" in d.field for d in warnings)


def test_unknown_assembly_field_is_warning():
    fn = _make_fn("PART", meta_hint={"assembly.future_field": "v2_value"})
    _, diags = _run([fn])
    errors = [d for d in diags if d.level == "error"]
    warnings = [d for d in diags if d.level == "warning"]
    assert errors == [], errors
    assert any("future_field" in d.message or "future_field" in d.field for d in warnings)


def test_unknown_operation_field_is_warning():
    fn = _make_fn("CUTTER", meta_hint={
        "operation.kind": "subtract",
        "operation.future_option": True,
    })
    _, diags = _run([fn])
    errors = [d for d in diags if d.level == "error"]
    warnings = [d for d in diags if d.level == "warning"]
    assert errors == [], errors
    assert any("future_option" in d.message or "future_option" in d.field for d in warnings)


# ---------------------------------------------------------------------------
# 9. Assembly list-of-dict structural checks
# ---------------------------------------------------------------------------

def test_assembly_list_field_not_a_list_is_error():
    fn = _make_fn("PART", meta_hint={
        "assembly.bolt_patterns": {"id": "primary"},  # dict, not list of dicts
    })
    _, diags = _run([fn])
    errors = [d for d in diags if d.level == "error"]
    assert any("bolt_patterns" in d.field for d in errors), \
        f"Expected bolt_patterns structural error, got: {diags}"


def test_assembly_list_field_entry_not_dict_is_error():
    fn = _make_fn("PART", meta_hint={
        "assembly.datums": ["not_a_dict", 42],
    })
    _, diags = _run([fn])
    errors = [d for d in diags if d.level == "error"]
    assert any("datums" in d.field for d in errors)


# ---------------------------------------------------------------------------
# 10. Multiple functions — diagnostics keyed per-function
# ---------------------------------------------------------------------------

def test_multiple_fns_independent_diagnostics():
    good = _make_fn("GOOD", meta_hint={
        "operation.kind": "subtract",
        "operation.target_filter": ["x"],
    })
    bad = _make_fn("BAD", meta_hint={
        "operation.policy": "strict",  # operation.* but no kind → error
    })
    _, diags = _run([good, bad])

    error_cmds = {d.command for d in diags if d.level == "error"}
    assert "BAD" in error_cmds, f"Expected BAD to have errors, got: {error_cmds}"
    assert "GOOD" not in error_cmds, f"Expected GOOD to be clean, got: {error_cmds}"


# ---------------------------------------------------------------------------
# 11. Idempotency — running twice gives same result
# ---------------------------------------------------------------------------

def test_idempotent():
    fn = _make_fn("CUTTER", meta_hint={
        "operation.kind": "subtract",
        "operation.through": "true",   # string → will be coerced to bool
        "operation.target_filter": "forward_bulkhead",  # string → list
    })
    module = Module(span=_SPAN, name="test_module", functions=[fn])
    t = MetadataTransform()
    t.transform(module)
    hint_after_first = dict(fn.meta_hint)
    t.transform(module)
    hint_after_second = dict(fn.meta_hint)
    assert hint_after_first == hint_after_second


# ---------------------------------------------------------------------------
# 12. get_transform_diagnostics on un-transformed module → empty
# ---------------------------------------------------------------------------

def test_untransformed_module_returns_empty_diags():
    module = Module(span=_SPAN, name="test_module", functions=[])
    diags = get_transform_diagnostics(module)
    assert diags == []


# ---------------------------------------------------------------------------
# 13. Verify diagnostics carry command name and level correctly
# ---------------------------------------------------------------------------

def test_diagnostic_fields_are_populated():
    fn = _make_fn("MY_CUTTER", meta_hint={"operation.policy": "bad_val", "operation.kind": "subtract"})
    _, diags = _run([fn])
    for d in diags:
        assert isinstance(d, MetaDiagnostic)
        assert d.command == "MY_CUTTER"
        assert d.level in ("error", "warning", "info")
        assert isinstance(d.field, str) and d.field
        assert isinstance(d.message, str) and d.message
