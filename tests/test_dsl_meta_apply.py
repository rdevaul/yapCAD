"""Tests for yapcad.dsl.meta_apply — @meta hint → v1.1 namespace bridge."""

import pytest
from yapcad.dsl.meta_apply import apply_meta_hint, apply_meta_hint_to_raw
from yapcad.metadata import (
    get_solid_metadata,
    get_assembly_metadata,
    get_operation_metadata,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_solid():
    """Return a minimal yapCAD solid stub sufficient for metadata operations.

    A valid solid is ['solid', surfaces, xforms, ids, optional_meta_dict].
    We use empty lists for surfaces/xforms/ids and an empty dict for metadata.
    """
    return ['solid', [], [], [], {}]


def get_assembly(solid):
    meta = get_solid_metadata(solid)
    return get_assembly_metadata(meta)


def get_operation(solid):
    meta = get_solid_metadata(solid)
    return get_operation_metadata(meta)


# ---------------------------------------------------------------------------
# Basic routing
# ---------------------------------------------------------------------------

class TestApplyMetaHintAssembly:
    def test_joint_kind_axial(self):
        s = make_solid()
        apply_meta_hint(s, {"assembly.joint_kind": "axial"})
        assert get_assembly(s)["joint_kind"] == "axial"

    def test_joint_kind_radial(self):
        s = make_solid()
        apply_meta_hint(s, {"assembly.joint_kind": "radial"})
        assert get_assembly(s)["joint_kind"] == "radial"

    def test_no_cut_true(self):
        s = make_solid()
        apply_meta_hint(s, {"assembly.no_cut": True})
        assert get_assembly(s)["no_cut"] is True

    def test_no_cut_string_coercion(self):
        """String 'true' is coerced to bool True."""
        s = make_solid()
        apply_meta_hint(s, {"assembly.no_cut": "true"})
        assert get_assembly(s)["no_cut"] is True

    def test_assembly_unknown_field_stored_verbatim(self):
        """Unknown assembly fields are stored without error (forward-compat)."""
        s = make_solid()
        apply_meta_hint(s, {"assembly.future_field": "some_value"})
        assert get_assembly(s)["future_field"] == "some_value"

    def test_assembly_list_field_non_list_raises(self):
        """List-shaped assembly fields reject non-list values."""
        s = make_solid()
        with pytest.raises(TypeError, match="bolt_patterns"):
            apply_meta_hint(s, {"assembly.bolt_patterns": "not a list"})

    def test_assembly_list_field_non_dict_entry_raises(self):
        """Entries inside the list must be dicts."""
        s = make_solid()
        with pytest.raises(TypeError, match=r"datums\[0\]"):
            apply_meta_hint(s, {"assembly.datums": ["not a dict"]})


class TestApplyMetaHintAssemblyListFields:
    """Coverage for the four ``assembly.<list>`` list-of-dict fields.

    These exercise the dispatch from ``@meta(assembly.<field>=[{...}, ...])``
    through ``_apply_assembly_list_field`` to the dedicated
    ``yapcad.metadata`` helpers (``add_bolt_pattern`` / ``add_datum`` /
    ``add_surface`` / ``add_keepout``).
    """

    def test_bolt_patterns_minimal(self):
        s = make_solid()
        apply_meta_hint(
            s,
            {
                "assembly.bolt_patterns": [
                    {
                        "id": "primary",
                        "ring": "axial",
                        "R_mm": 149.4,
                        "z_mm": 317.3,
                        "count": 8,
                    }
                ]
            },
        )
        bolt_patterns = get_assembly(s)["bolt_patterns"]
        assert len(bolt_patterns) == 1
        entry = bolt_patterns[0]
        assert entry["id"] == "primary"
        assert entry["ring"] == "axial"
        assert entry["R_mm"] == 149.4
        assert entry["z_mm"] == 317.3
        assert entry["count"] == 8

    def test_bolt_patterns_with_fastener_dict(self):
        s = make_solid()
        apply_meta_hint(
            s,
            {
                "assembly.bolt_patterns": [
                    {
                        "id": "primary",
                        "ring": "axial",
                        "R_mm": 149.4,
                        "z_mm": 317.3,
                        "count": 8,
                        "theta0_deg": 22.5,
                        "fastener": {
                            "designation": "1/4-20 SHCS",
                            "thread": "heat-set",
                            "head": "SHCS",
                        },
                        "tolerance_mm": 0.05,
                    }
                ]
            },
        )
        entry = get_assembly(s)["bolt_patterns"][0]
        assert entry["theta0_deg"] == 22.5
        assert entry["fastener"]["designation"] == "1/4-20 SHCS"
        assert entry["fastener"]["thread"] == "heat-set"
        assert entry["fastener"]["head"] == "SHCS"
        assert entry["tolerance_mm"] == 0.05

    def test_bolt_patterns_multiple_entries_append_in_order(self):
        s = make_solid()
        apply_meta_hint(
            s,
            {
                "assembly.bolt_patterns": [
                    {
                        "id": "top_ring",
                        "ring": "axial",
                        "R_mm": 149.4,
                        "z_mm": 317.3,
                        "count": 8,
                    },
                    {
                        "id": "side_ring",
                        "ring": "radial",
                        "R_mm": 145.0,
                        "z_mm": 12.7,
                        "count": 8,
                    },
                ]
            },
        )
        patterns = get_assembly(s)["bolt_patterns"]
        assert [p["id"] for p in patterns] == ["top_ring", "side_ring"]

    def test_bolt_patterns_invalid_ring_raises(self):
        """Enum validation from the helper surfaces correctly."""
        s = make_solid()
        with pytest.raises(ValueError, match="bolt_pattern.ring"):
            apply_meta_hint(
                s,
                {
                    "assembly.bolt_patterns": [
                        {
                            "id": "bad",
                            "ring": "sideways",  # not in _BOLT_RINGS
                            "R_mm": 100.0,
                            "z_mm": 0.0,
                            "count": 4,
                        }
                    ]
                },
            )

    def test_bolt_patterns_missing_required_key_raises(self):
        s = make_solid()
        with pytest.raises(TypeError, match=r"bolt_patterns\[0\] rejected"):
            apply_meta_hint(
                s,
                {"assembly.bolt_patterns": [{"ring": "upper"}]},  # no id, etc.
            )

    def test_datums_minimal(self):
        s = make_solid()
        apply_meta_hint(
            s,
            {
                "assembly.datums": [
                    {"id": "bore_axis", "kind": "axis"},
                    {
                        "id": "neck_bolt_circle",
                        "kind": "bolt_circle",
                        "ring": "axial",
                        "R_mm": 151.9,
                        "z_mm": 330.0,
                        "direction": [0.0, 0.0, 1.0],
                    },
                ]
            },
        )
        datums = get_assembly(s)["datums"]
        assert len(datums) == 2
        assert datums[0]["id"] == "bore_axis"
        assert datums[1]["direction"] == [0.0, 0.0, 1.0]
        assert datums[1]["R_mm"] == 151.9
        assert datums[1]["kind"] == "bolt_circle"

    def test_datums_invalid_direction_length_raises(self):
        s = make_solid()
        with pytest.raises(ValueError, match="datum.direction"):
            apply_meta_hint(
                s,
                {
                    "assembly.datums": [
                        {"id": "bad", "kind": "axis", "direction": [0.0, 1.0]}
                    ]
                },
            )

    def test_surfaces_minimal(self):
        s = make_solid()
        apply_meta_hint(
            s,
            {
                "assembly.surfaces": [
                    {
                        "id": "mate_top",
                        "kind": "mating",
                        "mate_to": "nosecone.base_bore",
                        "tolerance_mm": 0.05,
                    }
                ]
            },
        )
        surface = get_assembly(s)["surfaces"][0]
        assert surface["id"] == "mate_top"
        assert surface["kind"] == "mating"
        assert surface["mate_to"] == "nosecone.base_bore"
        assert surface["tolerance_mm"] == 0.05

    def test_keepouts_with_shape(self):
        s = make_solid()
        apply_meta_hint(
            s,
            {
                "assembly.keepouts": [
                    {
                        "id": "approach",
                        "kind": "volume",
                        "shape": {"type": "cylinder", "R_mm": 160.0, "h_mm": 50.0},
                        "reason": "hand clearance during install",
                    }
                ]
            },
        )
        keepout = get_assembly(s)["keepouts"][0]
        assert keepout["id"] == "approach"
        assert keepout["shape"]["type"] == "cylinder"
        assert keepout["reason"].startswith("hand clearance")

    def test_mixed_list_and_scalar_in_single_hint(self):
        """List fields coexist with scalar fields in the same hint dict."""
        s = make_solid()
        apply_meta_hint(
            s,
            {
                "assembly.joint_kind": "axial",
                "assembly.datums": [{"id": "bore", "kind": "axis"}],
                "assembly.bolt_patterns": [
                    {
                        "id": "primary",
                        "ring": "axial",
                        "R_mm": 149.4,
                        "z_mm": 317.3,
                        "count": 8,
                    }
                ],
                "layer": "kinematics",
            },
        )
        asm = get_assembly(s)
        assert asm["joint_kind"] == "axial"
        assert asm["datums"][0]["id"] == "bore"
        assert asm["bolt_patterns"][0]["id"] == "primary"
        # Root field still landed correctly
        from yapcad.metadata import get_solid_metadata, _ensure_root
        meta = get_solid_metadata(s)
        root = _ensure_root(meta)
        assert root["layer"] == "kinematics"


class TestApplyMetaHintOperation:
    def test_operation_kind(self):
        s = make_solid()
        apply_meta_hint(s, {"operation.kind": "subtract"})
        assert get_operation(s)["kind"] == "subtract"

    def test_operation_feature_kind(self):
        s = make_solid()
        # kind must come first so set_operation has a valid required field
        apply_meta_hint(s, {"operation.kind": "subtract"})
        apply_meta_hint(s, {"operation.feature_kind": "pocket"})
        op = get_operation(s)
        assert op["kind"] == "subtract"
        assert op["feature_kind"] == "pocket"

    def test_operation_through_bool(self):
        s = make_solid()
        apply_meta_hint(s, {"operation.kind": "subtract"})
        apply_meta_hint(s, {"operation.through": True})
        assert get_operation(s)["through"] is True

    def test_operation_through_string_coercion(self):
        s = make_solid()
        apply_meta_hint(s, {"operation.kind": "subtract"})
        apply_meta_hint(s, {"operation.through": "true"})
        assert get_operation(s)["through"] is True

    def test_operation_target_filter_list(self):
        s = make_solid()
        apply_meta_hint(s, {"operation.kind": "subtract"})
        apply_meta_hint(s, {"operation.target_filter": ["tank", "bulkhead"]})
        assert get_operation(s)["target_filter"] == ["tank", "bulkhead"]

    def test_operation_target_filter_csv_string(self):
        """Comma-separated string is split into a list."""
        s = make_solid()
        apply_meta_hint(s, {"operation.kind": "subtract"})
        apply_meta_hint(s, {"operation.target_filter": "tank, bulkhead"})
        assert get_operation(s)["target_filter"] == ["tank", "bulkhead"]

    def test_operation_priority(self):
        s = make_solid()
        apply_meta_hint(s, {"operation.kind": "subtract"})
        apply_meta_hint(s, {"operation.priority": 2.5})
        assert get_operation(s)["priority"] == 2.5

    def test_operation_unknown_field_verbatim(self):
        s = make_solid()
        apply_meta_hint(s, {"operation.kind": "subtract"})
        apply_meta_hint(s, {"operation.future_key": "xyz"})
        assert get_operation(s)["future_key"] == "xyz"


class TestApplyMetaHintRoot:
    def test_layer(self):
        s = make_solid()
        apply_meta_hint(s, {"layer": "kinematics"})
        meta = get_solid_metadata(s)
        assert meta["layer"] == "kinematics"

    def test_tags_string(self):
        s = make_solid()
        apply_meta_hint(s, {"tags": "structural"})
        meta = get_solid_metadata(s)
        assert "structural" in meta["tags"]

    def test_tags_list(self):
        s = make_solid()
        apply_meta_hint(s, {"tags": ["structural", "load-bearing"]})
        meta = get_solid_metadata(s)
        assert "structural" in meta["tags"]
        assert "load-bearing" in meta["tags"]

    def test_free_form_root_field(self):
        s = make_solid()
        apply_meta_hint(s, {"description": "Main hinge bracket"})
        meta = get_solid_metadata(s)
        assert meta["description"] == "Main hinge bracket"

    def test_unknown_namespace_verbatim(self):
        s = make_solid()
        apply_meta_hint(s, {"future_ns.some_field": "value"})
        meta = get_solid_metadata(s)
        assert meta.get("future_ns", {}).get("some_field") == "value"


# ---------------------------------------------------------------------------
# Combined / integration
# ---------------------------------------------------------------------------

class TestApplyMetaHintCombined:
    def test_full_stacked_meta(self):
        """Mirror the typical stacked @meta use-case end-to-end."""
        s = make_solid()
        # Simulates:
        #   @meta(assembly.joint_kind="axial", layer="kinematics")
        #   @meta(operation.kind="subtract", operation.feature_kind="pocket")
        hint = {
            "assembly.joint_kind": "axial",
            "layer": "kinematics",
            "operation.kind": "subtract",
            "operation.feature_kind": "pocket",
        }
        apply_meta_hint(s, hint)

        assert get_assembly(s)["joint_kind"] == "axial"
        assert get_operation(s)["kind"] == "subtract"
        assert get_operation(s)["feature_kind"] == "pocket"
        meta = get_solid_metadata(s)
        assert meta["layer"] == "kinematics"

    def test_idempotent(self):
        """Applying the same hint twice does not corrupt metadata."""
        s = make_solid()
        hint = {"assembly.joint_kind": "axial", "operation.kind": "subtract"}
        apply_meta_hint(s, hint)
        apply_meta_hint(s, hint)
        assert get_assembly(s)["joint_kind"] == "axial"
        assert get_operation(s)["kind"] == "subtract"

    def test_empty_hint_is_noop(self):
        s = make_solid()
        apply_meta_hint(s, {})
        # No metadata sections created
        meta = get_solid_metadata(s)
        assert "assembly" not in meta
        assert "operation" not in meta

    def test_invalid_meta_hint_type_raises(self):
        s = make_solid()
        with pytest.raises(TypeError):
            apply_meta_hint(s, "not a dict")

    def test_invalid_assembly_enum_raises(self):
        """Invalid enum values propagate from the v1.1 helpers."""
        s = make_solid()
        with pytest.raises(ValueError, match="assembly.joint_kind"):
            apply_meta_hint(s, {"assembly.joint_kind": "revolute"})  # not in v1.1 vocab

    def test_invalid_operation_enum_raises(self):
        s = make_solid()
        with pytest.raises(ValueError, match="operation.kind"):
            apply_meta_hint(s, {"operation.kind": "cut"})  # valid word, wrong enum value


# ---------------------------------------------------------------------------
# apply_meta_hint_to_raw
# ---------------------------------------------------------------------------

class TestApplyMetaHintToRaw:
    def test_raw_dict_assembly(self):
        meta = {}
        apply_meta_hint_to_raw(meta, {"assembly.joint_kind": "axial"})
        assert get_assembly_metadata(meta)["joint_kind"] == "axial"

    def test_raw_dict_operation(self):
        meta = {}
        apply_meta_hint_to_raw(meta, {"operation.kind": "union"})
        assert get_operation_metadata(meta)["kind"] == "union"

    def test_raw_dict_invalid_type_raises(self):
        with pytest.raises(TypeError):
            apply_meta_hint_to_raw({}, 42)


# ---------------------------------------------------------------------------
# DSL round-trip: parse @meta → apply_meta_hint
# ---------------------------------------------------------------------------

class TestDslRoundTrip:
    def test_parse_and_apply(self):
        """Full round-trip: DSL source → parse → meta_hint → apply → assert."""
        import sys, pathlib
        sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))
        from yapcad.dsl.lexer import Lexer
        from yapcad.dsl.parser import Parser

        src = """
module test

@meta(assembly.joint_kind="axial", layer="kinematics")
@meta(operation.kind="subtract", operation.feature_kind="pocket")
command MAKE_POCKET(depth: float = 10.0) -> solid:
    emit box(1.0, 1.0, 1.0)
"""
        ast = Parser(Lexer(src).tokenize()).parse_module()
        fn = ast.functions[0]
        assert fn.meta_hint is not None

        s = make_solid()
        apply_meta_hint(s, fn.meta_hint)

        assert get_assembly(s)["joint_kind"] == "axial"
        assert get_operation(s)["kind"] == "subtract"
        assert get_operation(s)["feature_kind"] == "pocket"
        meta = get_solid_metadata(s)
        assert meta["layer"] == "kinematics"

    def test_parse_and_apply_list_fields(self):
        """DSL round-trip exercising @meta list-of-dict assembly fields.

        This is the headline new capability: ``@meta(assembly.datums=[...],
        assembly.bolt_patterns=[...])`` parses through dict literals to a
        flat hint dict, then ``apply_meta_hint`` dispatches each entry to
        the dedicated ``yapcad.metadata`` helper.
        """
        from yapcad.dsl.lexer import Lexer
        from yapcad.dsl.parser import Parser

        src = '''
module test

@meta(
    assembly.joint_kind="axial",
    assembly.datums=[
        {"id": "bore_axis", "kind": "axis"},
        {"id": "neck_circle", "kind": "bolt_circle", "R_mm": 151.9, "z_mm": 330.0},
    ],
    assembly.bolt_patterns=[
        {"id": "primary", "ring": "axial", "R_mm": 149.4, "z_mm": 317.3, "count": 8},
    ],
    assembly.surfaces=[
        {"id": "mate_top", "kind": "mating", "mate_to": "nosecone.base_bore"},
    ],
)
command FORWARD_BULKHEAD() -> solid:
    emit box(1.0, 1.0, 1.0)
'''
        ast = Parser(Lexer(src).tokenize()).parse_module()
        fn = ast.functions[0]
        assert fn.meta_hint is not None

        s = make_solid()
        apply_meta_hint(s, fn.meta_hint)

        asm = get_assembly(s)
        assert asm["joint_kind"] == "axial"

        datums = asm["datums"]
        assert [d["id"] for d in datums] == ["bore_axis", "neck_circle"]
        assert datums[1]["kind"] == "bolt_circle"
        assert datums[1]["R_mm"] == 151.9

        bp = asm["bolt_patterns"][0]
        assert bp["id"] == "primary"
        assert bp["ring"] == "axial"
        assert bp["count"] == 8

        surface = asm["surfaces"][0]
        assert surface["id"] == "mate_top"
        assert surface["kind"] == "mating"
        assert surface["mate_to"] == "nosecone.base_bore"
