"""Tests for the DSL assembly builtins (Phase 2 of yapcad-assembly-integration).

Exercises the five Phase-2 builtins ``assembly``, ``add_part``, ``add_mate``,
``validate_assembly``, ``assembly_report`` end-to-end through the DSL
interpreter, plus targeted unit tests for the metadata→PartDefinition
lifter that ``add_part`` uses internally.
"""

import pytest

from yapcad.assembly.assembly import Assembly
from yapcad.assembly.datum import Datum, DatumType, PartDefinition
from yapcad.assembly.mate import Mate, MateType
from yapcad.dsl.runtime.builtins import call_builtin, get_builtin_registry
from yapcad.dsl.runtime.values import (
    assembly_val,
    solid_val,
    string_val,
)
from yapcad.dsl.meta_apply import apply_meta_hint


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_solid_with_assembly_meta(meta_hint=None):
    """Return a solid stub with optional @meta-style metadata applied."""
    solid = ['solid', [], [], [], {}]
    if meta_hint:
        apply_meta_hint(solid, meta_hint)
    return solid


# ---------------------------------------------------------------------------
# 1. ``assembly(name)`` — constructor
# ---------------------------------------------------------------------------

class TestAssemblyBuiltin:
    def test_assembly_creates_handle(self):
        result = call_builtin("assembly", [string_val("forward_section")])
        assert result.type.name == "assembly"
        assert isinstance(result.data, Assembly)
        assert result.data.name == "forward_section"

    def test_assembly_name_coerced_to_string(self):
        # Pass a string explicitly; numeric coercion is the call site's job
        result = call_builtin("assembly", [string_val("123")])
        assert result.data.name == "123"


# ---------------------------------------------------------------------------
# 2. ``add_part(asm, solid, name)`` — instance + datum lifting
# ---------------------------------------------------------------------------

class TestAddPartBuiltin:
    def test_add_part_without_datums(self):
        """Bare solid with no metadata still adds a part instance."""
        asm = call_builtin("assembly", [string_val("test")])
        solid = make_solid_with_assembly_meta()
        result = call_builtin(
            "add_part",
            [asm, solid_val(solid), string_val("bulkhead")],
        )
        # Builtin returns the same assembly handle (chainable)
        assert result is asm
        assembly_obj = asm.data
        assert "bulkhead" in assembly_obj.parts
        part = assembly_obj.parts["bulkhead"]
        assert part.name == "bulkhead"
        assert part.datums == {}

    def test_add_part_lifts_axis_datum_from_metadata(self):
        """``@meta(assembly.datums=[{kind=axis, ...}])`` → Datum on PartDef."""
        solid = make_solid_with_assembly_meta({
            "assembly.datums": [
                {
                    "id": "bore_axis",
                    "kind": "axis",
                    "R_mm": 0.0,
                    "z_mm": 0.0,
                    "direction": [0.0, 0.0, 1.0],
                }
            ]
        })
        asm = call_builtin("assembly", [string_val("test")])
        call_builtin(
            "add_part",
            [asm, solid_val(solid), string_val("part_a")],
        )
        part = asm.data.parts["part_a"]
        assert "bore_axis" in part.datums
        datum = part.datums["bore_axis"]
        assert datum.datum_type == DatumType.AXIS
        assert datum.direction == [0.0, 0.0, 1.0, 0.0]

    def test_add_part_lifts_bolt_circle_datum(self):
        """``kind=bolt_circle`` lifts to a CIRCLE datum with radius+normal."""
        solid = make_solid_with_assembly_meta({
            "assembly.datums": [
                {
                    "id": "neck",
                    "kind": "bolt_circle",
                    "R_mm": 149.4,
                    "z_mm": 317.3,
                    "direction": [0.0, 0.0, 1.0],
                }
            ]
        })
        asm = call_builtin("assembly", [string_val("test")])
        call_builtin(
            "add_part",
            [asm, solid_val(solid), string_val("bulkhead")],
        )
        datum = asm.data.parts["bulkhead"].datums["neck"]
        assert datum.datum_type == DatumType.CIRCLE
        assert datum.radius == 149.4
        # origin is [R, 0, z, 1]
        assert datum.origin == [149.4, 0.0, 317.3, 1.0]
        assert datum.normal == [0.0, 0.0, 1.0, 0.0]

    def test_add_part_lifts_plane_datum(self):
        solid = make_solid_with_assembly_meta({
            "assembly.datums": [
                {
                    "id": "top_face",
                    "kind": "plane",
                    "R_mm": 0.0,
                    "z_mm": 330.0,
                    "direction": [0.0, 0.0, 1.0],
                }
            ]
        })
        asm = call_builtin("assembly", [string_val("test")])
        call_builtin(
            "add_part",
            [asm, solid_val(solid), string_val("p")],
        )
        datum = asm.data.parts["p"].datums["top_face"]
        assert datum.datum_type == DatumType.PLANE
        assert datum.normal == [0.0, 0.0, 1.0, 0.0]

    def test_add_part_skips_unsupported_datum_kind(self):
        """``edge`` has no Datum equivalent; the entry is silently skipped.

        Both ``apply_meta_hint`` (the @meta path) and the direct
        ``add_datum`` helper accept ``edge`` as a valid v1.1 kind, but
        Datum.DatumType has no edge variant. ``_part_def_from_solid``
        must handle this gracefully (skip the entry, keep the rest).
        """
        # Build metadata directly — apply_meta_hint also accepts edge, but
        # we want full control over the entries for this regression test.
        from yapcad.metadata import get_solid_metadata, add_datum
        solid = make_solid_with_assembly_meta()
        meta = get_solid_metadata(solid, create=True)
        add_datum(meta, id="edge_ref", kind="edge")
        add_datum(meta, id="ax", kind="axis")

        asm = call_builtin("assembly", [string_val("test")])
        call_builtin(
            "add_part",
            [asm, solid_val(solid), string_val("p")],
        )
        part = asm.data.parts["p"]
        assert "edge_ref" not in part.datums  # unsupported kind skipped
        assert "ax" in part.datums            # supported kind kept

    def test_add_part_skips_entries_without_id(self):
        """Entries lacking an ``id`` field are skipped by the lifter.

        This guards against partial/garbage data that may live in the
        metadata dict if a future v1.x writer left a sentinel entry.
        We inject the malformed entry directly because ``add_datum``
        and ``apply_meta_hint`` both reject id-less entries upfront.
        """
        from yapcad.metadata import get_solid_metadata, get_assembly_metadata, add_datum
        solid = make_solid_with_assembly_meta()
        meta = get_solid_metadata(solid, create=True)
        # First the well-formed one through the helper (validates)
        add_datum(meta, id="real", kind="axis")
        # Then the malformed one injected directly
        asm_meta = get_assembly_metadata(meta, create=True)
        asm_meta.setdefault("datums", []).append({"kind": "axis"})  # no id

        asm = call_builtin("assembly", [string_val("test")])
        call_builtin(
            "add_part",
            [asm, solid_val(solid), string_val("p")],
        )
        datums = asm.data.parts["p"].datums
        assert "real" in datums
        assert len(datums) == 1  # id-less entry skipped

    def test_add_part_duplicate_name_raises(self):
        """The underlying Assembly raises on duplicate part names."""
        from yapcad.assembly.assembly import AssemblyError
        solid = make_solid_with_assembly_meta()
        asm = call_builtin("assembly", [string_val("test")])
        call_builtin("add_part", [asm, solid_val(solid), string_val("p")])
        with pytest.raises(AssemblyError):
            call_builtin("add_part", [asm, solid_val(solid), string_val("p")])


# ---------------------------------------------------------------------------
# 3. ``add_mate(asm, kind, part_a, datum_a, part_b, datum_b)``
# ---------------------------------------------------------------------------

class TestAddMateBuiltin:
    def _setup_two_parts(self, kind_a="axis", kind_b="axis"):
        """Build an assembly with two parts, each carrying one datum."""
        asm = call_builtin("assembly", [string_val("test")])
        for name, kind in [("part_a", kind_a), ("part_b", kind_b)]:
            solid = make_solid_with_assembly_meta({
                "assembly.datums": [
                    {
                        "id": "shaft",
                        "kind": kind,
                        "direction": [0.0, 0.0, 1.0],
                    }
                ]
            })
            call_builtin(
                "add_part",
                [asm, solid_val(solid), string_val(name)],
            )
        return asm

    def test_add_mate_concentric(self):
        asm = self._setup_two_parts()
        result = call_builtin(
            "add_mate",
            [
                asm,
                string_val("concentric"),
                string_val("part_a"),
                string_val("shaft"),
                string_val("part_b"),
                string_val("shaft"),
            ],
        )
        # Chainable
        assert result is asm
        mates = asm.data.mates
        assert len(mates) == 1
        m = mates[0]
        assert m.mate_type == MateType.CONCENTRIC
        assert m.part_a == "part_a"
        assert m.datum_a == "shaft"
        assert m.part_b == "part_b"
        assert m.datum_b == "shaft"

    def test_add_mate_revolute(self):
        asm = self._setup_two_parts()
        call_builtin(
            "add_mate",
            [
                asm,
                string_val("revolute"),
                string_val("part_a"),
                string_val("shaft"),
                string_val("part_b"),
                string_val("shaft"),
            ],
        )
        assert asm.data.mates[0].mate_type == MateType.REVOLUTE

    def test_add_mate_invalid_kind_raises(self):
        asm = self._setup_two_parts()
        with pytest.raises(ValueError, match="add_mate: invalid kind"):
            call_builtin(
                "add_mate",
                [
                    asm,
                    string_val("not_a_real_mate_kind"),
                    string_val("part_a"),
                    string_val("shaft"),
                    string_val("part_b"),
                    string_val("shaft"),
                ],
            )

    def test_add_mate_includes_valid_kinds_in_error(self):
        """The error message lists valid kinds so DSL authors can recover."""
        asm = self._setup_two_parts()
        try:
            call_builtin(
                "add_mate",
                [
                    asm,
                    string_val("not_real"),
                    string_val("part_a"),
                    string_val("shaft"),
                    string_val("part_b"),
                    string_val("shaft"),
                ],
            )
        except ValueError as exc:
            msg = str(exc)
            # The message lists at least the basic mate kinds
            for expected in ("concentric", "revolute", "coincident"):
                assert expected in msg, f"missing {expected!r} in error: {msg}"


# ---------------------------------------------------------------------------
# 4. ``validate_assembly(asm)`` / 5. ``assembly_report(asm)``
# ---------------------------------------------------------------------------

class TestValidateAndReportBuiltins:
    def test_validate_empty_assembly(self):
        """An empty assembly is trivially valid."""
        asm = call_builtin("assembly", [string_val("empty")])
        result = call_builtin("validate_assembly", [asm])
        assert result.type.name == "bool"
        assert result.data is True

    def test_validate_single_part_assembly(self):
        asm = call_builtin("assembly", [string_val("solo")])
        solid = make_solid_with_assembly_meta()
        call_builtin(
            "add_part",
            [asm, solid_val(solid), string_val("only_part")],
        )
        result = call_builtin("validate_assembly", [asm])
        assert result.data is True

    def test_report_is_string_containing_assembly_name(self):
        asm = call_builtin("assembly", [string_val("forward_section")])
        result = call_builtin("assembly_report", [asm])
        assert result.type.name == "string"
        assert "forward_section" in result.data



# ---------------------------------------------------------------------------
# End-to-end: full assembly construction through builtins
# ---------------------------------------------------------------------------

class TestAssemblyBuiltinsEndToEnd:
    def test_two_part_concentric_assembly(self):
        """Build a two-part assembly with one CONCENTRIC mate, validate it.

        This is the workflow we want DSL authors to be able to write::

            let nosecone = NOSECONE()        # @meta(assembly.datums=[shaft axis])
            let bulkhead = FORWARD_BULKHEAD()  # @meta(assembly.datums=[shaft axis])
            let asm = assembly("forward_section")
            add_part(asm, nosecone, "nosecone")
            add_part(asm, bulkhead, "bulkhead")
            add_mate(asm, "concentric", "nosecone", "shaft", "bulkhead", "shaft")
            assert validate_assembly(asm)
        """
        def for_solid():
            return make_solid_with_assembly_meta({
                "assembly.datums": [
                    {"id": "shaft", "kind": "axis", "direction": [0.0, 0.0, 1.0]}
                ]
            })

        asm = call_builtin("assembly", [string_val("forward_section")])
        call_builtin(
            "add_part",
            [asm, solid_val(for_solid()), string_val("nosecone")],
        )
        call_builtin(
            "add_part",
            [asm, solid_val(for_solid()), string_val("bulkhead")],
        )
        call_builtin(
            "add_mate",
            [
                asm,
                string_val("concentric"),
                string_val("nosecone"),
                string_val("shaft"),
                string_val("bulkhead"),
                string_val("shaft"),
            ],
        )

        valid = call_builtin("validate_assembly", [asm])
        assert valid.data is True

        report = call_builtin("assembly_report", [asm]).data
        assert "forward_section" in report
        assert "nosecone" in report
        assert "bulkhead" in report

    def test_add_mate_fails_fast_on_unknown_part(self):
        """Adding a mate that references a missing part raises immediately.

        The Python ``Assembly.add_mate`` validates part names at add time,
        so dangling references are caught before they reach
        ``validate_assembly``. This is the right behavior: it gives DSL
        authors a clear error at the offending line instead of a deferred
        "validation failed" with no location.
        """
        from yapcad.assembly.assembly import AssemblyError
        solid = make_solid_with_assembly_meta({
            "assembly.datums": [
                {"id": "shaft", "kind": "axis", "direction": [0.0, 0.0, 1.0]}
            ]
        })
        asm = call_builtin("assembly", [string_val("broken")])
        call_builtin(
            "add_part",
            [asm, solid_val(solid), string_val("only_one")],
        )
        with pytest.raises(AssemblyError, match="phantom"):
            call_builtin(
                "add_mate",
                [
                    asm,
                    string_val("concentric"),
                    string_val("only_one"),
                    string_val("shaft"),
                    string_val("phantom"),
                    string_val("shaft"),
                ],
            )
