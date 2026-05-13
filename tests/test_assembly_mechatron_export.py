"""Tests for the Mechatron assembly-graph JSON exporter.

Two test families:

1. ``to_mechatron_snapshot`` / ``to_mechatron_json`` directly \u2014 unit tests
   for the converter, covering empty assemblies, parts with datums of
   each supported kind, all four mate-to-joint mappings, materials,
   and the round-trip shape.
2. The ``emit_assembly`` DSL builtin \u2014 exercises the integration
   path the way a DSL author would.
"""

import json

import pytest

from yapcad.assembly.assembly import Assembly
from yapcad.assembly.datum import Datum, DatumType, PartDefinition
from yapcad.assembly.mate import Mate, MateType
from yapcad.assembly.mechatron_export import (
    to_mechatron_json,
    to_mechatron_snapshot,
)
from yapcad.dsl.runtime.builtins import call_builtin
from yapcad.dsl.runtime.values import solid_val, string_val
from yapcad.dsl.meta_apply import apply_meta_hint


# ---------------------------------------------------------------------------
# Helpers (shared with test_dsl_assembly_builtins)
# ---------------------------------------------------------------------------

def make_solid_with_assembly_meta(meta_hint=None):
    solid = ['solid', [], [], [], {}]
    if meta_hint:
        apply_meta_hint(solid, meta_hint)
    return solid


def _bare_part(name: str) -> PartDefinition:
    return PartDefinition(name=name)


# ---------------------------------------------------------------------------
# 1. Direct converter tests
# ---------------------------------------------------------------------------

class TestEmptyAssembly:
    def test_empty_assembly_shape(self):
        asm = Assembly("empty")
        snap = to_mechatron_snapshot(asm)
        assert snap["parts"] == []
        assert snap["interfaces"] == []
        # Mechatron-only features come out empty
        assert snap["loop_closures"] == []
        assert snap["joint_couplings"] == []
        assert snap["keyframes"] == []
        assert snap["generators"] == []
        assert snap["design_constraints"] == []
        # Provenance hint carries the assembly name
        assert snap["last_script"] == "yapcad-dsl: empty"

    def test_json_is_pretty_printed(self):
        asm = Assembly("empty")
        json_str = to_mechatron_json(asm)
        # Pretty-printed JSON has newlines + indentation
        assert "\n" in json_str
        # Round-trips through json.loads to the same shape
        assert json.loads(json_str) == to_mechatron_snapshot(asm)


class TestPartExport:
    def test_part_with_no_datums(self):
        asm = Assembly("a")
        asm.add_part(_bare_part("solo"), name="solo")
        snap = to_mechatron_snapshot(asm)
        assert len(snap["parts"]) == 1
        part = snap["parts"][0]
        assert part["id"] == "solo"
        assert part["name"] == "solo"
        assert part["datums"] == []
        # Defaults
        assert part["material"] == "PETG"  # yapcad PartDefinition default
        assert part["process"] == "FDM"
        assert part["tags"] == []

    def test_part_material_mapping(self):
        asm = Assembly("a")
        part = _bare_part("p")
        part.material = "PETG-CF"
        asm.add_part(part, name="p")
        snap = to_mechatron_snapshot(asm)
        assert snap["parts"][0]["material"] == "PetgCf"

    def test_unknown_material_falls_back_to_petg(self):
        asm = Assembly("a")
        part = _bare_part("p")
        part.material = "Some Future Material 9000"
        asm.add_part(part, name="p")
        snap = to_mechatron_snapshot(asm)
        assert snap["parts"][0]["material"] == "PETG"

    def test_parts_emit_in_insertion_order(self):
        asm = Assembly("a")
        for name in ("nosecone", "bulkhead", "tank"):
            asm.add_part(_bare_part(name), name=name)
        ids = [p["id"] for p in to_mechatron_snapshot(asm)["parts"]]
        assert ids == ["nosecone", "bulkhead", "tank"]


class TestDatumExport:
    """Cover each DatumType → Mechatron datum_type PascalCase mapping."""

    def _build_part(self, datum: Datum) -> Assembly:
        asm = Assembly("a")
        part = _bare_part("p")
        part.add_datum(datum)
        asm.add_part(part, name="p")
        return asm

    def test_axis_datum(self):
        asm = self._build_part(Datum(
            name="shaft",
            datum_type=DatumType.AXIS,
            origin=[0.0, 0.0, 0.0, 1.0],
            direction=[0.0, 0.0, 1.0, 0.0],
        ))
        d = to_mechatron_snapshot(asm)["parts"][0]["datums"][0]
        assert d["name"] == "shaft"
        assert d["datum_type"] == "Axis"
        assert d["origin"] == {"x": 0.0, "y": 0.0, "z": 0.0}
        assert d["direction"] == {"x": 0.0, "y": 0.0, "z": 1.0}
        assert "normal" not in d
        assert "radius" not in d

    def test_plane_datum(self):
        asm = self._build_part(Datum(
            name="top",
            datum_type=DatumType.PLANE,
            origin=[0.0, 0.0, 330.0, 1.0],
            normal=[0.0, 0.0, 1.0, 0.0],
        ))
        d = to_mechatron_snapshot(asm)["parts"][0]["datums"][0]
        assert d["datum_type"] == "Plane"
        assert d["origin"] == {"x": 0.0, "y": 0.0, "z": 330.0}
        assert d["normal"] == {"x": 0.0, "y": 0.0, "z": 1.0}
        assert "direction" not in d

    def test_circle_datum_carries_radius(self):
        asm = self._build_part(Datum(
            name="bolt_circle",
            datum_type=DatumType.CIRCLE,
            origin=[149.4, 0.0, 317.3, 1.0],
            normal=[0.0, 0.0, 1.0, 0.0],
            radius=149.4,
        ))
        d = to_mechatron_snapshot(asm)["parts"][0]["datums"][0]
        assert d["datum_type"] == "Circle"
        assert d["radius"] == 149.4

    def test_point_datum(self):
        asm = self._build_part(Datum(
            name="origin",
            datum_type=DatumType.POINT,
            origin=[10.0, 20.0, 30.0, 1.0],
        ))
        d = to_mechatron_snapshot(asm)["parts"][0]["datums"][0]
        assert d["datum_type"] == "Point"
        assert d["origin"] == {"x": 10.0, "y": 20.0, "z": 30.0}

    def test_homogeneous_w_component_is_dropped(self):
        """yapCAD's 4th component (w=1 for points, w=0 for vectors) drops out."""
        asm = self._build_part(Datum(
            name="ax",
            datum_type=DatumType.AXIS,
            origin=[1.0, 2.0, 3.0, 1.0],
            direction=[0.0, 0.0, 1.0, 0.0],
        ))
        d = to_mechatron_snapshot(asm)["parts"][0]["datums"][0]
        # Point3D / Vec3D objects have only x/y/z keys, no w
        assert set(d["origin"].keys()) == {"x", "y", "z"}
        assert set(d["direction"].keys()) == {"x", "y", "z"}


class TestMateToJointMapping:
    """Cover each yapCAD MateType → Mechatron Joint variant mapping."""

    def _two_part_assembly_with_mate(self, mate_type: MateType) -> Assembly:
        asm = Assembly("a")
        for name in ("part_a", "part_b"):
            part = _bare_part(name)
            part.add_datum(Datum(
                name="shaft",
                datum_type=DatumType.AXIS,
                origin=[0.0, 0.0, 0.0, 1.0],
                direction=[0.0, 0.0, 1.0, 0.0],
            ))
            asm.add_part(part, name=name)
        asm.add_mate(Mate(
            name="test",
            mate_type=mate_type,
            part_a="part_a",
            datum_a="shaft",
            part_b="part_b",
            datum_b="shaft",
        ))
        return asm

    def test_rigid_maps_to_fixed(self):
        asm = self._two_part_assembly_with_mate(MateType.RIGID)
        iface = to_mechatron_snapshot(asm)["interfaces"][0]
        assert iface["joint"] == {"type": "Fixed"}
        # Rigid is a motion mate captured by Joint::Fixed; no mate_constraint
        assert iface["mate_constraints"] == []

    def test_revolute_maps_to_revolute_with_axis(self):
        asm = self._two_part_assembly_with_mate(MateType.REVOLUTE)
        iface = to_mechatron_snapshot(asm)["interfaces"][0]
        assert iface["joint"]["type"] == "Revolute"
        # Placeholder +Z axis until the Mechatron solver overrides it
        assert iface["joint"]["axis"] == [0.0, 0.0, 1.0]

    def test_prismatic_maps_to_prismatic_with_axis(self):
        asm = self._two_part_assembly_with_mate(MateType.PRISMATIC)
        iface = to_mechatron_snapshot(asm)["interfaces"][0]
        assert iface["joint"]["type"] == "Prismatic"
        assert iface["joint"]["axis"] == [0.0, 0.0, 1.0]

    def test_spherical_maps_to_ball(self):
        asm = self._two_part_assembly_with_mate(MateType.SPHERICAL)
        iface = to_mechatron_snapshot(asm)["interfaces"][0]
        assert iface["joint"] == {"type": "Ball"}

    def test_concentric_emits_axis_coincident(self):
        # yapCAD CONCENTRIC maps to Mechatron AxisCoincident (not
        # Concentric) so the validator's composite rule
        # "AxisCoincident + Coincident -> axial rotation locked" can fire
        # and Fixed mounts come out fully constrained.
        asm = self._two_part_assembly_with_mate(MateType.CONCENTRIC)
        iface = to_mechatron_snapshot(asm)["interfaces"][0]
        assert iface["joint"] == {"type": "Fixed"}
        c = iface["mate_constraints"][0]
        assert c["type"] == "AxisCoincident"
        assert c["parent_datum"] == "shaft"
        assert c["child_datum"] == "shaft"

    def test_coincident_emits_mate_constraint(self):
        asm = self._two_part_assembly_with_mate(MateType.COINCIDENT)
        iface = to_mechatron_snapshot(asm)["interfaces"][0]
        assert iface["joint"] == {"type": "Fixed"}
        c = iface["mate_constraints"][0]
        assert c["type"] == "Coincident"
        assert c["parent_datum"] == "shaft"
        assert c["child_datum"] == "shaft"

    def test_tangent_emits_mate_constraint(self):
        # TANGENT requires CIRCLE or PLANE datums per yapcad's validator,
        # so build the fixture with CIRCLE datums for this case only.
        asm = Assembly("a")
        for name in ("part_a", "part_b"):
            part = _bare_part(name)
            part.add_datum(Datum(
                name="face",
                datum_type=DatumType.CIRCLE,
                origin=[0.0, 0.0, 0.0, 1.0],
                normal=[0.0, 0.0, 1.0, 0.0],
                radius=10.0,
            ))
            asm.add_part(part, name=name)
        asm.add_mate(Mate(
            name="test", mate_type=MateType.TANGENT,
            part_a="part_a", datum_a="face",
            part_b="part_b", datum_b="face",
        ))
        c = to_mechatron_snapshot(asm)["interfaces"][0]["mate_constraints"][0]
        assert c["type"] == "Tangent"
        assert c["parent_datum"] == "face"
        assert c["child_datum"] == "face"

    def test_parallel_maps_to_axis_align(self):
        """PARALLEL is now mapped to Mechatron's AxisAlign (unilateral)."""
        asm = self._two_part_assembly_with_mate(MateType.PARALLEL)
        iface = to_mechatron_snapshot(asm)["interfaces"][0]
        assert iface["joint"] == {"type": "Fixed"}
        c = iface["mate_constraints"][0]
        assert c["type"] == "AxisAlign"
        # AxisAlign has a different shape: child_datum + target_axis
        # (no parent_datum). target_axis comes from datum_a's direction.
        assert c["child_datum"] == "shaft"
        assert c["target_axis"] == [0.0, 0.0, 1.0]
        assert "parent_datum" not in c

    def test_perpendicular_drops_constraint(self):
        """PERPENDICULAR still has no Mechatron equivalent and is dropped."""
        asm = self._two_part_assembly_with_mate(MateType.PERPENDICULAR)
        iface = to_mechatron_snapshot(asm)["interfaces"][0]
        assert iface["joint"] == {"type": "Fixed"}
        assert iface["mate_constraints"] == []

    def test_two_mates_same_pair_collapse_to_one_interface(self):
        """yapCAD authors typically stack Concentric+Coincident to fully
        constrain a rigid mount. Mechatron treats them as one Interface."""
        asm = Assembly("a")
        for name in ("part_a", "part_b"):
            part = _bare_part(name)
            part.add_datum(Datum(
                name="centerline", datum_type=DatumType.AXIS,
                origin=[0.0, 0.0, 0.0, 1.0],
                direction=[0.0, 0.0, 1.0, 0.0],
            ))
            part.add_datum(Datum(
                name="face", datum_type=DatumType.PLANE,
                origin=[0.0, 0.0, 0.0, 1.0],
                normal=[0.0, 0.0, 1.0, 0.0],
            ))
            asm.add_part(part, name=name)
        asm.add_mate(Mate(
            name="m1", mate_type=MateType.CONCENTRIC,
            part_a="part_a", datum_a="centerline",
            part_b="part_b", datum_b="centerline",
        ))
        asm.add_mate(Mate(
            name="m2", mate_type=MateType.COINCIDENT,
            part_a="part_a", datum_a="face",
            part_b="part_b", datum_b="face",
        ))
        snap = to_mechatron_snapshot(asm)
        # Two mates -> one Interface, two constraints
        assert len(snap["interfaces"]) == 1
        iface = snap["interfaces"][0]
        assert iface["joint"] == {"type": "Fixed"}
        kinds = [c["type"] for c in iface["mate_constraints"]]
        assert kinds == ["AxisCoincident", "Coincident"]

    def test_motion_and_constraint_mates_share_one_interface(self):
        """A Revolute joint plus a Coincident positioning constraint on the
        same pair must produce one Interface with Joint::Revolute AND
        the Coincident mate_constraint."""
        asm = Assembly("a")
        for name in ("part_a", "part_b"):
            part = _bare_part(name)
            part.add_datum(Datum(
                name="shaft", datum_type=DatumType.AXIS,
                origin=[0.0, 0.0, 0.0, 1.0],
                direction=[1.0, 0.0, 0.0, 0.0],
            ))
            part.add_datum(Datum(
                name="face", datum_type=DatumType.PLANE,
                origin=[0.0, 0.0, 0.0, 1.0],
                normal=[0.0, 0.0, 1.0, 0.0],
            ))
            asm.add_part(part, name=name)
        asm.add_mate(Mate(
            name="contact", mate_type=MateType.COINCIDENT,
            part_a="part_a", datum_a="face",
            part_b="part_b", datum_b="face",
        ))
        asm.add_mate(Mate(
            name="hinge", mate_type=MateType.REVOLUTE,
            part_a="part_a", datum_a="shaft",
            part_b="part_b", datum_b="shaft",
        ))
        snap = to_mechatron_snapshot(asm)
        assert len(snap["interfaces"]) == 1
        iface = snap["interfaces"][0]
        # Motion mate wins the joint slot...
        assert iface["joint"]["type"] == "Revolute"
        assert iface["joint"]["axis"] == [1.0, 0.0, 0.0]
        # ...and the Coincident shows up as a mate_constraint
        assert len(iface["mate_constraints"]) == 1
        assert iface["mate_constraints"][0]["type"] == "Coincident"

    def test_insertion_vector_from_plane_datum(self):
        """insertion_vector should derive from the parent-side PLANE datum
        normal, negated (= direction the child moves to seat)."""
        asm = Assembly("a")
        for name in ("part_a", "part_b"):
            part = _bare_part(name)
            part.add_datum(Datum(
                name="face", datum_type=DatumType.PLANE,
                origin=[0.0, 0.0, 0.0, 1.0],
                normal=[0.0, 0.0, 1.0, 0.0],  # parent face points +Z
            ))
            asm.add_part(part, name=name)
        asm.add_mate(Mate(
            name="contact", mate_type=MateType.COINCIDENT,
            part_a="part_a", datum_a="face",
            part_b="part_b", datum_b="face",
        ))
        iface = to_mechatron_snapshot(asm)["interfaces"][0]
        # Child seats by moving along -Z (into the parent face)
        assert iface["insertion_vector"] == [0.0, 0.0, -1.0]

    def test_no_insertion_vector_when_no_plane_or_circle(self):
        """AXIS-only interfaces leave insertion_vector unset."""
        asm = Assembly("a")
        for name in ("part_a", "part_b"):
            part = _bare_part(name)
            part.add_datum(Datum(
                name="shaft", datum_type=DatumType.AXIS,
                origin=[0.0, 0.0, 0.0, 1.0],
                direction=[0.0, 0.0, 1.0, 0.0],
            ))
            asm.add_part(part, name=name)
        asm.add_mate(Mate(
            name="c", mate_type=MateType.CONCENTRIC,
            part_a="part_a", datum_a="shaft",
            part_b="part_b", datum_b="shaft",
        ))
        iface = to_mechatron_snapshot(asm)["interfaces"][0]
        assert "insertion_vector" not in iface

    def test_interface_ids_are_sequential(self):
        asm = Assembly("a")
        for name in ("part_a", "part_b", "part_c"):
            part = _bare_part(name)
            part.add_datum(Datum(
                name="shaft",
                datum_type=DatumType.AXIS,
                origin=[0.0, 0.0, 0.0, 1.0],
                direction=[0.0, 0.0, 1.0, 0.0],
            ))
            asm.add_part(part, name=name)
        asm.add_mate(Mate(
            name="m1", mate_type=MateType.CONCENTRIC,
            part_a="part_a", datum_a="shaft",
            part_b="part_b", datum_b="shaft",
        ))
        asm.add_mate(Mate(
            name="m2", mate_type=MateType.RIGID,
            part_a="part_b", datum_a="shaft",
            part_b="part_c", datum_b="shaft",
        ))
        ids = [i["id"] for i in to_mechatron_snapshot(asm)["interfaces"]]
        assert ids == ["iface_0", "iface_1"]


# ---------------------------------------------------------------------------
# 2. DSL builtin tests
# ---------------------------------------------------------------------------

class TestEmitAssemblyBuiltin:
    def test_emit_returns_pretty_json_string(self):
        asm = call_builtin("assembly", [string_val("forward_section")])
        result = call_builtin("emit_assembly", [asm])
        assert result.type.name == "string"
        # Returns valid JSON
        snap = json.loads(result.data)
        assert snap["parts"] == []
        assert snap["interfaces"] == []
        assert snap["last_script"] == "yapcad-dsl: forward_section"

    def test_emit_includes_parts_with_datums_from_meta(self):
        """End-to-end: build via builtins, emit, validate JSON shape."""
        solid = make_solid_with_assembly_meta({
            "assembly.datums": [
                {
                    "id": "shaft",
                    "kind": "axis",
                    "direction": [0.0, 0.0, 1.0],
                },
                {
                    "id": "neck",
                    "kind": "bolt_circle",
                    "R_mm": 149.4,
                    "z_mm": 317.3,
                    "direction": [0.0, 0.0, 1.0],
                },
            ],
        })
        asm = call_builtin("assembly", [string_val("test")])
        call_builtin(
            "add_part",
            [asm, solid_val(solid), string_val("bulkhead")],
        )
        snap = json.loads(call_builtin("emit_assembly", [asm]).data)
        assert len(snap["parts"]) == 1
        part = snap["parts"][0]
        assert part["id"] == "bulkhead"
        datums_by_name = {d["name"]: d for d in part["datums"]}
        assert datums_by_name["shaft"]["datum_type"] == "Axis"
        assert datums_by_name["neck"]["datum_type"] == "Circle"
        assert datums_by_name["neck"]["radius"] == 149.4

    def test_emit_full_two_part_assembly_via_builtins(self):
        """Run the entire workflow through DSL builtins, then verify JSON."""
        def for_solid():
            return make_solid_with_assembly_meta({
                "assembly.datums": [
                    {"id": "shaft", "kind": "axis", "direction": [0.0, 0.0, 1.0]},
                ],
            })

        asm = call_builtin("assembly", [string_val("forward_section")])
        call_builtin("add_part", [asm, solid_val(for_solid()), string_val("nosecone")])
        call_builtin("add_part", [asm, solid_val(for_solid()), string_val("bulkhead")])
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

        snap = json.loads(call_builtin("emit_assembly", [asm]).data)
        assert [p["id"] for p in snap["parts"]] == ["nosecone", "bulkhead"]
        iface = snap["interfaces"][0]
        assert iface["parent_part"] == "nosecone"
        assert iface["child_part"] == "bulkhead"
        assert iface["joint"] == {"type": "Fixed"}
        c = iface["mate_constraints"][0]
        # CONCENTRIC -> AxisCoincident (see _MATE_TO_CONSTRAINT_VARIANT)
        assert c["type"] == "AxisCoincident"
        assert c["parent_datum"] == "shaft"
        assert c["child_datum"] == "shaft"

    def test_emit_output_is_deterministic(self):
        """Same input → byte-identical output (no dict-iter randomness)."""
        asm = call_builtin("assembly", [string_val("deterministic")])
        solid = make_solid_with_assembly_meta({
            "assembly.datums": [
                {"id": "shaft", "kind": "axis", "direction": [0.0, 0.0, 1.0]},
            ],
        })
        for name in ("a", "b", "c"):
            call_builtin("add_part", [asm, solid_val(solid), string_val(name)])
        first = call_builtin("emit_assembly", [asm]).data
        second = call_builtin("emit_assembly", [asm]).data
        assert first == second
