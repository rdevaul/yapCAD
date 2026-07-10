"""Tests for LoadCase, BoltPattern, and the Mechatron round-trip emission.

History: added 2026-05-20 alongside `src/yapcad/assembly/load_case.py` to
close the FEA setup gap exposed by the Agentic-1 chute-deploy iterations
(v4-v7). See Whiteboard message 20260520201704-46e1d134 for context.
"""
import json

import pytest

from yapcad.assembly import (
    Assembly,
    AssemblyError,
    PartDefinition,
    LoadCase,
    LoadAttach,
    BoltPattern,
    COORDINATE_FRAME_VEHICLE,
    COORDINATE_FRAME_MESH,
    STATUS_CONFIRMED,
    GROUP_RECOVERY,
    GROUP_INERTIAL,
)
from yapcad.assembly.mechatron_export import (
    to_mechatron_snapshot,
    to_mechatron_json,
)


# ---------------------------------------------------------------------------
# BoltPattern dataclass
# ---------------------------------------------------------------------------

class TestBoltPattern:
    def test_valid_construction(self):
        bp = BoltPattern(pcd_mm=276.0, n_bolts=12, bolt_spec="1-4-20-x-1in-button-head-shcs",
                         access_direction=[0, 0, -1])
        assert bp.pcd_mm == 276.0
        assert bp.n_bolts == 12

    def test_rejects_zero_pcd(self):
        with pytest.raises(ValueError, match="pcd_mm must be positive"):
            BoltPattern(pcd_mm=0, n_bolts=12, bolt_spec="x", access_direction=[0, 0, 1])

    def test_rejects_zero_n_bolts(self):
        with pytest.raises(ValueError, match="n_bolts must be positive"):
            BoltPattern(pcd_mm=100, n_bolts=0, bolt_spec="x", access_direction=[0, 0, 1])

    def test_rejects_bad_access_direction(self):
        with pytest.raises(ValueError, match="access_direction must be 3-vector"):
            BoltPattern(pcd_mm=100, n_bolts=4, bolt_spec="x", access_direction=[0, 0])

    def test_to_dict_roundtrip(self):
        bp = BoltPattern(pcd_mm=276.0, n_bolts=12, bolt_spec="1-4-20-x-1in-button-head-shcs",
                         access_direction=[0, 0, -1])
        d = bp.to_dict()
        assert d == {
            "pcd_mm": 276.0,
            "n_bolts": 12,
            "bolt_spec": "1-4-20-x-1in-button-head-shcs",
            "access_direction": [0, 0, -1],
        }
        bp2 = BoltPattern.from_dict(d)
        assert bp2.pcd_mm == bp.pcd_mm
        assert bp2.n_bolts == bp.n_bolts


# ---------------------------------------------------------------------------
# LoadCase dataclass
# ---------------------------------------------------------------------------

class TestLoadCase:
    def test_valid_construction(self):
        lc = LoadCase(
            id="LC-004a",
            name="Y-bridle leg (LOX aft skirt, clock 90°)",
            group=GROUP_RECOVERY,
            direction=[0, 0, 1],
            attach=LoadAttach(part="oxidizer-tank", position=[1664.6, 0, 145],
                              interface="iface_intertank_lox_tank", clock_deg=90.0),
            magnitude_n=6000.0,
            coordinate_frame=COORDINATE_FRAME_VEHICLE,
            source="dashboard LOAD_CASES",
            status=STATUS_CONFIRMED,
        )
        assert lc.id == "LC-004a"
        assert lc.attach.part == "oxidizer-tank"
        assert lc.magnitude_n == 6000.0

    def test_rejects_unknown_group(self):
        with pytest.raises(ValueError, match="must be one of"):
            LoadCase(
                id="X", name="X", group="bogus", direction=[1, 0, 0],
                attach=LoadAttach(part="p"), magnitude_n=100,
            )

    def test_rejects_no_magnitude(self):
        with pytest.raises(ValueError, match="at least one of magnitude_n / magnitude_g"):
            LoadCase(
                id="X", name="X", group=GROUP_RECOVERY, direction=[1, 0, 0],
                attach=LoadAttach(part="p"),
            )

    def test_rejects_bad_direction(self):
        with pytest.raises(ValueError, match="direction must be 3-vector"):
            LoadCase(
                id="X", name="X", group=GROUP_RECOVERY, direction=[1, 0],
                attach=LoadAttach(part="p"), magnitude_n=100,
            )

    def test_to_dict_roundtrip(self):
        lc = LoadCase(
            id="LC-004a", name="Y-bridle", group=GROUP_RECOVERY,
            direction=[0, 0, 1],
            attach=LoadAttach(part="oxidizer-tank", clock_deg=90.0),
            magnitude_n=6000,
            warn=True,
        )
        d = lc.to_dict()
        assert d["id"] == "LC-004a"
        assert d["attach"]["part"] == "oxidizer-tank"
        assert d["warn"] is True
        lc2 = LoadCase.from_dict(d)
        assert lc2.id == lc.id
        assert lc2.attach.part == lc.attach.part


# ---------------------------------------------------------------------------
# Assembly load_case / bolt_pattern registration
# ---------------------------------------------------------------------------

def _make_minimal_assembly():
    asm = Assembly("test")
    asm.add_part(PartDefinition(name="boattail"))
    asm.add_part(PartDefinition(name="upper-thrust"))
    asm.add_part(PartDefinition(name="oxidizer-tank"))
    return asm


class TestAssemblyLoadCaseRegistration:
    def test_add_load_case_succeeds_for_valid_part(self):
        asm = _make_minimal_assembly()
        lc = LoadCase(
            id="LC-004a", name="Y-bridle", group=GROUP_RECOVERY,
            direction=[0, 0, 1],
            attach=LoadAttach(part="oxidizer-tank", clock_deg=90.0),
            magnitude_n=6000,
        )
        asm.add_load_case(lc)
        assert "LC-004a" in asm.load_cases
        assert asm.get_load_case("LC-004a") is lc

    def test_add_load_case_rejects_unknown_part(self):
        asm = _make_minimal_assembly()
        lc = LoadCase(
            id="X", name="X", group=GROUP_RECOVERY, direction=[1, 0, 0],
            attach=LoadAttach(part="nonexistent"), magnitude_n=100,
        )
        with pytest.raises(AssemblyError, match="not in the assembly"):
            asm.add_load_case(lc)

    def test_add_load_case_allows_assembly_sentinel(self):
        asm = _make_minimal_assembly()
        # "_assembly" is the special whole-stack body force sentinel
        lc = LoadCase(
            id="LC-003", name="Burnout 6g", group=GROUP_INERTIAL,
            direction=[1, 0, 0],
            attach=LoadAttach(part="_assembly", position=[300, 0, 0]),
            magnitude_g=6, magnitude_n=6500,
        )
        asm.add_load_case(lc)
        assert "LC-003" in asm.load_cases

    def test_add_load_case_rejects_duplicate_id(self):
        asm = _make_minimal_assembly()
        lc = LoadCase(
            id="LC-1", name="A", group=GROUP_RECOVERY, direction=[1, 0, 0],
            attach=LoadAttach(part="oxidizer-tank"), magnitude_n=100,
        )
        asm.add_load_case(lc)
        with pytest.raises(AssemblyError, match="already registered"):
            asm.add_load_case(lc)


class TestAssemblyBoltPatternRegistration:
    def test_add_bolt_pattern(self):
        asm = _make_minimal_assembly()
        bp = BoltPattern(pcd_mm=284, n_bolts=12, bolt_spec="1-4-20",
                         access_direction=[0, 0, 1])
        asm.add_bolt_pattern("boattail", "upper-thrust", bp)
        assert asm.get_bolt_pattern("boattail", "upper-thrust") is bp

    def test_rejects_unknown_part(self):
        asm = _make_minimal_assembly()
        bp = BoltPattern(pcd_mm=100, n_bolts=4, bolt_spec="x", access_direction=[0, 0, 1])
        with pytest.raises(AssemblyError, match="not in assembly"):
            asm.add_bolt_pattern("boattail", "nonexistent", bp)

    def test_rejects_duplicate(self):
        asm = _make_minimal_assembly()
        bp = BoltPattern(pcd_mm=100, n_bolts=4, bolt_spec="x", access_direction=[0, 0, 1])
        asm.add_bolt_pattern("boattail", "upper-thrust", bp)
        with pytest.raises(AssemblyError, match="already registered"):
            asm.add_bolt_pattern("boattail", "upper-thrust", bp)


# ---------------------------------------------------------------------------
# Mechatron snapshot emission
# ---------------------------------------------------------------------------

class TestMechatronSnapshotEmission:
    def test_bolt_pattern_emitted_on_interface(self):
        """Verify that BoltPattern registered on an assembly is emitted onto
        the corresponding Mechatron Interface in the snapshot output."""
        from yapcad.assembly.mate import Mate, MateType
        asm = _make_minimal_assembly()
        bp = BoltPattern(pcd_mm=284, n_bolts=12, bolt_spec="1-4-20-x-1in-button-head-shcs",
                         access_direction=[0, 0, 1])
        asm.add_bolt_pattern("boattail", "upper-thrust", bp)
        # Add a minimal mate so an interface is created
        mate = Mate(name="m1", mate_type=MateType.RIGID,
                    part_a="boattail", part_b="upper-thrust",
                    datum_a="dummy_a", datum_b="dummy_b")
        # We bypass the datum check by appending directly; the exporter only
        # needs the part_a/part_b to be set for interface grouping.
        asm.mates.append(mate)
        snap = to_mechatron_snapshot(asm)
        assert len(snap["interfaces"]) == 1
        iface = snap["interfaces"][0]
        assert iface["parent_part"] == "boattail"
        assert iface["child_part"] == "upper-thrust"
        assert "bolt_pattern" in iface
        assert iface["bolt_pattern"]["pcd_mm"] == 284
        assert iface["bolt_pattern"]["n_bolts"] == 12
        assert iface["bolt_pattern"]["bolt_spec"] == "1-4-20-x-1in-button-head-shcs"

    def test_load_cases_emitted_on_snapshot(self):
        asm = _make_minimal_assembly()
        lc = LoadCase(
            id="LC-004a", name="Y-bridle", group=GROUP_RECOVERY,
            direction=[0, 0, 1],
            attach=LoadAttach(part="oxidizer-tank", clock_deg=90.0,
                              interface="iface_intertank_lox_tank"),
            magnitude_n=6000,
            warn=True,
        )
        asm.add_load_case(lc)
        snap = to_mechatron_snapshot(asm)
        assert "load_cases" in snap
        assert len(snap["load_cases"]) == 1
        emitted = snap["load_cases"][0]
        assert emitted["id"] == "LC-004a"
        assert emitted["attach"]["part"] == "oxidizer-tank"
        assert emitted["attach"]["clock_deg"] == 90.0
        assert emitted["warn"] is True

    def test_no_load_cases_means_no_field(self):
        """Backward-compat: assemblies without load_cases shouldn't have
        an empty load_cases key in the snapshot (avoids polluting older
        consumers that don't expect it)."""
        asm = _make_minimal_assembly()
        snap = to_mechatron_snapshot(asm)
        assert "load_cases" not in snap

    def test_no_bolt_pattern_means_no_field(self):
        from yapcad.assembly.mate import Mate, MateType
        asm = _make_minimal_assembly()
        mate = Mate(name="m1", mate_type=MateType.RIGID,
                    part_a="boattail", part_b="upper-thrust",
                    datum_a="x", datum_b="y")
        asm.mates.append(mate)
        snap = to_mechatron_snapshot(asm)
        iface = snap["interfaces"][0]
        assert "bolt_pattern" not in iface

    def test_full_snapshot_json_roundtrip(self):
        """End-to-end: build an assembly with LoadCase + BoltPattern, emit
        JSON, parse it back, verify shape."""
        from yapcad.assembly.mate import Mate, MateType
        asm = _make_minimal_assembly()
        bp = BoltPattern(pcd_mm=276, n_bolts=12, bolt_spec="1-4-20-x-1in-button-head-shcs",
                         access_direction=[0, 0, -1])
        asm.add_bolt_pattern("oxidizer-tank", "boattail", bp)  # arbitrary pair
        mate = Mate(name="m1", mate_type=MateType.RIGID,
                    part_a="oxidizer-tank", part_b="boattail",
                    datum_a="a", datum_b="b")
        asm.mates.append(mate)
        lc = LoadCase(
            id="LC-004a", name="Y-bridle", group=GROUP_RECOVERY,
            direction=[0, 0, 1],
            attach=LoadAttach(part="oxidizer-tank"),
            magnitude_n=6000,
        )
        asm.add_load_case(lc)
        s = to_mechatron_json(asm)
        parsed = json.loads(s)
        assert parsed["interfaces"][0]["bolt_pattern"]["pcd_mm"] == 276
        assert parsed["load_cases"][0]["id"] == "LC-004a"
