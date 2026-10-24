"""Tests for yapcad.rom package.

Covers:
- DSL annotation parsing (@meta(kinematic.*))
- Tree building from compiled DSL
- Clearance evaluation between adjacent links
- ROM sweep orchestration
"""

from __future__ import annotations

import pytest
import numpy as np

from yapcad.rom import (
    RomTree,
    RomLink,
    sweep_kinematics,
    SweepResult,
    SweepConfig,
    InterferenceEvent,
    dof_schedule_from_tree,
    interference_summary,
)
from yapcad.rom.annotations import apply_kinematic_hint
from yapcad.rom.clearance import clearance_sdf
from yapcad.kinematics import Transform


# ===========================================================================
# Annotations
# ===========================================================================

class TestAnnotations:
    """Test kinematic hint routing into metadata."""

    def test_basic_hint_routing(self):
        """Apply a kinematic hint dict and verify it lands in metadata."""
        meta = {}
        apply_kinematic_hint(meta, "joint_kind", "revolute")
        apply_kinematic_hint(meta, "parent", "BASE")
        apply_kinematic_hint(meta, "axis", [0, 0, 1])
        apply_kinematic_hint(meta, "limits", [-90, 90])
        apply_kinematic_hint(meta, "zero_pose", 0.0)

        kin = meta.get("kinematic", {})
        assert kin["joint_kind"] == "revolute"
        assert kin["parent"] == "BASE"
        assert kin["axis"] == [0.0, 0.0, 1.0]
        assert kin["limits"] == [-90.0, 90.0]
        assert kin["zero_pose"] == 0.0

    def test_printable_bool(self):
        meta = {}
        apply_kinematic_hint(meta, "printable", "true")
        assert meta["kinematic"]["printable"] is True
        apply_kinematic_hint(meta, "printable", False)
        assert meta["kinematic"]["printable"] is False


# ===========================================================================
# RomLink
# ===========================================================================

class TestRomLink:
    def test_from_meta_dict(self):
        meta = {
            "parent": "PELVIS",
            "parent_frame": "HIP_YAW",
            "joint_kind": "revolute",
            "axis": [0, 0, 1],
            "limits": [-90, 90],
            "zero_pose": 0.0,
            "description": "Hip yaw joint",
        }
        link = RomLink.from_command_meta("HIP_YAW", meta)
        assert link.name == "HIP_YAW"
        assert link.parent == "PELVIS"
        assert link.parent_frame == "HIP_YAW"
        assert link.joint_kind == "revolute"
        assert link.axis == (0, 0, 1)
        assert link.limits == (-90.0, 90.0)
        assert link.description == "Hip yaw joint"


# ===========================================================================
# RomTree
# ===========================================================================

class TestRomTree:
    def test_build_from_links(self):
        """Build a tree from explicit RomLink objects."""
        links = [
            RomLink(name="BASE", parent=None, joint_kind="fixed"),
            RomLink(
                name="KNEE", parent="BASE",
                joint_kind="revolute", axis=(1, 0, 0),
                limits=[0, 135], zero_pose=0.0,
            ),
        ]
        tree = RomTree.from_links("test", links)
        assert set(tree.links.keys()) == {"BASE", "KNEE"}
        pairs = tree.adjacent_pairs()
        assert ("BASE", "KNEE") in pairs
        assert tree.chain.get_joint_value("KNEE") == 0.0

    def test_set_joints(self):
        links = [
            RomLink(name="A", parent=None, joint_kind="fixed"),
            RomLink(
                name="B", parent="A",
                joint_kind="revolute", axis=(0, 0, 1),
                limits=[-90, 90], zero_pose=0.0,
            ),
        ]
        tree = RomTree.from_links("test", links)
        tree.set_joints({"B": 45.0})
        assert tree.chain.get_joint_value("B") == 45.0
        tf = tree.get_world_transform("B")
        assert tf.translation == (0.0, 0.0, 0.0)

    def test_clear_joints(self):
        links = [
            RomLink(name="A", parent=None, joint_kind="fixed"),
            RomLink(
                name="B", parent="A",
                joint_kind="revolute", axis=(0, 0, 1),
                limits=[-90, 90], zero_pose=15.0,
            ),
        ]
        tree = RomTree.from_links("test", links)
        tree.set_joints({"B": 45.0})
        assert tree.chain.get_joint_value("B") == 45.0
        tree.clear_joints()
        assert tree.chain.get_joint_value("B") == 15.0


# ===========================================================================
# Clearance
# ===========================================================================

class TestClearance:
    def test_clearance_overlapping_boxes(self):
        """Two boxes that overlap should report negative clearance."""
        from yapcad.sdf import box, translate

        node_a = box((10, 10, 10))
        node_b = box((10, 10, 10))

        # B is shifted by 8 mm in X → 2 mm overlap
        from yapcad.kinematics import Transform
        tf_b = Transform.from_translation(8.0, 0.0, 0.0)
        tf_a = Transform.identity()

        result = clearance_sdf(node_a, node_b, tf_a, tf_b, resolution=20)
        assert result["min_distance"] < 0, \
            f"Expected overlap, got min_distance={result['min_distance']:.2f}"
        assert result["penetration_depth"] > 0

    def test_clearance_separated_boxes(self):
        """Two boxes far apart should report positive clearance."""
        from yapcad.sdf import box
        from yapcad.kinematics import Transform

        node_a = box((10, 10, 10))
        node_b = box((10, 10, 10))

        tf_a = Transform.identity()
        tf_b = Transform.from_translation(50.0, 0.0, 0.0)

        result = clearance_sdf(node_a, node_b, tf_a, tf_b, resolution=20)
        assert result["min_distance"] > 0, \
            f"Expected separation, got min_distance={result['min_distance']:.2f}"

    def test_clearance_transformed_box(self):
        """B shifted by 12 mm → 2 mm gap (box edge to edge at X=5 and X=7)."""
        from yapcad.sdf import box
        from yapcad.kinematics import Transform

        node_a = box((10, 10, 10))
        node_b = box((10, 10, 10))

        tf_a = Transform.identity()
        tf_b = Transform.from_translation(12.0, 0.0, 0.0)

        result = clearance_sdf(node_a, node_b, tf_a, tf_b, resolution=24)
        # 12 mm separation, each box is 10 mm → 2 mm gap at surface
        # At resolution=24, the min_distance should be ≥ 0 (no overlap)
        assert result["min_distance"] >= 0
        assert result["min_distance"] < 4  # should be close to 2 mm


# ===========================================================================
# Sweep
# ===========================================================================

class TestSweep:
    def test_dof_schedule_from_tree(self):
        links = [
            RomLink(name="A", parent=None, joint_kind="fixed"),
            RomLink(
                name="KNEE", parent="A",
                joint_kind="revolute", axis=(1, 0, 0),
                limits=[0, 135], zero_pose=0.0,
            ),
        ]
        tree = RomTree.from_links("test", links)
        sched = dof_schedule_from_tree(tree, step_deg=10.0)
        assert "KNEE" in sched
        lo, hi, step = sched["KNEE"]
        assert lo == 0.0
        assert step == 10.0
        # hi should be 135 + 10 = 145 (so arange goes up to 140)
        assert hi == 145.0

    def test_sweep_no_sdf_skips_clearance(self):
        """Sweep with no SDF nodes should still produce a result with 0 pairs checked."""
        links = [
            RomLink(name="A", parent=None, joint_kind="fixed"),
            RomLink(
                name="KNEE", parent="A",
                joint_kind="revolute", axis=(1, 0, 0),
                limits=[0, 90], zero_pose=0.0,
            ),
        ]
        tree = RomTree.from_links("test", links)
        result = sweep_kinematics(
            tree,
            {"KNEE": (0, 91, 15)},
        )
        assert result.n_configurations == 7  # 0, 15, 30, 45, 60, 75, 90
        assert len(result.events) == 0
        assert "KNEE" not in str(result.summary()) or "min" not in str(result.summary())

    def test_sweep_with_sdf_overlapping(self):
        """Sweep with overlapping SDF boxes should detect interference."""
        from yapcad.sdf import box, translate

        node_a = box((10, 10, 10))
        node_b = box((10, 10, 10))

        links = [
            RomLink(
                name="A", parent=None, joint_kind="fixed",
                sdf_node=node_a,
            ),
            RomLink(
                name="B", parent="A",
                joint_kind="revolute", axis=(0, 0, 1),
                limits=[-30, 30], zero_pose=0.0,
                sdf_node=node_b,
            ),
        ]
        tree = RomTree.from_links("overlap_test", links)

        result = sweep_kinematics(
            tree,
            {"B": (-30, 31, 10)},
            config=SweepConfig(sdf_resolution=20, sdf_padding=1.0),
        )
        assert result.n_configurations > 0
        # B and A overlap (they're both centered at origin)
        assert len(result.events) > 0, \
            f"Expected interference events, got {len(result.events)}"

    def test_interference_event_data(self):
        """Verify that interference events carry useful data."""
        from yapcad.sdf import box

        node_a = box((10, 10, 10))
        node_b = box((10, 10, 10))

        links = [
            RomLink(name="A", parent=None, joint_kind="fixed", sdf_node=node_a),
            RomLink(
                name="B", parent="A",
                joint_kind="revolute", axis=(0, 0, 1),
                limits=[-10, 10], zero_pose=0.0,
                sdf_node=node_b,
            ),
        ]
        tree = RomTree.from_links("evt_test", links)

        result = sweep_kinematics(
            tree,
            {"B": (-10, 11, 10)},
            config=SweepConfig(sdf_resolution=16, sdf_padding=0.0),
        )

        if result.events:
            ev = result.events[0]
            assert isinstance(ev, InterferenceEvent)
            assert isinstance(ev.config, dict)
            assert ev.pair == ("A", "B")
            assert ev.min_distance < 0
            assert ev.penetration_depth > 0
            assert "min_distance" in ev.pairwise_results

    def test_summary_output(self):
        result = SweepResult(
            tree=RomTree("summary_test"),
            config=SweepConfig(step_deg=5.0),
            n_configurations=42,
            events=[
                InterferenceEvent(
                    config={"J1": 45.0, "J2": 90.0},
                    pair=("A", "B"),
                    min_distance=-2.5,
                    penetration_depth=2.5,
                ),
            ],
        )
        txt = result.summary()
        assert "42" in txt
        assert "1" in txt  # interference count
        assert "-2.5" in txt

    def test_report_output(self):
        from yapcad.sdf import box

        node_a = box((10, 10, 10))
        node_b = box((10, 10, 10))
        links = [
            RomLink(name="A", parent=None, joint_kind="fixed", sdf_node=node_a),
            RomLink(
                name="B", parent="A",
                joint_kind="revolute", axis=(0, 0, 1),
                limits=[-10, 10], zero_pose=0.0,
                sdf_node=node_b,
            ),
        ]
        tree = RomTree.from_links("report_test", links)

        result = sweep_kinematics(
            tree,
            {"B": (-10, 11, 10)},
            config=SweepConfig(sdf_resolution=16, sdf_padding=0.0),
        )

        report = interference_summary(result)
        assert "ROM Sweep Report" in report
        assert "report_test" in report
        # Should mention pair clearance minima
        assert "Pair" in report or "Clearance" in report


# ===========================================================================
# DSL integration (optional, requires full interpreter)
# ===========================================================================

class TestDSLIntegration:
    """End-to-end test: DSL → RomTree → sweep."""

    def test_compile_and_build_tree(self):
        """Verify that @meta(kinematic.*) annotations are preserved and parsed."""
        import yapcad.rom.annotations  # noqa: trigger patch
        from yapcad.dsl import compile_and_run
        from yapcad.metadata import get_solid_metadata

        source = '''
module test_mod

@meta(kinematic.joint_kind=\"revolute\")
@meta(kinematic.parent=\"BASE\")
@meta(kinematic.axis=[0, 0, 1])
@meta(kinematic.limits=[-90, 90])
@meta(kinematic.zero_pose=0.0)
command TEST_JOINT(x: float = 0.0, y: float = 0.0, z: float = 0.0) -> solid:
    emit box(10, 10, 10)

@meta(kinematic.joint_kind=\"fixed\")
command BASE() -> solid:
    emit box(20, 20, 20)
'''
        # The type checker should now accept kinematic namespace
        from yapcad.dsl import tokenize, parse, check
        tokens = tokenize(source)
        mod = parse(tokens)
        chk = check(mod)
        errors = [d for d in chk.diagnostics if d.severity.name == 'error']
        assert not errors, f"Type check errors: {errors}"

        # compile_and_run should attach kinematic metadata
        res = compile_and_run(source, "TEST_JOINT", {"x": 0.0, "y": 0.0, "z": 0.0})
        assert res.success
        meta = get_solid_metadata(res.geometry)
        assert "kinematic" in meta
        assert meta["kinematic"]["joint_kind"] == "revolute"
        assert meta["kinematic"]["parent"] == "BASE"
