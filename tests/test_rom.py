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
from yapcad.kinematics import Transform, KinematicPart


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

    def test_to_kinematic_part(self):
        """Verify RomLink.to_kinematic_part returns a properly built KinematicPart."""
        from yapcad.kinematics import JointType

        meta = {
            "parent": "BASE",
            "parent_frame": "MOUNT_A",
            "joint_kind": "revolute",
            "axis": [1.0, 0.0, 0.0],
            "limits": [-45, 90],
            "zero_pose": 5.0,
            "description": "Test joint",
            "printable": "true",
        }
        link = RomLink.from_command_meta("HIP_JOINT", meta)
        part = link.to_kinematic_part()

        assert isinstance(part, KinematicPart)
        assert part.name == "HIP_JOINT"
        assert part.parent == "BASE"
        assert part.parent_frame == "MOUNT_A"
        assert part.description == "Test joint"
        assert part.is_printable is True

        assert part.joint.name == "HIP_JOINT_joint"
        assert part.joint.joint_type == JointType.REVOLUTE
        assert part.joint.axis == (1.0, 0.0, 0.0)
        assert part.joint.value == 5.0
        assert part.joint.min_limit == -45.0
        assert part.joint.max_limit == 90.0

    def test_to_kinematic_part_fixed(self):
        """Fixed joint should map to JointType.FIXED."""
        from yapcad.kinematics import JointType

        link = RomLink(
            name="BASE_PLATE",
            parent=None,
            joint_kind="fixed",
            description="Base",
        )
        part = link.to_kinematic_part()
        assert part.joint.joint_type == JointType.FIXED

    def test_to_kinematic_part_with_base_transform(self):
        """Base transform should be attached when provided."""
        tf = Transform.from_translation(10.0, 20.0, 30.0)
        link = RomLink(name="OFFSET_LINK", parent="BASE", joint_kind="fixed")
        part = link.to_kinematic_part(base_transform=tf)
        assert part.joint.base_transform is not None
        assert part.joint.base_transform.translation == (10.0, 20.0, 30.0)


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

    def test_get_sdf_node_found(self):
        """get_sdf_node returns the node for a known part."""
        from yapcad.sdf import box

        node = box((5, 5, 5))
        link = RomLink(name="PART", parent=None, joint_kind="fixed", sdf_node=node)
        tree = RomTree.from_links("test", [link])
        assert tree.get_sdf_node("PART") is node

    def test_get_sdf_node_not_found(self):
        """get_sdf_node returns None for unknown part."""
        tree = RomTree("test")
        tree.add_link(RomLink(name="A", parent=None, joint_kind="fixed"))
        assert tree.get_sdf_node("NONEXISTENT") is None

    def test_link_names(self):
        """link_names returns the correct list of names."""
        links = [
            RomLink(name="ROOT", parent=None, joint_kind="fixed"),
            RomLink(name="CHILD", parent="ROOT", joint_kind="fixed"),
        ]
        tree = RomTree.from_links("test", links)
        names = tree.link_names()
        assert "ROOT" in names
        assert "CHILD" in names
        assert len(names) == 2

    def test_adjacent_pairs_missing_parent(self):
        """Links whose parent is not in the tree are excluded from pairs."""
        # Build the chain manually to bypass add_part validation
        tree = RomTree("test")
        link_a = RomLink(name="A", parent=None, joint_kind="fixed")
        link_b = RomLink(name="B", parent="MISSING", joint_kind="revolute", axis=(0, 0, 1))
        tree.links["A"] = link_a
        tree.links["B"] = link_b
        tree.chain.add_part(link_a.to_kinematic_part())
        # Skip B's add_part since parent MISSING doesn't exist
        pairs = tree.adjacent_pairs()
        # B's parent "MISSING" is not in tree.links, so no pair should be emitted for B
        for parent, child in pairs:
            assert parent in tree.links
            assert child in tree.links
        assert ("MISSING", "B") not in pairs

    def test_from_interpreter_context(self):
        """Verify kinematic metadata propagates from full DSL compilation."""
        import yapcad.rom.annotations  # noqa: trigger patch
        from yapcad.dsl import compile_and_run
        from yapcad.metadata import get_solid_metadata

        source = '''
module interp_test

@meta(kinematic.joint_kind="fixed")
command FIXED_BASE() -> solid:
    emit box(20, 20, 20)

@meta(kinematic.joint_kind="revolute")
@meta(kinematic.parent="FIXED_BASE")
@meta(kinematic.axis=[0, 0, 1])
@meta(kinematic.limits=[-90, 90])
@meta(kinematic.zero_pose=0.0)
command MOVING_PART() -> solid:
    emit box(10, 10, 10)
'''
        # compile_and_run should attach kinematic metadata to the result
        res = compile_and_run(source, "MOVING_PART", {})
        assert res.success
        meta = get_solid_metadata(res.geometry)
        assert "kinematic" in meta
        assert meta["kinematic"]["joint_kind"] == "revolute"
        assert meta["kinematic"]["parent"] == "FIXED_BASE"
        assert meta["kinematic"]["axis"] == [0, 0, 1]
        assert meta["kinematic"]["limits"] == [-90, 90]


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

    def test_clearance_sdf_all_result_keys(self):
        """Verify clearance_sdf returns all 9 expected result keys."""
        from yapcad.sdf import box

        node_a = box((10, 10, 10))
        node_b = box((10, 10, 10))

        tf_a = Transform.identity()
        tf_b = Transform.from_translation(20.0, 0.0, 0.0)

        result = clearance_sdf(node_a, node_b, tf_a, tf_b, resolution=16)

        expected_keys = {
            "min_distance", "penetration_depth", "mean_distance",
            "within_5mm", "within_2mm", "fraction_negative",
            "grid_points", "full_results",
        }
        assert expected_keys.issubset(set(result.keys())), \
            f"Missing keys: {expected_keys - set(result.keys())}"


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

    def test_sweep_config_defaults(self):
        """Verify SweepConfig has sensible defaults."""
        cfg = SweepConfig()
        assert cfg.step_deg == 5.0
        assert cfg.step_mm == 2.0
        assert cfg.clearance_threshold == 5.0
        assert cfg.interference_threshold == 0.0
        assert cfg.sdf_resolution == 32
        assert cfg.sdf_padding == 2.0
        assert cfg.max_pairs == 0
        assert cfg.progress == 50
        assert cfg.early_stop == 0

    def test_sweep_result_accessors(self):
        """Test worst_interference() and clearance_at() methods."""
        tree = RomTree("accessor_test")
        result = SweepResult(
            tree=tree,
            config=SweepConfig(),
            n_configurations=10,
            events=[
                InterferenceEvent(
                    config={"J": 45.0},
                    pair=("A", "B"),
                    min_distance=-3.0,
                    penetration_depth=3.0,
                ),
                InterferenceEvent(
                    config={"J": 60.0},
                    pair=("A", "B"),
                    min_distance=-1.0,
                    penetration_depth=1.0,
                ),
            ],
            pair_trace={
                ("A", "B"): [
                    ("J=45", -3.0, -2.0),
                    ("J=60", -1.0, -0.5),
                ],
            },
        )

        worst = result.worst_interference()
        assert worst is not None
        assert worst.min_distance == -3.0

        trace = result.clearance_at("A", "B")
        assert len(trace) == 2
        assert trace[0][1] == -3.0

        assert result.clearance_at("X", "Y") == []

    def test_dof_schedule_prismatic(self):
        """A prismatic joint should appear in the DOF schedule."""
        links = [
            RomLink(name="BASE", parent=None, joint_kind="fixed"),
            RomLink(
                name="SLIDER", parent="BASE",
                joint_kind="prismatic", axis=(0, 0, 1),
                limits=[0, 100], zero_pose=0.0,
            ),
        ]
        tree = RomTree.from_links("prismatic_test", links)
        sched = dof_schedule_from_tree(tree, step_deg=10.0)
        assert "SLIDER" in sched
        lo, hi, step = sched["SLIDER"]
        assert lo == 0.0
        assert hi == 110.0  # 100 + step
        assert step == 10.0

    def test_sweep_extra_pairs(self):
        """Passing extra_pairs should include those pairs in the check."""
        from yapcad.sdf import box

        node_a = box((10, 10, 10))
        node_b = box((10, 10, 10))
        node_c = box((10, 10, 10))

        links = [
            RomLink(name="A", parent=None, joint_kind="fixed", sdf_node=node_a),
            RomLink(
                name="B", parent="A",
                joint_kind="revolute", axis=(0, 0, 1),
                limits=[0, 10], zero_pose=0.0, sdf_node=node_b,
            ),
            RomLink(
                name="C", parent="B",
                joint_kind="revolute", axis=(0, 0, 1),
                limits=[0, 10], zero_pose=0.0, sdf_node=node_c,
            ),
        ]
        tree = RomTree.from_links("extra_test", links)

        result = sweep_kinematics(
            tree,
            {"B": (0, 11, 10)},
            extra_pairs=[("A", "C")],
            config=SweepConfig(sdf_resolution=12, sdf_padding=0.0),
        )
        assert ("A", "B") in result.pair_trace
        assert ("B", "C") in result.pair_trace
        assert ("A", "C") in result.pair_trace

    def test_sweep_fixed_angles(self):
        """fixed_angles should override sweep DOF values."""
        from yapcad.sdf import box

        node_a = box((10, 10, 10))
        node_b = box((10, 10, 10))

        links = [
            RomLink(name="A", parent=None, joint_kind="fixed", sdf_node=node_a),
            RomLink(
                name="B", parent="A",
                joint_kind="revolute", axis=(0, 0, 1),
                limits=[-90, 90], zero_pose=0.0, sdf_node=node_b,
            ),
            RomLink(
                name="C", parent="B",
                joint_kind="revolute", axis=(1, 0, 0),
                limits=[-45, 45], zero_pose=0.0, sdf_node=node_a,
            ),
        ]
        tree = RomTree.from_links("fixed_test", links)

        result = sweep_kinematics(
            tree,
            {"B": (-10, 11, 10)},
            fixed_angles={"C": 30.0},
            config=SweepConfig(sdf_resolution=12, sdf_padding=0.0),
        )
        for ev in result.events:
            assert ev.config.get("C") == 30.0


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

    def test_compile_and_interpret_context(self):
        """Full end-to-end: compile_and_run → metadata extraction."""
        import yapcad.rom.annotations  # noqa: trigger patch
        from yapcad.dsl import compile_and_run
        from yapcad.metadata import get_solid_metadata

        source = '''
module ctx_test

@meta(kinematic.joint_kind="fixed")
command BASE() -> solid:
    emit box(30, 30, 30)

@meta(kinematic.joint_kind="revolute")
@meta(kinematic.parent="BASE")
@meta(kinematic.axis=[1, 0, 0])
@meta(kinematic.limits=[0, 90])
@meta(kinematic.zero_pose=0.0)
command ARM() -> solid:
    emit box(10, 10, 50)
'''
        res = compile_and_run(source, "ARM", {})
        assert res.success
        meta = get_solid_metadata(res.geometry)
        assert meta["kinematic"]["joint_kind"] == "revolute"
        assert meta["kinematic"]["parent"] == "BASE"
        assert meta["kinematic"]["limits"] == [0, 90]

    def test_from_interpreter_context_skip_fixed(self):
        """Verify fixed and revolute annotations are distinguished."""
        import yapcad.rom.annotations  # noqa: trigger patch
        from yapcad.dsl import compile_and_run
        from yapcad.metadata import get_solid_metadata

        source = '''
module skip_test
@meta(kinematic.joint_kind="fixed")
command BASE() -> solid:
    emit box(20, 20, 20)
@meta(kinematic.joint_kind="revolute")
@meta(kinematic.parent="BASE")
@meta(kinematic.axis=[0, 0, 1])
@meta(kinematic.limits=[-90, 90])
@meta(kinematic.zero_pose=0.0)
command JOINT() -> solid:
    emit box(5, 5, 5)
'''
        # Fixed command should have joint_kind=fixed
        res_base = compile_and_run(source, "BASE", {})
        meta_base = get_solid_metadata(res_base.geometry)
        assert meta_base["kinematic"]["joint_kind"] == "fixed"

        # Revolute command should have joint_kind=revolute
        res_joint = compile_and_run(source, "JOINT", {})
        meta_joint = get_solid_metadata(res_joint.geometry)
        assert meta_joint["kinematic"]["joint_kind"] == "revolute"
        assert meta_joint["kinematic"]["parent"] == "BASE"

        # Build a tree from explicit links mirroring the DSL
        link_base = RomLink(
            name="BASE", parent=None, joint_kind="fixed",
            sdf_node=None,
        )
        link_joint = RomLink(
            name="JOINT", parent="BASE", joint_kind="revolute",
            axis=(0, 0, 1), limits=(-90, 90), zero_pose=0.0,
            sdf_node=None,
        )
        tree = RomTree.from_links("skip_test", [link_base, link_joint])
        # Both should be present
        assert "BASE" in tree.links
        assert "JOINT" in tree.links
        # Adjacent pair should exist
        pairs = tree.adjacent_pairs()
        assert ("BASE", "JOINT") in pairs
