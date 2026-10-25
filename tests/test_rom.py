"""Range-of-motion sweeps: clearance measurement and assembly sweeps.

The clearance tests check measured numbers against analytic answers,
because a clearance tool that is merely self-consistent is worthless: an
earlier version of this code reported every distance short by exactly its
sampling margin, and passed its own tests.

The sweep tests drive a three-part fixture whose arm swings towards a stop
it shares no mate with -- the pair where interference actually happens.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from yapcad.dsl import compile_and_run
from yapcad.rom import (
    Clearance,
    RomError,
    RomJoint,
    RomLink,
    RomModel,
    SweepConfig,
    assembly_from_result,
    bounds_clearance,
    clearance_sdf,
    interference_summary,
    radian_range,
    schedule_from_model,
    sweep_rom,
    world_bounds,
)
from yapcad.sdf import box, half_space, sphere


def tf(x: float = 0.0, y: float = 0.0, z: float = 0.0) -> np.ndarray:
    matrix = np.eye(4)
    matrix[:3, 3] = (x, y, z)
    return matrix


def rot_z(degrees: float) -> np.ndarray:
    angle = math.radians(degrees)
    matrix = np.eye(4)
    matrix[0, 0] = math.cos(angle)
    matrix[0, 1] = -math.sin(angle)
    matrix[1, 0] = math.sin(angle)
    matrix[1, 1] = math.cos(angle)
    return matrix


# ===========================================================================
# Clearance, against answers worked out by hand
# ===========================================================================

class TestClearanceAccuracy:
    """Two 10 mm cubes, one offset along X: the gap is offset - 10."""

    @pytest.mark.parametrize("offset,expected", [
        (20.0, 10.0),
        (12.0, 2.0),
        (11.0, 1.0),
        (10.5, 0.5),
        (10.0, 0.0),
        (9.5, -0.5),
        (9.0, -1.0),
        (8.0, -2.0),
    ])
    def test_distance_matches_the_analytic_gap(self, offset, expected):
        result = clearance_sdf(box((10, 10, 10)), box((10, 10, 10)),
                               tf(), tf(offset), pitch_mm=0.5)
        assert result.min_distance == pytest.approx(
            expected, abs=result.tolerance_mm
        )

    @pytest.mark.parametrize("pitch", [2.0, 1.0, 0.5, 0.25])
    def test_the_pitch_does_not_bias_the_result(self, pitch):
        """The regression test for this module's worst bug.

        Sampling used to run over the part's bounding box grown by a
        margin, and kept the points in that margin -- so every clearance
        came back short by exactly the margin, at every resolution, and
        a pair 1 mm apart was reported as interfering.
        """
        result = clearance_sdf(box((10, 10, 10)), box((10, 10, 10)),
                               tf(), tf(12.0), pitch_mm=pitch)
        assert result.min_distance == pytest.approx(
            2.0, abs=result.tolerance_mm
        )
        assert result.min_distance > 0.0, "a 2 mm gap is not interference"

    def test_a_pair_one_millimetre_apart_does_not_interfere(self):
        result = clearance_sdf(box((10, 10, 10)), box((10, 10, 10)),
                               tf(), tf(11.0))
        assert not result.interferes
        assert result.penetration == 0.0

    def test_a_one_millimetre_overlap_does_interfere(self):
        result = clearance_sdf(box((10, 10, 10)), box((10, 10, 10)),
                               tf(), tf(9.0), pitch_mm=0.25)
        assert result.interferes
        assert result.penetration == pytest.approx(
            1.0, abs=result.tolerance_mm
        )

    def test_spheres_converge_on_the_analytic_gap(self):
        """Curved surfaces, offset off-axis so no sample lands on the
        closest point: the error has to shrink with the pitch."""
        step = 12.0 / math.sqrt(3.0)
        apart = tf(step, step, step)          # centres 12 mm apart
        coarse = clearance_sdf(sphere(5.0), sphere(5.0), tf(), apart,
                               pitch_mm=1.0)
        fine = clearance_sdf(sphere(5.0), sphere(5.0), tf(), apart,
                             pitch_mm=0.2)
        assert abs(coarse.min_distance - 2.0) <= coarse.tolerance_mm
        assert abs(fine.min_distance - 2.0) <= fine.tolerance_mm
        assert fine.tolerance_mm < coarse.tolerance_mm
        assert (abs(fine.min_distance - 2.0)
                < abs(coarse.min_distance - 2.0) + 1e-12)

    def test_the_tolerance_is_proportional_to_the_pitch(self):
        half = clearance_sdf(box((10, 10, 10)), box((10, 10, 10)),
                             tf(), tf(12.0), pitch_mm=0.5)
        quarter = clearance_sdf(box((10, 10, 10)), box((10, 10, 10)),
                                tf(), tf(12.0), pitch_mm=0.25)
        assert half.tolerance_mm == pytest.approx(
            0.5 * math.sqrt(3) / 2, rel=0.2
        )
        assert quarter.tolerance_mm < half.tolerance_mm

    def test_moving_both_parts_together_changes_nothing(self):
        """The measurement is of the pair, not of the world frame."""
        plain = clearance_sdf(box((10, 10, 10)), box((20, 10, 10)),
                              tf(), tf(20.0), pitch_mm=0.5)
        shift = rot_z(37.0) @ tf(13.0, -4.0, 8.0)
        moved = clearance_sdf(box((10, 10, 10)), box((20, 10, 10)),
                              shift, shift @ tf(20.0), pitch_mm=0.5)
        assert moved.min_distance == pytest.approx(
            plain.min_distance, abs=1e-6
        )

    def test_a_rotated_part_is_measured_in_its_rotated_place(self):
        """A 40 mm bar end-on is 15 mm away; turned broadside, 5 mm."""
        bar = box((40, 10, 10))
        end_on = clearance_sdf(box((10, 10, 10)), bar, tf(), tf(40.0),
                               pitch_mm=0.5)
        broadside = clearance_sdf(box((10, 10, 10)), bar,
                                  tf(), tf(40.0) @ rot_z(90.0),
                                  pitch_mm=0.5)
        assert end_on.min_distance == pytest.approx(
            15.0, abs=end_on.tolerance_mm
        )
        assert broadside.min_distance == pytest.approx(
            30.0, abs=broadside.tolerance_mm
        )

    def test_a_part_thinner_than_the_pitch_is_still_measured(self):
        plate = box((10, 10, 0.2))
        result = clearance_sdf(box((10, 10, 10)), plate,
                               tf(), tf(0.0, 0.0, 7.1), pitch_mm=1.0)
        assert result.samples > 0
        assert result.min_distance == pytest.approx(
            2.0, abs=result.tolerance_mm
        )


class TestClearanceLimits:
    def test_an_unbounded_field_is_refused(self):
        with pytest.raises(RomError, match="unbounded"):
            clearance_sdf(box((10, 10, 10)), half_space((0, 0, 1)))

    def test_a_non_positive_pitch_is_refused(self):
        with pytest.raises(RomError, match="pitch_mm"):
            clearance_sdf(box((10, 10, 10)), box((10, 10, 10)),
                          pitch_mm=0.0)

    def test_a_bad_transform_is_refused(self):
        with pytest.raises(RomError, match="4x4"):
            clearance_sdf(box((10, 10, 10)), box((10, 10, 10)),
                          tf(), np.eye(3))

    def test_too_many_samples_coarsens_the_pitch_and_reports_it(self):
        result = clearance_sdf(box((10, 10, 10)), box((10, 10, 10)),
                               tf(), tf(12.0),
                               pitch_mm=0.01, max_samples=20_000)
        assert result.samples <= 20_000
        assert result.pitch_mm > 0.01, "the coarser pitch must be reported"
        assert result.min_distance == pytest.approx(
            2.0, abs=result.tolerance_mm
        )

    def test_a_result_smaller_than_its_tolerance_is_not_resolved(self):
        result = clearance_sdf(box((10, 10, 10)), box((10, 10, 10)),
                               tf(), tf(10.1), pitch_mm=2.0)
        assert not result.resolved
        assert "+/-" in result.describe()


class TestBounds:
    def test_world_bounds_follow_the_pose(self):
        lo, hi = world_bounds(box((10, 10, 10)), tf(5.0, 0.0, 0.0))
        assert lo == pytest.approx((0.0, -5.0, -5.0))
        assert hi == pytest.approx((10.0, 5.0, 5.0))

    def test_bounds_clearance_is_a_lower_bound(self):
        for offset in (12.0, 20.0, 50.0):
            bound = bounds_clearance(box((10, 10, 10)), box((10, 10, 10)),
                                     tf(), tf(offset))
            measured = clearance_sdf(box((10, 10, 10)), box((10, 10, 10)),
                                     tf(), tf(offset), pitch_mm=0.5)
            assert bound <= measured.min_distance + 1e-9
            assert bound == pytest.approx(offset - 10.0, abs=1e-9)

    def test_overlapping_boxes_bound_nothing(self):
        assert bounds_clearance(box((10, 10, 10)), box((10, 10, 10)),
                                tf(), tf(5.0)) == 0.0


# ===========================================================================
# The model, read from a DSL assembly
# ===========================================================================

ARM_DESIGN = """module rom_arm

@meta(assembly.datums=[
    {"id": "pivot", "kind": "axis",
     "origin_mm": [0.0, 0.0, 10.0], "direction": [0.0, 0.0, 1.0]},
    {"id": "stop_mount", "kind": "axis",
     "origin_mm": [0.0, 25.0, 15.0], "direction": [0.0, 0.0, 1.0]}
])
command BASE() -> solid:
    emit box(40.0, 40.0, 20.0)

@meta(assembly.datums=[
    {"id": "hub", "kind": "axis",
     "origin_mm": [0.0, 0.0, 0.0], "direction": [0.0, 0.0, 1.0]}
])
command ARM() -> solid:
    emit translate(box(60.0, 10.0, 10.0), 30.0, 0.0, 15.0)

@meta(assembly.datums=[
    {"id": "mount", "kind": "axis",
     "origin_mm": [0.0, 0.0, 0.0], "direction": [0.0, 0.0, 1.0]}
])
command STOP() -> solid:
    emit box(10.0, 20.0, 30.0)

command BUILD() -> solid:
    let a: assembly = assembly("rom_arm")
    add_part(a, BASE(), "base")
    add_part(a, ARM(), "arm")
    add_part(a, STOP(), "stop")
    add_named_mate(a, "elbow", "revolute", "base", "pivot", "arm", "hub")
    add_named_mate(a, "stop_mount", "rigid",
                   "base", "stop_mount", "stop", "mount")
    set_mate_limits(a, "elbow", radians(-90.0), radians(90.0))
    solve_assembly(a, "base")
    emit assembly_compound(a)
"""

LIMITS_LINE = (
    '    set_mate_limits(a, "elbow", radians(-90.0), radians(90.0))\n'
)

PLAIN_DESIGN = """module rom_plain

command PART() -> solid:
    emit box(10.0, 10.0, 10.0)
"""


def build(source: str, command: str, representation: str, tmp_path):
    path = tmp_path / "design.dsl"
    path.write_text(source)
    return compile_and_run(source, command, {},
                           representation=representation,
                           source_path=str(path))


@pytest.fixture(scope="module")
def field_result(tmp_path_factory):
    """The arm design built as fields -- the supported way to sweep."""
    return build(ARM_DESIGN, "BUILD", "sdf",
                 tmp_path_factory.mktemp("rom_sdf"))


@pytest.fixture
def model(field_result):
    return RomModel.from_result(field_result)


class TestModel:
    def test_parts_fields_joints_and_limits_are_read(self, model):
        assert sorted(model.parts) == ["arm", "base", "stop"]
        assert model.root == "base"
        assert model.parts_without_fields() == {}
        elbow = model.joint("elbow")
        assert elbow.movable and elbow.bounded
        assert elbow.lower == pytest.approx(math.radians(-90.0))
        assert elbow.upper == pytest.approx(math.radians(90.0))
        assert [j.name for j in model.movable_joints()] == ["elbow"]
        assert model.joint("stop_mount").kind == "rigid"

    def test_the_field_is_the_part_not_its_bounding_box(self, model):
        """An L-shaped part must not answer as the box around it."""
        field = model.field("arm")
        lo, hi = field.bounds
        assert lo == pytest.approx((0.0, -5.0, 10.0))
        assert hi == pytest.approx((60.0, 5.0, 20.0))

    def test_pose_rotates_about_the_datum_not_the_origin(self, model):
        """The whole point of reading the datums: the joint is at z=10."""
        straight = model.pose({"elbow": 0.0})
        turned = model.pose({"elbow": math.radians(90.0)})
        assert np.asarray(straight["arm"])[:3, 3] == pytest.approx(
            (0.0, 0.0, 10.0)
        )
        assert np.asarray(turned["arm"])[:3, 3] == pytest.approx(
            (0.0, 0.0, 10.0)
        )
        # +X in the arm's frame points along +Y once it has turned 90 deg.
        assert np.asarray(turned["arm"])[:3, 0] == pytest.approx(
            (0.0, 1.0, 0.0), abs=1e-9
        )
        # The stop is rigid: it does not move with the joint.
        assert np.asarray(straight["stop"])[:3, 3] == pytest.approx(
            np.asarray(turned["stop"])[:3, 3]
        )

    def test_moving_pairs_are_the_ones_the_joint_separates(self, model):
        assert model.moving_pairs() == [("base", "arm"), ("arm", "stop")]
        # arm and stop share no mate at all, which is where the
        # interference is; base and stop never move relative to each other.
        assert model.mate_pairs() == [("base", "arm"), ("base", "stop")]
        assert len(model.all_pairs()) == 3
        assert model.mate_chains() == {
            "base": (), "arm": ("elbow",), "stop": ("stop_mount",),
        }

    def test_an_unknown_joint_is_an_error(self, model):
        with pytest.raises(RomError, match="No joint named 'knee'"):
            model.joint("knee")

    def test_an_unknown_part_is_an_error(self, model):
        with pytest.raises(RomError, match="No part named 'ghost'"):
            model.field("ghost")

    def test_a_value_beyond_the_declared_limit_is_refused(self, model):
        with pytest.raises(RomError, match="above its limit"):
            model.pose({"elbow": math.radians(120.0)})
        with pytest.raises(RomError, match="below its limit"):
            model.pose({"elbow": math.radians(-120.0)})

    def test_a_rigid_joint_cannot_be_driven(self, model):
        with pytest.raises(RomError, match="not revolute"):
            model.pose({"stop_mount": 0.1})

    def test_describe_names_the_joints_and_their_limits(self, model):
        text = model.describe()
        assert "rooted at 'base'" in text
        assert "elbow: revolute, -90.0 deg to 90.0 deg" in text
        assert "stop_mount: rigid" in text


class TestModelWithoutFields:
    def test_a_mesh_build_says_what_to_do_about_it(self, tmp_path):
        result = build(ARM_DESIGN, "BUILD", "mesh", tmp_path)
        model = RomModel.from_result(result)
        missing = model.parts_without_fields()
        assert set(missing) == {"arm", "base", "stop"}
        assert "representation='sdf'" in missing["arm"]
        with pytest.raises(RomError, match="no field"):
            model.field("arm")

    def test_a_result_with_no_assembly_is_an_error(self, tmp_path):
        result = build(PLAIN_DESIGN, "PART", "sdf", tmp_path)
        with pytest.raises(RomError, match="assembly_compound"):
            RomModel.from_result(result)

    def test_a_failed_run_is_an_error(self):
        result = compile_and_run("module broken\n", "NOPE", {})
        with pytest.raises(RomError, match="DSL run failed"):
            RomModel.from_result(result)


class TestValueObjects:
    def test_a_joint_describes_its_range_in_degrees(self):
        joint = RomJoint(name="elbow", parent="base", child="arm",
                         kind="revolute", lower=math.radians(-30.0),
                         upper=math.radians(45.0))
        assert joint.movable and joint.bounded
        assert "elbow: revolute, -30.0 deg to 45.0 deg" in joint.describe()

    def test_an_unbounded_joint_says_so(self):
        joint = RomJoint(name="spin", parent="a", child="b",
                         kind="revolute")
        assert joint.movable and not joint.bounded
        assert "unbounded" in joint.describe()

    def test_a_rigid_joint_is_not_movable(self):
        joint = RomJoint(name="bolted", parent="a", child="b", kind="rigid")
        assert not joint.movable
        assert "rigid" in joint.describe()

    def test_a_link_knows_whether_it_has_a_field(self):
        assert RomLink(name="a", field=box((1, 1, 1))).has_field
        assert not RomLink(name="b", reason="imported mesh").has_field

    def test_assembly_from_result_returns_the_same_object(self,
                                                          field_result):
        assembly = assembly_from_result(field_result)
        assert assembly is RomModel.from_result(field_result).assembly
        assert assembly.name == "rom_arm"

    def test_from_assembly_needs_a_root_when_none_was_solved(self,
                                                             field_result):
        assembly = assembly_from_result(field_result)
        assert RomModel.from_assembly(assembly, root="arm").root == "arm"
        with pytest.raises(RomError, match="not in the assembly"):
            RomModel.from_assembly(assembly, root="ghost")

    def test_check_joint_values_accepts_the_declared_range(self, model):
        model.check_joint_values({"elbow": math.radians(90.0)})
        model.check_joint_values({})
        with pytest.raises(RomError, match="was given"):
            model.check_joint_values({"elbow": float("nan")})


# ===========================================================================
# Sweeps
# ===========================================================================

class TestSchedules:
    def test_radian_range_includes_the_stop(self):
        values = radian_range(0.0, 90.0, 30.0)
        assert [round(math.degrees(v), 6) for v in values] == [
            0.0, 30.0, 60.0, 90.0
        ]

    def test_radian_range_adds_the_stop_when_the_step_misses_it(self):
        values = radian_range(0.0, 100.0, 30.0)
        assert round(math.degrees(values[-1]), 6) == 100.0

    def test_radian_range_rejects_a_backwards_range(self):
        with pytest.raises(RomError, match="below start_deg"):
            radian_range(10.0, 0.0, 1.0)

    def test_schedule_from_model_spans_the_declared_limits(self, model):
        schedule = schedule_from_model(model, step_deg=45.0)
        assert list(schedule) == ["elbow"]
        degrees = [round(math.degrees(v), 6) for v in schedule["elbow"]]
        assert degrees == [-90.0, -45.0, 0.0, 45.0, 90.0]

    def test_schedule_from_model_refuses_an_unbounded_joint(self, tmp_path):
        source = ARM_DESIGN.replace(
            LIMITS_LINE, "",
        )
        result = build(source, "BUILD", "sdf", tmp_path)
        model = RomModel.from_result(result)
        with pytest.raises(RomError, match="set_mate_limits"):
            schedule_from_model(model)

    def test_schedule_from_model_refuses_a_rigid_joint(self, model):
        with pytest.raises(RomError, match="cannot be driven"):
            schedule_from_model(model, joints=["stop_mount"])


class TestSweep:
    def test_the_arm_clears_the_stop_then_hits_it(self, model):
        sweep = sweep_rom(
            model, {"elbow": radian_range(0.0, 90.0, 30.0)},
            config=SweepConfig(pitch_mm=0.5),
        )
        assert sweep.configurations == 4
        assert sweep.pairs == (("base", "arm"), ("arm", "stop"))

        trace = {round(math.degrees(s.config["elbow"])): s
                 for s in sweep.trace("arm", "stop")}
        assert trace[0].min_distance == pytest.approx(
            10.0, abs=max(trace[0].clearance.tolerance_mm, 1e-9)
        )
        assert not trace[0].clearance.interferes
        assert trace[90].clearance.interferes
        assert trace[90].min_distance == pytest.approx(
            -5.0, abs=trace[90].clearance.tolerance_mm
        )

        worst = sweep.worst()
        assert worst.pair == ("arm", "stop")
        assert round(math.degrees(worst.config["elbow"])) == 90

    def test_a_pair_that_never_moves_is_not_measured_by_default(self, model):
        sweep = sweep_rom(model, {"elbow": [0.0]})
        assert ("base", "stop") not in sweep.pairs

    def test_every_pair_can_be_asked_for(self, model):
        sweep = sweep_rom(model, {"elbow": [0.0]}, pairs="all")
        assert len(sweep.pairs) == 3

    def test_mate_pairs_can_be_asked_for(self, model):
        sweep = sweep_rom(model, {"elbow": [0.0]}, pairs="mates")
        assert set(sweep.pairs) == {("base", "arm"), ("base", "stop")}

    def test_an_unknown_pair_selector_is_an_error(self, model):
        with pytest.raises(RomError, match="pairs must be"):
            sweep_rom(model, {"elbow": [0.0]}, pairs="adjacent")

    def test_extra_pairs_are_added_and_validated(self, model):
        sweep = sweep_rom(model, {"elbow": [0.0]},
                          extra_pairs=[("base", "stop")])
        assert ("base", "stop") in sweep.pairs
        with pytest.raises(RomError, match="No part named 'ghost'"):
            sweep_rom(model, {"elbow": [0.0]},
                      extra_pairs=[("base", "ghost")])

    def test_distant_pairs_are_resolved_from_their_bounds(self, model):
        sweep = sweep_rom(model, {"elbow": [0.0]},
                          config=SweepConfig(clearance_threshold=5.0))
        methods = {s.pair: s.clearance.method for s in sweep.samples}
        # base and arm are 10 mm apart: no sampling needed to know that.
        assert methods[("base", "arm")] == "bounds"
        assert sweep.from_bounds >= 1

    def test_the_smaller_part_of_a_pair_is_the_one_sampled(self, model):
        """Sampling the 40 mm base at a 0.2 mm pitch would need 8 million
        points and be coarsened; the 10 mm stop needs 125 thousand."""
        sweep = sweep_rom(
            model, {"elbow": [math.radians(90.0)]},
            pairs=[("base", "stop")],
            config=SweepConfig(pitch_mm=0.2, clearance_threshold=5.0,
                               max_samples=1_000_000),
        )
        sample = sweep.samples[0]
        assert sample.clearance.method == "field"
        assert sample.clearance.pitch_mm == pytest.approx(0.2)
        assert sample.clearance.samples < 1_000_000

    def test_thresholds_are_views_over_the_samples(self, model):
        sweep = sweep_rom(
            model, {"elbow": radian_range(0.0, 90.0, 30.0)},
            config=SweepConfig(pitch_mm=1.0, clearance_threshold=12.0),
        )
        assert len(sweep.events) == 2
        assert len(sweep.near_misses) + len(sweep.events) == len(
            sweep.samples
        )
        minima = sweep.pair_minima()
        assert set(minima) == {("base", "arm"), ("arm", "stop")}
        assert minima[("arm", "stop")].clearance.interferes

    def test_an_unknown_joint_is_refused_before_any_measurement(self, model):
        with pytest.raises(RomError, match="No joint named 'knee'"):
            sweep_rom(model, {"knee": [0.0]})

    def test_a_value_beyond_the_limits_is_refused(self, model):
        with pytest.raises(RomError, match="above its limit"):
            sweep_rom(model, {"elbow": [math.radians(100.0)]})

    def test_a_joint_cannot_be_both_swept_and_held(self, model):
        with pytest.raises(RomError, match="swept and fixed"):
            sweep_rom(model, {"elbow": [0.0]}, fixed={"elbow": 0.1})

    def test_a_schedule_with_no_joints_is_refused(self, model):
        with pytest.raises(RomError, match="names no joints"):
            sweep_rom(model, {})

    def test_an_empty_schedule_is_refused(self, model):
        with pytest.raises(RomError, match="schedule for 'elbow' is empty"):
            sweep_rom(model, {"elbow": []})

    def test_early_stop_says_the_sweep_is_incomplete(self, model):
        sweep = sweep_rom(
            model, {"elbow": radian_range(0.0, 90.0, 10.0)},
            config=SweepConfig(pitch_mm=1.0, early_stop=1),
        )
        assert sweep.stopped_early
        assert "incomplete" in sweep.summary()
        assert sweep.configurations == 10

    def test_a_sweep_restores_the_pose_it_found(self, model):
        model.pose({"elbow": math.radians(30.0)})
        sweep_rom(model, {"elbow": radian_range(0.0, 90.0, 45.0)})
        assert model.joint_values()["elbow"] == pytest.approx(
            math.radians(30.0)
        )

    def test_the_default_schedule_sweeps_every_movable_joint(self, model):
        sweep = sweep_rom(model, config=SweepConfig(pitch_mm=2.0))
        assert list(sweep.schedule) == ["elbow"]
        assert sweep.configurations == 37  # -90..90 at the 5 deg default


class TestSweepSilenceIsImpossible:
    """A sweep that measured nothing must not read as a clean result."""

    @pytest.fixture
    def mesh_model(self, tmp_path):
        return RomModel.from_result(build(ARM_DESIGN, "BUILD", "mesh",
                                          tmp_path))

    def test_a_pair_without_fields_stops_the_sweep(self, mesh_model):
        with pytest.raises(RomError, match="Cannot measure"):
            sweep_rom(mesh_model, {"elbow": [0.0]})

    def test_skipping_is_allowed_but_always_reported(self, mesh_model):
        sweep = sweep_rom(mesh_model, {"elbow": [0.0]},
                          config=SweepConfig(skip_unmeasurable=True))
        assert sweep.samples == []
        assert len(sweep.skipped) == 2
        summary = sweep.summary()
        assert "NOTHING WAS MEASURED" in summary
        assert "not a clean result" in summary
        assert "pairs NOT measured" in summary
        assert "representation='sdf'" in summary
        report = interference_summary(sweep)
        assert "NOT MEASURED" in report
        assert "says nothing about the" in report

    def test_a_partial_skip_still_measures_the_rest(self, tmp_path,
                                                    field_result):
        """One part without a field must not hide the others."""
        model = RomModel.from_result(field_result)
        model.links["stop"] = type(model.links["stop"])(
            name="stop", field=None, reason="pretend it was imported",
        )
        sweep = sweep_rom(model, {"elbow": radian_range(0.0, 90.0, 45.0)},
                          config=SweepConfig(skip_unmeasurable=True))
        assert sweep.pairs == (("base", "arm"),)
        assert ("arm", "stop") in sweep.skipped
        assert sweep.samples
        assert "pairs NOT measured: 1" in sweep.summary()

    def test_a_sweep_with_no_movable_joint_is_refused(self, tmp_path):
        source = ARM_DESIGN.replace(
            '"elbow", "revolute"', '"elbow", "rigid"'
        ).replace(LIMITS_LINE, "")
        result = build(source, "BUILD", "sdf", tmp_path)
        model = RomModel.from_result(result)
        with pytest.raises(RomError, match="no movable joints"):
            sweep_rom(model)


# ===========================================================================
# Reports
# ===========================================================================

class TestReport:
    def test_the_report_names_the_worst_pair_and_the_tolerance(self, model):
        sweep = sweep_rom(
            model, {"elbow": radian_range(0.0, 90.0, 45.0)},
            config=SweepConfig(pitch_mm=1.0),
        )
        report = interference_summary(sweep)
        assert "rom_arm" in report
        assert "arm - stop" in report
        assert "Interference" in report
        assert "sampling pitch" in report
        assert "elbow=90.0 deg" in report

    def test_a_clean_sweep_says_so_without_crashing(self, model):
        sweep = sweep_rom(
            model, {"elbow": radian_range(0.0, 30.0, 30.0)},
            config=SweepConfig(pitch_mm=1.0),
        )
        report = interference_summary(sweep)
        assert "No interference found." in report
        assert "Worst configuration per pair" in report

    def test_the_report_survives_a_hand_built_result(self, model):
        """Formatting used to crash on a result it had not measured."""
        from yapcad.rom.sweep import PairSample, SweepResult

        clearance = Clearance(
            min_distance=-2.5, penetration=2.5, mean_distance=4.0,
            fraction_inside=0.1, samples=10, pitch_mm=1.0,
            tolerance_mm=0.9, method="field",
        )
        result = SweepResult(
            model=model, config=SweepConfig(), configurations=1,
            pairs=(("arm", "stop"),),
            samples=[PairSample(pair=("arm", "stop"),
                                config={"elbow": 0.5},
                                clearance=clearance)],
        )
        report = interference_summary(result)
        assert "-2.500" in report
        assert "penetration 2.500" in report
        assert result.summary().count("arm - stop") >= 1
