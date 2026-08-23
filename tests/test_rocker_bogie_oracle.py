"""Regression tests for the independent YapRover terrain/contact oracle."""

import math

import pytest

from yapcad.kinematics.rocker_bogie import (
    RockerBogieGeometry,
    solve_terrain_pose,
    wheel_centers,
)


def _degrees(value):
    return math.degrees(value)


def _assert_pose(pose, height, pitch, roll, left_rocker,
                 left_bogie=None, right_bogie=None):
    assert pose.success, pose.errors
    assert pose.chassis_height == pytest.approx(height, abs=0.25)
    assert _degrees(pose.chassis_pitch) == pytest.approx(pitch, abs=0.25)
    assert _degrees(pose.chassis_roll) == pytest.approx(roll, abs=0.25)
    assert _degrees(pose.left_rocker) == pytest.approx(left_rocker, abs=0.25)
    if left_bogie is not None:
        assert _degrees(pose.left_bogie) == pytest.approx(left_bogie, abs=0.25)
    if right_bogie is not None:
        assert _degrees(pose.right_bogie) == pytest.approx(right_bogie, abs=0.25)
    assert pose.max_contact_residual <= 1e-7
    assert pose.closure_residual == pytest.approx(0.0, abs=1e-12)
    assert pose.left_rocker + pose.right_rocker == pytest.approx(0.0, abs=1e-12)


def test_flat_pose_matches_nominal_geometry():
    pose = solve_terrain_pose({})

    _assert_pose(pose, 140.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    assert all(center[1] == pytest.approx(65.0)
               for center in pose.wheel_centers.values())


@pytest.mark.parametrize(
    "terrain, expected",
    [
        ({"lf": 80.0}, (152.4, 7.99, 4.59, 7.99, -15.98, 0.0)),
        ({"lf": 80.0, "rf": 80.0},
         (164.8, 15.98, 0.0, 0.0, -15.98, -15.98)),
        ({"lm": 80.0}, (151.9, -3.69, 4.40, -3.69, 32.28, 0.0)),
        ({"lr": 80.0}, (151.9, -3.69, 4.40, -3.69, -17.52, 0.0)),
        ({"lf": 80.0, "rr": 80.0},
         (164.3, 4.30, 0.19, 11.68, -15.98, -17.52)),
    ],
)
def test_eighty_mm_terrain_regression_poses(terrain, expected):
    pose = solve_terrain_pose(terrain)

    _assert_pose(pose, *expected)


def test_wheel_centers_uses_parent_relative_joint_angles():
    centers = wheel_centers(
        chassis_height=140.0,
        chassis_pitch=math.radians(10.0),
        chassis_roll=0.0,
        left_rocker=math.radians(-10.0),
        left_bogie=0.0,
        right_rocker=math.radians(-10.0),
        right_bogie=0.0,
    )

    # Absolute link angles are zero, so nominal X/Z coordinates are retained.
    assert centers["lf"] == pytest.approx((190.0, 65.0))
    assert centers["lm"] == pytest.approx((0.0, 65.0))
    assert centers["lr"] == pytest.approx((-190.0, 65.0))


def test_middle_wheel_at_100_mm_is_rejected_by_bogie_limit():
    pose = solve_terrain_pose({"lm": 100.0})

    assert not pose.success
    assert pose.max_contact_residual <= 1e-7
    assert _degrees(pose.left_bogie) == pytest.approx(40.7, abs=0.25)
    assert any("left bogie" in error and "limit" in error
               for error in pose.errors)


def test_unknown_wheel_name_is_rejected():
    pose = solve_terrain_pose({"left_front": 80.0})

    assert not pose.success
    assert any("unknown wheel" in error.lower() for error in pose.errors)


@pytest.mark.parametrize("height", [math.nan, math.inf, -math.inf])
def test_non_finite_terrain_height_is_rejected(height):
    pose = solve_terrain_pose({"lf": height})

    assert not pose.success
    assert any("finite" in error.lower() for error in pose.errors)


def test_geometry_defaults_are_the_approved_baseline():
    geometry = RockerBogieGeometry()

    assert geometry.center_to_front == (190.0, -75.0)
    assert geometry.center_to_bogie == (-95.0, -35.0)
    assert geometry.bogie_to_middle == (95.0, -40.0)
    assert geometry.bogie_to_rear == (-95.0, -40.0)
    assert geometry.wheel_radius == 65.0
    assert geometry.half_track == 155.0
