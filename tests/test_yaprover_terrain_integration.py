"""End-to-end acceptance of oracle-driven YapRover assembly poses."""

from pathlib import Path

import numpy as np
import pytest

from yapcad.dsl import compile_and_run
from yapcad.kinematics.rocker_bogie import apply_terrain_pose

from test_yaprover_suspension import _make_assembly


EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "yaprover_suspension.dsl"
SHOWCASE_TERRAINS = (
    {},
    {"lf": 80.0},
    {"lf": 80.0, "rf": 80.0},
    {"lm": 80.0},
    {"lr": 80.0},
    {"lf": 80.0, "rr": 80.0},
)
WHEEL_DATUMS = {
    "lf": ("left_rocker", "front_axle"),
    "lm": ("left_bogie", "middle_axle"),
    "lr": ("left_bogie", "rear_axle"),
    "rf": ("right_rocker", "front_axle"),
    "rm": ("right_bogie", "middle_axle"),
    "rr": ("right_bogie", "rear_axle"),
}


@pytest.fixture(scope="module")
def rover_assembly():
    source = EXAMPLE.read_text(encoding="utf-8")
    solids = {}
    for command in ("CHASSIS", "ROCKER", "BOGIE", "WHEEL"):
        result = compile_and_run(source, command, {})
        assert result.success, result.error_message
        solids[command] = result.geometry
    return _make_assembly(solids).data


def _contact_residuals(assembly, terrain):
    return {
        wheel: assembly.get_transformed_datum(*datum).origin[2]
        - (65.0 + terrain.get(wheel, 0.0))
        for wheel, datum in WHEEL_DATUMS.items()
    }


@pytest.mark.parametrize("terrain", SHOWCASE_TERRAINS)
def test_showcase_terrain_pose_matches_all_six_contacts(rover_assembly, terrain):
    pose, result = apply_terrain_pose(rover_assembly, terrain)

    assert pose.success, pose.errors
    assert result.success, result.errors
    assert max(abs(value) for value in _contact_residuals(
        rover_assembly, terrain,
    ).values()) <= 0.10
    assert result.coupling_residuals["rocker_differential"] <= np.deg2rad(0.10)
    assert max(result.residuals.values()) <= 0.05


@pytest.mark.parametrize("target", SHOWCASE_TERRAINS[1:])
def test_every_showcase_terrain_has_a_continuous_21_step_path(
        rover_assembly, target):
    previous = None
    for fraction in np.linspace(0.0, 1.0, 21):
        terrain = {wheel: height * fraction for wheel, height in target.items()}
        pose, result = apply_terrain_pose(rover_assembly, terrain)
        assert pose.success, (fraction, pose.errors)
        assert result.success, (fraction, result.errors)
        assert max(abs(value) for value in _contact_residuals(
            rover_assembly, terrain,
        ).values()) <= 0.10
        coordinates = np.array([
            pose.chassis_height, pose.chassis_pitch, pose.chassis_roll,
            pose.left_rocker, pose.left_bogie, pose.right_bogie,
        ])
        if previous is not None:
            # A single interpolation interval cannot jump to another branch.
            assert abs(coordinates[0] - previous[0]) < 10.0
            assert np.max(np.abs(coordinates[1:] - previous[1:])) < 0.10
        previous = coordinates


def test_invalid_middle_wheel_pose_does_not_replace_last_valid_assembly_pose(
        rover_assembly):
    pose, result = apply_terrain_pose(rover_assembly, {"lm": 80.0})
    assert pose.success and result.success
    before = {
        name: transform.copy()
        for name, transform in rover_assembly.transforms.items()
    }

    invalid_pose, invalid_result = apply_terrain_pose(
        rover_assembly, {"lm": 100.0},
    )
    assert not invalid_pose.success
    assert invalid_result is None
    for name, transform in before.items():
        np.testing.assert_allclose(rover_assembly.transforms[name], transform)


def test_coupling_mismatch_restores_complete_previous_assembly_state(
        rover_assembly):
    pose, result = apply_terrain_pose(rover_assembly, {"lf": 80.0})
    assert pose.success and result.success
    before_transforms = {
        name: transform.copy()
        for name, transform in rover_assembly.transforms.items()
    }
    before_joint_values = rover_assembly._joint_values.copy()
    before_prescribed = rover_assembly._prescribed_joint_values.copy()
    before_root = rover_assembly._root_part
    before_solved = rover_assembly._solved

    with pytest.raises(ValueError, match="coupling"):
        apply_terrain_pose(
            rover_assembly,
            {"lm": 80.0},
            right_rocker_joint="left_rocker_pivot",
        )

    assert rover_assembly._joint_values == before_joint_values
    assert rover_assembly._prescribed_joint_values == before_prescribed
    assert rover_assembly._root_part == before_root
    assert rover_assembly._solved is before_solved
    for name, transform in before_transforms.items():
        np.testing.assert_allclose(rover_assembly.transforms[name], transform)
