"""Contract tests for rooted assembly graph solving and articulation.

These tests intentionally exercise :meth:`Assembly.solve` rather than the
older, ordered ``solve_mate_chain`` helper.  An assembly solve is a graph
operation: mate declaration order must not matter, each non-root part has one
placement parent, and the returned world transforms must form an independent
snapshot of the requested pose.
"""

import math

import numpy as np
import pytest

from yapcad.assembly import (
    Assembly,
    Datum,
    DatumType,
    Mate,
    MateLimits,
    MateType,
    PartDefinition,
)


IDENTITY = np.eye(4)


def _axis(name, origin=(0.0, 0.0, 0.0), direction=(0.0, 0.0, 1.0)):
    return Datum(name, DatumType.AXIS, origin=origin, direction=direction)


def _frame(name, origin=(0.0, 0.0, 0.0)):
    return Datum(
        name,
        DatumType.FRAME,
        origin=origin,
        x_axis=(1.0, 0.0, 0.0),
        y_axis=(0.0, 1.0, 0.0),
    )


def _part(name, *datums):
    part = PartDefinition(name)
    for datum in datums:
        part.add_datum(datum)
    return part


def _mate(name, mate_type, parent, parent_datum, child, child_datum,
          limits=None):
    return Mate(
        name=name,
        mate_type=mate_type,
        part_a=parent,
        datum_a=parent_datum,
        part_b=child,
        datum_b=child_datum,
        limits=limits,
    )


def _errors(result):
    """Flatten diagnostics without prescribing list versus mapping storage."""
    errors = result.errors
    if isinstance(errors, dict):
        errors = errors.values()
    return " ".join(str(error) for error in errors).lower()


def _assert_rigid(transform):
    assert transform.shape == (4, 4)
    np.testing.assert_allclose(transform[3], [0.0, 0.0, 0.0, 1.0])
    rotation = transform[:3, :3]
    np.testing.assert_allclose(rotation.T @ rotation, np.eye(3), atol=1e-10)
    assert np.linalg.det(rotation) == pytest.approx(1.0, abs=1e-10)


def _simple_rigid_assembly():
    assembly = Assembly("rigid")
    assembly.add_part(_part("base", _frame("mount", (10.0, 0.0, 0.0))))
    assembly.add_part(_part("body", _frame("plug", (2.0, 0.0, 0.0))))
    assembly.add_mate(
        _mate("base_body", MateType.RIGID, "base", "mount", "body", "plug")
    )
    return assembly


def _branching_assembly(mate_order=("left", "right")):
    assembly = Assembly("branch")
    assembly.add_part(
        _part(
            "base",
            _axis("left", (10.0, 0.0, 0.0)),
            _axis("right", (0.0, 20.0, 0.0)),
        )
    )
    assembly.add_part(_part("left_link", _axis("pivot")))
    assembly.add_part(_part("right_link", _axis("pivot")))
    mates = {
        "left": _mate(
            "left_joint", MateType.REVOLUTE,
            "base", "left", "left_link", "pivot",
        ),
        "right": _mate(
            "right_joint", MateType.REVOLUTE,
            "base", "right", "right_link", "pivot",
        ),
    }
    for key in mate_order:
        assembly.add_mate(mates[key])
    return assembly


def test_solve_rigid_mate_positions_child_and_updates_assembly():
    assembly = _simple_rigid_assembly()

    result = assembly.solve(root_part="base")

    assert result.success
    np.testing.assert_allclose(result.transforms["base"], IDENTITY)
    np.testing.assert_allclose(result.transforms["body"][:3, 3], [8.0, 0.0, 0.0])
    np.testing.assert_allclose(assembly.transforms["body"], result.transforms["body"])


def test_solve_uses_root_instance_transform_as_world_anchor():
    assembly = _simple_rigid_assembly()
    root_transform = np.eye(4)
    root_transform[:3, :3] = [[0.0, -1.0, 0.0],
                              [1.0, 0.0, 0.0],
                              [0.0, 0.0, 1.0]]
    root_transform[:3, 3] = [100.0, 200.0, 5.0]
    assembly.transforms["base"] = root_transform

    result = assembly.solve(root_part="base")

    assert result.success
    np.testing.assert_allclose(result.transforms["base"], root_transform)
    np.testing.assert_allclose(result.transforms["body"][:3, 3], [100.0, 208.0, 5.0])


def test_branching_solve_is_independent_of_mate_input_order():
    forward = _branching_assembly(("left", "right")).solve(
        root_part="base",
        joint_values={"left_joint": math.pi / 3, "right_joint": -math.pi / 4},
    )
    reversed_order = _branching_assembly(("right", "left")).solve(
        root_part="base",
        joint_values={"left_joint": math.pi / 3, "right_joint": -math.pi / 4},
    )

    assert forward.success and reversed_order.success
    assert forward.transforms.keys() == reversed_order.transforms.keys()
    for part_name in forward.transforms:
        np.testing.assert_allclose(
            forward.transforms[part_name], reversed_order.transforms[part_name]
        )


def test_branch_children_receive_their_own_parent_datum_transforms():
    result = _branching_assembly().solve(root_part="base")

    assert result.success
    np.testing.assert_allclose(result.transforms["left_link"][:3, 3], [10.0, 0.0, 0.0])
    np.testing.assert_allclose(result.transforms["right_link"][:3, 3], [0.0, 20.0, 0.0])


def test_revolute_joint_value_rotates_child_about_mated_axis():
    result = _branching_assembly().solve(
        root_part="base", joint_values={"left_joint": math.pi / 2}
    )

    assert result.success
    expected_rotation = np.array([[0.0, -1.0, 0.0],
                                  [1.0, 0.0, 0.0],
                                  [0.0, 0.0, 1.0]])
    np.testing.assert_allclose(
        result.transforms["left_link"][:3, :3], expected_rotation, atol=1e-12
    )
    np.testing.assert_allclose(result.transforms["left_link"][:3, 3], [10.0, 0.0, 0.0])


def test_revolute_rotation_propagates_to_descendant_branch():
    assembly = Assembly("chain")
    assembly.add_part(_part("base", _axis("pivot")))
    assembly.add_part(
        _part("rocker", _axis("root"), _frame("axle", (10.0, 0.0, 0.0)))
    )
    assembly.add_part(_part("wheel", _frame("mount")))
    assembly.add_mate(
        _mate("rocker_joint", MateType.REVOLUTE, "base", "pivot", "rocker", "root")
    )
    assembly.add_mate(
        _mate("wheel_mount", MateType.RIGID, "rocker", "axle", "wheel", "mount")
    )

    result = assembly.solve(
        root_part="base", joint_values={"rocker_joint": math.pi / 2}
    )

    assert result.success
    np.testing.assert_allclose(result.transforms["wheel"][:3, 3], [0.0, 10.0, 0.0], atol=1e-12)
    np.testing.assert_allclose(
        result.transforms["wheel"][:3, :3], result.transforms["rocker"][:3, :3]
    )


@pytest.mark.parametrize("angle", [-math.pi / 4, math.pi / 4])
def test_revolute_joint_accepts_inclusive_limit_boundaries(angle):
    assembly = _branching_assembly()
    assembly.mates[0].limits = MateLimits(
        min_value=-math.pi / 4, max_value=math.pi / 4
    )

    result = assembly.solve(root_part="base", joint_values={"left_joint": angle})

    assert result.success


def test_revolute_joint_rejects_value_outside_limits_without_mutation():
    assembly = _branching_assembly()
    assembly.mates[0].limits = MateLimits(min_value=-0.5, max_value=0.5)
    before = {name: transform.copy() for name, transform in assembly.transforms.items()}

    result = assembly.solve(root_part="base", joint_values={"left_joint": 0.5001})

    assert not result.success
    assert "left_joint" in _errors(result)
    assert "limit" in _errors(result)
    for part_name, transform in before.items():
        np.testing.assert_allclose(assembly.transforms[part_name], transform)


def test_successful_solve_reports_small_residuals_and_rigid_transforms():
    result = _branching_assembly().solve(
        root_part="base",
        joint_values={"left_joint": 0.37, "right_joint": -0.22},
    )

    assert result.success
    assert set(result.residuals) == {"left_joint", "right_joint"}
    assert all(residual <= 1e-9 for residual in result.residuals.values())
    for transform in result.transforms.values():
        _assert_rigid(transform)


def test_repeated_solve_recomputes_pose_and_returns_independent_snapshots():
    assembly = _branching_assembly()
    first = assembly.solve(
        root_part="base", joint_values={"left_joint": 0.0}
    )
    first_snapshot = first.transforms["left_link"].copy()

    second = assembly.solve(
        root_part="base", joint_values={"left_joint": math.pi / 2}
    )

    assert first.success and second.success
    assert not np.allclose(first.transforms["left_link"], second.transforms["left_link"])
    np.testing.assert_allclose(first.transforms["left_link"], first_snapshot)
    assert first.transforms["left_link"] is not second.transforms["left_link"]
    assert second.transforms["left_link"] is not assembly.transforms["left_link"]


def test_missing_root_returns_diagnostic_without_mutating_transforms():
    assembly = _simple_rigid_assembly()
    before = {name: transform.copy() for name, transform in assembly.transforms.items()}

    result = assembly.solve(root_part="not_a_part")

    assert not result.success
    assert "root" in _errors(result)
    assert "not_a_part" in _errors(result)
    for part_name, transform in before.items():
        np.testing.assert_allclose(assembly.transforms[part_name], transform)


def test_disconnected_part_returns_diagnostic_and_preserves_last_good_pose():
    assembly = _simple_rigid_assembly()
    good = assembly.solve(root_part="base")
    assert good.success
    last_good = {name: transform.copy() for name, transform in assembly.transforms.items()}
    assembly.add_part(_part("orphan", _frame("mount")))
    last_good["orphan"] = assembly.transforms["orphan"].copy()

    result = assembly.solve(root_part="base")

    assert not result.success
    assert "orphan" in _errors(result)
    assert "disconnect" in _errors(result) or "unreachable" in _errors(result)
    for part_name, transform in last_good.items():
        np.testing.assert_allclose(assembly.transforms[part_name], transform)


def test_cycle_returns_explicit_graph_diagnostic():
    assembly = Assembly("cycle")
    assembly.add_part(_part("a", _frame("mount")))
    assembly.add_part(_part("b", _frame("mount")))
    assembly.add_mate(_mate("a_b", MateType.RIGID, "a", "mount", "b", "mount"))
    assembly.add_mate(_mate("b_a", MateType.RIGID, "b", "mount", "a", "mount"))

    result = assembly.solve(root_part="a")

    assert not result.success
    assert "cycle" in _errors(result)


def test_duplicate_placement_parent_returns_child_and_mates_in_diagnostic():
    assembly = Assembly("duplicate-parent")
    assembly.add_part(_part("root", _frame("a"), _frame("b"), _frame("c")))
    assembly.add_part(_part("parent_2", _frame("root"), _frame("child")))
    assembly.add_part(_part("child", _frame("root")))
    assembly.add_mate(
        _mate("attach_parent_2", MateType.RIGID, "root", "a", "parent_2", "root")
    )
    assembly.add_mate(
        _mate("direct_child", MateType.RIGID, "root", "b", "child", "root")
    )
    assembly.add_mate(
        _mate("indirect_child", MateType.RIGID, "parent_2", "child", "child", "root")
    )

    result = assembly.solve(root_part="root")

    assert not result.success
    diagnostic = _errors(result)
    assert "child" in diagnostic
    assert "parent" in diagnostic
    assert "direct_child" in diagnostic
    assert "indirect_child" in diagnostic


def test_unsupported_placement_mate_returns_actionable_diagnostic():
    assembly = Assembly("unsupported")
    assembly.add_part(_part("base", _axis("rail")))
    assembly.add_part(_part("carriage", _axis("rail")))
    assembly.add_mate(
        _mate(
            "linear_slide", MateType.PRISMATIC,
            "base", "rail", "carriage", "rail",
        )
    )

    result = assembly.solve(root_part="base")

    assert not result.success
    diagnostic = _errors(result)
    assert "linear_slide" in diagnostic
    assert "prismatic" in diagnostic
    assert "unsupported" in diagnostic


def test_unknown_joint_value_name_is_rejected_as_likely_typo():
    assembly = _branching_assembly()

    result = assembly.solve(
        root_part="base", joint_values={"left_jiont": 0.25}
    )

    assert not result.success
    diagnostic = _errors(result)
    assert "left_jiont" in diagnostic
    assert "unknown" in diagnostic or "not found" in diagnostic
