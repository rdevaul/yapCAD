"""Contract tests for deterministic affine coupling of assembly joints."""

import math

import numpy as np
import pytest

from yapcad.assembly import (
    Assembly,
    AssemblyError,
    Datum,
    DatumType,
    LinearJointCoupling,
    Mate,
    MateLimits,
    MateType,
    PartDefinition,
)
from yapcad.assembly.mechatron_export import to_mechatron_snapshot


def _axis(name):
    return Datum(name, DatumType.AXIS, origin=(0, 0, 0), direction=(0, 0, 1))


def _part(name):
    part = PartDefinition(name)
    part.add_datum(_axis("pivot"))
    return part


def _joint(name, parent, child, limits=None):
    return Mate(
        name=name,
        mate_type=MateType.REVOLUTE,
        part_a=parent,
        datum_a="pivot",
        part_b=child,
        datum_b="pivot",
        limits=limits,
    )


def _differential():
    assembly = Assembly("differential")
    for name in ("chassis", "left_rocker", "right_rocker"):
        assembly.add_part(_part(name))
    assembly.add_mate(_joint("left", "chassis", "left_rocker"))
    assembly.add_mate(_joint("right", "chassis", "right_rocker"))
    assembly.add_joint_coupling(LinearJointCoupling(
        name="rocker_differential",
        dependent_joint="right",
        driver_coefficients={"left": -1.0},
    ))
    return assembly


def _rotation_angle(transform):
    return math.atan2(transform[1, 0], transform[0, 0])


def _errors(result):
    return " ".join(result.errors).lower()


def test_rover_differential_derives_equal_and_opposite_rocker_angles():
    assembly = _differential()
    result = assembly.solve("chassis", {"left": math.radians(12)})

    assert result.success, result.errors
    assert result.joint_values["left"] == pytest.approx(math.radians(12))
    assert result.joint_values["right"] == pytest.approx(math.radians(-12))
    assert _rotation_angle(result.transforms["left_rocker"]) == pytest.approx(
        math.radians(12)
    )
    assert _rotation_angle(result.transforms["right_rocker"]) == pytest.approx(
        math.radians(-12)
    )
    assert result.coupling_residuals["rocker_differential"] < 1e-12


def test_zero_pose_is_well_defined_when_no_driver_is_prescribed():
    result = _differential().solve("chassis")
    assert result.success
    assert result.joint_values == {"left": 0.0, "right": 0.0}


def test_affine_offset_and_multiple_drivers_are_supported():
    assembly = Assembly("multi")
    for name in ("base", "a", "b", "carrier"):
        assembly.add_part(_part(name))
    assembly.add_mate(_joint("left", "base", "a"))
    assembly.add_mate(_joint("right", "base", "b"))
    assembly.add_mate(_joint("average", "base", "carrier"))
    assembly.add_joint_coupling(LinearJointCoupling(
        name="carrier_average",
        dependent_joint="average",
        driver_coefficients={"left": 0.5, "right": 0.5},
        offset=0.1,
    ))

    result = assembly.solve("base", {"left": 0.4, "right": -0.2})
    assert result.success
    assert result.joint_values["average"] == pytest.approx(0.2)


def test_chained_couplings_resolve_independent_of_declaration_order():
    assembly = Assembly("chain")
    for name in ("base", "a", "b", "c"):
        assembly.add_part(_part(name))
    for joint, child in (("a", "a"), ("b", "b"), ("c", "c")):
        assembly.add_mate(_joint(joint, "base", child))
    assembly.add_joint_coupling(LinearJointCoupling(
        "second", "c", {"b": 3.0}, offset=0.25,
    ))
    assembly.add_joint_coupling(LinearJointCoupling(
        "first", "b", {"a": 2.0}, offset=-0.5,
    ))

    result = assembly.solve("base", {"a": 0.4})
    assert result.success
    assert result.joint_values["b"] == pytest.approx(0.3)
    assert result.joint_values["c"] == pytest.approx(1.15)


def test_prescribing_a_derived_joint_is_rejected_as_ambiguous():
    result = _differential().solve("chassis", {"right": 0.2})
    assert not result.success
    assert "dependent" in _errors(result) and "right" in _errors(result)


def test_set_joint_position_reposes_driver_and_all_dependents():
    assembly = _differential()
    assert assembly.solve("chassis").success
    result = assembly.set_joint_position("left", 0.35)
    assert result.joint_values["right"] == pytest.approx(-0.35)
    assert assembly._joint_values["right"] == pytest.approx(-0.35)


def test_set_joint_position_rejects_a_dependent_joint():
    assembly = _differential()
    assert assembly.solve("chassis").success
    with pytest.raises(AssemblyError, match="dependent.*right|right.*dependent"):
        assembly.set_joint_position("right", 0.2)


def test_derived_value_is_checked_against_joint_limits_transactionally():
    assembly = _differential()
    right = next(m for m in assembly.mates if m.name == "right")
    right.limits = MateLimits(min_value=-0.25, max_value=0.25)
    assert assembly.solve("chassis", {"left": 0.1}).success
    before = {name: value.copy() for name, value in assembly.transforms.items()}

    result = assembly.solve("chassis", {"left": 0.4})
    assert not result.success
    assert "limit" in _errors(result) and "right" in _errors(result)
    for name in before:
        np.testing.assert_allclose(assembly.transforms[name], before[name])


@pytest.mark.parametrize(
    "coupling, message",
    [
        (LinearJointCoupling("bad", "ghost", {"left": 1}), "ghost"),
        (LinearJointCoupling("bad", "right", {"ghost": 1}), "ghost"),
    ],
)
def test_invalid_coupling_references_are_diagnostic(coupling, message):
    assembly = _differential()
    assembly.joint_couplings.clear()
    with pytest.raises(AssemblyError, match=message):
        assembly.add_joint_coupling(coupling)


@pytest.mark.parametrize(
    "factory, message",
    [
        (lambda: LinearJointCoupling("bad", "right", {"right": 1}), "itself|self"),
        (lambda: LinearJointCoupling("bad", "right", {}), "driver"),
    ],
)
def test_structurally_invalid_couplings_fail_at_construction(factory, message):
    with pytest.raises(ValueError, match=message):
        factory()


def test_only_revolute_joint_coordinates_are_supported_initially():
    assembly = _differential()
    assembly.joint_couplings.clear()
    next(m for m in assembly.mates if m.name == "right").mate_type = MateType.RIGID
    with pytest.raises(AssemblyError, match="revolute"):
        assembly.add_joint_coupling(LinearJointCoupling("bad", "right", {"left": 1}))


def test_duplicate_dependent_joint_is_rejected():
    assembly = _differential()
    with pytest.raises(AssemblyError, match="dependent|already"):
        assembly.add_joint_coupling(LinearJointCoupling(
            "second_differential", "right", {"left": 2},
        ))


def test_coupling_dependency_cycle_is_reported_without_recursion_failure():
    assembly = _differential()
    assembly.joint_couplings.clear()
    assembly.add_joint_coupling(LinearJointCoupling("left_from_right", "left", {"right": -1}))
    assembly.add_joint_coupling(LinearJointCoupling("right_from_left", "right", {"left": -1}))

    result = assembly.solve("chassis")
    assert not result.success
    assert "cycle" in _errors(result)


def test_coupling_contract_is_preserved_in_semantic_graph_export():
    coupling = to_mechatron_snapshot(_differential())["joint_couplings"][0]
    assert coupling == {
        "id": "rocker_differential",
        "type": "Affine",
        "dependent_joint": "right",
        "driver_coefficients": {"left": -1.0},
        "offset": 0.0,
        "tolerance": 1e-9,
    }
