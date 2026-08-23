"""Integration contract for the annotated YapRover suspension mini-chain."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from yapcad.brep import has_brep_data, occ_available
from yapcad.dsl import compile_and_run
from yapcad.dsl.runtime.builtins import call_builtin
from yapcad.dsl.runtime.values import (
    float_val,
    list_val,
    solid_val,
    string_val,
)
from yapcad.dsl.types import FLOAT, STRING
from yapcad.geom3d import issolid, solidbbox, volumeof
from yapcad.io.step import write_step_analytic


EXAMPLE = (
    Path(__file__).resolve().parents[1]
    / "examples"
    / "yaprover_suspension.dsl"
)

PART_INSTANCES = (
    ("CHASSIS", "chassis"),
    ("ROCKER", "left_rocker"),
    ("ROCKER", "right_rocker"),
    ("BOGIE", "left_bogie"),
    ("BOGIE", "right_bogie"),
    ("WHEEL", "left_front_wheel"),
    ("WHEEL", "left_middle_wheel"),
    ("WHEEL", "left_rear_wheel"),
    ("WHEEL", "right_front_wheel"),
    ("WHEEL", "right_middle_wheel"),
    ("WHEEL", "right_rear_wheel"),
)

MATES = (
    ("left_rocker_pivot", "chassis", "left_rocker",
     "left_rocker", "chassis_pivot"),
    ("right_rocker_pivot", "chassis", "right_rocker",
     "right_rocker", "chassis_pivot"),
    ("left_bogie_pivot", "left_rocker", "left_bogie",
     "bogie_pivot", "rocker_pivot"),
    ("right_bogie_pivot", "right_rocker", "right_bogie",
     "bogie_pivot", "rocker_pivot"),
    ("left_front_axle", "left_rocker", "left_front_wheel",
     "front_axle", "axle"),
    ("left_middle_axle", "left_bogie", "left_middle_wheel",
     "middle_axle", "axle"),
    ("left_rear_axle", "left_bogie", "left_rear_wheel",
     "rear_axle", "axle"),
    ("right_front_axle", "right_rocker", "right_front_wheel",
     "front_axle", "axle"),
    ("right_middle_axle", "right_bogie", "right_middle_wheel",
     "middle_axle", "axle"),
    ("right_rear_axle", "right_bogie", "right_rear_wheel",
     "rear_axle", "axle"),
)


@pytest.fixture(scope="module")
def source():
    return EXAMPLE.read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def proxy_solids(source):
    solids = {}
    for command in ("CHASSIS", "ROCKER", "BOGIE", "WHEEL"):
        result = compile_and_run(source, command, {})
        assert result.success, result.error_message
        assert issolid(result.geometry)
        solids[command] = result.geometry
    return solids


def _add_named_mate(asm, spec):
    name, parent, child, parent_datum, child_datum = spec
    call_builtin("add_named_mate", [
        asm,
        string_val(name),
        string_val("revolute"),
        string_val(parent),
        string_val(parent_datum),
        string_val(child),
        string_val(child_datum),
    ])


def _make_assembly(proxy_solids):
    asm = call_builtin("assembly", [string_val("yaprover_suspension")])
    for command, instance in PART_INSTANCES:
        call_builtin("add_part", [
            asm, solid_val(proxy_solids[command]), string_val(instance),
        ])
    for mate in MATES:
        _add_named_mate(asm, mate)
    call_builtin("add_joint_coupling", [
        asm,
        string_val("rocker_differential"),
        string_val("right_rocker_pivot"),
        list_val([string_val("left_rocker_pivot")], STRING),
        list_val([float_val(-1.0)], FLOAT),
        float_val(0.0),
    ])
    for name, minimum, maximum in (
        ("left_rocker_pivot", -18.0, 18.0),
        ("right_rocker_pivot", -18.0, 18.0),
        ("left_bogie_pivot", -35.0, 38.0),
        ("right_bogie_pivot", -35.0, 38.0),
    ):
        call_builtin("set_mate_limits", [
            asm, string_val(name),
            float_val(math.radians(minimum)),
            float_val(math.radians(maximum)),
        ])
    return asm


def _extent(solid):
    bounds = solidbbox(solid)
    return np.asarray(bounds[1][:3]) - np.asarray(bounds[0][:3])


def test_proxy_parts_are_annotated_solids_within_fdm_build_volume(proxy_solids):
    """Every unique proxy fits a 220 x 220 x 250 mm printer."""
    for name, solid in proxy_solids.items():
        assert issolid(solid), name
        assert np.all(_extent(solid) <= [220.0, 220.0, 250.0]), name


def test_full_rover_part_set_and_positive_y_axis_annotations(proxy_solids):
    asm = _make_assembly(proxy_solids)

    assert list(asm.data.parts) == [instance for _, instance in PART_INSTANCES]
    assert len(asm.data.parts) == 11
    for part in asm.data.parts.values():
        assert part.datums
        for datum in part.datums.values():
            np.testing.assert_allclose(datum.direction[:3], [0.0, 1.0, 0.0])

    chassis = asm.data.parts["chassis"]
    assert chassis.datums["left_rocker"].origin[:3] == [0.0, 155.0, 0.0]
    assert chassis.datums["right_rocker"].origin[:3] == [0.0, -155.0, 0.0]
    rocker = asm.data.parts["left_rocker"]
    assert rocker.datums["front_axle"].origin[:3] == [190.0, 0.0, -75.0]
    assert rocker.datums["bogie_pivot"].origin[:3] == [-95.0, 0.0, -35.0]
    bogie = asm.data.parts["left_bogie"]
    assert bogie.datums["middle_axle"].origin[:3] == [95.0, 0.0, -40.0]
    assert bogie.datums["rear_axle"].origin[:3] == [-95.0, 0.0, -40.0]


def test_full_mate_tree_is_connected_and_acyclic(proxy_solids):
    asm = _make_assembly(proxy_solids)

    assert len(asm.data.mates) == len(asm.data.parts) - 1 == 10
    assert {mate.name for mate in asm.data.mates} == {
        spec[0] for spec in MATES
    }
    result = call_builtin("solve_assembly", [asm, string_val("chassis")])
    assert result is asm
    assert set(asm.data.transforms) == set(asm.data.parts)
    for mate in asm.data.mates:
        parent = asm.data.get_transformed_datum(mate.part_a, mate.datum_a)
        child = asm.data.get_transformed_datum(mate.part_b, mate.datum_b)
        np.testing.assert_allclose(parent.origin, child.origin, atol=1e-8)


def test_affine_differential_poses_rockers_equally_and_oppositely(proxy_solids):
    asm = _make_assembly(proxy_solids)
    call_builtin("solve_assembly", [asm, string_val("chassis")])
    call_builtin("set_joint_position", [
        asm, string_val("left_rocker_pivot"), float_val(math.radians(18.0)),
    ])

    assert asm.data._joint_values["left_rocker_pivot"] == pytest.approx(
        math.radians(18.0)
    )
    assert asm.data._joint_values["right_rocker_pivot"] == pytest.approx(
        math.radians(-18.0)
    )
    left_fixed = asm.data.get_transformed_datum("chassis", "left_rocker")
    left_moving = asm.data.get_transformed_datum(
        "left_rocker", "chassis_pivot"
    )
    right_fixed = asm.data.get_transformed_datum("chassis", "right_rocker")
    right_moving = asm.data.get_transformed_datum(
        "right_rocker", "chassis_pivot"
    )
    np.testing.assert_allclose(left_fixed.origin, left_moving.origin, atol=1e-8)
    np.testing.assert_allclose(right_fixed.origin, right_moving.origin, atol=1e-8)


def test_all_six_wheel_axes_remain_coincident_after_articulation(proxy_solids):
    asm = _make_assembly(proxy_solids)
    call_builtin("solve_assembly", [asm, string_val("chassis")])
    call_builtin("set_joint_position", [
        asm, string_val("left_rocker_pivot"), float_val(0.25),
    ])

    by_name = {mate.name: mate for mate in asm.data.mates}
    for name in (
        "left_front_axle", "left_middle_axle", "left_rear_axle",
        "right_front_axle", "right_middle_axle", "right_rear_axle",
    ):
        mate = by_name[name]
        parent = asm.data.get_transformed_datum(mate.part_a, mate.datum_a)
        wheel = asm.data.get_transformed_datum(mate.part_b, mate.datum_b)
        np.testing.assert_allclose(parent.origin, wheel.origin, atol=1e-8)
        np.testing.assert_allclose(parent.direction, wheel.direction, atol=1e-8)


def test_positioned_compound_contains_all_eleven_proxy_bodies(proxy_solids):
    asm = _make_assembly(proxy_solids)
    call_builtin("solve_assembly", [asm, string_val("chassis")])
    call_builtin("set_joint_position", [
        asm, string_val("left_rocker_pivot"), float_val(0.2),
    ])

    compound = call_builtin("assembly_compound", [asm]).data
    positioned = asm.data.positioned_parts()

    assert issolid(compound)
    assert len(positioned) == 11
    expected_volume = sum(volumeof(solid) for solid in positioned.values())
    assert volumeof(compound) == pytest.approx(expected_volume, rel=1e-8)
    # The assembled tabletop rover remains compact enough for a workbench.
    # Articulation increases the axis-aligned X envelope beyond the 510 mm
    # level-pose target while remaining a compact workbench assembly.
    assert np.all(_extent(compound) <= [550.0, 346.0, 250.0])


def test_complete_suspension_is_authored_and_evaluated_in_dsl(source):
    result = compile_and_run(
        source, "BUILD_SUSPENSION", {"rocker_angle": math.radians(12.0)}
    )

    assert result.success, result.error_message
    assert issolid(result.geometry)
    assert volumeof(result.geometry) > 0.0


@pytest.mark.requires_occ
@pytest.mark.skipif(not occ_available(), reason="pythonocc-core is not available")
def test_rover_compound_supports_strict_analytic_step(
        tmp_path, source):
    result = compile_and_run(
        source, "BUILD_SUSPENSION", {"rocker_angle": math.radians(12.0)}
    )
    assert result.success, result.error_message
    assert has_brep_data(result.geometry)

    output = tmp_path / "yaprover_suspension.step"
    assert write_step_analytic(
        result.geometry, str(output), fallback_to_faceted=False,
    ) is True
    assert output.exists() and output.stat().st_size > 0
