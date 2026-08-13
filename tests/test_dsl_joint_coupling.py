"""DSL contract tests for affine joint coupling declarations."""

import math
import textwrap

import pytest

from yapcad.dsl import compile_and_run
from yapcad.dsl.runtime.builtins import call_builtin
from yapcad.dsl.runtime.values import float_val, list_val, string_val
from yapcad.dsl.types import FLOAT, STRING
from yapcad.geom3d import issolid

from test_dsl_assembly_solving import _add_named_mate, _add_part, _new_assembly


def _coupled_assembly():
    asm = _new_assembly("differential")
    _add_part(asm, "BASE", "chassis")
    _add_part(asm, "LINK", "left_rocker")
    _add_part(asm, "LINK", "right_rocker")
    _add_named_mate(
        asm, "left", "revolute", "chassis", "joint", "left_rocker", "joint"
    )
    _add_named_mate(
        asm, "right", "revolute", "chassis", "joint", "right_rocker", "joint"
    )
    call_builtin("add_joint_coupling", [
        asm,
        string_val("rocker_differential"),
        string_val("right"),
        list_val([string_val("left")], STRING),
        list_val([float_val(-1.0)], FLOAT),
        float_val(0.0),
    ])
    return asm


def test_dsl_builtin_declares_and_poses_a_differential():
    asm = _coupled_assembly()
    call_builtin("solve_assembly", [asm, string_val("chassis")])
    call_builtin("set_joint_position", [asm, string_val("left"), float_val(0.3)])
    assert asm.data._joint_values["left"] == pytest.approx(0.3)
    assert asm.data._joint_values["right"] == pytest.approx(-0.3)


def test_dsl_builtin_rejects_mismatched_joint_and_coefficient_lists():
    asm = _coupled_assembly()
    with pytest.raises(ValueError, match="same length"):
        call_builtin("add_joint_coupling", [
            asm,
            string_val("bad"),
            string_val("left"),
            list_val([string_val("right")], STRING),
            list_val([], FLOAT),
            float_val(0.0),
        ])


def test_complete_coupled_mechanism_can_be_authored_in_dsl():
    from test_dsl_assembly_solving import PART_SOURCE

    source = PART_SOURCE + """

command BUILD_DIFFERENTIAL(angle: float) -> solid:
    let mechanism: assembly = assembly("dsl_differential")
    add_part(mechanism, BASE(), "chassis")
    add_part(mechanism, LINK(), "left_rocker")
    add_part(mechanism, LINK(), "right_rocker")
    add_named_mate(mechanism, "left", "revolute",
                   "chassis", "joint", "left_rocker", "joint")
    add_named_mate(mechanism, "right", "revolute",
                   "chassis", "joint", "right_rocker", "joint")
    add_joint_coupling(mechanism, "rocker_differential", "right",
                       ["left"], [-1.0], 0.0)
    solve_assembly(mechanism, "chassis")
    set_joint_position(mechanism, "left", angle)
    emit assembly_compound(mechanism)
"""
    result = compile_and_run(
        textwrap.dedent(source), "BUILD_DIFFERENTIAL", {"angle": math.pi / 8}
    )
    assert result.success, result.error_message
    assert issolid(result.geometry)

