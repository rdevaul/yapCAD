"""DSL contract tests for revolute mate position limits.

``set_mate_limits`` is additive: mate creation remains unchanged and limits
are attached by stable mate name.  Revolute values are expressed in radians,
matching ``set_joint_position`` and the assembly solver.
"""

import math
import textwrap

import pytest

from yapcad.assembly.assembly import AssemblyError
from yapcad.dsl import compile_and_run
from yapcad.dsl.runtime.builtins import call_builtin
from yapcad.dsl.runtime.values import (
    float_val,
    list_val,
    solid_val,
    string_val,
)
from yapcad.dsl.types import FLOAT, STRING
from yapcad.geom3d import issolid


PART_SOURCE = """
module mate_limit_parts

@meta(assembly.datums=[{
    "id": "pivot",
    "kind": "axis",
    "origin_mm": [0.0, 0.0, 0.0],
    "direction": [0.0, 0.0, 1.0],
}])
command PART() -> solid:
    emit box(4.0, 2.0, 1.0)
"""


def _solid():
    result = compile_and_run(textwrap.dedent(PART_SOURCE), "PART", {})
    assert result.success, result.error_message
    return result.geometry


def _assembly_with_joint(kind="revolute", name="pivot_joint"):
    asm = call_builtin("assembly", [string_val("limited_mechanism")])
    for part_name in ("base", "link"):
        call_builtin(
            "add_part", [asm, solid_val(_solid()), string_val(part_name)]
        )
    call_builtin(
        "add_named_mate",
        [
            asm,
            string_val(name),
            string_val(kind),
            string_val("base"),
            string_val("pivot"),
            string_val("link"),
            string_val("pivot"),
        ],
    )
    return asm


def _set_limits(asm, name, minimum, maximum):
    return call_builtin(
        "set_mate_limits",
        [
            asm,
            string_val(name),
            float_val(minimum),
            float_val(maximum),
        ],
    )


class TestSetMateLimitsBuiltin:
    def test_attaches_radian_limits_and_is_chainable(self):
        asm = _assembly_with_joint()

        result = _set_limits(asm, "pivot_joint", -math.pi / 4, math.pi / 3)

        assert result is asm
        assert asm.data.mates[0].limits.min_value == pytest.approx(-math.pi / 4)
        assert asm.data.mates[0].limits.max_value == pytest.approx(math.pi / 3)

    def test_inclusive_boundaries_are_accepted(self):
        asm = _assembly_with_joint()
        _set_limits(asm, "pivot_joint", -0.5, 0.75)
        call_builtin("solve_assembly", [asm, string_val("base")])

        for angle in (-0.5, 0.75):
            result = call_builtin(
                "set_joint_position",
                [asm, string_val("pivot_joint"), float_val(angle)],
            )
            assert result is asm

    @pytest.mark.parametrize("angle", [-0.500001, 0.750001])
    def test_values_outside_limits_are_rejected(self, angle):
        asm = _assembly_with_joint()
        _set_limits(asm, "pivot_joint", -0.5, 0.75)
        call_builtin("solve_assembly", [asm, string_val("base")])

        with pytest.raises(AssemblyError, match=r"pivot_joint.*limit"):
            call_builtin(
                "set_joint_position",
                [asm, string_val("pivot_joint"), float_val(angle)],
            )

    def test_invalid_ordering_is_rejected_without_mutating_existing_limits(self):
        asm = _assembly_with_joint()
        _set_limits(asm, "pivot_joint", -1.0, 1.0)

        with pytest.raises(ValueError, match=r"minimum.*maximum"):
            _set_limits(asm, "pivot_joint", 2.0, -2.0)

        limits = asm.data.mates[0].limits
        assert (limits.min_value, limits.max_value) == (-1.0, 1.0)

    @pytest.mark.parametrize(
        "minimum, maximum",
        [(math.nan, 1.0), (-1.0, math.inf), (-math.inf, 1.0)],
    )
    def test_non_finite_limits_are_rejected_without_mutation(
        self, minimum, maximum,
    ):
        asm = _assembly_with_joint()
        _set_limits(asm, "pivot_joint", -1.0, 1.0)

        with pytest.raises(ValueError, match="finite"):
            _set_limits(asm, "pivot_joint", minimum, maximum)

        limits = asm.data.mates[0].limits
        assert (limits.min_value, limits.max_value) == (-1.0, 1.0)

    def test_unknown_mate_is_rejected(self):
        asm = _assembly_with_joint()
        with pytest.raises(AssemblyError, match=r"missing.*not found"):
            _set_limits(asm, "missing", -1.0, 1.0)

    def test_nonrevolute_mate_is_rejected(self):
        asm = _assembly_with_joint(kind="rigid", name="fixed_mount")
        with pytest.raises(AssemblyError, match=r"fixed_mount.*rigid.*revolute"):
            _set_limits(asm, "fixed_mount", -1.0, 1.0)


class TestCoupledMateLimits:
    def test_limit_is_enforced_on_derived_dependent_value(self):
        asm = call_builtin("assembly", [string_val("differential")])
        for part_name in ("base", "left", "right"):
            call_builtin(
                "add_part", [asm, solid_val(_solid()), string_val(part_name)]
            )
        for mate_name, child in (("left_joint", "left"), ("right_joint", "right")):
            call_builtin(
                "add_named_mate",
                [
                    asm,
                    string_val(mate_name),
                    string_val("revolute"),
                    string_val("base"),
                    string_val("pivot"),
                    string_val(child),
                    string_val("pivot"),
                ],
            )
        call_builtin(
            "add_joint_coupling",
            [
                asm,
                string_val("differential"),
                string_val("right_joint"),
                list_val([string_val("left_joint")], STRING),
                list_val([float_val(-2.0)], FLOAT),
                float_val(0.0),
            ],
        )
        _set_limits(asm, "right_joint", -0.8, 0.8)
        call_builtin("solve_assembly", [asm, string_val("base")])

        # left=0.4 derives right=-0.8 and is accepted inclusively.
        call_builtin(
            "set_joint_position",
            [asm, string_val("left_joint"), float_val(0.4)],
        )
        assert asm.data._joint_values["right_joint"] == pytest.approx(-0.8)

        # left=0.41 derives right=-0.82, so the dependent joint's limit wins.
        with pytest.raises(AssemblyError, match=r"right_joint.*limit"):
            call_builtin(
                "set_joint_position",
                [asm, string_val("left_joint"), float_val(0.41)],
            )


class TestParsedDslMateLimits:
    def test_limits_can_be_declared_and_exercised_in_parsed_dsl(self):
        source = PART_SOURCE + """

command BUILD(angle: float) -> solid:
    let mechanism: assembly = assembly("limited")
    let base: solid = PART()
    let link: solid = PART()
    add_part(mechanism, base, "base")
    add_part(mechanism, link, "link")
    add_named_mate(mechanism, "pivot_joint", "revolute",
                   "base", "pivot", "link", "pivot")
    set_mate_limits(mechanism, "pivot_joint", -0.75, 0.75)
    solve_assembly(mechanism, "base")
    set_joint_position(mechanism, "pivot_joint", angle)
    emit assembly_compound(mechanism)
"""

        result = compile_and_run(textwrap.dedent(source), "BUILD", {"angle": 0.75})

        assert result.success, result.error_message
        assert issolid(result.geometry)

    def test_parsed_dsl_reports_out_of_range_pose(self):
        source = PART_SOURCE + """

command BUILD(angle: float) -> solid:
    let mechanism: assembly = assembly("limited")
    let base: solid = PART()
    let link: solid = PART()
    add_part(mechanism, base, "base")
    add_part(mechanism, link, "link")
    add_named_mate(mechanism, "pivot_joint", "revolute",
                   "base", "pivot", "link", "pivot")
    set_mate_limits(mechanism, "pivot_joint", -0.75, 0.75)
    solve_assembly(mechanism, "base")
    set_joint_position(mechanism, "pivot_joint", angle)
    emit assembly_compound(mechanism)
"""

        with pytest.raises(AssemblyError, match=r"pivot_joint.*limit"):
            compile_and_run(textwrap.dedent(source), "BUILD", {"angle": 0.8})


def test_legacy_and_named_mate_builtins_remain_registered():
    """Adding limit configuration must not replace either mate constructor."""
    legacy = _assembly_with_joint()
    assert legacy.data.mates[0].name == "pivot_joint"

    asm = call_builtin("assembly", [string_val("legacy")])
    for name in ("a", "b"):
        call_builtin("add_part", [asm, solid_val(_solid()), string_val(name)])
    result = call_builtin(
        "add_mate",
        [
            asm,
            string_val("revolute"),
            string_val("a"),
            string_val("pivot"),
            string_val("b"),
            string_val("pivot"),
        ],
    )
    assert result is asm
