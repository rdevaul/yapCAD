"""Executable contract for solving and posing assemblies built by the DSL.

These tests intentionally describe the next assembly-integration increment.
They use real parsed/interpreted DSL solids and the public builtin dispatcher;
production internals are not mocked.  The proposed additive builtin surface is:

``add_named_mate(asm, name, kind, part_a, datum_a, part_b, datum_b)``
    Add an addressable mate.  The existing six-argument ``add_mate`` remains
    valid and continues to synthesize a name.
``solve_assembly(asm, root_part)``
    Solve a rooted, acyclic mate graph and return the assembly handle.
``set_joint_position(asm, mate_name, value)``
    Set a revolute angle (radians) or prismatic distance (millimetres).
``part_transform(asm, part_name)``
    Return a DSL ``transform`` compatible with ``apply``.
``assembly_compound(asm)``
    Return one positioned yapCAD solid containing every retained part solid.

The deliberately small API separates stable identity (named mates) from pose,
keeps solve explicit, and does not overload the legacy ``add_mate`` signature.
"""

import math
import textwrap

import numpy as np
import pytest

from yapcad.assembly.assembly import AssemblyError
from yapcad.assembly.mate import MateType
from yapcad.dsl import compile_and_run
from yapcad.dsl.runtime.builtins import call_builtin
from yapcad.dsl.runtime.values import float_val, solid_val, string_val
from yapcad.geom3d import issolid, solidbbox, volumeof


PART_SOURCE = """
module assembly_solver_parts

@meta(assembly.datums=[{
    "id": "joint",
    "kind": "axis",
    "origin_mm": [12.5, -7.0, 4.0],
    "direction": [0.0, 0.0, 1.0],
}])
command BASE() -> solid:
    emit box(20.0, 16.0, 8.0)

@meta(assembly.datums=[{
    "id": "joint",
    "kind": "axis",
    "origin_mm": [2.5, 3.0, 1.0],
    "direction": [0.0, 0.0, 1.0],
}])
command LINK() -> solid:
    emit box(10.0, 4.0, 3.0)

@meta(assembly.datums=[{
    "id": "tip",
    "kind": "axis",
    "origin_mm": [10.0, 0.0, 0.0],
    "direction": [0.0, 0.0, 1.0],
}])
command TIP() -> solid:
    emit box(3.0, 3.0, 3.0)
"""


def _dsl_solid(command):
    result = compile_and_run(textwrap.dedent(PART_SOURCE), command, {})
    assert result.success, result.error_message
    assert issolid(result.geometry)
    return result.geometry


def _new_assembly(name="mechanism"):
    return call_builtin("assembly", [string_val(name)])


def _add_part(asm, command, name):
    return call_builtin(
        "add_part", [asm, solid_val(_dsl_solid(command)), string_val(name)]
    )


def _add_named_mate(asm, name, kind, part_a, datum_a, part_b, datum_b):
    return call_builtin(
        "add_named_mate",
        [
            asm,
            string_val(name),
            string_val(kind),
            string_val(part_a),
            string_val(datum_a),
            string_val(part_b),
            string_val(datum_b),
        ],
    )


def _matrix_data(transform_value):
    """Normalize the public transform payload without prescribing storage."""
    matrix = transform_value.data
    if hasattr(matrix, "m"):
        matrix = matrix.m
    return np.asarray(matrix, dtype=float)


def _two_part_assembly(kind="rigid", mate_name="base_to_link"):
    asm = _new_assembly()
    _add_part(asm, "BASE", "base")
    _add_part(asm, "LINK", "link")
    _add_named_mate(
        asm, mate_name, kind, "base", "joint", "link", "joint"
    )
    return asm


class TestPartGeometryAndDatumBridge:
    def test_real_dsl_lifts_arbitrary_3d_datum(self):
        asm = _new_assembly()
        _add_part(asm, "BASE", "base")

        datum = asm.data.parts["base"].datums["joint"]
        assert datum.origin == [12.5, -7.0, 4.0, 1.0]
        assert datum.direction == [0.0, 0.0, 1.0, 0.0]

    def test_add_part_retains_the_interpreted_solid(self):
        asm = _new_assembly()
        source_solid = _dsl_solid("BASE")
        call_builtin(
            "add_part", [asm, solid_val(source_solid), string_val("base")]
        )

        # Geometry must remain available for positioned output.  A defensive
        # copy is acceptable, but a metadata-only placeholder is not.
        retained = asm.data.parts["base"].geometry
        assert issolid(retained)
        assert volumeof(retained) == pytest.approx(volumeof(source_solid))


class TestNamedMateContract:
    def test_six_argument_add_mate_remains_compatible(self):
        asm = _new_assembly()
        _add_part(asm, "BASE", "base")
        _add_part(asm, "LINK", "link")

        result = call_builtin(
            "add_mate",
            [
                asm,
                string_val("revolute"),
                string_val("base"),
                string_val("joint"),
                string_val("link"),
                string_val("joint"),
            ],
        )

        assert result is asm
        assert asm.data.mates[0].mate_type == MateType.REVOLUTE

    @pytest.mark.parametrize("kind", ["rigid", "revolute"])
    def test_named_rigid_and_revolute_mates(self, kind):
        asm = _two_part_assembly(kind=kind, mate_name="suspension_pivot")

        assert asm.data.mates[0].name == "suspension_pivot"
        assert asm.data.mates[0].mate_type == MateType(kind)

    def test_duplicate_mate_name_is_rejected(self):
        asm = _two_part_assembly()
        with pytest.raises(AssemblyError, match=r"base_to_link.*already exists"):
            _add_named_mate(
                asm,
                "base_to_link",
                "rigid",
                "base",
                "joint",
                "link",
                "joint",
            )


class TestRootedSolveAndTransformAccess:
    def test_rigid_solve_places_arbitrary_local_datums_together(self):
        asm = _two_part_assembly()
        result = call_builtin("solve_assembly", [asm, string_val("base")])

        assert result is asm
        np.testing.assert_allclose(asm.data.transforms["base"], np.eye(4))
        base_datum = asm.data.get_transformed_datum("base", "joint")
        link_datum = asm.data.get_transformed_datum("link", "joint")
        np.testing.assert_allclose(link_datum.origin, base_datum.origin, atol=1e-8)
        np.testing.assert_allclose(
            link_datum.direction, base_datum.direction, atol=1e-8
        )

    def test_solve_propagates_over_a_three_part_tree(self):
        asm = _two_part_assembly()
        _add_part(asm, "TIP", "tip")
        _add_named_mate(
            asm, "link_to_tip", "rigid", "link", "joint", "tip", "tip"
        )

        call_builtin("solve_assembly", [asm, string_val("base")])

        link_joint = asm.data.get_transformed_datum("link", "joint")
        tip_joint = asm.data.get_transformed_datum("tip", "tip")
        np.testing.assert_allclose(tip_joint.origin, link_joint.origin, atol=1e-8)

    def test_either_endpoint_can_be_selected_as_root(self):
        asm = _two_part_assembly()

        call_builtin("solve_assembly", [asm, string_val("link")])

        np.testing.assert_allclose(asm.data.transforms["link"], np.eye(4))
        base_joint = asm.data.get_transformed_datum("base", "joint")
        link_joint = asm.data.get_transformed_datum("link", "joint")
        np.testing.assert_allclose(base_joint.origin, link_joint.origin, atol=1e-8)

    def test_part_transform_is_a_dsl_transform_usable_by_apply(self):
        asm = _two_part_assembly()
        call_builtin("solve_assembly", [asm, string_val("base")])

        transform = call_builtin("part_transform", [asm, string_val("link")])

        assert transform.type.name == "transform"
        assert _matrix_data(transform).shape == (4, 4)
        positioned = call_builtin("apply", [transform, solid_val(_dsl_solid("LINK"))])
        assert positioned.type.name == "solid"
        assert issolid(positioned.data)


class TestRevolutePose:
    def test_joint_position_rotates_child_about_the_solved_world_axis(self):
        asm = _two_part_assembly(kind="revolute", mate_name="rocker_pivot")
        call_builtin("solve_assembly", [asm, string_val("base")])
        zero = _matrix_data(
            call_builtin("part_transform", [asm, string_val("link")])
        )

        result = call_builtin(
            "set_joint_position",
            [asm, string_val("rocker_pivot"), float_val(math.pi / 2.0)],
        )
        posed = _matrix_data(
            call_builtin("part_transform", [asm, string_val("link")])
        )

        assert result is asm
        assert not np.allclose(posed, zero)
        # The mate datum stays coincident while the link's orientation changes.
        base_joint = asm.data.get_transformed_datum("base", "joint")
        link_joint = asm.data.get_transformed_datum("link", "joint")
        np.testing.assert_allclose(link_joint.origin, base_joint.origin, atol=1e-8)
        np.testing.assert_allclose(
            posed[:3, :3],
            [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
            atol=1e-8,
        )

    def test_joint_positions_are_absolute_not_accumulated_deltas(self):
        sequential = _two_part_assembly(
            kind="revolute", mate_name="rocker_pivot"
        )
        call_builtin("solve_assembly", [sequential, string_val("base")])
        call_builtin(
            "set_joint_position",
            [sequential, string_val("rocker_pivot"), float_val(math.pi / 4.0)],
        )
        call_builtin(
            "set_joint_position",
            [sequential, string_val("rocker_pivot"), float_val(math.pi / 2.0)],
        )

        direct = _two_part_assembly(kind="revolute", mate_name="rocker_pivot")
        call_builtin("solve_assembly", [direct, string_val("base")])
        call_builtin(
            "set_joint_position",
            [direct, string_val("rocker_pivot"), float_val(math.pi / 2.0)],
        )

        sequential_tf = _matrix_data(
            call_builtin("part_transform", [sequential, string_val("link")])
        )
        direct_tf = _matrix_data(
            call_builtin("part_transform", [direct, string_val("link")])
        )
        np.testing.assert_allclose(sequential_tf, direct_tf, atol=1e-8)

    def test_unknown_or_non_revolute_joint_cannot_be_posed(self):
        rigid = _two_part_assembly(kind="rigid")
        call_builtin("solve_assembly", [rigid, string_val("base")])

        with pytest.raises(AssemblyError, match=r"missing_joint.*not found"):
            call_builtin(
                "set_joint_position",
                [rigid, string_val("missing_joint"), float_val(0.1)],
            )
        with pytest.raises(AssemblyError, match=r"base_to_link.*rigid"):
            call_builtin(
                "set_joint_position",
                [rigid, string_val("base_to_link"), float_val(0.1)],
            )


class TestAssemblyCompound:
    def test_compound_contains_positioned_retained_part_geometry(self):
        asm = _two_part_assembly()
        call_builtin("solve_assembly", [asm, string_val("base")])

        compound = call_builtin("assembly_compound", [asm])

        assert compound.type.name == "solid"
        assert issolid(compound.data)
        expected_volume = (
            volumeof(_dsl_solid("BASE")) + volumeof(_dsl_solid("LINK"))
        )
        assert volumeof(compound.data) == pytest.approx(expected_volume, rel=1e-8)
        bbox = solidbbox(compound.data)
        assert bbox and len(bbox) == 2

    def test_compound_requires_a_solved_assembly(self):
        asm = _two_part_assembly()
        with pytest.raises(AssemblyError, match=r"solve.*assembly"):
            call_builtin("assembly_compound", [asm])


class TestSolveValidationFailures:
    def test_unknown_root_is_rejected(self):
        asm = _two_part_assembly()
        with pytest.raises(AssemblyError, match=r"root.*ghost.*not found"):
            call_builtin("solve_assembly", [asm, string_val("ghost")])

    def test_disconnected_part_is_reported(self):
        asm = _two_part_assembly()
        _add_part(asm, "TIP", "orphan")
        with pytest.raises(AssemblyError, match=r"[Dd]isconnected.*orphan"):
            call_builtin("solve_assembly", [asm, string_val("base")])

    def test_cycle_is_reported_instead_of_silently_overwriting_a_pose(self):
        asm = _two_part_assembly()
        _add_named_mate(
            asm, "return_edge", "rigid", "link", "joint", "base", "joint"
        )
        with pytest.raises(AssemblyError, match=r"[Cc]ycle|over.?constrained"):
            call_builtin("solve_assembly", [asm, string_val("base")])


class TestParsedDslAssemblyWorkflow:
    def test_complete_workflow_can_be_authored_in_dsl(self):
        source = PART_SOURCE + """

command BUILD_MECHANISM(angle: float) -> solid:
    let base: solid = BASE()
    let link: solid = LINK()
    let mechanism: assembly = assembly("dsl_mechanism")
    add_part(mechanism, base, "base")
    add_part(mechanism, link, "link")
    add_named_mate(mechanism, "rocker_pivot", "revolute",
                   "base", "joint", "link", "joint")
    solve_assembly(mechanism, "base")
    set_joint_position(mechanism, "rocker_pivot", angle)
    emit assembly_compound(mechanism)
"""

        result = compile_and_run(
            textwrap.dedent(source), "BUILD_MECHANISM", {"angle": math.pi / 4.0}
        )

        assert result.success, result.error_message
        assert issolid(result.geometry)
        assert volumeof(result.geometry) > 0.0
