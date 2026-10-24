"""
yapcad.rom — Range-of-Motion (kinematic sweep) package.

Overview
========
This package lets you annotate DSL commands with kinematic metadata
(``@meta(kinematic.*)``), build a kinematic tree, and sweep the DOF space
to find pairwise clearance between adjacent links — all without meshing,
using SDF field evaluation on sampled grids.

Quick start
===========
In your .dsl file::

    @meta(kinematic.joint_kind="revolute")
    @meta(kinematic.axis=(0, 0, 1))
    @meta(kinematic.limits=[-90, 90])
    @meta(kinematic.parent="PELVIS")
    @meta(kinematic.parent_frame="HIP_YAW_MOUNT")
    command HIP_YAW_JOINT() -> solid:
        srv = translate(servo.XL430_ENVELOPE(), 0, ho*s, 0)
        emit srv

From Python::

    from yapcad.rom import RomTree, sweep_kinematics
    from yapcad.dsl import compile_to_ast, interpret_module

    ast = compile_to_ast(source)
    ctx = interpret_module(ast)
    tree = RomTree.from_interpreter_context("my_robot", ctx)

    result = sweep_kinematics(
        tree,
        {"HIP_YAW": (-90, 91, 5)},   # (start, stop, step)
        {"KNEE_PITCH": (0, 136, 5)},
    )

    print(result.summary())
    # → "KNEE_PITCH @ 112°: -0.2 mm (thigh_link vs shin_link)"
    # → "HIP_YAW @ ±45..90°: 1.0 mm (pelvis vs hip_bracket)"
"""

from .tree import RomTree, RomLink
from .clearance import clearance_sdf, clearance_mesh
from .sweep import (
    sweep_kinematics,
    SweepResult,
    SweepConfig,
    InterferenceEvent,
    dof_schedule_from_tree,
)
from .report import interference_summary

__all__ = [
    "RomTree",
    "RomLink",
    "clearance_sdf",
    "clearance_mesh",
    "sweep_kinematics",
    "SweepResult",
    "SweepConfig",
    "InterferenceEvent",
    "dof_schedule_from_tree",
    "interference_summary",
]