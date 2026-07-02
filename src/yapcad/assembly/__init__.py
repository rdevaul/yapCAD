"""Constraint-based assembly system for yapCAD.

This package provides a declarative assembly system inspired by modern CAD
tools (SolidWorks mates, Fusion 360 joints) but tailored for programmatic
3D generative design workflows. It enables definition and validation of
assembly relationships to prevent common orientation errors.

Key Concepts:
    Datum Features:
        Named geometric references (points, axes, planes, frames, circles)
        defined on parts in local coordinates. Datums transform with parts.

    Assembly Mates:
        Relationships between datum features that position parts relative
        to each other. Types include coincident, concentric, parallel,
        perpendicular, tangent, and angle mates.

    Design Constraints:
        High-level rules that validate assembly intent beyond positioning.
        Examples: "motor axis must be tangent to chassis circle",
        "stator face must point inward toward mounting surface".

    Constraint Solver:
        Computes transforms from mates and validates all constraints,
        providing clear error messages when design intent is violated.

Design Principles:
    - Constraints are first-class executable objects, not documentation
    - Fail-fast validation detects invalid assemblies at definition time
    - Separation of concerns: mates position, constraints verify intent
    - Composability: complex assemblies built from validated subassemblies
    - Integration with existing kinematic chain transform system

Quick Start:
    >>> from yapcad.assembly import (
    ...     Datum, DatumType, PartDefinition,
    ...     Mate, MateType,
    ...     Constraint, ConstraintType,
    ...     Assembly
    ... )

    >>> # Define a part with datum features
    >>> motor = PartDefinition("motor")
    >>> motor.add_datum(Datum(
    ...     "axis", DatumType.AXIS,
    ...     origin=(0, 0, 0), direction=(0, 1, 0),
    ...     description="Motor rotation axis"
    ... ))
    >>> motor.add_datum(Datum(
    ...     "stator_face", DatumType.PLANE,
    ...     origin=(0, 5, 0), normal=(0, 1, 0),
    ...     description="Mounting surface"
    ... ))

    >>> # Define assembly with mates
    >>> assembly = Assembly("wheel_assembly")
    >>> assembly.add_part(motor)
    >>> assembly.add_mate(Mate(
    ...     MateType.CONCENTRIC,
    ...     "motor", "axis",
    ...     "bracket", "hole_axis"
    ... ))

    >>> # Add design constraints
    >>> assembly.add_constraint(Constraint(
    ...     ConstraintType.TANGENT_TO_CIRCLE,
    ...     "motor", "axis",
    ...     "chassis", "wheel_orbit",
    ...     description="Motor axis tangent for proper wheel rolling"
    ... ))

    >>> # Solve and validate
    >>> result = assembly.solve()
    >>> if result.is_valid:
    ...     transforms = result.transforms
    >>> else:
    ...     print(result.errors)

Available Classes:
    Datum Types:
        Datum - Named geometric reference feature on a part
        DatumType - Enum: POINT, AXIS, PLANE, FRAME, CIRCLE
        PartDefinition - Part with datum features and metadata

    Mate System:
        Mate - Relationship between two datum features
        MateType - Enum: COINCIDENT, CONCENTRIC, PARALLEL, PERPENDICULAR,
                  TANGENT, ANGLE, ALIGNED, OPPOSED

    Constraint System:
        Constraint - Design rule that validates assembly intent
        ConstraintType - Enum: TANGENT_TO_CIRCLE, POINTS_TOWARD,
                        POINTS_AWAY, AXIS_PARALLEL_TO_PLANE,
                        AXIS_PERPENDICULAR_TO_PLANE, CUSTOM

    Assembly:
        Assembly - Container for parts, mates, and constraints
        AssemblyResult - Validation result with transforms or errors

See Also:
    - Examples: examples/assembly_demo.py
    - Integration: src/yapcad/assembly/solver.py
"""

from __future__ import annotations

from .datum import Datum, DatumType, PartDefinition
from .constraint import (
    Constraint,
    ConstraintType,
    ConstraintResult,
    angle_between_vectors,
    distance_to_point,
    dot_product,
    is_tangent_to_circle,
    is_radial_from_center,
)
from .mate import (
    Mate,
    MateType,
    MateLimits,
    MateDynamics,
    CoincidentResult,
    create_revolute_mate,
    create_prismatic_mate,
    create_gear_mate,
    create_coincident_mate,
    evaluate_coincident,
    check_bolt_circle_alignment,
)
from .assembly import (
    Assembly,
    AssemblyError,
    AssemblyValidationResult,
)
from .kinematic_integration import (
    KinematicConstraint,
    KinematicConstraintType,
    ConstraintEvaluationResult,
    ValidationReport,
    AssemblyValidator,
    validate_assembly,
    create_tangent_constraint,
    create_parallel_constraint,
    create_coincident_constraint,
    create_z_stack_clearance_constraint,
    create_no_overlap_constraint,
    create_min_distance_constraint,
)

# Declarative assembly intent system
from .intent import (
    # Reference Geometry
    GeometryType,
    ReferenceGeometry,
    # Functional Requirements
    FunctionalRequirement,
    ContactRequirement,
    RollRequirement,
    AxisOrientationRequirement,
    ParallelAxesRequirement,
    ReachRequirement,
    # Connections and Clearances
    Connection,
    Clearance,
    ClearanceResult,
    # Solve Result
    SolveResult,
    # Main Class
    AssemblyIntent,
    # Example Factories
    create_wheel_assembly_intent,
    create_scara_arm_intent,
)

# Datum registry and face-to-face mate system
from .datum_registry import (
    DatumRegistry,
    DatumSource,
    datum_to_transform_matrix,
)
from .face_mate import (
    FaceToFaceMate,
    MateValidationResult,
    compute_flush_transform,
    compute_concentric_transform,
    compute_parallel_transform,
    create_face_mate,
    create_axis_mate,
)
from .load_case import (
    LoadCase,
    LoadAttach,
    BoltPattern,
    COORDINATE_FRAME_VEHICLE,
    COORDINATE_FRAME_MESH,
    STATUS_CONFIRMED,
    STATUS_ESTIMATED,
    STATUS_TRIBAL,
    STATUS_PLACEHOLDER,
    GROUP_INERTIAL,
    GROUP_RECOVERY,
    GROUP_PRESSURE,
)

from .resolver import (
    AssemblyResolver,
    BasicResolver,
    StagedResolver,
    KinematicResolver,
    ResolvedPart,
    ResolveResult,
    RESOLVER_REGISTRY,
    get_resolver,
    resolve_assembly,
)

# Note: Additional modules will be implemented in separate tasks:
# - solver.py: Full constraint solver

__all__ = [
    # Datum system (implemented)
    "Datum",
    "DatumType",
    "PartDefinition",
    # Constraint system (implemented)
    "Constraint",
    "ConstraintType",
    "ConstraintResult",
    "angle_between_vectors",
    "distance_to_point",
    "dot_product",
    "is_tangent_to_circle",
    "is_radial_from_center",
    # Mate system (implemented)
    "Mate",
    "MateType",
    "MateLimits",
    "MateDynamics",
    "CoincidentResult",
    "create_revolute_mate",
    "create_prismatic_mate",
    "create_gear_mate",
    "create_coincident_mate",
    "evaluate_coincident",
    "check_bolt_circle_alignment",
    # Assembly system (implemented)
    "Assembly",
    "AssemblyError",
    "AssemblyValidationResult",
    # Kinematic chain integration
    "KinematicConstraint",
    "KinematicConstraintType",
    "ConstraintEvaluationResult",
    "ValidationReport",
    "AssemblyValidator",
    "validate_assembly",
    "create_tangent_constraint",
    "create_parallel_constraint",
    "create_coincident_constraint",
    "create_z_stack_clearance_constraint",
    "create_no_overlap_constraint",
    "create_min_distance_constraint",
    # Declarative assembly intent system
    "GeometryType",
    "ReferenceGeometry",
    "FunctionalRequirement",
    "ContactRequirement",
    "RollRequirement",
    "AxisOrientationRequirement",
    "ParallelAxesRequirement",
    "ReachRequirement",
    "Connection",
    "Clearance",
    "ClearanceResult",
    "SolveResult",
    "AssemblyIntent",
    "create_wheel_assembly_intent",
    "create_scara_arm_intent",
    # Datum registry and face-to-face mates
    "DatumRegistry",
    "DatumSource",
    "datum_to_transform_matrix",
    "FaceToFaceMate",
    "MateValidationResult",
    # Load case + bolt pattern system (added 2026-05-20)
    "LoadCase",
    "LoadAttach",
    "BoltPattern",
    "COORDINATE_FRAME_VEHICLE",
    "COORDINATE_FRAME_MESH",
    "STATUS_CONFIRMED",
    "STATUS_ESTIMATED",
    "STATUS_TRIBAL",
    "STATUS_PLACEHOLDER",
    "GROUP_INERTIAL",
    "GROUP_RECOVERY",
    "GROUP_PRESSURE",
    "compute_flush_transform",
    "compute_concentric_transform",
    "compute_parallel_transform",
    "create_face_mate",
    "create_axis_mate",
    # Assembly operation resolver (Step 6 — v1.1)
    "AssemblyResolver",
    "BasicResolver",
    "StagedResolver",
    "KinematicResolver",
    "ResolvedPart",
    "ResolveResult",
    "RESOLVER_REGISTRY",
    "get_resolver",
    "resolve_assembly",
]

__version__ = "0.1.0"
