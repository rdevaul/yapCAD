"""
Assembly mate/constraint system for yapCAD.

This module provides kinematic constraints for defining relationships between
components in an assembly. Mates define both static positioning and allowed
motion, supporting mechanical design, animation export, and physics simulation.

Based on industry standards from SolidWorks, Fusion 360, CATIA, Siemens NX,
and PTC Creo, adapted for yapCAD's programmatic workflow.

Key concepts:
    - **Mate**: Geometric relationship that constrains degrees of freedom (DOF)
    - **DOF Removal**: Each mate removes one or more of 6 DOF (3 trans + 3 rot)
    - **Limits**: Min/max position or angle constraints for joints
    - **Dynamics**: Friction, damping, stiffness for motion simulation
    - **Coupling**: Gear ratios and other motion relationships

Example usage:

    from yapcad.assembly.mate import Mate, MateType, MateLimits
    from yapcad.geom import point, vect

    # Define a revolute joint (hinge) for a robot arm
    shoulder = Mate(
        name="shoulder_pitch",
        mate_type=MateType.REVOLUTE,
        part_a="base",
        datum_a="shoulder_axis",
        part_b="upper_arm",
        datum_b="arm_root_axis",
        offset=0.0,
        angle=0.0,
        limits=MateLimits(
            min_value=-1.57,  # -90 degrees
            max_value=1.57,   # +90 degrees
            max_velocity=2.0  # rad/s
        )
    )

    # Check degrees of freedom
    dof = shoulder.degrees_of_freedom()  # Returns 1 (rotation only)

    # Evaluate constraint satisfaction
    result = shoulder.evaluate()
    print(f"DOF remaining: {result['dof_remaining']}")
    print(f"Constraint error: {result['error']}")

For animation and simulation export, see the CAD_MATE_SYSTEMS_RESEARCH.md
document for URDF, SDF, and Blender armature generation strategies.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional, Dict, Any, List, Tuple, TYPE_CHECKING
import math

if TYPE_CHECKING:
    from .datum import Datum, DatumType

# Import numpy with fallback
try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None


class MateType(Enum):
    """
    Kinematic mate types based on industry CAD standards.

    Each mate type removes specific degrees of freedom (DOF) from a rigid body.
    An unconstrained body has 6 DOF: 3 translational + 3 rotational.

    Basic Mates (geometric constraints):
    ------------------------------------
    COINCIDENT: Points, edges, or faces share same location (1-3 DOF removed)
    CONCENTRIC: Cylindrical axes are colinear (4 DOF removed: 2t + 2r)
    PARALLEL: Directions remain parallel (2 DOF removed)
    PERPENDICULAR: Directions maintain 90-degree angle (1 DOF removed)
    TANGENT: Surfaces remain tangent (1 DOF removed)
    DISTANCE: Fixed offset between datums (1 DOF removed)
    ANGLE: Fixed angular relationship (1 DOF removed)

    Standard Joints (motion primitives):
    ------------------------------------
    RIGID: No relative motion (6 DOF removed, 0 remaining)
    REVOLUTE: Rotation about single axis (5 DOF removed, 1 remaining)
    PRISMATIC: Translation along single axis (5 DOF removed, 1 remaining)
    CYLINDRICAL: Rotation + translation on same axis (4 DOF removed, 2 remaining)
    SPHERICAL: Ball-and-socket, rotation only (3 DOF removed, 3 remaining)
    PLANAR: Two translation + one rotation in plane (3 DOF removed, 3 remaining)

    Compound Joints:
    ----------------
    PIN_SLOT: Translation along slot + rotation about pin (4 DOF removed, 2 remaining)
    UNIVERSAL: Two perpendicular rotation axes (4 DOF removed, 2 remaining)
    SCREW: Coupled rotation and translation (5 DOF removed, 1 coupled)

    Coupled Mates (motion relationships):
    -------------------------------------
    GEAR: Coupled rotation with ratio (no DOF removed, creates relationship)
    RACK_PINION: Couples rotation to translation (no DOF removed, creates relationship)
    CAM: Follower constrained to cam profile path (path-following)
    SLOT: Linear sliding constraint along trajectory (varies by slot type)

    References:
        - SolidWorks Mates Overview (2025)
        - Fusion 360 Joint Types
        - CATIA DMU Kinematics
        - Siemens NX Motion Joints
        - PTC Creo Mechanism Design
    """

    # Geometric constraints (standard mates)
    COINCIDENT = "coincident"
    CONCENTRIC = "concentric"
    PARALLEL = "parallel"
    PERPENDICULAR = "perpendicular"
    TANGENT = "tangent"
    DISTANCE = "distance"
    ANGLE = "angle"

    # Standard joints (motion primitives)
    RIGID = "rigid"
    REVOLUTE = "revolute"
    PRISMATIC = "prismatic"
    CYLINDRICAL = "cylindrical"
    SPHERICAL = "spherical"
    PLANAR = "planar"

    # Compound joints
    PIN_SLOT = "pin_slot"
    UNIVERSAL = "universal"
    SCREW = "screw"

    # Coupled mates (motion relationships)
    GEAR = "gear"
    RACK_PINION = "rack_pinion"
    CAM = "cam"
    SLOT = "slot"


@dataclass
class MateLimits:
    """
    Position, velocity, and effort limits for mate degrees of freedom.

    Used for both mechanical design (hard stops) and motion simulation
    (soft limits, compliance). All limits are optional; None means unlimited.

    For revolute joints:
        - min_value/max_value in radians
        - max_velocity in rad/s
        - max_effort in N*m (torque)

    For prismatic joints:
        - min_value/max_value in mm
        - max_velocity in mm/s
        - max_effort in N (force)

    For multi-DOF joints (cylindrical, spherical, planar):
        - Use min_value/max_value for primary DOF
        - Use min_secondary/max_secondary for secondary DOF

    Attributes:
        min_value: Minimum position (mm) or angle (radians) for primary DOF
        max_value: Maximum position (mm) or angle (radians) for primary DOF
        min_velocity: Minimum velocity (can be negative for reversal)
        max_velocity: Maximum velocity (mm/s or rad/s)
        min_effort: Minimum force/torque (N or N*m)
        max_effort: Maximum force/torque (N or N*m)
        min_secondary: Minimum value for secondary DOF (multi-DOF joints)
        max_secondary: Maximum value for secondary DOF (multi-DOF joints)
        limit_stiffness: Stiffness of soft limit (N/m or N*m/rad), 1e6 = hard
        limit_damping: Damping at limit stop (N*s/m or N*m*s/rad)
        restitution: Bounce coefficient at limits (0.0 = no bounce, 1.0 = perfect)

    Example:
        # Revolute joint with ±90 degree limits
        limits = MateLimits(
            min_value=-math.pi/2,
            max_value=math.pi/2,
            max_velocity=2.0,      # 2 rad/s
            max_effort=50.0        # 50 N*m torque limit
        )

        # Prismatic joint (linear slide) with hard stops
        limits = MateLimits(
            min_value=0.0,
            max_value=100.0,       # 100mm travel
            max_velocity=50.0,     # 50 mm/s
            limit_stiffness=1e6    # Very stiff stop
        )
    """

    min_value: Optional[float] = None
    max_value: Optional[float] = None
    min_velocity: Optional[float] = None
    max_velocity: Optional[float] = None
    min_effort: Optional[float] = None
    max_effort: Optional[float] = None
    min_secondary: Optional[float] = None
    max_secondary: Optional[float] = None
    limit_stiffness: float = 1e6
    limit_damping: float = 1e3
    restitution: float = 0.0


@dataclass
class MateDynamics:
    """
    Physical dynamics parameters for motion simulation.

    Defines friction, damping, and compliance for realistic motion behavior.
    Used for physics simulation, animation, and digital twin applications.

    Friction types:
        - Static (Coulomb): Resistance to start motion
        - Kinetic (dynamic): Resistance during motion
        - Viscous: Velocity-dependent resistance (linear with speed)

    Typical values:
        - Lubricated metal bearings: friction_static=0.01-0.05
        - Dry metal: friction_static=0.15-0.30
        - Typical joint damping: 0.1-1.0
        - Stiffness: 0 for rigid, >0 for compliant/springy

    Attributes:
        friction_static: Static (Coulomb) friction coefficient
        friction_kinetic: Kinetic/dynamic friction coefficient
        friction_viscous: Viscous friction coefficient (velocity-dependent)
        damping: Viscous damping coefficient (energy dissipation)
        stiffness: Spring stiffness for compliant joints (N/m or N*m/rad)
        rest_position: Equilibrium position for spring (mm or radians)

    Example:
        # Low-friction bearing with slight damping
        dynamics = MateDynamics(
            friction_static=0.02,
            friction_kinetic=0.015,
            damping=0.5
        )

        # Compliant joint with spring return
        dynamics = MateDynamics(
            stiffness=1000.0,      # Spring constant
            rest_position=0.0,     # Returns to center
            damping=50.0           # Energy dissipation
        )
    """

    friction_static: float = 0.0
    friction_kinetic: float = 0.0
    friction_viscous: float = 0.0
    damping: float = 0.0
    stiffness: float = 0.0
    rest_position: float = 0.0


@dataclass
class Mate:
    """
    Kinematic constraint defining relationship between two assembly components.

    A mate constrains the relative position and/or orientation of two parts
    by removing degrees of freedom (DOF). Mates can be purely geometric
    (COINCIDENT, PARALLEL) or define motion primitives (REVOLUTE, PRISMATIC).

    The mate references two parts via named datum features (points, axes,
    planes, or surfaces). The constraint solver uses these datums to compute
    the relative transformation that satisfies the mate.

    Attributes:
        name: Human-readable identifier (e.g., "shoulder_pitch", "wheel_axle")
        mate_type: Type of constraint from MateType enum
        part_a: Identifier of first (parent/base) component
        datum_a: Named datum feature on part_a (e.g., "mounting_axis")
        part_b: Identifier of second (child/moving) component
        datum_b: Named datum feature on part_b (e.g., "joint_axis")
        offset: Distance offset between datums (mm, used by DISTANCE mate)
        angle: Angular offset between datums (radians, used by ANGLE mate)
        axis: Primary motion axis for joints (e.g., [0,0,1,0] for Z-axis)
        secondary_axis: Secondary reference axis for compound joints
        limits: Optional position/velocity/effort limits
        dynamics: Optional friction/damping/stiffness parameters
        coupling_ratio: Motion ratio for coupled mates (GEAR, SCREW)
        coupling_offset: Phase offset for coupled motion (radians or mm)
        coupling_reverse: Reverse direction of coupled motion
        coupling_pitch: Thread pitch for SCREW mates (mm per revolution)
        metadata: Additional application-specific data

    Example:
        # Revolute joint for robot shoulder
        shoulder_mate = Mate(
            name="shoulder_pitch",
            mate_type=MateType.REVOLUTE,
            part_a="robot_base",
            datum_a="shoulder_mount_axis",
            part_b="upper_arm",
            datum_b="arm_root_axis",
            offset=0.0,
            angle=0.0,
            axis=[0, 1, 0, 0],  # Y-axis rotation
            limits=MateLimits(
                min_value=-math.pi/2,
                max_value=math.pi/2,
                max_velocity=1.5,
                max_effort=100.0
            ),
            dynamics=MateDynamics(
                friction_static=0.05,
                damping=0.2
            )
        )

        # Gear coupling between two shafts
        gear_mate = Mate(
            name="gear_1_to_2",
            mate_type=MateType.GEAR,
            part_a="gear_1",
            datum_a="gear_1_axis",
            part_b="gear_2",
            datum_b="gear_2_axis",
            coupling_ratio=2.5,  # gear_2 rotates 2.5x for each rotation of gear_1
            coupling_reverse=True  # Opposite rotation direction
        )

        # Check constraint properties
        dof = shoulder_mate.degrees_of_freedom()  # Returns 1
        result = shoulder_mate.evaluate()
    """

    name: str
    mate_type: MateType
    part_a: str
    datum_a: str
    part_b: str
    datum_b: str
    offset: float = 0.0
    angle: float = 0.0
    axis: List[float] = field(default_factory=lambda: [0, 0, 1, 0])
    secondary_axis: List[float] = field(default_factory=lambda: [1, 0, 0, 0])
    limits: Optional[MateLimits] = None
    dynamics: Optional[MateDynamics] = None
    coupling_ratio: float = 1.0
    coupling_offset: float = 0.0
    coupling_reverse: bool = False
    coupling_pitch: Optional[float] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    # Properties for alternative naming (part1/part2/datum1/datum2)
    # These provide compatibility with code that uses the alternate naming convention
    @property
    def part1(self) -> str:
        """Alias for part_a (first/parent component)."""
        return self.part_a

    @property
    def part2(self) -> str:
        """Alias for part_b (second/child component)."""
        return self.part_b

    @property
    def datum1(self) -> str:
        """Alias for datum_a (datum on first component)."""
        return self.datum_a

    @property
    def datum2(self) -> str:
        """Alias for datum_b (datum on second component)."""
        return self.datum_b

    def validate(self, datum_a: 'Datum', datum_b: 'Datum') -> List[str]:
        """Validate that this mate is compatible with the given datums.

        Checks that the mate type is appropriate for the datum types.
        For example, CONCENTRIC requires two AXIS or CIRCLE datums,
        COINCIDENT can work with POINT, PLANE, or CIRCLE datums.

        Args:
            datum_a: First datum feature
            datum_b: Second datum feature

        Returns:
            List of validation error messages (empty if valid)

        Example:
            >>> mate = Mate("test", MateType.CONCENTRIC, ...)
            >>> issues = mate.validate(axis_datum, point_datum)
            >>> if issues:
            ...     print(f"Invalid: {issues}")
        """
        from .datum import DatumType
        issues = []

        # Define valid datum type combinations for each mate type
        if self.mate_type == MateType.COINCIDENT:
            # COINCIDENT works with:
            # - Point-to-Point: centers coincide
            # - Point-to-Plane: point lies on plane
            # - Plane-to-Plane: planes are coplanar
            # - Circle-to-Circle: centers coincide (bolt circle alignment)
            # - Axis-to-Axis: axes are colinear (special case)
            valid_a = {DatumType.POINT, DatumType.PLANE, DatumType.CIRCLE,
                       DatumType.AXIS, DatumType.FRAME}
            valid_b = {DatumType.POINT, DatumType.PLANE, DatumType.CIRCLE,
                       DatumType.AXIS, DatumType.FRAME}
            if datum_a.datum_type not in valid_a:
                issues.append(
                    f"COINCIDENT datum_a must be POINT, PLANE, CIRCLE, AXIS, or FRAME, "
                    f"got {datum_a.datum_type.value}"
                )
            if datum_b.datum_type not in valid_b:
                issues.append(
                    f"COINCIDENT datum_b must be POINT, PLANE, CIRCLE, AXIS, or FRAME, "
                    f"got {datum_b.datum_type.value}"
                )

        elif self.mate_type == MateType.CONCENTRIC:
            # CONCENTRIC requires AXIS or CIRCLE datums
            valid = {DatumType.AXIS, DatumType.CIRCLE}
            if datum_a.datum_type not in valid:
                issues.append(
                    f"CONCENTRIC datum_a must be AXIS or CIRCLE, "
                    f"got {datum_a.datum_type.value}"
                )
            if datum_b.datum_type not in valid:
                issues.append(
                    f"CONCENTRIC datum_b must be AXIS or CIRCLE, "
                    f"got {datum_b.datum_type.value}"
                )

        elif self.mate_type in (MateType.PARALLEL, MateType.PERPENDICULAR):
            # Require datums with direction (AXIS, PLANE)
            valid = {DatumType.AXIS, DatumType.PLANE}
            if datum_a.datum_type not in valid:
                issues.append(
                    f"{self.mate_type.value.upper()} datum_a must be AXIS or PLANE, "
                    f"got {datum_a.datum_type.value}"
                )
            if datum_b.datum_type not in valid:
                issues.append(
                    f"{self.mate_type.value.upper()} datum_b must be AXIS or PLANE, "
                    f"got {datum_b.datum_type.value}"
                )

        elif self.mate_type == MateType.TANGENT:
            # TANGENT typically requires at least one PLANE or curved surface
            valid = {DatumType.PLANE, DatumType.CIRCLE}
            if datum_a.datum_type not in valid and datum_b.datum_type not in valid:
                issues.append(
                    f"TANGENT requires at least one PLANE or CIRCLE datum"
                )

        return issues

    def degrees_of_freedom(self) -> int:
        """
        Return the number of degrees of freedom (DOF) remaining after this mate.

        An unconstrained rigid body has 6 DOF (3 translational + 3 rotational).
        Each mate removes one or more DOF. This method returns how many DOF
        remain after applying this mate constraint.

        Returns:
            Number of remaining DOF (0-6)

        DOF by mate type:
            RIGID: 0 DOF (fully constrained)
            REVOLUTE, PRISMATIC, SCREW: 1 DOF
            CYLINDRICAL, PIN_SLOT, UNIVERSAL: 2 DOF
            SPHERICAL, PLANAR: 3 DOF
            COINCIDENT: depends on geometry (1-3 DOF removed)
            CONCENTRIC: 2 DOF (translation along + rotation about axis)
            DISTANCE, ANGLE, PARALLEL, PERPENDICULAR, TANGENT: varies
            GEAR, RACK_PINION, CAM: coupled (creates relationship, not DOF)

        Example:
            shoulder = Mate(name="shoulder", mate_type=MateType.REVOLUTE, ...)
            dof = shoulder.degrees_of_freedom()  # Returns 1
        """
        dof_map = {
            # Basic joints
            MateType.RIGID: 0,
            MateType.REVOLUTE: 1,
            MateType.PRISMATIC: 1,
            MateType.CYLINDRICAL: 2,
            MateType.SPHERICAL: 3,
            MateType.PLANAR: 3,

            # Compound joints
            MateType.PIN_SLOT: 2,
            MateType.UNIVERSAL: 2,
            MateType.SCREW: 1,

            # Geometric constraints (approximate - depends on combination)
            MateType.COINCIDENT: 3,  # Typically removes 3 DOF
            MateType.CONCENTRIC: 2,  # Translation + rotation on axis
            MateType.PARALLEL: 4,    # 2 rotational DOF removed
            MateType.PERPENDICULAR: 5,  # 1 rotational DOF removed
            MateType.TANGENT: 5,     # 1 DOF removed
            MateType.DISTANCE: 5,    # 1 translational DOF removed
            MateType.ANGLE: 5,       # 1 rotational DOF removed

            # Coupled mates (these create relationships, not DOF reduction)
            MateType.GEAR: -1,       # Special: coupled motion
            MateType.RACK_PINION: -1,
            MateType.CAM: -1,
            MateType.SLOT: 4,        # Point constrained to curve
        }

        return dof_map.get(self.mate_type, 6)

    def evaluate(self, current_position: Optional[float] = None,
                 current_velocity: Optional[float] = None,
                 current_effort: Optional[float] = None) -> Dict[str, Any]:
        """
        Evaluate constraint satisfaction and compute error metrics.

        This method checks if the mate is satisfied given current state,
        computes constraint violation errors, and validates limits.

        Args:
            current_position: Current position (mm) or angle (radians) of primary DOF
            current_velocity: Current velocity (mm/s or rad/s)
            current_effort: Current force (N) or torque (N*m)

        Returns:
            Dictionary containing:
                - dof_remaining: Number of DOF after constraint (int)
                - error: Constraint violation error magnitude (float)
                - satisfied: True if constraint is satisfied within tolerance (bool)
                - limit_violations: List of violated limits (List[str])
                - dynamics_active: True if dynamics parameters are defined (bool)

        Example:
            result = mate.evaluate(current_position=0.5, current_velocity=1.2)
            if not result['satisfied']:
                print(f"Constraint error: {result['error']}")
            if result['limit_violations']:
                print(f"Limit violations: {result['limit_violations']}")
        """
        result = {
            'dof_remaining': self.degrees_of_freedom(),
            'error': 0.0,
            'satisfied': True,
            'limit_violations': [],
            'dynamics_active': self.dynamics is not None
        }

        # Check position limits
        if current_position is not None and self.limits is not None:
            if self.limits.min_value is not None and current_position < self.limits.min_value:
                violation = self.limits.min_value - current_position
                result['error'] += abs(violation)
                result['satisfied'] = False
                result['limit_violations'].append(
                    f"min_value: {current_position:.4f} < {self.limits.min_value:.4f}"
                )

            if self.limits.max_value is not None and current_position > self.limits.max_value:
                violation = current_position - self.limits.max_value
                result['error'] += abs(violation)
                result['satisfied'] = False
                result['limit_violations'].append(
                    f"max_value: {current_position:.4f} > {self.limits.max_value:.4f}"
                )

        # Check velocity limits
        if current_velocity is not None and self.limits is not None:
            if self.limits.min_velocity is not None and current_velocity < self.limits.min_velocity:
                result['limit_violations'].append(
                    f"min_velocity: {current_velocity:.4f} < {self.limits.min_velocity:.4f}"
                )

            if self.limits.max_velocity is not None and abs(current_velocity) > self.limits.max_velocity:
                result['limit_violations'].append(
                    f"max_velocity: |{current_velocity:.4f}| > {self.limits.max_velocity:.4f}"
                )

        # Check effort limits
        if current_effort is not None and self.limits is not None:
            if self.limits.min_effort is not None and current_effort < self.limits.min_effort:
                result['limit_violations'].append(
                    f"min_effort: {current_effort:.4f} < {self.limits.min_effort:.4f}"
                )

            if self.limits.max_effort is not None and abs(current_effort) > self.limits.max_effort:
                result['limit_violations'].append(
                    f"max_effort: |{current_effort:.4f}| > {self.limits.max_effort:.4f}"
                )

        return result

    def is_coupled(self) -> bool:
        """
        Return True if this mate defines a coupled motion relationship.

        Coupled mates (GEAR, RACK_PINION, SCREW, CAM) create motion
        relationships between components rather than removing DOF.

        Returns:
            True if mate is a coupled motion type

        Example:
            if mate.is_coupled():
                print(f"Coupling ratio: {mate.coupling_ratio}")
        """
        return self.mate_type in (
            MateType.GEAR,
            MateType.RACK_PINION,
            MateType.CAM,
            MateType.SCREW
        )

    def compute_coupled_motion(self, driver_position: float) -> float:
        """
        Compute driven position from driver position for coupled mates.

        For coupled mates (GEAR, SCREW, RACK_PINION), compute the position
        of the driven component given the position of the driving component.

        Args:
            driver_position: Position of driving component (mm or radians)

        Returns:
            Position of driven component (mm or radians)

        Raises:
            ValueError: If mate is not a coupled type

        Example:
            # Gear mate with 3:1 ratio
            gear_mate = Mate(mate_type=MateType.GEAR, coupling_ratio=3.0, ...)
            driven_angle = gear_mate.compute_coupled_motion(1.0)  # Returns 3.0

            # Screw mate with 2mm pitch
            screw_mate = Mate(mate_type=MateType.SCREW, coupling_pitch=2.0, ...)
            linear_pos = screw_mate.compute_coupled_motion(math.pi)  # Returns 2.0
        """
        if not self.is_coupled():
            raise ValueError(f"Mate {self.name} is not a coupled type")

        # Apply coupling ratio and offset
        driven = driver_position * self.coupling_ratio + self.coupling_offset

        # Apply reversal if specified
        if self.coupling_reverse:
            driven = -driven

        # For SCREW mates, use pitch instead of ratio
        if self.mate_type == MateType.SCREW and self.coupling_pitch is not None:
            # Position = rotations * pitch
            # driver_position is in radians, convert to rotations
            rotations = driver_position / (2.0 * math.pi)
            driven = rotations * self.coupling_pitch

        return driven

    def to_dict(self) -> Dict[str, Any]:
        """
        Serialize mate to dictionary for export to JSON, YAML, or other formats.

        Returns:
            Dictionary representation of mate with all parameters

        Example:
            mate_dict = mate.to_dict()
            import json
            json.dump(mate_dict, f, indent=2)
        """
        return {
            'name': self.name,
            'mate_type': self.mate_type.value,
            'part_a': self.part_a,
            'datum_a': self.datum_a,
            'part_b': self.part_b,
            'datum_b': self.datum_b,
            'offset': self.offset,
            'angle': self.angle,
            'axis': self.axis,
            'secondary_axis': self.secondary_axis,
            'limits': {
                'min_value': self.limits.min_value,
                'max_value': self.limits.max_value,
                'min_velocity': self.limits.min_velocity,
                'max_velocity': self.limits.max_velocity,
                'min_effort': self.limits.min_effort,
                'max_effort': self.limits.max_effort,
                'min_secondary': self.limits.min_secondary,
                'max_secondary': self.limits.max_secondary,
                'limit_stiffness': self.limits.limit_stiffness,
                'limit_damping': self.limits.limit_damping,
                'restitution': self.limits.restitution,
            } if self.limits else None,
            'dynamics': {
                'friction_static': self.dynamics.friction_static,
                'friction_kinetic': self.dynamics.friction_kinetic,
                'friction_viscous': self.dynamics.friction_viscous,
                'damping': self.dynamics.damping,
                'stiffness': self.dynamics.stiffness,
                'rest_position': self.dynamics.rest_position,
            } if self.dynamics else None,
            'coupling_ratio': self.coupling_ratio,
            'coupling_offset': self.coupling_offset,
            'coupling_reverse': self.coupling_reverse,
            'coupling_pitch': self.coupling_pitch,
            'metadata': self.metadata,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'Mate':
        """
        Deserialize mate from dictionary.

        Args:
            data: Dictionary representation of mate

        Returns:
            Mate instance

        Example:
            import json
            mate_dict = json.load(f)
            mate = Mate.from_dict(mate_dict)
        """
        # Convert mate_type string to enum
        mate_type = MateType(data['mate_type'])

        # Reconstruct limits if present
        limits = None
        if data.get('limits'):
            limits = MateLimits(**data['limits'])

        # Reconstruct dynamics if present
        dynamics = None
        if data.get('dynamics'):
            dynamics = MateDynamics(**data['dynamics'])

        return cls(
            name=data['name'],
            mate_type=mate_type,
            part_a=data['part_a'],
            datum_a=data['datum_a'],
            part_b=data['part_b'],
            datum_b=data['datum_b'],
            offset=data.get('offset', 0.0),
            angle=data.get('angle', 0.0),
            axis=data.get('axis', [0, 0, 1, 0]),
            secondary_axis=data.get('secondary_axis', [1, 0, 0, 0]),
            limits=limits,
            dynamics=dynamics,
            coupling_ratio=data.get('coupling_ratio', 1.0),
            coupling_offset=data.get('coupling_offset', 0.0),
            coupling_reverse=data.get('coupling_reverse', False),
            coupling_pitch=data.get('coupling_pitch'),
            metadata=data.get('metadata', {}),
        )


def create_revolute_mate(
    name: str,
    part_a: str,
    datum_a: str,
    part_b: str,
    datum_b: str,
    axis: Optional[List[float]] = None,
    min_angle: Optional[float] = None,
    max_angle: Optional[float] = None,
    max_velocity: Optional[float] = None,
    max_torque: Optional[float] = None,
    friction: float = 0.0,
    damping: float = 0.0
) -> Mate:
    """
    Convenience function to create a revolute (hinge) joint mate.

    A revolute joint allows rotation about a single axis. This is the most
    common joint type for mechanisms, robots, and articulated assemblies.

    Args:
        name: Descriptive name (e.g., "shoulder_pitch", "door_hinge")
        part_a: Parent/base component identifier
        datum_a: Axis datum on parent component
        part_b: Child/moving component identifier
        datum_b: Axis datum on child component
        axis: Rotation axis direction [x,y,z,w], defaults to [0,0,1,0] (Z-axis)
        min_angle: Minimum rotation angle in radians (None = unlimited)
        max_angle: Maximum rotation angle in radians (None = unlimited)
        max_velocity: Maximum angular velocity in rad/s (None = unlimited)
        max_torque: Maximum torque in N*m (None = unlimited)
        friction: Static friction coefficient (0.0 = frictionless)
        damping: Viscous damping coefficient (0.0 = no damping)

    Returns:
        Mate configured as revolute joint

    Example:
        # Robot elbow with 180-degree range
        elbow = create_revolute_mate(
            name="elbow_flex",
            part_a="upper_arm",
            datum_a="elbow_axis",
            part_b="forearm",
            datum_b="forearm_root_axis",
            min_angle=0.0,
            max_angle=math.pi,
            max_velocity=2.0,
            friction=0.02,
            damping=0.1
        )
    """
    limits = None
    if any(x is not None for x in [min_angle, max_angle, max_velocity, max_torque]):
        limits = MateLimits(
            min_value=min_angle,
            max_value=max_angle,
            max_velocity=max_velocity,
            max_effort=max_torque
        )

    dynamics = None
    if friction > 0.0 or damping > 0.0:
        dynamics = MateDynamics(
            friction_static=friction,
            friction_kinetic=friction * 0.8,  # Kinetic typically 80% of static
            damping=damping
        )

    return Mate(
        name=name,
        mate_type=MateType.REVOLUTE,
        part_a=part_a,
        datum_a=datum_a,
        part_b=part_b,
        datum_b=datum_b,
        axis=axis or [0, 0, 1, 0],
        limits=limits,
        dynamics=dynamics
    )


def create_prismatic_mate(
    name: str,
    part_a: str,
    datum_a: str,
    part_b: str,
    datum_b: str,
    axis: Optional[List[float]] = None,
    min_position: Optional[float] = None,
    max_position: Optional[float] = None,
    max_velocity: Optional[float] = None,
    max_force: Optional[float] = None,
    friction: float = 0.0,
    damping: float = 0.0
) -> Mate:
    """
    Convenience function to create a prismatic (slider) joint mate.

    A prismatic joint allows translation along a single axis with no rotation.
    Common for linear actuators, pistons, and drawer slides.

    Args:
        name: Descriptive name (e.g., "piston_stroke", "drawer_slide")
        part_a: Parent/base component identifier
        datum_a: Axis datum on parent component
        part_b: Child/moving component identifier
        datum_b: Axis datum on child component
        axis: Translation axis direction [x,y,z,w], defaults to [0,0,1,0] (Z-axis)
        min_position: Minimum position in mm (None = unlimited)
        max_position: Maximum position in mm (None = unlimited)
        max_velocity: Maximum velocity in mm/s (None = unlimited)
        max_force: Maximum force in N (None = unlimited)
        friction: Static friction coefficient (0.0 = frictionless)
        damping: Viscous damping coefficient (0.0 = no damping)

    Returns:
        Mate configured as prismatic joint

    Example:
        # Linear actuator with 100mm stroke
        actuator = create_prismatic_mate(
            name="z_axis_slide",
            part_a="base",
            datum_a="rail_axis",
            part_b="carriage",
            datum_b="slider_axis",
            min_position=0.0,
            max_position=100.0,
            max_velocity=50.0,
            friction=0.05,
            damping=5.0
        )
    """
    limits = None
    if any(x is not None for x in [min_position, max_position, max_velocity, max_force]):
        limits = MateLimits(
            min_value=min_position,
            max_value=max_position,
            max_velocity=max_velocity,
            max_effort=max_force
        )

    dynamics = None
    if friction > 0.0 or damping > 0.0:
        dynamics = MateDynamics(
            friction_static=friction,
            friction_kinetic=friction * 0.8,
            damping=damping
        )

    return Mate(
        name=name,
        mate_type=MateType.PRISMATIC,
        part_a=part_a,
        datum_a=datum_a,
        part_b=part_b,
        datum_b=datum_b,
        axis=axis or [0, 0, 1, 0],
        limits=limits,
        dynamics=dynamics
    )


def create_gear_mate(
    name: str,
    part_a: str,
    datum_a: str,
    part_b: str,
    datum_b: str,
    ratio: float,
    reverse: bool = False
) -> Mate:
    """
    Convenience function to create a gear coupling mate.

    Couples rotation of two components with a fixed ratio. The ratio is
    defined as: ratio = driven_rotations / driver_rotations

    Args:
        name: Descriptive name (e.g., "gear_1_to_2", "pulley_coupling")
        part_a: Driver (input) component identifier
        datum_a: Rotation axis on driver
        part_b: Driven (output) component identifier
        datum_b: Rotation axis on driven
        ratio: Motion ratio (driven_rotations / driver_rotations)
        reverse: If True, gears rotate in opposite directions

    Returns:
        Mate configured as gear coupling

    Example:
        # 3:1 reduction gearbox (output rotates 1/3 speed of input)
        gearbox = create_gear_mate(
            name="motor_to_output",
            part_a="motor_shaft",
            datum_a="motor_axis",
            part_b="output_shaft",
            datum_b="output_axis",
            ratio=1.0/3.0,
            reverse=False
        )
    """
    return Mate(
        name=name,
        mate_type=MateType.GEAR,
        part_a=part_a,
        datum_a=datum_a,
        part_b=part_b,
        datum_b=datum_b,
        coupling_ratio=ratio,
        coupling_reverse=reverse
    )


# =============================================================================
# COINCIDENT Constraint Evaluation Functions
# =============================================================================

@dataclass
class CoincidentResult:
    """Result of evaluating a COINCIDENT constraint.

    Attributes:
        satisfied: True if constraint is satisfied within tolerance
        error_distance: Distance error in mm (for point/origin alignment)
        error_angle: Angular error in degrees (for plane/axis alignment)
        error_message: Human-readable description of the result
        details: Additional information about the evaluation
    """
    satisfied: bool
    error_distance: float = 0.0
    error_angle: float = 0.0
    error_message: str = ""
    details: Dict[str, Any] = field(default_factory=dict)

    def __str__(self) -> str:
        """Format result as a readable string."""
        status = "PASS" if self.satisfied else "FAIL"
        errors = []
        if self.error_distance > 0:
            errors.append(f"distance: {self.error_distance:.3f}mm")
        if self.error_angle > 0:
            errors.append(f"angle: {self.error_angle:.2f}deg")
        error_str = ", ".join(errors) if errors else "none"
        return f"[{status}] {self.error_message} (error: {error_str})"


def evaluate_coincident(
    datum_a: 'Datum',
    datum_b: 'Datum',
    tolerance_mm: float = 0.1,
    tolerance_deg: float = 1.0
) -> CoincidentResult:
    """Evaluate a COINCIDENT constraint between two datums.

    COINCIDENT means two geometric features occupy the same location/orientation:
    - Point-to-Point: Origins coincide
    - Point-to-Plane: Point lies on plane
    - Plane-to-Plane: Planes are coplanar (same origin + parallel normals)
    - Circle-to-Circle: Centers coincide (for bolt circle alignment)
    - Axis-to-Axis: Axes are colinear

    Args:
        datum_a: First datum in world coordinates
        datum_b: Second datum in world coordinates
        tolerance_mm: Linear tolerance in millimeters (default: 0.1mm)
        tolerance_deg: Angular tolerance in degrees (default: 1.0 deg)

    Returns:
        CoincidentResult with evaluation status and error metrics

    Example:
        >>> # Check if bolt holes align
        >>> result = evaluate_coincident(
        ...     horn_bolt_circle,
        ...     sun_gear_bolt_circle,
        ...     tolerance_mm=0.05
        ... )
        >>> if result.satisfied:
        ...     print("Bolt holes align!")
    """
    if not HAS_NUMPY:
        return CoincidentResult(
            satisfied=False,
            error_message="NumPy required for constraint evaluation"
        )

    from .datum import DatumType

    # Get origins as numpy arrays
    origin_a = np.array(datum_a.origin[:3])
    origin_b = np.array(datum_b.origin[:3])

    # Calculate distance between origins
    distance = float(np.linalg.norm(origin_b - origin_a))

    # Handle different datum type combinations
    type_a = datum_a.datum_type
    type_b = datum_b.datum_type

    # Point-to-Point coincidence
    if type_a == DatumType.POINT and type_b == DatumType.POINT:
        satisfied = distance <= tolerance_mm
        return CoincidentResult(
            satisfied=satisfied,
            error_distance=distance,
            error_message=f"Point-to-Point: distance={distance:.4f}mm",
            details={
                "origin_a": origin_a.tolist(),
                "origin_b": origin_b.tolist(),
                "type": "point_to_point"
            }
        )

    # Circle-to-Circle coincidence (bolt circles)
    if type_a == DatumType.CIRCLE and type_b == DatumType.CIRCLE:
        # Check center alignment
        satisfied_distance = distance <= tolerance_mm

        # Also check normal alignment for coplanar circles
        normal_a = np.array(datum_a.normal[:3])
        normal_b = np.array(datum_b.normal[:3])
        normal_a = normal_a / np.linalg.norm(normal_a)
        normal_b = normal_b / np.linalg.norm(normal_b)

        # Normals can be parallel or antiparallel for coplanar
        dot = abs(float(np.dot(normal_a, normal_b)))
        angle_error = math.degrees(math.acos(min(1.0, dot)))
        satisfied_angle = angle_error <= tolerance_deg

        satisfied = satisfied_distance and satisfied_angle

        # Check radius match (informational)
        radius_diff = abs(datum_a.radius - datum_b.radius) if (
            datum_a.radius and datum_b.radius
        ) else 0.0

        return CoincidentResult(
            satisfied=satisfied,
            error_distance=distance,
            error_angle=angle_error,
            error_message=f"Circle-to-Circle: center_dist={distance:.4f}mm, "
                         f"angle_err={angle_error:.2f}deg, "
                         f"radius_diff={radius_diff:.3f}mm",
            details={
                "center_a": origin_a.tolist(),
                "center_b": origin_b.tolist(),
                "normal_a": normal_a.tolist(),
                "normal_b": normal_b.tolist(),
                "radius_a": datum_a.radius,
                "radius_b": datum_b.radius,
                "radius_difference": radius_diff,
                "type": "circle_to_circle"
            }
        )

    # Plane-to-Plane coincidence (coplanar)
    if type_a == DatumType.PLANE and type_b == DatumType.PLANE:
        normal_a = np.array(datum_a.normal[:3])
        normal_b = np.array(datum_b.normal[:3])
        normal_a = normal_a / np.linalg.norm(normal_a)
        normal_b = normal_b / np.linalg.norm(normal_b)

        # Check parallel normals (can be same or opposite direction)
        dot = abs(float(np.dot(normal_a, normal_b)))
        angle_error = math.degrees(math.acos(min(1.0, dot)))
        satisfied_angle = angle_error <= tolerance_deg

        # Check that origin_b lies on plane_a (distance to plane)
        # Distance from point to plane: |n . (p - p0)|
        plane_distance = abs(float(np.dot(normal_a, origin_b - origin_a)))
        satisfied_distance = plane_distance <= tolerance_mm

        satisfied = satisfied_distance and satisfied_angle

        return CoincidentResult(
            satisfied=satisfied,
            error_distance=plane_distance,
            error_angle=angle_error,
            error_message=f"Plane-to-Plane: plane_dist={plane_distance:.4f}mm, "
                         f"angle_err={angle_error:.2f}deg",
            details={
                "origin_a": origin_a.tolist(),
                "origin_b": origin_b.tolist(),
                "normal_a": normal_a.tolist(),
                "normal_b": normal_b.tolist(),
                "plane_distance": plane_distance,
                "type": "plane_to_plane"
            }
        )

    # Axis-to-Axis coincidence (colinear)
    if type_a == DatumType.AXIS and type_b == DatumType.AXIS:
        dir_a = np.array(datum_a.direction[:3])
        dir_b = np.array(datum_b.direction[:3])
        dir_a = dir_a / np.linalg.norm(dir_a)
        dir_b = dir_b / np.linalg.norm(dir_b)

        # Check parallel directions
        dot = abs(float(np.dot(dir_a, dir_b)))
        angle_error = math.degrees(math.acos(min(1.0, dot)))
        satisfied_angle = angle_error <= tolerance_deg

        # Check that origin_b lies on axis_a (perpendicular distance)
        # Vector from origin_a to origin_b
        v = origin_b - origin_a
        # Component along axis
        along = np.dot(v, dir_a) * dir_a
        # Perpendicular component
        perp = v - along
        perp_distance = float(np.linalg.norm(perp))
        satisfied_distance = perp_distance <= tolerance_mm

        satisfied = satisfied_distance and satisfied_angle

        return CoincidentResult(
            satisfied=satisfied,
            error_distance=perp_distance,
            error_angle=angle_error,
            error_message=f"Axis-to-Axis: perp_dist={perp_distance:.4f}mm, "
                         f"angle_err={angle_error:.2f}deg",
            details={
                "origin_a": origin_a.tolist(),
                "origin_b": origin_b.tolist(),
                "direction_a": dir_a.tolist(),
                "direction_b": dir_b.tolist(),
                "perpendicular_distance": perp_distance,
                "type": "axis_to_axis"
            }
        )

    # Point-to-Plane coincidence (point lies on plane)
    if (type_a == DatumType.POINT and type_b == DatumType.PLANE):
        normal = np.array(datum_b.normal[:3])
        normal = normal / np.linalg.norm(normal)
        plane_origin = origin_b
        point = origin_a

        # Distance from point to plane
        plane_distance = abs(float(np.dot(normal, point - plane_origin)))
        satisfied = plane_distance <= tolerance_mm

        return CoincidentResult(
            satisfied=satisfied,
            error_distance=plane_distance,
            error_message=f"Point-to-Plane: distance={plane_distance:.4f}mm",
            details={
                "point": point.tolist(),
                "plane_origin": plane_origin.tolist(),
                "plane_normal": normal.tolist(),
                "type": "point_to_plane"
            }
        )

    if (type_a == DatumType.PLANE and type_b == DatumType.POINT):
        normal = np.array(datum_a.normal[:3])
        normal = normal / np.linalg.norm(normal)
        plane_origin = origin_a
        point = origin_b

        # Distance from point to plane
        plane_distance = abs(float(np.dot(normal, point - plane_origin)))
        satisfied = plane_distance <= tolerance_mm

        return CoincidentResult(
            satisfied=satisfied,
            error_distance=plane_distance,
            error_message=f"Plane-to-Point: distance={plane_distance:.4f}mm",
            details={
                "point": point.tolist(),
                "plane_origin": plane_origin.tolist(),
                "plane_normal": normal.tolist(),
                "type": "plane_to_point"
            }
        )

    # Default: just check origin coincidence
    satisfied = distance <= tolerance_mm
    return CoincidentResult(
        satisfied=satisfied,
        error_distance=distance,
        error_message=f"{type_a.value}-to-{type_b.value}: distance={distance:.4f}mm",
        details={
            "origin_a": origin_a.tolist(),
            "origin_b": origin_b.tolist(),
            "type": f"{type_a.value}_to_{type_b.value}"
        }
    )


def check_bolt_circle_alignment(
    bolt_circle_a: 'Datum',
    bolt_circle_b: 'Datum',
    hole_count: int,
    tolerance_mm: float = 0.1,
    angular_offset_deg: float = 0.0
) -> CoincidentResult:
    """Check if two bolt circles align for proper mating.

    This validates that mounting hole patterns match for servo horn to
    gear hub interfaces, motor mounts, etc.

    Args:
        bolt_circle_a: First bolt circle datum (e.g., servo horn)
        bolt_circle_b: Second bolt circle datum (e.g., sun gear hub)
        hole_count: Number of holes in each pattern
        tolerance_mm: Position tolerance for hole centers
        angular_offset_deg: Expected angular offset between patterns (e.g., 45 degrees)

    Returns:
        CoincidentResult with alignment status

    Example:
        >>> # XH540 servo horn (4 holes at 45 deg offset) to sun gear
        >>> result = check_bolt_circle_alignment(
        ...     horn_bolt_circle,
        ...     sun_gear_bolt_circle,
        ...     hole_count=4,
        ...     angular_offset_deg=45.0
        ... )
    """
    if not HAS_NUMPY:
        return CoincidentResult(
            satisfied=False,
            error_message="NumPy required for constraint evaluation"
        )

    from .datum import DatumType

    # Verify both are circle datums
    if bolt_circle_a.datum_type != DatumType.CIRCLE:
        return CoincidentResult(
            satisfied=False,
            error_message=f"First datum must be CIRCLE, got {bolt_circle_a.datum_type.value}"
        )
    if bolt_circle_b.datum_type != DatumType.CIRCLE:
        return CoincidentResult(
            satisfied=False,
            error_message=f"Second datum must be CIRCLE, got {bolt_circle_b.datum_type.value}"
        )

    # First check basic coincidence (centers align, normals parallel)
    basic_result = evaluate_coincident(bolt_circle_a, bolt_circle_b, tolerance_mm)

    if not basic_result.satisfied:
        return basic_result

    # Check radius match
    radius_a = bolt_circle_a.radius
    radius_b = bolt_circle_b.radius
    radius_diff = abs(radius_a - radius_b)

    if radius_diff > tolerance_mm:
        return CoincidentResult(
            satisfied=False,
            error_distance=radius_diff,
            error_message=f"Bolt circle radii don't match: {radius_a:.3f}mm vs {radius_b:.3f}mm",
            details={
                "radius_a": radius_a,
                "radius_b": radius_b,
                "radius_difference": radius_diff,
                "type": "radius_mismatch"
            }
        )

    # Calculate hole positions for both circles
    center_a = np.array(bolt_circle_a.origin[:3])
    center_b = np.array(bolt_circle_b.origin[:3])
    normal_a = np.array(bolt_circle_a.normal[:3])
    normal_a = normal_a / np.linalg.norm(normal_a)

    # Generate basis vectors in the plane of the bolt circle
    # Find a vector perpendicular to the normal
    if abs(normal_a[0]) < 0.9:
        ref = np.array([1, 0, 0])
    else:
        ref = np.array([0, 1, 0])

    u = np.cross(normal_a, ref)
    u = u / np.linalg.norm(u)
    v = np.cross(normal_a, u)
    v = v / np.linalg.norm(v)

    # Calculate hole positions
    hole_spacing = 360.0 / hole_count
    offset_rad = math.radians(angular_offset_deg)

    max_hole_error = 0.0
    hole_errors = []

    for i in range(hole_count):
        # Hole position in circle A
        angle_a = math.radians(i * hole_spacing)
        pos_a = center_a + radius_a * (math.cos(angle_a) * u + math.sin(angle_a) * v)

        # Corresponding hole position in circle B (with angular offset)
        angle_b = angle_a + offset_rad
        pos_b = center_b + radius_b * (math.cos(angle_b) * u + math.sin(angle_b) * v)

        # Calculate error
        error = float(np.linalg.norm(pos_b - pos_a))
        hole_errors.append(error)
        max_hole_error = max(max_hole_error, error)

    satisfied = max_hole_error <= tolerance_mm

    return CoincidentResult(
        satisfied=satisfied,
        error_distance=max_hole_error,
        error_message=f"Bolt circle alignment: max_error={max_hole_error:.4f}mm "
                     f"({hole_count} holes, {angular_offset_deg}deg offset)",
        details={
            "center_a": center_a.tolist(),
            "center_b": center_b.tolist(),
            "radius_a": radius_a,
            "radius_b": radius_b,
            "hole_count": hole_count,
            "angular_offset_deg": angular_offset_deg,
            "hole_errors": hole_errors,
            "max_error": max_hole_error,
            "type": "bolt_circle_alignment"
        }
    )


def create_coincident_mate(
    name: str,
    part_a: str,
    datum_a: str,
    part_b: str,
    datum_b: str,
    offset: float = 0.0
) -> Mate:
    """Convenience function to create a COINCIDENT mate.

    A COINCIDENT mate constrains two features to occupy the same location:
    - Points share the same position
    - Planes are coplanar
    - Bolt circles align (centers coincide)

    Args:
        name: Descriptive name (e.g., "horn_to_sun_gear")
        part_a: First component identifier
        datum_a: Datum on first component
        part_b: Second component identifier
        datum_b: Datum on second component
        offset: Optional offset distance along normal (default: 0)

    Returns:
        Mate configured as COINCIDENT constraint

    Example:
        >>> # Align servo horn bolt holes with sun gear bolt holes
        >>> mate = create_coincident_mate(
        ...     name="horn_bolt_alignment",
        ...     part_a="xh540_servo",
        ...     datum_a="horn_bolt_circle",
        ...     part_b="axis1_sun_gear",
        ...     datum_b="hub_bolt_circle"
        ... )
    """
    return Mate(
        name=name,
        mate_type=MateType.COINCIDENT,
        part_a=part_a,
        datum_a=datum_a,
        part_b=part_b,
        datum_b=datum_b,
        offset=offset
    )
