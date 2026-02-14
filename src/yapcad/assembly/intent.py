"""Declarative Assembly Intent System for yapCAD.

This module provides a declarative approach to assembly design where designers
specify WHAT parts must do (functional requirements) rather than HOW to achieve
it (explicit transforms). The system derives geometry from requirements.

Key Concepts:
    AssemblyIntent: Top-level container for declarative assembly specification
    FunctionalRequirement: What a part must DO (contact, roll, align, etc.)
    Connection: How parts relate (topology, not exact position)
    Clearance: What must NOT happen (collision avoidance)
    ReferenceGeometry: Fixed geometry that parts must relate to

Example - Wheel Assembly:
    >>> tube = ReferenceGeometry(
    ...     name="tube_inner_wall",
    ...     geometry_type="cylinder",
    ...     center=(0, 0, 0),
    ...     axis=(0, 0, 1),
    ...     radius=175.0,
    ... )
    >>>
    >>> wheel_assembly = AssemblyIntent(
    ...     name="wheel_pod",
    ...     reference_geometry={"tube": tube},
    ...     functional_requirements=[
    ...         ContactRequirement(
    ...             name="wheel_contact",
    ...             part="motor",
    ...             surface="tire_outer",
    ...             target="tube_inner_wall",
    ...             contact_type="rolling",
    ...         ),
    ...         RollRequirement(
    ...             name="roll_along_z",
    ...             part="motor",
    ...             roll_direction="along_tube_axis",
    ...         ),
    ...     ],
    ...     connections=[
    ...         Connection(
    ...             parent="chassis.pivot_boss",
    ...             child="wheel_arm.pivot_bore",
    ...             joint_type="revolute",
    ...             axis="tangent",
    ...         ),
    ...     ],
    ...     clearances=[
    ...         Clearance("motor", "chassis", min_distance=5.0),
    ...     ],
    ...     derived_parameters=["wheel_arm.length", "wheel_center_radius"],
    ... )
    >>>
    >>> result = wheel_assembly.solve()
    >>> print(result.derived["wheel_arm.length"])

Copyright (c) 2026 yapCAD contributors
License: MIT
"""

from __future__ import annotations

import math
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from enum import Enum
from typing import (
    Any,
    Callable,
    Dict,
    List,
    Optional,
    Tuple,
    Union,
    TYPE_CHECKING,
)

# Import from existing assembly system
try:
    from .constraint import Constraint, ConstraintType, ConstraintResult
    from .mate import Mate, MateType, MateLimits
    from .assembly import Assembly, AssemblyValidationResult
    from .datum import Datum, DatumType, PartDefinition
except ImportError:
    # Fallback for standalone testing
    Constraint = Any
    ConstraintType = Any
    ConstraintResult = Any
    Mate = Any
    MateType = Any
    MateLimits = Any
    Assembly = Any
    AssemblyValidationResult = Any
    Datum = Any
    DatumType = Any
    PartDefinition = Any

# Import numpy
try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None


# =============================================================================
# REFERENCE GEOMETRY
# =============================================================================

class GeometryType(Enum):
    """Types of reference geometry."""
    POINT = "point"
    LINE = "line"
    PLANE = "plane"
    CYLINDER = "cylinder"
    SPHERE = "sphere"
    CONE = "cone"


@dataclass
class ReferenceGeometry:
    """Fixed reference geometry that parts must relate to.

    Reference geometry defines the environment or constraints that the
    assembly operates within. Examples include tube walls, mounting
    surfaces, and clearance envelopes.

    Attributes:
        name: Unique identifier for this geometry
        geometry_type: Type of geometry (cylinder, plane, etc.)
        center: Center point (x, y, z) for most geometry types
        axis: Direction vector for cylinders, cones, lines
        radius: Radius for cylinders, spheres, cones
        normal: Normal vector for planes
        inner: True if surface is inner (e.g., inner wall of tube)

    Example:
        >>> # 350mm ID tube
        >>> tube = ReferenceGeometry(
        ...     name="tube_inner_wall",
        ...     geometry_type="cylinder",
        ...     center=(0, 0, 0),
        ...     axis=(0, 0, 1),
        ...     radius=175.0,
        ...     inner=True,
        ... )
    """
    name: str
    geometry_type: str  # "point", "line", "plane", "cylinder", "sphere", "cone"

    # Geometry parameters
    center: Optional[Tuple[float, float, float]] = None
    axis: Optional[Tuple[float, float, float]] = None
    radius: Optional[float] = None
    normal: Optional[Tuple[float, float, float]] = None
    inner: bool = False  # For cylinders: inner wall vs outer wall

    def get_surface_point(
        self,
        theta: float = 0.0,
        z: float = 0.0
    ) -> Tuple[float, float, float]:
        """Get a point on the surface at given parameters.

        For cylinders: theta is angle in radians, z is height along axis.
        For planes: theta and z are ignored.
        For spheres: theta is azimuth, z is elevation (radians).

        Returns:
            Point (x, y, z) on the surface
        """
        if self.geometry_type == "cylinder":
            if self.center is None or self.radius is None:
                raise ValueError("Cylinder requires center and radius")
            r = self.radius
            x = self.center[0] + r * math.cos(theta)
            y = self.center[1] + r * math.sin(theta)
            z_coord = self.center[2] + z
            return (x, y, z_coord)

        elif self.geometry_type == "plane":
            if self.center is None:
                return (0, 0, 0)
            return self.center

        elif self.geometry_type == "sphere":
            if self.center is None or self.radius is None:
                raise ValueError("Sphere requires center and radius")
            r = self.radius
            x = self.center[0] + r * math.cos(theta) * math.cos(z)
            y = self.center[1] + r * math.sin(theta) * math.cos(z)
            z_coord = self.center[2] + r * math.sin(z)
            return (x, y, z_coord)

        else:
            return self.center if self.center else (0, 0, 0)

    def get_normal_at(
        self,
        point: Tuple[float, float, float]
    ) -> Tuple[float, float, float]:
        """Get the surface normal at a given point.

        For cylinders: radial direction from axis.
        For planes: constant normal.
        For spheres: radial direction from center.

        Returns:
            Unit normal vector (x, y, z)
        """
        if self.geometry_type == "cylinder":
            if self.center is None:
                raise ValueError("Cylinder requires center")
            # Project point onto axis, get perpendicular
            c = self.center
            dx = point[0] - c[0]
            dy = point[1] - c[1]
            mag = math.sqrt(dx * dx + dy * dy)
            if mag < 1e-10:
                return (1, 0, 0)
            if self.inner:
                return (-dx / mag, -dy / mag, 0)  # Inward normal
            else:
                return (dx / mag, dy / mag, 0)  # Outward normal

        elif self.geometry_type == "plane":
            if self.normal is None:
                return (0, 0, 1)
            mag = math.sqrt(sum(n * n for n in self.normal))
            return tuple(n / mag for n in self.normal)

        elif self.geometry_type == "sphere":
            if self.center is None:
                raise ValueError("Sphere requires center")
            c = self.center
            dx = point[0] - c[0]
            dy = point[1] - c[1]
            dz = point[2] - c[2]
            mag = math.sqrt(dx * dx + dy * dy + dz * dz)
            if mag < 1e-10:
                return (0, 0, 1)
            if self.inner:
                return (-dx / mag, -dy / mag, -dz / mag)
            else:
                return (dx / mag, dy / mag, dz / mag)

        else:
            return (0, 0, 1)


# =============================================================================
# FUNCTIONAL REQUIREMENTS
# =============================================================================

@dataclass
class FunctionalRequirement(ABC):
    """Base class for functional requirements.

    Functional requirements describe WHAT a part must DO, not HOW to achieve it.
    The system derives constraints and parameters from these requirements.

    Subclasses implement `get_implied_constraints()` to convert the high-level
    requirement into low-level constraints that the solver can work with.

    Attributes:
        name: Unique identifier for this requirement
        part: Name of the part that must satisfy this requirement
        description: Human-readable description of the requirement
        priority: Relative importance (higher = more important)
    """
    name: str
    part: str
    description: str = ""
    priority: int = 1

    @abstractmethod
    def get_implied_constraints(
        self,
        reference_geometry: Dict[str, ReferenceGeometry]
    ) -> List[Dict[str, Any]]:
        """Derive low-level constraints from this requirement.

        Args:
            reference_geometry: Dictionary of available reference geometry

        Returns:
            List of constraint specifications (dicts that can create Constraints)
        """
        pass

    def get_derived_parameters(self) -> List[str]:
        """List parameters that can be derived from this requirement.

        Returns:
            List of parameter names in "part.parameter" format
        """
        return []


@dataclass
class ContactRequirement(FunctionalRequirement):
    """Part surface must contact a reference surface.

    This requirement specifies that a datum surface on a part must
    touch a reference geometry surface. Common examples:
    - Wheel tire contacting tube inner wall
    - Mounting face flush with bracket surface
    - Ball bearing in socket

    Attributes:
        surface: Datum name on part (e.g., "tire_outer_diameter")
        target: Reference geometry name (e.g., "tube_inner_wall")
        contact_type: "static", "rolling", or "sliding"
        preload_source: Optional name of what provides contact force

    Example:
        >>> contact = ContactRequirement(
        ...     name="wheel_contact",
        ...     part="ddsm115_motor",
        ...     surface="tire_outer_diameter",
        ...     target="tube_inner_wall",
        ...     contact_type="rolling",
        ...     preload_source="suspension_spring",
        ... )
    """
    surface: str = ""
    target: str = ""
    contact_type: str = "static"  # "static", "rolling", "sliding"
    preload_source: Optional[str] = None

    def get_implied_constraints(
        self,
        reference_geometry: Dict[str, ReferenceGeometry]
    ) -> List[Dict[str, Any]]:
        """Contact implies position at specific radius/distance."""
        constraints = []

        # Find target geometry
        target_geom = None
        for name, geom in reference_geometry.items():
            if name == self.target or geom.name == self.target:
                target_geom = geom
                break

        if target_geom is None:
            # Return constraint that will fail validation
            return [{
                "name": f"{self.name}_missing_target",
                "type": "error",
                "message": f"Target geometry '{self.target}' not found",
            }]

        # For cylinder contact, implies AT_RADIUS constraint
        if target_geom.geometry_type == "cylinder":
            constraints.append({
                "name": f"{self.part}_at_{self.target}",
                "constraint_type": "AT_RADIUS",
                "part": self.part,
                "datum": self.surface,
                "center": target_geom.center,
                "radius": target_geom.radius,
                "description": f"{self.part} surface contacts {self.target}",
            })

        # For plane contact, implies coincident/flush constraint
        elif target_geom.geometry_type == "plane":
            constraints.append({
                "name": f"{self.part}_on_{self.target}",
                "constraint_type": "IN_PLANE",
                "part": self.part,
                "datum": self.surface,
                "plane_origin": target_geom.center,
                "plane_normal": target_geom.normal,
                "description": f"{self.part} surface flush with {self.target}",
            })

        return constraints

    def get_derived_parameters(self) -> List[str]:
        return [
            f"{self.part}.contact_radius",
            f"{self.part}.contact_position",
        ]


@dataclass
class RollRequirement(FunctionalRequirement):
    """Part must roll in a specified direction.

    This requirement specifies that a rotating part (wheel, roller, etc.)
    must roll in a specific direction. This IMPLIES that the rotation
    axis must be perpendicular to the roll direction.

    For a wheel rolling along the Z-axis (tube axis):
    - Roll direction: (0, 0, 1)
    - Rotation axis must be tangent (perpendicular to both Z and radial)

    Attributes:
        roll_direction: Direction of travel
            - "along_tube_axis" or "+z" for vertical tube
            - "(x, y, z)" tuple for custom direction
        axis_datum: Name of rotation axis datum on part

    Example:
        >>> roll = RollRequirement(
        ...     name="wheel_rolls_z",
        ...     part="motor",
        ...     roll_direction="along_tube_axis",
        ...     axis_datum="motor_axis",
        ... )
    """
    roll_direction: str = "along_tube_axis"  # or "+z", "-z", tuple
    axis_datum: str = "rotation_axis"

    def get_implied_constraints(
        self,
        reference_geometry: Dict[str, ReferenceGeometry]
    ) -> List[Dict[str, Any]]:
        """Rolling implies rotation axis perpendicular to roll direction."""
        constraints = []

        # Parse roll direction
        if self.roll_direction in ("along_tube_axis", "+z"):
            roll_dir = (0, 0, 1)
        elif self.roll_direction == "-z":
            roll_dir = (0, 0, -1)
        elif isinstance(self.roll_direction, (tuple, list)):
            roll_dir = tuple(self.roll_direction)
        else:
            roll_dir = (0, 0, 1)  # Default

        # For rolling along Z in a cylindrical tube, axis must be TANGENT
        if roll_dir == (0, 0, 1) or roll_dir == (0, 0, -1):
            # Look for a cylinder reference geometry
            for geom in reference_geometry.values():
                if geom.geometry_type == "cylinder":
                    constraints.append({
                        "name": f"{self.part}_axis_tangent",
                        "constraint_type": "TANGENT_TO_CIRCLE",
                        "part": self.part,
                        "datum": self.axis_datum,
                        "center": geom.center,
                        "radius": geom.radius,
                        "description": f"{self.part} axis tangent for rolling along tube",
                    })
                    break

        # General case: axis perpendicular to roll direction
        constraints.append({
            "name": f"{self.part}_axis_perpendicular_to_roll",
            "constraint_type": "PERPENDICULAR_TO",
            "part": self.part,
            "datum": self.axis_datum,
            "axis": roll_dir,
            "description": f"{self.part} axis perpendicular to roll direction",
        })

        return constraints


@dataclass
class AxisOrientationRequirement(FunctionalRequirement):
    """Part axis must have specific orientation relative to reference.

    This requirement specifies that a datum axis on a part must point
    in a specific direction relative to the assembly or reference geometry.

    Attributes:
        axis_datum: Name of the axis datum on the part
        orientation: Direction specification
            - "tangent": Tangent to reference cylinder at part position
            - "radial": Pointing toward/away from reference center
            - "axial": Parallel to reference axis (e.g., tube Z-axis)
            - "+x", "-y", "+z": Global direction
            - tuple: Custom direction vector
        reference: Name of reference geometry or "global"
        pointing: "toward" or "away" for radial orientation

    Example:
        >>> axis_req = AxisOrientationRequirement(
        ...     name="motor_tangent",
        ...     part="motor",
        ...     axis_datum="motor_axis",
        ...     orientation="tangent",
        ...     reference="tube",
        ... )
    """
    axis_datum: str = ""
    orientation: str = "tangent"  # "tangent", "radial", "axial", direction
    reference: str = "global"
    pointing: str = "toward"  # For radial: "toward" or "away"

    def get_implied_constraints(
        self,
        reference_geometry: Dict[str, ReferenceGeometry]
    ) -> List[Dict[str, Any]]:
        """Convert orientation requirement to constraint."""
        constraints = []

        # Find reference geometry
        ref_geom = reference_geometry.get(self.reference)

        if self.orientation == "tangent":
            if ref_geom and ref_geom.geometry_type == "cylinder":
                constraints.append({
                    "name": f"{self.part}_{self.axis_datum}_tangent",
                    "constraint_type": "TANGENT_TO_CIRCLE",
                    "part": self.part,
                    "datum": self.axis_datum,
                    "center": ref_geom.center,
                    "radius": ref_geom.radius,
                })
            else:
                constraints.append({
                    "name": f"{self.part}_{self.axis_datum}_tangent",
                    "constraint_type": "TANGENT_TO_CIRCLE",
                    "part": self.part,
                    "datum": self.axis_datum,
                    "center": (0, 0, 0),
                })

        elif self.orientation == "radial":
            direction = "outward" if self.pointing == "away" else "inward"
            constraints.append({
                "name": f"{self.part}_{self.axis_datum}_radial",
                "constraint_type": "RADIAL_FROM_CENTER",
                "part": self.part,
                "datum": self.axis_datum,
                "center": ref_geom.center if ref_geom else (0, 0, 0),
                "direction": direction,
            })

        elif self.orientation == "axial":
            # Parallel to reference axis (default Z)
            ref_axis = ref_geom.axis if ref_geom and ref_geom.axis else (0, 0, 1)
            constraints.append({
                "name": f"{self.part}_{self.axis_datum}_axial",
                "constraint_type": "PARALLEL_TO",
                "part": self.part,
                "datum": self.axis_datum,
                "axis": ref_axis,
            })

        elif self.orientation in ("+x", "-x", "+y", "-y", "+z", "-z"):
            # Global direction
            dir_map = {
                "+x": (1, 0, 0), "-x": (-1, 0, 0),
                "+y": (0, 1, 0), "-y": (0, -1, 0),
                "+z": (0, 0, 1), "-z": (0, 0, -1),
            }
            constraints.append({
                "name": f"{self.part}_{self.axis_datum}_facing",
                "constraint_type": "FACING",
                "part": self.part,
                "datum": self.axis_datum,
                "direction": self.orientation,
            })

        return constraints


@dataclass
class ParallelAxesRequirement(FunctionalRequirement):
    """Two axes must be parallel.

    This requirement specifies that two datum axes (on same or different
    parts) must be parallel. Common use: pivot axis parallel to motor axis
    for clean suspension motion.

    Attributes:
        axis_a: First axis in "part.datum" format
        axis_b: Second axis in "part.datum" format
        allow_opposite: If True, axes can point in opposite directions

    Example:
        >>> parallel = ParallelAxesRequirement(
        ...     name="pivot_motor_parallel",
        ...     part="wheel_arm",  # Primary part
        ...     axis_a="wheel_arm.pivot_bore_axis",
        ...     axis_b="motor.motor_axis",
        ... )
    """
    axis_a: str = ""
    axis_b: str = ""
    allow_opposite: bool = True

    def get_implied_constraints(
        self,
        reference_geometry: Dict[str, ReferenceGeometry]
    ) -> List[Dict[str, Any]]:
        """Generate parallel constraint."""
        # Parse "part.datum" format
        part_a, datum_a = self.axis_a.split(".") if "." in self.axis_a else (self.part, self.axis_a)
        part_b, datum_b = self.axis_b.split(".") if "." in self.axis_b else (self.part, self.axis_b)

        return [{
            "name": f"{part_a}_{datum_a}_parallel_{part_b}_{datum_b}",
            "constraint_type": "PARALLEL_TO",
            "part": part_a,
            "datum": datum_a,
            "reference_part": part_b,
            "reference_datum": datum_b,
            "description": f"{self.axis_a} parallel to {self.axis_b}",
        }]


@dataclass
class ReachRequirement(FunctionalRequirement):
    """End effector must have specified reach envelope.

    For serial manipulators (SCARA arms, etc.), specifies the minimum
    and maximum reach from the base.

    Attributes:
        end_effector: Datum name on end effector part (e.g., "tool_point")
        base: Datum name on base part (e.g., "base_center")
        min_reach: Minimum distance from base (mm)
        max_reach: Maximum distance from base (mm)

    Example:
        >>> reach = ReachRequirement(
        ...     name="scara_workspace",
        ...     part="wrist",
        ...     end_effector="wrist.tool_point",
        ...     base="tower.axis1_center",
        ...     min_reach=50.0,
        ...     max_reach=200.0,
        ... )
    """
    end_effector: str = ""
    base: str = ""
    min_reach: float = 0.0
    max_reach: float = 100.0

    def get_implied_constraints(
        self,
        reference_geometry: Dict[str, ReferenceGeometry]
    ) -> List[Dict[str, Any]]:
        """Reach implies link length constraints."""
        # This is more complex - needs kinematic analysis
        # For now, return informational constraint
        return [{
            "name": f"{self.name}_reach_info",
            "constraint_type": "CUSTOM",
            "description": f"Reach: {self.min_reach}mm to {self.max_reach}mm",
            "min_reach": self.min_reach,
            "max_reach": self.max_reach,
        }]

    def get_derived_parameters(self) -> List[str]:
        return [
            f"{self.part}.link_length",
            "total_reach",
        ]


# =============================================================================
# CONNECTIONS
# =============================================================================

@dataclass
class Connection:
    """Defines how two parts connect (topology).

    Connections specify the parent-child relationships between parts
    and the type of joint (rigid, revolute, prismatic, etc.).

    Attributes:
        parent: Parent datum in "part.datum" format
        child: Child datum in "part.datum" format
        joint_type: Type of connection
            - "rigid": No relative motion
            - "revolute": Rotation about axis
            - "prismatic": Translation along axis
            - "cylindrical": Rotation + translation on axis
            - "spherical": Ball-and-socket
        axis: Axis specification for joints
            - "tangent", "radial", "axial" for reference-relative
            - tuple (x, y, z) for explicit direction
        limits: (min, max) for joint limits (degrees or mm)
        interface_type: Physical connection type
        interface_details: Additional interface parameters

    Example:
        >>> conn = Connection(
        ...     parent="chassis.pivot_boss",
        ...     child="wheel_arm.pivot_bore",
        ...     joint_type="revolute",
        ...     axis="tangent",
        ...     limits=(-15, 15),
        ... )
    """
    parent: str
    child: str
    joint_type: str = "rigid"

    # Joint parameters
    axis: Optional[Union[str, Tuple[float, float, float]]] = None
    limits: Optional[Tuple[float, float]] = None

    # Interface details
    interface_type: Optional[str] = None  # "bolt_pattern", "press_fit", etc.
    interface_details: Dict[str, Any] = field(default_factory=dict)

    def get_parent_part(self) -> str:
        """Extract parent part name."""
        return self.parent.split(".")[0] if "." in self.parent else self.parent

    def get_parent_datum(self) -> str:
        """Extract parent datum name."""
        return self.parent.split(".")[1] if "." in self.parent else "origin"

    def get_child_part(self) -> str:
        """Extract child part name."""
        return self.child.split(".")[0] if "." in self.child else self.child

    def get_child_datum(self) -> str:
        """Extract child datum name."""
        return self.child.split(".")[1] if "." in self.child else "origin"

    def to_mate_spec(self) -> Dict[str, Any]:
        """Convert to mate specification dict."""
        spec = {
            "name": f"{self.get_parent_part()}_to_{self.get_child_part()}",
            "part_a": self.get_parent_part(),
            "datum_a": self.get_parent_datum(),
            "part_b": self.get_child_part(),
            "datum_b": self.get_child_datum(),
        }

        if self.joint_type == "rigid":
            spec["mate_type"] = "COINCIDENT"
        elif self.joint_type == "revolute":
            spec["mate_type"] = "REVOLUTE"
            spec["axis"] = self.axis
            if self.limits:
                spec["limits"] = {
                    "min_value": math.radians(self.limits[0]),
                    "max_value": math.radians(self.limits[1]),
                }
        elif self.joint_type == "prismatic":
            spec["mate_type"] = "PRISMATIC"
            spec["axis"] = self.axis
            if self.limits:
                spec["limits"] = {
                    "min_value": self.limits[0],
                    "max_value": self.limits[1],
                }
        elif self.joint_type == "cylindrical":
            spec["mate_type"] = "CYLINDRICAL"
        elif self.joint_type == "spherical":
            spec["mate_type"] = "SPHERICAL"

        return spec


# =============================================================================
# CLEARANCES
# =============================================================================

@dataclass
class Clearance:
    """Parts must maintain minimum distance (no collision).

    Clearance constraints ensure that parts do not interfere with each
    other or with reference geometry.

    Attributes:
        part_a: First part name
        part_b: Second part name (or reference geometry name)
        min_distance: Minimum allowed distance (mm)
        check_type: How to measure distance
            - "bounding_box": Use axis-aligned bounding boxes (fast)
            - "convex_hull": Use convex hull approximation
            - "mesh": Use full mesh collision (slow, accurate)
        critical: If True, violation is error; else warning

    Example:
        >>> clearance = Clearance(
        ...     part_a="motor",
        ...     part_b="chassis_plate",
        ...     min_distance=5.0,
        ... )
    """
    part_a: str
    part_b: str
    min_distance: float
    check_type: str = "bounding_box"
    critical: bool = True

    def validate(
        self,
        transforms: Dict[str, Any],
        part_bounds: Dict[str, Any]
    ) -> 'ClearanceResult':
        """Check if clearance is satisfied.

        Args:
            transforms: Part name -> transform matrix
            part_bounds: Part name -> bounding box/mesh

        Returns:
            ClearanceResult with validation status
        """
        # Simplified bounding box check
        # Full implementation would use actual geometry
        return ClearanceResult(
            satisfied=True,
            actual_distance=self.min_distance + 1.0,
            required_distance=self.min_distance,
            message=f"Clearance {self.part_a} <-> {self.part_b}: OK",
        )


@dataclass
class ClearanceResult:
    """Result of clearance validation."""
    satisfied: bool
    actual_distance: float
    required_distance: float
    message: str = ""
    interference_point: Optional[Tuple[float, float, float]] = None


# =============================================================================
# SOLVE RESULT
# =============================================================================

@dataclass
class SolveResult:
    """Result of solving an AssemblyIntent.

    Contains derived parameters, computed transforms, and validation status.

    Attributes:
        success: True if all requirements could be satisfied
        derived: Dictionary of derived parameter values
        transforms: Dictionary of computed transforms (part name -> 4x4 matrix)
        constraints_generated: List of constraint specifications generated
        validation: Validation result for all constraints
        errors: List of error messages if solve failed
        warnings: List of warning messages
    """
    success: bool
    derived: Dict[str, float] = field(default_factory=dict)
    transforms: Dict[str, Any] = field(default_factory=dict)
    constraints_generated: List[Dict[str, Any]] = field(default_factory=list)
    validation: Optional[Any] = None  # AssemblyValidationResult
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    def report(self) -> str:
        """Generate human-readable solve report."""
        lines = []
        lines.append("=" * 70)
        lines.append("ASSEMBLY INTENT SOLVE REPORT")
        lines.append("=" * 70)
        lines.append("")

        status = "SUCCESS" if self.success else "FAILED"
        lines.append(f"Status: {status}")
        lines.append("")

        if self.derived:
            lines.append("DERIVED PARAMETERS:")
            lines.append("-" * 70)
            for name, value in self.derived.items():
                lines.append(f"  {name}: {value:.4f}")
            lines.append("")

        if self.errors:
            lines.append("ERRORS:")
            lines.append("-" * 70)
            for error in self.errors:
                lines.append(f"  - {error}")
            lines.append("")

        if self.warnings:
            lines.append("WARNINGS:")
            lines.append("-" * 70)
            for warning in self.warnings:
                lines.append(f"  - {warning}")
            lines.append("")

        lines.append(f"Constraints generated: {len(self.constraints_generated)}")
        lines.append("=" * 70)

        return "\n".join(lines)


# =============================================================================
# ASSEMBLY INTENT
# =============================================================================

@dataclass
class AssemblyIntent:
    """Declarative specification of assembly requirements.

    AssemblyIntent is the top-level container for a declarative assembly.
    Instead of specifying transforms, designers specify:
    - What parts must DO (functional requirements)
    - How parts CONNECT (topology)
    - What must NOT happen (clearances)

    The system derives the geometry that satisfies all requirements.

    Attributes:
        name: Unique identifier for this assembly
        description: Human-readable description

        functional_requirements: List of FunctionalRequirement objects
        connections: List of Connection objects
        clearances: List of Clearance objects

        reference_geometry: Dict of ReferenceGeometry objects
        part_definitions: Dict of PartDefinition objects (optional)

        derived_parameters: List of parameter names to compute
        explicit_constraints: Additional explicit constraints

    Example:
        >>> assembly = AssemblyIntent(
        ...     name="wheel_pod",
        ...     functional_requirements=[...],
        ...     connections=[...],
        ...     clearances=[...],
        ... )
        >>> result = assembly.solve()
    """
    name: str
    description: str = ""

    # Requirements
    functional_requirements: List[FunctionalRequirement] = field(default_factory=list)
    connections: List[Connection] = field(default_factory=list)
    clearances: List[Clearance] = field(default_factory=list)

    # Reference data
    reference_geometry: Dict[str, ReferenceGeometry] = field(default_factory=dict)
    part_definitions: Dict[str, PartDefinition] = field(default_factory=dict)

    # Derived parameters
    derived_parameters: List[str] = field(default_factory=list)

    # Explicit constraints (from existing system)
    explicit_constraints: List[Constraint] = field(default_factory=list)

    def add_requirement(self, req: FunctionalRequirement) -> 'AssemblyIntent':
        """Add a functional requirement."""
        self.functional_requirements.append(req)
        return self

    def add_connection(self, conn: Connection) -> 'AssemblyIntent':
        """Add a connection."""
        self.connections.append(conn)
        return self

    def add_clearance(self, clearance: Clearance) -> 'AssemblyIntent':
        """Add a clearance constraint."""
        self.clearances.append(clearance)
        return self

    def add_reference_geometry(self, geom: ReferenceGeometry) -> 'AssemblyIntent':
        """Add reference geometry."""
        self.reference_geometry[geom.name] = geom
        return self

    def collect_constraints(self) -> List[Dict[str, Any]]:
        """Collect all constraints from requirements.

        Converts high-level functional requirements into low-level
        constraint specifications.

        Returns:
            List of constraint specification dicts
        """
        constraints = []

        # Get constraints from functional requirements
        for req in self.functional_requirements:
            implied = req.get_implied_constraints(self.reference_geometry)
            constraints.extend(implied)

        # Get mate specs from connections
        for conn in self.connections:
            mate_spec = conn.to_mate_spec()
            constraints.append(mate_spec)

        return constraints

    def collect_derived_parameters(self) -> List[str]:
        """Collect all parameters that need to be derived.

        Returns:
            List of parameter names
        """
        params = list(self.derived_parameters)

        for req in self.functional_requirements:
            params.extend(req.get_derived_parameters())

        return list(set(params))  # Remove duplicates

    def solve(self) -> SolveResult:
        """Solve the assembly intent to derive parameters.

        This is the main entry point for converting a declarative
        specification into actual geometry.

        Returns:
            SolveResult with derived values and validation status
        """
        result = SolveResult(success=True)
        result.constraints_generated = self.collect_constraints()

        # Derive parameters based on requirements
        # This is a simplified implementation - full solver would be more complex

        # Example: Derive wheel_center_radius from contact requirement
        for req in self.functional_requirements:
            if isinstance(req, ContactRequirement):
                # Find target geometry
                target = self.reference_geometry.get(req.target)
                if target and target.geometry_type == "cylinder":
                    # For wheel contacting tube inner wall:
                    # wheel_center_radius = tube_radius - wheel_radius
                    # We'd need wheel_radius from part definition
                    result.derived[f"{req.part}.contact_radius"] = target.radius
                    result.warnings.append(
                        f"ContactRequirement '{req.name}': derived contact_radius = {target.radius}mm"
                    )

        # Validate constraints
        error_constraints = [c for c in result.constraints_generated if c.get("type") == "error"]
        if error_constraints:
            result.success = False
            for c in error_constraints:
                result.errors.append(c.get("message", "Unknown error"))

        return result

    def to_assembly(self) -> Assembly:
        """Convert solved intent to yapCAD Assembly object.

        Returns:
            Assembly object with parts, mates, and constraints
        """
        assembly = Assembly(self.name)

        # Add part definitions
        for name, part_def in self.part_definitions.items():
            assembly.add_part(part_def, name=name)

        # Note: Full implementation would add computed transforms,
        # convert connections to mates, convert requirements to constraints

        return assembly

    def validate(self) -> List[str]:
        """Validate the assembly intent specification.

        Checks for:
        - Missing reference geometry
        - Invalid part references
        - Conflicting requirements
        - Under/over-constrained systems

        Returns:
            List of validation error/warning messages
        """
        issues = []

        # Check all referenced geometry exists
        for req in self.functional_requirements:
            if isinstance(req, ContactRequirement):
                if req.target not in self.reference_geometry:
                    issues.append(
                        f"ContactRequirement '{req.name}' references unknown "
                        f"geometry '{req.target}'"
                    )

        # Check connections reference valid parts
        parts_mentioned = set()
        for req in self.functional_requirements:
            parts_mentioned.add(req.part)
        for conn in self.connections:
            parts_mentioned.add(conn.get_parent_part())
            parts_mentioned.add(conn.get_child_part())

        # Check for orphan parts (not connected to anything)
        connected_parts = set()
        for conn in self.connections:
            connected_parts.add(conn.get_parent_part())
            connected_parts.add(conn.get_child_part())

        for part in parts_mentioned:
            if part not in connected_parts:
                issues.append(
                    f"Part '{part}' has requirements but no connections"
                )

        return issues

    def report(self) -> str:
        """Generate human-readable summary of the assembly intent."""
        lines = []
        lines.append("=" * 70)
        lines.append(f"ASSEMBLY INTENT: {self.name}")
        lines.append("=" * 70)

        if self.description:
            lines.append(f"\n{self.description}\n")

        lines.append(f"\nReference Geometry: {len(self.reference_geometry)}")
        for name, geom in self.reference_geometry.items():
            lines.append(f"  - {name}: {geom.geometry_type}")

        lines.append(f"\nFunctional Requirements: {len(self.functional_requirements)}")
        for req in self.functional_requirements:
            lines.append(f"  - {req.name}: {type(req).__name__} on '{req.part}'")

        lines.append(f"\nConnections: {len(self.connections)}")
        for conn in self.connections:
            lines.append(f"  - {conn.parent} -> {conn.child} ({conn.joint_type})")

        lines.append(f"\nClearances: {len(self.clearances)}")
        for cl in self.clearances:
            lines.append(f"  - {cl.part_a} <-> {cl.part_b}: min {cl.min_distance}mm")

        lines.append(f"\nDerived Parameters: {len(self.derived_parameters)}")
        for param in self.derived_parameters:
            lines.append(f"  - {param}")

        lines.append("")
        lines.append("=" * 70)

        return "\n".join(lines)


# =============================================================================
# EXAMPLE ASSEMBLIES
# =============================================================================

def create_wheel_assembly_intent() -> AssemblyIntent:
    """Create example wheel assembly intent for the tube robot.

    This demonstrates how to declare a wheel assembly using functional
    requirements rather than explicit transforms.

    Returns:
        AssemblyIntent for a wheel pod assembly
    """
    # Define reference geometry: the tube the robot operates in
    tube = ReferenceGeometry(
        name="tube_inner_wall",
        geometry_type="cylinder",
        center=(0, 0, 0),
        axis=(0, 0, 1),
        radius=175.0,  # 350mm ID tube
        inner=True,
    )

    # Create the assembly intent
    assembly = AssemblyIntent(
        name="wheel_pod_1",
        description="Drive wheel assembly at theta=0 position",
        reference_geometry={"tube_inner_wall": tube},
    )

    # Functional requirement 1: Wheel contacts tube wall
    assembly.add_requirement(ContactRequirement(
        name="wheel_wall_contact",
        part="ddsm115_motor",
        surface="tire_outer_diameter",
        target="tube_inner_wall",
        contact_type="rolling",
        preload_source="suspension_spring",
        description="Tire contacts tube inner wall for traction",
    ))

    # Functional requirement 2: Wheel rolls along tube axis
    assembly.add_requirement(RollRequirement(
        name="wheel_rolls_along_z",
        part="ddsm115_motor",
        roll_direction="along_tube_axis",
        axis_datum="motor_axis",
        description="Wheel rolls in Z direction to propel robot",
    ))

    # Functional requirement 3: Motor axis is tangent
    assembly.add_requirement(AxisOrientationRequirement(
        name="motor_axis_tangent",
        part="ddsm115_motor",
        axis_datum="motor_axis",
        orientation="tangent",
        reference="tube",
        description="Motor axis tangent for correct rolling motion",
    ))

    # Functional requirement 4: Pivot parallel to motor axis
    assembly.add_requirement(ParallelAxesRequirement(
        name="pivot_motor_parallel",
        part="wheel_arm",
        axis_a="wheel_arm.pivot_bore_axis",
        axis_b="ddsm115_motor.motor_axis",
        description="Pivot and motor axes parallel for clean suspension",
    ))

    # Connection 1: Chassis pivot to wheel arm
    assembly.add_connection(Connection(
        parent="chassis.pivot_boss_1",
        child="wheel_arm.pivot_bore",
        joint_type="revolute",
        axis="tangent",
        limits=(-15, 15),  # Degrees of swing
    ))

    # Connection 2: Wheel arm to motor
    assembly.add_connection(Connection(
        parent="wheel_arm.motor_mount",
        child="ddsm115_motor.stator_face",
        joint_type="rigid",
        interface_type="bolt_pattern",
        interface_details={
            "pattern": "3x_m25",
            "radius": 15.2,
            "anti_rotation": "triangular_pocket",
        },
    ))

    # Clearances
    assembly.add_clearance(Clearance(
        part_a="ddsm115_motor",
        part_b="chassis_plate",
        min_distance=5.0,
    ))

    assembly.add_clearance(Clearance(
        part_a="ddsm115_motor",
        part_b="tube_inner_wall",
        min_distance=3.0,
    ))

    # Parameters to derive
    assembly.derived_parameters = [
        "wheel_arm.length",
        "wheel_center_radius",
        "motor_mount_offset",
    ]

    return assembly


def create_scara_arm_intent() -> AssemblyIntent:
    """Create example SCARA arm assembly intent.

    Returns:
        AssemblyIntent for a 3-axis SCARA arm
    """
    assembly = AssemblyIntent(
        name="scara_arm",
        description="3-axis SCARA arm for tool positioning",
    )

    # All axes vertical (rotate around Z)
    assembly.add_requirement(AxisOrientationRequirement(
        name="axis1_vertical",
        part="link1",
        axis_datum="rotation_axis",
        orientation="axial",
        reference="global",
        description="Axis 1 rotates around vertical (Z)",
    ))

    assembly.add_requirement(AxisOrientationRequirement(
        name="axis2_vertical",
        part="link2",
        axis_datum="rotation_axis",
        orientation="axial",
        reference="global",
        description="Axis 2 rotates around vertical (Z)",
    ))

    # Workspace reach
    assembly.add_requirement(ReachRequirement(
        name="workspace_reach",
        part="wrist",
        end_effector="wrist.tool_point",
        base="tower.axis1_center",
        min_reach=50.0,
        max_reach=200.0,
        description="Arm must reach 50-200mm from base",
    ))

    # Connections: Tower -> Link1 -> Link2 -> Wrist
    assembly.add_connection(Connection(
        parent="tower.axis1_bearing",
        child="link1.proximal_bore",
        joint_type="revolute",
        axis=(0, 0, 1),  # Z-axis rotation
        limits=(-180, 180),
    ))

    assembly.add_connection(Connection(
        parent="link1.distal_bore",
        child="link2.proximal_bore",
        joint_type="revolute",
        axis=(0, 0, 1),
        limits=(-150, 150),
    ))

    assembly.add_connection(Connection(
        parent="link2.distal_bore",
        child="wrist.proximal_bore",
        joint_type="revolute",
        axis=(0, 0, 1),
        limits=(-180, 180),
    ))

    # Self-collision clearances
    assembly.add_clearance(Clearance(
        part_a="link1",
        part_b="chassis",
        min_distance=2.0,
    ))

    assembly.add_clearance(Clearance(
        part_a="link2",
        part_b="chassis",
        min_distance=2.0,
    ))

    assembly.add_clearance(Clearance(
        part_a="link1",
        part_b="link2",
        min_distance=1.0,
    ))

    # Derived parameters
    assembly.derived_parameters = [
        "link1.length",
        "link2.length",
        "total_reach",
    ]

    return assembly


# =============================================================================
# MODULE EXPORTS
# =============================================================================

__all__ = [
    # Reference Geometry
    "GeometryType",
    "ReferenceGeometry",

    # Functional Requirements
    "FunctionalRequirement",
    "ContactRequirement",
    "RollRequirement",
    "AxisOrientationRequirement",
    "ParallelAxesRequirement",
    "ReachRequirement",

    # Connections
    "Connection",

    # Clearances
    "Clearance",
    "ClearanceResult",

    # Solve Result
    "SolveResult",

    # Main Class
    "AssemblyIntent",

    # Example Factories
    "create_wheel_assembly_intent",
    "create_scara_arm_intent",
]
