"""Assembly constraint validation for yapCAD.

This module provides high-level design constraints that validate assembly intent
after mates have positioned parts. Constraints check that the resulting assembly
meets design requirements like tangency, radial orientation, and clearances.

Constraints differ from mates:
- Mates POSITION parts (determine transforms)
- Constraints VALIDATE the result (check design intent)

Example:
    from yapcad.assembly.constraint import (
        Constraint, ConstraintType, ConstraintResult
    )

    # Motor axis must be tangent to chassis (not radial)
    constraint = Constraint(
        name="wheel_axis_tangent",
        constraint_type=ConstraintType.TANGENT_TO_CIRCLE,
        part="DDSM115_MOTOR_1",
        datum="motor_axis",
        center=(0.0, 0.0, 0.0),
        radius=124.5,
        description="Motor axis tangent to wheel path for rolling motion"
    )

    # Evaluate against a datum in world coordinates
    result = constraint.evaluate(datum_world)
    if not result.passed:
        print(f"Constraint violated: {result.error_message}")
        print(f"Error: {result.error_value:.2f}°")
"""

from dataclasses import dataclass, field
from typing import Optional, Tuple, Callable, Any, Dict, TYPE_CHECKING
from enum import Enum
import math

if TYPE_CHECKING:
    from .datum import Datum, DatumType

# Import numpy with fallback
try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    # Provide minimal fallback for type hints
    np = None


class ConstraintType(Enum):
    """High-level design constraints that validate assembly intent.

    These constraint types capture common design requirements in mechanical
    assemblies that are difficult to express with low-level mates alone.

    Attributes:
        TANGENT_TO_CIRCLE: Axis is tangent (perpendicular to radial) to a circle.
            Used for wheels that must roll along a circular path.
        RADIAL_FROM_CENTER: Datum normal points toward (inward) or away from
            (outward) a center point. Used for features that must face the center
            or periphery of an assembly.
        FACING: Datum normal points in a specified global direction (+x, -y, +z, etc.).
        AT_RADIUS: Datum origin is at a specific radius from a center point.
        PARALLEL_TO: Datum direction is parallel to a reference direction.
        PERPENDICULAR_TO: Datum direction is perpendicular to a reference direction.
    """

    # Directional constraints
    FACING = "facing"
    TANGENT_TO_CIRCLE = "tangent_to_circle"
    RADIAL_FROM_CENTER = "radial_from_center"
    PARALLEL_TO = "parallel_to"
    PERPENDICULAR_TO = "perpendicular_to"

    # Positional constraints
    AT_RADIUS = "at_radius"
    ON_CIRCLE = "on_circle"
    IN_PLANE = "in_plane"

    # Clearance constraints
    MIN_DISTANCE = "min_distance"
    NO_INTERFERENCE = "no_interference"

    # Custom constraints
    CUSTOM = "custom"


@dataclass
class ConstraintResult:
    """Result of evaluating a design constraint.

    Attributes:
        passed: True if constraint is satisfied within tolerance
        error_value: Numeric measure of constraint violation (degrees or mm)
        error_message: Human-readable description of the result
        details: Additional information about the evaluation
    """

    passed: bool
    error_value: float = 0.0
    error_message: str = ""
    details: Dict[str, Any] = field(default_factory=dict)

    def __str__(self) -> str:
        """Format result as a readable string."""
        status = "PASS" if self.passed else "FAIL"
        if self.error_value > 0:
            return f"[{status}] {self.error_message} (error: {self.error_value:.3f})"
        return f"[{status}] {self.error_message}"


@dataclass
class Constraint:
    """A design constraint that validates assembly intent.

    Constraints are evaluated after mates have positioned all parts.
    They verify that the resulting assembly meets design requirements.

    Unlike mates (which determine part positions), constraints validate
    that the positions satisfy high-level design intent.

    Attributes:
        name: Unique identifier for this constraint
        constraint_type: Type of constraint to evaluate
        part: Name of the part containing the datum to check
        datum: Name of the datum feature to evaluate
        center: Reference center point for circular/radial constraints
        radius: Reference radius for circular constraints (mm)
        direction: Reference direction ("inward", "outward", "+x", "-z", etc.)
        axis: Reference axis vector for parallel/perpendicular constraints
        reference_datum: Reference datum for relative constraints
        plane_normal: Normal vector for plane constraints
        min_value: Minimum allowed value for range constraints
        max_value: Maximum allowed value for range constraints
        validator: Custom validation function for CUSTOM constraint type
        tolerance_deg: Angular tolerance in degrees (default: 1.0°)
        tolerance_mm: Linear tolerance in millimeters (default: 0.1mm)
        description: Human-readable description of design intent
        severity: "error", "warning", or "info"

    Examples:
        # Motor axis tangent to wheel path (perpendicular to radial)
        >>> tangent = Constraint(
        ...     name="wheel_axis_tangent",
        ...     constraint_type=ConstraintType.TANGENT_TO_CIRCLE,
        ...     part="DDSM115_MOTOR_1",
        ...     datum="motor_axis",
        ...     center=(0.0, 0.0, 0.0),
        ...     radius=124.5,
        ...     tolerance_deg=2.0,
        ...     description="Motor axis must be tangent for rolling motion"
        ... )

        # Tire tread must face outward toward tube wall
        >>> radial = Constraint(
        ...     name="tire_faces_outward",
        ...     constraint_type=ConstraintType.RADIAL_FROM_CENTER,
        ...     part="DDSM115_MOTOR_1",
        ...     datum="tire_face",
        ...     center=(0.0, 0.0, 0.0),
        ...     direction="outward",
        ...     tolerance_deg=5.0,
        ...     description="Tire tread faces tube wall for traction"
        ... )

        # Mounting face must face upward
        >>> facing = Constraint(
        ...     name="mount_faces_up",
        ...     constraint_type=ConstraintType.FACING,
        ...     part="BRACKET",
        ...     datum="mounting_face",
        ...     direction="+z",
        ...     tolerance_deg=1.0,
        ...     description="Mounting face must be horizontal"
        ... )

        # Part at specific radius from center
        >>> at_radius = Constraint(
        ...     name="motor_at_wheel_radius",
        ...     constraint_type=ConstraintType.AT_RADIUS,
        ...     part="MOTOR",
        ...     datum="motor_center",
        ...     center=(0.0, 0.0, 0.0),
        ...     radius=124.5,
        ...     tolerance_mm=0.5,
        ...     description="Motor positioned on wheel path circle"
        ... )
    """

    name: str
    constraint_type: ConstraintType

    # Primary datum reference
    part: str
    datum: str

    # Constraint-specific parameters
    center: Optional[Tuple[float, float, float]] = None
    radius: Optional[float] = None
    direction: Optional[str] = None
    axis: Optional[Tuple[float, float, float]] = None
    reference_datum: Optional[str] = None
    plane_normal: Optional[Tuple[float, float, float]] = None
    min_value: Optional[float] = None
    max_value: Optional[float] = None

    # For CUSTOM constraints
    validator: Optional[Callable[[Any], bool]] = None

    # Tolerances
    tolerance_deg: float = 1.0
    tolerance_mm: float = 0.1

    # Metadata
    description: str = ""
    severity: str = "error"

    def evaluate(self, datum_world: 'Datum') -> ConstraintResult:
        """Evaluate this constraint against a datum in world coordinates.

        Args:
            datum_world: The datum feature transformed to world coordinates

        Returns:
            ConstraintResult indicating whether constraint is satisfied

        Raises:
            ValueError: If required parameters are missing
            ImportError: If numpy is required but not available
        """
        if not HAS_NUMPY:
            return ConstraintResult(
                passed=False,
                error_message="NumPy required for constraint evaluation"
            )

        if self.constraint_type == ConstraintType.TANGENT_TO_CIRCLE:
            return self._check_tangent_to_circle(datum_world)
        elif self.constraint_type == ConstraintType.RADIAL_FROM_CENTER:
            return self._check_radial(datum_world)
        elif self.constraint_type == ConstraintType.FACING:
            return self._check_facing(datum_world)
        elif self.constraint_type == ConstraintType.AT_RADIUS:
            return self._check_at_radius(datum_world)
        elif self.constraint_type == ConstraintType.PARALLEL_TO:
            return self._check_parallel(datum_world)
        elif self.constraint_type == ConstraintType.PERPENDICULAR_TO:
            return self._check_perpendicular(datum_world)
        elif self.constraint_type == ConstraintType.CUSTOM:
            if self.validator:
                try:
                    passed = self.validator(datum_world)
                    return ConstraintResult(
                        passed=passed,
                        error_message=f"Custom constraint '{self.name}'"
                    )
                except Exception as e:
                    return ConstraintResult(
                        passed=False,
                        error_message=f"Custom validator error: {e}"
                    )
            else:
                return ConstraintResult(
                    passed=False,
                    error_message="CUSTOM constraint requires validator function"
                )

        return ConstraintResult(
            passed=False,
            error_message=f"Unknown constraint type: {self.constraint_type}"
        )

    def _check_tangent_to_circle(self, datum: 'Datum') -> ConstraintResult:
        """Check if an axis is tangent to a circle at a given radius.

        An axis is tangent when it is perpendicular to the radial direction
        from the center to the datum origin. This is used for wheels that must
        roll along a circular path rather than pointing radially inward/outward.

        Args:
            datum: Datum in world coordinates (must have direction or normal)

        Returns:
            ConstraintResult with angular error from tangency
        """
        if self.center is None or self.radius is None:
            return ConstraintResult(
                passed=False,
                error_message="TANGENT_TO_CIRCLE requires center and radius"
            )

        center = np.array(self.center)
        origin = np.array(datum.origin[:3])

        # Vector from center to datum origin (project to XY plane for Z-axis circle)
        radial = origin - center
        radial[2] = 0.0  # Project to XY plane
        radial_norm = np.linalg.norm(radial)

        if radial_norm < 1e-10:
            return ConstraintResult(
                passed=False,
                error_message="Datum at center - cannot determine tangent direction"
            )

        radial_unit = radial / radial_norm

        # Get the direction to check
        if hasattr(datum, 'direction') and datum.direction is not None:
            check_dir = np.array(datum.direction[:3])
        elif hasattr(datum, 'normal') and datum.normal is not None:
            check_dir = np.array(datum.normal[:3])
        else:
            return ConstraintResult(
                passed=False,
                error_message="Datum has no direction or normal vector"
            )

        # Project direction to XY plane
        check_dir_xy = check_dir.copy()
        check_dir_xy[2] = 0.0
        dir_norm = np.linalg.norm(check_dir_xy)

        if dir_norm < 1e-10:
            # Direction is purely vertical - perpendicular to XY, so tangent to any XY circle
            return ConstraintResult(
                passed=True,
                error_value=0.0,
                error_message="Direction is vertical (tangent to horizontal circle)",
                details={"radial_direction": radial_unit.tolist()}
            )

        check_dir_unit = check_dir_xy / dir_norm

        # Tangent means perpendicular to radial (dot product = 0)
        dot = abs(np.dot(radial_unit, check_dir_unit))
        angle_from_perpendicular = math.degrees(math.asin(min(1.0, dot)))

        passed = angle_from_perpendicular <= self.tolerance_deg

        return ConstraintResult(
            passed=passed,
            error_value=angle_from_perpendicular,
            error_message=f"Tangent angle error: {angle_from_perpendicular:.2f}° "
                         f"(tolerance: {self.tolerance_deg}°)",
            details={
                "radial_direction": radial_unit.tolist(),
                "axis_direction": check_dir_unit.tolist(),
                "dot_product": float(dot)
            }
        )

    def _check_radial(self, datum: 'Datum') -> ConstraintResult:
        """Check if a datum normal faces radially inward or outward from center.

        This validates that a feature (like a tire tread or mounting face)
        points toward or away from a center point in the XY plane.

        Args:
            datum: Datum in world coordinates (must have normal or direction)

        Returns:
            ConstraintResult with angular error from radial alignment
        """
        if self.center is None:
            return ConstraintResult(
                passed=False,
                error_message="RADIAL_FROM_CENTER requires center parameter"
            )

        if self.direction not in ("inward", "outward"):
            return ConstraintResult(
                passed=False,
                error_message=f"direction must be 'inward' or 'outward', got '{self.direction}'"
            )

        center = np.array(self.center)
        origin = np.array(datum.origin[:3])

        # Radial direction (outward from center)
        radial = origin - center
        radial[2] = 0.0  # Project to XY plane for cylindrical radial
        radial_norm = np.linalg.norm(radial)

        if radial_norm < 1e-10:
            return ConstraintResult(
                passed=False,
                error_message="Datum at center - cannot determine radial direction"
            )

        radial_unit = radial / radial_norm

        # Get the direction to check
        if hasattr(datum, 'normal') and datum.normal is not None:
            check_dir = np.array(datum.normal[:3])
        elif hasattr(datum, 'direction') and datum.direction is not None:
            check_dir = np.array(datum.direction[:3])
        else:
            return ConstraintResult(
                passed=False,
                error_message="Datum has no normal or direction vector"
            )

        check_dir = check_dir / np.linalg.norm(check_dir)

        # Check alignment
        dot = np.dot(radial_unit, check_dir)

        if self.direction == "outward":
            # Normal should point same direction as radial (dot > 0)
            cos_tolerance = math.cos(math.radians(self.tolerance_deg))
            satisfied = dot > cos_tolerance
            expected = "outward (+radial)"
        else:  # inward
            # Normal should point opposite to radial (dot < 0)
            cos_tolerance = math.cos(math.radians(self.tolerance_deg))
            satisfied = dot < -cos_tolerance
            expected = "inward (-radial)"

        angle_error = math.degrees(math.acos(min(1.0, abs(dot))))
        if not satisfied:
            # Calculate actual angle from expected direction
            if self.direction == "outward" and dot < 0:
                angle_error = 180.0 - angle_error
            elif self.direction == "inward" and dot > 0:
                angle_error = 180.0 - angle_error

        return ConstraintResult(
            passed=satisfied,
            error_value=angle_error if not satisfied else 0.0,
            error_message=f"Radial alignment: dot={dot:.3f}, expected {expected}, "
                         f"angle_error={angle_error:.2f}°",
            details={
                "radial_direction": radial_unit.tolist(),
                "datum_direction": check_dir.tolist(),
                "dot_product": float(dot),
                "expected_direction": self.direction
            }
        )

    def _check_facing(self, datum: 'Datum') -> ConstraintResult:
        """Check if a datum normal faces a specified global direction.

        Validates that a plane or surface normal points in a cardinal direction
        like +z (up), -y (back), etc.

        Args:
            datum: Datum in world coordinates (must have normal)

        Returns:
            ConstraintResult with angular error from target direction
        """
        if not hasattr(datum, 'normal') or datum.normal is None:
            return ConstraintResult(
                passed=False,
                error_message="FACING requires datum with normal vector"
            )

        normal = np.array(datum.normal[:3])
        normal = normal / np.linalg.norm(normal)

        # Parse direction specification
        target = self._parse_direction(self.direction)
        if target is None:
            return ConstraintResult(
                passed=False,
                error_message=f"Invalid direction: {self.direction}"
            )

        dot = np.dot(normal, target)
        angle = math.degrees(math.acos(min(1.0, max(-1.0, dot))))

        satisfied = angle <= self.tolerance_deg

        return ConstraintResult(
            passed=satisfied,
            error_value=angle,
            error_message=f"Facing angle error: {angle:.2f}° from {self.direction} "
                         f"(tolerance: {self.tolerance_deg}°)",
            details={
                "datum_normal": normal.tolist(),
                "target_direction": target.tolist(),
                "dot_product": float(dot)
            }
        )

    def _check_at_radius(self, datum: 'Datum') -> ConstraintResult:
        """Check if a datum origin is at a specified radius from center.

        Validates radial positioning in the XY plane (cylindrical coordinates).

        Args:
            datum: Datum in world coordinates

        Returns:
            ConstraintResult with distance error
        """
        if self.center is None or self.radius is None:
            return ConstraintResult(
                passed=False,
                error_message="AT_RADIUS requires center and radius parameters"
            )

        center = np.array(self.center)
        origin = np.array(datum.origin[:3])

        # Distance in XY plane (cylindrical radius)
        delta = origin - center
        delta[2] = 0.0
        actual_radius = np.linalg.norm(delta)

        error = abs(actual_radius - self.radius)
        satisfied = error <= self.tolerance_mm

        return ConstraintResult(
            passed=satisfied,
            error_value=error,
            error_message=f"Radius: {actual_radius:.3f}mm "
                         f"(expected {self.radius:.3f}mm, error {error:.3f}mm, "
                         f"tolerance: {self.tolerance_mm}mm)",
            details={
                "actual_radius": float(actual_radius),
                "expected_radius": float(self.radius),
                "center": self.center
            }
        )

    def _check_parallel(self, datum: 'Datum') -> ConstraintResult:
        """Check if a datum direction is parallel to a reference axis.

        Args:
            datum: Datum in world coordinates (must have direction or normal)

        Returns:
            ConstraintResult with angular error from parallel
        """
        if self.axis is None:
            return ConstraintResult(
                passed=False,
                error_message="PARALLEL_TO requires axis parameter"
            )

        ref_axis = np.array(self.axis[:3])
        ref_axis = ref_axis / np.linalg.norm(ref_axis)

        if hasattr(datum, 'direction') and datum.direction is not None:
            check_dir = np.array(datum.direction[:3])
        elif hasattr(datum, 'normal') and datum.normal is not None:
            check_dir = np.array(datum.normal[:3])
        else:
            return ConstraintResult(
                passed=False,
                error_message="Datum has no direction or normal vector"
            )

        check_dir = check_dir / np.linalg.norm(check_dir)

        # Parallel means dot product is +1 or -1
        dot = abs(np.dot(ref_axis, check_dir))
        angle = math.degrees(math.acos(min(1.0, dot)))

        satisfied = angle <= self.tolerance_deg

        return ConstraintResult(
            passed=satisfied,
            error_value=angle,
            error_message=f"Parallel angle error: {angle:.2f}° "
                         f"(tolerance: {self.tolerance_deg}°)",
            details={
                "reference_axis": ref_axis.tolist(),
                "datum_direction": check_dir.tolist(),
                "dot_product": float(dot)
            }
        )

    def _check_perpendicular(self, datum: 'Datum') -> ConstraintResult:
        """Check if a datum direction is perpendicular to a reference axis.

        Args:
            datum: Datum in world coordinates (must have direction or normal)

        Returns:
            ConstraintResult with angular error from perpendicular
        """
        if self.axis is None:
            return ConstraintResult(
                passed=False,
                error_message="PERPENDICULAR_TO requires axis parameter"
            )

        ref_axis = np.array(self.axis[:3])
        ref_axis = ref_axis / np.linalg.norm(ref_axis)

        if hasattr(datum, 'direction') and datum.direction is not None:
            check_dir = np.array(datum.direction[:3])
        elif hasattr(datum, 'normal') and datum.normal is not None:
            check_dir = np.array(datum.normal[:3])
        else:
            return ConstraintResult(
                passed=False,
                error_message="Datum has no direction or normal vector"
            )

        check_dir = check_dir / np.linalg.norm(check_dir)

        # Perpendicular means dot product is 0
        dot = abs(np.dot(ref_axis, check_dir))
        angle_from_perpendicular = math.degrees(math.asin(min(1.0, dot)))

        satisfied = angle_from_perpendicular <= self.tolerance_deg

        return ConstraintResult(
            passed=satisfied,
            error_value=angle_from_perpendicular,
            error_message=f"Perpendicular angle error: {angle_from_perpendicular:.2f}° "
                         f"(tolerance: {self.tolerance_deg}°)",
            details={
                "reference_axis": ref_axis.tolist(),
                "datum_direction": check_dir.tolist(),
                "dot_product": float(dot)
            }
        )

    def _parse_direction(self, direction: Optional[str]) -> Optional[np.ndarray]:
        """Parse a direction string into a unit vector.

        Args:
            direction: String like "+x", "-y", "+z", etc.

        Returns:
            Unit vector as numpy array, or None if invalid
        """
        if direction is None:
            return None

        directions = {
            "+x": np.array([1.0, 0.0, 0.0]),
            "-x": np.array([-1.0, 0.0, 0.0]),
            "+y": np.array([0.0, 1.0, 0.0]),
            "-y": np.array([0.0, -1.0, 0.0]),
            "+z": np.array([0.0, 0.0, 1.0]),
            "-z": np.array([0.0, 0.0, -1.0]),
        }
        return directions.get(direction.lower())


# Helper functions for geometric calculations

def angle_between_vectors(v1: Tuple[float, float, float],
                         v2: Tuple[float, float, float]) -> float:
    """Calculate angle between two vectors in degrees.

    Args:
        v1: First vector (x, y, z)
        v2: Second vector (x, y, z)

    Returns:
        Angle in degrees [0, 180]

    Example:
        >>> angle_between_vectors((1, 0, 0), (0, 1, 0))
        90.0
        >>> angle_between_vectors((1, 0, 0), (1, 0, 0))
        0.0
    """
    if not HAS_NUMPY:
        raise ImportError("NumPy required for angle_between_vectors")

    a = np.array(v1)
    b = np.array(v2)

    a_norm = np.linalg.norm(a)
    b_norm = np.linalg.norm(b)

    if a_norm < 1e-10 or b_norm < 1e-10:
        return 0.0

    a_unit = a / a_norm
    b_unit = b / b_norm

    dot = np.dot(a_unit, b_unit)
    # Clamp to [-1, 1] to handle numerical errors
    dot = max(-1.0, min(1.0, dot))

    return math.degrees(math.acos(dot))


def distance_to_point(p1: Tuple[float, float, float],
                     p2: Tuple[float, float, float]) -> float:
    """Calculate Euclidean distance between two points.

    Args:
        p1: First point (x, y, z)
        p2: Second point (x, y, z)

    Returns:
        Distance in same units as input

    Example:
        >>> distance_to_point((0, 0, 0), (3, 4, 0))
        5.0
    """
    if not HAS_NUMPY:
        raise ImportError("NumPy required for distance_to_point")

    a = np.array(p1)
    b = np.array(p2)
    return float(np.linalg.norm(b - a))


def dot_product(v1: Tuple[float, float, float],
               v2: Tuple[float, float, float]) -> float:
    """Calculate dot product of two vectors.

    Args:
        v1: First vector (x, y, z)
        v2: Second vector (x, y, z)

    Returns:
        Dot product (scalar)

    Example:
        >>> dot_product((1, 0, 0), (0, 1, 0))
        0.0
        >>> dot_product((1, 2, 3), (1, 2, 3))
        14.0
    """
    if not HAS_NUMPY:
        raise ImportError("NumPy required for dot_product")

    a = np.array(v1)
    b = np.array(v2)
    return float(np.dot(a, b))


def is_tangent_to_circle(axis_origin: Tuple[float, float, float],
                        axis_direction: Tuple[float, float, float],
                        circle_center: Tuple[float, float, float],
                        circle_radius: float,
                        tolerance_deg: float = 1.0) -> bool:
    """Check if an axis is tangent to a circle in the XY plane.

    An axis is tangent when its direction is perpendicular to the radial
    vector from the circle center to the axis origin.

    Args:
        axis_origin: Point on the axis (x, y, z)
        axis_direction: Direction vector of the axis (x, y, z)
        circle_center: Center of the circle (x, y, z)
        circle_radius: Radius of the circle
        tolerance_deg: Angular tolerance in degrees

    Returns:
        True if axis is tangent within tolerance

    Example:
        >>> # Axis at (10, 0, 0) pointing in +y direction, tangent to origin circle
        >>> is_tangent_to_circle((10, 0, 0), (0, 1, 0), (0, 0, 0), 10.0)
        True
    """
    if not HAS_NUMPY:
        raise ImportError("NumPy required for is_tangent_to_circle")

    center = np.array(circle_center)
    origin = np.array(axis_origin)
    direction = np.array(axis_direction)

    # Radial vector from center to axis origin (XY plane only)
    radial = origin - center
    radial[2] = 0.0
    radial_norm = np.linalg.norm(radial)

    if radial_norm < 1e-10:
        return False  # Axis passes through center

    radial_unit = radial / radial_norm

    # Project axis direction to XY plane
    direction_xy = direction.copy()
    direction_xy[2] = 0.0
    dir_norm = np.linalg.norm(direction_xy)

    if dir_norm < 1e-10:
        return True  # Vertical axis is tangent to horizontal circle

    direction_unit = direction_xy / dir_norm

    # Tangent means perpendicular to radial (dot product ≈ 0)
    dot = abs(np.dot(radial_unit, direction_unit))
    angle_from_perpendicular = math.degrees(math.asin(min(1.0, dot)))

    return angle_from_perpendicular <= tolerance_deg


def is_radial_from_center(point: Tuple[float, float, float],
                         normal: Tuple[float, float, float],
                         center: Tuple[float, float, float],
                         direction: str = "outward",
                         tolerance_deg: float = 1.0) -> bool:
    """Check if a normal vector points radially from a center point.

    Args:
        point: Location of the datum (x, y, z)
        normal: Normal vector to check (x, y, z)
        center: Center point (x, y, z)
        direction: "outward" or "inward"
        tolerance_deg: Angular tolerance in degrees

    Returns:
        True if normal is radial within tolerance

    Example:
        >>> # Point at (10, 0, 0) with normal pointing in +x (outward from origin)
        >>> is_radial_from_center((10, 0, 0), (1, 0, 0), (0, 0, 0), "outward")
        True
    """
    if not HAS_NUMPY:
        raise ImportError("NumPy required for is_radial_from_center")

    if direction not in ("inward", "outward"):
        raise ValueError("direction must be 'inward' or 'outward'")

    c = np.array(center)
    p = np.array(point)
    n = np.array(normal)

    # Radial vector (outward)
    radial = p - c
    radial[2] = 0.0  # Project to XY plane
    radial_norm = np.linalg.norm(radial)

    if radial_norm < 1e-10:
        return False  # Point at center

    radial_unit = radial / radial_norm
    normal_unit = n / np.linalg.norm(n)

    dot = np.dot(radial_unit, normal_unit)

    cos_tolerance = math.cos(math.radians(tolerance_deg))

    if direction == "outward":
        return dot > cos_tolerance
    else:  # inward
        return dot < -cos_tolerance
