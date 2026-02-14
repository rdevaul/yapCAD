"""Face-to-face mate system with automatic transform computation.

This module provides the FaceToFaceMate class which computes transforms
automatically from datum face constraints, replacing hardcoded transform
values with explicit geometric relationships.

Key Features:
    - Automatic transform computation from datum geometry
    - Support for FLUSH, ALIGNED, OFFSET, and CONCENTRIC constraints
    - Validation of mate satisfaction with tolerances
    - Integration with DatumRegistry for cross-DSL datum lookup

Example usage::

    from yapcad.assembly.face_mate import FaceToFaceMate
    from yapcad.assembly.mate import MateType

    # Define a mate between servo stator and link mounting face
    servo_mate = FaceToFaceMate(
        name="axis3_servo_to_link",
        parent_part="LINK_2_3",
        parent_datum="servo_mount_face",
        parent_source="scara_arm/scara_arm.dsl",
        child_part="AXIS3_SERVO_XH430",
        child_datum="stator_mounting_face",
        child_source="cots/xh430_surrogate.json",
        constraint=MateType.COINCIDENT,
    )

    # Compute the transform automatically
    transform_matrix = servo_mate.compute_transform()

Copyright (c) 2026 yapCAD contributors
License: MIT
"""

from __future__ import annotations

import math
import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any, Tuple, TYPE_CHECKING

from .datum import Datum, DatumType
from .mate import MateType

if TYPE_CHECKING:
    import numpy as np

logger = logging.getLogger(__name__)

# Check for numpy
try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None


@dataclass
class MateValidationResult:
    """Result of validating a mate constraint.

    Attributes:
        satisfied: True if constraint is satisfied within tolerance
        error_distance: Position error in mm
        error_angle: Angular error in degrees
        error_message: Human-readable description
        details: Additional validation information
    """
    satisfied: bool
    error_distance: float = 0.0
    error_angle: float = 0.0
    error_message: str = ""
    details: Dict[str, Any] = field(default_factory=dict)

    def __str__(self) -> str:
        status = "PASS" if self.satisfied else "FAIL"
        errors = []
        if self.error_distance > 0:
            errors.append(f"distance: {self.error_distance:.3f}mm")
        if self.error_angle > 0:
            errors.append(f"angle: {self.error_angle:.2f}deg")
        error_str = ", ".join(errors) if errors else "none"
        return f"[{status}] {self.error_message} (error: {error_str})"


@dataclass
class FaceToFaceMate:
    """Defines how two parts mate via their datum faces.

    This is the key abstraction that replaces hardcoded transforms with
    explicit geometric constraints. The transform is computed automatically
    from the datum geometry.

    Constraint Types:
        - COINCIDENT/FLUSH: Faces touch, normals point in opposite directions
        - CONCENTRIC: Axes are colinear (for cylindrical features)
        - PARALLEL: Faces parallel with optional offset
        - ALIGNED: Faces touch, normals point in same direction

    Attributes:
        name: Descriptive name for this mate
        parent_part: Name of the parent (fixed) part
        parent_datum: Datum name on parent part
        parent_source: Source file for parent datums (DSL or surrogate JSON)
        child_part: Name of the child (moving) part
        child_datum: Datum name on child part
        child_source: Source file for child datums
        constraint: Constraint type (MateType enum)
        offset_distance: Optional offset along normal (for OFFSET constraint)
        rotation_about_normal: Optional rotation about mating axis (degrees)
    """
    name: str

    # Parent (fixed) part and its datum (required fields first)
    parent_part: str
    parent_datum: str

    # Child (moving) part and its datum (required fields)
    child_part: str
    child_datum: str

    # Optional source files (default = use part name as source)
    parent_source: Optional[str] = None
    child_source: Optional[str] = None

    # Constraint type
    constraint: MateType = MateType.COINCIDENT

    # Optional offset along normal
    offset_distance: float = 0.0

    # Optional rotation about the mating axis (degrees)
    rotation_about_normal: float = 0.0

    # Cached computed transform
    _computed_transform: Optional[Any] = field(default=None, repr=False)

    # Cached datum references
    _parent_datum_obj: Optional[Datum] = field(default=None, repr=False)
    _child_datum_obj: Optional[Datum] = field(default=None, repr=False)

    def compute_transform(self,
                          parent_world_transform: Optional[Any] = None
                          ) -> Optional[Any]:
        """Compute the transform that places child relative to parent.

        This is the core algorithm that replaces hardcoded transforms.

        Args:
            parent_world_transform: Optional 4x4 world transform of parent
                                    (for computing absolute child position)

        Returns:
            4x4 numpy array transform, or None if numpy not available

        Raises:
            KeyError: If datums cannot be found
            ValueError: If constraint type not supported
        """
        if not HAS_NUMPY:
            logger.error("NumPy required for transform computation")
            return None

        if self._computed_transform is not None:
            return self._computed_transform

        # Load datums from registry
        from .datum_registry import DatumRegistry

        parent_source = self.parent_source or self.parent_part
        child_source = self.child_source or self.child_part

        parent_datum = DatumRegistry.get_datum(parent_source, self.parent_datum)
        child_datum = DatumRegistry.get_datum(child_source, self.child_datum)

        # Cache for later use
        self._parent_datum_obj = parent_datum
        self._child_datum_obj = child_datum

        # Compute transform based on constraint type
        if self.constraint in (MateType.COINCIDENT, MateType.ALIGNED):
            flip_normal = (self.constraint == MateType.COINCIDENT)
            tf = compute_flush_transform(
                parent_datum, child_datum,
                offset=self.offset_distance,
                rotation_deg=self.rotation_about_normal,
                flip_normal=flip_normal
            )
        elif self.constraint == MateType.CONCENTRIC:
            tf = compute_concentric_transform(
                parent_datum, child_datum,
                offset=self.offset_distance
            )
        elif self.constraint == MateType.PARALLEL:
            tf = compute_parallel_transform(
                parent_datum, child_datum,
                offset=self.offset_distance
            )
        else:
            raise ValueError(f"Unsupported constraint type: {self.constraint}")

        # Apply parent world transform if provided
        if parent_world_transform is not None:
            tf = parent_world_transform @ tf

        self._computed_transform = tf
        return tf

    def validate(self,
                 tolerance_mm: float = 0.1,
                 tolerance_deg: float = 1.0) -> MateValidationResult:
        """Validate that the mate constraint is satisfied.

        Args:
            tolerance_mm: Linear tolerance in millimeters
            tolerance_deg: Angular tolerance in degrees

        Returns:
            MateValidationResult with satisfaction status and errors
        """
        if not HAS_NUMPY:
            return MateValidationResult(
                satisfied=False,
                error_message="NumPy required for validation"
            )

        # Ensure transform is computed
        if self._computed_transform is None:
            try:
                self.compute_transform()
            except (KeyError, ValueError) as e:
                return MateValidationResult(
                    satisfied=False,
                    error_message=f"Failed to compute transform: {e}"
                )

        if self._parent_datum_obj is None or self._child_datum_obj is None:
            return MateValidationResult(
                satisfied=False,
                error_message="Datums not loaded"
            )

        # Transform child datum to parent frame
        child_origin = np.array(self._child_datum_obj.origin[:3])
        child_transformed = self._computed_transform[:3, :3] @ child_origin + \
                           self._computed_transform[:3, 3]

        parent_origin = np.array(self._parent_datum_obj.origin[:3])

        # Check distance
        distance = float(np.linalg.norm(child_transformed - parent_origin))
        distance_ok = distance <= tolerance_mm

        # Check angle (for PLANE/AXIS datums)
        angle_error = 0.0
        angle_ok = True

        if self._parent_datum_obj.datum_type in (DatumType.PLANE, DatumType.CIRCLE):
            parent_normal = np.array(self._parent_datum_obj.normal[:3])
            parent_normal = parent_normal / np.linalg.norm(parent_normal)

            child_normal = np.array(self._child_datum_obj.normal[:3])
            child_normal_transformed = self._computed_transform[:3, :3] @ child_normal
            child_normal_transformed = child_normal_transformed / \
                                       np.linalg.norm(child_normal_transformed)

            # For COINCIDENT, normals should be opposite
            if self.constraint == MateType.COINCIDENT:
                target = -parent_normal
            else:
                target = parent_normal

            dot = float(np.dot(child_normal_transformed, target))
            dot = max(-1.0, min(1.0, dot))  # Clamp for numerical stability
            angle_error = math.degrees(math.acos(dot))
            angle_ok = angle_error <= tolerance_deg

        satisfied = distance_ok and angle_ok
        status = "satisfied" if satisfied else "not satisfied"

        return MateValidationResult(
            satisfied=satisfied,
            error_distance=distance,
            error_angle=angle_error,
            error_message=f"Mate '{self.name}' {status}",
            details={
                "parent_origin": parent_origin.tolist(),
                "child_transformed": child_transformed.tolist(),
                "constraint": self.constraint.value
            }
        )

    def invalidate_cache(self) -> None:
        """Clear cached transform (call if datums change)."""
        self._computed_transform = None
        self._parent_datum_obj = None
        self._child_datum_obj = None


# =============================================================================
# Transform Computation Functions
# =============================================================================

def compute_flush_transform(parent_datum: Datum,
                            child_datum: Datum,
                            offset: float = 0.0,
                            rotation_deg: float = 0.0,
                            flip_normal: bool = True) -> Any:
    """Compute transform to make child flush with parent.

    The result is a transform T such that:
    - child_datum.origin moves to parent_datum.origin (+ offset along normal)
    - child_datum.normal aligns with parent_datum.normal
      (opposite if flip_normal=True for face-to-face contact)

    Uses Rodrigues' rotation formula for efficient axis-angle rotation.

    Args:
        parent_datum: Parent datum (target)
        child_datum: Child datum (to be transformed)
        offset: Optional offset along parent normal
        rotation_deg: Optional rotation about mating axis
        flip_normal: If True, child normal opposes parent (face-to-face)

    Returns:
        4x4 numpy transform matrix
    """
    if not HAS_NUMPY:
        raise RuntimeError("NumPy required for transform computation")

    # Extract vectors
    P_origin = np.array(parent_datum.origin[:3])
    C_origin = np.array(child_datum.origin[:3])

    # Get normals based on datum type
    if parent_datum.datum_type == DatumType.AXIS:
        P_normal = np.array(parent_datum.direction[:3])
    elif parent_datum.normal is not None:
        P_normal = np.array(parent_datum.normal[:3])
    else:
        P_normal = np.array([0, 0, 1])

    if child_datum.datum_type == DatumType.AXIS:
        C_normal = np.array(child_datum.direction[:3])
    elif child_datum.normal is not None:
        C_normal = np.array(child_datum.normal[:3])
    else:
        C_normal = np.array([0, 0, 1])

    # Normalize
    P_normal = P_normal / np.linalg.norm(P_normal)
    C_normal = C_normal / np.linalg.norm(C_normal)

    # Target: child normal should align with (or oppose) parent normal
    target_normal = -P_normal if flip_normal else P_normal

    # Step 1: Rotation to align normals using Rodrigues' formula
    R = _rotation_between_vectors(C_normal, target_normal)

    # Step 2: Apply rotation to child origin
    C_origin_rotated = R @ C_origin

    # Step 3: Compute translation
    t = P_origin - C_origin_rotated

    # Add offset along parent normal
    if offset != 0:
        t += offset * P_normal

    # Step 4: Build 4x4 transform
    T = np.eye(4)
    T[:3, :3] = R
    T[:3, 3] = t

    # Step 5: Optional rotation about mating axis
    if rotation_deg != 0:
        T = _apply_rotation_about_axis(T, P_normal, P_origin, rotation_deg)

    return T


def compute_concentric_transform(parent_datum: Datum,
                                 child_datum: Datum,
                                 offset: float = 0.0) -> Any:
    """Compute transform to make child axis colinear with parent axis.

    For CONCENTRIC mates, the axes must be colinear but origins can differ.

    Args:
        parent_datum: Parent axis datum
        child_datum: Child axis datum
        offset: Offset along axis direction

    Returns:
        4x4 numpy transform matrix
    """
    if not HAS_NUMPY:
        raise RuntimeError("NumPy required for transform computation")

    # Get axis directions
    if parent_datum.direction is not None:
        P_axis = np.array(parent_datum.direction[:3])
    elif parent_datum.normal is not None:
        P_axis = np.array(parent_datum.normal[:3])
    else:
        P_axis = np.array([0, 0, 1])

    if child_datum.direction is not None:
        C_axis = np.array(child_datum.direction[:3])
    elif child_datum.normal is not None:
        C_axis = np.array(child_datum.normal[:3])
    else:
        C_axis = np.array([0, 0, 1])

    P_axis = P_axis / np.linalg.norm(P_axis)
    C_axis = C_axis / np.linalg.norm(C_axis)

    # Rotation to align axes
    R = _rotation_between_vectors(C_axis, P_axis)

    # Project child origin onto parent axis
    P_origin = np.array(parent_datum.origin[:3])
    C_origin = np.array(child_datum.origin[:3])
    C_origin_rotated = R @ C_origin

    # Translation: move child origin to parent axis line
    # Find closest point on parent axis to rotated child origin
    v = C_origin_rotated - P_origin
    along_axis = np.dot(v, P_axis)
    closest_point = P_origin + along_axis * P_axis

    # For concentric, we want origins to coincide (or be offset along axis)
    t = P_origin - C_origin_rotated + offset * P_axis

    T = np.eye(4)
    T[:3, :3] = R
    T[:3, 3] = t

    return T


def compute_parallel_transform(parent_datum: Datum,
                               child_datum: Datum,
                               offset: float = 0.0) -> Any:
    """Compute transform to make child parallel to parent.

    For PARALLEL mates, normals align but origins don't necessarily coincide.

    Args:
        parent_datum: Parent plane datum
        child_datum: Child plane datum
        offset: Offset along normal direction

    Returns:
        4x4 numpy transform matrix
    """
    if not HAS_NUMPY:
        raise RuntimeError("NumPy required for transform computation")

    # Get normals
    P_normal = np.array(parent_datum.normal[:3]) if parent_datum.normal else np.array([0, 0, 1])
    C_normal = np.array(child_datum.normal[:3]) if child_datum.normal else np.array([0, 0, 1])

    P_normal = P_normal / np.linalg.norm(P_normal)
    C_normal = C_normal / np.linalg.norm(C_normal)

    # Rotation to make parallel
    R = _rotation_between_vectors(C_normal, P_normal)

    # For parallel, we maintain the relative XY position
    P_origin = np.array(parent_datum.origin[:3])
    C_origin = np.array(child_datum.origin[:3])
    C_origin_rotated = R @ C_origin

    # Only adjust Z (along normal) by offset
    t = C_origin - C_origin_rotated + offset * P_normal

    T = np.eye(4)
    T[:3, :3] = R
    T[:3, 3] = t

    return T


# =============================================================================
# Helper Functions
# =============================================================================

def _rotation_between_vectors(v1: Any, v2: Any) -> Any:
    """Compute rotation matrix to rotate v1 to v2 using Rodrigues' formula.

    Args:
        v1: Source unit vector (3,)
        v2: Target unit vector (3,)

    Returns:
        3x3 rotation matrix
    """
    v1 = v1 / np.linalg.norm(v1)
    v2 = v2 / np.linalg.norm(v2)

    # Cross product gives axis of rotation
    v = np.cross(v1, v2)
    s = np.linalg.norm(v)  # sin(angle)
    c = np.dot(v1, v2)     # cos(angle)

    if s < 1e-10:
        # Vectors are parallel or anti-parallel
        if c > 0:
            # Same direction - no rotation needed
            return np.eye(3)
        else:
            # Opposite direction - 180 degree rotation
            # Find any perpendicular axis
            if abs(v1[0]) < 0.9:
                perp = np.cross(v1, np.array([1, 0, 0]))
            else:
                perp = np.cross(v1, np.array([0, 1, 0]))
            perp = perp / np.linalg.norm(perp)
            # Rodrigues for 180 degrees: R = 2*perp*perp^T - I
            return 2 * np.outer(perp, perp) - np.eye(3)

    # General case - Rodrigues' formula
    # R = I + sin(θ)[v]× + (1-cos(θ))[v]×²
    v_normalized = v / s
    vx = np.array([
        [0, -v_normalized[2], v_normalized[1]],
        [v_normalized[2], 0, -v_normalized[0]],
        [-v_normalized[1], v_normalized[0], 0]
    ])
    R = np.eye(3) + s * vx + (1 - c) * (vx @ vx)

    return R


def _apply_rotation_about_axis(T: Any, axis: Any, point: Any,
                               angle_deg: float) -> Any:
    """Apply additional rotation about an axis through a point.

    Args:
        T: 4x4 transform to modify
        axis: Rotation axis (unit vector)
        point: Point on the axis
        angle_deg: Rotation angle in degrees

    Returns:
        Modified 4x4 transform
    """
    theta = np.radians(angle_deg)
    c, s = np.cos(theta), np.sin(theta)
    n = axis / np.linalg.norm(axis)

    # Rodrigues for rotation about arbitrary axis
    nx = np.array([
        [0, -n[2], n[1]],
        [n[2], 0, -n[0]],
        [-n[1], n[0], 0]
    ])
    R_about = c * np.eye(3) + s * nx + (1 - c) * np.outer(n, n)

    # Build transform that rotates about 'point'
    T_rot = np.eye(4)
    T_rot[:3, :3] = R_about

    # Translate to origin, rotate, translate back
    T_to_origin = np.eye(4)
    T_to_origin[:3, 3] = -point

    T_from_origin = np.eye(4)
    T_from_origin[:3, 3] = point

    return T_from_origin @ T_rot @ T_to_origin @ T


# =============================================================================
# Convenience Functions
# =============================================================================

def create_face_mate(name: str,
                     parent_part: str,
                     parent_datum: str,
                     child_part: str,
                     child_datum: str,
                     parent_source: Optional[str] = None,
                     child_source: Optional[str] = None,
                     offset: float = 0.0,
                     rotation: float = 0.0) -> FaceToFaceMate:
    """Create a face-to-face (COINCIDENT) mate.

    Convenience function for the most common mate type.

    Args:
        name: Mate name
        parent_part: Parent part name
        parent_datum: Parent datum name
        child_part: Child part name
        child_datum: Child datum name
        parent_source: Optional source file for parent
        child_source: Optional source file for child
        offset: Optional offset along normal
        rotation: Optional rotation about normal (degrees)

    Returns:
        Configured FaceToFaceMate
    """
    return FaceToFaceMate(
        name=name,
        parent_part=parent_part,
        parent_datum=parent_datum,
        parent_source=parent_source,
        child_part=child_part,
        child_datum=child_datum,
        child_source=child_source,
        constraint=MateType.COINCIDENT,
        offset_distance=offset,
        rotation_about_normal=rotation
    )


def create_axis_mate(name: str,
                     parent_part: str,
                     parent_datum: str,
                     child_part: str,
                     child_datum: str,
                     parent_source: Optional[str] = None,
                     child_source: Optional[str] = None,
                     offset: float = 0.0) -> FaceToFaceMate:
    """Create a concentric (axis-aligned) mate.

    Convenience function for coaxial alignment.

    Args:
        name: Mate name
        parent_part: Parent part name
        parent_datum: Parent axis datum name
        child_part: Child part name
        child_datum: Child axis datum name
        parent_source: Optional source file for parent
        child_source: Optional source file for child
        offset: Optional offset along axis

    Returns:
        Configured FaceToFaceMate
    """
    return FaceToFaceMate(
        name=name,
        parent_part=parent_part,
        parent_datum=parent_datum,
        parent_source=parent_source,
        child_part=child_part,
        child_datum=child_datum,
        child_source=child_source,
        constraint=MateType.CONCENTRIC,
        offset_distance=offset
    )
