"""Generic Mate Constraint Solver for Multi-Body Assemblies.

====================
OVERVIEW
====================

This module provides a generic constraint solver that computes 6DOF transforms
from mate constraints between parts. It serves as the core computation engine
for datum-driven assembly positioning, eliminating hardcoded transform values.

Key Concepts:
    - **Mate Constraint**: Geometric relationship between two datum features
    - **Datum**: Named geometric reference (point, axis, plane, frame) on a part
    - **Transform**: 4x4 homogeneous transformation matrix positioning a part
    - **Mate Axis**: Coordinate remapping for rotated reference frames

Quick Start::

    from yapcad.assembly.solver import MateConstraintSolver, MateMateSolveResult
    from yapcad.assembly.mate import Mate, MateType
    from yapcad.assembly.datum_registry import DatumRegistry

    # Register datums for parts
    DatumRegistry.register_source("servo", servo_datums)
    DatumRegistry.register_source("bracket", bracket_datums)

    # Create mate constraint
    mate = Mate(
        name="servo_to_bracket",
        mate_type=MateType.COINCIDENT,
        part_a="bracket",
        datum_a="servo_mount_face",
        part_b="servo",
        datum_b="stator_face",
    )

    # Solve for child transform
    solver = MateConstraintSolver()
    result = solver.solve_mate(mate)

    if result.success:
        child_transform = result.transform  # 4x4 numpy array

Note:
    This module's MateMateSolveResult is distinct from yapcad.assembly.intent.MateSolveResult.

Mate constraint solver contributed by Jeremy Mika.

Copyright (c) 2026 yapCAD contributors
License: MIT
"""

from __future__ import annotations

import math
import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any, Tuple, Union

from .datum import Datum, DatumType
from .mate import Mate, MateType

logger = logging.getLogger(__name__)

# =============================================================================
# Optional Dependencies
# =============================================================================

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None

try:
    from .datum_registry import DatumRegistry
    HAS_REGISTRY = True
except ImportError:
    HAS_REGISTRY = False
    DatumRegistry = None

try:
    from .face_mate import compute_flush_transform, compute_concentric_transform
    HAS_FACE_MATE = True
except ImportError:
    HAS_FACE_MATE = False
    compute_flush_transform = None
    compute_concentric_transform = None


# =============================================================================
# Constants
# =============================================================================

# Mate axis remapping rotations (3x3 rotation matrices)
# These transform local Z axis to the specified mate axis
MATE_AXIS_ROTATIONS: Dict[str, Any] = {}

if HAS_NUMPY:
    MATE_AXIS_ROTATIONS = {
        "X": np.array([  # Ry(-90): Z -> X
            [0, 0, -1],
            [0, 1, 0],
            [1, 0, 0]
        ], dtype=float),
        "Y": np.array([  # Rx(90): Z -> Y
            [1, 0, 0],
            [0, 0, -1],
            [0, 1, 0]
        ], dtype=float),
        "Z": np.eye(3, dtype=float),  # Identity
        "-X": np.array([  # Ry(90): Z -> -X
            [0, 0, 1],
            [0, 1, 0],
            [-1, 0, 0]
        ], dtype=float),
        "-Y": np.array([  # Rx(-90): Z -> -Y
            [1, 0, 0],
            [0, 0, 1],
            [0, -1, 0]
        ], dtype=float),
        "-Z": np.array([  # Rx(180): Z -> -Z
            [1, 0, 0],
            [0, -1, 0],
            [0, 0, -1]
        ], dtype=float),
    }


# =============================================================================
# Result Types
# =============================================================================

@dataclass
class MateMateSolveResult:
    """Result of solving a mate constraint.

    Attributes:
        success: True if constraint was successfully solved
        transform: 4x4 homogeneous transform matrix (child in parent frame)
        error_message: Human-readable error description if failed
        residual: Constraint satisfaction error (0.0 = perfect)
        details: Additional solver information
    """
    success: bool
    transform: Optional[Any] = None  # np.ndarray
    error_message: str = ""
    residual: float = 0.0
    details: Dict[str, Any] = field(default_factory=dict)

    def __str__(self) -> str:
        status = "OK" if self.success else "FAIL"
        if self.success:
            return f"[{status}] residual={self.residual:.6f}"
        return f"[{status}] {self.error_message}"


@dataclass
class ValidationResult:
    """Result of validating a transform matrix.

    Attributes:
        valid: True if transform is valid (orthonormal, proper shape)
        is_rigid: True if rotation is orthonormal with det=+1
        is_orthonormal: True if rotation columns are unit length
        position_error: Maximum deviation from orthonormality
        orientation_error: Deviation from det=1
        error_messages: List of validation failure reasons
    """
    valid: bool
    is_rigid: bool = True
    is_orthonormal: bool = True
    position_error: float = 0.0
    orientation_error: float = 0.0
    error_messages: List[str] = field(default_factory=list)

    def __str__(self) -> str:
        if self.valid:
            return "[VALID] Transform is rigid and orthonormal"
        return f"[INVALID] {'; '.join(self.error_messages)}"


# =============================================================================
# MateConstraintSolver
# =============================================================================

class MateConstraintSolver:
    """Generic constraint solver for computing transforms from mates.

    This solver takes Mate constraint definitions and computes the 4x4
    homogeneous transform that positions the child part relative to the
    parent part such that the constraint is satisfied.

    Supported constraint types:
        - COINCIDENT: Face-to-face contact (normals opposed)
        - CONCENTRIC: Axes colinear
        - PARALLEL: Directions aligned
        - PERPENDICULAR: Directions at 90 degrees
        - DISTANCE: Fixed offset along normal
        - ANGLE: Fixed angular relationship
        - REVOLUTE: Rotation about axis (returns base transform)

    Example::

        solver = MateConstraintSolver()

        # Solve single mate
        result = solver.solve_mate(my_mate)

        # Solve assembly
        results = solver.solve_all([mate1, mate2, mate3], base_transforms)

    Attributes:
        cache: Dict mapping mate names to cached transform results
        tolerance: Position tolerance in mm (default 0.001)
        angle_tolerance: Angular tolerance in degrees (default 0.1)
    """

    def __init__(
        self,
        tolerance: float = 0.001,
        angle_tolerance: float = 0.1,
        use_cache: bool = True,
    ):
        """Initialize the constraint solver.

        :param tolerance: Position tolerance in mm
        :param angle_tolerance: Angular tolerance in degrees
        :param use_cache: Whether to cache solved transforms
        """
        if not HAS_NUMPY:
            raise RuntimeError(
                "numpy is required for MateConstraintSolver. "
                "Install with: pip install numpy"
            )

        self.tolerance = tolerance
        self.angle_tolerance = angle_tolerance
        self._use_cache = use_cache
        self._cache: Dict[str, MateMateSolveResult] = {}

    def solve_mate(
        self,
        mate: Mate,
        parent_world_transform: Optional[Any] = None,
        mate_axis: str = "Z",
    ) -> MateMateSolveResult:
        """Solve a mate constraint to compute the child transform.

        :param mate: Mate constraint definition
        :param parent_world_transform: Optional 4x4 world transform of parent
        :param mate_axis: Coordinate axis for mate ("X", "Y", "Z", "-Z")
        :returns: MateSolveResult with transform or error

        The returned transform positions the child part origin in the parent
        part's coordinate frame such that the mate constraint is satisfied.
        """
        # Check cache
        cache_key = f"{mate.name}:{mate_axis}"
        if self._use_cache and cache_key in self._cache:
            return self._cache[cache_key]

        try:
            # Get datums from registry
            datum_a = self._get_datum(mate.part_a, mate.datum_a)
            datum_b = self._get_datum(mate.part_b, mate.datum_b)

            if datum_a is None:
                return MateMateSolveResult(
                    success=False,
                    error_message=f"Datum not found: {mate.part_a}.{mate.datum_a}"
                )
            if datum_b is None:
                return MateMateSolveResult(
                    success=False,
                    error_message=f"Datum not found: {mate.part_b}.{mate.datum_b}"
                )

            # Dispatch to constraint-specific solver
            transform = self._solve_constraint(mate.mate_type, datum_a, datum_b, mate)

            # Apply mate axis rotation if needed
            if mate_axis != "Z":
                transform = self._apply_mate_axis_rotation(transform, mate_axis)

            # Apply parent world transform if provided
            if parent_world_transform is not None:
                transform = parent_world_transform @ transform

            result = MateMateSolveResult(
                success=True,
                transform=transform,
                residual=0.0,
                details={
                    "mate_type": mate.mate_type.value,
                    "mate_axis": mate_axis,
                }
            )

            # Cache result
            if self._use_cache:
                self._cache[cache_key] = result

            return result

        except Exception as e:
            logger.exception(f"Error solving mate {mate.name}")
            return MateMateSolveResult(
                success=False,
                error_message=str(e)
            )

    def solve_all(
        self,
        mates: List[Mate],
        parent_transforms: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, MateMateSolveResult]:
        """Solve multiple mate constraints.

        :param mates: List of Mate constraints to solve
        :param parent_transforms: Optional dict of {part_name: 4x4 transform}
        :returns: Dict mapping mate names to MateMateSolveResults
        """
        results = {}
        parent_transforms = parent_transforms or {}

        for mate in mates:
            parent_tf = parent_transforms.get(mate.part_a)
            results[mate.name] = self.solve_mate(mate, parent_tf)

        return results

    def validate_transform(
        self,
        transform: Any,
        tolerance_det: float = 1e-6,
    ) -> ValidationResult:
        """Validate a 4x4 transformation matrix.

        :param transform: 4x4 numpy array to validate
        :param tolerance_det: Tolerance for determinant check
        :returns: ValidationResult with validation details

        Checks:
            - Shape is (4, 4)
            - Bottom row is [0, 0, 0, 1]
            - Rotation part is orthonormal (R^T R = I)
            - Determinant is +1 (proper rotation, not reflection)
        """
        errors = []

        # Shape check
        if transform.shape != (4, 4):
            return ValidationResult(
                valid=False,
                error_messages=[f"Invalid shape: {transform.shape}, expected (4, 4)"]
            )

        # Bottom row check
        bottom = transform[3, :]
        if not np.allclose(bottom, [0, 0, 0, 1], atol=1e-10):
            errors.append(f"Invalid bottom row: {bottom}")

        # Extract rotation
        R = transform[:3, :3]

        # Orthonormality check: R^T @ R should be identity
        RtR = R.T @ R
        ortho_error = np.max(np.abs(RtR - np.eye(3)))
        is_orthonormal = ortho_error < 1e-6

        if not is_orthonormal:
            errors.append(f"Rotation not orthonormal: max error {ortho_error:.2e}")

        # Determinant check
        det = np.linalg.det(R)
        is_rigid = abs(det - 1.0) < tolerance_det

        if not is_rigid:
            errors.append(f"Determinant {det:.6f} != 1 (reflection or scaling)")

        return ValidationResult(
            valid=len(errors) == 0,
            is_rigid=is_rigid,
            is_orthonormal=is_orthonormal,
            orientation_error=abs(det - 1.0),
            position_error=ortho_error,
            error_messages=errors,
        )

    def invalidate_cache(self, mate_name: Optional[str] = None) -> None:
        """Invalidate cached transforms.

        :param mate_name: Specific mate to invalidate, or None for all
        """
        if mate_name is None:
            self._cache.clear()
        else:
            # Remove all cache entries starting with mate_name
            keys_to_remove = [k for k in self._cache if k.startswith(mate_name)]
            for key in keys_to_remove:
                del self._cache[key]

    # =========================================================================
    # Private: Constraint Solvers
    # =========================================================================

    def _solve_constraint(
        self,
        mate_type: MateType,
        datum_a: Datum,
        datum_b: Datum,
        mate: Mate,
    ) -> Any:
        """Dispatch to constraint-specific solver.

        :param mate_type: Type of constraint
        :param datum_a: Parent datum
        :param datum_b: Child datum
        :param mate: Full mate definition (for offset/angle params)
        :returns: 4x4 transform matrix
        """
        # Use face_mate functions if available
        if HAS_FACE_MATE:
            if mate_type in (MateType.COINCIDENT, MateType.RIGID):
                return compute_flush_transform(datum_a, datum_b)
            if mate_type == MateType.CONCENTRIC:
                return compute_concentric_transform(datum_a, datum_b)

        # Fallback implementations
        if mate_type in (MateType.COINCIDENT, MateType.RIGID):
            return self._solve_coincident(datum_a, datum_b)
        elif mate_type == MateType.CONCENTRIC:
            return self._solve_concentric(datum_a, datum_b)
        elif mate_type == MateType.PARALLEL:
            return self._solve_parallel(datum_a, datum_b)
        elif mate_type == MateType.PERPENDICULAR:
            return self._solve_perpendicular(datum_a, datum_b)
        elif mate_type == MateType.DISTANCE:
            return self._solve_distance(datum_a, datum_b, mate.offset)
        elif mate_type == MateType.ANGLE:
            return self._solve_angle(datum_a, datum_b, mate.angle)
        elif mate_type == MateType.REVOLUTE:
            return self._solve_revolute(datum_a, datum_b)
        else:
            raise ValueError(f"Unsupported mate type: {mate_type}")

    def _solve_coincident(self, datum_a: Datum, datum_b: Datum) -> Any:
        """Solve COINCIDENT constraint (face-to-face contact).

        The child is positioned so its datum origin coincides with the
        parent datum origin, and datum normals point in opposite directions.
        """
        # Get origins (homogeneous: [x, y, z, 1])
        origin_a = np.array(datum_a.origin[:3])
        origin_b = np.array(datum_b.origin[:3])

        # Get normals
        normal_a = self._get_datum_direction(datum_a)
        normal_b = self._get_datum_direction(datum_b)

        # Child normal should oppose parent normal
        target_normal = -normal_a

        # Compute rotation to align normal_b with target_normal
        R = self._rotation_between_vectors(normal_b, target_normal)

        # Build transform: rotate child, then translate to parent origin
        # The child datum origin moves to parent datum origin
        transform = np.eye(4)
        transform[:3, :3] = R

        # Translation: move rotated child datum origin to parent datum origin
        rotated_origin_b = R @ origin_b
        transform[:3, 3] = origin_a - rotated_origin_b

        return transform

    def _solve_concentric(self, datum_a: Datum, datum_b: Datum) -> Any:
        """Solve CONCENTRIC constraint (axes colinear)."""
        origin_a = np.array(datum_a.origin[:3])
        origin_b = np.array(datum_b.origin[:3])

        axis_a = self._get_datum_direction(datum_a)
        axis_b = self._get_datum_direction(datum_b)

        # Align axes (same direction)
        R = self._rotation_between_vectors(axis_b, axis_a)

        transform = np.eye(4)
        transform[:3, :3] = R

        # Position: project parent origin onto child axis, translate
        rotated_origin_b = R @ origin_b
        transform[:3, 3] = origin_a - rotated_origin_b

        return transform

    def _solve_parallel(self, datum_a: Datum, datum_b: Datum) -> Any:
        """Solve PARALLEL constraint (directions aligned)."""
        dir_a = self._get_datum_direction(datum_a)
        dir_b = self._get_datum_direction(datum_b)

        R = self._rotation_between_vectors(dir_b, dir_a)

        transform = np.eye(4)
        transform[:3, :3] = R
        # No translation constraint for PARALLEL
        return transform

    def _solve_perpendicular(self, datum_a: Datum, datum_b: Datum) -> Any:
        """Solve PERPENDICULAR constraint (90 degree angle)."""
        dir_a = self._get_datum_direction(datum_a)
        dir_b = self._get_datum_direction(datum_b)

        # Find perpendicular direction
        # Use cross product to find rotation axis
        cross = np.cross(dir_b, dir_a)
        if np.linalg.norm(cross) < 1e-10:
            # Already perpendicular or parallel - find arbitrary perpendicular
            if abs(dir_a[0]) < 0.9:
                perp = np.cross(dir_a, [1, 0, 0])
            else:
                perp = np.cross(dir_a, [0, 1, 0])
            perp = perp / np.linalg.norm(perp)
        else:
            perp = cross / np.linalg.norm(cross)

        # Rotate dir_b to be perpendicular to dir_a
        R = self._rotation_between_vectors(dir_b, perp)

        transform = np.eye(4)
        transform[:3, :3] = R
        return transform

    def _solve_distance(self, datum_a: Datum, datum_b: Datum, offset: float) -> Any:
        """Solve DISTANCE constraint (fixed offset along normal)."""
        origin_a = np.array(datum_a.origin[:3])
        normal_a = self._get_datum_direction(datum_a)

        # Position child at offset distance from parent along normal
        transform = np.eye(4)
        transform[:3, 3] = origin_a + normal_a * offset

        return transform

    def _solve_angle(
        self, datum_a: Datum, datum_b: Datum, angle_deg: float
    ) -> Any:
        """Solve ANGLE constraint (fixed angular relationship)."""
        axis_a = self._get_datum_direction(datum_a)

        # Rotate around parent axis by specified angle
        R = self._rodrigues_rotation(axis_a, math.radians(angle_deg))

        transform = np.eye(4)
        transform[:3, :3] = R
        return transform

    def _solve_revolute(self, datum_a: Datum, datum_b: Datum) -> Any:
        """Solve REVOLUTE joint (returns base transform, rotation is free)."""
        # For revolute joints, we solve as CONCENTRIC for the base transform
        return self._solve_concentric(datum_a, datum_b)

    # =========================================================================
    # Private: Helper Functions
    # =========================================================================

    def _get_datum(self, part_id: str, datum_name: str) -> Optional[Datum]:
        """Retrieve a datum from the registry.

        :param part_id: Part identifier (source ID in registry)
        :param datum_name: Name of datum on the part
        :returns: Datum object or None if not found
        """
        if not HAS_REGISTRY or DatumRegistry is None:
            logger.warning("DatumRegistry not available")
            return None

        return DatumRegistry.get_datum(part_id, datum_name)

    def _get_datum_direction(self, datum: Datum) -> Any:
        """Extract direction vector from a datum.

        For PLANE datums, returns the normal.
        For AXIS datums, returns the axis direction.
        For POINT datums, returns [0, 0, 1] (default Z).
        """
        if datum.normal is not None:
            n = np.array(datum.normal[:3])
            norm = np.linalg.norm(n)
            return n / norm if norm > 1e-10 else np.array([0.0, 0.0, 1.0])

        if datum.direction is not None:
            d = np.array(datum.direction[:3])
            norm = np.linalg.norm(d)
            return d / norm if norm > 1e-10 else np.array([0.0, 0.0, 1.0])

        # Default to Z axis
        return np.array([0.0, 0.0, 1.0])

    def _rotation_between_vectors(self, v1: Any, v2: Any) -> Any:
        """Compute rotation matrix that rotates v1 to v2.

        Uses Rodrigues' rotation formula.

        :param v1: Source unit vector
        :param v2: Target unit vector
        :returns: 3x3 rotation matrix
        """
        v1 = v1 / np.linalg.norm(v1)
        v2 = v2 / np.linalg.norm(v2)

        # Check for parallel vectors
        dot = np.dot(v1, v2)
        if dot > 0.9999:
            return np.eye(3)  # Already aligned
        if dot < -0.9999:
            # Opposite directions - rotate 180 around any perpendicular axis
            if abs(v1[0]) < 0.9:
                perp = np.cross(v1, [1, 0, 0])
            else:
                perp = np.cross(v1, [0, 1, 0])
            perp = perp / np.linalg.norm(perp)
            return self._rodrigues_rotation(perp, math.pi)

        # Rodrigues formula
        axis = np.cross(v1, v2)
        axis = axis / np.linalg.norm(axis)
        angle = math.acos(np.clip(dot, -1.0, 1.0))

        return self._rodrigues_rotation(axis, angle)

    def _rodrigues_rotation(self, axis: Any, angle: float) -> Any:
        """Compute rotation matrix from axis-angle using Rodrigues' formula.

        :param axis: Unit rotation axis
        :param angle: Rotation angle in radians
        :returns: 3x3 rotation matrix

        R = I + sin(θ)K + (1-cos(θ))K²
        where K is the skew-symmetric cross-product matrix of axis.
        """
        axis = axis / np.linalg.norm(axis)
        K = np.array([
            [0, -axis[2], axis[1]],
            [axis[2], 0, -axis[0]],
            [-axis[1], axis[0], 0]
        ])
        return np.eye(3) + math.sin(angle) * K + (1 - math.cos(angle)) * (K @ K)

    def _apply_mate_axis_rotation(self, transform: Any, mate_axis: str) -> Any:
        """Apply mate axis coordinate remapping.

        :param transform: 4x4 transform to modify
        :param mate_axis: Target axis ("X", "Y", "Z", "-X", "-Y", "-Z")
        :returns: Modified 4x4 transform
        """
        if mate_axis not in MATE_AXIS_ROTATIONS:
            logger.warning(f"Unknown mate_axis: {mate_axis}, using Z")
            return transform

        R_axis = MATE_AXIS_ROTATIONS[mate_axis]

        result = transform.copy()
        result[:3, :3] = transform[:3, :3] @ R_axis
        return result


# =============================================================================
# Convenience Functions
# =============================================================================

def apply_mate_axis_rotation(transform: Any, mate_axis: str) -> Any:
    """Apply mate axis coordinate remapping to a transform.

    :param transform: 4x4 numpy array
    :param mate_axis: Target axis ("X", "Y", "Z", "-X", "-Y", "-Z")
    :returns: Modified 4x4 transform

    Example::

        # Remap Z axis to Y axis
        new_tf = apply_mate_axis_rotation(transform, "Y")
    """
    if mate_axis not in MATE_AXIS_ROTATIONS:
        return transform

    R_axis = MATE_AXIS_ROTATIONS[mate_axis]
    result = transform.copy()
    result[:3, :3] = transform[:3, :3] @ R_axis
    return result


def solve_mate_chain(
    mates: List[Mate],
    base_transform: Optional[Any] = None,
    solver: Optional[MateConstraintSolver] = None,
) -> Dict[str, Any]:
    """Solve a chain of mate constraints sequentially.

    :param mates: Ordered list of mates (parent -> child order)
    :param base_transform: Optional world transform for first parent
    :param solver: Optional solver instance (creates one if not provided)
    :returns: Dict mapping part names to world transforms

    Example::

        world_transforms = solve_mate_chain(
            [base_to_link1, link1_to_link2, link2_to_tool],
            base_transform=np.eye(4)
        )
    """
    if solver is None:
        solver = MateConstraintSolver()

    if base_transform is None:
        base_transform = np.eye(4)

    transforms: Dict[str, Any] = {}
    current_transform = base_transform

    for mate in mates:
        result = solver.solve_mate(mate, parent_world_transform=current_transform)
        if result.success:
            transforms[mate.part_b] = result.transform
            current_transform = result.transform
        else:
            logger.error(f"Failed to solve mate {mate.name}: {result.error_message}")
            break

    return transforms


# =============================================================================
# Exports
# =============================================================================

__all__ = [
    # Classes
    "MateConstraintSolver",
    "MateMateSolveResult",
    "ValidationResult",
    # Constants
    "MATE_AXIS_ROTATIONS",
    # Functions
    "apply_mate_axis_rotation",
    "solve_mate_chain",
]
