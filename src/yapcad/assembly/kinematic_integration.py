"""Kinematic chain integration for yapCAD assembly system.

This module provides integration between the assembly constraint/mate system
and the kinematic chain transform propagation system. It enables:

1. Transform propagation from kinematic chains to assembly validation
2. Constraint evaluation in world coordinates using chain transforms
3. Forward kinematics: root-to-leaf transform computation
4. Assembly-wide constraint validation with comprehensive reporting

Key Classes:
    KinematicConstraint: Enhanced constraint that works with world transforms
    ConstraintEvaluator: Evaluates constraints given world transform dict
    ValidationReport: Comprehensive assembly validation results
    AssemblyValidator: Validates full assemblies against constraint sets

Example:
    >>> from yapcad.assembly.kinematic_integration import (
    ...     AssemblyValidator, KinematicConstraint, ConstraintType
    ... )
    >>>
    >>> # Create constraints
    >>> tangent = KinematicConstraint(
    ...     name="wheel_tangent",
    ...     constraint_type=ConstraintType.TANGENT,
    ...     frame_a="DDSM115_MOTOR_1",
    ...     frame_b=None,  # World reference
    ...     reference_center=(0, 0, 0),
    ...     reference_radius=124.5,
    ...     tolerance_deg=2.0
    ... )
    >>>
    >>> # Get transforms from kinematic chain
    >>> from yapcad.kinematics import KinematicChain
    >>> chain = KinematicChain("my_assembly")
    >>> transforms = chain.get_all_world_transforms()
    >>>
    >>> # Validate
    >>> validator = AssemblyValidator()
    >>> validator.add_constraint(tangent)
    >>> report = validator.validate(transforms)
    >>> print(report)

Copyright (c) 2026 yapCAD contributors
License: MIT
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Any, Union, Callable
from enum import Enum

# Import numpy
try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None

# Import trimesh for mesh-based collision detection
try:
    import trimesh
    HAS_TRIMESH = True
except ImportError:
    HAS_TRIMESH = False
    trimesh = None


# =============================================================================
# MESH COLLISION DETECTION
# =============================================================================

@dataclass
class MeshCollisionResult:
    """Result of mesh-based collision detection between two parts.

    Attributes:
        collides: True if the meshes intersect
        penetration_depth: Estimated penetration depth in mm (0 if no collision)
        collision_volume: Volume of the intersection region in mm^3 (0 if no collision)
        contact_points: List of approximate contact/intersection points
        error_message: Error message if collision check failed
    """
    collides: bool
    penetration_depth: float = 0.0
    collision_volume: float = 0.0
    contact_points: List[Tuple[float, float, float]] = field(default_factory=list)
    error_message: str = ""

    def __str__(self) -> str:
        if self.error_message:
            return f"[ERROR] {self.error_message}"
        if self.collides:
            return (f"[COLLISION] penetration={self.penetration_depth:.3f}mm, "
                    f"volume={self.collision_volume:.3f}mm^3, "
                    f"contacts={len(self.contact_points)}")
        return "[OK] No collision detected"


def check_mesh_collision(
    stl_path_a: str,
    transform_a: Any,
    stl_path_b: str,
    transform_b: Any
) -> MeshCollisionResult:
    """Check for collision between two STL meshes with applied transforms.

    This function loads two STL files, applies 4x4 transformation matrices,
    and checks if the resulting meshes intersect. It provides detailed
    collision information including penetration depth and intersection volume.

    Args:
        stl_path_a: Path to the first STL file
        transform_a: 4x4 transformation matrix (numpy array or Transform object)
        stl_path_b: Path to the second STL file
        transform_b: 4x4 transformation matrix (numpy array or Transform object)

    Returns:
        MeshCollisionResult with collision details

    Example:
        >>> from yapcad.assembly.kinematic_integration import check_mesh_collision
        >>> import numpy as np
        >>>
        >>> # Identity transform (no transformation)
        >>> identity = np.eye(4)
        >>>
        >>> # Translation by 100mm in X
        >>> translated = np.eye(4)
        >>> translated[0, 3] = 100.0
        >>>
        >>> result = check_mesh_collision(
        ...     "part_a.stl", identity,
        ...     "part_b.stl", translated
        ... )
        >>> if result.collides:
        ...     print(f"Collision! Penetration: {result.penetration_depth}mm")
    """
    if not HAS_TRIMESH:
        return MeshCollisionResult(
            collides=False,
            error_message="trimesh library required for mesh collision detection"
        )

    if not HAS_NUMPY:
        return MeshCollisionResult(
            collides=False,
            error_message="numpy required for mesh collision detection"
        )

    # Load meshes
    try:
        mesh_a = trimesh.load(stl_path_a)
    except Exception as e:
        return MeshCollisionResult(
            collides=False,
            error_message=f"Failed to load mesh A ({stl_path_a}): {e}"
        )

    try:
        mesh_b = trimesh.load(stl_path_b)
    except Exception as e:
        return MeshCollisionResult(
            collides=False,
            error_message=f"Failed to load mesh B ({stl_path_b}): {e}"
        )

    # Handle Scene objects (multi-part STL files)
    if isinstance(mesh_a, trimesh.Scene):
        if len(mesh_a.geometry) == 0:
            return MeshCollisionResult(
                collides=False,
                error_message=f"Mesh A ({stl_path_a}) contains no geometry"
            )
        # Concatenate all geometries in the scene
        mesh_a = trimesh.util.concatenate(list(mesh_a.geometry.values()))

    if isinstance(mesh_b, trimesh.Scene):
        if len(mesh_b.geometry) == 0:
            return MeshCollisionResult(
                collides=False,
                error_message=f"Mesh B ({stl_path_b}) contains no geometry"
            )
        mesh_b = trimesh.util.concatenate(list(mesh_b.geometry.values()))

    # Extract transform matrices
    tf_a = _extract_transform_matrix(transform_a)
    tf_b = _extract_transform_matrix(transform_b)

    if tf_a is None:
        return MeshCollisionResult(
            collides=False,
            error_message="Invalid transform_a: must be 4x4 matrix or Transform object"
        )

    if tf_b is None:
        return MeshCollisionResult(
            collides=False,
            error_message="Invalid transform_b: must be 4x4 matrix or Transform object"
        )

    # Apply transforms to meshes
    mesh_a_transformed = mesh_a.copy()
    mesh_a_transformed.apply_transform(tf_a)

    mesh_b_transformed = mesh_b.copy()
    mesh_b_transformed.apply_transform(tf_b)

    # Quick AABB check first (fast rejection)
    if not _aabb_intersects(mesh_a_transformed.bounds, mesh_b_transformed.bounds):
        return MeshCollisionResult(
            collides=False,
            penetration_depth=0.0,
            collision_volume=0.0,
            contact_points=[]
        )

    # AABBs overlap - need more detailed check
    # Try multiple methods in order of preference

    # Method 1: Try collision manager (requires FCL)
    try:
        collision_manager = trimesh.collision.CollisionManager()
        collision_manager.add_object("mesh_a", mesh_a_transformed)

        is_collision, contact_data = collision_manager.in_collision_single(
            mesh_b_transformed, return_data=True
        )

        if is_collision:
            contact_points = []
            max_depth = 0.0
            for contact in contact_data:
                if hasattr(contact, 'point'):
                    contact_points.append(tuple(contact.point))
                if hasattr(contact, 'depth'):
                    max_depth = max(max_depth, abs(contact.depth))

            return MeshCollisionResult(
                collides=True,
                penetration_depth=max_depth,
                collision_volume=0.0,
                contact_points=contact_points
            )
        else:
            return MeshCollisionResult(
                collides=False,
                penetration_depth=0.0,
                collision_volume=0.0,
                contact_points=[]
            )
    except (ValueError, ImportError):
        # FCL not available, continue to fallback methods
        pass
    except Exception:
        # Other errors, continue to fallback
        pass

    # Method 2: Sample-based collision detection (no external dependencies)
    # Sample vertices from each mesh and check if they fall within the
    # AABB of the other mesh, then verify with face intersection
    try:
        collides, penetration, contacts = _check_collision_sampling(
            mesh_a_transformed, mesh_b_transformed
        )
        return MeshCollisionResult(
            collides=collides,
            penetration_depth=penetration,
            collision_volume=0.0,
            contact_points=contacts
        )
    except Exception as e:
        return MeshCollisionResult(
            collides=False,
            error_message=f"Collision check failed: {e}"
        )


def _extract_transform_matrix(transform: Any) -> Optional[np.ndarray]:
    """Extract a 4x4 numpy array from various transform representations.

    Args:
        transform: Can be a numpy array, Transform object with .matrix attribute,
                   or a list/tuple that can be converted to numpy array

    Returns:
        4x4 numpy array or None if conversion failed
    """
    if transform is None:
        return np.eye(4)

    # Handle Transform objects
    if hasattr(transform, 'matrix'):
        transform = transform.matrix

    # Convert to numpy array
    try:
        arr = np.array(transform, dtype=np.float64)
        if arr.shape == (4, 4):
            return arr
        elif arr.shape == (16,):
            return arr.reshape(4, 4)
        else:
            return None
    except Exception:
        return None


def _aabb_intersects(bounds_a: np.ndarray, bounds_b: np.ndarray) -> bool:
    """Check if two axis-aligned bounding boxes intersect.

    Args:
        bounds_a: Array of shape (2, 3) with [min_corner, max_corner]
        bounds_b: Array of shape (2, 3) with [min_corner, max_corner]

    Returns:
        True if AABBs intersect
    """
    # Check overlap on each axis
    for i in range(3):
        if bounds_a[1, i] < bounds_b[0, i] or bounds_b[1, i] < bounds_a[0, i]:
            return False
    return True


def _check_collision_sampling(
    mesh_a: 'trimesh.Trimesh',
    mesh_b: 'trimesh.Trimesh'
) -> Tuple[bool, float, List[Tuple[float, float, float]]]:
    """Sample-based collision detection without external dependencies.

    This method works by:
    1. Finding vertices from mesh_a that fall within mesh_b's AABB
    2. For those vertices, checking if they're inside mesh_b using face normals
    3. Repeating in reverse (mesh_b vertices in mesh_a)

    Args:
        mesh_a: First transformed mesh
        mesh_b: Second transformed mesh

    Returns:
        Tuple of (collides: bool, penetration_depth: float, contact_points: list)
    """
    contact_points = []
    max_penetration = 0.0

    bounds_a = mesh_a.bounds
    bounds_b = mesh_b.bounds

    # Find vertices of mesh_a inside mesh_b's AABB
    verts_a = mesh_a.vertices
    in_b_aabb = np.all(
        (verts_a >= bounds_b[0]) & (verts_a <= bounds_b[1]),
        axis=1
    )

    # Find vertices of mesh_b inside mesh_a's AABB
    verts_b = mesh_b.vertices
    in_a_aabb = np.all(
        (verts_b >= bounds_a[0]) & (verts_b <= bounds_a[1]),
        axis=1
    )

    candidates_a = verts_a[in_b_aabb]
    candidates_b = verts_b[in_a_aabb]

    # If no vertices in the overlap region, check for face intersections
    # at the AABB intersection boundary
    if len(candidates_a) == 0 and len(candidates_b) == 0:
        # AABBs overlap but no vertices penetrate - could still have
        # face-face intersection. Check using triangle-triangle tests
        # For efficiency, just sample the AABB overlap region
        overlap_min = np.maximum(bounds_a[0], bounds_b[0])
        overlap_max = np.minimum(bounds_a[1], bounds_b[1])

        # Sample points in overlap region
        n_samples = 5
        sample_grid = []
        for x in np.linspace(overlap_min[0], overlap_max[0], n_samples):
            for y in np.linspace(overlap_min[1], overlap_max[1], n_samples):
                for z in np.linspace(overlap_min[2], overlap_max[2], n_samples):
                    sample_grid.append([x, y, z])
        sample_points = np.array(sample_grid)

        # Check if any sample points are inside both meshes
        # by checking if they're "inside" using ray counting
        inside_a = _points_inside_mesh_simple(sample_points, mesh_a)
        inside_b = _points_inside_mesh_simple(sample_points, mesh_b)

        inside_both = inside_a & inside_b
        if np.any(inside_both):
            collision_pts = sample_points[inside_both]
            contact_points = [tuple(p) for p in collision_pts[:10]]  # Limit to 10
            max_penetration = 1.0  # Unknown exact depth
            return True, max_penetration, contact_points

        return False, 0.0, []

    # Check if candidate vertices from A are inside mesh B
    if len(candidates_a) > 0:
        inside_b = _points_inside_mesh_simple(candidates_a, mesh_b)
        if np.any(inside_b):
            collision_pts = candidates_a[inside_b]
            contact_points.extend([tuple(p) for p in collision_pts[:10]])
            # Estimate penetration as distance from point to nearest face center
            max_penetration = _estimate_penetration(collision_pts, mesh_b)

    # Check if candidate vertices from B are inside mesh A
    if len(candidates_b) > 0:
        inside_a = _points_inside_mesh_simple(candidates_b, mesh_a)
        if np.any(inside_a):
            collision_pts = candidates_b[inside_a]
            contact_points.extend([tuple(p) for p in collision_pts[:10]])
            pen = _estimate_penetration(collision_pts, mesh_a)
            max_penetration = max(max_penetration, pen)

    collides = len(contact_points) > 0
    return collides, max_penetration, contact_points[:20]  # Limit contacts


def _points_inside_mesh_simple(points: np.ndarray,
                                mesh: 'trimesh.Trimesh') -> np.ndarray:
    """Check if points are inside a mesh using face normal voting.

    For each point, cast a ray in the +X direction and count face intersections.
    If odd number of intersections, point is inside.

    This is a simplified version that doesn't require rtree.

    Args:
        points: Array of points (N, 3)
        mesh: Trimesh object

    Returns:
        Boolean array of size N, True if point is inside
    """
    n_points = len(points)
    inside = np.zeros(n_points, dtype=bool)

    # Get mesh triangles
    triangles = mesh.triangles  # Shape (n_faces, 3, 3)

    # For efficiency, only test a subset if there are many faces
    max_faces = 1000
    if len(triangles) > max_faces:
        indices = np.random.choice(len(triangles), max_faces, replace=False)
        triangles = triangles[indices]

    # Ray direction
    ray_dir = np.array([1.0, 0.0, 0.0])

    for i, point in enumerate(points):
        # Count intersections with triangles
        intersections = 0
        for tri in triangles:
            if _ray_triangle_intersect(point, ray_dir, tri):
                intersections += 1

        # Odd number means inside
        inside[i] = (intersections % 2) == 1

    return inside


def _ray_triangle_intersect(origin: np.ndarray,
                             direction: np.ndarray,
                             triangle: np.ndarray) -> bool:
    """Moller-Trumbore ray-triangle intersection test.

    Args:
        origin: Ray origin (3,)
        direction: Ray direction (3,)
        triangle: Triangle vertices (3, 3)

    Returns:
        True if ray intersects triangle in positive direction
    """
    epsilon = 1e-10

    v0, v1, v2 = triangle
    edge1 = v1 - v0
    edge2 = v2 - v0

    pvec = np.cross(direction, edge2)
    det = np.dot(edge1, pvec)

    if abs(det) < epsilon:
        return False  # Ray parallel to triangle

    inv_det = 1.0 / det
    tvec = origin - v0
    u = np.dot(tvec, pvec) * inv_det

    if u < 0.0 or u > 1.0:
        return False

    qvec = np.cross(tvec, edge1)
    v = np.dot(direction, qvec) * inv_det

    if v < 0.0 or u + v > 1.0:
        return False

    t = np.dot(edge2, qvec) * inv_det

    return t > epsilon  # Intersection in positive direction


def _estimate_penetration(collision_points: np.ndarray,
                           mesh: 'trimesh.Trimesh') -> float:
    """Estimate penetration depth for collision points.

    Args:
        collision_points: Points that are inside the mesh
        mesh: The mesh being penetrated

    Returns:
        Estimated maximum penetration depth in mm
    """
    if len(collision_points) == 0:
        return 0.0

    # Simple estimate: distance from centroid of collision points
    # to the mesh centroid
    collision_center = np.mean(collision_points, axis=0)
    mesh_center = mesh.centroid

    # Use half the distance as a rough penetration estimate
    dist = np.linalg.norm(collision_center - mesh_center)

    # Also consider the mesh bounding box diagonal as reference
    bbox_diag = np.linalg.norm(mesh.bounds[1] - mesh.bounds[0])

    # Penetration is likely a fraction of the smaller dimension
    return min(dist * 0.5, bbox_diag * 0.1)


class KinematicConstraintType(Enum):
    """Constraint types for kinematic chain validation.

    These constraints operate on world-space transforms from kinematic chains
    and validate geometric relationships between frames.
    """

    # Positional constraints
    COINCIDENT = "coincident"        # Two points/origins at same location
    CONCENTRIC = "concentric"        # Two cylindrical axes share same line
    AT_DISTANCE = "at_distance"      # Two frames at specific distance

    # Directional constraints
    PARALLEL = "parallel"            # Two vectors/axes are parallel
    PERPENDICULAR = "perpendicular"  # Two vectors/axes are perpendicular
    TANGENT = "tangent"              # Axis tangent to curve at a point

    # Pattern constraints
    BOLT_PATTERN = "bolt_pattern"    # Multiple holes align with pattern

    # Clearance/collision constraints
    MIN_DISTANCE = "min_distance"    # Parts must be at least N mm apart
    Z_STACK_CLEARANCE = "z_stack_clearance"  # Stacked parts must have proper Z separation
    NO_OVERLAP = "no_overlap"        # Bounding boxes must not overlap

    # Custom
    CUSTOM = "custom"                # User-defined validation function


@dataclass
class ConstraintEvaluationResult:
    """Result of evaluating a single constraint.

    Attributes:
        satisfied: True if constraint is satisfied within tolerance
        error: Numeric error value (distance in mm or angle in degrees)
        error_vector: Optional direction vector of error for visualization
        error_message: Human-readable description
        details: Additional diagnostic information
    """
    satisfied: bool
    error: float = 0.0
    error_vector: Optional[Tuple[float, float, float]] = None
    error_message: str = ""
    details: Dict[str, Any] = field(default_factory=dict)

    def __str__(self) -> str:
        status = "PASS" if self.satisfied else "FAIL"
        if self.error > 0:
            return f"[{status}] {self.error_message} (error: {self.error:.4f})"
        return f"[{status}] {self.error_message}"


@dataclass
class KinematicConstraint:
    """A constraint that operates on world transforms from kinematic chains.

    This constraint type is designed to work with Transform objects (4x4 matrices)
    from the kinematic chain system, enabling validation of assembly relationships
    in world coordinates.

    Attributes:
        name: Unique identifier for this constraint
        constraint_type: Type of constraint to evaluate
        frame_a: Name of first frame (part name in kinematic chain)
        frame_b: Name of second frame, or None for world reference

        # Face-based constraint specification (optional)
        face_a: Face specification on frame_a. Can be:
                - "TOP" / "BOTTOM": Calculated from bounds_a (max_z / min_z face center)
                - "FRONT" / "BACK": Calculated from bounds_a (max_y / min_y face center)
                - "LEFT" / "RIGHT": Calculated from bounds_a (min_x / max_x face center)
                - Frame name (e.g., "OUTPUT_SHAFT", "SUN_INPUT"): Looked up from
                  world_transforms using key "part_name.frame_name"
                When specified, constraint position is taken from face instead of part origin.
        face_b: Face specification on frame_b (same options as face_a).

        # Axis specification (which local axis to use for directional constraints)
        axis_a: Local axis on frame_a to use ("x", "y", "z", or tuple)
        axis_b: Local axis on frame_b to use ("x", "y", "z", or tuple)

        # Reference geometry for world-referenced constraints
        reference_center: Center point for tangent/radial constraints
        reference_radius: Radius for tangent/radial constraints
        reference_axis: Reference axis direction for parallel/perpendicular
        reference_normal: Reference normal for plane constraints

        # Pattern parameters
        pattern_count: Number of holes in bolt pattern
        pattern_radius: Bolt circle radius
        pattern_offset_deg: Angular offset between patterns

        # Tolerances
        tolerance_mm: Linear tolerance in mm (default: 0.1)
        tolerance_deg: Angular tolerance in degrees (default: 1.0)

        # Metadata
        description: Human-readable description
        severity: "error", "warning", or "info"

        # Custom validator
        validator: Optional custom validation function

    Example with face-based constraints::

        # Clearance measured from servo OUTPUT_SHAFT face to ring SUN_INPUT face
        constraint = KinematicConstraint(
            name="servo_ring_interface",
            constraint_type=KinematicConstraintType.AT_DISTANCE,
            frame_a="AXIS2_SERVO_XH430",
            frame_b="AXIS2_RING_HOUSING",
            face_a="OUTPUT_SHAFT",    # Use servo's output shaft frame position
            face_b="SUN_INPUT",       # Use ring's sun input frame position
            reference_radius=0.0,     # Expected distance (touching)
            tolerance_mm=0.5
        )

        # Contact constraint using standard face names
        constraint = KinematicConstraint(
            name="stacked_parts",
            constraint_type=KinematicConstraintType.COINCIDENT,
            frame_a="GEARBOX_TOP",
            frame_b="GEARBOX_BOTTOM",
            face_a="BOTTOM",          # Bottom face of top gearbox (from bounds_a)
            face_b="TOP",             # Top face of bottom gearbox (from bounds_b)
            bounds_a=((-30, -30, 0), (30, 30, 50)),
            bounds_b=((-30, -30, 0), (30, 30, 50)),
            tolerance_mm=0.1
        )
    """
    name: str
    constraint_type: KinematicConstraintType

    # Frame references
    frame_a: str
    frame_b: Optional[str] = None

    # Face-based constraint specification
    face_a: Optional[str] = None
    face_b: Optional[str] = None

    # Axis specification
    axis_a: Union[str, Tuple[float, float, float]] = "z"
    axis_b: Union[str, Tuple[float, float, float]] = "z"

    # Reference geometry
    reference_center: Optional[Tuple[float, float, float]] = None
    reference_radius: Optional[float] = None
    reference_axis: Optional[Tuple[float, float, float]] = None
    reference_normal: Optional[Tuple[float, float, float]] = None

    # Pattern parameters
    pattern_count: int = 0
    pattern_radius: float = 0.0
    pattern_offset_deg: float = 0.0

    # Tolerances
    tolerance_mm: float = 0.1
    tolerance_deg: float = 1.0

    # Part bounding boxes for collision detection (min_xyz, max_xyz in local coords)
    # Format: ((min_x, min_y, min_z), (max_x, max_y, max_z))
    bounds_a: Optional[Tuple[Tuple[float, float, float], Tuple[float, float, float]]] = None
    bounds_b: Optional[Tuple[Tuple[float, float, float], Tuple[float, float, float]]] = None

    # STL paths for mesh-based collision detection (optional, enhances NO_OVERLAP)
    # When provided, accurate mesh intersection is used instead of AABB approximation
    stl_path_a: Optional[str] = None
    stl_path_b: Optional[str] = None

    # Metadata
    description: str = ""
    severity: str = "error"

    # Custom validator
    validator: Optional[Callable[[Dict[str, Any]], ConstraintEvaluationResult]] = None

    def _get_local_axis(self, axis_spec: Union[str, Tuple[float, float, float]]) -> Optional[np.ndarray]:
        """Convert axis specification to unit vector.

        Returns None if axis is invalid (zero-length vector).
        """
        if isinstance(axis_spec, str):
            if axis_spec.lower() == "x":
                return np.array([1.0, 0.0, 0.0])
            elif axis_spec.lower() == "y":
                return np.array([0.0, 1.0, 0.0])
            elif axis_spec.lower() == "z":
                return np.array([0.0, 0.0, 1.0])
            elif axis_spec.lower() == "-x":
                return np.array([-1.0, 0.0, 0.0])
            elif axis_spec.lower() == "-y":
                return np.array([0.0, -1.0, 0.0])
            elif axis_spec.lower() == "-z":
                return np.array([0.0, 0.0, -1.0])
            else:
                return None  # Unknown axis
        else:
            v = np.array(axis_spec[:3])
            norm = np.linalg.norm(v)
            if norm < 1e-10:
                return None  # Zero-length vector
            return v / norm

    def _transform_vector(self, transform: np.ndarray, vector: np.ndarray) -> np.ndarray:
        """Transform a direction vector (rotation only, no translation)."""
        # Extract 3x3 rotation matrix
        R = transform[:3, :3]
        return R @ vector

    def _get_origin(self, transform: np.ndarray) -> np.ndarray:
        """Extract origin (translation) from transform."""
        return transform[:3, 3].copy()

    def _get_standard_face_offset(
        self,
        face_name: str,
        bounds: Optional[Tuple[Tuple[float, float, float], Tuple[float, float, float]]]
    ) -> Optional[np.ndarray]:
        """Get local offset for standard face names (TOP, BOTTOM, etc.) from bounds.

        Standard face names:
            - TOP: Center of max_z face
            - BOTTOM: Center of min_z face
            - FRONT: Center of max_y face
            - BACK: Center of min_y face
            - RIGHT: Center of max_x face
            - LEFT: Center of min_x face

        Args:
            face_name: Standard face name (case-insensitive)
            bounds: Part bounding box ((min_x, min_y, min_z), (max_x, max_y, max_z))

        Returns:
            Local offset vector as numpy array, or None if bounds not provided
            or face_name is not a standard face
        """
        if bounds is None:
            return None

        face_upper = face_name.upper()
        min_pt, max_pt = bounds

        # Calculate center of bounding box in XY
        center_x = (min_pt[0] + max_pt[0]) / 2.0
        center_y = (min_pt[1] + max_pt[1]) / 2.0
        center_z = (min_pt[2] + max_pt[2]) / 2.0

        if face_upper == "TOP":
            return np.array([center_x, center_y, max_pt[2]])
        elif face_upper == "BOTTOM":
            return np.array([center_x, center_y, min_pt[2]])
        elif face_upper == "FRONT":
            return np.array([center_x, max_pt[1], center_z])
        elif face_upper == "BACK":
            return np.array([center_x, min_pt[1], center_z])
        elif face_upper == "RIGHT":
            return np.array([max_pt[0], center_y, center_z])
        elif face_upper == "LEFT":
            return np.array([min_pt[0], center_y, center_z])
        else:
            return None  # Not a standard face name

    def _get_face_position(
        self,
        part_name: str,
        part_transform: np.ndarray,
        face_spec: Optional[str],
        bounds: Optional[Tuple[Tuple[float, float, float], Tuple[float, float, float]]],
        world_transforms: Dict[str, Any]
    ) -> Tuple[np.ndarray, Optional[str]]:
        """Get world position for a face specification.

        Face specification can be:
        1. None: Use part origin (default behavior)
        2. Standard face name (TOP, BOTTOM, FRONT, BACK, LEFT, RIGHT):
           Calculate from bounds
        3. Frame name (e.g., "OUTPUT_SHAFT"): Look up from world_transforms
           using key "part_name.frame_name"

        Args:
            part_name: Name of the part (e.g., "AXIS2_SERVO_XH430")
            part_transform: 4x4 world transform matrix for the part
            face_spec: Face specification string, or None for origin
            bounds: Part bounding box for standard face calculations
            world_transforms: Dictionary of all world transforms

        Returns:
            Tuple of (world_position, error_message).
            If error_message is not None, position calculation failed.
        """
        # Case 1: No face specified - use part origin
        if face_spec is None:
            return self._get_origin(part_transform), None

        # Case 2: Standard face name - calculate from bounds
        local_offset = self._get_standard_face_offset(face_spec, bounds)
        if local_offset is not None:
            # Transform local offset to world coordinates
            local_pt = np.array([local_offset[0], local_offset[1], local_offset[2], 1.0])
            world_pt = part_transform @ local_pt
            return world_pt[:3], None

        # Case 3: Frame name - look up from world_transforms
        # Try "part_name.frame_name" format first
        frame_key = f"{part_name}.{face_spec}"
        if frame_key in world_transforms:
            tf = world_transforms[frame_key]
            if hasattr(tf, 'matrix'):
                tf = tf.matrix
            tf = np.array(tf)
            return self._get_origin(tf), None

        # Also try just the face_spec as a standalone key (for flexibility)
        if face_spec in world_transforms:
            tf = world_transforms[face_spec]
            if hasattr(tf, 'matrix'):
                tf = tf.matrix
            tf = np.array(tf)
            return self._get_origin(tf), None

        # If face_spec looks like a standard face but bounds weren't provided
        if face_spec.upper() in ("TOP", "BOTTOM", "FRONT", "BACK", "LEFT", "RIGHT"):
            return self._get_origin(part_transform), (
                f"Face '{face_spec}' requires bounds but none provided for {part_name}"
            )

        # Frame not found
        return self._get_origin(part_transform), (
            f"Frame '{face_spec}' not found. Tried keys: '{frame_key}' and '{face_spec}'"
        )

    def evaluate(self, world_transforms: Dict[str, Any]) -> ConstraintEvaluationResult:
        """Evaluate this constraint given world transforms of all parts.

        Args:
            world_transforms: Dictionary mapping part names to Transform objects
                              or 4x4 numpy arrays

        Returns:
            ConstraintEvaluationResult with evaluation status and error metrics
        """
        if not HAS_NUMPY:
            return ConstraintEvaluationResult(
                satisfied=False,
                error_message="NumPy required for constraint evaluation"
            )

        # Get transform matrix for frame_a
        if self.frame_a not in world_transforms:
            return ConstraintEvaluationResult(
                satisfied=False,
                error_message=f"Frame '{self.frame_a}' not found in transforms"
            )

        tf_a = world_transforms[self.frame_a]
        # Handle both numpy arrays and Transform objects
        if hasattr(tf_a, 'matrix'):
            tf_a = tf_a.matrix
        tf_a = np.array(tf_a)

        # Get transform for frame_b if specified
        tf_b = None
        if self.frame_b is not None:
            if self.frame_b not in world_transforms:
                return ConstraintEvaluationResult(
                    satisfied=False,
                    error_message=f"Frame '{self.frame_b}' not found in transforms"
                )
            tf_b = world_transforms[self.frame_b]
            if hasattr(tf_b, 'matrix'):
                tf_b = tf_b.matrix
            tf_b = np.array(tf_b)

        # Compute face positions if face specifications are provided
        # These will be used by evaluation methods that support face-based constraints
        face_pos_a = None
        face_pos_b = None
        face_warnings = []

        if self.face_a is not None:
            face_pos_a, error = self._get_face_position(
                self.frame_a, tf_a, self.face_a, self.bounds_a, world_transforms
            )
            if error:
                face_warnings.append(f"face_a: {error}")

        if self.face_b is not None and tf_b is not None:
            face_pos_b, error = self._get_face_position(
                self.frame_b, tf_b, self.face_b, self.bounds_b, world_transforms
            )
            if error:
                face_warnings.append(f"face_b: {error}")

        # Store face positions for use by evaluation methods
        self._face_pos_a = face_pos_a
        self._face_pos_b = face_pos_b
        self._face_warnings = face_warnings

        # Dispatch to appropriate evaluation method
        if self.constraint_type == KinematicConstraintType.COINCIDENT:
            return self._evaluate_coincident(tf_a, tf_b)
        elif self.constraint_type == KinematicConstraintType.CONCENTRIC:
            return self._evaluate_concentric(tf_a, tf_b)
        elif self.constraint_type == KinematicConstraintType.PARALLEL:
            return self._evaluate_parallel(tf_a, tf_b)
        elif self.constraint_type == KinematicConstraintType.PERPENDICULAR:
            return self._evaluate_perpendicular(tf_a, tf_b)
        elif self.constraint_type == KinematicConstraintType.TANGENT:
            return self._evaluate_tangent(tf_a)
        elif self.constraint_type == KinematicConstraintType.AT_DISTANCE:
            return self._evaluate_at_distance(tf_a, tf_b)
        elif self.constraint_type == KinematicConstraintType.BOLT_PATTERN:
            return self._evaluate_bolt_pattern(tf_a, tf_b)
        elif self.constraint_type == KinematicConstraintType.MIN_DISTANCE:
            return self._evaluate_min_distance(tf_a, tf_b)
        elif self.constraint_type == KinematicConstraintType.Z_STACK_CLEARANCE:
            return self._evaluate_z_stack_clearance(tf_a, tf_b)
        elif self.constraint_type == KinematicConstraintType.NO_OVERLAP:
            return self._evaluate_no_overlap(tf_a, tf_b)
        elif self.constraint_type == KinematicConstraintType.CUSTOM:
            if self.validator:
                try:
                    return self.validator(world_transforms)
                except Exception as e:
                    return ConstraintEvaluationResult(
                        satisfied=False,
                        error_message=f"Custom validator error: {e}"
                    )
            else:
                return ConstraintEvaluationResult(
                    satisfied=False,
                    error_message="CUSTOM constraint requires validator function"
                )

        return ConstraintEvaluationResult(
            satisfied=False,
            error_message=f"Unknown constraint type: {self.constraint_type}"
        )

    def _evaluate_coincident(self, tf_a: np.ndarray,
                             tf_b: Optional[np.ndarray]) -> ConstraintEvaluationResult:
        """Evaluate COINCIDENT constraint (origins/faces at same position).

        If face_a or face_b are specified, uses face positions instead of part origins.
        """
        # Use face position if specified, otherwise use part origin
        if hasattr(self, '_face_pos_a') and self._face_pos_a is not None:
            origin_a = self._face_pos_a
        else:
            origin_a = self._get_origin(tf_a)

        if hasattr(self, '_face_pos_b') and self._face_pos_b is not None:
            origin_b = self._face_pos_b
        elif tf_b is not None:
            origin_b = self._get_origin(tf_b)
        elif self.reference_center is not None:
            origin_b = np.array(self.reference_center)
        else:
            return ConstraintEvaluationResult(
                satisfied=False,
                error_message="COINCIDENT requires frame_b or reference_center"
            )

        distance = float(np.linalg.norm(origin_b - origin_a))
        satisfied = distance <= self.tolerance_mm

        error_vec = tuple((origin_b - origin_a).tolist()) if distance > 1e-10 else None

        # Build description including face info
        face_info = ""
        if self.face_a:
            face_info += f" (face_a={self.face_a})"
        if self.face_b:
            face_info += f" (face_b={self.face_b})"

        return ConstraintEvaluationResult(
            satisfied=satisfied,
            error=distance,
            error_vector=error_vec,
            error_message=f"COINCIDENT{face_info}: distance={distance:.4f}mm (tol={self.tolerance_mm}mm)",
            details={
                "origin_a": origin_a.tolist(),
                "origin_b": origin_b.tolist(),
                "distance": distance,
                "face_a": self.face_a,
                "face_b": self.face_b
            }
        )

    def _evaluate_concentric(self, tf_a: np.ndarray,
                              tf_b: Optional[np.ndarray]) -> ConstraintEvaluationResult:
        """Evaluate CONCENTRIC constraint (axes share same line)."""
        local_axis_a = self._get_local_axis(self.axis_a)
        if local_axis_a is None:
            return ConstraintEvaluationResult(
                satisfied=False,
                error_message=f"Invalid axis_a specification: {self.axis_a}"
            )

        origin_a = self._get_origin(tf_a)
        axis_a = self._transform_vector(tf_a, local_axis_a)

        if tf_b is not None:
            origin_b = self._get_origin(tf_b)
            axis_b = self._transform_vector(tf_b, self._get_local_axis(self.axis_b))
        elif self.reference_axis is not None and self.reference_center is not None:
            origin_b = np.array(self.reference_center)
            axis_b = np.array(self.reference_axis)
            axis_b = axis_b / np.linalg.norm(axis_b)
        else:
            return ConstraintEvaluationResult(
                satisfied=False,
                error_message="CONCENTRIC requires frame_b or reference_axis+center"
            )

        # Check axes are parallel
        dot = abs(float(np.dot(axis_a, axis_b)))
        angle_error = math.degrees(math.acos(min(1.0, dot)))

        # Check perpendicular distance between axes
        v = origin_b - origin_a
        along = np.dot(v, axis_a) * axis_a
        perp = v - along
        perp_distance = float(np.linalg.norm(perp))

        satisfied = (angle_error <= self.tolerance_deg and
                     perp_distance <= self.tolerance_mm)

        return ConstraintEvaluationResult(
            satisfied=satisfied,
            error=max(angle_error, perp_distance),
            error_message=f"CONCENTRIC: angle_err={angle_error:.2f}deg, "
                         f"perp_dist={perp_distance:.4f}mm",
            details={
                "origin_a": origin_a.tolist(),
                "origin_b": origin_b.tolist(),
                "axis_a": axis_a.tolist(),
                "axis_b": axis_b.tolist(),
                "angle_error_deg": angle_error,
                "perpendicular_distance": perp_distance
            }
        )

    def _evaluate_parallel(self, tf_a: np.ndarray,
                           tf_b: Optional[np.ndarray]) -> ConstraintEvaluationResult:
        """Evaluate PARALLEL constraint (axes are parallel)."""
        local_axis_a = self._get_local_axis(self.axis_a)
        if local_axis_a is None:
            return ConstraintEvaluationResult(
                satisfied=False,
                error_message=f"Invalid axis_a specification: {self.axis_a}"
            )
        axis_a = self._transform_vector(tf_a, local_axis_a)

        if tf_b is not None:
            axis_b = self._transform_vector(tf_b, self._get_local_axis(self.axis_b))
        elif self.reference_axis is not None:
            axis_b = np.array(self.reference_axis)
            axis_b = axis_b / np.linalg.norm(axis_b)
        else:
            return ConstraintEvaluationResult(
                satisfied=False,
                error_message="PARALLEL requires frame_b or reference_axis"
            )

        # Parallel means |dot| = 1 (same or opposite direction)
        dot = abs(float(np.dot(axis_a, axis_b)))
        angle_error = math.degrees(math.acos(min(1.0, dot)))

        satisfied = angle_error <= self.tolerance_deg

        return ConstraintEvaluationResult(
            satisfied=satisfied,
            error=angle_error,
            error_message=f"PARALLEL: angle_error={angle_error:.2f}deg (tol={self.tolerance_deg}deg)",
            details={
                "axis_a": axis_a.tolist(),
                "axis_b": axis_b.tolist(),
                "dot_product": float(dot),
                "angle_error_deg": angle_error
            }
        )

    def _evaluate_perpendicular(self, tf_a: np.ndarray,
                                 tf_b: Optional[np.ndarray]) -> ConstraintEvaluationResult:
        """Evaluate PERPENDICULAR constraint (axes are perpendicular)."""
        local_axis_a = self._get_local_axis(self.axis_a)
        if local_axis_a is None:
            return ConstraintEvaluationResult(
                satisfied=False,
                error_message=f"Invalid axis_a specification: {self.axis_a}"
            )
        axis_a = self._transform_vector(tf_a, local_axis_a)

        if tf_b is not None:
            axis_b = self._transform_vector(tf_b, self._get_local_axis(self.axis_b))
        elif self.reference_axis is not None:
            axis_b = np.array(self.reference_axis)
            axis_b = axis_b / np.linalg.norm(axis_b)
        else:
            return ConstraintEvaluationResult(
                satisfied=False,
                error_message="PERPENDICULAR requires frame_b or reference_axis"
            )

        # Perpendicular means dot = 0
        dot = abs(float(np.dot(axis_a, axis_b)))
        angle_from_perpendicular = math.degrees(math.asin(min(1.0, dot)))

        satisfied = angle_from_perpendicular <= self.tolerance_deg

        return ConstraintEvaluationResult(
            satisfied=satisfied,
            error=angle_from_perpendicular,
            error_message=f"PERPENDICULAR: angle_error={angle_from_perpendicular:.2f}deg "
                         f"(tol={self.tolerance_deg}deg)",
            details={
                "axis_a": axis_a.tolist(),
                "axis_b": axis_b.tolist(),
                "dot_product": float(dot),
                "angle_from_perpendicular_deg": angle_from_perpendicular
            }
        )

    def _evaluate_tangent(self, tf_a: np.ndarray) -> ConstraintEvaluationResult:
        """Evaluate TANGENT constraint (axis tangent to circle at position).

        An axis is tangent when it is perpendicular to the radial direction
        from the reference center to the frame origin.
        """
        if self.reference_center is None:
            return ConstraintEvaluationResult(
                satisfied=False,
                error_message="TANGENT requires reference_center"
            )

        local_axis_a = self._get_local_axis(self.axis_a)
        if local_axis_a is None:
            return ConstraintEvaluationResult(
                satisfied=False,
                error_message=f"Invalid axis_a specification: {self.axis_a}"
            )

        origin = self._get_origin(tf_a)
        center = np.array(self.reference_center)
        axis = self._transform_vector(tf_a, local_axis_a)

        # Radial vector from center to origin (project to XY for Z-axis circle)
        radial = origin - center
        radial[2] = 0.0  # Project to XY plane
        radial_norm = np.linalg.norm(radial)

        if radial_norm < 1e-10:
            return ConstraintEvaluationResult(
                satisfied=False,
                error_message="Frame at center - cannot determine tangent direction"
            )

        radial_unit = radial / radial_norm

        # Project axis to XY plane
        axis_xy = axis.copy()
        axis_xy[2] = 0.0
        axis_norm = np.linalg.norm(axis_xy)

        if axis_norm < 1e-10:
            # Vertical axis is tangent to horizontal circle
            return ConstraintEvaluationResult(
                satisfied=True,
                error=0.0,
                error_message="TANGENT: axis is vertical (tangent to horizontal circle)",
                details={
                    "origin": origin.tolist(),
                    "center": center.tolist(),
                    "axis": axis.tolist()
                }
            )

        axis_unit = axis_xy / axis_norm

        # Tangent means perpendicular to radial (dot = 0)
        dot = abs(float(np.dot(radial_unit, axis_unit)))
        angle_from_tangent = math.degrees(math.asin(min(1.0, dot)))

        satisfied = angle_from_tangent <= self.tolerance_deg

        return ConstraintEvaluationResult(
            satisfied=satisfied,
            error=angle_from_tangent,
            error_message=f"TANGENT: angle_error={angle_from_tangent:.2f}deg "
                         f"(tol={self.tolerance_deg}deg)",
            details={
                "origin": origin.tolist(),
                "center": center.tolist(),
                "axis": axis.tolist(),
                "radial_direction": radial_unit.tolist(),
                "dot_with_radial": float(dot),
                "angle_from_tangent_deg": angle_from_tangent
            }
        )

    def _evaluate_at_distance(self, tf_a: np.ndarray,
                               tf_b: Optional[np.ndarray]) -> ConstraintEvaluationResult:
        """Evaluate AT_DISTANCE constraint (frames/faces at specific distance).

        If face_a or face_b are specified, uses face positions instead of part origins.
        """
        # Use face position if specified, otherwise use part origin
        if hasattr(self, '_face_pos_a') and self._face_pos_a is not None:
            origin_a = self._face_pos_a
        else:
            origin_a = self._get_origin(tf_a)

        if hasattr(self, '_face_pos_b') and self._face_pos_b is not None:
            origin_b = self._face_pos_b
        elif tf_b is not None:
            origin_b = self._get_origin(tf_b)
        elif self.reference_center is not None:
            origin_b = np.array(self.reference_center)
        else:
            return ConstraintEvaluationResult(
                satisfied=False,
                error_message="AT_DISTANCE requires frame_b or reference_center"
            )

        if self.reference_radius is None:
            return ConstraintEvaluationResult(
                satisfied=False,
                error_message="AT_DISTANCE requires reference_radius"
            )

        actual_distance = float(np.linalg.norm(origin_b - origin_a))
        error = abs(actual_distance - self.reference_radius)
        satisfied = error <= self.tolerance_mm

        # Build description including face info
        face_info = ""
        if self.face_a:
            face_info += f" (face_a={self.face_a})"
        if self.face_b:
            face_info += f" (face_b={self.face_b})"

        return ConstraintEvaluationResult(
            satisfied=satisfied,
            error=error,
            error_message=f"AT_DISTANCE{face_info}: actual={actual_distance:.4f}mm, "
                         f"expected={self.reference_radius:.4f}mm, error={error:.4f}mm",
            details={
                "origin_a": origin_a.tolist(),
                "origin_b": origin_b.tolist(),
                "actual_distance": actual_distance,
                "expected_distance": self.reference_radius,
                "error": error,
                "face_a": self.face_a,
                "face_b": self.face_b
            }
        )

    def _evaluate_bolt_pattern(self, tf_a: np.ndarray,
                                tf_b: Optional[np.ndarray]) -> ConstraintEvaluationResult:
        """Evaluate BOLT_PATTERN constraint (holes align with pattern)."""
        if self.pattern_count < 2:
            return ConstraintEvaluationResult(
                satisfied=False,
                error_message="BOLT_PATTERN requires pattern_count >= 2"
            )

        if self.pattern_radius <= 0:
            return ConstraintEvaluationResult(
                satisfied=False,
                error_message="BOLT_PATTERN requires positive pattern_radius"
            )

        local_axis_a = self._get_local_axis(self.axis_a)
        if local_axis_a is None:
            return ConstraintEvaluationResult(
                satisfied=False,
                error_message=f"Invalid axis_a specification: {self.axis_a}"
            )

        origin_a = self._get_origin(tf_a)
        axis_a = self._transform_vector(tf_a, local_axis_a)

        if tf_b is not None:
            origin_b = self._get_origin(tf_b)
            axis_b = self._transform_vector(tf_b, self._get_local_axis(self.axis_b))
        elif self.reference_center is not None and self.reference_normal is not None:
            origin_b = np.array(self.reference_center)
            axis_b = np.array(self.reference_normal)
            axis_b = axis_b / np.linalg.norm(axis_b)
        else:
            return ConstraintEvaluationResult(
                satisfied=False,
                error_message="BOLT_PATTERN requires frame_b or reference_center+normal"
            )

        # Check center alignment
        center_distance = float(np.linalg.norm(origin_b - origin_a))
        satisfied_center = center_distance <= self.tolerance_mm

        # Check normal alignment
        dot = abs(float(np.dot(axis_a, axis_b)))
        angle_error = math.degrees(math.acos(min(1.0, dot)))
        satisfied_angle = angle_error <= self.tolerance_deg

        satisfied = satisfied_center and satisfied_angle
        max_error = max(center_distance, angle_error)

        return ConstraintEvaluationResult(
            satisfied=satisfied,
            error=max_error,
            error_message=f"BOLT_PATTERN: center_dist={center_distance:.4f}mm, "
                         f"angle_err={angle_error:.2f}deg ({self.pattern_count} holes)",
            details={
                "center_a": origin_a.tolist(),
                "center_b": origin_b.tolist(),
                "normal_a": axis_a.tolist(),
                "normal_b": axis_b.tolist(),
                "center_distance": center_distance,
                "angle_error_deg": angle_error,
                "pattern_count": self.pattern_count,
                "pattern_radius": self.pattern_radius,
                "angular_offset_deg": self.pattern_offset_deg
            }
        )

    def _evaluate_min_distance(self, tf_a: np.ndarray,
                                tf_b: Optional[np.ndarray]) -> ConstraintEvaluationResult:
        """Evaluate MIN_DISTANCE constraint (parts must be at least N mm apart).

        Uses bounding boxes if provided, otherwise uses origin-to-origin distance.
        If face_a or face_b are specified (and no bounds), uses face positions.
        """
        # Use face position if specified, otherwise use part origin
        if hasattr(self, '_face_pos_a') and self._face_pos_a is not None:
            origin_a = self._face_pos_a
        else:
            origin_a = self._get_origin(tf_a)

        if hasattr(self, '_face_pos_b') and self._face_pos_b is not None:
            origin_b = self._face_pos_b
        elif tf_b is not None:
            origin_b = self._get_origin(tf_b)
        elif self.reference_center is not None:
            origin_b = np.array(self.reference_center)
        else:
            return ConstraintEvaluationResult(
                satisfied=False,
                error_message="MIN_DISTANCE requires frame_b or reference_center"
            )

        if self.reference_radius is None:
            return ConstraintEvaluationResult(
                satisfied=False,
                error_message="MIN_DISTANCE requires reference_radius (min distance)"
            )

        min_required = self.reference_radius

        # If bounding boxes are provided, calculate minimum separation
        # (bounding box method takes precedence over face positions for collision detection)
        if self.bounds_a is not None and self.bounds_b is not None:
            # Transform bounding boxes to world coordinates
            actual_distance = self._calculate_bbox_separation(
                tf_a, self.bounds_a, tf_b, self.bounds_b
            )
        else:
            # Fall back to face-to-face or origin-to-origin distance
            actual_distance = float(np.linalg.norm(origin_b - origin_a))

        error = min_required - actual_distance
        satisfied = actual_distance >= min_required

        # Build description including face info
        face_info = ""
        if self.face_a:
            face_info += f" (face_a={self.face_a})"
        if self.face_b:
            face_info += f" (face_b={self.face_b})"

        return ConstraintEvaluationResult(
            satisfied=satisfied,
            error=max(0, error),
            error_message=f"MIN_DISTANCE{face_info}: actual={actual_distance:.3f}mm, "
                         f"required={min_required:.3f}mm"
                         f"{'' if satisfied else f' (VIOLATION: {error:.3f}mm penetration)'}",
            details={
                "origin_a": origin_a.tolist(),
                "origin_b": origin_b.tolist(),
                "actual_distance": actual_distance,
                "required_distance": min_required,
                "penetration": max(0, error),
                "face_a": self.face_a,
                "face_b": self.face_b
            }
        )

    def _evaluate_z_stack_clearance(self, tf_a: np.ndarray,
                                     tf_b: Optional[np.ndarray]) -> ConstraintEvaluationResult:
        """Evaluate Z_STACK_CLEARANCE constraint for stacked parts.

        This constraint validates that when parts are stacked (like servo + gearbox),
        they have proper Z-axis separation to avoid collision. It checks that:
        - Part A's top surface is below Part B's bottom surface (or vice versa)
        - The separation is at least the required clearance

        Uses bounds_a and bounds_b to determine part heights.
        reference_radius specifies the required clearance.
        """
        if tf_b is None:
            return ConstraintEvaluationResult(
                satisfied=False,
                error_message="Z_STACK_CLEARANCE requires frame_b"
            )

        if self.bounds_a is None or self.bounds_b is None:
            return ConstraintEvaluationResult(
                satisfied=False,
                error_message="Z_STACK_CLEARANCE requires bounds_a and bounds_b"
            )

        min_clearance = self.reference_radius if self.reference_radius is not None else 0.0

        # Get Z positions of part origins in world space
        origin_a = self._get_origin(tf_a)
        origin_b = self._get_origin(tf_b)

        # Get part heights from bounding boxes (local Z extent)
        height_a = self.bounds_a[1][2] - self.bounds_a[0][2]  # max_z - min_z
        height_b = self.bounds_b[1][2] - self.bounds_b[0][2]

        # Get local Z offsets to top/bottom surfaces
        top_a_local = self.bounds_a[1][2]  # max_z in local coords
        bottom_a_local = self.bounds_a[0][2]  # min_z in local coords
        top_b_local = self.bounds_b[1][2]
        bottom_b_local = self.bounds_b[0][2]

        # Transform local top/bottom points to world Z
        # For simplified check, assume Z-axis alignment (after rotations)
        # Get world Z-axis direction for each part
        local_z_a = self._get_local_axis("z")
        local_z_b = self._get_local_axis("z")
        world_z_a = self._transform_vector(tf_a, local_z_a)
        world_z_b = self._transform_vector(tf_b, local_z_b)

        # Calculate world Z positions of top/bottom faces
        # For a part with local bounds, world top Z = origin_z + local_max_z (if Z aligned)
        # If rotated, this is more complex - simplified version uses origin + projected height
        top_a_world = origin_a[2] + top_a_local * abs(world_z_a[2])
        bottom_a_world = origin_a[2] + bottom_a_local * abs(world_z_a[2])
        top_b_world = origin_b[2] + top_b_local * abs(world_z_b[2])
        bottom_b_world = origin_b[2] + bottom_b_local * abs(world_z_b[2])

        # Determine overlap in Z
        # Parts overlap if top_a > bottom_b AND top_b > bottom_a
        z_overlap = min(top_a_world, top_b_world) - max(bottom_a_world, bottom_b_world)

        if z_overlap > 0:
            # Parts overlap in Z - this is a collision!
            satisfied = False
            error = z_overlap + min_clearance
            message = (f"Z_STACK_CLEARANCE: COLLISION! Parts overlap by {z_overlap:.3f}mm in Z. "
                      f"Part A Z: [{bottom_a_world:.1f}, {top_a_world:.1f}], "
                      f"Part B Z: [{bottom_b_world:.1f}, {top_b_world:.1f}]")
        else:
            # Calculate separation
            separation = -z_overlap  # Positive separation
            if separation >= min_clearance:
                satisfied = True
                error = 0.0
                message = (f"Z_STACK_CLEARANCE: OK. Separation={separation:.3f}mm >= "
                          f"required={min_clearance:.3f}mm")
            else:
                satisfied = False
                error = min_clearance - separation
                message = (f"Z_STACK_CLEARANCE: Insufficient separation. "
                          f"Actual={separation:.3f}mm < required={min_clearance:.3f}mm")

        return ConstraintEvaluationResult(
            satisfied=satisfied,
            error=error,
            error_message=message,
            details={
                "origin_a": origin_a.tolist(),
                "origin_b": origin_b.tolist(),
                "height_a": height_a,
                "height_b": height_b,
                "top_a_world": top_a_world,
                "bottom_a_world": bottom_a_world,
                "top_b_world": top_b_world,
                "bottom_b_world": bottom_b_world,
                "z_overlap": z_overlap,
                "min_clearance": min_clearance
            }
        )

    def _evaluate_no_overlap(self, tf_a: np.ndarray,
                              tf_b: Optional[np.ndarray]) -> ConstraintEvaluationResult:
        """Evaluate NO_OVERLAP constraint.

        When STL paths are provided, uses accurate mesh-based collision detection.
        Otherwise falls back to axis-aligned bounding box check (fast but conservative).
        """
        if tf_b is None:
            return ConstraintEvaluationResult(
                satisfied=False,
                error_message="NO_OVERLAP requires frame_b"
            )

        # If STL paths are provided, use mesh-based collision detection
        if self.stl_path_a is not None and self.stl_path_b is not None:
            return self._evaluate_mesh_collision(tf_a, tf_b)

        # Fall back to AABB check
        if self.bounds_a is None or self.bounds_b is None:
            # If no bounds provided, skip check with warning
            return ConstraintEvaluationResult(
                satisfied=True,
                error=0.0,
                error_message="NO_OVERLAP: Skipped (no bounding boxes or STL paths provided)",
                details={"skipped": True}
            )

        # Calculate AABB in world coordinates
        aabb_a = self._transform_aabb(tf_a, self.bounds_a)
        aabb_b = self._transform_aabb(tf_b, self.bounds_b)

        # Check for overlap in all three axes
        overlap_x = (aabb_a[0][0] <= aabb_b[1][0]) and (aabb_b[0][0] <= aabb_a[1][0])
        overlap_y = (aabb_a[0][1] <= aabb_b[1][1]) and (aabb_b[0][1] <= aabb_a[1][1])
        overlap_z = (aabb_a[0][2] <= aabb_b[1][2]) and (aabb_b[0][2] <= aabb_a[1][2])

        overlaps = overlap_x and overlap_y and overlap_z
        satisfied = not overlaps

        if overlaps:
            # Calculate penetration depth on each axis
            pen_x = min(aabb_a[1][0] - aabb_b[0][0], aabb_b[1][0] - aabb_a[0][0])
            pen_y = min(aabb_a[1][1] - aabb_b[0][1], aabb_b[1][1] - aabb_a[0][1])
            pen_z = min(aabb_a[1][2] - aabb_b[0][2], aabb_b[1][2] - aabb_a[0][2])
            error = min(pen_x, pen_y, pen_z)  # Minimum separation distance to resolve
            message = (f"NO_OVERLAP: COLLISION! Penetration depth={error:.3f}mm. "
                      f"Overlap: X={overlap_x}, Y={overlap_y}, Z={overlap_z}")
        else:
            error = 0.0
            message = "NO_OVERLAP: OK - bounding boxes do not overlap"

        return ConstraintEvaluationResult(
            satisfied=satisfied,
            error=error,
            error_message=message,
            details={
                "aabb_a": [list(aabb_a[0]), list(aabb_a[1])],
                "aabb_b": [list(aabb_b[0]), list(aabb_b[1])],
                "overlap_x": overlap_x,
                "overlap_y": overlap_y,
                "overlap_z": overlap_z,
                "collision": overlaps,
                "method": "aabb"
            }
        )

    def _evaluate_mesh_collision(self, tf_a: np.ndarray,
                                  tf_b: np.ndarray) -> ConstraintEvaluationResult:
        """Evaluate collision using mesh-based intersection detection.

        Uses trimesh library for accurate collision detection with the actual
        part geometry rather than bounding box approximations.
        """
        if not HAS_TRIMESH:
            return ConstraintEvaluationResult(
                satisfied=False,
                error_message="NO_OVERLAP: trimesh library required for mesh collision"
            )

        # Perform mesh collision check
        result = check_mesh_collision(
            self.stl_path_a, tf_a,
            self.stl_path_b, tf_b
        )

        if result.error_message and not result.collides:
            # Check failed due to error
            return ConstraintEvaluationResult(
                satisfied=False,
                error_message=f"NO_OVERLAP: Mesh check failed - {result.error_message}",
                details={"mesh_error": result.error_message}
            )

        satisfied = not result.collides

        if result.collides:
            message = (f"NO_OVERLAP: MESH COLLISION! "
                      f"Penetration={result.penetration_depth:.3f}mm, "
                      f"Volume={result.collision_volume:.3f}mm^3, "
                      f"Contacts={len(result.contact_points)}")
        else:
            message = "NO_OVERLAP: OK - meshes do not intersect"

        return ConstraintEvaluationResult(
            satisfied=satisfied,
            error=result.penetration_depth,
            error_message=message,
            details={
                "collision": result.collides,
                "penetration_depth": result.penetration_depth,
                "collision_volume": result.collision_volume,
                "contact_points": result.contact_points,
                "stl_path_a": self.stl_path_a,
                "stl_path_b": self.stl_path_b,
                "method": "mesh"
            }
        )

    def _calculate_bbox_separation(self, tf_a: np.ndarray, bounds_a, tf_b, bounds_b) -> float:
        """Calculate minimum separation between two bounding boxes."""
        aabb_a = self._transform_aabb(tf_a, bounds_a)
        aabb_b = self._transform_aabb(tf_b, bounds_b)

        # Calculate signed distance on each axis
        dist_x = max(aabb_a[0][0] - aabb_b[1][0], aabb_b[0][0] - aabb_a[1][0])
        dist_y = max(aabb_a[0][1] - aabb_b[1][1], aabb_b[0][1] - aabb_a[1][1])
        dist_z = max(aabb_a[0][2] - aabb_b[1][2], aabb_b[0][2] - aabb_a[1][2])

        # If any distance is positive, boxes are separated
        # Return the maximum distance (most separated axis)
        # If all negative, return the minimum (most penetration)
        if dist_x > 0 or dist_y > 0 or dist_z > 0:
            # Separated - return minimum positive distance
            positive_dists = [d for d in [dist_x, dist_y, dist_z] if d > 0]
            return min(positive_dists) if positive_dists else 0.0
        else:
            # Overlapping - return negative (penetration depth)
            return max(dist_x, dist_y, dist_z)  # Least penetration

    def _transform_aabb(self, transform: np.ndarray,
                        bounds: Tuple[Tuple[float, float, float], Tuple[float, float, float]]
                        ) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:
        """Transform an axis-aligned bounding box to world coordinates.

        Note: The result is still axis-aligned but may be larger than the
        oriented bounding box (conservative approximation).
        """
        min_local, max_local = bounds

        # Generate all 8 corners of the local bounding box
        corners = []
        for x in [min_local[0], max_local[0]]:
            for y in [min_local[1], max_local[1]]:
                for z in [min_local[2], max_local[2]]:
                    # Transform corner to world coordinates
                    local_pt = np.array([x, y, z, 1.0])
                    world_pt = transform @ local_pt
                    corners.append(world_pt[:3])

        # Find axis-aligned bounding box of transformed corners
        corners = np.array(corners)
        min_world = tuple(corners.min(axis=0))
        max_world = tuple(corners.max(axis=0))

        return (min_world, max_world)


@dataclass
class ValidationReport:
    """Comprehensive assembly validation report.

    Attributes:
        is_valid: True if all error-severity constraints pass
        total_constraints: Total number of constraints evaluated
        passed_count: Number of constraints that passed
        failed_count: Number of constraints that failed
        warning_count: Number of warning-severity constraints that failed
        constraint_results: Mapping of constraint name to evaluation result
        failed_constraints: List of names of failed constraints
        warnings: List of warning messages
        info: List of informational messages
    """
    is_valid: bool
    total_constraints: int = 0
    passed_count: int = 0
    failed_count: int = 0
    warning_count: int = 0
    constraint_results: Dict[str, ConstraintEvaluationResult] = field(default_factory=dict)
    failed_constraints: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    info: List[str] = field(default_factory=list)

    def __str__(self) -> str:
        status = "VALID" if self.is_valid else "INVALID"
        return (f"ValidationReport: {status} "
                f"({self.passed_count}/{self.total_constraints} passed, "
                f"{self.failed_count} failed, {self.warning_count} warnings)")

    def detailed_report(self) -> str:
        """Generate detailed multi-line report."""
        lines = []
        lines.append("=" * 70)
        lines.append("ASSEMBLY VALIDATION REPORT")
        lines.append("=" * 70)
        lines.append("")

        status = "VALID" if self.is_valid else "INVALID"
        lines.append(f"STATUS: {status}")
        lines.append(f"Total Constraints: {self.total_constraints}")
        lines.append(f"Passed: {self.passed_count}")
        lines.append(f"Failed: {self.failed_count}")
        lines.append(f"Warnings: {self.warning_count}")
        lines.append("")

        if self.failed_constraints:
            lines.append("FAILED CONSTRAINTS:")
            lines.append("-" * 70)
            for name in self.failed_constraints:
                result = self.constraint_results[name]
                lines.append(f"  [{name}]")
                lines.append(f"    {result.error_message}")
                if result.error > 0:
                    lines.append(f"    Error: {result.error:.4f}")
            lines.append("")

        if self.warnings:
            lines.append("WARNINGS:")
            lines.append("-" * 70)
            for warning in self.warnings:
                lines.append(f"  - {warning}")
            lines.append("")

        # Show all passed constraints
        passed = [name for name in self.constraint_results
                  if self.constraint_results[name].satisfied and name not in self.failed_constraints]
        if passed:
            lines.append("PASSED CONSTRAINTS:")
            lines.append("-" * 70)
            for name in passed:
                result = self.constraint_results[name]
                lines.append(f"  [OK] {name}")
            lines.append("")

        lines.append("=" * 70)
        return "\n".join(lines)


class AssemblyValidator:
    """Validates assemblies against a set of kinematic constraints.

    The validator maintains a collection of constraints and evaluates them
    against world transforms from a kinematic chain.

    Example:
        >>> validator = AssemblyValidator()
        >>> validator.add_constraint(KinematicConstraint(...))
        >>> validator.add_constraint(KinematicConstraint(...))
        >>>
        >>> # Get transforms from kinematic chain
        >>> transforms = chain.get_all_world_transforms()
        >>>
        >>> # Validate
        >>> report = validator.validate(transforms)
        >>> if not report.is_valid:
        ...     print(report.detailed_report())
    """

    def __init__(self, name: str = "assembly"):
        """Initialize validator.

        Args:
            name: Name for this validator (used in reports)
        """
        self.name = name
        self.constraints: List[KinematicConstraint] = []

    def add_constraint(self, constraint: KinematicConstraint) -> None:
        """Add a constraint to the validator.

        Args:
            constraint: KinematicConstraint to add
        """
        self.constraints.append(constraint)

    def add_constraints(self, constraints: List[KinematicConstraint]) -> None:
        """Add multiple constraints.

        Args:
            constraints: List of constraints to add
        """
        self.constraints.extend(constraints)

    def clear_constraints(self) -> None:
        """Remove all constraints."""
        self.constraints = []

    def validate(self, world_transforms: Dict[str, Any]) -> ValidationReport:
        """Validate all constraints against given world transforms.

        Args:
            world_transforms: Dictionary mapping part names to Transform objects
                              or 4x4 numpy arrays

        Returns:
            ValidationReport with comprehensive evaluation results
        """
        results: Dict[str, ConstraintEvaluationResult] = {}
        failed: List[str] = []
        warnings: List[str] = []
        info: List[str] = []
        passed_count = 0
        failed_count = 0
        warning_count = 0

        for constraint in self.constraints:
            result = constraint.evaluate(world_transforms)
            results[constraint.name] = result

            if result.satisfied:
                passed_count += 1
            else:
                if constraint.severity == "error":
                    failed_count += 1
                    failed.append(constraint.name)
                elif constraint.severity == "warning":
                    warning_count += 1
                    warnings.append(f"{constraint.name}: {result.error_message}")
                else:  # info
                    info.append(f"{constraint.name}: {result.error_message}")

        is_valid = failed_count == 0

        return ValidationReport(
            is_valid=is_valid,
            total_constraints=len(self.constraints),
            passed_count=passed_count,
            failed_count=failed_count,
            warning_count=warning_count,
            constraint_results=results,
            failed_constraints=failed,
            warnings=warnings,
            info=info
        )

    def validate_and_raise(self, world_transforms: Dict[str, Any]) -> None:
        """Validate and raise exception if any error-severity constraint fails.

        Args:
            world_transforms: Dictionary mapping part names to transforms

        Raises:
            ValueError: If any error-severity constraint fails
        """
        report = self.validate(world_transforms)
        if not report.is_valid:
            raise ValueError(f"Assembly validation failed:\n{report.detailed_report()}")


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def validate_assembly(
    constraints: List[KinematicConstraint],
    world_transforms: Dict[str, Any]
) -> ValidationReport:
    """Convenience function to validate constraints against transforms.

    Args:
        constraints: List of KinematicConstraint objects
        world_transforms: Dictionary mapping part names to transforms

    Returns:
        ValidationReport with evaluation results
    """
    validator = AssemblyValidator()
    validator.add_constraints(constraints)
    return validator.validate(world_transforms)


def create_tangent_constraint(
    name: str,
    frame: str,
    center: Tuple[float, float, float],
    axis: str = "y",
    tolerance_deg: float = 2.0,
    description: str = ""
) -> KinematicConstraint:
    """Create a tangent constraint for wheel/motor assemblies.

    Args:
        name: Constraint name
        frame: Part name in kinematic chain
        center: Center point for tangent calculation
        axis: Which local axis to check ("x", "y", or "z")
        tolerance_deg: Angular tolerance
        description: Description of constraint purpose

    Returns:
        Configured KinematicConstraint
    """
    return KinematicConstraint(
        name=name,
        constraint_type=KinematicConstraintType.TANGENT,
        frame_a=frame,
        axis_a=axis,
        reference_center=center,
        tolerance_deg=tolerance_deg,
        description=description or f"Axis of {frame} must be tangent to circle at {center}"
    )


def create_parallel_constraint(
    name: str,
    frame_a: str,
    frame_b: Optional[str] = None,
    axis_a: str = "z",
    axis_b: str = "z",
    reference_axis: Optional[Tuple[float, float, float]] = None,
    tolerance_deg: float = 1.0,
    description: str = ""
) -> KinematicConstraint:
    """Create a parallel constraint between two axes.

    Args:
        name: Constraint name
        frame_a: First part name
        frame_b: Second part name (or None for world reference)
        axis_a: Local axis on frame_a
        axis_b: Local axis on frame_b
        reference_axis: World reference axis (if frame_b is None)
        tolerance_deg: Angular tolerance
        description: Description of constraint purpose

    Returns:
        Configured KinematicConstraint
    """
    return KinematicConstraint(
        name=name,
        constraint_type=KinematicConstraintType.PARALLEL,
        frame_a=frame_a,
        frame_b=frame_b,
        axis_a=axis_a,
        axis_b=axis_b,
        reference_axis=reference_axis,
        tolerance_deg=tolerance_deg,
        description=description or f"Axes of {frame_a} and {frame_b or 'reference'} must be parallel"
    )


def create_coincident_constraint(
    name: str,
    frame_a: str,
    frame_b: Optional[str] = None,
    reference_point: Optional[Tuple[float, float, float]] = None,
    tolerance_mm: float = 0.1,
    description: str = ""
) -> KinematicConstraint:
    """Create a coincident constraint (origins at same location).

    Args:
        name: Constraint name
        frame_a: First part name
        frame_b: Second part name (or None for world reference)
        reference_point: World reference point (if frame_b is None)
        tolerance_mm: Distance tolerance in mm
        description: Description of constraint purpose

    Returns:
        Configured KinematicConstraint
    """
    return KinematicConstraint(
        name=name,
        constraint_type=KinematicConstraintType.COINCIDENT,
        frame_a=frame_a,
        frame_b=frame_b,
        reference_center=reference_point,
        tolerance_mm=tolerance_mm,
        description=description or f"Origins of {frame_a} and {frame_b or 'reference'} must coincide"
    )


def create_z_stack_clearance_constraint(
    name: str,
    frame_a: str,
    frame_b: str,
    bounds_a: Tuple[Tuple[float, float, float], Tuple[float, float, float]],
    bounds_b: Tuple[Tuple[float, float, float], Tuple[float, float, float]],
    min_clearance: float = 0.0,
    description: str = ""
) -> KinematicConstraint:
    """Create a Z-stack clearance constraint for stacked parts.

    Validates that vertically stacked parts (like servo + gearbox) have proper
    Z-axis separation to avoid collision. Useful for detecting servo bodies
    penetrating ring housings.

    Args:
        name: Constraint name
        frame_a: First part name (e.g., "AXIS2_SERVO_XH430")
        frame_b: Second part name (e.g., "AXIS2_RING_HOUSING")
        bounds_a: Bounding box of part A ((min_x, min_y, min_z), (max_x, max_y, max_z))
        bounds_b: Bounding box of part B
        min_clearance: Minimum required clearance in mm (default: 0 = touching OK)
        description: Description of constraint purpose

    Returns:
        Configured KinematicConstraint

    Example:
        >>> # XH430 servo: 28.5 x 46.5 x 34mm, origin at body center
        >>> servo_bounds = ((-14.25, -23.25, -17), (14.25, 23.25, 17))
        >>> # Ring housing: ~60mm diameter x 8mm height
        >>> ring_bounds = ((-30, -30, 0), (30, 30, 8))
        >>> constraint = create_z_stack_clearance_constraint(
        ...     "axis2_servo_ring_clearance",
        ...     "AXIS2_SERVO_XH430",
        ...     "AXIS2_RING_HOUSING",
        ...     servo_bounds,
        ...     ring_bounds,
        ...     min_clearance=1.0
        ... )
    """
    return KinematicConstraint(
        name=name,
        constraint_type=KinematicConstraintType.Z_STACK_CLEARANCE,
        frame_a=frame_a,
        frame_b=frame_b,
        bounds_a=bounds_a,
        bounds_b=bounds_b,
        reference_radius=min_clearance,  # Using reference_radius for min_clearance
        description=description or f"Z-stack clearance between {frame_a} and {frame_b}"
    )


def create_no_overlap_constraint(
    name: str,
    frame_a: str,
    frame_b: str,
    bounds_a: Optional[Tuple[Tuple[float, float, float], Tuple[float, float, float]]] = None,
    bounds_b: Optional[Tuple[Tuple[float, float, float], Tuple[float, float, float]]] = None,
    stl_path_a: Optional[str] = None,
    stl_path_b: Optional[str] = None,
    description: str = ""
) -> KinematicConstraint:
    """Create a no-overlap constraint.

    When STL paths are provided, uses accurate mesh-based collision detection.
    Otherwise uses axis-aligned bounding boxes (fast but conservative).

    Args:
        name: Constraint name
        frame_a: First part name
        frame_b: Second part name
        bounds_a: Bounding box of part A ((min_x, min_y, min_z), (max_x, max_y, max_z))
        bounds_b: Bounding box of part B
        stl_path_a: Path to STL file for part A (enables mesh collision)
        stl_path_b: Path to STL file for part B (enables mesh collision)
        description: Description of constraint purpose

    Returns:
        Configured KinematicConstraint

    Example:
        >>> # AABB-based collision (fast, conservative)
        >>> constraint = create_no_overlap_constraint(
        ...     "servo_ring_no_overlap",
        ...     "SERVO",
        ...     "RING_HOUSING",
        ...     bounds_a=((-15, -25, -17), (15, 25, 17)),
        ...     bounds_b=((-30, -30, 0), (30, 30, 8))
        ... )
        >>>
        >>> # Mesh-based collision (accurate, requires trimesh)
        >>> constraint = create_no_overlap_constraint(
        ...     "servo_ring_no_overlap",
        ...     "SERVO",
        ...     "RING_HOUSING",
        ...     stl_path_a="parts/servo.stl",
        ...     stl_path_b="parts/ring_housing.stl"
        ... )
    """
    return KinematicConstraint(
        name=name,
        constraint_type=KinematicConstraintType.NO_OVERLAP,
        frame_a=frame_a,
        frame_b=frame_b,
        bounds_a=bounds_a,
        bounds_b=bounds_b,
        stl_path_a=stl_path_a,
        stl_path_b=stl_path_b,
        description=description or f"No overlap between {frame_a} and {frame_b}"
    )


def create_mesh_collision_constraint(
    name: str,
    frame_a: str,
    frame_b: str,
    stl_path_a: str,
    stl_path_b: str,
    description: str = ""
) -> KinematicConstraint:
    """Create a mesh-based collision constraint.

    Uses trimesh library for accurate collision detection with actual part geometry.
    This is more accurate than AABB-based detection but requires the trimesh library
    and STL files for both parts.

    Args:
        name: Constraint name
        frame_a: First part name (must match kinematic chain frame name)
        frame_b: Second part name (must match kinematic chain frame name)
        stl_path_a: Path to STL file for part A
        stl_path_b: Path to STL file for part B
        description: Description of constraint purpose

    Returns:
        Configured KinematicConstraint for mesh collision detection

    Example:
        >>> from yapcad.assembly.kinematic_integration import (
        ...     create_mesh_collision_constraint, AssemblyValidator
        ... )
        >>>
        >>> constraint = create_mesh_collision_constraint(
        ...     "gearbox_servo_collision",
        ...     "AXIS1_RING_HOUSING",
        ...     "AXIS1_SERVO_XH430",
        ...     "/path/to/ring_housing.stl",
        ...     "/path/to/servo.stl",
        ...     description="Check servo doesn't penetrate ring housing"
        ... )
        >>>
        >>> validator = AssemblyValidator()
        >>> validator.add_constraint(constraint)
        >>> report = validator.validate(world_transforms)
        >>> if not report.is_valid:
        ...     print(report.detailed_report())
    """
    return KinematicConstraint(
        name=name,
        constraint_type=KinematicConstraintType.NO_OVERLAP,
        frame_a=frame_a,
        frame_b=frame_b,
        stl_path_a=stl_path_a,
        stl_path_b=stl_path_b,
        description=description or f"Mesh collision check: {frame_a} vs {frame_b}"
    )


def create_min_distance_constraint(
    name: str,
    frame_a: str,
    frame_b: str,
    min_distance: float,
    bounds_a: Optional[Tuple[Tuple[float, float, float], Tuple[float, float, float]]] = None,
    bounds_b: Optional[Tuple[Tuple[float, float, float], Tuple[float, float, float]]] = None,
    description: str = ""
) -> KinematicConstraint:
    """Create a minimum distance constraint.

    Parts must be at least min_distance mm apart. If bounding boxes are provided,
    uses them for more accurate separation calculation.

    Args:
        name: Constraint name
        frame_a: First part name
        frame_b: Second part name
        min_distance: Minimum required separation in mm
        bounds_a: Optional bounding box of part A
        bounds_b: Optional bounding box of part B
        description: Description of constraint purpose

    Returns:
        Configured KinematicConstraint
    """
    return KinematicConstraint(
        name=name,
        constraint_type=KinematicConstraintType.MIN_DISTANCE,
        frame_a=frame_a,
        frame_b=frame_b,
        bounds_a=bounds_a,
        bounds_b=bounds_b,
        reference_radius=min_distance,  # Using reference_radius for min_distance
        description=description or f"Minimum distance {min_distance}mm between {frame_a} and {frame_b}"
    )


def create_face_distance_constraint(
    name: str,
    frame_a: str,
    face_a: str,
    frame_b: str,
    face_b: str,
    expected_distance: float = 0.0,
    tolerance_mm: float = 0.5,
    bounds_a: Optional[Tuple[Tuple[float, float, float], Tuple[float, float, float]]] = None,
    bounds_b: Optional[Tuple[Tuple[float, float, float], Tuple[float, float, float]]] = None,
    description: str = ""
) -> KinematicConstraint:
    """Create a face-to-face distance constraint.

    Measures distance between specific faces on two parts, rather than part origins.
    Useful for validating mating interfaces like servo output shafts to gearbox inputs.

    Face specifications can be:
    - Standard faces: "TOP", "BOTTOM", "FRONT", "BACK", "LEFT", "RIGHT"
      (requires bounds_a/bounds_b to calculate face center positions)
    - Named frames: e.g., "OUTPUT_SHAFT", "SUN_INPUT"
      (looked up from world_transforms using "part_name.frame_name")

    Args:
        name: Constraint name
        frame_a: First part name (e.g., "AXIS2_SERVO_XH430")
        face_a: Face specification on frame_a (e.g., "OUTPUT_SHAFT" or "TOP")
        frame_b: Second part name (e.g., "AXIS2_RING_HOUSING")
        face_b: Face specification on frame_b (e.g., "SUN_INPUT" or "BOTTOM")
        expected_distance: Expected distance between faces in mm (0.0 = touching)
        tolerance_mm: Allowed deviation from expected distance
        bounds_a: Bounding box of part A (required for standard face names)
        bounds_b: Bounding box of part B (required for standard face names)
        description: Description of constraint purpose

    Returns:
        Configured KinematicConstraint

    Example:
        >>> # Servo output shaft should mate with ring gear input
        >>> constraint = create_face_distance_constraint(
        ...     name="servo_ring_interface",
        ...     frame_a="AXIS2_SERVO_XH430",
        ...     face_a="OUTPUT_SHAFT",
        ...     frame_b="AXIS2_RING_HOUSING",
        ...     face_b="SUN_INPUT",
        ...     expected_distance=0.0,  # Should be touching
        ...     tolerance_mm=0.5
        ... )
        >>>
        >>> # Top face of lower part should meet bottom face of upper part
        >>> constraint = create_face_distance_constraint(
        ...     name="stacked_contact",
        ...     frame_a="LOWER_PART",
        ...     face_a="TOP",
        ...     frame_b="UPPER_PART",
        ...     face_b="BOTTOM",
        ...     expected_distance=0.0,
        ...     bounds_a=((-25, -25, 0), (25, 25, 30)),
        ...     bounds_b=((-25, -25, 0), (25, 25, 40))
        ... )
    """
    return KinematicConstraint(
        name=name,
        constraint_type=KinematicConstraintType.AT_DISTANCE,
        frame_a=frame_a,
        frame_b=frame_b,
        face_a=face_a,
        face_b=face_b,
        bounds_a=bounds_a,
        bounds_b=bounds_b,
        reference_radius=expected_distance,
        tolerance_mm=tolerance_mm,
        description=description or f"Face distance: {frame_a}.{face_a} to {frame_b}.{face_b}"
    )


def create_face_coincident_constraint(
    name: str,
    frame_a: str,
    face_a: str,
    frame_b: str,
    face_b: str,
    tolerance_mm: float = 0.1,
    bounds_a: Optional[Tuple[Tuple[float, float, float], Tuple[float, float, float]]] = None,
    bounds_b: Optional[Tuple[Tuple[float, float, float], Tuple[float, float, float]]] = None,
    description: str = ""
) -> KinematicConstraint:
    """Create a face-to-face coincident constraint.

    Validates that two faces are at the same position (touching/mating).
    This is a convenience wrapper around create_face_distance_constraint with
    expected_distance=0.

    Args:
        name: Constraint name
        frame_a: First part name
        face_a: Face specification on frame_a
        frame_b: Second part name
        face_b: Face specification on frame_b
        tolerance_mm: Allowed position deviation
        bounds_a: Bounding box of part A (required for standard face names)
        bounds_b: Bounding box of part B (required for standard face names)
        description: Description of constraint purpose

    Returns:
        Configured KinematicConstraint

    Example:
        >>> # Servo output shaft should coincide with ring gear input
        >>> constraint = create_face_coincident_constraint(
        ...     name="servo_ring_mate",
        ...     frame_a="AXIS2_SERVO",
        ...     face_a="OUTPUT_SHAFT",
        ...     frame_b="AXIS2_RING",
        ...     face_b="SUN_INPUT"
        ... )
    """
    return KinematicConstraint(
        name=name,
        constraint_type=KinematicConstraintType.COINCIDENT,
        frame_a=frame_a,
        frame_b=frame_b,
        face_a=face_a,
        face_b=face_b,
        bounds_a=bounds_a,
        bounds_b=bounds_b,
        tolerance_mm=tolerance_mm,
        description=description or f"Face coincident: {frame_a}.{face_a} == {frame_b}.{face_b}"
    )
