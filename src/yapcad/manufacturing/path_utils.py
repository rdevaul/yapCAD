"""Path3D utilities for manufacturing post-processing.

This module provides functions for evaluating and manipulating path3d
objects, which are essential for beam segmentation operations.
"""

import math
from typing import Any, Dict, List, Optional, Tuple


def _normalize_vector(v: List[float]) -> List[float]:
    """Normalize a 3D vector to unit length."""
    length = math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    if length < 1e-10:
        return [0.0, 0.0, 1.0]  # Default up vector
    return [v[0]/length, v[1]/length, v[2]/length]


def _vector_subtract(a: List[float], b: List[float]) -> List[float]:
    """Subtract two 3D vectors."""
    return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]


def _vector_add(a: List[float], b: List[float]) -> List[float]:
    """Add two 3D vectors."""
    return [a[0] + b[0], a[1] + b[1], a[2] + b[2]]


def _vector_scale(v: List[float], s: float) -> List[float]:
    """Scale a 3D vector."""
    return [v[0] * s, v[1] * s, v[2] * s]


def _cross_product(a: List[float], b: List[float]) -> List[float]:
    """Cross product of two 3D vectors."""
    return [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    ]


def _lerp_point(a: List[float], b: List[float], t: float) -> List[float]:
    """Linear interpolation between two points."""
    return [
        a[0] + t * (b[0] - a[0]),
        a[1] + t * (b[1] - a[1]),
        a[2] + t * (b[2] - a[2])
    ]


def _evaluate_line_segment(
    start: List[float],
    end: List[float],
    t: float
) -> Tuple[List[float], List[float]]:
    """Evaluate a line segment at parameter t.

    Args:
        start: Start point [x, y, z]
        end: End point [x, y, z]
        t: Parameter in [0, 1]

    Returns:
        Tuple of (point, tangent) where tangent is unit vector
    """
    point = _lerp_point(start, end, t)
    tangent = _normalize_vector(_vector_subtract(end, start))
    return point, tangent


def _evaluate_arc_segment(
    center: List[float],
    start: List[float],
    end: List[float],
    normal: List[float],
    t: float
) -> Tuple[List[float], List[float]]:
    """Evaluate an arc segment at parameter t.

    Args:
        center: Arc center point
        start: Start point on arc
        end: End point on arc
        normal: Arc plane normal (defines rotation direction)
        t: Parameter in [0, 1]

    Returns:
        Tuple of (point, tangent) where tangent is unit vector
    """
    # Vector from center to start
    r_start = _vector_subtract(start, center)
    radius = math.sqrt(r_start[0]**2 + r_start[1]**2 + r_start[2]**2)

    # Vector from center to end
    r_end = _vector_subtract(end, center)

    # Normalize radius vectors
    r_start_norm = _normalize_vector(r_start)
    r_end_norm = _normalize_vector(r_end)

    # Calculate angle using dot product
    dot = r_start_norm[0]*r_end_norm[0] + r_start_norm[1]*r_end_norm[1] + r_start_norm[2]*r_end_norm[2]
    dot = max(-1.0, min(1.0, dot))
    total_angle = math.acos(dot)

    # Check if we need to go the "long way" around
    # Use cross product to determine direction
    cross = _cross_product(r_start_norm, r_end_norm)
    normal_norm = _normalize_vector(normal)
    sign = cross[0]*normal_norm[0] + cross[1]*normal_norm[1] + cross[2]*normal_norm[2]
    if sign < 0:
        total_angle = 2 * math.pi - total_angle

    # Interpolate angle
    angle = t * total_angle

    # Rodrigues' rotation formula to rotate r_start around normal
    cos_a = math.cos(angle)
    sin_a = math.sin(angle)

    # k x r (cross product of normal and radius vector)
    k_cross_r = _cross_product(normal_norm, r_start_norm)
    # k . r (dot product)
    k_dot_r = normal_norm[0]*r_start_norm[0] + normal_norm[1]*r_start_norm[1] + normal_norm[2]*r_start_norm[2]

    # Rotated radius direction: r*cos(a) + (k x r)*sin(a) + k*(k.r)*(1-cos(a))
    rotated = [
        r_start_norm[i] * cos_a + k_cross_r[i] * sin_a + normal_norm[i] * k_dot_r * (1 - cos_a)
        for i in range(3)
    ]

    # Point on arc
    point = _vector_add(center, _vector_scale(rotated, radius))

    # Tangent is perpendicular to radius, in direction of rotation
    # tangent = normal x rotated_radius_direction
    tangent = _normalize_vector(_cross_product(normal_norm, rotated))

    return point, tangent


def evaluate_path3d_at_t(
    path3d: Dict[str, Any],
    t: float
) -> Tuple[List[float], List[float]]:
    """Evaluate a path3d at a global parameter value.

    Args:
        path3d: Dict with 'segments' list of line/arc segments
        t: Global parameter in [0, 1] across entire path

    Returns:
        Tuple of (point, tangent) at parameter t
        - point: [x, y, z] position on path
        - tangent: [tx, ty, tz] unit tangent vector

    Raises:
        ValueError: If path has no segments or t is out of range
    """
    segments = path3d.get('segments', [])
    if not segments:
        raise ValueError("Path has no segments")

    # Clamp t to valid range
    t = max(0.0, min(1.0, t))

    # Convert global t to segment index and local t
    n_segments = len(segments)
    if t >= 1.0:
        # At the end, evaluate last segment at t=1
        seg_idx = n_segments - 1
        local_t = 1.0
    else:
        # Find which segment and local parameter
        scaled = t * n_segments
        seg_idx = int(scaled)
        local_t = scaled - seg_idx

    seg = segments[seg_idx]
    seg_type = seg.get('type', 'line')

    if seg_type == 'line':
        return _evaluate_line_segment(seg['start'], seg['end'], local_t)
    elif seg_type == 'arc':
        normal = seg.get('normal', [0, 0, 1])
        return _evaluate_arc_segment(
            seg['center'], seg['start'], seg['end'], normal, local_t
        )
    else:
        raise ValueError(f"Unknown segment type: {seg_type}")


def compute_cut_plane(
    path3d: Dict[str, Any],
    t: float
) -> Tuple[List[float], List[float]]:
    """Compute the cut plane at a parameter along the path.

    The cut plane is perpendicular to the path tangent at the given parameter.

    Args:
        path3d: Dict with 'segments' list
        t: Global parameter in [0, 1]

    Returns:
        Tuple of (point, normal) defining the plane
        - point: A point on the plane (the path point at t)
        - normal: Unit normal to the plane (the path tangent at t)
    """
    point, tangent = evaluate_path3d_at_t(path3d, t)
    return point, tangent


def extract_sub_path(
    path3d: Dict[str, Any],
    t_start: float,
    t_end: float
) -> Dict[str, Any]:
    """Extract a portion of a path3d between two parameters.

    Args:
        path3d: Original path dict
        t_start: Start parameter in [0, 1]
        t_end: End parameter in [0, 1], must be > t_start

    Returns:
        New path3d dict containing only the portion between t_start and t_end

    Raises:
        ValueError: If t_start >= t_end or parameters out of range
    """
    if t_start >= t_end:
        raise ValueError(f"t_start ({t_start}) must be less than t_end ({t_end})")

    t_start = max(0.0, min(1.0, t_start))
    t_end = max(0.0, min(1.0, t_end))

    segments = path3d.get('segments', [])
    if not segments:
        return {'segments': []}

    n_segments = len(segments)
    new_segments = []

    # Find start and end segment indices
    start_scaled = t_start * n_segments
    end_scaled = t_end * n_segments

    start_seg_idx = int(start_scaled)
    end_seg_idx = min(int(end_scaled), n_segments - 1)

    start_local_t = start_scaled - start_seg_idx
    end_local_t = end_scaled - end_seg_idx if end_seg_idx < n_segments else 1.0

    for seg_idx in range(start_seg_idx, end_seg_idx + 1):
        seg = segments[seg_idx]
        seg_type = seg.get('type', 'line')

        # Determine local t range for this segment
        if seg_idx == start_seg_idx:
            local_t_min = start_local_t
        else:
            local_t_min = 0.0

        if seg_idx == end_seg_idx:
            local_t_max = end_local_t
        else:
            local_t_max = 1.0

        # Skip if this segment contributes nothing
        if local_t_min >= local_t_max:
            continue

        # Extract the relevant portion
        if seg_type == 'line':
            new_start = _lerp_point(seg['start'], seg['end'], local_t_min)
            new_end = _lerp_point(seg['start'], seg['end'], local_t_max)
            new_segments.append({
                'type': 'line',
                'start': new_start,
                'end': new_end
            })
        elif seg_type == 'arc':
            # For arcs, we need to find the points at the local parameters
            normal = seg.get('normal', [0, 0, 1])
            new_start_pt, _ = _evaluate_arc_segment(
                seg['center'], seg['start'], seg['end'], normal, local_t_min
            )
            new_end_pt, _ = _evaluate_arc_segment(
                seg['center'], seg['start'], seg['end'], normal, local_t_max
            )
            new_segments.append({
                'type': 'arc',
                'center': seg['center'],
                'start': new_start_pt,
                'end': new_end_pt,
                'normal': normal
            })

    return {'segments': new_segments}


def path_length(path3d: Dict[str, Any]) -> float:
    """Compute the total arc length of a path3d.

    Args:
        path3d: Dict with 'segments' list

    Returns:
        Total length in same units as path coordinates
    """
    segments = path3d.get('segments', [])
    total = 0.0

    for seg in segments:
        seg_type = seg.get('type', 'line')

        if seg_type == 'line':
            diff = _vector_subtract(seg['end'], seg['start'])
            total += math.sqrt(diff[0]**2 + diff[1]**2 + diff[2]**2)
        elif seg_type == 'arc':
            # Arc length = radius * angle
            center = seg['center']
            start = seg['start']
            end = seg['end']
            normal = seg.get('normal', [0, 0, 1])

            r_start = _vector_subtract(start, center)
            radius = math.sqrt(r_start[0]**2 + r_start[1]**2 + r_start[2]**2)

            r_end = _vector_subtract(end, center)
            r_start_norm = _normalize_vector(r_start)
            r_end_norm = _normalize_vector(r_end)

            dot = r_start_norm[0]*r_end_norm[0] + r_start_norm[1]*r_end_norm[1] + r_start_norm[2]*r_end_norm[2]
            dot = max(-1.0, min(1.0, dot))
            angle = math.acos(dot)

            # Check direction
            cross = _cross_product(r_start_norm, r_end_norm)
            normal_norm = _normalize_vector(normal)
            sign = cross[0]*normal_norm[0] + cross[1]*normal_norm[1] + cross[2]*normal_norm[2]
            if sign < 0:
                angle = 2 * math.pi - angle

            total += radius * angle

    return total


def length_to_parameter(path3d: Dict[str, Any], arc_length: float) -> float:
    """Convert an arc length to a global parameter.

    Args:
        path3d: Dict with 'segments' list
        arc_length: Distance along the path from start

    Returns:
        Global parameter t in [0, 1]
    """
    total_length = path_length(path3d)
    if total_length < 1e-10:
        return 0.0
    return min(1.0, max(0.0, arc_length / total_length))


def parameter_to_length(path3d: Dict[str, Any], t: float) -> float:
    """Convert a global parameter to arc length.

    Args:
        path3d: Dict with 'segments' list
        t: Global parameter in [0, 1]

    Returns:
        Arc length from path start to parameter t
    """
    # Extract sub-path and compute its length
    sub_path = extract_sub_path(path3d, 0.0, t)
    return path_length(sub_path)
