"""Ring generation and female hole creation for manufacturing.

This module provides functions to create base and cradle rings
with female holes that receive terminal connectors from arcs.
"""

import math
from typing import Any, Dict, List, Optional, Tuple

from .data import SweptElementProvenance
from .connectors import (
    FIT_CLEARANCE,
    compute_connector_profile_dimensions,
    create_connector_region2d,
)


def create_ring_spine(
    radius: float,
    center: Tuple[float, float, float] = (0.0, 0.0, 0.0),
    tilt_angle_deg: float = 0.0,
    tilt_axis: str = "x",
    num_segments: int = 72,
) -> Dict[str, Any]:
    """Create a circular path3d for a ring.

    Args:
        radius: Radius of the ring (to centerline of profile)
        center: Center point of the ring (x, y, z)
        tilt_angle_deg: Tilt angle in degrees (0 = horizontal)
        tilt_axis: Axis to tilt around ("x", "y", or "z")
        num_segments: Number of line segments to approximate the circle

    Returns:
        yapCAD path3d dictionary with 'segments' key
    """
    # Build ring as a sequence of line segments in XY plane
    # This matches the path3d dict format used elsewhere in yapCAD
    segments = []
    angle_step = 2 * math.pi / num_segments

    for i in range(num_segments):
        angle_start = i * angle_step
        angle_end = (i + 1) * angle_step

        # Start and end points on circle in XY plane
        x0 = radius * math.cos(angle_start)
        y0 = radius * math.sin(angle_start)
        z0 = 0.0

        x1 = radius * math.cos(angle_end)
        y1 = radius * math.sin(angle_end)
        z1 = 0.0

        segments.append({
            'type': 'line',
            'start': [x0, y0, z0],
            'end': [x1, y1, z1]
        })

    ring_path = {'segments': segments}

    # Apply tilt rotation if needed
    if abs(tilt_angle_deg) > 0.001:
        ring_path = _rotate_path3d(ring_path, tilt_angle_deg, tilt_axis)

    # Apply translation to center
    if center != (0.0, 0.0, 0.0):
        ring_path = _translate_path3d(ring_path, center)

    return ring_path


def _rotate_path3d(
    path: Dict[str, Any],
    angle_deg: float,
    axis: str,
) -> Dict[str, Any]:
    """Rotate a path3d around an axis through the origin.

    Args:
        path: path3d dict with 'segments' key
        angle_deg: rotation angle in degrees
        axis: "x", "y", or "z"

    Returns:
        Rotated path3d dict
    """
    angle_rad = math.radians(angle_deg)
    c = math.cos(angle_rad)
    s = math.sin(angle_rad)

    def rotate_point(x, y, z):
        if axis.lower() == "x":
            return (x, y * c - z * s, y * s + z * c)
        elif axis.lower() == "y":
            return (x * c + z * s, y, -x * s + z * c)
        else:  # z axis
            return (x * c - y * s, x * s + y * c, z)

    new_segments = []
    for seg in path.get('segments', []):
        start = seg['start']
        end = seg['end']
        new_start = rotate_point(start[0], start[1], start[2])
        new_end = rotate_point(end[0], end[1], end[2])
        new_segments.append({
            'type': seg['type'],
            'start': list(new_start),
            'end': list(new_end)
        })

    return {'segments': new_segments}


def _translate_path3d(
    path: Dict[str, Any],
    offset: Tuple[float, float, float],
) -> Dict[str, Any]:
    """Translate a path3d by an offset.

    Args:
        path: path3d dict with 'segments' key
        offset: (dx, dy, dz) translation

    Returns:
        Translated path3d dict
    """
    dx, dy, dz = offset

    new_segments = []
    for seg in path.get('segments', []):
        start = seg['start']
        end = seg['end']
        new_segments.append({
            'type': seg['type'],
            'start': [start[0] + dx, start[1] + dy, start[2] + dz],
            'end': [end[0] + dx, end[1] + dy, end[2] + dz]
        })

    return {'segments': new_segments}


def create_ring_profile(
    outer_width: float,
    outer_height: float,
    wall_thickness: float,
) -> List:
    """Create a hollow rectangular profile (region2d) for a ring.

    Args:
        outer_width: Width of outer profile
        outer_height: Height of outer profile
        wall_thickness: Wall thickness

    Returns:
        yapCAD region2d (outer boundary, inner hole)
    """
    from yapcad.geom import line, point

    # Outer rectangle
    hw = outer_width / 2
    hh = outer_height / 2

    outer = [
        line(point(-hw, -hh, 0), point(hw, -hh, 0)),
        line(point(hw, -hh, 0), point(hw, hh, 0)),
        line(point(hw, hh, 0), point(-hw, hh, 0)),
        line(point(-hw, hh, 0), point(-hw, -hh, 0)),
    ]

    # Inner rectangle (hole)
    inner_hw = hw - wall_thickness
    inner_hh = hh - wall_thickness

    if inner_hw <= 0 or inner_hh <= 0:
        # Solid profile if wall too thick
        return [outer]

    # Inner winding is opposite (clockwise for hole)
    inner = [
        line(point(-inner_hw, -inner_hh, 0), point(-inner_hw, inner_hh, 0)),
        line(point(-inner_hw, inner_hh, 0), point(inner_hw, inner_hh, 0)),
        line(point(inner_hw, inner_hh, 0), point(inner_hw, -inner_hh, 0)),
        line(point(inner_hw, -inner_hh, 0), point(-inner_hw, -inner_hh, 0)),
    ]

    return [outer, inner]


def create_ring_solid(
    radius: float,
    outer_width: float,
    outer_height: float,
    wall_thickness: float,
    center: Tuple[float, float, float] = (0.0, 0.0, 0.0),
    tilt_angle_deg: float = 0.0,
    tilt_axis: str = "x",
    ring_id: str = "ring",
) -> Tuple[Any, SweptElementProvenance]:
    """Create a hollow box-beam ring with provenance tracking.

    Args:
        radius: Radius to centerline of profile
        outer_width: Width of beam profile
        outer_height: Height of beam profile
        wall_thickness: Wall thickness for hollow profile
        center: Center point of ring
        tilt_angle_deg: Tilt angle in degrees
        tilt_axis: Axis to tilt around
        ring_id: Identifier for the ring

    Returns:
        Tuple of (ring_solid, SweptElementProvenance)
    """
    from yapcad.geom3d_util import sweep_adaptive

    # Create ring spine
    spine = create_ring_spine(
        radius, center, tilt_angle_deg, tilt_axis
    )

    # Create hollow profile
    profile = create_ring_profile(outer_width, outer_height, wall_thickness)

    # Sweep profile along ring spine
    ring_solid = sweep_adaptive(
        profile[0],  # outer boundary
        spine,
        angle_threshold_deg=5.0,
        inner_profiles=[profile[1]] if len(profile) > 1 else None,
    )

    # Create provenance
    provenance = SweptElementProvenance(
        id=ring_id,
        operation="sweep_adaptive",
        outer_profile=profile,
        spine=spine,
        inner_profile=profile[1] if len(profile) > 1 else None,
        wall_thickness=wall_thickness,
        semantic_type="ring",
        metadata={
            'solid': ring_solid,
            'radius': radius,
            'center': center,
            'tilt_angle_deg': tilt_angle_deg,
        },
    )

    return ring_solid, provenance


def create_female_hole_solid(
    position: Tuple[float, float, float],
    direction: Tuple[float, float, float],
    outer_width: float,
    outer_height: float,
    wall_thickness: float,
    hole_depth: float,
    fit_clearance: float = FIT_CLEARANCE['press'],
) -> Any:
    """Create a solid to subtract from a ring for a female connector hole.

    The hole allows a male terminal connector to slot in.

    Args:
        position: Starting position of hole (on ring surface)
        direction: Direction hole extends (into ring)
        outer_width: Width of beam profile (for sizing hole)
        outer_height: Height of beam profile (for sizing hole)
        wall_thickness: Wall thickness of beam
        hole_depth: How deep the hole extends
        fit_clearance: Clearance to add for fit

    Returns:
        yapCAD solid representing the hole volume to subtract
    """
    from yapcad.geom3d_util import sweep_profile_along_path

    # Compute hole dimensions - should match connector profile but with
    # small additional clearance for the hole itself
    hole_clearance = fit_clearance / 2  # Extra clearance for hole
    conn_width, conn_height = compute_connector_profile_dimensions(
        outer_width, outer_height, wall_thickness, fit_clearance
    )
    # Add extra clearance for the female hole
    hole_width = conn_width + 2 * hole_clearance
    hole_height = conn_height + 2 * hole_clearance

    # Create hole profile
    hole_profile = create_connector_region2d(hole_width, hole_height)

    # Normalize direction
    dx, dy, dz = direction
    mag = math.sqrt(dx*dx + dy*dy + dz*dz)
    if mag > 0:
        dx, dy, dz = dx/mag, dy/mag, dz/mag

    # Create linear path for hole (using dict format)
    end_pos = (
        position[0] + dx * hole_depth,
        position[1] + dy * hole_depth,
        position[2] + dz * hole_depth,
    )

    hole_spine = {
        'segments': [{
            'type': 'line',
            'start': [position[0], position[1], position[2]],
            'end': [end_pos[0], end_pos[1], end_pos[2]]
        }]
    }

    # Sweep hole profile to create hole solid
    hole_solid = sweep_profile_along_path(
        hole_profile[0],  # outer boundary
        hole_spine,
    )

    return hole_solid


def compute_arc_attachment_point(
    ring_radius: float,
    attachment_angle_deg: float,
    ring_center: Tuple[float, float, float] = (0.0, 0.0, 0.0),
    ring_tilt_deg: float = 0.0,
    ring_tilt_axis: str = "x",
) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:
    """Compute where an arc attaches to a ring.

    Returns the position on the ring centerline and the radial direction
    pointing inward (toward ring center).

    Args:
        ring_radius: Radius of ring
        attachment_angle_deg: Angle around ring where arc attaches
        ring_center: Center of ring
        ring_tilt_deg: Ring tilt angle
        ring_tilt_axis: Axis ring is tilted around

    Returns:
        Tuple of (position, inward_direction)
    """
    # Compute position on untilted ring in XY plane
    angle_rad = math.radians(attachment_angle_deg)
    px = ring_radius * math.cos(angle_rad)
    py = ring_radius * math.sin(angle_rad)
    pz = 0.0

    # Radial direction (pointing inward, toward center)
    dx = -math.cos(angle_rad)
    dy = -math.sin(angle_rad)
    dz = 0.0

    # Apply tilt rotation
    if abs(ring_tilt_deg) > 0.001:
        px, py, pz, dx, dy, dz = _rotate_point_and_direction(
            px, py, pz, dx, dy, dz, ring_tilt_deg, ring_tilt_axis
        )

    # Apply translation
    position = (
        px + ring_center[0],
        py + ring_center[1],
        pz + ring_center[2],
    )

    return position, (dx, dy, dz)


def _rotate_point_and_direction(
    px: float, py: float, pz: float,
    dx: float, dy: float, dz: float,
    angle_deg: float,
    axis: str,
) -> Tuple[float, float, float, float, float, float]:
    """Rotate a point and direction vector around an axis."""
    angle_rad = math.radians(angle_deg)
    c = math.cos(angle_rad)
    s = math.sin(angle_rad)

    if axis.lower() == "x":
        # Rotate around X axis
        py_new = py * c - pz * s
        pz_new = py * s + pz * c
        py, pz = py_new, pz_new

        dy_new = dy * c - dz * s
        dz_new = dy * s + dz * c
        dy, dz = dy_new, dz_new
    elif axis.lower() == "y":
        # Rotate around Y axis
        px_new = px * c + pz * s
        pz_new = -px * s + pz * c
        px, pz = px_new, pz_new

        dx_new = dx * c + dz * s
        dz_new = -dx * s + dz * c
        dx, dz = dx_new, dz_new
    else:  # z axis
        # Rotate around Z axis
        px_new = px * c - py * s
        py_new = px * s + py * c
        px, py = px_new, py_new

        dx_new = dx * c - dy * s
        dy_new = dx * s + dy * c
        dx, dy = dx_new, dy_new

    return px, py, pz, dx, dy, dz


def add_female_holes_to_ring(
    ring_solid: Any,
    attachment_points: List[Tuple[Tuple[float, float, float], Tuple[float, float, float]]],
    outer_width: float,
    outer_height: float,
    wall_thickness: float,
    hole_depth: Optional[float] = None,
    fit_clearance: float = FIT_CLEARANCE['press'],
) -> Any:
    """Subtract female holes at arc attachment positions.

    Args:
        ring_solid: Ring solid to modify
        attachment_points: List of (position, direction) tuples
        outer_width: Beam profile width
        outer_height: Beam profile height
        wall_thickness: Wall thickness
        hole_depth: Depth of holes (auto-computed if None)
        fit_clearance: Fit clearance for holes

    Returns:
        Ring solid with female holes subtracted
    """
    from yapcad.geom3d import solid_boolean

    if hole_depth is None:
        # Default hole depth matches terminal connector length
        from .connectors import DEFAULT_LENGTH_FACTOR
        hole_depth = max(outer_width, outer_height) * DEFAULT_LENGTH_FACTOR / 2

    result = ring_solid

    for position, direction in attachment_points:
        hole_solid = create_female_hole_solid(
            position,
            direction,
            outer_width,
            outer_height,
            wall_thickness,
            hole_depth,
            fit_clearance,
        )

        try:
            result = solid_boolean(result, hole_solid, 'difference')
        except Exception as e:
            # Log warning but continue
            print(f"Warning: Failed to subtract hole at {position}: {e}")

    return result


def trim_segment_against_ring(
    segment_solid: Any,
    ring_solid: Any,
) -> Any:
    """Trim a segment by subtracting the ring solid from it.

    This creates a clean interface where the arc meets the ring,
    removing any overlap between the arc segment and ring geometry.

    Args:
        segment_solid: Arc segment to trim
        ring_solid: Ring solid to subtract

    Returns:
        Trimmed segment solid
    """
    from yapcad.geom3d import solid_boolean
    from yapcad.brep import brep_from_solid, attach_brep_to_solid

    result = solid_boolean(segment_solid, ring_solid, 'difference')

    # Ensure BREP data is attached to the result
    brep = brep_from_solid(result)
    if brep is not None:
        attach_brep_to_solid(result, brep)

    return result


def compute_ring_cuts_avoiding_holes(
    ring_circumference: float,
    max_segment_length: float,
    hole_angles_deg: List[float],
    min_distance_from_hole: float = 30.0,
    *,
    target_segments: Optional[int] = None,
) -> List[float]:
    """Compute cut parameters for ring that avoid female hole positions.

    Produces evenly-spaced segments with cuts placed to avoid holes.
    The algorithm finds cut positions that result in equal-length segments
    (relative to t=0) while keeping cuts away from hole locations.

    For equal segments on a closed ring, we need cuts at positions that
    divide [0, 1) into equal parts. With N segments, cuts should be at
    1/N, 2/N, ..., (N-1)/N. If holes are at these positions, we find
    alternate cut positions that still produce equal segments.

    Args:
        ring_circumference: Total circumference of ring in mm
        max_segment_length: Maximum segment length in mm
        hole_angles_deg: Angles where holes are located (0-360)
        min_distance_from_hole: Minimum distance from hole center in mm
        target_segments: If specified, use this many segments instead of
            computing from max_segment_length. Useful when build volume
            constraints allow fewer segments than the length calculation.

    Returns:
        List of cut parameters (0-1) that avoid hole locations
    """
    # Compute number of segments needed
    if target_segments is not None:
        num_segments = target_segments
    else:
        num_segments = math.ceil(ring_circumference / max_segment_length)

    if num_segments < 2:
        return []  # No cuts needed

    num_cuts = num_segments - 1

    # Convert hole angles to parameters (0-360 -> 0-1)
    hole_params = sorted([angle / 360.0 for angle in hole_angles_deg])

    # Convert min_distance to parameter space
    min_param_distance = min_distance_from_hole / ring_circumference

    # Build exclusion zones around holes
    exclusion_zones = []
    for hole_param in hole_params:
        zone_start = hole_param - min_param_distance
        zone_end = hole_param + min_param_distance
        exclusion_zones.append((zone_start, zone_end))

    def is_cut_valid(cut_param: float) -> bool:
        """Check if a cut position avoids all holes."""
        for zone_start, zone_end in exclusion_zones:
            # Handle wrap-around
            if zone_start < 0:
                # Zone wraps around 0
                if cut_param < zone_end or cut_param > (1.0 + zone_start):
                    return False
            elif zone_end > 1.0:
                # Zone wraps around 1
                if cut_param > zone_start or cut_param < (zone_end - 1.0):
                    return False
            else:
                # Normal case
                if zone_start <= cut_param <= zone_end:
                    return False
        return True

    def compute_segment_sizes(cuts: List[float]) -> List[float]:
        """Compute segment sizes from sorted cut positions."""
        if not cuts:
            return [1.0]
        sorted_cuts = sorted(cuts)
        sizes = []
        # First segment: from 0 to first cut
        sizes.append(sorted_cuts[0])
        # Middle segments
        for i in range(1, len(sorted_cuts)):
            sizes.append(sorted_cuts[i] - sorted_cuts[i-1])
        # Last segment: from last cut to 1.0
        sizes.append(1.0 - sorted_cuts[-1])
        return sizes

    def segment_size_variance(cuts: List[float]) -> float:
        """Compute variance of segment sizes (lower is better)."""
        sizes = compute_segment_sizes(cuts)
        mean_size = sum(sizes) / len(sizes)
        return sum((s - mean_size) ** 2 for s in sizes)

    # For equal segments, cuts MUST be at exactly k/num_segments for k=1..num_cuts
    # E.g., for 3 segments: cuts at 1/3 and 2/3 give equal segments [0,1/3], [1/3,2/3], [2/3,1]
    #
    # If holes block these positions, we find cuts as close as possible to minimize variance.

    segment_size = 1.0 / num_segments

    # Ideal cut positions for equal segments
    ideal_cuts = [segment_size * (i + 1) for i in range(num_cuts)]

    # First, check if ideal cuts are all valid
    if all(is_cut_valid(c) for c in ideal_cuts):
        return ideal_cuts

    # Ideal cuts are blocked by holes. Find cuts as close as possible to ideal positions.
    # Strategy: for each ideal cut position, find the nearest valid position.

    def find_nearest_valid_position(ideal_pos: float, search_range: float = 0.5) -> Optional[float]:
        """Find the nearest valid position to ideal_pos within search_range."""
        if is_cut_valid(ideal_pos):
            return ideal_pos

        # Search in both directions with fine resolution
        best_pos = None
        best_dist = float('inf')

        # Search 1000 positions within range
        for i in range(1000):
            offset = (i / 1000.0) * search_range

            # Try position above ideal
            pos_above = ideal_pos + offset
            if pos_above < 1.0 and is_cut_valid(pos_above):
                dist = abs(pos_above - ideal_pos)
                if dist < best_dist:
                    best_dist = dist
                    best_pos = pos_above
                    break  # Found closest valid above

            # Try position below ideal
            pos_below = ideal_pos - offset
            if pos_below > 0.0 and is_cut_valid(pos_below):
                dist = abs(pos_below - ideal_pos)
                if dist < best_dist:
                    best_dist = dist
                    best_pos = pos_below
                    break  # Found closest valid below

        return best_pos

    # Find nearest valid position for each ideal cut
    adjusted_cuts = []
    for ideal_pos in ideal_cuts:
        nearest = find_nearest_valid_position(ideal_pos)
        if nearest is not None:
            adjusted_cuts.append(nearest)
        else:
            # No valid position found for this cut - fall back to search
            break

    if len(adjusted_cuts) == num_cuts:
        # Verify no duplicates and proper ordering
        adjusted_cuts = sorted(set(adjusted_cuts))
        if len(adjusted_cuts) == num_cuts:
            return adjusted_cuts

    # Fallback: exhaustive search for minimum variance solution
    # Find valid regions (gaps between exclusion zones)
    exclusion_points = []
    for zone_start, zone_end in exclusion_zones:
        if zone_start < 0:
            exclusion_points.append((zone_start + 1.0, 1.0))
            exclusion_points.append((0.0, zone_end))
        elif zone_end > 1.0:
            exclusion_points.append((zone_start, 1.0))
            exclusion_points.append((0.0, zone_end - 1.0))
        else:
            exclusion_points.append((zone_start, zone_end))

    exclusion_points.sort()
    merged = []
    for start, end in exclusion_points:
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))

    valid_regions = []
    if not merged:
        valid_regions = [(0.0, 1.0)]
    else:
        if merged[0][0] > 0:
            valid_regions.append((0.0, merged[0][0]))
        for i in range(len(merged) - 1):
            valid_regions.append((merged[i][1], merged[i+1][0]))
        if merged[-1][1] < 1.0:
            valid_regions.append((merged[-1][1], 1.0))

    # Generate candidate positions: boundaries of valid regions and ideal cut positions
    candidate_cuts = []

    # Add boundary positions (just inside valid regions)
    for region_start, region_end in valid_regions:
        candidate_cuts.append(region_start + 0.001)
        candidate_cuts.append(region_end - 0.001)
        # Also add midpoint
        candidate_cuts.append((region_start + region_end) / 2)

    # Add positions near ideal cuts
    for ideal in ideal_cuts:
        for offset in [0.0, 0.01, 0.02, 0.05, 0.1, -0.01, -0.02, -0.05, -0.1]:
            pos = ideal + offset
            if 0 < pos < 1:
                candidate_cuts.append(pos)

    # Filter to valid positions only
    candidate_cuts = sorted(set(c for c in candidate_cuts if 0 < c < 1 and is_cut_valid(c)))

    if len(candidate_cuts) >= num_cuts:
        from itertools import combinations
        best_cuts = None
        best_variance = float('inf')

        for combo in combinations(candidate_cuts, num_cuts):
            cuts = sorted(combo)
            variance = segment_size_variance(cuts)
            if variance < best_variance:
                best_variance = variance
                best_cuts = cuts

        if best_cuts:
            return list(best_cuts)

    # Last resort: return evenly spaced cuts, even if they hit holes
    return [segment_size * (i + 1) for i in range(num_cuts)]
