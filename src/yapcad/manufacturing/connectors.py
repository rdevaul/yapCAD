"""Interior connector generation for beam segmentation.

This module provides functions to create interior connectors that join
segmented swept elements. Connectors fit inside hollow profiles and
provide structural continuity across cut planes.
"""

import math
from typing import Any, Dict, List, Optional, Tuple, TYPE_CHECKING

if TYPE_CHECKING:
    from .data import ConnectorSpec, SweptElementProvenance

from .path_utils import (
    extract_sub_path,
    path_length,
    length_to_parameter,
    parameter_to_length,
)


# Default fit clearance values (mm per side)
FIT_CLEARANCE = {
    'press': 0.18,      # Press-fit (structural)
    'slip': 0.30,       # Slip-fit (easy assembly)
    'loose': 0.45,      # Loose-fit (adjustable)
}

# Connector length factor (multiple of largest profile dimension)
DEFAULT_LENGTH_FACTOR = 3.0


def offset_rectangular_profile(
    width: float,
    height: float,
    clearance: float
) -> Tuple[float, float]:
    """Offset a rectangular profile inward by clearance.

    Args:
        width: Original profile width
        height: Original profile height
        clearance: Amount to shrink per side

    Returns:
        Tuple of (new_width, new_height)

    Raises:
        ValueError: If clearance would result in zero or negative dimensions
    """
    new_width = width - 2 * clearance
    new_height = height - 2 * clearance

    if new_width <= 0 or new_height <= 0:
        raise ValueError(
            f"Clearance {clearance} too large for profile {width}x{height}. "
            f"Would result in {new_width}x{new_height}"
        )

    return new_width, new_height


def compute_inner_profile_dimensions(
    outer_width: float,
    outer_height: float,
    wall_thickness: float
) -> Tuple[float, float]:
    """Compute inner void dimensions from outer profile and wall thickness.

    Args:
        outer_width: Outer profile width
        outer_height: Outer profile height
        wall_thickness: Wall thickness

    Returns:
        Tuple of (inner_width, inner_height)
    """
    inner_width = outer_width - 2 * wall_thickness
    inner_height = outer_height - 2 * wall_thickness

    if inner_width <= 0 or inner_height <= 0:
        raise ValueError(
            f"Wall thickness {wall_thickness} too large for profile "
            f"{outer_width}x{outer_height}"
        )

    return inner_width, inner_height


def compute_connector_profile_dimensions(
    outer_width: float,
    outer_height: float,
    wall_thickness: float,
    fit_clearance: float
) -> Tuple[float, float]:
    """Compute connector cross-section dimensions.

    The connector fits inside the hollow interior with specified clearance.

    Args:
        outer_width: Outer profile width
        outer_height: Outer profile height
        wall_thickness: Wall thickness of hollow profile
        fit_clearance: Clearance per side for desired fit

    Returns:
        Tuple of (connector_width, connector_height)
    """
    inner_w, inner_h = compute_inner_profile_dimensions(
        outer_width, outer_height, wall_thickness
    )
    return offset_rectangular_profile(inner_w, inner_h, fit_clearance)


def compute_connector_length(
    profile_width: float,
    profile_height: float,
    spine: Dict[str, Any],
    center_parameter: float,
    *,
    length_factor: float = DEFAULT_LENGTH_FACTOR,
    min_arc_degrees: float = 15.0
) -> float:
    """Compute the appropriate connector length.

    For straight sections, length is based on profile size.
    For curved sections, length may be extended to span a minimum arc.

    Args:
        profile_width: Profile width
        profile_height: Profile height
        spine: Path3d dict
        center_parameter: Where connector is centered (t in [0, 1])
        length_factor: Multiple of largest profile dimension
        min_arc_degrees: Minimum arc span for curved sections

    Returns:
        Connector length in mm
    """
    max_dim = max(profile_width, profile_height)
    base_length = length_factor * max_dim

    # Check curvature at cut point
    # For now, use the base length
    # Future: analyze local curvature and extend for tight curves

    return base_length


def create_connector_region2d(
    width: float,
    height: float,
    *,
    corner_radius: float = 0.0
) -> List:
    """Create a rectangular region2d for the connector profile.

    Args:
        width: Connector width
        height: Connector height
        corner_radius: Optional fillet radius for corners

    Returns:
        yapCAD region2d (list of polylines with proper winding)
    """
    from yapcad.geom import line, point

    hw = width / 2
    hh = height / 2

    if corner_radius > 0:
        # Rectangle with filleted corners - use RoundRect from poly
        from yapcad.poly import RoundRect
        poly = RoundRect(width, height, corner_radius * 2)  # chamf = diameter
        outline = poly.geom  # geom is a property that returns the geometry list
    else:
        # Simple rectangle centered at origin
        # yapCAD points have format [x, y, z, 1]
        outline = [
            line(point(-hw, -hh, 0), point(hw, -hh, 0)),
            line(point(hw, -hh, 0), point(hw, hh, 0)),
            line(point(hw, hh, 0), point(-hw, hh, 0)),
            line(point(-hw, hh, 0), point(-hw, -hh, 0)),
        ]

    return [outline]


def create_interior_connector(
    outer_profile_width: float,
    outer_profile_height: float,
    spine: Dict[str, Any],
    center_parameter: float,
    *,
    wall_thickness: float,
    connector_length: Optional[float] = None,
    fit_clearance: float = FIT_CLEARANCE['press'],
    corner_radius: float = 0.0,
) -> Any:
    """Create an interior connector solid.

    The connector fits inside the hollow interior of a swept beam element,
    spanning across a cut plane to join two segments.

    Args:
        outer_profile_width: Width of outer beam profile
        outer_profile_height: Height of outer beam profile
        spine: Path3d that the original beam follows
        center_parameter: Location along spine (t in [0, 1]) where cut occurs
        wall_thickness: Wall thickness of hollow beam
        connector_length: Length of connector (auto-computed if None)
        fit_clearance: Clearance for desired fit (mm per side)
        corner_radius: Optional fillet radius for connector corners

    Returns:
        yapCAD solid representing the connector
    """
    from yapcad.geom3d_util import sweep_adaptive

    # Compute connector profile dimensions
    conn_width, conn_height = compute_connector_profile_dimensions(
        outer_profile_width, outer_profile_height,
        wall_thickness, fit_clearance
    )

    # Compute connector length if not specified
    if connector_length is None:
        connector_length = compute_connector_length(
            outer_profile_width, outer_profile_height,
            spine, center_parameter
        )

    # Create connector profile (region2d)
    connector_profile = create_connector_region2d(
        conn_width, conn_height, corner_radius=corner_radius
    )

    # Extract spine segment for connector
    # Connector extends half-length on each side of center
    total_spine_length = path_length(spine)
    half_connector_length = connector_length / 2

    # Convert lengths to parameters
    center_length = parameter_to_length(spine, center_parameter)
    start_length = max(0, center_length - half_connector_length)
    end_length = min(total_spine_length, center_length + half_connector_length)

    start_t = length_to_parameter(spine, start_length)
    end_t = length_to_parameter(spine, end_length)

    # Extract sub-spine
    connector_spine = extract_sub_path(spine, start_t, end_t)

    # Sweep connector profile along spine segment
    # Note: sweep_adaptive expects the outer boundary polyline, not the full region2d
    connector_solid = sweep_adaptive(
        connector_profile[0],  # Extract outer boundary polyline from region2d
        connector_spine,
        angle_threshold_deg=5.0
    )

    return connector_solid


def compute_connector_spec(
    element_id: str,
    cut_parameter: float,
    outer_width: float,
    outer_height: float,
    wall_thickness: float,
    spine: Dict[str, Any],
    *,
    fit_clearance: float = FIT_CLEARANCE['press'],
    connector_id: Optional[str] = None,
) -> 'ConnectorSpec':
    """Compute full connector specification for a cut point.

    Args:
        element_id: ID of parent swept element
        cut_parameter: Where cut occurs (t in [0, 1])
        outer_width: Outer profile width
        outer_height: Outer profile height
        wall_thickness: Wall thickness
        spine: Path3d of parent element
        fit_clearance: Desired fit clearance
        connector_id: Optional ID (auto-generated if None)

    Returns:
        ConnectorSpec object with all computed values
    """
    from .data import ConnectorSpec

    length = compute_connector_length(
        outer_width, outer_height, spine, cut_parameter
    )

    if connector_id is None:
        connector_id = f"{element_id}_conn_{int(cut_parameter * 100)}"

    return ConnectorSpec(
        id=connector_id,
        parent_element_id=element_id,
        center_parameter=cut_parameter,
        length=length,
        fit_clearance=fit_clearance,
        profile_type="box"
    )


def create_terminal_connector(
    outer_profile_width: float,
    outer_profile_height: float,
    spine: Dict[str, Any],
    end: str,  # "start" or "end"
    *,
    wall_thickness: float,
    connector_length: Optional[float] = None,
    fit_clearance: float = FIT_CLEARANCE['press'],
    corner_radius: float = 0.0,
) -> Any:
    """Create a terminal connector tab at the start or end of a spine.

    Terminal connectors extend outward from arc endpoints and are designed
    to slot into female holes in rings or other structures.

    Args:
        outer_profile_width: Width of outer beam profile
        outer_profile_height: Height of outer beam profile
        spine: Path3d that the beam follows
        end: "start" for t=0 end, "end" for t=1 end
        wall_thickness: Wall thickness of hollow beam
        connector_length: Length of connector tab (auto-computed if None)
        fit_clearance: Clearance for desired fit (mm per side)
        corner_radius: Optional fillet radius for connector corners

    Returns:
        yapCAD solid representing the terminal connector tab
    """
    from yapcad.geom3d_util import sweep_profile_along_path

    if end not in ("start", "end"):
        raise ValueError(f"end must be 'start' or 'end', got '{end}'")

    # Compute connector profile dimensions (same as interior connector)
    conn_width, conn_height = compute_connector_profile_dimensions(
        outer_profile_width, outer_profile_height,
        wall_thickness, fit_clearance
    )

    # Compute connector length if not specified
    if connector_length is None:
        # For terminal connectors, use half the interior connector length
        # since they only extend one direction
        connector_length = compute_connector_length(
            outer_profile_width, outer_profile_height,
            spine, 0.5  # Use midpoint for length calculation
        ) / 2

    # Create connector profile
    connector_profile = create_connector_region2d(
        conn_width, conn_height, corner_radius=corner_radius
    )

    # Extract endpoint position and tangent from spine
    endpoint, tangent = _get_spine_endpoint_and_tangent(spine, end)

    # Create linear spine extending outward from endpoint
    # Tangent points along spine direction, so:
    # - At "start" (t=0), tangent points toward t=1, we want to extend OPPOSITE
    # - At "end" (t=1), tangent points toward t=1 (continuation), we want same dir
    if end == "start":
        # Extend in opposite direction of tangent (outward from start)
        direction = [-tangent[0], -tangent[1], -tangent[2]]
    else:
        # Extend in same direction as tangent (outward from end)
        direction = tangent

    # Create linear path from endpoint extending outward
    end_point = [
        endpoint[0] + direction[0] * connector_length,
        endpoint[1] + direction[1] * connector_length,
        endpoint[2] + direction[2] * connector_length,
    ]

    # Build linear path3d as dict (same format as rings.py)
    connector_spine = {
        'segments': [{
            'type': 'line',
            'start': [endpoint[0], endpoint[1], endpoint[2]],
            'end': [end_point[0], end_point[1], end_point[2]]
        }]
    }

    # Sweep connector profile along linear spine
    connector_solid = sweep_profile_along_path(
        connector_profile[0],  # Outer boundary polyline
        connector_spine,
    )

    return connector_solid


def _get_spine_endpoint_and_tangent(
    spine: Dict[str, Any],
    end: str,
) -> Tuple[List[float], List[float]]:
    """Extract endpoint position and tangent direction from spine.

    Args:
        spine: Path3d dict
        end: "start" for t=0, "end" for t=1

    Returns:
        Tuple of (position [x,y,z], tangent [dx,dy,dz] normalized)
    """
    from .path_utils import evaluate_path3d_at_t

    if end == "start":
        t = 0.0
    else:
        t = 1.0

    # Get position and tangent at endpoint
    pos, tangent = evaluate_path3d_at_t(spine, t)

    # Normalize tangent (should already be normalized, but just in case)
    mag = math.sqrt(tangent[0]**2 + tangent[1]**2 + tangent[2]**2)
    if mag > 0:
        tangent = [tangent[0]/mag, tangent[1]/mag, tangent[2]/mag]

    return pos, tangent


def add_terminal_connectors_to_segment(
    segment_solid: Any,
    provenance: 'SweptElementProvenance',
    *,
    add_start: bool = False,
    add_end: bool = False,
    connector_length: Optional[float] = None,
    fit_clearance: float = FIT_CLEARANCE['press'],
) -> Any:
    """Union terminal connector tabs with a segment's endpoints.

    Use this to add male tabs to the ends of arc segments that will
    slot into female holes in rings or other structures.

    Args:
        segment_solid: The segment solid to modify
        provenance: SweptElementProvenance with profile and spine data
        add_start: Add terminal tab at start (t=0) of spine
        add_end: Add terminal tab at end (t=1) of spine
        connector_length: Length of terminal tabs (auto-computed if None)
        fit_clearance: Clearance for fit

    Returns:
        Modified segment solid with terminal tabs unioned

    Raises:
        ValueError: If provenance lacks required data
        RuntimeError: If boolean union fails
    """
    from yapcad.geom3d import solid_boolean
    from yapcad.geom import geomlistbbox

    if not add_start and not add_end:
        return segment_solid  # Nothing to do

    # Get profile dimensions
    outer_profile = provenance.outer_profile
    if not outer_profile or not isinstance(outer_profile, list):
        raise ValueError("Provenance must have valid outer_profile")

    outer = outer_profile[0]  # First element is outer boundary
    bbox = geomlistbbox(outer)
    if not bbox or len(bbox) != 2:
        raise ValueError("Cannot extract dimensions from outer_profile")

    outer_w = bbox[1][0] - bbox[0][0]
    outer_h = bbox[1][1] - bbox[0][1]

    wall_thickness = provenance.wall_thickness
    if wall_thickness is None:
        raise ValueError("Wall thickness required for terminal connectors")

    result = segment_solid

    if add_start:
        start_tab = create_terminal_connector(
            outer_w, outer_h,
            provenance.spine,
            "start",
            wall_thickness=wall_thickness,
            connector_length=connector_length,
            fit_clearance=fit_clearance,
        )
        try:
            result = solid_boolean(result, start_tab, 'union')
        except Exception as e:
            raise RuntimeError(f"Terminal connector union (start) failed: {e}") from e

    if add_end:
        end_tab = create_terminal_connector(
            outer_w, outer_h,
            provenance.spine,
            "end",
            wall_thickness=wall_thickness,
            connector_length=connector_length,
            fit_clearance=fit_clearance,
        )
        try:
            result = solid_boolean(result, end_tab, 'union')
        except Exception as e:
            raise RuntimeError(f"Terminal connector union (end) failed: {e}") from e

    return result
