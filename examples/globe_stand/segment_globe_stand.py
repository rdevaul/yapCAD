#!/usr/bin/env python3
"""Globe Stand Segmentation Example

This script demonstrates the yapCAD manufacturing post-processing module by
segmenting the support arcs from the Mars globe stand design into printable
pieces with interior connectors for reassembly.

The globe stand has three support arcs that span ~300-450mm in arc length.
For 3D printers with limited build volume (e.g., Bambu X1C: 256mm), these
arcs need to be segmented into smaller pieces. The manufacturing module:

1. Analyzes swept element geometry and provenance
2. Computes optimal cut locations based on max segment length
3. Splits solids at cut planes using OCC boolean operations
4. Generates interior connectors that follow the original spine curve
5. Exports all geometry as STEP files for CAM/slicing software

Usage:
    python segment_globe_stand.py [OPTIONS] [max_segment_length]

Options:
    --no-step       Skip STEP file export
    --package       Create yapCAD package with assembled visualization
    --help, -h      Show help message

Arguments:
    max_segment_length: Maximum length of each segment in mm (default: 200)

Output:
    Creates output/ directory with:
    - Individual STEP files for each segment and connector
    - assembly_manifest.json with part relationships and file paths
    - Assembly instructions printed to console
    - Optional: assembled.ycpkg package for visualization

Example:
    python segment_globe_stand.py 150        # Segment at 150mm max length
    python segment_globe_stand.py --no-step  # Skip STEP export for testing
    python segment_globe_stand.py --package  # Create assembled visualization package
"""

import json
import math
import os
import sys
from pathlib import Path

# Add src to path for development
sys.path.insert(0, str(Path(__file__).parent.parent.parent / 'src'))

from yapcad.geom import line, point
from yapcad.geom3d_util import sweep_adaptive
from yapcad.manufacturing import (
    SweptElementProvenance,
    CutPoint,
    compute_optimal_cuts,
    segment_swept_element,
    segment_closed_ring,
    build_connector_solids,
    path_length,
    FIT_CLEARANCE,
    # Terminal connector functions
    create_terminal_connector,
    add_terminal_connectors_to_segment,
    # Ring functions
    create_ring_solid,
    compute_arc_attachment_point,
    add_female_holes_to_ring,
    trim_segment_against_ring,
    compute_ring_cuts_avoiding_holes,
)
from yapcad.geom3d import solid_boolean
from yapcad.brep import brep_from_solid
from yapcad.package import create_package_from_entities
from yapcad.metadata import get_solid_metadata


def export_solid_to_step(solid, filepath: str) -> bool:
    """Export a yapCAD solid to a STEP file.

    Args:
        solid: yapCAD solid with BREP data
        filepath: Output STEP file path

    Returns:
        True if export succeeded, False otherwise
    """
    try:
        from OCC.Core.STEPControl import STEPControl_Writer, STEPControl_AsIs
        from OCC.Core.IFSelect import IFSelect_RetDone

        brep = brep_from_solid(solid)
        if brep is None or brep.shape is None:
            print(f"    Warning: No BREP data for export to {filepath}")
            return False

        # Create STEP writer
        writer = STEPControl_Writer()
        writer.Transfer(brep.shape, STEPControl_AsIs)

        # Write file
        status = writer.Write(filepath)
        if status == IFSelect_RetDone:
            return True
        else:
            print(f"    Warning: STEP write failed for {filepath}")
            return False

    except ImportError as e:
        print(f"    Warning: OCC not available for STEP export: {e}")
        return False
    except Exception as e:
        print(f"    Warning: STEP export failed: {e}")
        return False


def create_assembled_package(
    all_solids: list,
    output_path: Path,
    max_segment_length: float,
) -> bool:
    """Create a yapCAD package with all segments and connectors in assembled position.

    Args:
        all_solids: List of (id, solid, part_type) tuples
        output_path: Directory containing output files
        max_segment_length: Max segment length used for segmentation

    Returns:
        True if package was created successfully
    """
    if not all_solids:
        print("  No solids available for package creation")
        return False

    # Collect solids and tag them with metadata
    entities = []
    for part_id, solid, part_type in all_solids:
        if solid is not None:
            # Update metadata for the solid (preserve existing BREP data)
            meta = get_solid_metadata(solid, create=True)
            meta['name'] = part_id
            meta['tags'] = [part_type]
            entities.append(solid)

    if not entities:
        print("  No valid solids for package creation")
        return False

    package_path = output_path / "assembled.ycpkg"

    try:
        manifest = create_package_from_entities(
            entities,
            package_path,
            name="Mars Globe Stand - Segmented Assembly",
            version="1.0.0",
            description=f"Segmented globe stand support arcs with {max_segment_length}mm max segment length. "
                       f"Contains {len(entities)} parts (segments + connectors) in assembled position.",
            units="mm",
            overwrite=True,
        )
        print(f"\n  Created assembly package: {package_path}")
        print(f"    Parts included: {len(entities)}")
        return True
    except Exception as e:
        print(f"  Warning: Package creation failed: {e}")
        return False


# ============================================================================
# Globe Stand Parameters (from globe_stand_v5.dsl)
# ============================================================================

GLOBE_DIAMETER = 304.8      # 12 inches in mm
BASE_DIAMETER = 400.0
TILT_ANGLE = 25.2           # degrees
CRADLE_LATITUDE = -35.0     # degrees (Southern hemisphere)
WAIST_RATIO = 0.7
TWIST_DEG = 60.0
BEAM_OUTER = 10.0           # mm
BEAM_WALL = 2.0             # mm
MARS_OBLATENESS = 0.00648

# Ring parameters
BASE_RING_RADIUS = BASE_DIAMETER / 2.0  # Same as arc base radius
CRADLE_RING_RADIUS = None  # Computed based on cradle latitude and globe size


# ============================================================================
# Profile Generation
# ============================================================================

def box_beam_profile(outer_size: float, wall_thickness: float) -> list:
    """Create outer profile for hollow box beam (region2d format).

    Returns a region2d as list containing one polyline (the outer boundary).
    """
    half = outer_size / 2.0
    outline = [
        line(point(-half, -half, 0), point(half, -half, 0)),
        line(point(half, -half, 0), point(half, half, 0)),
        line(point(half, half, 0), point(-half, half, 0)),
        line(point(-half, half, 0), point(-half, -half, 0)),
    ]
    return [outline]


def box_beam_inner(outer_size: float, wall_thickness: float) -> list:
    """Create inner void profile for hollow box beam (region2d format)."""
    inner_half = outer_size / 2.0 - wall_thickness
    outline = [
        line(point(-inner_half, -inner_half, 0), point(inner_half, -inner_half, 0)),
        line(point(inner_half, -inner_half, 0), point(inner_half, inner_half, 0)),
        line(point(inner_half, inner_half, 0), point(-inner_half, inner_half, 0)),
        line(point(-inner_half, inner_half, 0), point(-inner_half, -inner_half, 0)),
    ]
    return [outline]


# ============================================================================
# Support Arc Path Generation
# ============================================================================

def deg2rad(degrees: float) -> float:
    return degrees * math.pi / 180.0


def support_arc_point(
    base_x: float, base_y: float, base_z: float,
    top_x: float, top_y: float, top_z: float,
    waist_ratio: float,
    a_clear: float, c_clear: float,
    cos_tilt: float, sin_tilt: float,
    globe_center_z: float,
    sigmoid_k: float,
    t: float
) -> tuple:
    """Compute a single point on the support arc with sigmoid blending.

    Returns (x, y, z) tuple.
    """
    # Parabolic factor peaks at t=0.5
    para = 4.0 * t * (1.0 - t)

    # Linear interpolation
    lin_x = base_x + (top_x - base_x) * t
    lin_y = base_y + (top_y - base_y) * t
    lin_z = base_z + (top_z - base_z) * t

    # Apply waist (radial contraction)
    lin_r = math.sqrt(lin_x * lin_x + lin_y * lin_y)
    if lin_r < 1e-10:
        lin_r = 1e-10
    waist_offset = (1.0 - waist_ratio) * lin_r * para
    scale = (lin_r - waist_offset) / lin_r
    raw_x = lin_x * scale
    raw_y = lin_y * scale
    raw_z = lin_z

    # Transform to local (untilted) frame to check against clearance ellipsoid
    loc_x = raw_x
    loc_y = raw_y * cos_tilt + (raw_z - globe_center_z) * sin_tilt
    loc_z = -raw_y * sin_tilt + (raw_z - globe_center_z) * cos_tilt

    # Check if inside clearance ellipsoid
    a_clear2 = a_clear * a_clear
    c_clear2 = c_clear * c_clear
    check = (loc_x*loc_x + loc_y*loc_y) / a_clear2 + loc_z*loc_z / c_clear2

    if check < 1e-10:
        check = 1e-10

    # Sigmoid blend factor: 1 when inside (check<1), 0 when outside
    d = check - 1.0
    s = 1.0 / (1.0 + math.exp(sigmoid_k * d))

    # Project to clearance surface
    proj_t = 1.0 / math.sqrt(check)
    proj_loc_x = loc_x * proj_t
    proj_loc_y = loc_y * proj_t
    proj_loc_z = loc_z * proj_t

    # Transform projected point back to world frame
    proj_world_x = proj_loc_x
    proj_world_y = proj_loc_y * cos_tilt - proj_loc_z * sin_tilt
    proj_world_z = proj_loc_y * sin_tilt + proj_loc_z * cos_tilt + globe_center_z

    # Blend between raw and projected based on sigmoid
    x = (1.0 - s) * raw_x + s * proj_world_x
    y = (1.0 - s) * raw_y + s * proj_world_y
    z = (1.0 - s) * raw_z + s * proj_world_z

    return (x, y, z)


def wrapped_support_arc_path(
    base_radius: float,
    globe_diameter: float,
    oblateness: float,
    latitude_deg: float,
    tilt_deg: float,
    globe_center_z: float,
    beam_offset: float,
    waist_ratio: float,
    base_angle_deg: float,
    top_longitude_deg: float,
    num_segments: int = 17
) -> dict:
    """Generate the 3D path for a support arc.

    Returns path3d dict compatible with manufacturing module.
    """
    globe_radius = globe_diameter / 2.0
    polar_radius = globe_radius * (1.0 - oblateness)

    # Clearance ellipsoid = globe + beam_offset
    a_clear = globe_radius + beam_offset
    c_clear = polar_radius + beam_offset

    sigmoid_k = 20.0

    lat_rad = deg2rad(latitude_deg)
    tilt_rad = deg2rad(tilt_deg)
    base_angle_rad = deg2rad(base_angle_deg)
    top_lon_rad = deg2rad(top_longitude_deg)

    cos_lat = math.cos(lat_rad)
    sin_lat = math.sin(lat_rad)
    cos_tilt = math.cos(tilt_rad)
    sin_tilt = math.sin(tilt_rad)

    a2 = globe_radius * globe_radius
    c2 = polar_radius * polar_radius

    # Base point (on flat ring at z=0)
    base_x = base_radius * math.cos(base_angle_rad)
    base_y = base_radius * math.sin(base_angle_rad)
    base_z = 0.0

    # Top point: on latitude circle with beam offset
    tx = globe_radius * cos_lat * math.cos(top_lon_rad)
    ty = globe_radius * cos_lat * math.sin(top_lon_rad)
    tz = polar_radius * sin_lat
    tnx = tx / a2
    tny = ty / a2
    tnz = tz / c2
    tnlen = math.sqrt(tnx*tnx + tny*tny + tnz*tnz)
    if tnlen < 1e-10:
        tnlen = 1e-10
    tpx = tx + beam_offset * tnx / tnlen
    tpy = ty + beam_offset * tny / tnlen
    tpz = tz + beam_offset * tnz / tnlen
    top_x = tpx
    top_y = tpy * cos_tilt - tpz * sin_tilt
    top_z = tpy * sin_tilt + tpz * cos_tilt + globe_center_z

    # Generate intermediate points
    # t values: 1/(num_segments), 2/(num_segments), ..., (num_segments-1)/(num_segments)
    arc_points = []
    for i in range(1, num_segments):
        t = i / num_segments
        pt = support_arc_point(
            base_x, base_y, base_z, top_x, top_y, top_z,
            waist_ratio, a_clear, c_clear,
            cos_tilt, sin_tilt, globe_center_z, sigmoid_k, t
        )
        arc_points.append(pt)

    # Build path segments
    segments = []

    # First segment: base -> first intermediate
    segments.append({
        'type': 'line',
        'start': [base_x, base_y, base_z],
        'end': list(arc_points[0])
    })

    # Middle segments
    for i in range(len(arc_points) - 1):
        segments.append({
            'type': 'line',
            'start': list(arc_points[i]),
            'end': list(arc_points[i + 1])
        })

    # Last segment: last intermediate -> top
    segments.append({
        'type': 'line',
        'start': list(arc_points[-1]),
        'end': [top_x, top_y, top_z]
    })

    return {'segments': segments}


def compute_globe_center_z(
    globe_diameter: float,
    oblateness: float,
    latitude_deg: float,
    tilt_deg: float,
    beam_offset: float,
    target_cradle_low_z: float = 300.0
) -> float:
    """Calculate globe center Z so lowest cradle point is at target height."""
    globe_radius = globe_diameter / 2.0
    polar_radius = globe_radius * (1.0 - oblateness)

    lat_rad = deg2rad(latitude_deg)
    tilt_rad = deg2rad(tilt_deg)

    cos_lat = math.cos(lat_rad)
    sin_lat = math.sin(lat_rad)
    cos_tilt = math.cos(tilt_rad)
    sin_tilt = math.sin(tilt_rad)

    a2 = globe_radius * globe_radius
    c2 = polar_radius * polar_radius

    # Point at 270° longitude (negative Y) is the lowest
    x270 = 0.0
    y270 = -globe_radius * cos_lat
    z270 = polar_radius * sin_lat
    nx270 = x270 / a2
    ny270 = y270 / a2
    nz270 = z270 / c2
    nlen270 = math.sqrt(nx270*nx270 + ny270*ny270 + nz270*nz270)
    if nlen270 < 1e-10:
        nlen270 = 1e-10
    px270 = x270 + beam_offset * nx270 / nlen270
    py270 = y270 + beam_offset * ny270 / nlen270
    pz270 = z270 + beam_offset * nz270 / nlen270
    rz270 = py270 * sin_tilt + pz270 * cos_tilt

    return target_cradle_low_z - rz270


def cradle_point(
    globe_radius: float,
    polar_radius: float,
    cos_lat: float,
    sin_lat: float,
    cos_tilt: float,
    sin_tilt: float,
    beam_offset: float,
    globe_center_z: float,
    longitude_deg: float,
) -> tuple:
    """Compute a single point on the latitude circle (matching DSL cradle_point).

    Returns world-space (x, y, z) after ellipsoid offset and tilt rotation.
    """
    lon_rad = deg2rad(longitude_deg)
    cos_lon = math.cos(lon_rad)
    sin_lon = math.sin(lon_rad)

    # Point on latitude circle (before tilt)
    x = globe_radius * cos_lat * cos_lon
    y = globe_radius * cos_lat * sin_lon
    z = polar_radius * sin_lat

    # Normal vector for offset (gradient of ellipsoid)
    a2 = globe_radius * globe_radius
    c2 = polar_radius * polar_radius
    nx = x / a2
    ny = y / a2
    nz = z / c2
    nlen = math.sqrt(nx*nx + ny*ny + nz*nz)
    if nlen < 1e-10:
        nlen = 1e-10

    # Offset point along normal
    px = x + beam_offset * nx / nlen
    py = y + beam_offset * ny / nlen
    pz = z + beam_offset * nz / nlen

    # Apply tilt rotation around X axis, then translate
    rx = px
    ry = py * cos_tilt - pz * sin_tilt
    rz = py * sin_tilt + pz * cos_tilt + globe_center_z

    return (rx, ry, rz)


def latitude_cradle_path(
    globe_diameter: float,
    oblateness: float,
    latitude_deg: float,
    tilt_deg: float,
    globe_center_z: float,
    beam_offset: float,
    num_segments: int = 32,
) -> dict:
    """Generate the 3D path for a latitude cradle ring (matching DSL).

    Returns path3d dict compatible with manufacturing module.
    """
    globe_radius = globe_diameter / 2.0
    polar_radius = globe_radius * (1.0 - oblateness)

    lat_rad = deg2rad(latitude_deg)
    tilt_rad = deg2rad(tilt_deg)

    cos_lat = math.cos(lat_rad)
    sin_lat = math.sin(lat_rad)
    cos_tilt = math.cos(tilt_rad)
    sin_tilt = math.sin(tilt_rad)

    # Generate points around the latitude circle
    angle_step = 360.0 / num_segments
    points = []
    for i in range(num_segments):
        longitude = i * angle_step
        pt = cradle_point(
            globe_radius, polar_radius, cos_lat, sin_lat,
            cos_tilt, sin_tilt, beam_offset, globe_center_z,
            longitude
        )
        points.append(pt)

    # Build path segments connecting consecutive points (closing the loop)
    segments = []
    for i in range(num_segments):
        next_i = (i + 1) % num_segments
        segments.append({
            'type': 'line',
            'start': list(points[i]),
            'end': list(points[next_i])
        })

    return {'segments': segments}


def create_latitude_cradle_ring(
    globe_diameter: float,
    oblateness: float,
    latitude_deg: float,
    tilt_deg: float,
    globe_center_z: float,
    beam_outer: float,
    beam_wall: float,
    ring_id: str = "cradle_ring",
    num_path_segments: int = 32,
) -> tuple:
    """Create a cradle ring solid following the latitude circle on tilted oblate spheroid.

    Returns:
        Tuple of (solid, SweptElementProvenance)
    """
    from yapcad.geom3d_util import sweep_adaptive
    from yapcad.manufacturing import SweptElementProvenance, path_length

    beam_offset = beam_outer / 2.0

    # Generate the latitude path
    spine = latitude_cradle_path(
        globe_diameter, oblateness, latitude_deg, tilt_deg,
        globe_center_z, beam_offset, num_path_segments
    )

    # Create hollow box profile
    outer_profile = box_beam_profile(beam_outer, beam_wall)
    inner_profile = box_beam_inner(beam_outer, beam_wall)

    # Sweep the profile along the path
    ring_solid = sweep_adaptive(
        outer_profile[0],
        spine,
        angle_threshold_deg=5.0,
        inner_profiles=[inner_profile[0]],
    )

    # Create provenance
    provenance = SweptElementProvenance(
        id=ring_id,
        operation="sweep_adaptive",
        outer_profile=outer_profile,
        spine=spine,
        wall_thickness=beam_wall,
        metadata={'solid': ring_solid}
    )

    return ring_solid, provenance


def compute_cradle_ring_params(
    globe_diameter: float,
    oblateness: float,
    latitude_deg: float,
    tilt_deg: float,
    beam_offset: float,
    globe_center_z: float,
) -> tuple:
    """Compute cradle ring approximate radius and center Z position.

    NOTE: This returns approximate values for display purposes only.
    The actual cradle ring is a 3D path, not a simple tilted circle.

    Returns:
        Tuple of (approximate_radius, approximate_center_z)
    """
    globe_radius = globe_diameter / 2.0
    polar_radius = globe_radius * (1.0 - oblateness)

    lat_rad = deg2rad(latitude_deg)
    tilt_rad = deg2rad(tilt_deg)

    cos_lat = math.cos(lat_rad)
    sin_lat = math.sin(lat_rad)
    cos_tilt = math.cos(tilt_rad)
    sin_tilt = math.sin(tilt_rad)

    # Compute a representative point on the latitude circle (at 0° longitude)
    pt0 = cradle_point(
        globe_radius, polar_radius, cos_lat, sin_lat,
        cos_tilt, sin_tilt, beam_offset, globe_center_z,
        0.0
    )

    # Approximate radius from the XY plane distance
    approx_radius = math.sqrt(pt0[0]**2 + pt0[1]**2)

    # Use average Z from 0° and 180° longitude points
    pt180 = cradle_point(
        globe_radius, polar_radius, cos_lat, sin_lat,
        cos_tilt, sin_tilt, beam_offset, globe_center_z,
        180.0
    )
    approx_z = (pt0[2] + pt180[2]) / 2.0

    return approx_radius, approx_z


def compute_arc_endpoint_for_ring(
    base_radius: float,
    base_angle_deg: float,
    top_longitude_deg: float,
    globe_diameter: float,
    oblateness: float,
    latitude_deg: float,
    tilt_deg: float,
    beam_offset: float,
    globe_center_z: float,
    end: str,  # "start" or "end"
) -> tuple:
    """Compute the position and direction for arc terminal connector.

    Args:
        end: "start" for base ring attachment, "end" for cradle ring

    Returns:
        Tuple of (position, inward_direction)
    """
    globe_radius = globe_diameter / 2.0
    polar_radius = globe_radius * (1.0 - oblateness)

    lat_rad = deg2rad(latitude_deg)
    tilt_rad = deg2rad(tilt_deg)
    base_angle_rad = deg2rad(base_angle_deg)
    top_lon_rad = deg2rad(top_longitude_deg)

    cos_lat = math.cos(lat_rad)
    sin_lat = math.sin(lat_rad)
    cos_tilt = math.cos(tilt_rad)
    sin_tilt = math.sin(tilt_rad)

    if end == "start":
        # Base attachment point - on flat ring at z=0
        x = base_radius * math.cos(base_angle_rad)
        y = base_radius * math.sin(base_angle_rad)
        z = 0.0
        # Direction: pointing into the ring (toward center)
        dx = -math.cos(base_angle_rad)
        dy = -math.sin(base_angle_rad)
        dz = 0.0
    else:
        # Cradle attachment point - on tilted latitude circle
        a2 = globe_radius * globe_radius
        c2 = polar_radius * polar_radius

        tx = globe_radius * cos_lat * math.cos(top_lon_rad)
        ty = globe_radius * cos_lat * math.sin(top_lon_rad)
        tz = polar_radius * sin_lat

        # Normal at this point
        tnx = tx / a2
        tny = ty / a2
        tnz = tz / c2
        tnlen = math.sqrt(tnx*tnx + tny*tny + tnz*tnz)
        if tnlen < 1e-10:
            tnlen = 1e-10

        # Offset by beam_offset along normal
        x = tx + beam_offset * tnx / tnlen
        y = ty + beam_offset * tny / tnlen
        z = tz + beam_offset * tnz / tnlen

        # Transform to world coords (apply tilt)
        world_x = x
        world_y = y * cos_tilt - z * sin_tilt
        world_z = y * sin_tilt + z * cos_tilt + globe_center_z

        # Direction: tangent to cradle circle, pointing inward
        # Simplified: radial inward in the tilted cradle plane
        dx = -math.cos(top_lon_rad)
        dy = -math.sin(top_lon_rad)
        dz = 0.0
        # Transform direction
        world_dx = dx
        world_dy = dy * cos_tilt
        world_dz = dy * sin_tilt

        return (world_x, world_y, world_z), (world_dx, world_dy, world_dz)

    return (x, y, z), (dx, dy, dz)


# ============================================================================
# Main Segmentation Logic
# ============================================================================

def create_support_arc_provenance(
    arc_index: int,
    base_angle_deg: float,
    top_longitude_deg: float,
    globe_center_z: float
) -> SweptElementProvenance:
    """Create a support arc with full provenance tracking."""
    base_radius = BASE_DIAMETER / 2.0
    beam_offset = BEAM_OUTER / 2.0

    # Create profiles
    outer_profile = box_beam_profile(BEAM_OUTER, BEAM_WALL)
    inner_profile = box_beam_inner(BEAM_OUTER, BEAM_WALL)

    # Create path
    spine = wrapped_support_arc_path(
        base_radius, GLOBE_DIAMETER, MARS_OBLATENESS, CRADLE_LATITUDE,
        TILT_ANGLE, globe_center_z, beam_offset,
        WAIST_RATIO, base_angle_deg, top_longitude_deg
    )

    # Build swept solid
    print(f"  Building support arc {arc_index} solid...")
    try:
        arc_solid = sweep_adaptive(
            outer_profile[0],  # outer boundary polyline
            spine,
            inner_profiles=[inner_profile[0]],  # inner void polyline
            angle_threshold_deg=5.0
        )
    except Exception as e:
        print(f"  Warning: Could not build solid for arc {arc_index}: {e}")
        arc_solid = None

    # Create provenance
    return SweptElementProvenance(
        id=f"support_arc_{arc_index}",
        operation="sweep_adaptive_hollow",
        outer_profile=outer_profile,
        spine=spine,
        inner_profile=inner_profile,
        wall_thickness=BEAM_WALL,
        metadata={'solid': arc_solid}
    )


def segment_globe_stand(max_segment_length: float = 200.0, output_dir: str = "output",
                        export_step: bool = True, create_package: bool = False,
                        include_rings: bool = False):
    """Segment the globe stand support arcs and generate outputs.

    Args:
        max_segment_length: Maximum length of each segment in mm (default: 200)
        output_dir: Output directory name (default: "output")
        export_step: Whether to export STEP files (default: True)
        create_package: Whether to create a yapCAD package for visualization (default: False)
        include_rings: Whether to generate base and cradle rings (default: False)
    """

    print("=" * 60)
    print("Globe Stand Segmentation Example")
    print("=" * 60)
    print()
    print(f"Parameters:")
    print(f"  Globe diameter: {GLOBE_DIAMETER} mm")
    print(f"  Base diameter: {BASE_DIAMETER} mm")
    print(f"  Beam profile: {BEAM_OUTER}x{BEAM_OUTER} mm, wall: {BEAM_WALL} mm")
    print(f"  Max segment length: {max_segment_length} mm")
    print(f"  STEP export: {export_step}")
    print(f"  Create package: {create_package}")
    print(f"  Include rings: {include_rings}")
    print()

    # Calculate globe center Z
    beam_offset = BEAM_OUTER / 2.0
    globe_center_z = compute_globe_center_z(
        GLOBE_DIAMETER, MARS_OBLATENESS, CRADLE_LATITUDE,
        TILT_ANGLE, beam_offset
    )
    print(f"  Globe center Z: {globe_center_z:.1f} mm")
    print()

    # Create output directory
    output_path = Path(__file__).parent / output_dir
    output_path.mkdir(exist_ok=True)

    # Define the three support arcs (at 0°, 120°, 240°)
    arc_configs = [
        (0, 0.0, TWIST_DEG),
        (1, 120.0, 120.0 + TWIST_DEG),
        (2, 240.0, 240.0 + TWIST_DEG),
    ]

    all_results = []
    all_solids = []  # Collect (id, solid, type) tuples for package creation
    terminal_segments = {}  # Track terminal segments for ring integration: {arc_id: (first_seg, last_seg, prov)}
    arc_solids = {}  # Store full arc solids before segmentation: {arc_idx: solid}
    assembly_manifest = {
        'project': 'Mars Globe Stand',
        'max_segment_length': max_segment_length,
        'beam_profile': f'{BEAM_OUTER}x{BEAM_OUTER}mm hollow, {BEAM_WALL}mm wall',
        'fit_clearance': FIT_CLEARANCE['press'],
        'parts': []
    }

    print("Building and segmenting support arcs...")
    print("-" * 40)

    for arc_idx, base_angle, top_lon in arc_configs:
        print(f"\nArc {arc_idx}: base={base_angle}°, top={top_lon}°")

        # Create provenance
        prov = create_support_arc_provenance(arc_idx, base_angle, top_lon, globe_center_z)

        # Calculate path length
        arc_length = path_length(prov.spine)
        print(f"  Path length: {arc_length:.1f} mm")

        # Compute optimal cuts
        cuts = compute_optimal_cuts(prov, max_segment_length)
        num_cuts = len(cuts)
        print(f"  Cuts needed: {num_cuts}")

        if cuts:
            for i, cut in enumerate(cuts):
                print(f"    Cut {i}: t={cut.parameter:.3f}")

        # Skip segmentation if no solid was built
        if prov.metadata.get('solid') is None:
            print(f"  Skipping segmentation (no solid)")
            assembly_manifest['parts'].append({
                'id': prov.id,
                'type': 'support_arc',
                'status': 'not_segmented',
                'reason': 'solid_build_failed',
                'path_length_mm': arc_length
            })
            continue

        # Store full arc solid for cutting ring holes (before segmentation)
        arc_solids[arc_idx] = prov.metadata['solid']

        # Perform segmentation
        if cuts:
            try:
                result = segment_swept_element(prov, cuts)
                print(f"  Created {result.segment_count} segments")

                # Record and export segments
                for seg in result.segments:
                    step_file = None
                    if export_step and seg.solid is not None:
                        step_filename = f"{seg.id}.step"
                        step_filepath = str(output_path / step_filename)
                        if export_solid_to_step(seg.solid, step_filepath):
                            step_file = step_filename
                            print(f"    Exported: {step_filename}")

                    # Collect solid for package creation
                    if seg.solid is not None:
                        all_solids.append((seg.id, seg.solid, 'segment'))

                    part_info = {
                        'id': seg.id,
                        'type': 'arc_segment',
                        'parent': prov.id,
                        'parameter_range': list(seg.parameter_range),
                        'mates_with': seg.mates_with,
                        'has_connector_tab': seg.has_connector_tab,
                        'connector_type': seg.connector_type,
                        'step_file': step_file
                    }
                    assembly_manifest['parts'].append(part_info)

                # Build and export connectors
                # Build connector solids if we need them for STEP export or package
                if (export_step or create_package) and result.connectors:
                    print(f"  Building connector solids...")
                    try:
                        connector_solids = build_connector_solids(prov, result.connectors)
                        for conn_spec, conn_solid in connector_solids:
                            # Collect solid for package creation
                            all_solids.append((conn_spec.id, conn_solid, 'connector'))

                            step_filename = None
                            if export_step:
                                step_filename = f"{conn_spec.id}.step"
                                step_filepath = str(output_path / step_filename)
                                if export_solid_to_step(conn_solid, step_filepath):
                                    print(f"    Exported: {step_filename}")
                                else:
                                    step_filename = None

                            conn_info = {
                                'id': conn_spec.id,
                                'type': 'interior_connector',
                                'parent': prov.id,
                                'center_parameter': conn_spec.center_parameter,
                                'length_mm': conn_spec.length,
                                'fit_clearance_mm': conn_spec.fit_clearance,
                                'step_file': step_filename
                            }
                            assembly_manifest['parts'].append(conn_info)
                    except Exception as e:
                        print(f"  Warning: Connector build failed: {e}")
                        # Still record connector specs without solid export
                        for conn in result.connectors:
                            conn_info = {
                                'id': conn.id,
                                'type': 'interior_connector',
                                'parent': prov.id,
                                'center_parameter': conn.center_parameter,
                                'length_mm': conn.length,
                                'fit_clearance_mm': conn.fit_clearance,
                                'step_file': None
                            }
                            assembly_manifest['parts'].append(conn_info)
                elif result.connectors:
                    # Record connectors without building solids
                    for conn in result.connectors:
                        conn_info = {
                            'id': conn.id,
                            'type': 'interior_connector',
                            'parent': prov.id,
                            'center_parameter': conn.center_parameter,
                            'length_mm': conn.length,
                            'fit_clearance_mm': conn.fit_clearance,
                            'step_file': None
                        }
                        assembly_manifest['parts'].append(conn_info)

                all_results.append((prov, result))

                # Track terminal segments for ring integration
                if result.segments:
                    first_seg = result.segments[0]
                    last_seg = result.segments[-1]
                    terminal_segments[prov.id] = (first_seg, last_seg, prov)

                # Print assembly instructions
                if result.assembly_instructions:
                    print(f"\n  Assembly instructions for arc {arc_idx}:")
                    for line_text in result.assembly_instructions.split('\n')[:10]:
                        print(f"    {line_text}")

            except Exception as e:
                print(f"  Segmentation failed: {e}")
                assembly_manifest['parts'].append({
                    'id': prov.id,
                    'type': 'support_arc',
                    'status': 'segmentation_failed',
                    'error': str(e)
                })
        else:
            # No cuts needed, record as single piece
            assembly_manifest['parts'].append({
                'id': prov.id,
                'type': 'support_arc',
                'status': 'no_segmentation_needed',
                'path_length_mm': arc_length
            })

    # ========================================================================
    # Ring Generation (if --rings flag is set)
    # ========================================================================
    if include_rings:
        print("\n" + "-" * 40)
        print("Generating rings with female holes...")
        print("-" * 40)

        # Compute cradle ring parameters
        cradle_radius, cradle_z = compute_cradle_ring_params(
            GLOBE_DIAMETER, MARS_OBLATENESS, CRADLE_LATITUDE,
            TILT_ANGLE, beam_offset, globe_center_z
        )
        print(f"  Cradle ring: radius={cradle_radius:.1f}mm, z={cradle_z:.1f}mm")

        # Base ring attachment angles (same as arc base angles)
        base_attachment_angles = [0.0, 120.0, 240.0]
        # Cradle ring attachment angles (twisted by TWIST_DEG)
        cradle_attachment_angles = [TWIST_DEG, 120.0 + TWIST_DEG, 240.0 + TWIST_DEG]

        # Create base ring
        print("\n  Creating base ring...")
        base_ring_solid, base_ring_prov = create_ring_solid(
            radius=BASE_RING_RADIUS,
            outer_width=BEAM_OUTER,
            outer_height=BEAM_OUTER,
            wall_thickness=BEAM_WALL,
            center=(0.0, 0.0, 0.0),
            tilt_angle_deg=0.0,
            tilt_axis="x",
            ring_id="base_ring"
        )

        # Cut female holes using arc solids (properly aligned to arc direction)
        print("    Cutting holes using arc geometry...")
        base_ring_with_holes = base_ring_solid
        for arc_idx, arc_solid in arc_solids.items():
            try:
                base_ring_with_holes = solid_boolean(
                    base_ring_with_holes, arc_solid, 'difference'
                )
                print(f"      Cut hole for arc {arc_idx} in base ring")
            except Exception as e:
                print(f"      Warning: Failed to cut hole for arc {arc_idx}: {e}")
        base_ring_prov.metadata['solid'] = base_ring_with_holes

        # Create cradle ring following the actual latitude path on tilted oblate spheroid
        print("\n  Creating cradle ring (latitude path)...")
        cradle_ring_solid, cradle_ring_prov = create_latitude_cradle_ring(
            globe_diameter=GLOBE_DIAMETER,
            oblateness=MARS_OBLATENESS,
            latitude_deg=CRADLE_LATITUDE,
            tilt_deg=TILT_ANGLE,
            globe_center_z=globe_center_z,
            beam_outer=BEAM_OUTER,
            beam_wall=BEAM_WALL,
            ring_id="cradle_ring"
        )
        print(f"    Cradle ring: approx radius={cradle_radius:.1f}mm, approx z={cradle_z:.1f}mm")

        # Cut female holes using arc solids (properly aligned to arc direction)
        print("    Cutting holes using arc geometry...")
        cradle_ring_with_holes = cradle_ring_solid
        for arc_idx, arc_solid in arc_solids.items():
            try:
                cradle_ring_with_holes = solid_boolean(
                    cradle_ring_with_holes, arc_solid, 'difference'
                )
                print(f"      Cut hole for arc {arc_idx} in cradle ring")
            except Exception as e:
                print(f"      Warning: Failed to cut hole for arc {arc_idx}: {e}")
        cradle_ring_prov.metadata['solid'] = cradle_ring_with_holes

        # ====================================================================
        # Arc-Ring Integration: Trim terminal segments and add terminal tabs
        # ====================================================================
        if terminal_segments:
            print("\n  Integrating arc segments with rings...")

            for arc_id, (first_seg, last_seg, prov) in terminal_segments.items():
                print(f"    Processing {arc_id}...")

                # Find and update first segment (base ring connection) in all_solids
                if first_seg.solid is not None:
                    # Trim against base ring
                    try:
                        trimmed_first = trim_segment_against_ring(
                            first_seg.solid, base_ring_solid
                        )
                        # Update segment solid reference
                        first_seg.solid = trimmed_first

                        # Update in all_solids list
                        for i, (sid, solid, stype) in enumerate(all_solids):
                            if sid == first_seg.id:
                                all_solids[i] = (sid, trimmed_first, stype)
                                break

                        print(f"      Trimmed {first_seg.id} against base ring")
                    except Exception as e:
                        print(f"      Warning: Failed to trim {first_seg.id}: {e}")

                # Find and update last segment (cradle ring connection) in all_solids
                if last_seg.solid is not None and last_seg.id != first_seg.id:
                    # Trim against cradle ring
                    try:
                        trimmed_last = trim_segment_against_ring(
                            last_seg.solid, cradle_ring_solid
                        )
                        # Update segment solid reference
                        last_seg.solid = trimmed_last

                        # Update in all_solids list
                        for i, (sid, solid, stype) in enumerate(all_solids):
                            if sid == last_seg.id:
                                all_solids[i] = (sid, trimmed_last, stype)
                                break

                        print(f"      Trimmed {last_seg.id} against cradle ring")
                    except Exception as e:
                        print(f"      Warning: Failed to trim {last_seg.id}: {e}")

                # Now add terminal tabs to trimmed segments
                # First segment gets a start tab (connects to base ring)
                if first_seg.solid is not None:
                    try:
                        seg_with_tab = add_terminal_connectors_to_segment(
                            first_seg.solid, prov,
                            add_start=True, add_end=False,
                            fit_clearance=FIT_CLEARANCE['press']
                        )
                        first_seg.solid = seg_with_tab

                        # Update in all_solids list
                        for i, (sid, solid, stype) in enumerate(all_solids):
                            if sid == first_seg.id:
                                all_solids[i] = (sid, seg_with_tab, stype)
                                break

                        print(f"      Added terminal tab to {first_seg.id} (start)")
                    except Exception as e:
                        print(f"      Warning: Failed to add tab to {first_seg.id}: {e}")

                # Last segment gets an end tab (connects to cradle ring)
                if last_seg.solid is not None and last_seg.id != first_seg.id:
                    try:
                        seg_with_tab = add_terminal_connectors_to_segment(
                            last_seg.solid, prov,
                            add_start=False, add_end=True,
                            fit_clearance=FIT_CLEARANCE['press']
                        )
                        last_seg.solid = seg_with_tab

                        # Update in all_solids list
                        for i, (sid, solid, stype) in enumerate(all_solids):
                            if sid == last_seg.id:
                                all_solids[i] = (sid, seg_with_tab, stype)
                                break

                        print(f"      Added terminal tab to {last_seg.id} (end)")
                    except Exception as e:
                        print(f"      Warning: Failed to add tab to {last_seg.id}: {e}")

        # Segment the rings
        # Use 3 segments per ring to fit within Bambu Lab H2S build volume
        # when laid flat (ring diameter ~400mm fits within 450x450mm build area)
        RING_SEGMENTS = 3
        print("\n  Segmenting base ring...")
        base_circumference = 2 * math.pi * BASE_RING_RADIUS
        base_cut_params = compute_ring_cuts_avoiding_holes(
            base_circumference, max_segment_length,
            base_attachment_angles, min_distance_from_hole=30.0,
            target_segments=RING_SEGMENTS
        )
        print(f"    Circumference: {base_circumference:.1f}mm, cuts: {len(base_cut_params)}, target segments: {RING_SEGMENTS}")

        if base_cut_params:
            # Use segment_closed_ring for proper closed ring segmentation
            try:
                base_ring_segments = segment_closed_ring(
                    base_ring_with_holes,
                    base_ring_prov.spine,
                    base_cut_params
                )
                print(f"    Created {len(base_ring_segments)} base ring segments")

                # Build segment ranges from cut params
                all_cuts = [0.0] + sorted(base_cut_params) + [1.0]
                for i, seg_solid in enumerate(base_ring_segments):
                    seg_id = f"base_ring_seg_{i}"
                    param_range = (all_cuts[i], all_cuts[i + 1])

                    if seg_solid is not None:
                        all_solids.append((seg_id, seg_solid, 'ring_segment'))

                        if export_step:
                            step_filename = f"{seg_id}.step"
                            step_filepath = str(output_path / step_filename)
                            if export_solid_to_step(seg_solid, step_filepath):
                                print(f"      Exported: {step_filename}")

                        # Determine mating relationships
                        mates = []
                        if i > 0:
                            mates.append(f"base_ring_seg_{i-1}")
                        if i < len(base_ring_segments) - 1:
                            mates.append(f"base_ring_seg_{i+1}")

                        # Connector type: male for all but last
                        connector_type = "male" if i < len(base_ring_segments) - 1 else "female"

                        assembly_manifest['parts'].append({
                            'id': seg_id,
                            'type': 'ring_segment',
                            'parent': 'base_ring',
                            'parameter_range': list(param_range),
                            'mates_with': mates,
                            'connector_type': connector_type,
                        })
            except Exception as e:
                print(f"    Base ring segmentation failed: {e}")
                import traceback
                traceback.print_exc()
                # Fall back to unsegmented ring
                all_solids.append(('base_ring', base_ring_with_holes, 'ring'))
                assembly_manifest['parts'].append({
                    'id': 'base_ring',
                    'type': 'ring',
                    'segmented': False,
                })
        else:
            # No segmentation needed
            all_solids.append(('base_ring', base_ring_with_holes, 'ring'))
            assembly_manifest['parts'].append({
                'id': 'base_ring',
                'type': 'ring',
                'segmented': False,
            })

        print("\n  Segmenting cradle ring...")
        cradle_circumference = 2 * math.pi * cradle_radius
        cradle_cut_params = compute_ring_cuts_avoiding_holes(
            cradle_circumference, max_segment_length,
            cradle_attachment_angles, min_distance_from_hole=30.0,
            target_segments=RING_SEGMENTS
        )
        print(f"    Circumference: {cradle_circumference:.1f}mm, cuts: {len(cradle_cut_params)}, target segments: {RING_SEGMENTS}")

        if cradle_cut_params:
            # Use segment_closed_ring for proper closed ring segmentation
            try:
                cradle_ring_segments = segment_closed_ring(
                    cradle_ring_with_holes,
                    cradle_ring_prov.spine,
                    cradle_cut_params
                )
                print(f"    Created {len(cradle_ring_segments)} cradle ring segments")

                # Build segment ranges from cut params
                all_cuts = [0.0] + sorted(cradle_cut_params) + [1.0]
                for i, seg_solid in enumerate(cradle_ring_segments):
                    seg_id = f"cradle_ring_seg_{i}"
                    param_range = (all_cuts[i], all_cuts[i + 1])

                    if seg_solid is not None:
                        all_solids.append((seg_id, seg_solid, 'ring_segment'))

                        if export_step:
                            step_filename = f"{seg_id}.step"
                            step_filepath = str(output_path / step_filename)
                            if export_solid_to_step(seg_solid, step_filepath):
                                print(f"      Exported: {step_filename}")

                        # Determine mating relationships
                        mates = []
                        if i > 0:
                            mates.append(f"cradle_ring_seg_{i-1}")
                        if i < len(cradle_ring_segments) - 1:
                            mates.append(f"cradle_ring_seg_{i+1}")

                        # Connector type: male for all but last
                        connector_type = "male" if i < len(cradle_ring_segments) - 1 else "female"

                        assembly_manifest['parts'].append({
                            'id': seg_id,
                            'type': 'ring_segment',
                            'parent': 'cradle_ring',
                            'parameter_range': list(param_range),
                            'mates_with': mates,
                            'connector_type': connector_type,
                        })
            except Exception as e:
                print(f"    Cradle ring segmentation failed: {e}")
                import traceback
                traceback.print_exc()
                all_solids.append(('cradle_ring', cradle_ring_with_holes, 'ring'))
                assembly_manifest['parts'].append({
                    'id': 'cradle_ring',
                    'type': 'ring',
                    'segmented': False,
                })
        else:
            all_solids.append(('cradle_ring', cradle_ring_with_holes, 'ring'))
            assembly_manifest['parts'].append({
                'id': 'cradle_ring',
                'type': 'ring',
                'segmented': False,
            })

    # Write assembly manifest
    manifest_path = output_path / "assembly_manifest.json"
    with open(manifest_path, 'w') as f:
        json.dump(assembly_manifest, f, indent=2)
    print(f"\nWrote assembly manifest to: {manifest_path}")

    # Create assembled package if requested
    package_created = False
    if create_package:
        package_created = create_assembled_package(
            all_solids, output_path, max_segment_length
        )

    # Count exported files
    step_files = [p.get('step_file') for p in assembly_manifest['parts'] if p.get('step_file')]

    # Summary
    print()
    print("=" * 60)
    print("Summary")
    print("=" * 60)
    total_segments = sum(r.segment_count for _, r in all_results)
    total_connectors = sum(len(r.connectors) for _, r in all_results)
    print(f"  Total segments: {total_segments}")
    print(f"  Total connectors: {total_connectors}")
    print(f"  STEP files exported: {len(step_files)}")
    if package_created:
        print(f"  Assembly package: assembled.ycpkg")
    print(f"  Output directory: {output_path}")
    print()
    print("Output files:")
    print(f"  - assembly_manifest.json")
    for f in step_files:
        print(f"  - {f}")
    if package_created:
        print(f"  - assembled.ycpkg/")

    return all_results, assembly_manifest


# ============================================================================
# Entry Point
# ============================================================================

def print_usage():
    """Print usage information."""
    print("Usage: python segment_globe_stand.py [OPTIONS] [max_segment_length]")
    print()
    print("Options:")
    print("  --no-step      Skip STEP file export")
    print("  --package      Create yapCAD package with assembled visualization")
    print("  --rings        Include base and cradle rings with female holes")
    print("  --help, -h     Show this help message")
    print()
    print("Arguments:")
    print("  max_segment_length  Maximum segment length in mm (default: 200)")
    print()
    print("Examples:")
    print("  python segment_globe_stand.py              # Default: 200mm segments, STEP export")
    print("  python segment_globe_stand.py 150          # 150mm segments, STEP export")
    print("  python segment_globe_stand.py --no-step    # 200mm segments, no STEP export")
    print("  python segment_globe_stand.py --package    # Create assembled visualization package")
    print("  python segment_globe_stand.py --rings      # Include rings with female holes")


if __name__ == '__main__':
    max_length = 200.0
    export_step = True
    create_package = False
    include_rings = False

    args = sys.argv[1:]

    # Parse arguments
    for arg in args:
        if arg in ('--help', '-h'):
            print_usage()
            sys.exit(0)
        elif arg == '--no-step':
            export_step = False
        elif arg == '--package':
            create_package = True
        elif arg == '--rings':
            include_rings = True
        else:
            try:
                max_length = float(arg)
            except ValueError:
                print(f"Error: Invalid argument: {arg}")
                print_usage()
                sys.exit(1)

    segment_globe_stand(max_length, export_step=export_step, create_package=create_package,
                        include_rings=include_rings)
