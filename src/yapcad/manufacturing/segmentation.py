"""Beam segmentation for manufacturing post-processing.

This module provides functions to segment swept elements (beams, pipes)
at specified cut planes, creating printable segments with interior
connectors for reassembly.
"""

from typing import Any, Dict, List, Optional, Tuple

from .data import (
    CutPoint,
    Segment,
    SegmentationResult,
    ConnectorSpec,
    SweptElementProvenance,
)
from .path_utils import (
    compute_cut_plane,
    extract_sub_path,
    path_length,
)
from .connectors import (
    FIT_CLEARANCE,
    compute_connector_spec,
    create_interior_connector,
)


def _fix_shape(shape):
    """Apply OCC shape fixing to correct orientation and other issues.

    Ensures solids have outward-facing normals (matter inside).

    Args:
        shape: OCC TopoDS_Shape to fix

    Returns:
        Fixed shape with correct orientation
    """
    from OCC.Core.ShapeFix import ShapeFix_Shape, ShapeFix_Solid
    from OCC.Core.TopAbs import TopAbs_SOLID
    from OCC.Core.BRepLib import breplib
    from OCC.Core.TopoDS import topods

    # Apply general shape fixing first
    fixer = ShapeFix_Shape(shape)
    fixer.SetPrecision(1e-6)
    fixer.SetMaxTolerance(1e-4)
    fixer.SetMinTolerance(1e-8)
    fixer.Perform()
    fixed = fixer.Shape()

    # If it's a solid, ensure proper orientation (outward-facing normals)
    if fixed.ShapeType() == TopAbs_SOLID:
        try:
            solid = topods.Solid(fixed)
            # OrientClosedSolid orients the shell so matter is inside
            # (i.e., normals point outward)
            breplib.OrientClosedSolid(solid)
            fixed = solid
        except Exception:
            # Fall back to ShapeFix_Solid if OrientClosedSolid fails
            try:
                solid_fixer = ShapeFix_Solid(fixed)
                solid_fixer.SetPrecision(1e-6)
                solid_fixer.Perform()
                fixed = solid_fixer.Solid()
            except Exception:
                pass  # Keep the shape-fixed result if all else fails

    return fixed


def extract_segment_between_planes(
    solid: Any,
    plane1_point: List[float],
    plane1_normal: List[float],
    plane2_point: List[float],
    plane2_normal: List[float],
) -> Any:
    """Extract a segment of a solid between two cutting planes.

    This is designed for closed rings where a single half-space cut doesn't
    correctly identify the segment. We cut with both planes to isolate the
    segment between them.

    The normals should point OUTWARD from the segment (toward the parts to remove).

    Args:
        solid: yapCAD solid (typically a closed ring)
        plane1_point: Point on first cutting plane
        plane1_normal: Normal pointing away from segment (toward part to remove)
        plane2_point: Point on second cutting plane
        plane2_normal: Normal pointing away from segment (toward part to remove)

    Returns:
        The segment between the two planes
    """
    try:
        from OCC.Core.gp import gp_Pnt, gp_Dir, gp_Pln
        from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeHalfSpace
        from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut
        from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
        from yapcad.geom3d import solid as make_solid
        from yapcad.brep import brep_from_solid, BrepSolid, attach_brep_to_solid
    except ImportError as e:
        raise RuntimeError(
            "OCC (pythonocc-core) required for solid splitting"
        ) from e

    brep = brep_from_solid(solid)
    if brep is None or brep.shape is None:
        raise ValueError("Could not extract BREP from solid")

    current_shape = brep.shape

    # Cut with first plane (remove the part where normal1 points)
    pnt1 = gp_Pnt(plane1_point[0], plane1_point[1], plane1_point[2])
    dir1 = gp_Dir(plane1_normal[0], plane1_normal[1], plane1_normal[2])
    plane1 = gp_Pln(pnt1, dir1)
    face1 = BRepBuilderAPI_MakeFace(plane1).Face()

    # Half-space on the side where normal points (the part to remove)
    hs1 = BRepPrimAPI_MakeHalfSpace(face1, gp_Pnt(
        plane1_point[0] + plane1_normal[0],
        plane1_point[1] + plane1_normal[1],
        plane1_point[2] + plane1_normal[2]
    ))
    if not hs1.IsDone():
        raise RuntimeError("Failed to create half-space 1")

    cut1 = BRepAlgoAPI_Cut(current_shape, hs1.Solid())
    if not cut1.IsDone():
        raise RuntimeError("First plane cut failed")
    current_shape = cut1.Shape()

    # Cut with second plane (remove the part where normal2 points)
    pnt2 = gp_Pnt(plane2_point[0], plane2_point[1], plane2_point[2])
    dir2 = gp_Dir(plane2_normal[0], plane2_normal[1], plane2_normal[2])
    plane2 = gp_Pln(pnt2, dir2)
    face2 = BRepBuilderAPI_MakeFace(plane2).Face()

    hs2 = BRepPrimAPI_MakeHalfSpace(face2, gp_Pnt(
        plane2_point[0] + plane2_normal[0],
        plane2_point[1] + plane2_normal[1],
        plane2_point[2] + plane2_normal[2]
    ))
    if not hs2.IsDone():
        raise RuntimeError("Failed to create half-space 2")

    cut2 = BRepAlgoAPI_Cut(current_shape, hs2.Solid())
    if not cut2.IsDone():
        raise RuntimeError("Second plane cut failed")
    result_shape = cut2.Shape()

    # Fix the resulting shape
    result_shape = _fix_shape(result_shape)

    # Wrap as yapCAD solid
    result_brep = BrepSolid(result_shape)
    result_solid = make_solid([], [])
    attach_brep_to_solid(result_solid, result_brep)

    return result_solid


def segment_closed_ring(
    solid: Any,
    spine: Dict[str, Any],
    cut_parameters: List[float],
) -> List[Any]:
    """Segment a closed ring using two-plane extraction.

    For closed rings, single half-space cuts don't work correctly because
    the geometry wraps around. This function uses two cutting planes per
    segment to properly isolate each piece.

    Args:
        solid: The complete ring solid
        spine: Path3d representing the ring's sweep path
        cut_parameters: Sorted list of t values where cuts occur (0 < t < 1)

    Returns:
        List of segment solids in parameter order (segment 0 is [0, t0], etc.)

    Raises:
        RuntimeError: If extraction fails
        ValueError: If cut_parameters is empty or invalid
    """
    if not cut_parameters:
        raise ValueError("At least one cut parameter required")

    # Validate and sort cut parameters
    cuts = sorted(cut_parameters)
    for t in cuts:
        if t <= 0 or t >= 1:
            raise ValueError(f"Cut parameter {t} must be in range (0, 1)")

    segments = []
    n_cuts = len(cuts)

    # For a closed ring with cuts at [t1, t2, ...], we create segments:
    # Segment 0: from t=0 to t=cuts[0]
    # Segment 1: from t=cuts[0] to t=cuts[1]
    # ...
    # Segment n: from t=cuts[n-1] to t=1.0 (which wraps to 0)

    # Build list of segment boundaries
    # Each segment is defined by [start_t, end_t]
    segment_ranges = []
    segment_ranges.append((0.0, cuts[0]))  # First segment
    for i in range(n_cuts - 1):
        segment_ranges.append((cuts[i], cuts[i + 1]))
    segment_ranges.append((cuts[-1], 1.0))  # Last segment (wraps to start)

    for start_t, end_t in segment_ranges:
        # Get cut planes
        # For each segment, we need two planes that bound it
        # The normal should point AWAY from the segment (toward part to remove)

        # Plane at start_t: normal points backward (toward lower t, away from segment)
        start_point, start_normal = compute_cut_plane(spine, start_t)
        # Flip normal to point backward (toward lower t values)
        start_normal_out = [-start_normal[0], -start_normal[1], -start_normal[2]]

        # Plane at end_t: normal points forward (toward higher t, away from segment)
        end_point, end_normal = compute_cut_plane(spine, end_t)
        # Normal already points forward (toward higher t values)
        end_normal_out = end_normal

        # Special handling for first segment (start_t = 0)
        # and last segment (end_t = 1.0 which equals 0 for closed ring)
        # For these, we still use the computed planes, but the geometry
        # naturally handles the wrap-around

        try:
            segment = extract_segment_between_planes(
                solid,
                start_point, start_normal_out,
                end_point, end_normal_out,
            )
            segments.append(segment)
        except RuntimeError as e:
            raise RuntimeError(
                f"Failed to extract segment [{start_t}, {end_t}]: {e}"
            ) from e

    return segments


def split_solid_at_plane(
    solid: Any,
    plane_point: List[float],
    plane_normal: List[float],
) -> Tuple[Any, Any]:
    """Split a solid into two parts at a plane.

    Uses OCC's boolean operations to cut the solid with a half-space
    defined by the plane.

    NOTE: This works well for open paths (arcs, beams). For closed rings,
    use segment_closed_ring() or extract_segment_between_planes() instead.

    Args:
        solid: yapCAD solid to split
        plane_point: [x, y, z] point on the cutting plane
        plane_normal: [nx, ny, nz] normal vector (points toward "B" side)

    Returns:
        Tuple of (solid_a, solid_b) where:
        - solid_a: portion on negative side of plane (lower parameter)
        - solid_b: portion on positive side of plane (higher parameter)

    Raises:
        RuntimeError: If OCC is not available or splitting fails
    """
    try:
        from OCC.Core.gp import gp_Pnt, gp_Dir, gp_Pln
        from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeHalfSpace
        from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut
        from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
        from yapcad.geom3d import solid as make_solid
        from yapcad.brep import brep_from_solid, BrepSolid, attach_brep_to_solid
    except ImportError as e:
        raise RuntimeError(
            "OCC (pythonocc-core) required for solid splitting"
        ) from e

    # Get BREP representation of input solid
    brep = brep_from_solid(solid)
    if brep is None or brep.shape is None:
        raise ValueError("Could not extract BREP from solid")

    # Create cutting plane
    pnt = gp_Pnt(plane_point[0], plane_point[1], plane_point[2])
    direction = gp_Dir(plane_normal[0], plane_normal[1], plane_normal[2])
    plane = gp_Pln(pnt, direction)

    # Create a face on the plane (needed for half-space)
    face_builder = BRepBuilderAPI_MakeFace(plane)
    if not face_builder.IsDone():
        raise RuntimeError("Failed to create cutting plane face")
    plane_face = face_builder.Face()

    # Create half-spaces for both sides
    # Half-space on positive side of plane (where normal points)
    half_space_pos = BRepPrimAPI_MakeHalfSpace(plane_face, gp_Pnt(
        plane_point[0] + plane_normal[0],
        plane_point[1] + plane_normal[1],
        plane_point[2] + plane_normal[2]
    ))
    if not half_space_pos.IsDone():
        raise RuntimeError("Failed to create positive half-space")

    # Half-space on negative side of plane
    half_space_neg = BRepPrimAPI_MakeHalfSpace(plane_face, gp_Pnt(
        plane_point[0] - plane_normal[0],
        plane_point[1] - plane_normal[1],
        plane_point[2] - plane_normal[2]
    ))
    if not half_space_neg.IsDone():
        raise RuntimeError("Failed to create negative half-space")

    # Cut solid with each half-space to get the two parts
    # solid_a = solid ∩ negative_half_space (cut away positive side)
    cut_a = BRepAlgoAPI_Cut(brep.shape, half_space_pos.Solid())
    if not cut_a.IsDone():
        raise RuntimeError("Boolean cut for segment A failed")
    shape_a = cut_a.Shape()

    # solid_b = solid ∩ positive_half_space (cut away negative side)
    cut_b = BRepAlgoAPI_Cut(brep.shape, half_space_neg.Solid())
    if not cut_b.IsDone():
        raise RuntimeError("Boolean cut for segment B failed")
    shape_b = cut_b.Shape()

    # Fix shapes to correct face orientations after boolean operations
    shape_a = _fix_shape(shape_a)
    shape_b = _fix_shape(shape_b)

    # Wrap results as yapCAD solids with BREP using standard storage method
    # This ensures entityId is created and BREP is properly serialized/cached
    brep_a = BrepSolid(shape_a)
    solid_a = make_solid([], [])
    attach_brep_to_solid(solid_a, brep_a)

    brep_b = BrepSolid(shape_b)
    solid_b = make_solid([], [])
    attach_brep_to_solid(solid_b, brep_b)

    return solid_a, solid_b


def _union_connector_with_segment(
    segment_solid: Any,
    connector_solid: Any,
) -> Any:
    """Union a connector with a segment solid to create a male joint.

    Args:
        segment_solid: The segment to union connector with
        connector_solid: The interior connector solid

    Returns:
        Resulting solid with connector unioned

    Raises:
        RuntimeError: If union operation fails
    """
    from yapcad.geom3d import solid_boolean
    from yapcad.brep import brep_from_solid, BrepSolid, attach_brep_to_solid
    from yapcad.geom3d import solid as make_solid

    try:
        result = solid_boolean(segment_solid, connector_solid, 'union')

        # Apply shape fixing to the result to ensure correct face orientations
        brep = brep_from_solid(result)
        if brep is not None and brep.shape is not None:
            fixed_shape = _fix_shape(brep.shape)
            fixed_brep = BrepSolid(fixed_shape)
            fixed_solid = make_solid([], [])
            attach_brep_to_solid(fixed_solid, fixed_brep)
            return fixed_solid

        return result
    except Exception as e:
        raise RuntimeError(f"Connector union failed: {e}") from e


def segment_swept_element(
    provenance: SweptElementProvenance,
    cut_points: List[CutPoint],
    *,
    build_connectors: bool = True,
    union_connectors: bool = True,
) -> SegmentationResult:
    """Segment a swept element at the specified cut points.

    Takes a swept element (beam, pipe, etc.) with its provenance data
    and produces segments that can be individually manufactured and
    reassembled.

    Args:
        provenance: SweptElementProvenance with profile, spine, and metadata
        cut_points: List of CutPoint specifying where to cut
        build_connectors: Whether to generate interior connector solids
        union_connectors: Whether to union connectors with segments (male/female)
            When True, connectors are integrated into segments based on
            CutPoint.union_connector_with. When False, connectors are separate.

    Returns:
        SegmentationResult with segments, connectors, and assembly info
    """
    # Validate inputs
    if not cut_points:
        raise ValueError("At least one cut point required")

    # Sort cut points by parameter
    sorted_cuts = sorted(cut_points, key=lambda cp: cp.parameter)

    # Validate cut points are for this element
    for cp in sorted_cuts:
        if cp.element_id != provenance.id:
            raise ValueError(
                f"Cut point element_id '{cp.element_id}' doesn't match "
                f"provenance id '{provenance.id}'"
            )

    # We need the solid - rebuild it from provenance if needed
    # For now, assume provenance has metadata['solid'] from original creation
    solid = provenance.metadata.get('solid')
    if solid is None:
        raise ValueError(
            "Provenance must include 'solid' in metadata for segmentation"
        )

    # Get profile dimensions upfront if building connectors
    outer_w, outer_h = None, None
    if build_connectors and provenance.wall_thickness is not None:
        outer_w, outer_h = _extract_profile_dimensions(provenance.outer_profile)

    # Generate segments by iterating through cuts
    segments: List[Segment] = []
    connectors: List[ConnectorSpec] = []
    assembly_graph: Dict[str, List[str]] = {}
    warnings: List[str] = []

    # Track remaining solid as we make cuts
    remaining_solid = solid
    prev_segment_id: Optional[str] = None
    segment_start_t = 0.0

    # Track connector type for segments
    # When union_connector_with=="b", the next segment gets the connector
    next_segment_has_connector = False

    for i, cut_point in enumerate(sorted_cuts):
        # Compute cut plane from spine
        plane_point, plane_normal = compute_cut_plane(
            provenance.spine, cut_point.parameter
        )

        # Split the remaining solid
        try:
            segment_solid, remaining_solid = split_solid_at_plane(
                remaining_solid, plane_point, plane_normal
            )
        except RuntimeError as e:
            warnings.append(f"Cut {i} at t={cut_point.parameter} failed: {e}")
            continue

        # Determine this segment's connector type based on previous cut's decision
        # and build connector for this cut
        segment_connector_type = "none"
        segment_has_tab = False

        # If previous cut designated "b" (remaining), this segment has the connector
        if next_segment_has_connector:
            segment_connector_type = "male"
            segment_has_tab = True
            next_segment_has_connector = False

        # Build and integrate connector for THIS cut
        if build_connectors and outer_w and outer_h and provenance.wall_thickness:
            # Build connector specification
            conn_spec = compute_connector_spec(
                provenance.id,
                cut_point.parameter,
                outer_w,
                outer_h,
                provenance.wall_thickness,
                provenance.spine,
                fit_clearance=cut_point.fit_clearance,
            )

            if union_connectors and cut_point.union_connector_with != "none":
                # Build actual connector solid
                connector_solid = create_interior_connector(
                    outer_w,
                    outer_h,
                    provenance.spine,
                    cut_point.parameter,
                    wall_thickness=provenance.wall_thickness,
                    connector_length=conn_spec.length,
                    fit_clearance=cut_point.fit_clearance,
                )

                if cut_point.union_connector_with == "a":
                    # Union connector with segment A (this segment)
                    try:
                        segment_solid = _union_connector_with_segment(
                            segment_solid, connector_solid
                        )
                        segment_connector_type = "male"
                        segment_has_tab = True
                    except RuntimeError as e:
                        warnings.append(
                            f"Connector union at t={cut_point.parameter} failed: {e}"
                        )
                        # Fall back to separate connector
                        connectors.append(conn_spec)

                elif cut_point.union_connector_with == "b":
                    # Union connector with segment B (remaining/next segment)
                    try:
                        remaining_solid = _union_connector_with_segment(
                            remaining_solid, connector_solid
                        )
                        next_segment_has_connector = True
                        # This segment (A) is female - receives the connector from B
                        if segment_connector_type == "none":
                            segment_connector_type = "female"
                    except RuntimeError as e:
                        warnings.append(
                            f"Connector union at t={cut_point.parameter} failed: {e}"
                        )
                        connectors.append(conn_spec)
            else:
                # Keep connector as separate piece
                connectors.append(conn_spec)

        # Create segment
        segment_id = f"{provenance.id}_seg_{i}"
        segment = Segment(
            id=segment_id,
            solid=segment_solid,
            parent_element_id=provenance.id,
            parameter_range=(segment_start_t, cut_point.parameter),
            has_connector_tab=segment_has_tab,
            connector_type=segment_connector_type,
        )

        # Track mating relationships
        if prev_segment_id:
            segment.mates_with.append(prev_segment_id)
            # Update previous segment's mates
            prev_seg = next(
                (s for s in segments if s.id == prev_segment_id), None
            )
            if prev_seg:
                prev_seg.mates_with.append(segment_id)
                # If this segment is male and prev segment doesn't have a tab,
                # prev segment becomes female at this joint
                if segment_connector_type == "male" and prev_seg.connector_type == "none":
                    prev_seg.connector_type = "female"

        segments.append(segment)
        assembly_graph[segment_id] = list(segment.mates_with)

        prev_segment_id = segment_id
        segment_start_t = cut_point.parameter

    # Final segment (from last cut to end)
    final_segment_id = f"{provenance.id}_seg_{len(sorted_cuts)}"

    # Check if final segment has connector from last cut's "b" designation
    final_has_tab = next_segment_has_connector

    # Determine final segment's connector type
    if next_segment_has_connector:
        # Last cut used union_connector_with="b", so final segment is male
        final_connector_type = "male"
    elif prev_segment_id:
        # Check if previous segment is male - if so, final segment receives its tab
        prev_seg = next((s for s in segments if s.id == prev_segment_id), None)
        if prev_seg and prev_seg.connector_type == "male":
            final_connector_type = "female"
        else:
            final_connector_type = "none"
    else:
        final_connector_type = "none"

    final_segment = Segment(
        id=final_segment_id,
        solid=remaining_solid,
        parent_element_id=provenance.id,
        parameter_range=(segment_start_t, 1.0),
        has_connector_tab=final_has_tab,
        connector_type=final_connector_type,
    )

    if prev_segment_id:
        final_segment.mates_with.append(prev_segment_id)
        prev_seg = next((s for s in segments if s.id == prev_segment_id), None)
        if prev_seg:
            prev_seg.mates_with.append(final_segment_id)

    segments.append(final_segment)
    assembly_graph[final_segment_id] = list(final_segment.mates_with)

    # Generate assembly instructions
    instructions = _generate_assembly_instructions(segments, connectors)

    return SegmentationResult(
        segments=segments,
        connectors=connectors,
        assembly_graph=assembly_graph,
        build_volume_ok=True,  # TODO: validate against target volume
        warnings=warnings,
        assembly_instructions=instructions,
    )


def _extract_profile_dimensions(
    profile: Any,
) -> Tuple[Optional[float], Optional[float]]:
    """Extract width and height from a rectangular region2d profile.

    Args:
        profile: yapCAD region2d

    Returns:
        Tuple of (width, height) or (None, None) if not extractable
    """
    if not profile:
        return None, None

    try:
        from yapcad.geom import geomlistbbox

        # Get bounding box of outer profile boundary
        if isinstance(profile, list) and len(profile) > 0:
            outer = profile[0]  # First element is outer boundary
            bbox = geomlistbbox(outer)
            if bbox and len(bbox) == 2:
                width = bbox[1][0] - bbox[0][0]
                height = bbox[1][1] - bbox[0][1]
                return width, height
    except Exception:
        pass

    return None, None


def _generate_assembly_instructions(
    segments: List[Segment],
    connectors: List[ConnectorSpec],
) -> str:
    """Generate human-readable assembly instructions.

    Args:
        segments: List of segment objects
        connectors: List of connector specifications (may be empty if integrated)

    Returns:
        Multi-line string with assembly steps
    """
    lines = ["Assembly Instructions", "=" * 21, ""]

    lines.append(f"Total segments: {len(segments)}")

    # Count male/female segments
    male_count = sum(1 for s in segments if s.connector_type == "male")
    female_count = sum(1 for s in segments if s.connector_type == "female")

    if male_count > 0 or female_count > 0:
        lines.append(f"Male segments (with integrated connectors): {male_count}")
        lines.append(f"Female segments (receive connectors): {female_count}")
    if connectors:
        lines.append(f"Separate connectors: {len(connectors)}")
    lines.append("")

    lines.append("Assembly Sequence:")
    lines.append("-" * 18)

    for i, seg in enumerate(segments):
        step_num = i + 1
        lines.append(f"{step_num}. Place segment '{seg.id}'")

        # Show connector type if applicable
        if seg.connector_type == "male":
            lines.append("   Type: MALE (has integrated connector tab)")
        elif seg.connector_type == "female":
            lines.append("   Type: FEMALE (receives connector from mate)")

        if seg.mates_with:
            mates = ", ".join(seg.mates_with)
            lines.append(f"   Connects to: {mates}")

            # Describe mating
            for mate_id in seg.mates_with:
                mate_seg = next((s for s in segments if s.id == mate_id), None)
                if mate_seg:
                    if seg.connector_type == "male" and mate_seg.connector_type == "female":
                        lines.append(f"   -> Insert tab into '{mate_id}'")
                    elif seg.connector_type == "female" and mate_seg.connector_type == "male":
                        lines.append(f"   -> Receive tab from '{mate_id}'")

        # Check for separate connectors at this joint
        for conn in connectors:
            seg_end = seg.parameter_range[1]
            if abs(conn.center_parameter - seg_end) < 0.01:
                lines.append(f"   Insert connector '{conn.id}' before next segment")
                lines.append(f"   Fit type: {_fit_type_description(conn.fit_clearance)}")

        lines.append("")

    return "\n".join(lines)


def _fit_type_description(clearance: float) -> str:
    """Convert clearance value to human-readable fit type."""
    if clearance <= FIT_CLEARANCE['press']:
        return "press-fit (structural)"
    elif clearance <= FIT_CLEARANCE['slip']:
        return "slip-fit (easy assembly)"
    else:
        return "loose-fit (adjustable)"


def compute_optimal_cuts(
    provenance: SweptElementProvenance,
    max_segment_length: float,
    *,
    prefer_straight_cuts: bool = True,
) -> List[CutPoint]:
    """Compute optimal cut locations for a given max segment length.

    Analyzes the spine and proposes cut locations that:
    - Keep segments under max_segment_length
    - Prefer cuts at straight sections (not mid-curve)
    - Maintain structural integrity

    Args:
        provenance: SweptElementProvenance with spine data
        max_segment_length: Maximum length of any single segment
        prefer_straight_cuts: Try to cut at straight sections when possible

    Returns:
        List of CutPoint at optimal locations
    """
    spine = provenance.spine
    total_length = path_length(spine)

    if total_length <= max_segment_length:
        # No cuts needed
        return []

    # Simple approach: evenly spaced cuts
    # TODO: Enhanced version would analyze curvature and avoid mid-arc cuts
    num_segments = int(total_length / max_segment_length) + 1
    segment_length = total_length / num_segments

    cut_points = []
    for i in range(1, num_segments):
        # Convert length to parameter
        cut_length = i * segment_length
        cut_t = cut_length / total_length

        cut_points.append(CutPoint(
            element_id=provenance.id,
            parameter=cut_t,
            fit_clearance=FIT_CLEARANCE['press'],
        ))

    return cut_points


def build_connector_solids(
    provenance: SweptElementProvenance,
    connector_specs: List[ConnectorSpec],
) -> List[Tuple[ConnectorSpec, Any]]:
    """Build actual connector solids from specifications.

    Args:
        provenance: SweptElementProvenance with profile and spine data
        connector_specs: List of ConnectorSpec from segmentation

    Returns:
        List of (ConnectorSpec, solid) tuples
    """
    results = []

    # Get profile dimensions
    outer_w, outer_h = _extract_profile_dimensions(provenance.outer_profile)
    if outer_w is None or outer_h is None:
        raise ValueError("Cannot extract profile dimensions from provenance")

    wall_thickness = provenance.wall_thickness
    if wall_thickness is None:
        raise ValueError("Wall thickness required for connector generation")

    for spec in connector_specs:
        connector_solid = create_interior_connector(
            outer_w,
            outer_h,
            provenance.spine,
            spec.center_parameter,
            wall_thickness=wall_thickness,
            connector_length=spec.length,
            fit_clearance=spec.fit_clearance,
        )
        results.append((spec, connector_solid))

    return results
