# Pipe Bending Jig - Sweep-based implementation
# Uses sweep operation instead of boolean unions for cleaner geometry
#
# This jig verifies 11.5 degree deflection bends in 1/2 inch OD copper pipe.
# Two bends, both left-hand turns, creating a shallow arc.
# Total length: approximately 13 inches (330mm)

module pipe_bending_jig_sweep

# Helper: create the outer box girder cross-section profile
command make_outer_profile(
    outer_width: float,
    outer_height: float
) -> region2d:
    profile: region2d = rectangle(outer_width, outer_height, point2d(0.0, 0.0))
    emit profile

# Helper: create the inner (hollow) box girder cross-section profile
# Note: bottom is solid, so inner profile is offset upward
command make_inner_profile(
    outer_width: float,
    outer_height: float,
    wall_thickness: float,
    bottom_thickness: float
) -> region2d:
    # Inner dimensions (walls on sides, solid bottom)
    inner_width: float = outer_width - 2.0 * wall_thickness
    inner_height: float = outer_height - wall_thickness - bottom_thickness

    # Center of inner rectangle is shifted up by (bottom - top_wall) / 2
    # Since we have solid bottom and wall on top:
    # bottom edge of inner = -outer_height/2 + bottom_thickness
    # top edge of inner = outer_height/2 - wall_thickness
    # center = (bottom_edge + top_edge) / 2
    center_offset: float = (bottom_thickness - wall_thickness) / 2.0

    profile: region2d = rectangle(inner_width, inner_height, point2d(0.0, center_offset))
    emit profile

# Helper: compute path points for the bent girder
# Returns the 4 vertices of the path (start, bend1, bend2, end)
# Path is in XY plane with Z=0
command make_girder_path(
    end_segment_length: float,
    middle_segment_length: float,
    deflection_angle: float
) -> path3d:
    # Convert angle to radians for trig
    angle_rad: float = radians(deflection_angle)
    double_angle_rad: float = radians(2.0 * deflection_angle)

    # Path starts at origin, goes along +Y
    p0_x: float = 0.0
    p0_y: float = 0.0
    p0_z: float = 0.0

    # First segment end (straight along Y)
    p1_x: float = 0.0
    p1_y: float = end_segment_length
    p1_z: float = 0.0

    # Direction after first bend (11.5 degrees from +Y toward +X)
    dir1_x: float = sin(angle_rad)
    dir1_y: float = cos(angle_rad)

    # Second bend point (after middle segment)
    p2_x: float = p1_x + middle_segment_length * dir1_x
    p2_y: float = p1_y + middle_segment_length * dir1_y
    p2_z: float = 0.0

    # Direction after second bend (23 degrees from +Y)
    dir2_x: float = sin(double_angle_rad)
    dir2_y: float = cos(double_angle_rad)

    # End point (after final segment)
    p3_x: float = p2_x + end_segment_length * dir2_x
    p3_y: float = p2_y + end_segment_length * dir2_y
    p3_z: float = 0.0

    # Create path segments
    seg1: path3d = path3d_line(point(p0_x, p0_y, p0_z), point(p1_x, p1_y, p1_z))
    seg2: path3d = path3d_line(point(p1_x, p1_y, p1_z), point(p2_x, p2_y, p2_z))
    seg3: path3d = path3d_line(point(p2_x, p2_y, p2_z), point(p3_x, p3_y, p3_z))

    # make_path3d is variadic - pass segments directly, not as list
    spine: path3d = make_path3d(seg1, seg2, seg3)

    emit spine

# Main command: create the pipe bending jig girder using sweep
# Default segment lengths adjusted so Y-projection = 330.2mm (13 inches)
command MAKE_GIRDER_SWEEP(
    end_segment_length: float = 52.2,
    middle_segment_length: float = 234.7,
    deflection_angle: float = 11.5,
    girder_width: float = 12.7,
    girder_height: float = 12.7,
    wall_thickness: float = 2.54,
    bottom_thickness: float = 2.54
) -> solid:
    # Create the outer cross-section profile
    outer_profile: region2d = make_outer_profile(girder_width, girder_height)

    # Create the inner cross-section profile (for hollow girder)
    inner_profile: region2d = make_inner_profile(girder_width, girder_height, wall_thickness, bottom_thickness)

    # Create the sweep path
    spine: path3d = make_girder_path(end_segment_length, middle_segment_length, deflection_angle)

    # Sweep hollow profile along path to create box girder
    girder: solid = sweep_hollow(outer_profile, inner_profile, spine)

    emit girder

# Adaptive sweep version: profile tracks tangent with 5 degree threshold
command MAKE_GIRDER_SWEEP_ADAPTIVE(
    end_segment_length: float = 52.2,
    middle_segment_length: float = 234.7,
    deflection_angle: float = 11.5,
    girder_width: float = 12.7,
    girder_height: float = 12.7,
    angle_threshold: float = 5.0
) -> solid:
    # Create the cross-section profile
    profile: region2d = make_outer_profile(girder_width, girder_height)

    # Create the sweep path
    spine: path3d = make_girder_path(end_segment_length, middle_segment_length, deflection_angle)

    # Sweep with adaptive tangent tracking
    # Profile rotates to stay perpendicular to path tangent
    girder: solid = sweep_adaptive(profile, spine, angle_threshold)

    emit girder

# Adaptive sweep version with hollow profile (box girder)
command MAKE_GIRDER_SWEEP_ADAPTIVE_HOLLOW(
    end_segment_length: float = 52.2,
    middle_segment_length: float = 234.7,
    deflection_angle: float = 11.5,
    girder_width: float = 12.7,
    girder_height: float = 12.7,
    wall_thickness: float = 2.54,
    bottom_thickness: float = 2.54,
    angle_threshold: float = 5.0
) -> solid:
    # Create the outer cross-section profile
    outer_profile: region2d = make_outer_profile(girder_width, girder_height)

    # Create the inner cross-section profile (for hollow girder)
    inner_profile: region2d = make_inner_profile(girder_width, girder_height, wall_thickness, bottom_thickness)

    # Create the sweep path
    spine: path3d = make_girder_path(end_segment_length, middle_segment_length, deflection_angle)

    # Sweep hollow profile with adaptive tangent tracking
    # Profile rotates to stay perpendicular to path tangent
    girder: solid = sweep_adaptive_hollow(outer_profile, inner_profile, spine, angle_threshold)

    emit girder

# Helper: create a single pipe cradle (120 degree C-shape)
# The cradle is oriented with the opening facing +Z (upward)
# and the pipe runs along the Y axis
command make_cradle(
    cradle_length: float,
    cradle_diameter: float
) -> solid:
    # Create the cradle block
    block_width: float = cradle_diameter + 4.0
    block_height: float = cradle_diameter / 2.0 + 2.0

    cradle_block: solid = box(block_width, cradle_length, block_height)

    # Create cylinder to subtract (along Y axis for pipe)
    # Cylinder is longer than block to ensure clean cut
    cutter_length: float = cradle_length + 2.0
    pipe_cutter: solid = cylinder(cradle_diameter / 2.0, cutter_length)
    # Rotate to align along Y axis (cylinder default is Z)
    pipe_rotated: solid = rotate(pipe_cutter, 90.0, 0.0, 0.0)
    # Center along Y
    half_cutter: float = cutter_length / 2.0
    pipe_centered: solid = translate(pipe_rotated, 0.0, half_cutter, 0.0)

    # Position cylinder so 120 degrees is captured
    cylinder_z_offset: float = block_height / 2.0 - (cradle_diameter / 2.0) * 0.5
    pipe_positioned: solid = translate(pipe_centered, 0.0, 0.0, cylinder_z_offset)

    # Subtract cylinder from block
    cradle: solid = difference(cradle_block, pipe_positioned)

    emit cradle

# Full jig with adaptive hollow girder and 7 cradles positioned along the path
# Cradle positions: 2 on first segment, 3 on middle segment, 2 on last segment
# Positions are computed directly along the path geometry
command MAKE_JIG_SWEEP_ADAPTIVE(
    end_segment_length: float = 52.2,
    middle_segment_length: float = 234.7,
    deflection_angle: float = 11.5,
    girder_width: float = 12.7,
    girder_height: float = 12.7,
    wall_thickness: float = 2.54,
    bottom_thickness: float = 2.54,
    cradle_length: float = 12.7,
    cradle_diameter: float = 12.8016,
    angle_threshold: float = 5.0
) -> solid:
    # Create the girder using adaptive hollow sweep
    outer_profile: region2d = make_outer_profile(girder_width, girder_height)
    inner_profile: region2d = make_inner_profile(girder_width, girder_height, wall_thickness, bottom_thickness)
    spine: path3d = make_girder_path(end_segment_length, middle_segment_length, deflection_angle)
    girder: solid = sweep_adaptive_hollow(outer_profile, inner_profile, spine, angle_threshold)

    # Create base cradle
    base_cradle: solid = make_cradle(cradle_length, cradle_diameter)

    # Cradle Z offset (on top of girder)
    cradle_block_height: float = cradle_diameter / 2.0 + 2.0
    z_offset: float = girder_height / 2.0 + cradle_block_height / 2.0 - 1.0

    # Path geometry calculations (same as make_girder_path)
    angle_rad: float = radians(deflection_angle)
    double_angle_rad: float = radians(2.0 * deflection_angle)

    # Path vertices
    # p0 = (0, 0) - start
    # p1 = (0, end_segment_length) - first bend
    p1_x: float = 0.0
    p1_y: float = end_segment_length

    # Direction after first bend
    dir1_x: float = sin(angle_rad)
    dir1_y: float = cos(angle_rad)

    # p2 = second bend
    p2_x: float = p1_x + middle_segment_length * dir1_x
    p2_y: float = p1_y + middle_segment_length * dir1_y

    # Direction after second bend
    dir2_x: float = sin(double_angle_rad)
    dir2_y: float = cos(double_angle_rad)

    # Cradle positions along each segment
    # Segment 1: cradles at 25% and 75% of end_segment_length
    # These are on the first straight segment (Y axis, angle=0)
    # Cradle groove runs along Y, so no rotation needed for segment 1
    c1_y: float = end_segment_length * 0.25
    c2_y: float = end_segment_length * 0.75

    c1: solid = translate(base_cradle, 0.0, c1_y, z_offset)
    c2: solid = translate(base_cradle, 0.0, c2_y, z_offset)

    # Segment 2: cradles at 25%, 50%, 75% of middle_segment_length
    # These are on the angled middle segment
    # Path tangent is at deflection_angle from +Y toward +X
    # Rotate cradle around Z by -deflection_angle to align groove with tangent
    # (negative because rotate() uses right-hand rule, and we want clockwise when viewed from +Z)
    neg_angle1: float = 0.0 - deflection_angle

    d3: float = middle_segment_length * 0.25
    c3_x: float = p1_x + d3 * dir1_x
    c3_y: float = p1_y + d3 * dir1_y

    d4: float = middle_segment_length * 0.5
    c4_x: float = p1_x + d4 * dir1_x
    c4_y: float = p1_y + d4 * dir1_y

    d5: float = middle_segment_length * 0.75
    c5_x: float = p1_x + d5 * dir1_x
    c5_y: float = p1_y + d5 * dir1_y

    c3_rot: solid = rotate(base_cradle, 0.0, 0.0, neg_angle1)
    c3: solid = translate(c3_rot, c3_x, c3_y, z_offset)

    c4_rot: solid = rotate(base_cradle, 0.0, 0.0, neg_angle1)
    c4: solid = translate(c4_rot, c4_x, c4_y, z_offset)

    c5_rot: solid = rotate(base_cradle, 0.0, 0.0, neg_angle1)
    c5: solid = translate(c5_rot, c5_x, c5_y, z_offset)

    # Segment 3: cradles at 25% and 75% of end_segment_length
    # These are on the final segment (angle = 2 * deflection_angle from +Y)
    neg_angle2: float = 0.0 - 2.0 * deflection_angle

    d6: float = end_segment_length * 0.25
    c6_x: float = p2_x + d6 * dir2_x
    c6_y: float = p2_y + d6 * dir2_y

    d7: float = end_segment_length * 0.75
    c7_x: float = p2_x + d7 * dir2_x
    c7_y: float = p2_y + d7 * dir2_y

    c6_rot: solid = rotate(base_cradle, 0.0, 0.0, neg_angle2)
    c6: solid = translate(c6_rot, c6_x, c6_y, z_offset)

    c7_rot: solid = rotate(base_cradle, 0.0, 0.0, neg_angle2)
    c7: solid = translate(c7_rot, c7_x, c7_y, z_offset)

    # Union all parts using union_all
    parts: list<solid> = [girder, c1, c2, c3, c4, c5, c6, c7]
    jig: solid = union_all(parts)

    emit jig
