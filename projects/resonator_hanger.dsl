module resonator_hanger

# All dimensions are in millimeters.

command gusset_base(width: float, height: float) -> region2d:
    require width > 0.0
    require height > 0.0

    p0: point = point(0.0,0.0)
    p1: point = point(0.0,-height)
    p2: point = point(width,0)

    pts: list<point> = [p0, p1, p2]
    base: region2d = polygon(pts)

    emit base

command gusset(length: float, width: float, height: float) -> solid:
    require height > 0.0

    base: region2d = gusset_base(width,height)

    # Extrude along Z, then rotate to get the gusset oriented correctly
    gusset0: solid = extrude(base, length)
    gusset1: solid = translate(gusset0, 0.0, 0.0, -length / 2.0)
    gusset2: solid = rotate(gusset1, -90.0, 0.0, -90.0)
    # gusset3: solid = rotate(gusset1, 0.0, 0.0, 0.0)

    emit gusset2

command GUSSET_BASE() -> region2d:
    base: region2d = gusset_base(5.0,7.0)
    emit base

command GUSSET(length: float, width: float, height: float) -> solid:
    emit gusset(length, width, height)


command box_tube_x(length: float, outer: float, wall: float) -> solid:
    require length > 0.0
    require outer > 0.0
    require wall > 0.0
    require outer > 2.0 * wall

    outer_box: solid = box(length, outer, outer)
    inner_size: float = outer - 2.0 * wall
    inner_box: solid = box(length + 2.0, inner_size, inner_size)
    emit difference(outer_box, inner_box)


command box_tube_z(length: float, outer: float, wall: float) -> solid:
    require length > 0.0
    require outer > 0.0
    require wall > 0.0
    require outer > 2.0 * wall

    outer_box: solid = box(outer, outer, length)
    inner_size: float = outer - 2.0 * wall
    inner_box: solid = box(inner_size, inner_size, length + 2.0)
    emit difference(outer_box, inner_box)


command cyl_x(radius: float, length: float) -> solid:
    require radius > 0.0
    require length > 0.0
    c0: solid = cylinder(radius, length)
    c1: solid = rotate(c0, 0.0, 90.0, 0.0)
    c2: solid = translate(c1, -length / 2.0, 0.0, 0.0)
    emit c2


command cyl_y(radius: float, length: float) -> solid:
    require radius > 0.0
    require length > 0.0
    c0: solid = cylinder(radius, length)
    c1: solid = rotate(c0, 90.0, 0.0, 0.0)
    c2: solid = translate(c1, 0.0, -length / 2.0, 0.0)
    emit c2


command add_reinforced_through_hole_y(
    base: solid,
    x: float,
    z: float,
    depth: float,
    hole_diameter: float,
    boss_wall: float
) -> solid:
    require depth > 0.0
    require hole_diameter > 0.0
    require boss_wall > 0.0

    hole_r: float = hole_diameter / 2.0
    boss_r: float = hole_r + boss_wall
    boss_len: float = depth + 4.0
    hole_len: float = depth + 6.0

    hole0: solid = cyl_y(hole_r, hole_len)
    # hole1: solid = translate(hole0, x, 0.0, z)
    hole1: solid = translate(hole0, x, hole_len , z)

    boss_outer0: solid = cyl_y(boss_r, boss_len)
    # boss_outer1: solid = translate(boss_outer0, x, 0.0, z)
    boss_outer1: solid = translate(boss_outer0, x, boss_len, z)
    boss_ring: solid = difference(boss_outer1, hole1)

    base_hole: solid = difference(base, hole1)
    out0: solid = union(base_hole, boss_ring)
    emit out0


command add_reinforced_hole_pair_y(
    base: solid,
    x_spacing: float,
    z: float,
    depth: float,
    hole_diameter: float,
    boss_wall: float
) -> solid:
    require x_spacing >= 0.0
    x0: float = 0.0 - x_spacing / 2.0
    x1: float = x_spacing / 2.0
    out0: solid = add_reinforced_through_hole_y(base, x0, z, depth, hole_diameter, boss_wall)
    out1: solid = add_reinforced_through_hole_y(out0, x1, z, depth, hole_diameter, boss_wall)
    emit out1


command MAKE_HANGER(
    overall_length: float = 160.0,
    clip_width: float = 30.0,
    clip_outer: float = 30.0,
    rail_length: float = 130.0,
    rail_outer: float = 30.0,
    wall: float = 3.0,
    riser_height: float = 10.0,
    # Plumbing 1/2" pipe OD is ~21.34mm (0.840"). Add clearance for sliding.
    pipe_od: float = 21.34,
    pipe_clearance: float = 0.3,
    pipe_top_wall: float = 3.0,
    pipe_opening: float = 14.0,
    # clip_tilt_deg: float = -45.0,
    clip_tilt_deg: float = 45.0,
    mount_x_spacing: float = 50.0,
    peg_od: float = 9.5,
    peg_clearance: float = 0.3,
    hole_boss_wall: float = 3.0
) -> solid:
    require overall_length > 2.0 * clip_width
    require rail_length > 0.0
    require rail_outer > 2.0 * wall
    require wall > 0.0
    require riser_height >= 0.0
    require pipe_od > 0.0
    require pipe_clearance >= 0.0
    require pipe_opening >= 0.0
    require pipe_opening <= clip_outer
    require peg_od > 0.0
    require peg_clearance >= 0.0

    hole_d: float = peg_od + peg_clearance
    hole_z: float = rail_outer / 2.0
    rail_depth: float = rail_outer

    rail0: solid = box_tube_x(rail_length, rail_outer, wall)
    rail1: solid = translate(rail0, 0.0, 0.0, hole_z)
    rail2: solid = add_reinforced_hole_pair_y(rail1, mount_x_spacing, hole_z, rail_depth, hole_d, hole_boss_wall)

    clip_bottom_z: float = rail_outer + riser_height
    clip_center_z: float = clip_bottom_z + clip_outer / 2.0

    clip_offset_x: float = (overall_length / 2.0) - (clip_width / 2.0)
    left_clip0: solid = box(clip_width, clip_outer, clip_outer)
    left_clip1: solid = translate(left_clip0, 0.0 - clip_offset_x, 0.0, clip_center_z)
    left_clip2: solid = translate(left_clip1, 0.0, 0.0, 0.0 - clip_center_z)
    left_clip3: solid = rotate(left_clip2, clip_tilt_deg, 0.0, 0.0)
    left_clip4: solid = translate(left_clip3, 0.0, 0.0, clip_center_z)

    right_clip0: solid = box(clip_width, clip_outer, clip_outer)
    right_clip1: solid = translate(right_clip0, clip_offset_x, 0.0, clip_center_z)
    right_clip2: solid = translate(right_clip1, 0.0, 0.0, 0.0 - clip_center_z)
    right_clip3: solid = rotate(right_clip2, clip_tilt_deg, 0.0, 0.0)
    right_clip4: solid = translate(right_clip3, 0.0, 0.0, clip_center_z)

    pipe_r: float = (pipe_od + pipe_clearance) / 2.0
    pipe_center_from_bottom: float = clip_outer - pipe_top_wall - pipe_r
    pipe_center_z: float = clip_bottom_z + pipe_center_from_bottom

    pipe_cut0: solid = cyl_x(pipe_r, clip_width + 4.0)
    pipe_cut_left: solid = translate(pipe_cut0, 0.0 - clip_offset_x, 0.0, pipe_center_z)
    left_clip5: solid = difference(left_clip4, pipe_cut_left)
    pipe_cut_right: solid = translate(pipe_cut0, clip_offset_x, 0.0, pipe_center_z)
    right_clip5: solid = difference(right_clip4, pipe_cut_right)

    # Side-entry opening: remove a slab from the -Y face to expose the pipe channel.
    # Reduce pipe_opening if you want more capture; increase it for easier insertion.
    opening_cut0: solid = box(clip_width + 6.0, pipe_opening, clip_outer + 6.0)
    opening_y: float = (0.0 - clip_outer / 2.0) + (pipe_opening / 2.0)

    opening_left0: solid = translate(opening_cut0, 0.0 - clip_offset_x, opening_y, clip_center_z)
    opening_left1: solid = translate(opening_left0, 0.0, 0.0, 0.0 - clip_center_z)
    opening_left2: solid = rotate(opening_left1, clip_tilt_deg, 0.0, 0.0)
    opening_left: solid = translate(opening_left2, 0.0, 0.0, clip_center_z)
    left_clip6: solid = difference(left_clip5, opening_left)

    opening_right0: solid = translate(opening_cut0, clip_offset_x, opening_y, clip_center_z)
    opening_right1: solid = translate(opening_right0, 0.0, 0.0, 0.0 - clip_center_z)
    opening_right2: solid = rotate(opening_right1, clip_tilt_deg, 0.0, 0.0)
    opening_right: solid = translate(opening_right2, 0.0, 0.0, clip_center_z)
    right_clip6: solid = difference(right_clip5, opening_right)

    riser0: solid = box(clip_width, rail_outer, riser_height)
    left_riser: solid = translate(riser0, 0.0 - clip_offset_x, 0.0, rail_outer + (riser_height / 2.0))
    right_riser: solid = translate(riser0, clip_offset_x, 0.0, rail_outer + (riser_height / 2.0))

    # Simple gusset block to better support the tilted clips (intentionally overbuilt for v1).
#    gusset0: solid = box(clip_width, rail_outer, clip_outer)
    gusset0: solid = gusset(clip_width, rail_outer/1.9, rail_outer/1.9)
    gusset1: solid = translate(gusset0,0,rail_outer/2,-rail_outer/6)
    left_gusset: solid = translate(gusset1, 0.0 - clip_offset_x, 0.0, rail_outer + (clip_outer / 2.0))
    right_gusset: solid = translate(gusset1, clip_offset_x, 0.0, rail_outer + (clip_outer / 2.0))

    # Union all parts using union_all
    parts: list<solid> = [rail2, left_riser, right_riser, left_gusset, right_gusset, left_clip6, right_clip6]
    result: solid = union_all(parts)
    emit result


command MAKE_C_RAIL(
    rail_length: float = 90.0,
    rail_outer: float = 30.0,
    wall: float = 3.0,
    span_between_rails: float = 306.0,
    mount_x_spacing: float = 50.0,
    peg_od: float = 9.5,
    peg_clearance: float = 0.3,
    hole_boss_wall: float = 3.0
) -> solid:
    require rail_length > 0.0
    require rail_outer > 2.0 * wall
    require span_between_rails >= 0.0

    hole_d: float = peg_od + peg_clearance
    depth: float = rail_outer

    # Build a "C" shape with the rails centered at X=0; the vertical leg sits on the left.
    x_leg_center: float = 0.0 - (rail_length / 2.0 - rail_outer / 2.0)

    bottom_rail0: solid = box_tube_x(rail_length, rail_outer, wall)
    bottom_rail1: solid = translate(bottom_rail0, 0.0, 0.0, rail_outer / 2.0)

    top_rail_z: float = rail_outer + span_between_rails + rail_outer / 2.0
    top_rail0: solid = box_tube_x(rail_length, rail_outer, wall)
    top_rail1: solid = translate(top_rail0, 0.0, 0.0, top_rail_z)

    vertical_len: float = 2.0 * rail_outer + span_between_rails
    vertical_leg0: solid = box_tube_z(vertical_len, rail_outer, wall)
    vertical_leg1: solid = translate(vertical_leg0, x_leg_center, 0.0, vertical_len / 2.0)

    # Union all rail parts using union_all
    rail_parts: list<solid> = [bottom_rail1, top_rail1, vertical_leg1]
    c_rail: solid = union_all(rail_parts)
    result: solid = add_reinforced_hole_pair_y(c_rail, mount_x_spacing, top_rail_z, depth, hole_d, hole_boss_wall)

    emit result


command MAKE_CLAMP(
    block_width: float = 15.0,
    block_height: float = 50.0,
    block_thickness: float = 10.0,
    peg_spacing_z: float = 30.0,
    peg_od: float = 9.5,
    peg_length: float = 14.0,
    taper_length: float = 5.0,
    taper_tip_od: float = 7.5
) -> solid:
    require block_width > 0.0
    require block_height > 0.0
    require block_thickness > 0.0
    require peg_spacing_z >= 0.0
    require peg_od > 0.0
    require peg_length > 0.0
    require taper_length > 0.0
    require taper_length <= peg_length
    require taper_tip_od > 0.0
    require taper_tip_od <= peg_od

    body0: solid = box(block_width, block_thickness, block_height)
    body: solid = translate(body0, 0.0, 0.0, block_height / 2.0)

    peg_r: float = peg_od / 2.0
    tip_r: float = taper_tip_od / 2.0
    straight_len: float = peg_length - taper_length

    peg0_z: float = (block_height / 2.0) - (peg_spacing_z / 2.0)
    peg1_z: float = (block_height / 2.0) + (peg_spacing_z / 2.0)

    peg_straight0: solid = cyl_y(peg_r, straight_len*1.1)
    # peg_straight_y: float = (block_thickness / 2.0) + (straight_len / 2.0)
    peg_straight_y: float = block_thickness + straight_len
    # peg_straight: solid = translate(peg_straight0, 0.0, peg_straight_y, peg0_z)
    peg_straight: solid = translate(peg_straight0, 0.0, peg_straight_y , peg0_z)

    peg_tip0: solid = cone(peg_r, tip_r, taper_length)
    peg_tip1: solid = rotate(peg_tip0, 90.0, 0.0, 180.0)
    peg_tip_y: float = (block_thickness / 2.0) + straight_len
    peg_tip: solid = translate(peg_tip1, 0.0, peg_tip_y, peg0_z)

    peg0: solid = union(peg_straight, peg_tip)
    peg1: solid = translate(peg0, 0.0, 0.0, peg1_z - peg0_z)

    # Union all parts using union_all
    parts: list<solid> = [body, peg0, peg1]
    result: solid = union_all(parts)
    emit result
