module globe_stand

# Globe Stand with Box Beam Construction
#
# A stand for a 30.48cm (12 inch) diameter globe featuring:
# - Flat base ring made from box beam
# - Tilted cradle ring to hold the globe (Mars axial tilt: 25.2 degrees)
# - Three parabolic arc supports connecting base to cradle
# - All beams are hollow square profile (10mm outer, 2mm wall)

# Pi constant for degree-to-radian conversion
# Note: cos/sin in DSL use radians


# ============================================================================
# Helper: degree to radian conversion
# ============================================================================

command deg2rad(degrees: float) -> float:
    pi: float = 3.14159265359
    emit degrees * pi / 180.0


# ============================================================================
# Profile definitions
# ============================================================================

# Hollow square profile for box beam construction
command box_beam_profile(outer_size: float, wall_thickness: float) -> region2d:
    require outer_size > 0.0
    require wall_thickness > 0.0
    require wall_thickness < outer_size / 2.0

    half: float = outer_size / 2.0

    p0: point = point(-half, -half)
    p1: point = point(half, -half)
    p2: point = point(half, half)
    p3: point = point(-half, half)

    pts: list<point> = [p0, p1, p2, p3]
    emit polygon(pts)


command box_beam_inner(outer_size: float, wall_thickness: float) -> region2d:
    require outer_size > 0.0
    require wall_thickness > 0.0
    require wall_thickness < outer_size / 2.0

    half: float = outer_size / 2.0
    inner_half: float = half - wall_thickness

    q0: point = point(-inner_half, -inner_half)
    q1: point = point(inner_half, -inner_half)
    q2: point = point(inner_half, inner_half)
    q3: point = point(-inner_half, inner_half)

    inner_pts: list<point> = [q0, q1, q2, q3]
    emit polygon(inner_pts)


# ============================================================================
# Ring path generation
# ============================================================================

# Create a circular ring path in the XY plane at height z
# Uses 4 arc segments (90 degrees each)
command ring_path(radius: float, z_height: float) -> path3d:
    require radius > 0.0

    pi: float = 3.14159265359

    # Four 90-degree arc segments
    # Segment 0: 0 to 90 degrees
    seg0: path3d = path3d_arc_auto(
        point(0.0, 0.0, z_height),
        point(radius, 0.0, z_height),
        point(0.0, radius, z_height),
        false
    )

    # Segment 1: 90 to 180 degrees
    seg1: path3d = path3d_arc_auto(
        point(0.0, 0.0, z_height),
        point(0.0, radius, z_height),
        point(-radius, 0.0, z_height),
        false
    )

    # Segment 2: 180 to 270 degrees
    seg2: path3d = path3d_arc_auto(
        point(0.0, 0.0, z_height),
        point(-radius, 0.0, z_height),
        point(0.0, -radius, z_height),
        false
    )

    # Segment 3: 270 to 360 degrees (closes the ring)
    seg3: path3d = path3d_arc_auto(
        point(0.0, 0.0, z_height),
        point(0.0, -radius, z_height),
        point(radius, 0.0, z_height),
        false
    )

    full_ring: path3d = make_path3d(seg0, seg1, seg2, seg3)
    emit full_ring


# Create a tilted circular ring path
# The ring is tilted around the X axis by tilt_deg degrees
# center_z is the Z height of the ring center
command tilted_ring_path(
    radius: float,
    center_z: float,
    tilt_deg: float
) -> path3d:
    require radius > 0.0

    pi: float = 3.14159265359


# ============================================================================
# Latitude-based cradle ring path
# ============================================================================

# Create a cradle ring path that follows a parallel of latitude on a tilted oblate spheroid
# The path is offset outward from the globe surface by beam_offset (typically beam_outer/2)
# globe_center_z is the Z position of the globe center
#
# This generates the ring as a series of line segments (16 segments for smooth approximation)
command latitude_cradle_path(
    globe_diameter: float,
    oblateness: float,
    latitude_deg: float,
    tilt_deg: float,
    globe_center_z: float,
    beam_offset: float
) -> path3d:
    require globe_diameter > 0.0
    require oblateness >= 0.0
    require oblateness < 1.0

    pi: float = 3.14159265359

    globe_radius: float = globe_diameter / 2.0
    polar_radius: float = globe_radius * (1.0 - oblateness)

    lat_rad: float = latitude_deg * pi / 180.0
    tilt_rad: float = tilt_deg * pi / 180.0

    cos_lat: float = cos(lat_rad)
    sin_lat: float = sin(lat_rad)
    cos_tilt: float = cos(tilt_rad)
    sin_tilt: float = sin(tilt_rad)

    # For normalization of surface normals
    a2: float = globe_radius * globe_radius
    c2: float = polar_radius * polar_radius

    # Generate 16 points around the latitude circle (every 22.5 degrees)
    # Point at longitude 0 degrees
    x0: float = globe_radius * cos_lat * 1.0  # cos(0) = 1
    y0: float = globe_radius * cos_lat * 0.0  # sin(0) = 0
    z0: float = polar_radius * sin_lat
    nx0: float = x0 / a2
    ny0: float = y0 / a2
    nz0: float = z0 / c2
    nlen0: float = sqrt(nx0*nx0 + ny0*ny0 + nz0*nz0)
    px0: float = x0 + beam_offset * nx0 / nlen0
    py0: float = y0 + beam_offset * ny0 / nlen0
    pz0: float = z0 + beam_offset * nz0 / nlen0
    # Apply tilt and translate
    rx0: float = px0
    ry0: float = py0 * cos_tilt - pz0 * sin_tilt
    rz0: float = py0 * sin_tilt + pz0 * cos_tilt + globe_center_z

    # Point at longitude 22.5 degrees
    cos22: float = 0.92387953251  # cos(22.5°)
    sin22: float = 0.38268343236  # sin(22.5°)
    x1: float = globe_radius * cos_lat * cos22
    y1: float = globe_radius * cos_lat * sin22
    z1: float = polar_radius * sin_lat
    nx1: float = x1 / a2
    ny1: float = y1 / a2
    nz1: float = z1 / c2
    nlen1: float = sqrt(nx1*nx1 + ny1*ny1 + nz1*nz1)
    px1: float = x1 + beam_offset * nx1 / nlen1
    py1: float = y1 + beam_offset * ny1 / nlen1
    pz1: float = z1 + beam_offset * nz1 / nlen1
    rx1: float = px1
    ry1: float = py1 * cos_tilt - pz1 * sin_tilt
    rz1: float = py1 * sin_tilt + pz1 * cos_tilt + globe_center_z

    # Point at longitude 45 degrees
    cos45: float = 0.70710678118  # cos(45°)
    sin45: float = 0.70710678118  # sin(45°)
    x2: float = globe_radius * cos_lat * cos45
    y2: float = globe_radius * cos_lat * sin45
    z2: float = polar_radius * sin_lat
    nx2: float = x2 / a2
    ny2: float = y2 / a2
    nz2: float = z2 / c2
    nlen2: float = sqrt(nx2*nx2 + ny2*ny2 + nz2*nz2)
    px2: float = x2 + beam_offset * nx2 / nlen2
    py2: float = y2 + beam_offset * ny2 / nlen2
    pz2: float = z2 + beam_offset * nz2 / nlen2
    rx2: float = px2
    ry2: float = py2 * cos_tilt - pz2 * sin_tilt
    rz2: float = py2 * sin_tilt + pz2 * cos_tilt + globe_center_z

    # Point at longitude 67.5 degrees
    cos67: float = 0.38268343236  # cos(67.5°)
    sin67: float = 0.92387953251  # sin(67.5°)
    x3: float = globe_radius * cos_lat * cos67
    y3: float = globe_radius * cos_lat * sin67
    z3: float = polar_radius * sin_lat
    nx3: float = x3 / a2
    ny3: float = y3 / a2
    nz3: float = z3 / c2
    nlen3: float = sqrt(nx3*nx3 + ny3*ny3 + nz3*nz3)
    px3: float = x3 + beam_offset * nx3 / nlen3
    py3: float = y3 + beam_offset * ny3 / nlen3
    pz3: float = z3 + beam_offset * nz3 / nlen3
    rx3: float = px3
    ry3: float = py3 * cos_tilt - pz3 * sin_tilt
    rz3: float = py3 * sin_tilt + pz3 * cos_tilt + globe_center_z

    # Point at longitude 90 degrees
    x4: float = globe_radius * cos_lat * 0.0  # cos(90) = 0
    y4: float = globe_radius * cos_lat * 1.0  # sin(90) = 1
    z4: float = polar_radius * sin_lat
    nx4: float = x4 / a2
    ny4: float = y4 / a2
    nz4: float = z4 / c2
    nlen4: float = sqrt(nx4*nx4 + ny4*ny4 + nz4*nz4)
    px4: float = x4 + beam_offset * nx4 / nlen4
    py4: float = y4 + beam_offset * ny4 / nlen4
    pz4: float = z4 + beam_offset * nz4 / nlen4
    rx4: float = px4
    ry4: float = py4 * cos_tilt - pz4 * sin_tilt
    rz4: float = py4 * sin_tilt + pz4 * cos_tilt + globe_center_z

    # Point at longitude 112.5 degrees
    cos112: float = -0.38268343236  # cos(112.5°)
    sin112: float = 0.92387953251   # sin(112.5°)
    x5: float = globe_radius * cos_lat * cos112
    y5: float = globe_radius * cos_lat * sin112
    z5: float = polar_radius * sin_lat
    nx5: float = x5 / a2
    ny5: float = y5 / a2
    nz5: float = z5 / c2
    nlen5: float = sqrt(nx5*nx5 + ny5*ny5 + nz5*nz5)
    px5: float = x5 + beam_offset * nx5 / nlen5
    py5: float = y5 + beam_offset * ny5 / nlen5
    pz5: float = z5 + beam_offset * nz5 / nlen5
    rx5: float = px5
    ry5: float = py5 * cos_tilt - pz5 * sin_tilt
    rz5: float = py5 * sin_tilt + pz5 * cos_tilt + globe_center_z

    # Point at longitude 135 degrees
    cos135: float = -0.70710678118  # cos(135°)
    sin135: float = 0.70710678118   # sin(135°)
    x6: float = globe_radius * cos_lat * cos135
    y6: float = globe_radius * cos_lat * sin135
    z6: float = polar_radius * sin_lat
    nx6: float = x6 / a2
    ny6: float = y6 / a2
    nz6: float = z6 / c2
    nlen6: float = sqrt(nx6*nx6 + ny6*ny6 + nz6*nz6)
    px6: float = x6 + beam_offset * nx6 / nlen6
    py6: float = y6 + beam_offset * ny6 / nlen6
    pz6: float = z6 + beam_offset * nz6 / nlen6
    rx6: float = px6
    ry6: float = py6 * cos_tilt - pz6 * sin_tilt
    rz6: float = py6 * sin_tilt + pz6 * cos_tilt + globe_center_z

    # Point at longitude 157.5 degrees
    cos157: float = -0.92387953251  # cos(157.5°)
    sin157: float = 0.38268343236   # sin(157.5°)
    x7: float = globe_radius * cos_lat * cos157
    y7: float = globe_radius * cos_lat * sin157
    z7: float = polar_radius * sin_lat
    nx7: float = x7 / a2
    ny7: float = y7 / a2
    nz7: float = z7 / c2
    nlen7: float = sqrt(nx7*nx7 + ny7*ny7 + nz7*nz7)
    px7: float = x7 + beam_offset * nx7 / nlen7
    py7: float = y7 + beam_offset * ny7 / nlen7
    pz7: float = z7 + beam_offset * nz7 / nlen7
    rx7: float = px7
    ry7: float = py7 * cos_tilt - pz7 * sin_tilt
    rz7: float = py7 * sin_tilt + pz7 * cos_tilt + globe_center_z

    # Point at longitude 180 degrees
    x8: float = globe_radius * cos_lat * -1.0  # cos(180) = -1
    y8: float = globe_radius * cos_lat * 0.0   # sin(180) = 0
    z8: float = polar_radius * sin_lat
    nx8: float = x8 / a2
    ny8: float = y8 / a2
    nz8: float = z8 / c2
    nlen8: float = sqrt(nx8*nx8 + ny8*ny8 + nz8*nz8)
    px8: float = x8 + beam_offset * nx8 / nlen8
    py8: float = y8 + beam_offset * ny8 / nlen8
    pz8: float = z8 + beam_offset * nz8 / nlen8
    rx8: float = px8
    ry8: float = py8 * cos_tilt - pz8 * sin_tilt
    rz8: float = py8 * sin_tilt + pz8 * cos_tilt + globe_center_z

    # Point at longitude 202.5 degrees
    cos202: float = -0.92387953251  # cos(202.5°)
    sin202: float = -0.38268343236  # sin(202.5°)
    x9: float = globe_radius * cos_lat * cos202
    y9: float = globe_radius * cos_lat * sin202
    z9: float = polar_radius * sin_lat
    nx9: float = x9 / a2
    ny9: float = y9 / a2
    nz9: float = z9 / c2
    nlen9: float = sqrt(nx9*nx9 + ny9*ny9 + nz9*nz9)
    px9: float = x9 + beam_offset * nx9 / nlen9
    py9: float = y9 + beam_offset * ny9 / nlen9
    pz9: float = z9 + beam_offset * nz9 / nlen9
    rx9: float = px9
    ry9: float = py9 * cos_tilt - pz9 * sin_tilt
    rz9: float = py9 * sin_tilt + pz9 * cos_tilt + globe_center_z

    # Point at longitude 225 degrees
    cos225: float = -0.70710678118  # cos(225°)
    sin225: float = -0.70710678118  # sin(225°)
    x10: float = globe_radius * cos_lat * cos225
    y10: float = globe_radius * cos_lat * sin225
    z10: float = polar_radius * sin_lat
    nx10: float = x10 / a2
    ny10: float = y10 / a2
    nz10: float = z10 / c2
    nlen10: float = sqrt(nx10*nx10 + ny10*ny10 + nz10*nz10)
    px10: float = x10 + beam_offset * nx10 / nlen10
    py10: float = y10 + beam_offset * ny10 / nlen10
    pz10: float = z10 + beam_offset * nz10 / nlen10
    rx10: float = px10
    ry10: float = py10 * cos_tilt - pz10 * sin_tilt
    rz10: float = py10 * sin_tilt + pz10 * cos_tilt + globe_center_z

    # Point at longitude 247.5 degrees
    cos247: float = -0.38268343236  # cos(247.5°)
    sin247: float = -0.92387953251  # sin(247.5°)
    x11: float = globe_radius * cos_lat * cos247
    y11: float = globe_radius * cos_lat * sin247
    z11: float = polar_radius * sin_lat
    nx11: float = x11 / a2
    ny11: float = y11 / a2
    nz11: float = z11 / c2
    nlen11: float = sqrt(nx11*nx11 + ny11*ny11 + nz11*nz11)
    px11: float = x11 + beam_offset * nx11 / nlen11
    py11: float = y11 + beam_offset * ny11 / nlen11
    pz11: float = z11 + beam_offset * nz11 / nlen11
    rx11: float = px11
    ry11: float = py11 * cos_tilt - pz11 * sin_tilt
    rz11: float = py11 * sin_tilt + pz11 * cos_tilt + globe_center_z

    # Point at longitude 270 degrees
    x12: float = globe_radius * cos_lat * 0.0   # cos(270) = 0
    y12: float = globe_radius * cos_lat * -1.0  # sin(270) = -1
    z12: float = polar_radius * sin_lat
    nx12: float = x12 / a2
    ny12: float = y12 / a2
    nz12: float = z12 / c2
    nlen12: float = sqrt(nx12*nx12 + ny12*ny12 + nz12*nz12)
    px12: float = x12 + beam_offset * nx12 / nlen12
    py12: float = y12 + beam_offset * ny12 / nlen12
    pz12: float = z12 + beam_offset * nz12 / nlen12
    rx12: float = px12
    ry12: float = py12 * cos_tilt - pz12 * sin_tilt
    rz12: float = py12 * sin_tilt + pz12 * cos_tilt + globe_center_z

    # Point at longitude 292.5 degrees
    cos292: float = 0.38268343236   # cos(292.5°)
    sin292: float = -0.92387953251  # sin(292.5°)
    x13: float = globe_radius * cos_lat * cos292
    y13: float = globe_radius * cos_lat * sin292
    z13: float = polar_radius * sin_lat
    nx13: float = x13 / a2
    ny13: float = y13 / a2
    nz13: float = z13 / c2
    nlen13: float = sqrt(nx13*nx13 + ny13*ny13 + nz13*nz13)
    px13: float = x13 + beam_offset * nx13 / nlen13
    py13: float = y13 + beam_offset * ny13 / nlen13
    pz13: float = z13 + beam_offset * nz13 / nlen13
    rx13: float = px13
    ry13: float = py13 * cos_tilt - pz13 * sin_tilt
    rz13: float = py13 * sin_tilt + pz13 * cos_tilt + globe_center_z

    # Point at longitude 315 degrees
    cos315: float = 0.70710678118   # cos(315°)
    sin315: float = -0.70710678118  # sin(315°)
    x14: float = globe_radius * cos_lat * cos315
    y14: float = globe_radius * cos_lat * sin315
    z14: float = polar_radius * sin_lat
    nx14: float = x14 / a2
    ny14: float = y14 / a2
    nz14: float = z14 / c2
    nlen14: float = sqrt(nx14*nx14 + ny14*ny14 + nz14*nz14)
    px14: float = x14 + beam_offset * nx14 / nlen14
    py14: float = y14 + beam_offset * ny14 / nlen14
    pz14: float = z14 + beam_offset * nz14 / nlen14
    rx14: float = px14
    ry14: float = py14 * cos_tilt - pz14 * sin_tilt
    rz14: float = py14 * sin_tilt + pz14 * cos_tilt + globe_center_z

    # Point at longitude 337.5 degrees
    cos337: float = 0.92387953251   # cos(337.5°)
    sin337: float = -0.38268343236  # sin(337.5°)
    x15: float = globe_radius * cos_lat * cos337
    y15: float = globe_radius * cos_lat * sin337
    z15: float = polar_radius * sin_lat
    nx15: float = x15 / a2
    ny15: float = y15 / a2
    nz15: float = z15 / c2
    nlen15: float = sqrt(nx15*nx15 + ny15*ny15 + nz15*nz15)
    px15: float = x15 + beam_offset * nx15 / nlen15
    py15: float = y15 + beam_offset * ny15 / nlen15
    pz15: float = z15 + beam_offset * nz15 / nlen15
    rx15: float = px15
    ry15: float = py15 * cos_tilt - pz15 * sin_tilt
    rz15: float = py15 * sin_tilt + pz15 * cos_tilt + globe_center_z

    # Build path segments (16 line segments forming closed loop)
    seg0: path3d = path3d_line(point(rx0, ry0, rz0), point(rx1, ry1, rz1))
    seg1: path3d = path3d_line(point(rx1, ry1, rz1), point(rx2, ry2, rz2))
    seg2: path3d = path3d_line(point(rx2, ry2, rz2), point(rx3, ry3, rz3))
    seg3: path3d = path3d_line(point(rx3, ry3, rz3), point(rx4, ry4, rz4))
    seg4: path3d = path3d_line(point(rx4, ry4, rz4), point(rx5, ry5, rz5))
    seg5: path3d = path3d_line(point(rx5, ry5, rz5), point(rx6, ry6, rz6))
    seg6: path3d = path3d_line(point(rx6, ry6, rz6), point(rx7, ry7, rz7))
    seg7: path3d = path3d_line(point(rx7, ry7, rz7), point(rx8, ry8, rz8))
    seg8: path3d = path3d_line(point(rx8, ry8, rz8), point(rx9, ry9, rz9))
    seg9: path3d = path3d_line(point(rx9, ry9, rz9), point(rx10, ry10, rz10))
    seg10: path3d = path3d_line(point(rx10, ry10, rz10), point(rx11, ry11, rz11))
    seg11: path3d = path3d_line(point(rx11, ry11, rz11), point(rx12, ry12, rz12))
    seg12: path3d = path3d_line(point(rx12, ry12, rz12), point(rx13, ry13, rz13))
    seg13: path3d = path3d_line(point(rx13, ry13, rz13), point(rx14, ry14, rz14))
    seg14: path3d = path3d_line(point(rx14, ry14, rz14), point(rx15, ry15, rz15))
    seg15: path3d = path3d_line(point(rx15, ry15, rz15), point(rx0, ry0, rz0))

    full_path: path3d = make_path3d(seg0, seg1, seg2, seg3, seg4, seg5, seg6, seg7,
                                     seg8, seg9, seg10, seg11, seg12, seg13, seg14, seg15)
    emit full_path


# Cradle ring using latitude-based path
command latitude_cradle_ring(
    globe_diameter: float,
    oblateness: float,
    latitude_deg: float,
    tilt_deg: float,
    globe_center_z: float,
    beam_outer: float,
    beam_wall: float,
    angle_threshold: float
) -> solid:
    outer: region2d = box_beam_profile(beam_outer, beam_wall)
    inner: region2d = box_beam_inner(beam_outer, beam_wall)
    beam_offset: float = beam_outer / 2.0
    path: path3d = latitude_cradle_path(globe_diameter, oblateness, latitude_deg,
                                         tilt_deg, globe_center_z, beam_offset)

    result: solid = sweep_adaptive_hollow(outer, inner, path, angle_threshold)
    emit result


# ============================================================================
# Parabolic support arc generation
# ============================================================================

# Create a parabolic arc path from base to cradle
# base_angle_deg: angle on base ring (degrees)
# top_angle_deg: angle on cradle ring (degrees)
# The arc bows inward with parabolic profile, creating a waist at mid-height
#
# The path interpolates from base point to top point with:
# - Linear interpolation of XYZ position (straight line in 3D)
# - Parabolic inward offset (waist at t=0.5)
# - The waist_ratio controls the minimum radius relative to the linear path
command parabolic_support_arc(
    base_radius: float,
    top_radius: float,
    top_center_z: float,
    top_tilt_deg: float,
    waist_ratio: float,
    base_angle_deg: float,
    top_angle_deg: float
) -> path3d:
    require base_radius > 0.0
    require top_radius > 0.0
    require waist_ratio > 0.0
    require waist_ratio < 1.0

    pi: float = 3.14159265359

    # Convert angles to radians
    base_angle_rad: float = base_angle_deg * pi / 180.0
    top_angle_rad: float = top_angle_deg * pi / 180.0
    tilt_rad: float = top_tilt_deg * pi / 180.0

    cos_tilt: float = cos(tilt_rad)
    sin_tilt: float = sin(tilt_rad)

    # Base point (on flat ring at z=0)
    base_x: float = base_radius * cos(base_angle_rad)
    base_y: float = base_radius * sin(base_angle_rad)
    base_z: float = 0.0

    # Top point (on tilted cradle ring)
    # Point on tilted ring at angle theta:
    # x = r * cos(theta)
    # y = r * sin(theta) * cos(tilt)
    # z = r * sin(theta) * sin(tilt) + center_z
    top_x: float = top_radius * cos(top_angle_rad)
    top_y: float = top_radius * sin(top_angle_rad) * cos_tilt
    top_z: float = top_radius * sin(top_angle_rad) * sin_tilt + top_center_z

    # Calculate the inward offset direction at each point
    # The waist pulls the path inward (toward z-axis) by a parabolic amount
    # Maximum inward pull at t=0.5
    #
    # For each point at parameter t:
    # 1. Linear position: lerp(base, top, t)
    # 2. Linear radius at that position: lerp(base_radius, top_radius, t)
    # 3. Parabolic factor: 4*t*(1-t), peaks at 1.0 when t=0.5
    # 4. Waist offset: (1 - waist_ratio) * linear_radius * parabolic_factor
    # 5. Final radius: linear_radius - waist_offset

    # Point at t=1/8
    t1: float = 0.125
    para1: float = 4.0 * t1 * (1.0 - t1)
    lin_x1: float = base_x + (top_x - base_x) * t1
    lin_y1: float = base_y + (top_y - base_y) * t1
    lin_z1: float = base_z + (top_z - base_z) * t1
    lin_r1: float = sqrt(lin_x1 * lin_x1 + lin_y1 * lin_y1)
    waist_offset1: float = (1.0 - waist_ratio) * lin_r1 * para1
    scale1: float = (lin_r1 - waist_offset1) / lin_r1
    x1: float = lin_x1 * scale1
    y1: float = lin_y1 * scale1
    z1: float = lin_z1

    # Point at t=2/8
    t2: float = 0.25
    para2: float = 4.0 * t2 * (1.0 - t2)
    lin_x2: float = base_x + (top_x - base_x) * t2
    lin_y2: float = base_y + (top_y - base_y) * t2
    lin_z2: float = base_z + (top_z - base_z) * t2
    lin_r2: float = sqrt(lin_x2 * lin_x2 + lin_y2 * lin_y2)
    waist_offset2: float = (1.0 - waist_ratio) * lin_r2 * para2
    scale2: float = (lin_r2 - waist_offset2) / lin_r2
    x2: float = lin_x2 * scale2
    y2: float = lin_y2 * scale2
    z2: float = lin_z2

    # Point at t=3/8
    t3: float = 0.375
    para3: float = 4.0 * t3 * (1.0 - t3)
    lin_x3: float = base_x + (top_x - base_x) * t3
    lin_y3: float = base_y + (top_y - base_y) * t3
    lin_z3: float = base_z + (top_z - base_z) * t3
    lin_r3: float = sqrt(lin_x3 * lin_x3 + lin_y3 * lin_y3)
    waist_offset3: float = (1.0 - waist_ratio) * lin_r3 * para3
    scale3: float = (lin_r3 - waist_offset3) / lin_r3
    x3: float = lin_x3 * scale3
    y3: float = lin_y3 * scale3
    z3: float = lin_z3

    # Point at t=4/8 (waist - maximum inward pull)
    t4: float = 0.5
    para4: float = 4.0 * t4 * (1.0 - t4)
    lin_x4: float = base_x + (top_x - base_x) * t4
    lin_y4: float = base_y + (top_y - base_y) * t4
    lin_z4: float = base_z + (top_z - base_z) * t4
    lin_r4: float = sqrt(lin_x4 * lin_x4 + lin_y4 * lin_y4)
    waist_offset4: float = (1.0 - waist_ratio) * lin_r4 * para4
    scale4: float = (lin_r4 - waist_offset4) / lin_r4
    x4: float = lin_x4 * scale4
    y4: float = lin_y4 * scale4
    z4: float = lin_z4

    # Point at t=5/8
    t5: float = 0.625
    para5: float = 4.0 * t5 * (1.0 - t5)
    lin_x5: float = base_x + (top_x - base_x) * t5
    lin_y5: float = base_y + (top_y - base_y) * t5
    lin_z5: float = base_z + (top_z - base_z) * t5
    lin_r5: float = sqrt(lin_x5 * lin_x5 + lin_y5 * lin_y5)
    waist_offset5: float = (1.0 - waist_ratio) * lin_r5 * para5
    scale5: float = (lin_r5 - waist_offset5) / lin_r5
    x5: float = lin_x5 * scale5
    y5: float = lin_y5 * scale5
    z5: float = lin_z5

    # Point at t=6/8
    t6: float = 0.75
    para6: float = 4.0 * t6 * (1.0 - t6)
    lin_x6: float = base_x + (top_x - base_x) * t6
    lin_y6: float = base_y + (top_y - base_y) * t6
    lin_z6: float = base_z + (top_z - base_z) * t6
    lin_r6: float = sqrt(lin_x6 * lin_x6 + lin_y6 * lin_y6)
    waist_offset6: float = (1.0 - waist_ratio) * lin_r6 * para6
    scale6: float = (lin_r6 - waist_offset6) / lin_r6
    x6: float = lin_x6 * scale6
    y6: float = lin_y6 * scale6
    z6: float = lin_z6

    # Point at t=7/8
    t7: float = 0.875
    para7: float = 4.0 * t7 * (1.0 - t7)
    lin_x7: float = base_x + (top_x - base_x) * t7
    lin_y7: float = base_y + (top_y - base_y) * t7
    lin_z7: float = base_z + (top_z - base_z) * t7
    lin_r7: float = sqrt(lin_x7 * lin_x7 + lin_y7 * lin_y7)
    waist_offset7: float = (1.0 - waist_ratio) * lin_r7 * para7
    scale7: float = (lin_r7 - waist_offset7) / lin_r7
    x7: float = lin_x7 * scale7
    y7: float = lin_y7 * scale7
    z7: float = lin_z7

    # Build path segments
    # Start at base, end exactly at top (on the cradle ring)
    seg0: path3d = path3d_line(point(base_x, base_y, base_z), point(x1, y1, z1))
    seg1: path3d = path3d_line(point(x1, y1, z1), point(x2, y2, z2))
    seg2: path3d = path3d_line(point(x2, y2, z2), point(x3, y3, z3))
    seg3: path3d = path3d_line(point(x3, y3, z3), point(x4, y4, z4))
    seg4: path3d = path3d_line(point(x4, y4, z4), point(x5, y5, z5))
    seg5: path3d = path3d_line(point(x5, y5, z5), point(x6, y6, z6))
    seg6: path3d = path3d_line(point(x6, y6, z6), point(x7, y7, z7))
    seg7: path3d = path3d_line(point(x7, y7, z7), point(top_x, top_y, top_z))

    full_arc: path3d = make_path3d(seg0, seg1, seg2, seg3, seg4, seg5, seg6, seg7)
    emit full_arc


# ============================================================================
# Main stand components
# ============================================================================

# Base ring - flat on the table
command base_ring(
    radius: float,
    beam_outer: float,
    beam_wall: float,
    angle_threshold: float
) -> solid:
    outer: region2d = box_beam_profile(beam_outer, beam_wall)
    inner: region2d = box_beam_inner(beam_outer, beam_wall)
    path: path3d = ring_path(radius, 0.0)

    result: solid = sweep_adaptive_hollow(outer, inner, path, angle_threshold)
    emit result


# Cradle ring - tilted to hold the globe
command cradle_ring(
    radius: float,
    center_z: float,
    tilt_deg: float,
    beam_outer: float,
    beam_wall: float,
    angle_threshold: float
) -> solid:
    outer: region2d = box_beam_profile(beam_outer, beam_wall)
    inner: region2d = box_beam_inner(beam_outer, beam_wall)
    path: path3d = tilted_ring_path(radius, center_z, tilt_deg)

    result: solid = sweep_adaptive_hollow(outer, inner, path, angle_threshold)
    emit result


# Single support arc
command support_arc(
    base_radius: float,
    top_radius: float,
    top_center_z: float,
    top_tilt_deg: float,
    waist_ratio: float,
    base_angle_deg: float,
    top_angle_deg: float,
    beam_outer: float,
    beam_wall: float,
    angle_threshold: float
) -> solid:
    outer: region2d = box_beam_profile(beam_outer, beam_wall)
    inner: region2d = box_beam_inner(beam_outer, beam_wall)
    path: path3d = parabolic_support_arc(
        base_radius,
        top_radius, top_center_z, top_tilt_deg,
        waist_ratio,
        base_angle_deg, top_angle_deg
    )

    result: solid = sweep_adaptive_hollow(outer, inner, path, angle_threshold)
    emit result


# ============================================================================
# Latitude-based support arc
# ============================================================================

# Create a parabolic support arc from base to a point on the latitude-based cradle
# The top point is computed as a point on the tilted oblate spheroid's latitude circle
command latitude_support_arc_path(
    base_radius: float,
    globe_diameter: float,
    oblateness: float,
    latitude_deg: float,
    tilt_deg: float,
    globe_center_z: float,
    beam_offset: float,
    waist_ratio: float,
    base_angle_deg: float,
    top_longitude_deg: float
) -> path3d:
    require base_radius > 0.0
    require globe_diameter > 0.0
    require waist_ratio > 0.0
    require waist_ratio < 1.0

    pi: float = 3.14159265359

    globe_radius: float = globe_diameter / 2.0
    polar_radius: float = globe_radius * (1.0 - oblateness)

    lat_rad: float = latitude_deg * pi / 180.0
    tilt_rad: float = tilt_deg * pi / 180.0
    base_angle_rad: float = base_angle_deg * pi / 180.0
    top_lon_rad: float = top_longitude_deg * pi / 180.0

    cos_lat: float = cos(lat_rad)
    sin_lat: float = sin(lat_rad)
    cos_tilt: float = cos(tilt_rad)
    sin_tilt: float = sin(tilt_rad)

    a2: float = globe_radius * globe_radius
    c2: float = polar_radius * polar_radius

    # Base point (on flat ring at z=0)
    base_x: float = base_radius * cos(base_angle_rad)
    base_y: float = base_radius * sin(base_angle_rad)
    base_z: float = 0.0

    # Top point: on latitude circle of tilted oblate spheroid, offset by beam radius
    # First compute point on ellipsoid surface
    tx: float = globe_radius * cos_lat * cos(top_lon_rad)
    ty: float = globe_radius * cos_lat * sin(top_lon_rad)
    tz: float = polar_radius * sin_lat

    # Compute outward normal and offset
    tnx: float = tx / a2
    tny: float = ty / a2
    tnz: float = tz / c2
    tnlen: float = sqrt(tnx*tnx + tny*tny + tnz*tnz)
    tpx: float = tx + beam_offset * tnx / tnlen
    tpy: float = ty + beam_offset * tny / tnlen
    tpz: float = tz + beam_offset * tnz / tnlen

    # Apply tilt rotation around X axis and translate to globe center
    top_x: float = tpx
    top_y: float = tpy * cos_tilt - tpz * sin_tilt
    top_z: float = tpy * sin_tilt + tpz * cos_tilt + globe_center_z

    # Now generate parabolic arc from base to top with waist
    # Point at t=1/8
    t1: float = 0.125
    para1: float = 4.0 * t1 * (1.0 - t1)
    lin_x1: float = base_x + (top_x - base_x) * t1
    lin_y1: float = base_y + (top_y - base_y) * t1
    lin_z1: float = base_z + (top_z - base_z) * t1
    lin_r1: float = sqrt(lin_x1 * lin_x1 + lin_y1 * lin_y1)
    waist_offset1: float = (1.0 - waist_ratio) * lin_r1 * para1
    scale1: float = (lin_r1 - waist_offset1) / lin_r1
    x1: float = lin_x1 * scale1
    y1: float = lin_y1 * scale1
    z1: float = lin_z1

    # Point at t=2/8
    t2: float = 0.25
    para2: float = 4.0 * t2 * (1.0 - t2)
    lin_x2: float = base_x + (top_x - base_x) * t2
    lin_y2: float = base_y + (top_y - base_y) * t2
    lin_z2: float = base_z + (top_z - base_z) * t2
    lin_r2: float = sqrt(lin_x2 * lin_x2 + lin_y2 * lin_y2)
    waist_offset2: float = (1.0 - waist_ratio) * lin_r2 * para2
    scale2: float = (lin_r2 - waist_offset2) / lin_r2
    x2: float = lin_x2 * scale2
    y2: float = lin_y2 * scale2
    z2: float = lin_z2

    # Point at t=3/8
    t3: float = 0.375
    para3: float = 4.0 * t3 * (1.0 - t3)
    lin_x3: float = base_x + (top_x - base_x) * t3
    lin_y3: float = base_y + (top_y - base_y) * t3
    lin_z3: float = base_z + (top_z - base_z) * t3
    lin_r3: float = sqrt(lin_x3 * lin_x3 + lin_y3 * lin_y3)
    waist_offset3: float = (1.0 - waist_ratio) * lin_r3 * para3
    scale3: float = (lin_r3 - waist_offset3) / lin_r3
    x3: float = lin_x3 * scale3
    y3: float = lin_y3 * scale3
    z3: float = lin_z3

    # Point at t=4/8 (waist)
    t4: float = 0.5
    para4: float = 4.0 * t4 * (1.0 - t4)
    lin_x4: float = base_x + (top_x - base_x) * t4
    lin_y4: float = base_y + (top_y - base_y) * t4
    lin_z4: float = base_z + (top_z - base_z) * t4
    lin_r4: float = sqrt(lin_x4 * lin_x4 + lin_y4 * lin_y4)
    waist_offset4: float = (1.0 - waist_ratio) * lin_r4 * para4
    scale4: float = (lin_r4 - waist_offset4) / lin_r4
    x4: float = lin_x4 * scale4
    y4: float = lin_y4 * scale4
    z4: float = lin_z4

    # Point at t=5/8
    t5: float = 0.625
    para5: float = 4.0 * t5 * (1.0 - t5)
    lin_x5: float = base_x + (top_x - base_x) * t5
    lin_y5: float = base_y + (top_y - base_y) * t5
    lin_z5: float = base_z + (top_z - base_z) * t5
    lin_r5: float = sqrt(lin_x5 * lin_x5 + lin_y5 * lin_y5)
    waist_offset5: float = (1.0 - waist_ratio) * lin_r5 * para5
    scale5: float = (lin_r5 - waist_offset5) / lin_r5
    x5: float = lin_x5 * scale5
    y5: float = lin_y5 * scale5
    z5: float = lin_z5

    # Point at t=6/8
    t6: float = 0.75
    para6: float = 4.0 * t6 * (1.0 - t6)
    lin_x6: float = base_x + (top_x - base_x) * t6
    lin_y6: float = base_y + (top_y - base_y) * t6
    lin_z6: float = base_z + (top_z - base_z) * t6
    lin_r6: float = sqrt(lin_x6 * lin_x6 + lin_y6 * lin_y6)
    waist_offset6: float = (1.0 - waist_ratio) * lin_r6 * para6
    scale6: float = (lin_r6 - waist_offset6) / lin_r6
    x6: float = lin_x6 * scale6
    y6: float = lin_y6 * scale6
    z6: float = lin_z6

    # Point at t=7/8
    t7: float = 0.875
    para7: float = 4.0 * t7 * (1.0 - t7)
    lin_x7: float = base_x + (top_x - base_x) * t7
    lin_y7: float = base_y + (top_y - base_y) * t7
    lin_z7: float = base_z + (top_z - base_z) * t7
    lin_r7: float = sqrt(lin_x7 * lin_x7 + lin_y7 * lin_y7)
    waist_offset7: float = (1.0 - waist_ratio) * lin_r7 * para7
    scale7: float = (lin_r7 - waist_offset7) / lin_r7
    x7: float = lin_x7 * scale7
    y7: float = lin_y7 * scale7
    z7: float = lin_z7

    # Build path segments
    seg0: path3d = path3d_line(point(base_x, base_y, base_z), point(x1, y1, z1))
    seg1: path3d = path3d_line(point(x1, y1, z1), point(x2, y2, z2))
    seg2: path3d = path3d_line(point(x2, y2, z2), point(x3, y3, z3))
    seg3: path3d = path3d_line(point(x3, y3, z3), point(x4, y4, z4))
    seg4: path3d = path3d_line(point(x4, y4, z4), point(x5, y5, z5))
    seg5: path3d = path3d_line(point(x5, y5, z5), point(x6, y6, z6))
    seg6: path3d = path3d_line(point(x6, y6, z6), point(x7, y7, z7))
    seg7: path3d = path3d_line(point(x7, y7, z7), point(top_x, top_y, top_z))

    full_arc: path3d = make_path3d(seg0, seg1, seg2, seg3, seg4, seg5, seg6, seg7)
    emit full_arc


# Support arc solid using latitude-based path
command latitude_support_arc(
    base_radius: float,
    globe_diameter: float,
    oblateness: float,
    latitude_deg: float,
    tilt_deg: float,
    globe_center_z: float,
    waist_ratio: float,
    base_angle_deg: float,
    top_longitude_deg: float,
    beam_outer: float,
    beam_wall: float,
    angle_threshold: float
) -> solid:
    outer: region2d = box_beam_profile(beam_outer, beam_wall)
    inner: region2d = box_beam_inner(beam_outer, beam_wall)
    beam_offset: float = beam_outer / 2.0
    path: path3d = latitude_support_arc_path(
        base_radius, globe_diameter, oblateness, latitude_deg,
        tilt_deg, globe_center_z, beam_offset,
        waist_ratio, base_angle_deg, top_longitude_deg
    )

    result: solid = sweep_adaptive_hollow(outer, inner, path, angle_threshold)
    emit result


# ============================================================================
# Wrapped support arc - projects arc points onto globe surface where they intersect
# ============================================================================

# This implementation uses a smooth sigmoid blend between the original parabolic
# arc path and the projected (surface-hugging) path.
#
# For each arc sample point:
# 1. Compute original parabolic arc point a1 (raw_x, raw_y, raw_z)
# 2. Transform to globe-local coords and compute implicit distance:
#    d = check - 1, where check = (x² + y²)/a² + z²/c² for clearance ellipsoid
#    d < 0: inside, d = 0: on surface, d > 0: outside
# 3. Compute projected point a2 (raycast to surface)
# 4. Blend: A = (1-s)*a1 + s*a2 where s = 1/(1 + exp(k*d))
#    - k controls transition sharpness (k=20 gives sharp transition)
#    - When d << 0 (inside): s ≈ 1, use projected
#    - When d >> 0 (outside): s ≈ 0, use original
#    - Smooth transition near surface

command wrapped_support_arc_path(
    base_radius: float,
    globe_diameter: float,
    oblateness: float,
    latitude_deg: float,
    tilt_deg: float,
    globe_center_z: float,
    beam_offset: float,
    waist_ratio: float,
    base_angle_deg: float,
    top_longitude_deg: float
) -> path3d:
    require base_radius > 0.0
    require globe_diameter > 0.0
    require waist_ratio > 0.0
    require waist_ratio < 1.0

    pi: float = 3.14159265359

    globe_radius: float = globe_diameter / 2.0
    polar_radius: float = globe_radius * (1.0 - oblateness)

    # Clearance ellipsoid = globe + beam_offset
    a_clear: float = globe_radius + beam_offset
    c_clear: float = polar_radius + beam_offset
    a_clear2: float = a_clear * a_clear
    c_clear2: float = c_clear * c_clear

    # Sigmoid sharpness parameter (k=20 gives transition over ~0.1 units of d)
    sigmoid_k: float = 20.0

    lat_rad: float = latitude_deg * pi / 180.0
    tilt_rad: float = tilt_deg * pi / 180.0
    base_angle_rad: float = base_angle_deg * pi / 180.0
    top_lon_rad: float = top_longitude_deg * pi / 180.0

    cos_lat: float = cos(lat_rad)
    sin_lat: float = sin(lat_rad)
    cos_tilt: float = cos(tilt_rad)
    sin_tilt: float = sin(tilt_rad)

    a2: float = globe_radius * globe_radius
    c2: float = polar_radius * polar_radius

    # Base point (on flat ring at z=0) - never needs wrapping
    base_x: float = base_radius * cos(base_angle_rad)
    base_y: float = base_radius * sin(base_angle_rad)
    base_z: float = 0.0

    # Top point: on latitude circle (already on clearance surface, no wrapping needed)
    tx: float = globe_radius * cos_lat * cos(top_lon_rad)
    ty: float = globe_radius * cos_lat * sin(top_lon_rad)
    tz: float = polar_radius * sin_lat
    tnx: float = tx / a2
    tny: float = ty / a2
    tnz: float = tz / c2
    tnlen: float = sqrt(tnx*tnx + tny*tny + tnz*tnz)
    tpx: float = tx + beam_offset * tnx / tnlen
    tpy: float = ty + beam_offset * tny / tnlen
    tpz: float = tz + beam_offset * tnz / tnlen
    top_x: float = tpx
    top_y: float = tpy * cos_tilt - tpz * sin_tilt
    top_z: float = tpy * sin_tilt + tpz * cos_tilt + globe_center_z

    # Generate parabolic arc points with sigmoid-blended wrapping

    # Point at t=1/8
    t1: float = 0.125
    para1: float = 4.0 * t1 * (1.0 - t1)
    lin_x1: float = base_x + (top_x - base_x) * t1
    lin_y1: float = base_y + (top_y - base_y) * t1
    lin_z1: float = base_z + (top_z - base_z) * t1
    lin_r1: float = sqrt(lin_x1 * lin_x1 + lin_y1 * lin_y1)
    waist_offset1: float = (1.0 - waist_ratio) * lin_r1 * para1
    scale1: float = (lin_r1 - waist_offset1) / lin_r1
    raw_x1: float = lin_x1 * scale1
    raw_y1: float = lin_y1 * scale1
    raw_z1: float = lin_z1
    # Transform to globe-local: translate then inverse-tilt
    loc_x1: float = raw_x1
    loc_y1: float = raw_y1 * cos_tilt + (raw_z1 - globe_center_z) * sin_tilt
    loc_z1: float = -raw_y1 * sin_tilt + (raw_z1 - globe_center_z) * cos_tilt
    # Compute implicit distance d = check - 1 (negative inside, positive outside)
    check1: float = (loc_x1*loc_x1 + loc_y1*loc_y1) / a_clear2 + loc_z1*loc_z1 / c_clear2
    d1: float = check1 - 1.0
    # Sigmoid blend factor: s = 1/(1 + exp(k*d))
    s1: float = 1.0 / (1.0 + exp(sigmoid_k * d1))
    # Compute projected point (raycast to surface)
    proj_t1: float = 1.0 / sqrt(check1)
    proj_loc_x1: float = loc_x1 * proj_t1
    proj_loc_y1: float = loc_y1 * proj_t1
    proj_loc_z1: float = loc_z1 * proj_t1
    # Transform projected point back to world coords
    proj_world_x1: float = proj_loc_x1
    proj_world_y1: float = proj_loc_y1 * cos_tilt - proj_loc_z1 * sin_tilt
    proj_world_z1: float = proj_loc_y1 * sin_tilt + proj_loc_z1 * cos_tilt + globe_center_z
    # Blend: final = (1-s)*original + s*projected
    x1: float = (1.0 - s1) * raw_x1 + s1 * proj_world_x1
    y1: float = (1.0 - s1) * raw_y1 + s1 * proj_world_y1
    z1: float = (1.0 - s1) * raw_z1 + s1 * proj_world_z1

    # Point at t=2/8
    t2: float = 0.25
    para2: float = 4.0 * t2 * (1.0 - t2)
    lin_x2: float = base_x + (top_x - base_x) * t2
    lin_y2: float = base_y + (top_y - base_y) * t2
    lin_z2: float = base_z + (top_z - base_z) * t2
    lin_r2: float = sqrt(lin_x2 * lin_x2 + lin_y2 * lin_y2)
    waist_offset2: float = (1.0 - waist_ratio) * lin_r2 * para2
    scale2: float = (lin_r2 - waist_offset2) / lin_r2
    raw_x2: float = lin_x2 * scale2
    raw_y2: float = lin_y2 * scale2
    raw_z2: float = lin_z2
    loc_x2: float = raw_x2
    loc_y2: float = raw_y2 * cos_tilt + (raw_z2 - globe_center_z) * sin_tilt
    loc_z2: float = -raw_y2 * sin_tilt + (raw_z2 - globe_center_z) * cos_tilt
    check2: float = (loc_x2*loc_x2 + loc_y2*loc_y2) / a_clear2 + loc_z2*loc_z2 / c_clear2
    d2: float = check2 - 1.0
    s2: float = 1.0 / (1.0 + exp(sigmoid_k * d2))
    proj_t2: float = 1.0 / sqrt(check2)
    proj_loc_x2: float = loc_x2 * proj_t2
    proj_loc_y2: float = loc_y2 * proj_t2
    proj_loc_z2: float = loc_z2 * proj_t2
    proj_world_x2: float = proj_loc_x2
    proj_world_y2: float = proj_loc_y2 * cos_tilt - proj_loc_z2 * sin_tilt
    proj_world_z2: float = proj_loc_y2 * sin_tilt + proj_loc_z2 * cos_tilt + globe_center_z
    x2: float = (1.0 - s2) * raw_x2 + s2 * proj_world_x2
    y2: float = (1.0 - s2) * raw_y2 + s2 * proj_world_y2
    z2: float = (1.0 - s2) * raw_z2 + s2 * proj_world_z2

    # Point at t=3/8
    t3: float = 0.375
    para3: float = 4.0 * t3 * (1.0 - t3)
    lin_x3: float = base_x + (top_x - base_x) * t3
    lin_y3: float = base_y + (top_y - base_y) * t3
    lin_z3: float = base_z + (top_z - base_z) * t3
    lin_r3: float = sqrt(lin_x3 * lin_x3 + lin_y3 * lin_y3)
    waist_offset3: float = (1.0 - waist_ratio) * lin_r3 * para3
    scale3: float = (lin_r3 - waist_offset3) / lin_r3
    raw_x3: float = lin_x3 * scale3
    raw_y3: float = lin_y3 * scale3
    raw_z3: float = lin_z3
    loc_x3: float = raw_x3
    loc_y3: float = raw_y3 * cos_tilt + (raw_z3 - globe_center_z) * sin_tilt
    loc_z3: float = -raw_y3 * sin_tilt + (raw_z3 - globe_center_z) * cos_tilt
    check3: float = (loc_x3*loc_x3 + loc_y3*loc_y3) / a_clear2 + loc_z3*loc_z3 / c_clear2
    d3: float = check3 - 1.0
    s3: float = 1.0 / (1.0 + exp(sigmoid_k * d3))
    proj_t3: float = 1.0 / sqrt(check3)
    proj_loc_x3: float = loc_x3 * proj_t3
    proj_loc_y3: float = loc_y3 * proj_t3
    proj_loc_z3: float = loc_z3 * proj_t3
    proj_world_x3: float = proj_loc_x3
    proj_world_y3: float = proj_loc_y3 * cos_tilt - proj_loc_z3 * sin_tilt
    proj_world_z3: float = proj_loc_y3 * sin_tilt + proj_loc_z3 * cos_tilt + globe_center_z
    x3: float = (1.0 - s3) * raw_x3 + s3 * proj_world_x3
    y3: float = (1.0 - s3) * raw_y3 + s3 * proj_world_y3
    z3: float = (1.0 - s3) * raw_z3 + s3 * proj_world_z3

    # Point at t=4/8 (waist)
    t4: float = 0.5
    para4: float = 4.0 * t4 * (1.0 - t4)
    lin_x4: float = base_x + (top_x - base_x) * t4
    lin_y4: float = base_y + (top_y - base_y) * t4
    lin_z4: float = base_z + (top_z - base_z) * t4
    lin_r4: float = sqrt(lin_x4 * lin_x4 + lin_y4 * lin_y4)
    waist_offset4: float = (1.0 - waist_ratio) * lin_r4 * para4
    scale4: float = (lin_r4 - waist_offset4) / lin_r4
    raw_x4: float = lin_x4 * scale4
    raw_y4: float = lin_y4 * scale4
    raw_z4: float = lin_z4
    loc_x4: float = raw_x4
    loc_y4: float = raw_y4 * cos_tilt + (raw_z4 - globe_center_z) * sin_tilt
    loc_z4: float = -raw_y4 * sin_tilt + (raw_z4 - globe_center_z) * cos_tilt
    check4: float = (loc_x4*loc_x4 + loc_y4*loc_y4) / a_clear2 + loc_z4*loc_z4 / c_clear2
    d4: float = check4 - 1.0
    s4: float = 1.0 / (1.0 + exp(sigmoid_k * d4))
    proj_t4: float = 1.0 / sqrt(check4)
    proj_loc_x4: float = loc_x4 * proj_t4
    proj_loc_y4: float = loc_y4 * proj_t4
    proj_loc_z4: float = loc_z4 * proj_t4
    proj_world_x4: float = proj_loc_x4
    proj_world_y4: float = proj_loc_y4 * cos_tilt - proj_loc_z4 * sin_tilt
    proj_world_z4: float = proj_loc_y4 * sin_tilt + proj_loc_z4 * cos_tilt + globe_center_z
    x4: float = (1.0 - s4) * raw_x4 + s4 * proj_world_x4
    y4: float = (1.0 - s4) * raw_y4 + s4 * proj_world_y4
    z4: float = (1.0 - s4) * raw_z4 + s4 * proj_world_z4

    # Point at t=5/8
    t5: float = 0.625
    para5: float = 4.0 * t5 * (1.0 - t5)
    lin_x5: float = base_x + (top_x - base_x) * t5
    lin_y5: float = base_y + (top_y - base_y) * t5
    lin_z5: float = base_z + (top_z - base_z) * t5
    lin_r5: float = sqrt(lin_x5 * lin_x5 + lin_y5 * lin_y5)
    waist_offset5: float = (1.0 - waist_ratio) * lin_r5 * para5
    scale5: float = (lin_r5 - waist_offset5) / lin_r5
    raw_x5: float = lin_x5 * scale5
    raw_y5: float = lin_y5 * scale5
    raw_z5: float = lin_z5
    loc_x5: float = raw_x5
    loc_y5: float = raw_y5 * cos_tilt + (raw_z5 - globe_center_z) * sin_tilt
    loc_z5: float = -raw_y5 * sin_tilt + (raw_z5 - globe_center_z) * cos_tilt
    check5: float = (loc_x5*loc_x5 + loc_y5*loc_y5) / a_clear2 + loc_z5*loc_z5 / c_clear2
    d5: float = check5 - 1.0
    s5: float = 1.0 / (1.0 + exp(sigmoid_k * d5))
    proj_t5: float = 1.0 / sqrt(check5)
    proj_loc_x5: float = loc_x5 * proj_t5
    proj_loc_y5: float = loc_y5 * proj_t5
    proj_loc_z5: float = loc_z5 * proj_t5
    proj_world_x5: float = proj_loc_x5
    proj_world_y5: float = proj_loc_y5 * cos_tilt - proj_loc_z5 * sin_tilt
    proj_world_z5: float = proj_loc_y5 * sin_tilt + proj_loc_z5 * cos_tilt + globe_center_z
    x5: float = (1.0 - s5) * raw_x5 + s5 * proj_world_x5
    y5: float = (1.0 - s5) * raw_y5 + s5 * proj_world_y5
    z5: float = (1.0 - s5) * raw_z5 + s5 * proj_world_z5

    # Point at t=6/8
    t6: float = 0.75
    para6: float = 4.0 * t6 * (1.0 - t6)
    lin_x6: float = base_x + (top_x - base_x) * t6
    lin_y6: float = base_y + (top_y - base_y) * t6
    lin_z6: float = base_z + (top_z - base_z) * t6
    lin_r6: float = sqrt(lin_x6 * lin_x6 + lin_y6 * lin_y6)
    waist_offset6: float = (1.0 - waist_ratio) * lin_r6 * para6
    scale6: float = (lin_r6 - waist_offset6) / lin_r6
    raw_x6: float = lin_x6 * scale6
    raw_y6: float = lin_y6 * scale6
    raw_z6: float = lin_z6
    loc_x6: float = raw_x6
    loc_y6: float = raw_y6 * cos_tilt + (raw_z6 - globe_center_z) * sin_tilt
    loc_z6: float = -raw_y6 * sin_tilt + (raw_z6 - globe_center_z) * cos_tilt
    check6: float = (loc_x6*loc_x6 + loc_y6*loc_y6) / a_clear2 + loc_z6*loc_z6 / c_clear2
    d6: float = check6 - 1.0
    s6: float = 1.0 / (1.0 + exp(sigmoid_k * d6))
    proj_t6: float = 1.0 / sqrt(check6)
    proj_loc_x6: float = loc_x6 * proj_t6
    proj_loc_y6: float = loc_y6 * proj_t6
    proj_loc_z6: float = loc_z6 * proj_t6
    proj_world_x6: float = proj_loc_x6
    proj_world_y6: float = proj_loc_y6 * cos_tilt - proj_loc_z6 * sin_tilt
    proj_world_z6: float = proj_loc_y6 * sin_tilt + proj_loc_z6 * cos_tilt + globe_center_z
    x6: float = (1.0 - s6) * raw_x6 + s6 * proj_world_x6
    y6: float = (1.0 - s6) * raw_y6 + s6 * proj_world_y6
    z6: float = (1.0 - s6) * raw_z6 + s6 * proj_world_z6

    # Point at t=7/8
    t7: float = 0.875
    para7: float = 4.0 * t7 * (1.0 - t7)
    lin_x7: float = base_x + (top_x - base_x) * t7
    lin_y7: float = base_y + (top_y - base_y) * t7
    lin_z7: float = base_z + (top_z - base_z) * t7
    lin_r7: float = sqrt(lin_x7 * lin_x7 + lin_y7 * lin_y7)
    waist_offset7: float = (1.0 - waist_ratio) * lin_r7 * para7
    scale7: float = (lin_r7 - waist_offset7) / lin_r7
    raw_x7: float = lin_x7 * scale7
    raw_y7: float = lin_y7 * scale7
    raw_z7: float = lin_z7
    loc_x7: float = raw_x7
    loc_y7: float = raw_y7 * cos_tilt + (raw_z7 - globe_center_z) * sin_tilt
    loc_z7: float = -raw_y7 * sin_tilt + (raw_z7 - globe_center_z) * cos_tilt
    check7: float = (loc_x7*loc_x7 + loc_y7*loc_y7) / a_clear2 + loc_z7*loc_z7 / c_clear2
    d7: float = check7 - 1.0
    s7: float = 1.0 / (1.0 + exp(sigmoid_k * d7))
    proj_t7: float = 1.0 / sqrt(check7)
    proj_loc_x7: float = loc_x7 * proj_t7
    proj_loc_y7: float = loc_y7 * proj_t7
    proj_loc_z7: float = loc_z7 * proj_t7
    proj_world_x7: float = proj_loc_x7
    proj_world_y7: float = proj_loc_y7 * cos_tilt - proj_loc_z7 * sin_tilt
    proj_world_z7: float = proj_loc_y7 * sin_tilt + proj_loc_z7 * cos_tilt + globe_center_z
    x7: float = (1.0 - s7) * raw_x7 + s7 * proj_world_x7
    y7: float = (1.0 - s7) * raw_y7 + s7 * proj_world_y7
    z7: float = (1.0 - s7) * raw_z7 + s7 * proj_world_z7

    # Build path segments
    seg0: path3d = path3d_line(point(base_x, base_y, base_z), point(x1, y1, z1))
    seg1: path3d = path3d_line(point(x1, y1, z1), point(x2, y2, z2))
    seg2: path3d = path3d_line(point(x2, y2, z2), point(x3, y3, z3))
    seg3: path3d = path3d_line(point(x3, y3, z3), point(x4, y4, z4))
    seg4: path3d = path3d_line(point(x4, y4, z4), point(x5, y5, z5))
    seg5: path3d = path3d_line(point(x5, y5, z5), point(x6, y6, z6))
    seg6: path3d = path3d_line(point(x6, y6, z6), point(x7, y7, z7))
    seg7: path3d = path3d_line(point(x7, y7, z7), point(top_x, top_y, top_z))

    full_arc: path3d = make_path3d(seg0, seg1, seg2, seg3, seg4, seg5, seg6, seg7)
    emit full_arc


# Wrapped support arc solid
command wrapped_support_arc(
    base_radius: float,
    globe_diameter: float,
    oblateness: float,
    latitude_deg: float,
    tilt_deg: float,
    globe_center_z: float,
    waist_ratio: float,
    base_angle_deg: float,
    top_longitude_deg: float,
    beam_outer: float,
    beam_wall: float,
    angle_threshold: float
) -> solid:
    outer: region2d = box_beam_profile(beam_outer, beam_wall)
    inner: region2d = box_beam_inner(beam_outer, beam_wall)
    beam_offset: float = beam_outer / 2.0
    path: path3d = wrapped_support_arc_path(
        base_radius, globe_diameter, oblateness, latitude_deg,
        tilt_deg, globe_center_z, beam_offset,
        waist_ratio, base_angle_deg, top_longitude_deg
    )

    result: solid = sweep_adaptive_hollow(outer, inner, path, angle_threshold)
    emit result


# ============================================================================
# Top-level assembly commands
# ============================================================================

# Complete globe stand
command GLOBE_STAND(
    globe_diameter: float = 304.8,      # 30.48cm = 12 inches
    stand_height: float = 300.0,        # 30cm tall
    base_diameter: float = 400.0,       # 40cm base
    tilt_angle: float = 25.2,           # Mars axial tilt
    waist_ratio: float = 0.7,           # 70% waist
    twist_deg: float = 60.0,            # twist from base to cradle (degrees)
    beam_outer: float = 10.0,           # 10mm beam
    beam_wall: float = 2.0,             # 2mm wall
    angle_threshold: float = 5.0
) -> solid:
    # Calculate derived dimensions
    base_radius: float = base_diameter / 2.0
    globe_radius: float = globe_diameter / 2.0

    # Cradle ring radius - slightly smaller than globe to cradle it
    cradle_radius: float = globe_radius * 0.85

    # Build base ring
    base: solid = base_ring(base_radius, beam_outer, beam_wall, angle_threshold)

    # Build cradle ring at the top
    cradle: solid = cradle_ring(cradle_radius, stand_height, tilt_angle, beam_outer, beam_wall, angle_threshold)

    # Build three support arcs, offset by 120 degrees
    # Each arc twists by twist_deg from base to top
    arc1: solid = support_arc(
        base_radius, cradle_radius, stand_height, tilt_angle, waist_ratio,
        0.0, twist_deg,
        beam_outer, beam_wall, angle_threshold
    )

    arc2: solid = support_arc(
        base_radius, cradle_radius, stand_height, tilt_angle, waist_ratio,
        120.0, 120.0 + twist_deg,
        beam_outer, beam_wall, angle_threshold
    )

    arc3: solid = support_arc(
        base_radius, cradle_radius, stand_height, tilt_angle, waist_ratio,
        240.0, 240.0 + twist_deg,
        beam_outer, beam_wall, angle_threshold
    )

    # Union all parts
    stand: solid = union(base, cradle)
    stand2: solid = union(stand, arc1)
    stand3: solid = union(stand2, arc2)
    result: solid = union(stand3, arc3)

    emit result


# Individual components for testing

command BASE_RING(
    base_diameter: float = 400.0,
    beam_outer: float = 10.0,
    beam_wall: float = 2.0,
    angle_threshold: float = 5.0
) -> solid:
    base_radius: float = base_diameter / 2.0
    emit base_ring(base_radius, beam_outer, beam_wall, angle_threshold)


command CRADLE_RING(
    globe_diameter: float = 304.8,
    stand_height: float = 300.0,
    tilt_angle: float = 25.2,
    beam_outer: float = 10.0,
    beam_wall: float = 2.0,
    angle_threshold: float = 5.0
) -> solid:
    globe_radius: float = globe_diameter / 2.0
    cradle_radius: float = globe_radius * 0.85
    emit cradle_ring(cradle_radius, stand_height, tilt_angle, beam_outer, beam_wall, angle_threshold)


command SUPPORT_ARC(
    base_diameter: float = 400.0,
    globe_diameter: float = 304.8,
    stand_height: float = 300.0,
    tilt_angle: float = 25.2,
    waist_ratio: float = 0.7,
    base_angle: float = 0.0,
    top_angle: float = 60.0,
    beam_outer: float = 10.0,
    beam_wall: float = 2.0,
    angle_threshold: float = 5.0
) -> solid:
    base_radius: float = base_diameter / 2.0
    globe_radius: float = globe_diameter / 2.0
    cradle_radius: float = globe_radius * 0.85

    emit support_arc(
        base_radius, cradle_radius, stand_height, tilt_angle, waist_ratio,
        base_angle, top_angle,
        beam_outer, beam_wall, angle_threshold
    )


# Test: just base and one arc union
command BASE_AND_ARC(
    globe_diameter: float = 304.8,
    stand_height: float = 300.0,
    base_diameter: float = 400.0,
    tilt_angle: float = 25.2,
    waist_ratio: float = 0.7,
    beam_outer: float = 10.0,
    beam_wall: float = 2.0,
    angle_threshold: float = 5.0
) -> solid:
    base_radius: float = base_diameter / 2.0
    globe_radius: float = globe_diameter / 2.0
    cradle_radius: float = globe_radius * 0.85

    base: solid = base_ring(base_radius, beam_outer, beam_wall, angle_threshold)
    arc1: solid = support_arc(
        base_radius, cradle_radius, stand_height, tilt_angle, waist_ratio,
        0.0, 60.0,
        beam_outer, beam_wall, angle_threshold
    )

    result: solid = union(base, arc1)
    emit result


# Test: just base and cradle union
command BASE_AND_CRADLE(
    globe_diameter: float = 304.8,
    stand_height: float = 300.0,
    base_diameter: float = 400.0,
    tilt_angle: float = 25.2,
    beam_outer: float = 10.0,
    beam_wall: float = 2.0,
    angle_threshold: float = 5.0
) -> solid:
    base_radius: float = base_diameter / 2.0
    globe_radius: float = globe_diameter / 2.0
    cradle_radius: float = globe_radius * 0.85

    base: solid = base_ring(base_radius, beam_outer, beam_wall, angle_threshold)
    cradle: solid = cradle_ring(cradle_radius, stand_height, tilt_angle, beam_outer, beam_wall, angle_threshold)

    result: solid = union(base, cradle)
    emit result


# Assembly using variadic union
command GLOBE_STAND_V2(
    globe_diameter: float = 304.8,
    stand_height: float = 300.0,
    base_diameter: float = 400.0,
    tilt_angle: float = 25.2,
    waist_ratio: float = 0.7,
    beam_outer: float = 10.0,
    beam_wall: float = 2.0,
    angle_threshold: float = 5.0
) -> solid:
    base_radius: float = base_diameter / 2.0
    globe_radius: float = globe_diameter / 2.0
    cradle_radius: float = globe_radius * 0.85

    base: solid = base_ring(base_radius, beam_outer, beam_wall, angle_threshold)
    cradle: solid = cradle_ring(cradle_radius, stand_height, tilt_angle, beam_outer, beam_wall, angle_threshold)

    arc1: solid = support_arc(
        base_radius, cradle_radius, stand_height, tilt_angle, waist_ratio,
        0.0, 60.0,
        beam_outer, beam_wall, angle_threshold
    )

    arc2: solid = support_arc(
        base_radius, cradle_radius, stand_height, tilt_angle, waist_ratio,
        120.0, 180.0,
        beam_outer, beam_wall, angle_threshold
    )

    arc3: solid = support_arc(
        base_radius, cradle_radius, stand_height, tilt_angle, waist_ratio,
        240.0, 300.0,
        beam_outer, beam_wall, angle_threshold
    )

    # Try variadic union with all parts
    result: solid = union(base, cradle, arc1, arc2, arc3)
    emit result


# Test paths for debugging
command BASE_RING_PATH(radius: float = 200.0) -> path3d:
    emit ring_path(radius, 0.0)


command TILTED_RING_PATH(
    radius: float = 129.5,
    center_z: float = 300.0,
    tilt_deg: float = 25.2
) -> path3d:
    emit tilted_ring_path(radius, center_z, tilt_deg)


command SUPPORT_ARC_PATH(
    base_radius: float = 200.0,
    top_radius: float = 129.5,
    top_z: float = 300.0,
    tilt_deg: float = 25.2,
    waist_ratio: float = 0.7,
    base_angle: float = 0.0,
    top_angle: float = 60.0
) -> path3d:
    emit parabolic_support_arc(
        base_radius,
        top_radius, top_z, tilt_deg,
        waist_ratio,
        base_angle, top_angle
    )


# ============================================================================
# Mars Globe
# ============================================================================

# Create a Mars globe - an oblate spheroid with Mars' oblateness
# positioned and tilted to rest ON TOP of the cradle ring
#
# The globe center must be positioned so the tilted oblate spheroid is tangent
# to the OUTER top edge of the tilted cradle ring (globe resting on top).
#
# Geometry:
# - Globe is tilted at tilt_angle around X axis (same as cradle)
# - Cradle ring outer top edge at theta=90°:
#   x=0, y = outer_r * cos(tilt), z = stand_height + outer_r * sin(tilt) + beam_outer/2
#   where outer_r = cradle_radius + beam_outer/2
# - Globe surface in tilted frame: y'²/a² + z'²/c² = 1
#   where a = equatorial radius, c = polar radius
# - Globe center is ABOVE the contact point (resting on top)
command mars_globe(
    globe_diameter: float,
    stand_height: float,
    tilt_angle: float,
    beam_outer: float,
    cradle_radius: float,
    oblateness: float = 0.00648
) -> solid:
    # Mars' geometric oblateness is approximately 0.00648
    # This means polar radius = equatorial radius * (1 - 0.00648)

    pi: float = 3.14159265359
    tilt_rad: float = tilt_angle * pi / 180.0
    cos_tilt: float = cos(tilt_rad)
    sin_tilt: float = sin(tilt_rad)

    globe_radius: float = globe_diameter / 2.0
    polar_radius: float = globe_radius * (1.0 - oblateness)

    # OUTER edge of cradle ring (globe rests ON TOP)
    outer_ring_radius: float = cradle_radius + beam_outer / 2.0

    # Contact point on outer top edge of tilted ring (at theta=90°, highest point)
    # World coords: x=0, y = outer_r * cos(tilt), z = stand_height + outer_r * sin(tilt) + beam_outer/2
    k1: float = outer_ring_radius * cos_tilt
    k2: float = stand_height + outer_ring_radius * sin_tilt + beam_outer / 2.0

    # For the tilted globe surface, in the rotated frame the ellipsoid equation is:
    # y'²/a² + z'²/c² = 1
    #
    # Transform contact point to globe-local rotated frame:
    # Let d = z0 - k2 (globe center ABOVE contact point, so d > 0)
    # y' = k1*cos_tilt - d*sin_tilt
    # z' = -k1*sin_tilt - d*cos_tilt
    #
    # Ellipsoid constraint: y'²/a² + z'²/c² = 1
    # This gives quadratic: A*d² + B*d + C = 0
    # where:
    # A = sin²/a² + cos²/c²
    # B = 2*k1*sincos*(1/c² - 1/a²)
    # C = k1²*cos²/a² + k1²*sin²/c² - 1

    a: float = globe_radius
    c: float = polar_radius
    a2: float = a * a
    c2: float = c * c

    sin2: float = sin_tilt * sin_tilt
    cos2: float = cos_tilt * cos_tilt
    sincos: float = sin_tilt * cos_tilt

    A: float = sin2 / a2 + cos2 / c2
    B: float = 2.0 * k1 * sincos * (1.0 / c2 - 1.0 / a2)
    C: float = k1 * k1 * cos2 / a2 + k1 * k1 * sin2 / c2 - 1.0

    # Solve quadratic: d = (-B ± sqrt(B² - 4AC)) / 2A
    # We want the positive d (globe center ABOVE contact point)
    discriminant: float = B * B - 4.0 * A * C
    sqrt_disc: float = sqrt(discriminant)

    d1: float = (-B + sqrt_disc) / (2.0 * A)
    d2: float = (-B - sqrt_disc) / (2.0 * A)

    # Take the positive solution (d > 0 means globe center above contact)
    d_positive: float = d1

    # Globe center Z = contact Z + d
    globe_center_z: float = k2 + d_positive

    # Create oblate spheroid at origin with polar axis along Z
    globe: solid = oblate_spheroid(globe_diameter, oblateness)

    # First rotate around X axis by tilt_angle degrees
    tilted_globe: solid = rotate(globe, tilt_angle, 0.0, 0.0)

    # Then translate to the calculated center position
    positioned_globe: solid = translate(tilted_globe, 0.0, 0.0, globe_center_z)

    emit positioned_globe


# Standalone Mars globe (just the globe, for testing)
command MARS_GLOBE(
    globe_diameter: float = 304.8,
    stand_height: float = 300.0,
    tilt_angle: float = 25.2,
    beam_outer: float = 10.0,
    oblateness: float = 0.00648
) -> solid:
    globe_radius: float = globe_diameter / 2.0
    cradle_radius: float = globe_radius * 0.85
    emit mars_globe(globe_diameter, stand_height, tilt_angle, beam_outer, cradle_radius, oblateness)


# Complete globe stand with Mars globe - exports as separate solids (no union with globe)
# Use this to visually check clearance between globe and stand
command GLOBE_STAND_WITH_MARS(
    globe_diameter: float = 304.8,      # 30.48cm = 12 inches
    stand_height: float = 300.0,        # 30cm tall
    base_diameter: float = 400.0,       # 40cm base
    tilt_angle: float = 25.2,           # Mars axial tilt
    waist_ratio: float = 0.7,           # 70% waist
    twist_deg: float = 60.0,            # twist from base to cradle (degrees)
    beam_outer: float = 10.0,           # 10mm beam
    beam_wall: float = 2.0,             # 2mm wall
    angle_threshold: float = 5.0,
    mars_oblateness: float = 0.00648    # Mars geometric oblateness
) -> solid:
    # Calculate derived dimensions
    base_radius: float = base_diameter / 2.0
    globe_radius: float = globe_diameter / 2.0

    # Cradle ring radius - slightly smaller than globe to cradle it
    cradle_radius: float = globe_radius * 0.85

    # Build base ring
    base: solid = base_ring(base_radius, beam_outer, beam_wall, angle_threshold)

    # Build cradle ring at the top
    cradle: solid = cradle_ring(cradle_radius, stand_height, tilt_angle, beam_outer, beam_wall, angle_threshold)

    # Build three support arcs, offset by 120 degrees
    # Each arc twists by twist_deg from base to top
    arc1: solid = support_arc(
        base_radius, cradle_radius, stand_height, tilt_angle, waist_ratio,
        0.0, twist_deg,
        beam_outer, beam_wall, angle_threshold
    )

    arc2: solid = support_arc(
        base_radius, cradle_radius, stand_height, tilt_angle, waist_ratio,
        120.0, 120.0 + twist_deg,
        beam_outer, beam_wall, angle_threshold
    )

    arc3: solid = support_arc(
        base_radius, cradle_radius, stand_height, tilt_angle, waist_ratio,
        240.0, 240.0 + twist_deg,
        beam_outer, beam_wall, angle_threshold
    )

    # Build Mars globe with correct positioning
    mars: solid = mars_globe(globe_diameter, stand_height, tilt_angle, beam_outer, cradle_radius, mars_oblateness)

    # Union stand parts together
    stand: solid = union(base, cradle, arc1, arc2, arc3)

    # Combine stand and globe as separate solids (no boolean merge)
    # Using compound() keeps them as separate bodies for visual inspection
    result: solid = compound(stand, mars)

    emit result


# Debug: Export components separately (no union) for intersection testing
command GLOBE_STAND_PARTS(
    globe_diameter: float = 304.8,
    stand_height: float = 300.0,
    base_diameter: float = 400.0,
    tilt_angle: float = 25.2,
    waist_ratio: float = 0.7,
    twist_deg: float = 60.0,
    beam_outer: float = 10.0,
    beam_wall: float = 2.0,
    angle_threshold: float = 5.0,
    mars_oblateness: float = 0.00648
) -> solid:
    # This creates a compound of separate solids for intersection testing
    base_radius: float = base_diameter / 2.0
    globe_radius: float = globe_diameter / 2.0
    cradle_radius: float = globe_radius * 0.85

    base: solid = base_ring(base_radius, beam_outer, beam_wall, angle_threshold)
    cradle: solid = cradle_ring(cradle_radius, stand_height, tilt_angle, beam_outer, beam_wall, angle_threshold)

    arc1: solid = support_arc(
        base_radius, cradle_radius, stand_height, tilt_angle, waist_ratio,
        0.0, twist_deg,
        beam_outer, beam_wall, angle_threshold
    )
    arc2: solid = support_arc(
        base_radius, cradle_radius, stand_height, tilt_angle, waist_ratio,
        120.0, 120.0 + twist_deg,
        beam_outer, beam_wall, angle_threshold
    )
    arc3: solid = support_arc(
        base_radius, cradle_radius, stand_height, tilt_angle, waist_ratio,
        240.0, 240.0 + twist_deg,
        beam_outer, beam_wall, angle_threshold
    )

    mars: solid = mars_globe(globe_diameter, stand_height, tilt_angle, mars_oblateness)

    # Use variadic union which may preserve separate solids in compound
    result: solid = union(base, cradle, arc1, arc2, arc3, mars)
    emit result


# ============================================================================
# CENTERED GLOBE STAND - Globe CG centered over base
# ============================================================================

# This version uses a latitude-based cradle ring that follows the -35° parallel
# on the tilted oblate globe. The globe center is directly above the base center,
# providing better stability than the offset design.
#
# Key differences from GLOBE_STAND_WITH_MARS:
# - Globe is centered at (0, 0, globe_center_z) where globe_center_z is calculated
#   so the lowest point of the cradle ring is at a reasonable height
# - Cradle ring follows the -35° latitude parallel on the globe surface
# - Support arcs connect to specific longitudes on the latitude circle

# Simple Mars globe for centered design - just needs globe_center_z
command centered_mars_globe(
    globe_diameter: float,
    globe_center_z: float,
    tilt_angle: float,
    oblateness: float = 0.00648
) -> solid:
    # Create oblate spheroid at origin with polar axis along Z
    globe: solid = oblate_spheroid(globe_diameter, oblateness)

    # Rotate around X axis by tilt_angle degrees
    tilted_globe: solid = rotate(globe, tilt_angle, 0.0, 0.0)

    # Translate to the globe center position
    positioned_globe: solid = translate(tilted_globe, 0.0, 0.0, globe_center_z+5.0)

    emit positioned_globe


# Complete globe stand with centered globe
# The cradle follows the -35° latitude parallel on the tilted oblate spheroid
command CENTERED_GLOBE_STAND(
    globe_diameter: float = 304.8,      # 30.48cm = 12 inches
    base_diameter: float = 400.0,       # 40cm base
    tilt_angle: float = 25.2,           # Mars axial tilt
    cradle_latitude: float = -35.0,     # latitude for cradle ring (-35° = 35° south)
    waist_ratio: float = 0.7,           # 70% waist for support arcs
    twist_deg: float = 60.0,            # twist from base to cradle (degrees)
    beam_outer: float = 10.0,           # 10mm beam
    beam_wall: float = 2.0,             # 2mm wall
    angle_threshold: float = 5.0,
    mars_oblateness: float = 0.00648
) -> solid:
    pi: float = 3.14159265359

    base_radius: float = base_diameter / 2.0
    globe_radius: float = globe_diameter / 2.0
    polar_radius: float = globe_radius * (1.0 - mars_oblateness)
    beam_offset: float = beam_outer / 2.0

    # Calculate globe center height
    # We want the lowest point of the cradle ring (at longitude 270°, after tilt)
    # to be at a reasonable height above the base
    #
    # First, compute the lowest Z of the latitude circle (before globe translation)
    lat_rad: float = cradle_latitude * pi / 180.0
    tilt_rad: float = tilt_angle * pi / 180.0

    cos_lat: float = cos(lat_rad)
    sin_lat: float = sin(lat_rad)
    cos_tilt: float = cos(tilt_rad)
    sin_tilt: float = sin(tilt_rad)

    a2: float = globe_radius * globe_radius
    c2: float = polar_radius * polar_radius

    # Point at longitude 270° (lowest after tilt) on the latitude circle
    # Before tilt: x=0, y=-r*cos(lat), z=c*sin(lat) where r = a*cos(lat)
    x270: float = 0.0
    y270: float = -globe_radius * cos_lat
    z270: float = polar_radius * sin_lat

    # Compute normal and offset
    nx270: float = x270 / a2
    ny270: float = y270 / a2
    nz270: float = z270 / c2
    nlen270: float = sqrt(nx270*nx270 + ny270*ny270 + nz270*nz270)
    px270: float = x270 + beam_offset * nx270 / nlen270
    py270: float = y270 + beam_offset * ny270 / nlen270
    pz270: float = z270 + beam_offset * nz270 / nlen270

    # After tilt rotation around X
    rz270: float = py270 * sin_tilt + pz270 * cos_tilt

    # Set globe center so that the lowest cradle ring centerline point is at
    # a reasonable height (say, 300mm) - this can be adjusted
    target_cradle_low_z: float = 300.0
    globe_center_z: float = target_cradle_low_z - rz270

    # Build base ring
    base: solid = base_ring(base_radius, beam_outer, beam_wall, angle_threshold)

    # Build latitude-based cradle ring
    cradle: solid = latitude_cradle_ring(
        globe_diameter, mars_oblateness, cradle_latitude,
        tilt_angle, globe_center_z,
        beam_outer, beam_wall, angle_threshold
    )

    # Build three support arcs with globe wrapping
    # Each arc starts at base_angle on the base ring and connects to
    # base_angle + twist_deg longitude on the cradle latitude circle
    # The wrapped_support_arc projects any points inside the globe onto its surface
    arc1: solid = wrapped_support_arc(
        base_radius, globe_diameter, mars_oblateness, cradle_latitude,
        tilt_angle, globe_center_z, waist_ratio,
        0.0, twist_deg,
        beam_outer, beam_wall, angle_threshold
    )

    arc2: solid = wrapped_support_arc(
        base_radius, globe_diameter, mars_oblateness, cradle_latitude,
        tilt_angle, globe_center_z, waist_ratio,
        120.0, 120.0 + twist_deg,
        beam_outer, beam_wall, angle_threshold
    )

    arc3: solid = wrapped_support_arc(
        base_radius, globe_diameter, mars_oblateness, cradle_latitude,
        tilt_angle, globe_center_z, waist_ratio,
        240.0, 240.0 + twist_deg,
        beam_outer, beam_wall, angle_threshold
    )

    # Build Mars globe
    mars: solid = centered_mars_globe(globe_diameter, globe_center_z, tilt_angle, mars_oblateness)

    # Union stand parts together
    stand: solid = union(base, cradle, arc1, arc2, arc3)

    # Combine stand and globe as separate solids (no boolean merge)
    result: solid = compound(stand, mars)

    emit result


# Stand only version (no globe) for manufacturing
command CENTERED_GLOBE_STAND_ONLY(
    globe_diameter: float = 304.8,
    base_diameter: float = 400.0,
    tilt_angle: float = 25.2,
    cradle_latitude: float = -35.0,
    waist_ratio: float = 0.7,
    twist_deg: float = 60.0,
    beam_outer: float = 10.0,
    beam_wall: float = 2.0,
    angle_threshold: float = 5.0,
    mars_oblateness: float = 0.00648
) -> solid:
    pi: float = 3.14159265359

    base_radius: float = base_diameter / 2.0
    globe_radius: float = globe_diameter / 2.0
    polar_radius: float = globe_radius * (1.0 - mars_oblateness)
    beam_offset: float = beam_outer / 2.0

    lat_rad: float = cradle_latitude * pi / 180.0
    tilt_rad: float = tilt_angle * pi / 180.0

    cos_lat: float = cos(lat_rad)
    sin_lat: float = sin(lat_rad)
    cos_tilt: float = cos(tilt_rad)
    sin_tilt: float = sin(tilt_rad)

    a2: float = globe_radius * globe_radius
    c2: float = polar_radius * polar_radius

    x270: float = 0.0
    y270: float = -globe_radius * cos_lat
    z270: float = polar_radius * sin_lat

    nx270: float = x270 / a2
    ny270: float = y270 / a2
    nz270: float = z270 / c2
    nlen270: float = sqrt(nx270*nx270 + ny270*ny270 + nz270*nz270)
    px270: float = x270 + beam_offset * nx270 / nlen270
    py270: float = y270 + beam_offset * ny270 / nlen270
    pz270: float = z270 + beam_offset * nz270 / nlen270

    rz270: float = py270 * sin_tilt + pz270 * cos_tilt

    target_cradle_low_z: float = 300.0
    globe_center_z: float = target_cradle_low_z - rz270

    base: solid = base_ring(base_radius, beam_outer, beam_wall, angle_threshold)

    cradle: solid = latitude_cradle_ring(
        globe_diameter, mars_oblateness, cradle_latitude,
        tilt_angle, globe_center_z,
        beam_outer, beam_wall, angle_threshold
    )

    arc1: solid = wrapped_support_arc(
        base_radius, globe_diameter, mars_oblateness, cradle_latitude,
        tilt_angle, globe_center_z, waist_ratio,
        0.0, twist_deg,
        beam_outer, beam_wall, angle_threshold
    )

    arc2: solid = wrapped_support_arc(
        base_radius, globe_diameter, mars_oblateness, cradle_latitude,
        tilt_angle, globe_center_z, waist_ratio,
        120.0, 120.0 + twist_deg,
        beam_outer, beam_wall, angle_threshold
    )

    arc3: solid = wrapped_support_arc(
        base_radius, globe_diameter, mars_oblateness, cradle_latitude,
        tilt_angle, globe_center_z, waist_ratio,
        240.0, 240.0 + twist_deg,
        beam_outer, beam_wall, angle_threshold
    )

    result: solid = union(base, cradle, arc1, arc2, arc3)
    emit result


# ============================================================================
# HIGH-RESOLUTION VERSIONS (2x sampling frequency)
# ============================================================================

# High-resolution latitude cradle path with 32 segments (11.25° increments)
command latitude_cradle_path_hires(
    globe_diameter: float,
    oblateness: float,
    latitude_deg: float,
    tilt_deg: float,
    globe_center_z: float,
    beam_offset: float
) -> path3d:
    require globe_diameter > 0.0
    require oblateness >= 0.0
    require oblateness < 1.0

    pi: float = 3.14159265359

    globe_radius: float = globe_diameter / 2.0
    polar_radius: float = globe_radius * (1.0 - oblateness)

    lat_rad: float = latitude_deg * pi / 180.0
    tilt_rad: float = tilt_deg * pi / 180.0

    cos_lat: float = cos(lat_rad)
    sin_lat: float = sin(lat_rad)
    cos_tilt: float = cos(tilt_rad)
    sin_tilt: float = sin(tilt_rad)

    a2: float = globe_radius * globe_radius
    c2: float = polar_radius * polar_radius

    # Precomputed cos/sin values for 11.25° increments
    cos11: float = 0.98078528040
    sin11: float = 0.19509032201
    cos22: float = 0.92387953251
    sin22: float = 0.38268343236
    cos33: float = 0.83146961230
    sin33: float = 0.55557023302
    cos45: float = 0.70710678118
    sin45: float = 0.70710678118
    cos56: float = 0.55557023302
    sin56: float = 0.83146961230
    cos67: float = 0.38268343236
    sin67: float = 0.92387953251
    cos78: float = 0.19509032201
    sin78: float = 0.98078528040

    # Point 0: longitude 0°
    x0: float = globe_radius * cos_lat * 1.0
    y0: float = globe_radius * cos_lat * 0.0
    z0: float = polar_radius * sin_lat
    nx0: float = x0 / a2
    ny0: float = y0 / a2
    nz0: float = z0 / c2
    nlen0: float = sqrt(nx0*nx0 + ny0*ny0 + nz0*nz0)
    px0: float = x0 + beam_offset * nx0 / nlen0
    py0: float = y0 + beam_offset * ny0 / nlen0
    pz0: float = z0 + beam_offset * nz0 / nlen0
    rx0: float = px0
    ry0: float = py0 * cos_tilt - pz0 * sin_tilt
    rz0: float = py0 * sin_tilt + pz0 * cos_tilt + globe_center_z

    # Point 1: longitude 11.25°
    x1: float = globe_radius * cos_lat * cos11
    y1: float = globe_radius * cos_lat * sin11
    z1: float = polar_radius * sin_lat
    nx1: float = x1 / a2
    ny1: float = y1 / a2
    nz1: float = z1 / c2
    nlen1: float = sqrt(nx1*nx1 + ny1*ny1 + nz1*nz1)
    px1: float = x1 + beam_offset * nx1 / nlen1
    py1: float = y1 + beam_offset * ny1 / nlen1
    pz1: float = z1 + beam_offset * nz1 / nlen1
    rx1: float = px1
    ry1: float = py1 * cos_tilt - pz1 * sin_tilt
    rz1: float = py1 * sin_tilt + pz1 * cos_tilt + globe_center_z

    # Point 2: longitude 22.5°
    x2: float = globe_radius * cos_lat * cos22
    y2: float = globe_radius * cos_lat * sin22
    z2: float = polar_radius * sin_lat
    nx2: float = x2 / a2
    ny2: float = y2 / a2
    nz2: float = z2 / c2
    nlen2: float = sqrt(nx2*nx2 + ny2*ny2 + nz2*nz2)
    px2: float = x2 + beam_offset * nx2 / nlen2
    py2: float = y2 + beam_offset * ny2 / nlen2
    pz2: float = z2 + beam_offset * nz2 / nlen2
    rx2: float = px2
    ry2: float = py2 * cos_tilt - pz2 * sin_tilt
    rz2: float = py2 * sin_tilt + pz2 * cos_tilt + globe_center_z

    # Point 3: longitude 33.75°
    x3: float = globe_radius * cos_lat * cos33
    y3: float = globe_radius * cos_lat * sin33
    z3: float = polar_radius * sin_lat
    nx3: float = x3 / a2
    ny3: float = y3 / a2
    nz3: float = z3 / c2
    nlen3: float = sqrt(nx3*nx3 + ny3*ny3 + nz3*nz3)
    px3: float = x3 + beam_offset * nx3 / nlen3
    py3: float = y3 + beam_offset * ny3 / nlen3
    pz3: float = z3 + beam_offset * nz3 / nlen3
    rx3: float = px3
    ry3: float = py3 * cos_tilt - pz3 * sin_tilt
    rz3: float = py3 * sin_tilt + pz3 * cos_tilt + globe_center_z

    # Point 4: longitude 45°
    x4: float = globe_radius * cos_lat * cos45
    y4: float = globe_radius * cos_lat * sin45
    z4: float = polar_radius * sin_lat
    nx4: float = x4 / a2
    ny4: float = y4 / a2
    nz4: float = z4 / c2
    nlen4: float = sqrt(nx4*nx4 + ny4*ny4 + nz4*nz4)
    px4: float = x4 + beam_offset * nx4 / nlen4
    py4: float = y4 + beam_offset * ny4 / nlen4
    pz4: float = z4 + beam_offset * nz4 / nlen4
    rx4: float = px4
    ry4: float = py4 * cos_tilt - pz4 * sin_tilt
    rz4: float = py4 * sin_tilt + pz4 * cos_tilt + globe_center_z

    # Point 5: longitude 56.25°
    x5: float = globe_radius * cos_lat * cos56
    y5: float = globe_radius * cos_lat * sin56
    z5: float = polar_radius * sin_lat
    nx5: float = x5 / a2
    ny5: float = y5 / a2
    nz5: float = z5 / c2
    nlen5: float = sqrt(nx5*nx5 + ny5*ny5 + nz5*nz5)
    px5: float = x5 + beam_offset * nx5 / nlen5
    py5: float = y5 + beam_offset * ny5 / nlen5
    pz5: float = z5 + beam_offset * nz5 / nlen5
    rx5: float = px5
    ry5: float = py5 * cos_tilt - pz5 * sin_tilt
    rz5: float = py5 * sin_tilt + pz5 * cos_tilt + globe_center_z

    # Point 6: longitude 67.5°
    x6: float = globe_radius * cos_lat * cos67
    y6: float = globe_radius * cos_lat * sin67
    z6: float = polar_radius * sin_lat
    nx6: float = x6 / a2
    ny6: float = y6 / a2
    nz6: float = z6 / c2
    nlen6: float = sqrt(nx6*nx6 + ny6*ny6 + nz6*nz6)
    px6: float = x6 + beam_offset * nx6 / nlen6
    py6: float = y6 + beam_offset * ny6 / nlen6
    pz6: float = z6 + beam_offset * nz6 / nlen6
    rx6: float = px6
    ry6: float = py6 * cos_tilt - pz6 * sin_tilt
    rz6: float = py6 * sin_tilt + pz6 * cos_tilt + globe_center_z

    # Point 7: longitude 78.75°
    x7: float = globe_radius * cos_lat * cos78
    y7: float = globe_radius * cos_lat * sin78
    z7: float = polar_radius * sin_lat
    nx7: float = x7 / a2
    ny7: float = y7 / a2
    nz7: float = z7 / c2
    nlen7: float = sqrt(nx7*nx7 + ny7*ny7 + nz7*nz7)
    px7: float = x7 + beam_offset * nx7 / nlen7
    py7: float = y7 + beam_offset * ny7 / nlen7
    pz7: float = z7 + beam_offset * nz7 / nlen7
    rx7: float = px7
    ry7: float = py7 * cos_tilt - pz7 * sin_tilt
    rz7: float = py7 * sin_tilt + pz7 * cos_tilt + globe_center_z

    # Point 8: longitude 90°
    x8: float = globe_radius * cos_lat * 0.0
    y8: float = globe_radius * cos_lat * 1.0
    z8: float = polar_radius * sin_lat
    nx8: float = x8 / a2
    ny8: float = y8 / a2
    nz8: float = z8 / c2
    nlen8: float = sqrt(nx8*nx8 + ny8*ny8 + nz8*nz8)
    px8: float = x8 + beam_offset * nx8 / nlen8
    py8: float = y8 + beam_offset * ny8 / nlen8
    pz8: float = z8 + beam_offset * nz8 / nlen8
    rx8: float = px8
    ry8: float = py8 * cos_tilt - pz8 * sin_tilt
    rz8: float = py8 * sin_tilt + pz8 * cos_tilt + globe_center_z

    # Point 9: longitude 101.25°
    x9: float = globe_radius * cos_lat * (0.0 - sin11)
    y9: float = globe_radius * cos_lat * cos11
    z9: float = polar_radius * sin_lat
    nx9: float = x9 / a2
    ny9: float = y9 / a2
    nz9: float = z9 / c2
    nlen9: float = sqrt(nx9*nx9 + ny9*ny9 + nz9*nz9)
    px9: float = x9 + beam_offset * nx9 / nlen9
    py9: float = y9 + beam_offset * ny9 / nlen9
    pz9: float = z9 + beam_offset * nz9 / nlen9
    rx9: float = px9
    ry9: float = py9 * cos_tilt - pz9 * sin_tilt
    rz9: float = py9 * sin_tilt + pz9 * cos_tilt + globe_center_z

    # Point 10: longitude 112.5°
    x10: float = globe_radius * cos_lat * (0.0 - sin22)
    y10: float = globe_radius * cos_lat * cos22
    z10: float = polar_radius * sin_lat
    nx10: float = x10 / a2
    ny10: float = y10 / a2
    nz10: float = z10 / c2
    nlen10: float = sqrt(nx10*nx10 + ny10*ny10 + nz10*nz10)
    px10: float = x10 + beam_offset * nx10 / nlen10
    py10: float = y10 + beam_offset * ny10 / nlen10
    pz10: float = z10 + beam_offset * nz10 / nlen10
    rx10: float = px10
    ry10: float = py10 * cos_tilt - pz10 * sin_tilt
    rz10: float = py10 * sin_tilt + pz10 * cos_tilt + globe_center_z

    # Point 11: longitude 123.75°
    x11: float = globe_radius * cos_lat * (0.0 - sin33)
    y11: float = globe_radius * cos_lat * cos33
    z11: float = polar_radius * sin_lat
    nx11: float = x11 / a2
    ny11: float = y11 / a2
    nz11: float = z11 / c2
    nlen11: float = sqrt(nx11*nx11 + ny11*ny11 + nz11*nz11)
    px11: float = x11 + beam_offset * nx11 / nlen11
    py11: float = y11 + beam_offset * ny11 / nlen11
    pz11: float = z11 + beam_offset * nz11 / nlen11
    rx11: float = px11
    ry11: float = py11 * cos_tilt - pz11 * sin_tilt
    rz11: float = py11 * sin_tilt + pz11 * cos_tilt + globe_center_z

    # Point 12: longitude 135°
    x12: float = globe_radius * cos_lat * (0.0 - cos45)
    y12: float = globe_radius * cos_lat * cos45
    z12: float = polar_radius * sin_lat
    nx12: float = x12 / a2
    ny12: float = y12 / a2
    nz12: float = z12 / c2
    nlen12: float = sqrt(nx12*nx12 + ny12*ny12 + nz12*nz12)
    px12: float = x12 + beam_offset * nx12 / nlen12
    py12: float = y12 + beam_offset * ny12 / nlen12
    pz12: float = z12 + beam_offset * nz12 / nlen12
    rx12: float = px12
    ry12: float = py12 * cos_tilt - pz12 * sin_tilt
    rz12: float = py12 * sin_tilt + pz12 * cos_tilt + globe_center_z

    # Point 13: longitude 146.25°
    x13: float = globe_radius * cos_lat * (0.0 - sin56)
    y13: float = globe_radius * cos_lat * cos56
    z13: float = polar_radius * sin_lat
    nx13: float = x13 / a2
    ny13: float = y13 / a2
    nz13: float = z13 / c2
    nlen13: float = sqrt(nx13*nx13 + ny13*ny13 + nz13*nz13)
    px13: float = x13 + beam_offset * nx13 / nlen13
    py13: float = y13 + beam_offset * ny13 / nlen13
    pz13: float = z13 + beam_offset * nz13 / nlen13
    rx13: float = px13
    ry13: float = py13 * cos_tilt - pz13 * sin_tilt
    rz13: float = py13 * sin_tilt + pz13 * cos_tilt + globe_center_z

    # Point 14: longitude 157.5°
    x14: float = globe_radius * cos_lat * (0.0 - sin67)
    y14: float = globe_radius * cos_lat * cos67
    z14: float = polar_radius * sin_lat
    nx14: float = x14 / a2
    ny14: float = y14 / a2
    nz14: float = z14 / c2
    nlen14: float = sqrt(nx14*nx14 + ny14*ny14 + nz14*nz14)
    px14: float = x14 + beam_offset * nx14 / nlen14
    py14: float = y14 + beam_offset * ny14 / nlen14
    pz14: float = z14 + beam_offset * nz14 / nlen14
    rx14: float = px14
    ry14: float = py14 * cos_tilt - pz14 * sin_tilt
    rz14: float = py14 * sin_tilt + pz14 * cos_tilt + globe_center_z

    # Point 15: longitude 168.75°
    x15: float = globe_radius * cos_lat * (0.0 - sin78)
    y15: float = globe_radius * cos_lat * cos78
    z15: float = polar_radius * sin_lat
    nx15: float = x15 / a2
    ny15: float = y15 / a2
    nz15: float = z15 / c2
    nlen15: float = sqrt(nx15*nx15 + ny15*ny15 + nz15*nz15)
    px15: float = x15 + beam_offset * nx15 / nlen15
    py15: float = y15 + beam_offset * ny15 / nlen15
    pz15: float = z15 + beam_offset * nz15 / nlen15
    rx15: float = px15
    ry15: float = py15 * cos_tilt - pz15 * sin_tilt
    rz15: float = py15 * sin_tilt + pz15 * cos_tilt + globe_center_z

    # Point 16: longitude 180°
    x16: float = globe_radius * cos_lat * -1.0
    y16: float = globe_radius * cos_lat * 0.0
    z16: float = polar_radius * sin_lat
    nx16: float = x16 / a2
    ny16: float = y16 / a2
    nz16: float = z16 / c2
    nlen16: float = sqrt(nx16*nx16 + ny16*ny16 + nz16*nz16)
    px16: float = x16 + beam_offset * nx16 / nlen16
    py16: float = y16 + beam_offset * ny16 / nlen16
    pz16: float = z16 + beam_offset * nz16 / nlen16
    rx16: float = px16
    ry16: float = py16 * cos_tilt - pz16 * sin_tilt
    rz16: float = py16 * sin_tilt + pz16 * cos_tilt + globe_center_z

    # Point 17: longitude 191.25°
    x17: float = globe_radius * cos_lat * (0.0 - sin78)
    y17: float = globe_radius * cos_lat * (0.0 - cos78)
    z17: float = polar_radius * sin_lat
    nx17: float = x17 / a2
    ny17: float = y17 / a2
    nz17: float = z17 / c2
    nlen17: float = sqrt(nx17*nx17 + ny17*ny17 + nz17*nz17)
    px17: float = x17 + beam_offset * nx17 / nlen17
    py17: float = y17 + beam_offset * ny17 / nlen17
    pz17: float = z17 + beam_offset * nz17 / nlen17
    rx17: float = px17
    ry17: float = py17 * cos_tilt - pz17 * sin_tilt
    rz17: float = py17 * sin_tilt + pz17 * cos_tilt + globe_center_z

    # Point 18: longitude 202.5°
    x18: float = globe_radius * cos_lat * (0.0 - sin67)
    y18: float = globe_radius * cos_lat * (0.0 - cos67)
    z18: float = polar_radius * sin_lat
    nx18: float = x18 / a2
    ny18: float = y18 / a2
    nz18: float = z18 / c2
    nlen18: float = sqrt(nx18*nx18 + ny18*ny18 + nz18*nz18)
    px18: float = x18 + beam_offset * nx18 / nlen18
    py18: float = y18 + beam_offset * ny18 / nlen18
    pz18: float = z18 + beam_offset * nz18 / nlen18
    rx18: float = px18
    ry18: float = py18 * cos_tilt - pz18 * sin_tilt
    rz18: float = py18 * sin_tilt + pz18 * cos_tilt + globe_center_z

    # Point 19: longitude 213.75°
    x19: float = globe_radius * cos_lat * (0.0 - sin56)
    y19: float = globe_radius * cos_lat * (0.0 - cos56)
    z19: float = polar_radius * sin_lat
    nx19: float = x19 / a2
    ny19: float = y19 / a2
    nz19: float = z19 / c2
    nlen19: float = sqrt(nx19*nx19 + ny19*ny19 + nz19*nz19)
    px19: float = x19 + beam_offset * nx19 / nlen19
    py19: float = y19 + beam_offset * ny19 / nlen19
    pz19: float = z19 + beam_offset * nz19 / nlen19
    rx19: float = px19
    ry19: float = py19 * cos_tilt - pz19 * sin_tilt
    rz19: float = py19 * sin_tilt + pz19 * cos_tilt + globe_center_z

    # Point 20: longitude 225°
    x20: float = globe_radius * cos_lat * (0.0 - cos45)
    y20: float = globe_radius * cos_lat * (0.0 - cos45)
    z20: float = polar_radius * sin_lat
    nx20: float = x20 / a2
    ny20: float = y20 / a2
    nz20: float = z20 / c2
    nlen20: float = sqrt(nx20*nx20 + ny20*ny20 + nz20*nz20)
    px20: float = x20 + beam_offset * nx20 / nlen20
    py20: float = y20 + beam_offset * ny20 / nlen20
    pz20: float = z20 + beam_offset * nz20 / nlen20
    rx20: float = px20
    ry20: float = py20 * cos_tilt - pz20 * sin_tilt
    rz20: float = py20 * sin_tilt + pz20 * cos_tilt + globe_center_z

    # Point 21: longitude 236.25°
    x21: float = globe_radius * cos_lat * (0.0 - sin33)
    y21: float = globe_radius * cos_lat * (0.0 - cos33)
    z21: float = polar_radius * sin_lat
    nx21: float = x21 / a2
    ny21: float = y21 / a2
    nz21: float = z21 / c2
    nlen21: float = sqrt(nx21*nx21 + ny21*ny21 + nz21*nz21)
    px21: float = x21 + beam_offset * nx21 / nlen21
    py21: float = y21 + beam_offset * ny21 / nlen21
    pz21: float = z21 + beam_offset * nz21 / nlen21
    rx21: float = px21
    ry21: float = py21 * cos_tilt - pz21 * sin_tilt
    rz21: float = py21 * sin_tilt + pz21 * cos_tilt + globe_center_z

    # Point 22: longitude 247.5°
    x22: float = globe_radius * cos_lat * (0.0 - sin22)
    y22: float = globe_radius * cos_lat * (0.0 - cos22)
    z22: float = polar_radius * sin_lat
    nx22: float = x22 / a2
    ny22: float = y22 / a2
    nz22: float = z22 / c2
    nlen22: float = sqrt(nx22*nx22 + ny22*ny22 + nz22*nz22)
    px22: float = x22 + beam_offset * nx22 / nlen22
    py22: float = y22 + beam_offset * ny22 / nlen22
    pz22: float = z22 + beam_offset * nz22 / nlen22
    rx22: float = px22
    ry22: float = py22 * cos_tilt - pz22 * sin_tilt
    rz22: float = py22 * sin_tilt + pz22 * cos_tilt + globe_center_z

    # Point 23: longitude 258.75°
    x23: float = globe_radius * cos_lat * (0.0 - sin11)
    y23: float = globe_radius * cos_lat * (0.0 - cos11)
    z23: float = polar_radius * sin_lat
    nx23: float = x23 / a2
    ny23: float = y23 / a2
    nz23: float = z23 / c2
    nlen23: float = sqrt(nx23*nx23 + ny23*ny23 + nz23*nz23)
    px23: float = x23 + beam_offset * nx23 / nlen23
    py23: float = y23 + beam_offset * ny23 / nlen23
    pz23: float = z23 + beam_offset * nz23 / nlen23
    rx23: float = px23
    ry23: float = py23 * cos_tilt - pz23 * sin_tilt
    rz23: float = py23 * sin_tilt + pz23 * cos_tilt + globe_center_z

    # Point 24: longitude 270°
    x24: float = globe_radius * cos_lat * 0.0
    y24: float = globe_radius * cos_lat * -1.0
    z24: float = polar_radius * sin_lat
    nx24: float = x24 / a2
    ny24: float = y24 / a2
    nz24: float = z24 / c2
    nlen24: float = sqrt(nx24*nx24 + ny24*ny24 + nz24*nz24)
    px24: float = x24 + beam_offset * nx24 / nlen24
    py24: float = y24 + beam_offset * ny24 / nlen24
    pz24: float = z24 + beam_offset * nz24 / nlen24
    rx24: float = px24
    ry24: float = py24 * cos_tilt - pz24 * sin_tilt
    rz24: float = py24 * sin_tilt + pz24 * cos_tilt + globe_center_z

    # Point 25: longitude 281.25°
    x25: float = globe_radius * cos_lat * sin11
    y25: float = globe_radius * cos_lat * (0.0 - cos11)
    z25: float = polar_radius * sin_lat
    nx25: float = x25 / a2
    ny25: float = y25 / a2
    nz25: float = z25 / c2
    nlen25: float = sqrt(nx25*nx25 + ny25*ny25 + nz25*nz25)
    px25: float = x25 + beam_offset * nx25 / nlen25
    py25: float = y25 + beam_offset * ny25 / nlen25
    pz25: float = z25 + beam_offset * nz25 / nlen25
    rx25: float = px25
    ry25: float = py25 * cos_tilt - pz25 * sin_tilt
    rz25: float = py25 * sin_tilt + pz25 * cos_tilt + globe_center_z

    # Point 26: longitude 292.5°
    x26: float = globe_radius * cos_lat * sin22
    y26: float = globe_radius * cos_lat * (0.0 - cos22)
    z26: float = polar_radius * sin_lat
    nx26: float = x26 / a2
    ny26: float = y26 / a2
    nz26: float = z26 / c2
    nlen26: float = sqrt(nx26*nx26 + ny26*ny26 + nz26*nz26)
    px26: float = x26 + beam_offset * nx26 / nlen26
    py26: float = y26 + beam_offset * ny26 / nlen26
    pz26: float = z26 + beam_offset * nz26 / nlen26
    rx26: float = px26
    ry26: float = py26 * cos_tilt - pz26 * sin_tilt
    rz26: float = py26 * sin_tilt + pz26 * cos_tilt + globe_center_z

    # Point 27: longitude 303.75°
    x27: float = globe_radius * cos_lat * sin33
    y27: float = globe_radius * cos_lat * (0.0 - cos33)
    z27: float = polar_radius * sin_lat
    nx27: float = x27 / a2
    ny27: float = y27 / a2
    nz27: float = z27 / c2
    nlen27: float = sqrt(nx27*nx27 + ny27*ny27 + nz27*nz27)
    px27: float = x27 + beam_offset * nx27 / nlen27
    py27: float = y27 + beam_offset * ny27 / nlen27
    pz27: float = z27 + beam_offset * nz27 / nlen27
    rx27: float = px27
    ry27: float = py27 * cos_tilt - pz27 * sin_tilt
    rz27: float = py27 * sin_tilt + pz27 * cos_tilt + globe_center_z

    # Point 28: longitude 315°
    x28: float = globe_radius * cos_lat * cos45
    y28: float = globe_radius * cos_lat * (0.0 - cos45)
    z28: float = polar_radius * sin_lat
    nx28: float = x28 / a2
    ny28: float = y28 / a2
    nz28: float = z28 / c2
    nlen28: float = sqrt(nx28*nx28 + ny28*ny28 + nz28*nz28)
    px28: float = x28 + beam_offset * nx28 / nlen28
    py28: float = y28 + beam_offset * ny28 / nlen28
    pz28: float = z28 + beam_offset * nz28 / nlen28
    rx28: float = px28
    ry28: float = py28 * cos_tilt - pz28 * sin_tilt
    rz28: float = py28 * sin_tilt + pz28 * cos_tilt + globe_center_z

    # Point 29: longitude 326.25°
    x29: float = globe_radius * cos_lat * sin56
    y29: float = globe_radius * cos_lat * (0.0 - cos56)
    z29: float = polar_radius * sin_lat
    nx29: float = x29 / a2
    ny29: float = y29 / a2
    nz29: float = z29 / c2
    nlen29: float = sqrt(nx29*nx29 + ny29*ny29 + nz29*nz29)
    px29: float = x29 + beam_offset * nx29 / nlen29
    py29: float = y29 + beam_offset * ny29 / nlen29
    pz29: float = z29 + beam_offset * nz29 / nlen29
    rx29: float = px29
    ry29: float = py29 * cos_tilt - pz29 * sin_tilt
    rz29: float = py29 * sin_tilt + pz29 * cos_tilt + globe_center_z

    # Point 30: longitude 337.5°
    x30: float = globe_radius * cos_lat * sin67
    y30: float = globe_radius * cos_lat * (0.0 - cos67)
    z30: float = polar_radius * sin_lat
    nx30: float = x30 / a2
    ny30: float = y30 / a2
    nz30: float = z30 / c2
    nlen30: float = sqrt(nx30*nx30 + ny30*ny30 + nz30*nz30)
    px30: float = x30 + beam_offset * nx30 / nlen30
    py30: float = y30 + beam_offset * ny30 / nlen30
    pz30: float = z30 + beam_offset * nz30 / nlen30
    rx30: float = px30
    ry30: float = py30 * cos_tilt - pz30 * sin_tilt
    rz30: float = py30 * sin_tilt + pz30 * cos_tilt + globe_center_z

    # Point 31: longitude 348.75°
    x31: float = globe_radius * cos_lat * sin78
    y31: float = globe_radius * cos_lat * (0.0 - cos78)
    z31: float = polar_radius * sin_lat
    nx31: float = x31 / a2
    ny31: float = y31 / a2
    nz31: float = z31 / c2
    nlen31: float = sqrt(nx31*nx31 + ny31*ny31 + nz31*nz31)
    px31: float = x31 + beam_offset * nx31 / nlen31
    py31: float = y31 + beam_offset * ny31 / nlen31
    pz31: float = z31 + beam_offset * nz31 / nlen31
    rx31: float = px31
    ry31: float = py31 * cos_tilt - pz31 * sin_tilt
    rz31: float = py31 * sin_tilt + pz31 * cos_tilt + globe_center_z

    # Build 32 path segments (closing the loop)
    seg0: path3d = path3d_line(point(rx0, ry0, rz0), point(rx1, ry1, rz1))
    seg1: path3d = path3d_line(point(rx1, ry1, rz1), point(rx2, ry2, rz2))
    seg2: path3d = path3d_line(point(rx2, ry2, rz2), point(rx3, ry3, rz3))
    seg3: path3d = path3d_line(point(rx3, ry3, rz3), point(rx4, ry4, rz4))
    seg4: path3d = path3d_line(point(rx4, ry4, rz4), point(rx5, ry5, rz5))
    seg5: path3d = path3d_line(point(rx5, ry5, rz5), point(rx6, ry6, rz6))
    seg6: path3d = path3d_line(point(rx6, ry6, rz6), point(rx7, ry7, rz7))
    seg7: path3d = path3d_line(point(rx7, ry7, rz7), point(rx8, ry8, rz8))
    seg8: path3d = path3d_line(point(rx8, ry8, rz8), point(rx9, ry9, rz9))
    seg9: path3d = path3d_line(point(rx9, ry9, rz9), point(rx10, ry10, rz10))
    seg10: path3d = path3d_line(point(rx10, ry10, rz10), point(rx11, ry11, rz11))
    seg11: path3d = path3d_line(point(rx11, ry11, rz11), point(rx12, ry12, rz12))
    seg12: path3d = path3d_line(point(rx12, ry12, rz12), point(rx13, ry13, rz13))
    seg13: path3d = path3d_line(point(rx13, ry13, rz13), point(rx14, ry14, rz14))
    seg14: path3d = path3d_line(point(rx14, ry14, rz14), point(rx15, ry15, rz15))
    seg15: path3d = path3d_line(point(rx15, ry15, rz15), point(rx16, ry16, rz16))
    seg16: path3d = path3d_line(point(rx16, ry16, rz16), point(rx17, ry17, rz17))
    seg17: path3d = path3d_line(point(rx17, ry17, rz17), point(rx18, ry18, rz18))
    seg18: path3d = path3d_line(point(rx18, ry18, rz18), point(rx19, ry19, rz19))
    seg19: path3d = path3d_line(point(rx19, ry19, rz19), point(rx20, ry20, rz20))
    seg20: path3d = path3d_line(point(rx20, ry20, rz20), point(rx21, ry21, rz21))
    seg21: path3d = path3d_line(point(rx21, ry21, rz21), point(rx22, ry22, rz22))
    seg22: path3d = path3d_line(point(rx22, ry22, rz22), point(rx23, ry23, rz23))
    seg23: path3d = path3d_line(point(rx23, ry23, rz23), point(rx24, ry24, rz24))
    seg24: path3d = path3d_line(point(rx24, ry24, rz24), point(rx25, ry25, rz25))
    seg25: path3d = path3d_line(point(rx25, ry25, rz25), point(rx26, ry26, rz26))
    seg26: path3d = path3d_line(point(rx26, ry26, rz26), point(rx27, ry27, rz27))
    seg27: path3d = path3d_line(point(rx27, ry27, rz27), point(rx28, ry28, rz28))
    seg28: path3d = path3d_line(point(rx28, ry28, rz28), point(rx29, ry29, rz29))
    seg29: path3d = path3d_line(point(rx29, ry29, rz29), point(rx30, ry30, rz30))
    seg30: path3d = path3d_line(point(rx30, ry30, rz30), point(rx31, ry31, rz31))
    seg31: path3d = path3d_line(point(rx31, ry31, rz31), point(rx0, ry0, rz0))

    full_path: path3d = make_path3d(seg0, seg1, seg2, seg3, seg4, seg5, seg6, seg7,
                                    seg8, seg9, seg10, seg11, seg12, seg13, seg14, seg15,
                                    seg16, seg17, seg18, seg19, seg20, seg21, seg22, seg23,
                                    seg24, seg25, seg26, seg27, seg28, seg29, seg30, seg31)
    emit full_path


# High-resolution cradle ring solid (32 segments)
command latitude_cradle_ring_hires(
    globe_diameter: float,
    oblateness: float,
    latitude_deg: float,
    tilt_deg: float,
    globe_center_z: float,
    beam_outer: float,
    beam_wall: float,
    angle_threshold: float
) -> solid:
    outer: region2d = box_beam_profile(beam_outer, beam_wall)
    inner: region2d = box_beam_inner(beam_outer, beam_wall)
    beam_offset: float = beam_outer / 2.0
    path: path3d = latitude_cradle_path_hires(
        globe_diameter, oblateness, latitude_deg, tilt_deg, globe_center_z, beam_offset
    )

    result: solid = sweep_adaptive_hollow(outer, inner, path, angle_threshold)
    emit result


# High-resolution wrapped support arc path with 16 intermediate points and sigmoid blending
command wrapped_support_arc_path_hires(
    base_radius: float,
    globe_diameter: float,
    oblateness: float,
    latitude_deg: float,
    tilt_deg: float,
    globe_center_z: float,
    beam_offset: float,
    waist_ratio: float,
    base_angle_deg: float,
    top_longitude_deg: float
) -> path3d:
    require base_radius > 0.0
    require globe_diameter > 0.0
    require waist_ratio > 0.0
    require waist_ratio < 1.0

    pi: float = 3.14159265359

    globe_radius: float = globe_diameter / 2.0
    polar_radius: float = globe_radius * (1.0 - oblateness)

    # Clearance ellipsoid = globe + beam_offset
    a_clear: float = globe_radius + beam_offset
    c_clear: float = polar_radius + beam_offset
    a_clear2: float = a_clear * a_clear
    c_clear2: float = c_clear * c_clear

    # Sigmoid sharpness parameter
    sigmoid_k: float = 20.0

    lat_rad: float = latitude_deg * pi / 180.0
    tilt_rad: float = tilt_deg * pi / 180.0
    base_angle_rad: float = base_angle_deg * pi / 180.0
    top_lon_rad: float = top_longitude_deg * pi / 180.0

    cos_lat: float = cos(lat_rad)
    sin_lat: float = sin(lat_rad)
    cos_tilt: float = cos(tilt_rad)
    sin_tilt: float = sin(tilt_rad)

    a2: float = globe_radius * globe_radius
    c2: float = polar_radius * polar_radius

    # Base point (on flat ring at z=0) - never needs wrapping
    base_x: float = base_radius * cos(base_angle_rad)
    base_y: float = base_radius * sin(base_angle_rad)
    base_z: float = 0.0

    # Top point: on latitude circle (already on clearance surface)
    tx: float = globe_radius * cos_lat * cos(top_lon_rad)
    ty: float = globe_radius * cos_lat * sin(top_lon_rad)
    tz: float = polar_radius * sin_lat
    tnx: float = tx / a2
    tny: float = ty / a2
    tnz: float = tz / c2
    tnlen: float = sqrt(tnx*tnx + tny*tny + tnz*tnz)
    tpx: float = tx + beam_offset * tnx / tnlen
    tpy: float = ty + beam_offset * tny / tnlen
    tpz: float = tz + beam_offset * tnz / tnlen
    top_x: float = tpx
    top_y: float = tpy * cos_tilt - tpz * sin_tilt
    top_z: float = tpy * sin_tilt + tpz * cos_tilt + globe_center_z

    # Generate 16 intermediate points with sigmoid-blended wrapping
    # Point at t=1/17
    t1: float = 0.0588235
    para1: float = 4.0 * t1 * (1.0 - t1)
    lin_x1: float = base_x + (top_x - base_x) * t1
    lin_y1: float = base_y + (top_y - base_y) * t1
    lin_z1: float = base_z + (top_z - base_z) * t1
    lin_r1: float = sqrt(lin_x1 * lin_x1 + lin_y1 * lin_y1)
    waist_offset1: float = (1.0 - waist_ratio) * lin_r1 * para1
    scale1: float = (lin_r1 - waist_offset1) / lin_r1
    raw_x1: float = lin_x1 * scale1
    raw_y1: float = lin_y1 * scale1
    raw_z1: float = lin_z1
    loc_x1: float = raw_x1
    loc_y1: float = raw_y1 * cos_tilt + (raw_z1 - globe_center_z) * sin_tilt
    loc_z1: float = -raw_y1 * sin_tilt + (raw_z1 - globe_center_z) * cos_tilt
    check1: float = (loc_x1*loc_x1 + loc_y1*loc_y1) / a_clear2 + loc_z1*loc_z1 / c_clear2
    d1: float = check1 - 1.0
    s1: float = 1.0 / (1.0 + exp(sigmoid_k * d1))
    proj_t1: float = 1.0 / sqrt(check1)
    proj_loc_x1: float = loc_x1 * proj_t1
    proj_loc_y1: float = loc_y1 * proj_t1
    proj_loc_z1: float = loc_z1 * proj_t1
    proj_world_x1: float = proj_loc_x1
    proj_world_y1: float = proj_loc_y1 * cos_tilt - proj_loc_z1 * sin_tilt
    proj_world_z1: float = proj_loc_y1 * sin_tilt + proj_loc_z1 * cos_tilt + globe_center_z
    x1: float = (1.0 - s1) * raw_x1 + s1 * proj_world_x1
    y1: float = (1.0 - s1) * raw_y1 + s1 * proj_world_y1
    z1: float = (1.0 - s1) * raw_z1 + s1 * proj_world_z1

    # Point at t=2/17
    t2: float = 0.1176471
    para2: float = 4.0 * t2 * (1.0 - t2)
    lin_x2: float = base_x + (top_x - base_x) * t2
    lin_y2: float = base_y + (top_y - base_y) * t2
    lin_z2: float = base_z + (top_z - base_z) * t2
    lin_r2: float = sqrt(lin_x2 * lin_x2 + lin_y2 * lin_y2)
    waist_offset2: float = (1.0 - waist_ratio) * lin_r2 * para2
    scale2: float = (lin_r2 - waist_offset2) / lin_r2
    raw_x2: float = lin_x2 * scale2
    raw_y2: float = lin_y2 * scale2
    raw_z2: float = lin_z2
    loc_x2: float = raw_x2
    loc_y2: float = raw_y2 * cos_tilt + (raw_z2 - globe_center_z) * sin_tilt
    loc_z2: float = -raw_y2 * sin_tilt + (raw_z2 - globe_center_z) * cos_tilt
    check2: float = (loc_x2*loc_x2 + loc_y2*loc_y2) / a_clear2 + loc_z2*loc_z2 / c_clear2
    d2: float = check2 - 1.0
    s2: float = 1.0 / (1.0 + exp(sigmoid_k * d2))
    proj_t2: float = 1.0 / sqrt(check2)
    proj_loc_x2: float = loc_x2 * proj_t2
    proj_loc_y2: float = loc_y2 * proj_t2
    proj_loc_z2: float = loc_z2 * proj_t2
    proj_world_x2: float = proj_loc_x2
    proj_world_y2: float = proj_loc_y2 * cos_tilt - proj_loc_z2 * sin_tilt
    proj_world_z2: float = proj_loc_y2 * sin_tilt + proj_loc_z2 * cos_tilt + globe_center_z
    x2: float = (1.0 - s2) * raw_x2 + s2 * proj_world_x2
    y2: float = (1.0 - s2) * raw_y2 + s2 * proj_world_y2
    z2: float = (1.0 - s2) * raw_z2 + s2 * proj_world_z2

    # Point at t=3/17
    t3: float = 0.1764706
    para3: float = 4.0 * t3 * (1.0 - t3)
    lin_x3: float = base_x + (top_x - base_x) * t3
    lin_y3: float = base_y + (top_y - base_y) * t3
    lin_z3: float = base_z + (top_z - base_z) * t3
    lin_r3: float = sqrt(lin_x3 * lin_x3 + lin_y3 * lin_y3)
    waist_offset3: float = (1.0 - waist_ratio) * lin_r3 * para3
    scale3: float = (lin_r3 - waist_offset3) / lin_r3
    raw_x3: float = lin_x3 * scale3
    raw_y3: float = lin_y3 * scale3
    raw_z3: float = lin_z3
    loc_x3: float = raw_x3
    loc_y3: float = raw_y3 * cos_tilt + (raw_z3 - globe_center_z) * sin_tilt
    loc_z3: float = -raw_y3 * sin_tilt + (raw_z3 - globe_center_z) * cos_tilt
    check3: float = (loc_x3*loc_x3 + loc_y3*loc_y3) / a_clear2 + loc_z3*loc_z3 / c_clear2
    d3: float = check3 - 1.0
    s3: float = 1.0 / (1.0 + exp(sigmoid_k * d3))
    proj_t3: float = 1.0 / sqrt(check3)
    proj_loc_x3: float = loc_x3 * proj_t3
    proj_loc_y3: float = loc_y3 * proj_t3
    proj_loc_z3: float = loc_z3 * proj_t3
    proj_world_x3: float = proj_loc_x3
    proj_world_y3: float = proj_loc_y3 * cos_tilt - proj_loc_z3 * sin_tilt
    proj_world_z3: float = proj_loc_y3 * sin_tilt + proj_loc_z3 * cos_tilt + globe_center_z
    x3: float = (1.0 - s3) * raw_x3 + s3 * proj_world_x3
    y3: float = (1.0 - s3) * raw_y3 + s3 * proj_world_y3
    z3: float = (1.0 - s3) * raw_z3 + s3 * proj_world_z3

    # Point at t=4/17
    t4: float = 0.2352941
    para4: float = 4.0 * t4 * (1.0 - t4)
    lin_x4: float = base_x + (top_x - base_x) * t4
    lin_y4: float = base_y + (top_y - base_y) * t4
    lin_z4: float = base_z + (top_z - base_z) * t4
    lin_r4: float = sqrt(lin_x4 * lin_x4 + lin_y4 * lin_y4)
    waist_offset4: float = (1.0 - waist_ratio) * lin_r4 * para4
    scale4: float = (lin_r4 - waist_offset4) / lin_r4
    raw_x4: float = lin_x4 * scale4
    raw_y4: float = lin_y4 * scale4
    raw_z4: float = lin_z4
    loc_x4: float = raw_x4
    loc_y4: float = raw_y4 * cos_tilt + (raw_z4 - globe_center_z) * sin_tilt
    loc_z4: float = -raw_y4 * sin_tilt + (raw_z4 - globe_center_z) * cos_tilt
    check4: float = (loc_x4*loc_x4 + loc_y4*loc_y4) / a_clear2 + loc_z4*loc_z4 / c_clear2
    d4: float = check4 - 1.0
    s4: float = 1.0 / (1.0 + exp(sigmoid_k * d4))
    proj_t4: float = 1.0 / sqrt(check4)
    proj_loc_x4: float = loc_x4 * proj_t4
    proj_loc_y4: float = loc_y4 * proj_t4
    proj_loc_z4: float = loc_z4 * proj_t4
    proj_world_x4: float = proj_loc_x4
    proj_world_y4: float = proj_loc_y4 * cos_tilt - proj_loc_z4 * sin_tilt
    proj_world_z4: float = proj_loc_y4 * sin_tilt + proj_loc_z4 * cos_tilt + globe_center_z
    x4: float = (1.0 - s4) * raw_x4 + s4 * proj_world_x4
    y4: float = (1.0 - s4) * raw_y4 + s4 * proj_world_y4
    z4: float = (1.0 - s4) * raw_z4 + s4 * proj_world_z4

    # Point at t=5/17
    t5: float = 0.2941176
    para5: float = 4.0 * t5 * (1.0 - t5)
    lin_x5: float = base_x + (top_x - base_x) * t5
    lin_y5: float = base_y + (top_y - base_y) * t5
    lin_z5: float = base_z + (top_z - base_z) * t5
    lin_r5: float = sqrt(lin_x5 * lin_x5 + lin_y5 * lin_y5)
    waist_offset5: float = (1.0 - waist_ratio) * lin_r5 * para5
    scale5: float = (lin_r5 - waist_offset5) / lin_r5
    raw_x5: float = lin_x5 * scale5
    raw_y5: float = lin_y5 * scale5
    raw_z5: float = lin_z5
    loc_x5: float = raw_x5
    loc_y5: float = raw_y5 * cos_tilt + (raw_z5 - globe_center_z) * sin_tilt
    loc_z5: float = -raw_y5 * sin_tilt + (raw_z5 - globe_center_z) * cos_tilt
    check5: float = (loc_x5*loc_x5 + loc_y5*loc_y5) / a_clear2 + loc_z5*loc_z5 / c_clear2
    d5: float = check5 - 1.0
    s5: float = 1.0 / (1.0 + exp(sigmoid_k * d5))
    proj_t5: float = 1.0 / sqrt(check5)
    proj_loc_x5: float = loc_x5 * proj_t5
    proj_loc_y5: float = loc_y5 * proj_t5
    proj_loc_z5: float = loc_z5 * proj_t5
    proj_world_x5: float = proj_loc_x5
    proj_world_y5: float = proj_loc_y5 * cos_tilt - proj_loc_z5 * sin_tilt
    proj_world_z5: float = proj_loc_y5 * sin_tilt + proj_loc_z5 * cos_tilt + globe_center_z
    x5: float = (1.0 - s5) * raw_x5 + s5 * proj_world_x5
    y5: float = (1.0 - s5) * raw_y5 + s5 * proj_world_y5
    z5: float = (1.0 - s5) * raw_z5 + s5 * proj_world_z5

    # Point at t=6/17
    t6: float = 0.3529412
    para6: float = 4.0 * t6 * (1.0 - t6)
    lin_x6: float = base_x + (top_x - base_x) * t6
    lin_y6: float = base_y + (top_y - base_y) * t6
    lin_z6: float = base_z + (top_z - base_z) * t6
    lin_r6: float = sqrt(lin_x6 * lin_x6 + lin_y6 * lin_y6)
    waist_offset6: float = (1.0 - waist_ratio) * lin_r6 * para6
    scale6: float = (lin_r6 - waist_offset6) / lin_r6
    raw_x6: float = lin_x6 * scale6
    raw_y6: float = lin_y6 * scale6
    raw_z6: float = lin_z6
    loc_x6: float = raw_x6
    loc_y6: float = raw_y6 * cos_tilt + (raw_z6 - globe_center_z) * sin_tilt
    loc_z6: float = -raw_y6 * sin_tilt + (raw_z6 - globe_center_z) * cos_tilt
    check6: float = (loc_x6*loc_x6 + loc_y6*loc_y6) / a_clear2 + loc_z6*loc_z6 / c_clear2
    d6: float = check6 - 1.0
    s6: float = 1.0 / (1.0 + exp(sigmoid_k * d6))
    proj_t6: float = 1.0 / sqrt(check6)
    proj_loc_x6: float = loc_x6 * proj_t6
    proj_loc_y6: float = loc_y6 * proj_t6
    proj_loc_z6: float = loc_z6 * proj_t6
    proj_world_x6: float = proj_loc_x6
    proj_world_y6: float = proj_loc_y6 * cos_tilt - proj_loc_z6 * sin_tilt
    proj_world_z6: float = proj_loc_y6 * sin_tilt + proj_loc_z6 * cos_tilt + globe_center_z
    x6: float = (1.0 - s6) * raw_x6 + s6 * proj_world_x6
    y6: float = (1.0 - s6) * raw_y6 + s6 * proj_world_y6
    z6: float = (1.0 - s6) * raw_z6 + s6 * proj_world_z6

    # Point at t=7/17
    t7: float = 0.4117647
    para7: float = 4.0 * t7 * (1.0 - t7)
    lin_x7: float = base_x + (top_x - base_x) * t7
    lin_y7: float = base_y + (top_y - base_y) * t7
    lin_z7: float = base_z + (top_z - base_z) * t7
    lin_r7: float = sqrt(lin_x7 * lin_x7 + lin_y7 * lin_y7)
    waist_offset7: float = (1.0 - waist_ratio) * lin_r7 * para7
    scale7: float = (lin_r7 - waist_offset7) / lin_r7
    raw_x7: float = lin_x7 * scale7
    raw_y7: float = lin_y7 * scale7
    raw_z7: float = lin_z7
    loc_x7: float = raw_x7
    loc_y7: float = raw_y7 * cos_tilt + (raw_z7 - globe_center_z) * sin_tilt
    loc_z7: float = -raw_y7 * sin_tilt + (raw_z7 - globe_center_z) * cos_tilt
    check7: float = (loc_x7*loc_x7 + loc_y7*loc_y7) / a_clear2 + loc_z7*loc_z7 / c_clear2
    d7: float = check7 - 1.0
    s7: float = 1.0 / (1.0 + exp(sigmoid_k * d7))
    proj_t7: float = 1.0 / sqrt(check7)
    proj_loc_x7: float = loc_x7 * proj_t7
    proj_loc_y7: float = loc_y7 * proj_t7
    proj_loc_z7: float = loc_z7 * proj_t7
    proj_world_x7: float = proj_loc_x7
    proj_world_y7: float = proj_loc_y7 * cos_tilt - proj_loc_z7 * sin_tilt
    proj_world_z7: float = proj_loc_y7 * sin_tilt + proj_loc_z7 * cos_tilt + globe_center_z
    x7: float = (1.0 - s7) * raw_x7 + s7 * proj_world_x7
    y7: float = (1.0 - s7) * raw_y7 + s7 * proj_world_y7
    z7: float = (1.0 - s7) * raw_z7 + s7 * proj_world_z7

    # Point at t=8/17
    t8: float = 0.4705882
    para8: float = 4.0 * t8 * (1.0 - t8)
    lin_x8: float = base_x + (top_x - base_x) * t8
    lin_y8: float = base_y + (top_y - base_y) * t8
    lin_z8: float = base_z + (top_z - base_z) * t8
    lin_r8: float = sqrt(lin_x8 * lin_x8 + lin_y8 * lin_y8)
    waist_offset8: float = (1.0 - waist_ratio) * lin_r8 * para8
    scale8: float = (lin_r8 - waist_offset8) / lin_r8
    raw_x8: float = lin_x8 * scale8
    raw_y8: float = lin_y8 * scale8
    raw_z8: float = lin_z8
    loc_x8: float = raw_x8
    loc_y8: float = raw_y8 * cos_tilt + (raw_z8 - globe_center_z) * sin_tilt
    loc_z8: float = -raw_y8 * sin_tilt + (raw_z8 - globe_center_z) * cos_tilt
    check8: float = (loc_x8*loc_x8 + loc_y8*loc_y8) / a_clear2 + loc_z8*loc_z8 / c_clear2
    d8: float = check8 - 1.0
    s8: float = 1.0 / (1.0 + exp(sigmoid_k * d8))
    proj_t8: float = 1.0 / sqrt(check8)
    proj_loc_x8: float = loc_x8 * proj_t8
    proj_loc_y8: float = loc_y8 * proj_t8
    proj_loc_z8: float = loc_z8 * proj_t8
    proj_world_x8: float = proj_loc_x8
    proj_world_y8: float = proj_loc_y8 * cos_tilt - proj_loc_z8 * sin_tilt
    proj_world_z8: float = proj_loc_y8 * sin_tilt + proj_loc_z8 * cos_tilt + globe_center_z
    x8: float = (1.0 - s8) * raw_x8 + s8 * proj_world_x8
    y8: float = (1.0 - s8) * raw_y8 + s8 * proj_world_y8
    z8: float = (1.0 - s8) * raw_z8 + s8 * proj_world_z8

    # Point at t=9/17
    t9: float = 0.5294118
    para9: float = 4.0 * t9 * (1.0 - t9)
    lin_x9: float = base_x + (top_x - base_x) * t9
    lin_y9: float = base_y + (top_y - base_y) * t9
    lin_z9: float = base_z + (top_z - base_z) * t9
    lin_r9: float = sqrt(lin_x9 * lin_x9 + lin_y9 * lin_y9)
    waist_offset9: float = (1.0 - waist_ratio) * lin_r9 * para9
    scale9: float = (lin_r9 - waist_offset9) / lin_r9
    raw_x9: float = lin_x9 * scale9
    raw_y9: float = lin_y9 * scale9
    raw_z9: float = lin_z9
    loc_x9: float = raw_x9
    loc_y9: float = raw_y9 * cos_tilt + (raw_z9 - globe_center_z) * sin_tilt
    loc_z9: float = -raw_y9 * sin_tilt + (raw_z9 - globe_center_z) * cos_tilt
    check9: float = (loc_x9*loc_x9 + loc_y9*loc_y9) / a_clear2 + loc_z9*loc_z9 / c_clear2
    d9: float = check9 - 1.0
    s9: float = 1.0 / (1.0 + exp(sigmoid_k * d9))
    proj_t9: float = 1.0 / sqrt(check9)
    proj_loc_x9: float = loc_x9 * proj_t9
    proj_loc_y9: float = loc_y9 * proj_t9
    proj_loc_z9: float = loc_z9 * proj_t9
    proj_world_x9: float = proj_loc_x9
    proj_world_y9: float = proj_loc_y9 * cos_tilt - proj_loc_z9 * sin_tilt
    proj_world_z9: float = proj_loc_y9 * sin_tilt + proj_loc_z9 * cos_tilt + globe_center_z
    x9: float = (1.0 - s9) * raw_x9 + s9 * proj_world_x9
    y9: float = (1.0 - s9) * raw_y9 + s9 * proj_world_y9
    z9: float = (1.0 - s9) * raw_z9 + s9 * proj_world_z9

    # Point at t=10/17
    t10: float = 0.5882353
    para10: float = 4.0 * t10 * (1.0 - t10)
    lin_x10: float = base_x + (top_x - base_x) * t10
    lin_y10: float = base_y + (top_y - base_y) * t10
    lin_z10: float = base_z + (top_z - base_z) * t10
    lin_r10: float = sqrt(lin_x10 * lin_x10 + lin_y10 * lin_y10)
    waist_offset10: float = (1.0 - waist_ratio) * lin_r10 * para10
    scale10: float = (lin_r10 - waist_offset10) / lin_r10
    raw_x10: float = lin_x10 * scale10
    raw_y10: float = lin_y10 * scale10
    raw_z10: float = lin_z10
    loc_x10: float = raw_x10
    loc_y10: float = raw_y10 * cos_tilt + (raw_z10 - globe_center_z) * sin_tilt
    loc_z10: float = -raw_y10 * sin_tilt + (raw_z10 - globe_center_z) * cos_tilt
    check10: float = (loc_x10*loc_x10 + loc_y10*loc_y10) / a_clear2 + loc_z10*loc_z10 / c_clear2
    d10: float = check10 - 1.0
    s10: float = 1.0 / (1.0 + exp(sigmoid_k * d10))
    proj_t10: float = 1.0 / sqrt(check10)
    proj_loc_x10: float = loc_x10 * proj_t10
    proj_loc_y10: float = loc_y10 * proj_t10
    proj_loc_z10: float = loc_z10 * proj_t10
    proj_world_x10: float = proj_loc_x10
    proj_world_y10: float = proj_loc_y10 * cos_tilt - proj_loc_z10 * sin_tilt
    proj_world_z10: float = proj_loc_y10 * sin_tilt + proj_loc_z10 * cos_tilt + globe_center_z
    x10: float = (1.0 - s10) * raw_x10 + s10 * proj_world_x10
    y10: float = (1.0 - s10) * raw_y10 + s10 * proj_world_y10
    z10: float = (1.0 - s10) * raw_z10 + s10 * proj_world_z10

    # Point at t=11/17
    t11: float = 0.6470588
    para11: float = 4.0 * t11 * (1.0 - t11)
    lin_x11: float = base_x + (top_x - base_x) * t11
    lin_y11: float = base_y + (top_y - base_y) * t11
    lin_z11: float = base_z + (top_z - base_z) * t11
    lin_r11: float = sqrt(lin_x11 * lin_x11 + lin_y11 * lin_y11)
    waist_offset11: float = (1.0 - waist_ratio) * lin_r11 * para11
    scale11: float = (lin_r11 - waist_offset11) / lin_r11
    raw_x11: float = lin_x11 * scale11
    raw_y11: float = lin_y11 * scale11
    raw_z11: float = lin_z11
    loc_x11: float = raw_x11
    loc_y11: float = raw_y11 * cos_tilt + (raw_z11 - globe_center_z) * sin_tilt
    loc_z11: float = -raw_y11 * sin_tilt + (raw_z11 - globe_center_z) * cos_tilt
    check11: float = (loc_x11*loc_x11 + loc_y11*loc_y11) / a_clear2 + loc_z11*loc_z11 / c_clear2
    d11: float = check11 - 1.0
    s11: float = 1.0 / (1.0 + exp(sigmoid_k * d11))
    proj_t11: float = 1.0 / sqrt(check11)
    proj_loc_x11: float = loc_x11 * proj_t11
    proj_loc_y11: float = loc_y11 * proj_t11
    proj_loc_z11: float = loc_z11 * proj_t11
    proj_world_x11: float = proj_loc_x11
    proj_world_y11: float = proj_loc_y11 * cos_tilt - proj_loc_z11 * sin_tilt
    proj_world_z11: float = proj_loc_y11 * sin_tilt + proj_loc_z11 * cos_tilt + globe_center_z
    x11: float = (1.0 - s11) * raw_x11 + s11 * proj_world_x11
    y11: float = (1.0 - s11) * raw_y11 + s11 * proj_world_y11
    z11: float = (1.0 - s11) * raw_z11 + s11 * proj_world_z11

    # Point at t=12/17
    t12: float = 0.7058824
    para12: float = 4.0 * t12 * (1.0 - t12)
    lin_x12: float = base_x + (top_x - base_x) * t12
    lin_y12: float = base_y + (top_y - base_y) * t12
    lin_z12: float = base_z + (top_z - base_z) * t12
    lin_r12: float = sqrt(lin_x12 * lin_x12 + lin_y12 * lin_y12)
    waist_offset12: float = (1.0 - waist_ratio) * lin_r12 * para12
    scale12: float = (lin_r12 - waist_offset12) / lin_r12
    raw_x12: float = lin_x12 * scale12
    raw_y12: float = lin_y12 * scale12
    raw_z12: float = lin_z12
    loc_x12: float = raw_x12
    loc_y12: float = raw_y12 * cos_tilt + (raw_z12 - globe_center_z) * sin_tilt
    loc_z12: float = -raw_y12 * sin_tilt + (raw_z12 - globe_center_z) * cos_tilt
    check12: float = (loc_x12*loc_x12 + loc_y12*loc_y12) / a_clear2 + loc_z12*loc_z12 / c_clear2
    d12: float = check12 - 1.0
    s12: float = 1.0 / (1.0 + exp(sigmoid_k * d12))
    proj_t12: float = 1.0 / sqrt(check12)
    proj_loc_x12: float = loc_x12 * proj_t12
    proj_loc_y12: float = loc_y12 * proj_t12
    proj_loc_z12: float = loc_z12 * proj_t12
    proj_world_x12: float = proj_loc_x12
    proj_world_y12: float = proj_loc_y12 * cos_tilt - proj_loc_z12 * sin_tilt
    proj_world_z12: float = proj_loc_y12 * sin_tilt + proj_loc_z12 * cos_tilt + globe_center_z
    x12: float = (1.0 - s12) * raw_x12 + s12 * proj_world_x12
    y12: float = (1.0 - s12) * raw_y12 + s12 * proj_world_y12
    z12: float = (1.0 - s12) * raw_z12 + s12 * proj_world_z12

    # Point at t=13/17
    t13: float = 0.7647059
    para13: float = 4.0 * t13 * (1.0 - t13)
    lin_x13: float = base_x + (top_x - base_x) * t13
    lin_y13: float = base_y + (top_y - base_y) * t13
    lin_z13: float = base_z + (top_z - base_z) * t13
    lin_r13: float = sqrt(lin_x13 * lin_x13 + lin_y13 * lin_y13)
    waist_offset13: float = (1.0 - waist_ratio) * lin_r13 * para13
    scale13: float = (lin_r13 - waist_offset13) / lin_r13
    raw_x13: float = lin_x13 * scale13
    raw_y13: float = lin_y13 * scale13
    raw_z13: float = lin_z13
    loc_x13: float = raw_x13
    loc_y13: float = raw_y13 * cos_tilt + (raw_z13 - globe_center_z) * sin_tilt
    loc_z13: float = -raw_y13 * sin_tilt + (raw_z13 - globe_center_z) * cos_tilt
    check13: float = (loc_x13*loc_x13 + loc_y13*loc_y13) / a_clear2 + loc_z13*loc_z13 / c_clear2
    d13: float = check13 - 1.0
    s13: float = 1.0 / (1.0 + exp(sigmoid_k * d13))
    proj_t13: float = 1.0 / sqrt(check13)
    proj_loc_x13: float = loc_x13 * proj_t13
    proj_loc_y13: float = loc_y13 * proj_t13
    proj_loc_z13: float = loc_z13 * proj_t13
    proj_world_x13: float = proj_loc_x13
    proj_world_y13: float = proj_loc_y13 * cos_tilt - proj_loc_z13 * sin_tilt
    proj_world_z13: float = proj_loc_y13 * sin_tilt + proj_loc_z13 * cos_tilt + globe_center_z
    x13: float = (1.0 - s13) * raw_x13 + s13 * proj_world_x13
    y13: float = (1.0 - s13) * raw_y13 + s13 * proj_world_y13
    z13: float = (1.0 - s13) * raw_z13 + s13 * proj_world_z13

    # Point at t=14/17
    t14: float = 0.8235294
    para14: float = 4.0 * t14 * (1.0 - t14)
    lin_x14: float = base_x + (top_x - base_x) * t14
    lin_y14: float = base_y + (top_y - base_y) * t14
    lin_z14: float = base_z + (top_z - base_z) * t14
    lin_r14: float = sqrt(lin_x14 * lin_x14 + lin_y14 * lin_y14)
    waist_offset14: float = (1.0 - waist_ratio) * lin_r14 * para14
    scale14: float = (lin_r14 - waist_offset14) / lin_r14
    raw_x14: float = lin_x14 * scale14
    raw_y14: float = lin_y14 * scale14
    raw_z14: float = lin_z14
    loc_x14: float = raw_x14
    loc_y14: float = raw_y14 * cos_tilt + (raw_z14 - globe_center_z) * sin_tilt
    loc_z14: float = -raw_y14 * sin_tilt + (raw_z14 - globe_center_z) * cos_tilt
    check14: float = (loc_x14*loc_x14 + loc_y14*loc_y14) / a_clear2 + loc_z14*loc_z14 / c_clear2
    d14: float = check14 - 1.0
    s14: float = 1.0 / (1.0 + exp(sigmoid_k * d14))
    proj_t14: float = 1.0 / sqrt(check14)
    proj_loc_x14: float = loc_x14 * proj_t14
    proj_loc_y14: float = loc_y14 * proj_t14
    proj_loc_z14: float = loc_z14 * proj_t14
    proj_world_x14: float = proj_loc_x14
    proj_world_y14: float = proj_loc_y14 * cos_tilt - proj_loc_z14 * sin_tilt
    proj_world_z14: float = proj_loc_y14 * sin_tilt + proj_loc_z14 * cos_tilt + globe_center_z
    x14: float = (1.0 - s14) * raw_x14 + s14 * proj_world_x14
    y14: float = (1.0 - s14) * raw_y14 + s14 * proj_world_y14
    z14: float = (1.0 - s14) * raw_z14 + s14 * proj_world_z14

    # Point at t=15/17
    t15: float = 0.8823529
    para15: float = 4.0 * t15 * (1.0 - t15)
    lin_x15: float = base_x + (top_x - base_x) * t15
    lin_y15: float = base_y + (top_y - base_y) * t15
    lin_z15: float = base_z + (top_z - base_z) * t15
    lin_r15: float = sqrt(lin_x15 * lin_x15 + lin_y15 * lin_y15)
    waist_offset15: float = (1.0 - waist_ratio) * lin_r15 * para15
    scale15: float = (lin_r15 - waist_offset15) / lin_r15
    raw_x15: float = lin_x15 * scale15
    raw_y15: float = lin_y15 * scale15
    raw_z15: float = lin_z15
    loc_x15: float = raw_x15
    loc_y15: float = raw_y15 * cos_tilt + (raw_z15 - globe_center_z) * sin_tilt
    loc_z15: float = -raw_y15 * sin_tilt + (raw_z15 - globe_center_z) * cos_tilt
    check15: float = (loc_x15*loc_x15 + loc_y15*loc_y15) / a_clear2 + loc_z15*loc_z15 / c_clear2
    d15: float = check15 - 1.0
    s15: float = 1.0 / (1.0 + exp(sigmoid_k * d15))
    proj_t15: float = 1.0 / sqrt(check15)
    proj_loc_x15: float = loc_x15 * proj_t15
    proj_loc_y15: float = loc_y15 * proj_t15
    proj_loc_z15: float = loc_z15 * proj_t15
    proj_world_x15: float = proj_loc_x15
    proj_world_y15: float = proj_loc_y15 * cos_tilt - proj_loc_z15 * sin_tilt
    proj_world_z15: float = proj_loc_y15 * sin_tilt + proj_loc_z15 * cos_tilt + globe_center_z
    x15: float = (1.0 - s15) * raw_x15 + s15 * proj_world_x15
    y15: float = (1.0 - s15) * raw_y15 + s15 * proj_world_y15
    z15: float = (1.0 - s15) * raw_z15 + s15 * proj_world_z15

    # Point at t=16/17
    t16: float = 0.9411765
    para16: float = 4.0 * t16 * (1.0 - t16)
    lin_x16: float = base_x + (top_x - base_x) * t16
    lin_y16: float = base_y + (top_y - base_y) * t16
    lin_z16: float = base_z + (top_z - base_z) * t16
    lin_r16: float = sqrt(lin_x16 * lin_x16 + lin_y16 * lin_y16)
    waist_offset16: float = (1.0 - waist_ratio) * lin_r16 * para16
    scale16: float = (lin_r16 - waist_offset16) / lin_r16
    raw_x16: float = lin_x16 * scale16
    raw_y16: float = lin_y16 * scale16
    raw_z16: float = lin_z16
    loc_x16: float = raw_x16
    loc_y16: float = raw_y16 * cos_tilt + (raw_z16 - globe_center_z) * sin_tilt
    loc_z16: float = -raw_y16 * sin_tilt + (raw_z16 - globe_center_z) * cos_tilt
    check16: float = (loc_x16*loc_x16 + loc_y16*loc_y16) / a_clear2 + loc_z16*loc_z16 / c_clear2
    d16: float = check16 - 1.0
    s16: float = 1.0 / (1.0 + exp(sigmoid_k * d16))
    proj_t16: float = 1.0 / sqrt(check16)
    proj_loc_x16: float = loc_x16 * proj_t16
    proj_loc_y16: float = loc_y16 * proj_t16
    proj_loc_z16: float = loc_z16 * proj_t16
    proj_world_x16: float = proj_loc_x16
    proj_world_y16: float = proj_loc_y16 * cos_tilt - proj_loc_z16 * sin_tilt
    proj_world_z16: float = proj_loc_y16 * sin_tilt + proj_loc_z16 * cos_tilt + globe_center_z
    x16: float = (1.0 - s16) * raw_x16 + s16 * proj_world_x16
    y16: float = (1.0 - s16) * raw_y16 + s16 * proj_world_y16
    z16: float = (1.0 - s16) * raw_z16 + s16 * proj_world_z16

    # Build 17 path segments
    seg0: path3d = path3d_line(point(base_x, base_y, base_z), point(x1, y1, z1))
    seg1: path3d = path3d_line(point(x1, y1, z1), point(x2, y2, z2))
    seg2: path3d = path3d_line(point(x2, y2, z2), point(x3, y3, z3))
    seg3: path3d = path3d_line(point(x3, y3, z3), point(x4, y4, z4))
    seg4: path3d = path3d_line(point(x4, y4, z4), point(x5, y5, z5))
    seg5: path3d = path3d_line(point(x5, y5, z5), point(x6, y6, z6))
    seg6: path3d = path3d_line(point(x6, y6, z6), point(x7, y7, z7))
    seg7: path3d = path3d_line(point(x7, y7, z7), point(x8, y8, z8))
    seg8: path3d = path3d_line(point(x8, y8, z8), point(x9, y9, z9))
    seg9: path3d = path3d_line(point(x9, y9, z9), point(x10, y10, z10))
    seg10: path3d = path3d_line(point(x10, y10, z10), point(x11, y11, z11))
    seg11: path3d = path3d_line(point(x11, y11, z11), point(x12, y12, z12))
    seg12: path3d = path3d_line(point(x12, y12, z12), point(x13, y13, z13))
    seg13: path3d = path3d_line(point(x13, y13, z13), point(x14, y14, z14))
    seg14: path3d = path3d_line(point(x14, y14, z14), point(x15, y15, z15))
    seg15: path3d = path3d_line(point(x15, y15, z15), point(x16, y16, z16))
    seg16: path3d = path3d_line(point(x16, y16, z16), point(top_x, top_y, top_z))

    full_arc: path3d = make_path3d(seg0, seg1, seg2, seg3, seg4, seg5, seg6, seg7,
                                   seg8, seg9, seg10, seg11, seg12, seg13, seg14, seg15, seg16)
    emit full_arc


# High-resolution wrapped support arc solid
command wrapped_support_arc_hires(
    base_radius: float,
    globe_diameter: float,
    oblateness: float,
    latitude_deg: float,
    tilt_deg: float,
    globe_center_z: float,
    waist_ratio: float,
    base_angle_deg: float,
    top_longitude_deg: float,
    beam_outer: float,
    beam_wall: float,
    angle_threshold: float
) -> solid:
    outer: region2d = box_beam_profile(beam_outer, beam_wall)
    inner: region2d = box_beam_inner(beam_outer, beam_wall)
    beam_offset: float = beam_outer / 2.0
    path: path3d = wrapped_support_arc_path_hires(
        base_radius, globe_diameter, oblateness, latitude_deg,
        tilt_deg, globe_center_z, beam_offset,
        waist_ratio, base_angle_deg, top_longitude_deg
    )

    result: solid = sweep_adaptive_hollow(outer, inner, path, angle_threshold)
    emit result


# High-resolution globe stand with 2x sampling (32-seg cradle, 16-pt arcs)
command CENTERED_GLOBE_STAND_HIRES(
    globe_diameter: float = 304.8,
    base_diameter: float = 400.0,
    tilt_angle: float = 25.2,
    cradle_latitude: float = -35.0,
    waist_ratio: float = 0.7,
    twist_deg: float = 60.0,
    beam_outer: float = 10.0,
    beam_wall: float = 2.0,
    angle_threshold: float = 5.0,
    mars_oblateness: float = 0.00648
) -> solid:
    pi: float = 3.14159265359
    base_radius: float = base_diameter / 2.0
    globe_radius: float = globe_diameter / 2.0
    polar_radius: float = globe_radius * (1.0 - mars_oblateness)
    beam_offset: float = beam_outer / 2.0
    lat_rad: float = cradle_latitude * pi / 180.0
    tilt_rad: float = tilt_angle * pi / 180.0
    cos_lat: float = cos(lat_rad)
    sin_lat: float = sin(lat_rad)
    cos_tilt: float = cos(tilt_rad)
    sin_tilt: float = sin(tilt_rad)
    a2: float = globe_radius * globe_radius
    c2: float = polar_radius * polar_radius
    x270: float = 0.0
    y270: float = -globe_radius * cos_lat
    z270: float = polar_radius * sin_lat
    nx270: float = x270 / a2
    ny270: float = y270 / a2
    nz270: float = z270 / c2
    nlen270: float = sqrt(nx270*nx270 + ny270*ny270 + nz270*nz270)
    px270: float = x270 + beam_offset * nx270 / nlen270
    py270: float = y270 + beam_offset * ny270 / nlen270
    pz270: float = z270 + beam_offset * nz270 / nlen270
    rz270: float = py270 * sin_tilt + pz270 * cos_tilt
    target_cradle_low_z: float = 300.0
    globe_center_z: float = target_cradle_low_z - rz270
    base: solid = base_ring(base_radius, beam_outer, beam_wall, angle_threshold)
    cradle: solid = latitude_cradle_ring_hires(globe_diameter, mars_oblateness, cradle_latitude, tilt_angle, globe_center_z, beam_outer, beam_wall, angle_threshold)
    arc1: solid = wrapped_support_arc_hires(base_radius, globe_diameter, mars_oblateness, cradle_latitude, tilt_angle, globe_center_z, waist_ratio, 0.0, twist_deg, beam_outer, beam_wall, angle_threshold)
    arc2: solid = wrapped_support_arc_hires(base_radius, globe_diameter, mars_oblateness, cradle_latitude, tilt_angle, globe_center_z, waist_ratio, 120.0, 120.0 + twist_deg, beam_outer, beam_wall, angle_threshold)
    arc3: solid = wrapped_support_arc_hires(base_radius, globe_diameter, mars_oblateness, cradle_latitude, tilt_angle, globe_center_z, waist_ratio, 240.0, 240.0 + twist_deg, beam_outer, beam_wall, angle_threshold)
    mars: solid = centered_mars_globe(globe_diameter, globe_center_z, tilt_angle, mars_oblateness)
    stand: solid = union(base, cradle, arc1, arc2, arc3)
    result: solid = compound(stand, mars)
    # result: solid = stand
    emit result


# High-resolution with 90° twist
command CENTERED_GLOBE_STAND_HIRES_TWIST90(
    globe_diameter: float = 304.8,
    base_diameter: float = 400.0,
    tilt_angle: float = 25.2,
    cradle_latitude: float = -35.0,
    waist_ratio: float = 0.7,
    beam_outer: float = 10.0,
    beam_wall: float = 2.0,
    angle_threshold: float = 5.0,
    mars_oblateness: float = 0.00648
) -> solid:
    twist_deg: float = 90.0
    pi: float = 3.14159265359
    base_radius: float = base_diameter / 2.0
    globe_radius: float = globe_diameter / 2.0
    polar_radius: float = globe_radius * (1.0 - mars_oblateness)
    beam_offset: float = beam_outer / 2.0
    lat_rad: float = cradle_latitude * pi / 180.0
    tilt_rad: float = tilt_angle * pi / 180.0
    cos_lat: float = cos(lat_rad)
    sin_lat: float = sin(lat_rad)
    cos_tilt: float = cos(tilt_rad)
    sin_tilt: float = sin(tilt_rad)
    a2: float = globe_radius * globe_radius
    c2: float = polar_radius * polar_radius
    x270: float = 0.0
    y270: float = -globe_radius * cos_lat
    z270: float = polar_radius * sin_lat
    nx270: float = x270 / a2
    ny270: float = y270 / a2
    nz270: float = z270 / c2
    nlen270: float = sqrt(nx270*nx270 + ny270*ny270 + nz270*nz270)
    px270: float = x270 + beam_offset * nx270 / nlen270
    py270: float = y270 + beam_offset * ny270 / nlen270
    pz270: float = z270 + beam_offset * nz270 / nlen270
    rz270: float = py270 * sin_tilt + pz270 * cos_tilt
    target_cradle_low_z: float = 300.0
    globe_center_z: float = target_cradle_low_z - rz270
    base: solid = base_ring(base_radius, beam_outer, beam_wall, angle_threshold)
    cradle: solid = latitude_cradle_ring_hires(globe_diameter, mars_oblateness, cradle_latitude, tilt_angle, globe_center_z, beam_outer, beam_wall, angle_threshold)
    arc1: solid = wrapped_support_arc_hires(base_radius, globe_diameter, mars_oblateness, cradle_latitude, tilt_angle, globe_center_z, waist_ratio, 0.0, twist_deg, beam_outer, beam_wall, angle_threshold)
    arc2: solid = wrapped_support_arc_hires(base_radius, globe_diameter, mars_oblateness, cradle_latitude, tilt_angle, globe_center_z, waist_ratio, 120.0, 120.0 + twist_deg, beam_outer, beam_wall, angle_threshold)
    arc3: solid = wrapped_support_arc_hires(base_radius, globe_diameter, mars_oblateness, cradle_latitude, tilt_angle, globe_center_z, waist_ratio, 240.0, 240.0 + twist_deg, beam_outer, beam_wall, angle_threshold)
    mars: solid = centered_mars_globe(globe_diameter, globe_center_z, tilt_angle, mars_oblateness)
    stand: solid = union(base, cradle, arc1, arc2, arc3)
    result: solid = compound(stand, mars)
    emit result


# High-resolution with 120° twist
command CENTERED_GLOBE_STAND_HIRES_TWIST120(
    globe_diameter: float = 304.8,
    base_diameter: float = 400.0,
    tilt_angle: float = 25.2,
    cradle_latitude: float = -35.0,
    waist_ratio: float = 0.7,
    beam_outer: float = 10.0,
    beam_wall: float = 2.0,
    angle_threshold: float = 5.0,
    mars_oblateness: float = 0.00648
) -> solid:
    twist_deg: float = 120.0
    pi: float = 3.14159265359
    base_radius: float = base_diameter / 2.0
    globe_radius: float = globe_diameter / 2.0
    polar_radius: float = globe_radius * (1.0 - mars_oblateness)
    beam_offset: float = beam_outer / 2.0
    lat_rad: float = cradle_latitude * pi / 180.0
    tilt_rad: float = tilt_angle * pi / 180.0
    cos_lat: float = cos(lat_rad)
    sin_lat: float = sin(lat_rad)
    cos_tilt: float = cos(tilt_rad)
    sin_tilt: float = sin(tilt_rad)
    a2: float = globe_radius * globe_radius
    c2: float = polar_radius * polar_radius
    x270: float = 0.0
    y270: float = -globe_radius * cos_lat
    z270: float = polar_radius * sin_lat
    nx270: float = x270 / a2
    ny270: float = y270 / a2
    nz270: float = z270 / c2
    nlen270: float = sqrt(nx270*nx270 + ny270*ny270 + nz270*nz270)
    px270: float = x270 + beam_offset * nx270 / nlen270
    py270: float = y270 + beam_offset * ny270 / nlen270
    pz270: float = z270 + beam_offset * nz270 / nlen270
    rz270: float = py270 * sin_tilt + pz270 * cos_tilt
    target_cradle_low_z: float = 300.0
    globe_center_z: float = target_cradle_low_z - rz270
    base: solid = base_ring(base_radius, beam_outer, beam_wall, angle_threshold)
    cradle: solid = latitude_cradle_ring_hires(globe_diameter, mars_oblateness, cradle_latitude, tilt_angle, globe_center_z, beam_outer, beam_wall, angle_threshold)
    arc1: solid = wrapped_support_arc_hires(base_radius, globe_diameter, mars_oblateness, cradle_latitude, tilt_angle, globe_center_z, waist_ratio, 0.0, twist_deg, beam_outer, beam_wall, angle_threshold)
    arc2: solid = wrapped_support_arc_hires(base_radius, globe_diameter, mars_oblateness, cradle_latitude, tilt_angle, globe_center_z, waist_ratio, 120.0, 120.0 + twist_deg, beam_outer, beam_wall, angle_threshold)
    arc3: solid = wrapped_support_arc_hires(base_radius, globe_diameter, mars_oblateness, cradle_latitude, tilt_angle, globe_center_z, waist_ratio, 240.0, 240.0 + twist_deg, beam_outer, beam_wall, angle_threshold)
    mars: solid = centered_mars_globe(globe_diameter, globe_center_z, tilt_angle, mars_oblateness)
    stand: solid = union(base, cradle, arc1, arc2, arc3)
    result: solid = compound(stand, mars)
    emit result
