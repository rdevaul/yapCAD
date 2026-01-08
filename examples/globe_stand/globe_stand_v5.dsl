module globe_stand_v5

# Globe Stand v5 - Concise version using DSL Phase 1 & 2 features
#
# This is a functionally equivalent version of CENTERED_GLOBE_STAND_HIRES from
# globe_stand_v4.dsl, rewritten to use:
# - Conditional expressions (Phase 1): value if condition else other
# - List comprehensions (Phase 2): [expr for var in iterable]
# - Runtime trig computation instead of precomputed constants
#
# Original v4: ~1700 lines for these commands
# v5: ~250 lines (85% reduction)


# ============================================================================
# Constants
# ============================================================================

command pi() -> float:
    emit 3.14159265359

command deg2rad(degrees: float) -> float:
    emit degrees * pi() / 180.0


# ============================================================================
# Profile definitions
# ============================================================================

command box_beam_profile(outer_size: float, wall_thickness: float) -> region2d:
    require outer_size > 0.0
    require wall_thickness > 0.0
    require wall_thickness < outer_size / 2.0

    half: float = outer_size / 2.0
    pts: list<point> = [
        point(-half, -half),
        point(half, -half),
        point(half, half),
        point(-half, half)
    ]
    emit polygon(pts)


command box_beam_inner(outer_size: float, wall_thickness: float) -> region2d:
    require outer_size > 0.0
    require wall_thickness > 0.0
    require wall_thickness < outer_size / 2.0

    inner_half: float = outer_size / 2.0 - wall_thickness
    pts: list<point> = [
        point(-inner_half, -inner_half),
        point(inner_half, -inner_half),
        point(inner_half, inner_half),
        point(-inner_half, inner_half)
    ]
    emit polygon(pts)


# ============================================================================
# Ring path generation
# ============================================================================

command ring_path(radius: float, z_height: float) -> path3d:
    require radius > 0.0

    center: point = point(0.0, 0.0, z_height)

    # Four 90-degree arc segments
    seg0: path3d = path3d_arc_auto(center, point(radius, 0.0, z_height), point(0.0, radius, z_height), false)
    seg1: path3d = path3d_arc_auto(center, point(0.0, radius, z_height), point(-radius, 0.0, z_height), false)
    seg2: path3d = path3d_arc_auto(center, point(-radius, 0.0, z_height), point(0.0, -radius, z_height), false)
    seg3: path3d = path3d_arc_auto(center, point(0.0, -radius, z_height), point(radius, 0.0, z_height), false)

    emit make_path3d(seg0, seg1, seg2, seg3)


# ============================================================================
# Base ring
# ============================================================================

command base_ring(
    radius: float,
    beam_outer: float,
    beam_wall: float,
    angle_threshold: float
) -> solid:
    outer: region2d = box_beam_profile(beam_outer, beam_wall)
    inner: region2d = box_beam_inner(beam_outer, beam_wall)
    path: path3d = ring_path(radius, 0.0)

    emit sweep_adaptive_hollow(outer, inner, path, angle_threshold)


# ============================================================================
# Latitude cradle path - computed point on tilted oblate spheroid
# ============================================================================

# Helper: compute a single point on the latitude circle
# Returns the world-space coordinates after tilt rotation
command cradle_point(
    globe_radius: float,
    polar_radius: float,
    cos_lat: float,
    sin_lat: float,
    cos_tilt: float,
    sin_tilt: float,
    beam_offset: float,
    globe_center_z: float,
    longitude_deg: float
) -> point:
    # Compute longitude trig
    lon_rad: float = deg2rad(longitude_deg)
    cos_lon: float = cos(lon_rad)
    sin_lon: float = sin(lon_rad)

    # Point on latitude circle (before tilt)
    x: float = globe_radius * cos_lat * cos_lon
    y: float = globe_radius * cos_lat * sin_lon
    z: float = polar_radius * sin_lat

    # Normal vector for offset (gradient of ellipsoid)
    a2: float = globe_radius * globe_radius
    c2: float = polar_radius * polar_radius
    nx: float = x / a2
    ny: float = y / a2
    nz: float = z / c2
    nlen: float = sqrt(nx*nx + ny*ny + nz*nz)

    # Offset point along normal
    px: float = x + beam_offset * nx / nlen
    py: float = y + beam_offset * ny / nlen
    pz: float = z + beam_offset * nz / nlen

    # Apply tilt rotation around X axis, then translate
    rx: float = px
    ry: float = py * cos_tilt - pz * sin_tilt
    rz: float = py * sin_tilt + pz * cos_tilt + globe_center_z

    emit point(rx, ry, rz)


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

    globe_radius: float = globe_diameter / 2.0
    polar_radius: float = globe_radius * (1.0 - oblateness)

    lat_rad: float = deg2rad(latitude_deg)
    tilt_rad: float = deg2rad(tilt_deg)

    cos_lat: float = cos(lat_rad)
    sin_lat: float = sin(lat_rad)
    cos_tilt: float = cos(tilt_rad)
    sin_tilt: float = sin(tilt_rad)

    # Generate 32 points at 11.25° intervals using list comprehension
    # longitude = i * 11.25 for i in 0..31
    points: list<point> = [
        cradle_point(globe_radius, polar_radius, cos_lat, sin_lat,
                     cos_tilt, sin_tilt, beam_offset, globe_center_z,
                     i * 11.25)
        for i in range(32)
    ]

    # Build path segments connecting consecutive points (closing the loop)
    seg0: path3d = path3d_line(points[0], points[1])
    seg1: path3d = path3d_line(points[1], points[2])
    seg2: path3d = path3d_line(points[2], points[3])
    seg3: path3d = path3d_line(points[3], points[4])
    seg4: path3d = path3d_line(points[4], points[5])
    seg5: path3d = path3d_line(points[5], points[6])
    seg6: path3d = path3d_line(points[6], points[7])
    seg7: path3d = path3d_line(points[7], points[8])
    seg8: path3d = path3d_line(points[8], points[9])
    seg9: path3d = path3d_line(points[9], points[10])
    seg10: path3d = path3d_line(points[10], points[11])
    seg11: path3d = path3d_line(points[11], points[12])
    seg12: path3d = path3d_line(points[12], points[13])
    seg13: path3d = path3d_line(points[13], points[14])
    seg14: path3d = path3d_line(points[14], points[15])
    seg15: path3d = path3d_line(points[15], points[16])
    seg16: path3d = path3d_line(points[16], points[17])
    seg17: path3d = path3d_line(points[17], points[18])
    seg18: path3d = path3d_line(points[18], points[19])
    seg19: path3d = path3d_line(points[19], points[20])
    seg20: path3d = path3d_line(points[20], points[21])
    seg21: path3d = path3d_line(points[21], points[22])
    seg22: path3d = path3d_line(points[22], points[23])
    seg23: path3d = path3d_line(points[23], points[24])
    seg24: path3d = path3d_line(points[24], points[25])
    seg25: path3d = path3d_line(points[25], points[26])
    seg26: path3d = path3d_line(points[26], points[27])
    seg27: path3d = path3d_line(points[27], points[28])
    seg28: path3d = path3d_line(points[28], points[29])
    seg29: path3d = path3d_line(points[29], points[30])
    seg30: path3d = path3d_line(points[30], points[31])
    seg31: path3d = path3d_line(points[31], points[0])

    emit make_path3d(seg0, seg1, seg2, seg3, seg4, seg5, seg6, seg7,
                     seg8, seg9, seg10, seg11, seg12, seg13, seg14, seg15,
                     seg16, seg17, seg18, seg19, seg20, seg21, seg22, seg23,
                     seg24, seg25, seg26, seg27, seg28, seg29, seg30, seg31)


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

    emit sweep_adaptive_hollow(outer, inner, path, angle_threshold)


# ============================================================================
# Wrapped support arc with sigmoid blending
# ============================================================================

# Helper: compute a single intermediate point on the support arc
# with parabolic waist and sigmoid-blended wrapping around clearance ellipsoid
command support_arc_point(
    base_x: float, base_y: float, base_z: float,
    top_x: float, top_y: float, top_z: float,
    waist_ratio: float,
    a_clear: float, c_clear: float,
    cos_tilt: float, sin_tilt: float,
    globe_center_z: float,
    sigmoid_k: float,
    t: float
) -> point:
    # Parabolic factor peaks at t=0.5
    para: float = 4.0 * t * (1.0 - t)

    # Linear interpolation
    lin_x: float = base_x + (top_x - base_x) * t
    lin_y: float = base_y + (top_y - base_y) * t
    lin_z: float = base_z + (top_z - base_z) * t

    # Apply waist (radial contraction)
    lin_r: float = sqrt(lin_x * lin_x + lin_y * lin_y)
    waist_offset: float = (1.0 - waist_ratio) * lin_r * para
    scale: float = (lin_r - waist_offset) / lin_r
    raw_x: float = lin_x * scale
    raw_y: float = lin_y * scale
    raw_z: float = lin_z

    # Transform to local (untilted) frame to check against clearance ellipsoid
    loc_x: float = raw_x
    loc_y: float = raw_y * cos_tilt + (raw_z - globe_center_z) * sin_tilt
    loc_z: float = -raw_y * sin_tilt + (raw_z - globe_center_z) * cos_tilt

    # Check if inside clearance ellipsoid
    a_clear2: float = a_clear * a_clear
    c_clear2: float = c_clear * c_clear
    check: float = (loc_x*loc_x + loc_y*loc_y) / a_clear2 + loc_z*loc_z / c_clear2

    # Sigmoid blend factor: 1 when inside (check<1), 0 when outside
    d: float = check - 1.0
    s: float = 1.0 / (1.0 + exp(sigmoid_k * d))

    # Project to clearance surface
    proj_t: float = 1.0 / sqrt(check)
    proj_loc_x: float = loc_x * proj_t
    proj_loc_y: float = loc_y * proj_t
    proj_loc_z: float = loc_z * proj_t

    # Transform projected point back to world frame
    proj_world_x: float = proj_loc_x
    proj_world_y: float = proj_loc_y * cos_tilt - proj_loc_z * sin_tilt
    proj_world_z: float = proj_loc_y * sin_tilt + proj_loc_z * cos_tilt + globe_center_z

    # Blend between raw and projected based on sigmoid
    x: float = (1.0 - s) * raw_x + s * proj_world_x
    y: float = (1.0 - s) * raw_y + s * proj_world_y
    z: float = (1.0 - s) * raw_z + s * proj_world_z

    emit point(x, y, z)


# High-resolution wrapped support arc path with 16 intermediate points
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

    globe_radius: float = globe_diameter / 2.0
    polar_radius: float = globe_radius * (1.0 - oblateness)

    # Clearance ellipsoid = globe + beam_offset
    a_clear: float = globe_radius + beam_offset
    c_clear: float = polar_radius + beam_offset

    sigmoid_k: float = 20.0

    lat_rad: float = deg2rad(latitude_deg)
    tilt_rad: float = deg2rad(tilt_deg)
    base_angle_rad: float = deg2rad(base_angle_deg)
    top_lon_rad: float = deg2rad(top_longitude_deg)

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

    # Top point: on latitude circle with beam offset
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

    # Generate 16 intermediate points using list comprehension
    # t values: 1/17, 2/17, ..., 16/17
    arc_points: list<point> = [
        support_arc_point(base_x, base_y, base_z, top_x, top_y, top_z,
                          waist_ratio, a_clear, c_clear,
                          cos_tilt, sin_tilt, globe_center_z, sigmoid_k,
                          i * 1.0 / 17.0)
        for i in range(1, 17)
    ]

    # Build 17 path segments: base -> p1 -> p2 -> ... -> p16 -> top
    base_pt: point = point(base_x, base_y, base_z)
    top_pt: point = point(top_x, top_y, top_z)

    seg0: path3d = path3d_line(base_pt, arc_points[0])
    seg1: path3d = path3d_line(arc_points[0], arc_points[1])
    seg2: path3d = path3d_line(arc_points[1], arc_points[2])
    seg3: path3d = path3d_line(arc_points[2], arc_points[3])
    seg4: path3d = path3d_line(arc_points[3], arc_points[4])
    seg5: path3d = path3d_line(arc_points[4], arc_points[5])
    seg6: path3d = path3d_line(arc_points[5], arc_points[6])
    seg7: path3d = path3d_line(arc_points[6], arc_points[7])
    seg8: path3d = path3d_line(arc_points[7], arc_points[8])
    seg9: path3d = path3d_line(arc_points[8], arc_points[9])
    seg10: path3d = path3d_line(arc_points[9], arc_points[10])
    seg11: path3d = path3d_line(arc_points[10], arc_points[11])
    seg12: path3d = path3d_line(arc_points[11], arc_points[12])
    seg13: path3d = path3d_line(arc_points[12], arc_points[13])
    seg14: path3d = path3d_line(arc_points[13], arc_points[14])
    seg15: path3d = path3d_line(arc_points[14], arc_points[15])
    seg16: path3d = path3d_line(arc_points[15], top_pt)

    emit make_path3d(seg0, seg1, seg2, seg3, seg4, seg5, seg6, seg7,
                     seg8, seg9, seg10, seg11, seg12, seg13, seg14, seg15, seg16)


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

    emit sweep_adaptive_hollow(outer, inner, path, angle_threshold)


# ============================================================================
# Mars globe
# ============================================================================

command centered_mars_globe(
    globe_diameter: float,
    globe_center_z: float,
    tilt_angle: float,
    oblateness: float = 0.00648
) -> solid:
    globe: solid = oblate_spheroid(globe_diameter, oblateness)
    tilted_globe: solid = rotate(globe, tilt_angle, 0.0, 0.0)
    positioned_globe: solid = translate(tilted_globe, 0.0, 0.0, globe_center_z + 5.0)
    emit positioned_globe


# ============================================================================
# Main command: CENTERED_GLOBE_STAND_HIRES
# ============================================================================

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
    mars_oblateness: float = 0.00648,
    include_mars: bool = false
) -> solid:
    base_radius: float = base_diameter / 2.0
    globe_radius: float = globe_diameter / 2.0
    polar_radius: float = globe_radius * (1.0 - mars_oblateness)
    beam_offset: float = beam_outer / 2.0

    lat_rad: float = deg2rad(cradle_latitude)
    tilt_rad: float = deg2rad(tilt_angle)

    cos_lat: float = cos(lat_rad)
    sin_lat: float = sin(lat_rad)
    cos_tilt: float = cos(tilt_rad)
    sin_tilt: float = sin(tilt_rad)

    a2: float = globe_radius * globe_radius
    c2: float = polar_radius * polar_radius

    # Calculate globe center Z so lowest cradle point is at target height
    # Point at 270° longitude (negative Y) is the lowest
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

    # Build components
    base: solid = base_ring(base_radius, beam_outer, beam_wall, angle_threshold)

    cradle: solid = latitude_cradle_ring_hires(
        globe_diameter, mars_oblateness, cradle_latitude, tilt_angle,
        globe_center_z, beam_outer, beam_wall, angle_threshold
    )

    arc1: solid = wrapped_support_arc_hires(
        base_radius, globe_diameter, mars_oblateness, cradle_latitude,
        tilt_angle, globe_center_z, waist_ratio,
        0.0, twist_deg, beam_outer, beam_wall, angle_threshold
    )

    arc2: solid = wrapped_support_arc_hires(
        base_radius, globe_diameter, mars_oblateness, cradle_latitude,
        tilt_angle, globe_center_z, waist_ratio,
        120.0, 120.0 + twist_deg, beam_outer, beam_wall, angle_threshold
    )

    arc3: solid = wrapped_support_arc_hires(
        base_radius, globe_diameter, mars_oblateness, cradle_latitude,
        tilt_angle, globe_center_z, waist_ratio,
        240.0, 240.0 + twist_deg, beam_outer, beam_wall, angle_threshold
    )

    mars: solid = centered_mars_globe(globe_diameter, globe_center_z, tilt_angle, mars_oblateness)

    stand: solid = union(base, cradle, arc1, arc2, arc3)
    result: solid = compound(stand, mars) if include_mars else stand

    emit result
