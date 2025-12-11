module thrust_structure

# Thrust Structure for Liquid Fuel Rocket
# ========================================
#
# A 12" diameter aluminum thrust plate with:
# - Central motor mount hole (4" diameter)
# - Three stringer notches at 120° intervals
# - Optional lightening holes for mass optimization
#
# Design constraints:
# - Must support 1500 lbf thrust
# - Max deflection < 10mm
# - Trilateral symmetry enforced
#
# Material: 6061-T6 Aluminum
# - Density: 2700 kg/m³
# - Young's Modulus: 68.9 GPa

# =============================================================================
# Main Command: Create the thrust plate without lightening holes (baseline)
# =============================================================================

command MAKE_BASELINE_PLATE(
    outer_diameter_mm: float = 304.8,
    thickness_mm: float = 8.0,
    motor_mount_diameter_mm: float = 101.6,
    stringer_width_mm: float = 50.8,
    stringer_depth_mm: float = 25.4
) -> solid:
    require outer_diameter_mm > motor_mount_diameter_mm, "Outer diameter must exceed motor mount"
    require thickness_mm > 0.0, "Thickness must be positive"

    outer_radius: float = outer_diameter_mm / 2.0
    motor_radius: float = motor_mount_diameter_mm / 2.0

    # Create base disk
    base_disk: solid = cylinder(outer_radius, thickness_mm)

    # Create motor mount hole
    motor_hole: solid = cylinder(motor_radius, thickness_mm + 2.0)
    motor_hole_pos: solid = translate(motor_hole, 0.0, 0.0, -1.0)

    # Subtract motor mount
    plate1: solid = difference(base_disk, motor_hole_pos)

    # Create and subtract notch 1 at 0°
    notch1: solid = make_notch(outer_radius, stringer_width_mm, stringer_depth_mm, thickness_mm, 0.0)
    plate2: solid = difference(plate1, notch1)

    # Create and subtract notch 2 at 120°
    notch2: solid = make_notch(outer_radius, stringer_width_mm, stringer_depth_mm, thickness_mm, 120.0)
    plate3: solid = difference(plate2, notch2)

    # Create and subtract notch 3 at 240°
    notch3: solid = make_notch(outer_radius, stringer_width_mm, stringer_depth_mm, thickness_mm, 240.0)
    final_plate: solid = difference(plate3, notch3)

    emit final_plate


# =============================================================================
# Command: Create plate with 3 lightening holes per sector (9 total)
# =============================================================================

command MAKE_LIGHTENED_PLATE_3(
    outer_diameter_mm: float = 304.8,
    thickness_mm: float = 8.0,
    motor_mount_diameter_mm: float = 101.6,
    stringer_width_mm: float = 50.8,
    stringer_depth_mm: float = 25.4,
    hole_radius_mm: float = 20.0,
    hole_radial_fraction: float = 0.6
) -> solid:
    # Start with baseline
    base_plate: solid = MAKE_BASELINE_PLATE(
        outer_diameter_mm, thickness_mm, motor_mount_diameter_mm,
        stringer_width_mm, stringer_depth_mm
    )

    outer_radius: float = outer_diameter_mm / 2.0
    motor_radius: float = motor_mount_diameter_mm / 2.0

    # Calculate radial position for holes
    inner_bound: float = motor_radius + hole_radius_mm + 5.0
    outer_bound: float = outer_radius - stringer_depth_mm - hole_radius_mm - 5.0
    radial_dist: float = inner_bound + hole_radial_fraction * (outer_bound - inner_bound)

    # Sector 1 (centered at 60°): holes at 40°, 60°, 80°
    h1_1: solid = make_hole(radial_dist, hole_radius_mm, thickness_mm, 40.0)
    h1_2: solid = make_hole(radial_dist, hole_radius_mm, thickness_mm, 60.0)
    h1_3: solid = make_hole(radial_dist, hole_radius_mm, thickness_mm, 80.0)

    # Sector 2 (centered at 180°): holes at 160°, 180°, 200°
    h2_1: solid = make_hole(radial_dist, hole_radius_mm, thickness_mm, 160.0)
    h2_2: solid = make_hole(radial_dist, hole_radius_mm, thickness_mm, 180.0)
    h2_3: solid = make_hole(radial_dist, hole_radius_mm, thickness_mm, 200.0)

    # Sector 3 (centered at 300°): holes at 280°, 300°, 320°
    h3_1: solid = make_hole(radial_dist, hole_radius_mm, thickness_mm, 280.0)
    h3_2: solid = make_hole(radial_dist, hole_radius_mm, thickness_mm, 300.0)
    h3_3: solid = make_hole(radial_dist, hole_radius_mm, thickness_mm, 320.0)

    # Subtract all holes
    p1: solid = difference(base_plate, h1_1)
    p2: solid = difference(p1, h1_2)
    p3: solid = difference(p2, h1_3)
    p4: solid = difference(p3, h2_1)
    p5: solid = difference(p4, h2_2)
    p6: solid = difference(p5, h2_3)
    p7: solid = difference(p6, h3_1)
    p8: solid = difference(p7, h3_2)
    final: solid = difference(p8, h3_3)

    emit final


# =============================================================================
# Command: Create plate with 4 lightening holes per sector (12 total)
# =============================================================================

command MAKE_LIGHTENED_PLATE_4(
    outer_diameter_mm: float = 304.8,
    thickness_mm: float = 8.0,
    motor_mount_diameter_mm: float = 101.6,
    stringer_width_mm: float = 50.8,
    stringer_depth_mm: float = 25.4,
    hole_radius_mm: float = 20.0,
    hole_radial_fraction: float = 0.6
) -> solid:
    base_plate: solid = MAKE_BASELINE_PLATE(
        outer_diameter_mm, thickness_mm, motor_mount_diameter_mm,
        stringer_width_mm, stringer_depth_mm
    )

    outer_radius: float = outer_diameter_mm / 2.0
    motor_radius: float = motor_mount_diameter_mm / 2.0

    inner_bound: float = motor_radius + hole_radius_mm + 5.0
    outer_bound: float = outer_radius - stringer_depth_mm - hole_radius_mm - 5.0
    radial_dist: float = inner_bound + hole_radial_fraction * (outer_bound - inner_bound)

    # Sector 1 (60°): holes at 35°, 50°, 70°, 85°
    h1_1: solid = make_hole(radial_dist, hole_radius_mm, thickness_mm, 35.0)
    h1_2: solid = make_hole(radial_dist, hole_radius_mm, thickness_mm, 50.0)
    h1_3: solid = make_hole(radial_dist, hole_radius_mm, thickness_mm, 70.0)
    h1_4: solid = make_hole(radial_dist, hole_radius_mm, thickness_mm, 85.0)

    # Sector 2 (180°): holes at 155°, 170°, 190°, 205°
    h2_1: solid = make_hole(radial_dist, hole_radius_mm, thickness_mm, 155.0)
    h2_2: solid = make_hole(radial_dist, hole_radius_mm, thickness_mm, 170.0)
    h2_3: solid = make_hole(radial_dist, hole_radius_mm, thickness_mm, 190.0)
    h2_4: solid = make_hole(radial_dist, hole_radius_mm, thickness_mm, 205.0)

    # Sector 3 (300°): holes at 275°, 290°, 310°, 325°
    h3_1: solid = make_hole(radial_dist, hole_radius_mm, thickness_mm, 275.0)
    h3_2: solid = make_hole(radial_dist, hole_radius_mm, thickness_mm, 290.0)
    h3_3: solid = make_hole(radial_dist, hole_radius_mm, thickness_mm, 310.0)
    h3_4: solid = make_hole(radial_dist, hole_radius_mm, thickness_mm, 325.0)

    # Subtract all 12 holes
    p1: solid = difference(base_plate, h1_1)
    p2: solid = difference(p1, h1_2)
    p3: solid = difference(p2, h1_3)
    p4: solid = difference(p3, h1_4)
    p5: solid = difference(p4, h2_1)
    p6: solid = difference(p5, h2_2)
    p7: solid = difference(p6, h2_3)
    p8: solid = difference(p7, h2_4)
    p9: solid = difference(p8, h3_1)
    p10: solid = difference(p9, h3_2)
    p11: solid = difference(p10, h3_3)
    final: solid = difference(p11, h3_4)

    emit final


# =============================================================================
# Command: Optimized configuration (found via optimization)
# =============================================================================

command MAKE_OPTIMIZED_PLATE() -> solid:
    # Optimized parameters for minimum mass meeting 10mm deflection limit
    plate: solid = MAKE_LIGHTENED_PLATE_4(
        304.8,   # outer diameter mm
        8.0,     # thickness mm
        101.6,   # motor mount diameter mm
        50.8,    # stringer width mm
        25.4,    # stringer depth mm
        25.0,    # hole radius mm (larger holes = more weight savings)
        0.65     # radial fraction (position holes)
    )
    emit plate


# =============================================================================
# Helper: Create a single rectangular stringer notch
# =============================================================================

command make_notch(
    outer_radius: float,
    width_mm: float,
    depth_mm: float,
    thickness_mm: float,
    angle_deg: float
) -> solid:
    notch_box: solid = box(width_mm, depth_mm, thickness_mm + 2.0)
    radial_offset: float = outer_radius - depth_mm / 2.0
    notch_at_edge: solid = translate(notch_box, 0.0, radial_offset, -1.0)
    notch_final: solid = rotate(notch_at_edge, 0.0, 0.0, angle_deg)
    emit notch_final


# =============================================================================
# Helper: Create a single cylindrical lightening hole
# =============================================================================

command make_hole(
    radial_distance: float,
    radius: float,
    thickness_mm: float,
    angle_deg: float
) -> solid:
    hole_cyl: solid = cylinder(radius, thickness_mm + 2.0)
    hole_at_radius: solid = translate(hole_cyl, radial_distance, 0.0, -1.0)
    hole_final: solid = rotate(hole_at_radius, 0.0, 0.0, angle_deg)
    emit hole_final
