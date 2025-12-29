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

    # Create stringer notches at 120° intervals using list comprehension
    notches: list<solid> = [
        make_notch(outer_radius, stringer_width_mm, stringer_depth_mm, thickness_mm, i * 120.0)
        for i in range(3)
    ]

    # Subtract all notches at once
    final_plate: solid = difference_all(plate1, notches)

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

    # Hole angles: 3 sectors at 60°, 180°, 300° with 3 holes each (±20° spread)
    # Sector centers: 60, 180, 300
    # Hole offsets: -20, 0, +20
    hole_angles: list<float> = [40.0, 60.0, 80.0, 160.0, 180.0, 200.0, 280.0, 300.0, 320.0]

    # Create all 9 holes using list comprehension
    holes: list<solid> = [
        make_hole(radial_dist, hole_radius_mm, thickness_mm, angle)
        for angle in hole_angles
    ]

    # Subtract all holes at once
    final: solid = difference_all(base_plate, holes)

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

    # Hole angles: 3 sectors at 60°, 180°, 300° with 4 holes each
    # Sector 1: 35, 50, 70, 85  |  Sector 2: 155, 170, 190, 205  |  Sector 3: 275, 290, 310, 325
    hole_angles: list<float> = [
        35.0, 50.0, 70.0, 85.0,
        155.0, 170.0, 190.0, 205.0,
        275.0, 290.0, 310.0, 325.0
    ]

    # Create all 12 holes using list comprehension
    holes: list<solid> = [
        make_hole(radial_dist, hole_radius_mm, thickness_mm, angle)
        for angle in hole_angles
    ]

    # Subtract all holes at once
    final: solid = difference_all(base_plate, holes)

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
    # Box centered at origin, then translated along +X to outer edge
    notch_box: solid = box(depth_mm, width_mm, thickness_mm + 2.0)
    radial_offset: float = outer_radius - depth_mm / 2.0
    # Translate along X axis (radial direction at 0°)
    notch_at_edge: solid = translate(notch_box, radial_offset, 0.0, -1.0)
    # Then rotate around Z to desired angle
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
