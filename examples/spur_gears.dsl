module spur_gears

# Create a single involute spur gear
command MAKE_GEAR(
    teeth: int,
    module_mm: float,
    pressure_angle: float,
    face_width: float,
    hub_diameter: float,
    hub_height: float
) -> solid:
    # Create the gear body using proper involute profile
    gear_body: solid = involute_gear(teeth, module_mm, pressure_angle, face_width)

    # Calculate pitch diameter for hub positioning
    pitch_diameter: float = teeth * module_mm

    # Create a hub if specified
    hub: solid = cylinder(hub_diameter / 2.0, face_width + 2.0 * hub_height)
    hub_offset: float = -hub_height
    positioned_hub: solid = translate(hub, 0.0, 0.0, hub_offset)

    # Union gear body with hub
    final_gear: solid = union(gear_body, positioned_hub)

    emit final_gear


# Create a meshed gear pair
command MAKE_GEAR_PAIR(
    teeth1: int,
    teeth2: int,
    module_mm: float,
    pressure_angle: float,
    face_width: float,
    hub_diameter: float,
    hub_height: float
) -> solid:
    # Create both gears with proper involute profiles
    gear1: solid = involute_gear(teeth1, module_mm, pressure_angle, face_width)
    gear2: solid = involute_gear(teeth2, module_mm, pressure_angle, face_width)

    # Calculate center distance for meshing (sum of pitch radii)
    pitch_radius1: float = (teeth1 * module_mm) / 2.0
    pitch_radius2: float = (teeth2 * module_mm) / 2.0
    center_distance: float = pitch_radius1 + pitch_radius2

    # Rotate gear2 by half a tooth pitch to mesh properly
    # mesh_rotation: float = 180.0 / teeth2
    mesh_rotation: float = 360.0 / teeth2
    gear2_rotated: solid = rotate(gear2, 0.0, 0.0, mesh_rotation * 0.85)

    # Position gear2 at center distance
    gear2_positioned: solid = translate(gear2_rotated, center_distance, 0.0, 0.0)

    # Combine both gears
    gear_pair: solid = union(gear1, gear2_positioned)

    emit gear_pair


# Simple demo with a single gear
command DEMO_SINGLE_GEAR() -> solid:
    # 24-tooth gear, 2mm module, 20 degree pressure angle, 10mm face width
    gear: solid = involute_gear(24, 2.0, 20.0, 10.0)
    emit gear
