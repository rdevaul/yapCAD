# Involute gear generation using native Python block
# This demonstrates the native block feature for wrapping Python libraries

module figgear_demo;

# Native block wrapping the figgear module
native python {
from yapcad.contrib.figgear import make_gear_figure
from yapcad.geom3d import poly2surfaceXY
from yapcad.geom3d_util import extrude
from yapcad.geom import point

def create_gear_profile(teeth, module_mm, pressure_angle):
    """Create a 2D gear profile using figgear."""
    profile_points, blueprints = make_gear_figure(
        m=module_mm,
        z=int(teeth),
        alpha_deg=pressure_angle,
        bottom_type="spline"
    )
    # Convert to yapCAD points with z=0 for XY plane
    pts = [point(x, y, 0.0) for x, y in profile_points]
    return pts

def create_gear_solid(teeth, module_mm, pressure_angle, face_width):
    """Create a 3D gear solid by extruding the profile."""
    pts = create_gear_profile(teeth, module_mm, pressure_angle)
    # Convert points to surface in XY plane (returns tuple: surface, info)
    surface, _ = poly2surfaceXY(pts)
    # Extrude the surface into a solid
    return extrude(surface, face_width)

} exports {
    fn create_gear_profile(teeth: int, module_mm: float, pressure_angle: float) -> region2d;
    fn create_gear_solid(teeth: int, module_mm: float, pressure_angle: float, face_width: float) -> solid;
}


# Create a single gear using the native figgear wrapper
command MAKE_FIGGEAR(
    teeth: int,
    module_mm: float,
    pressure_angle: float,
    face_width: float
) -> solid {
    # Use the native function to create the gear solid
    let gear: solid = create_gear_solid(teeth, module_mm, pressure_angle, face_width);
    emit gear;
}


# Create a gear with a center bore
command MAKE_BORED_GEAR(
    teeth: int,
    module_mm: float,
    pressure_angle: float,
    face_width: float,
    bore_diameter: float
) -> solid {
    # Create the gear body
    let gear_body: solid = create_gear_solid(teeth, module_mm, pressure_angle, face_width);

    # Create bore cylinder
    let bore: solid = cylinder(bore_diameter / 2.0, face_width + 0.1);
    let bore_centered: solid = translate(bore, 0.0, 0.0, -0.05);

    # Subtract bore from gear
    let final_gear: solid = difference(gear_body, bore_centered);
    emit final_gear;
}


# Create meshed gear pair using figgear profiles
command MAKE_FIGGEAR_PAIR(
    teeth1: int,
    teeth2: int,
    module_mm: float,
    pressure_angle: float,
    face_width: float
) -> solid {
    # Create both gears
    let gear1: solid = create_gear_solid(teeth1, module_mm, pressure_angle, face_width);
    let gear2: solid = create_gear_solid(teeth2, module_mm, pressure_angle, face_width);

    # Calculate center distance for meshing
    let pitch_radius1: float = (teeth1 * module_mm) / 2.0;
    let pitch_radius2: float = (teeth2 * module_mm) / 2.0;
    let center_distance: float = pitch_radius1 + pitch_radius2;

    # Rotate gear2 to mesh properly (half tooth pitch)
    let mesh_rotation: float = 180.0 / teeth2;
    let gear2_rotated: solid = rotate(gear2, 0.0, 0.0, mesh_rotation);

    # Position gear2
    let gear2_positioned: solid = translate(gear2_rotated, center_distance, 0.0, 0.0);

    # Combine gears
    let gear_pair: solid = union(gear1, gear2_positioned);
    emit gear_pair;
}


# Demo command showing a simple gear
command DEMO() -> solid {
    # 20-tooth gear, 2mm module, 20 degree pressure angle, 8mm face width
    let gear: solid = create_gear_solid(20, 2.0, 20.0, 8.0);
    emit gear;
}
