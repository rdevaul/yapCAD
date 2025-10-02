"""Rocket internal layout demo using yapCAD parametric primitives.

Builds a single-stage liquid rocket cross-section with fuel/oxidizer tanks,
feed lines, turbopumps, and a simplified engine.  Outputs STL and STEP files
for inspection.  Geometry is generated headlessly; set ``YAPCAD_OUTPUT_DIR``
to control where artifacts land (defaults to the examples directory).
"""

from __future__ import annotations

import math
import os
from pathlib import Path

from yapcad.geom import epsilon, point, vect
from yapcad.geom3d import solid, solidbbox
from yapcad.geom3d_util import (
    conic_tube,
    rotatesolid,
    spherical_shell,
    stack_solids,
    translatesolid,
    tube,
)
from yapcad.io import write_step, write_stl


OUTPUT_DIR = Path(os.environ.get('YAPCAD_OUTPUT_DIR', 'examples'))
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def merge_solids(solids):
    surfaces = []
    for sld in solids:
        surfaces.extend(sld[1])
    return solid(surfaces, [], ['procedure', 'rocket_internal_demo'])


def cutaway_solid_x(sld, plane_x=0.0, eps=1e-6):
    """Return a copy of ``sld`` with triangles strictly on the +X side removed."""

    stripped = []
    for surf in sld[1]:
        faces = []
        for tri in surf[3]:
            xs = [surf[1][idx][0] for idx in tri]
            if max(xs) <= plane_x + eps:
                faces.append(tri)
        if faces:
            stripped.append(['surface', surf[1], surf[2], faces, surf[4], surf[5]])
    return solid(stripped, [], sld[3] if len(sld) > 3 else [])


def make_tank(z_start, length, outer_radius, wall):
    cylinder_length = max(length - 2 * outer_radius, wall)
    cylinder_base = z_start + outer_radius
    body = tube(outer_diameter=outer_radius * 2,
                wall_thickness=wall,
                length=cylinder_length,
                base_point=point(0, 0, cylinder_base),
                include_caps=False)

    front_cap_center = point(0, 0, z_start + outer_radius)
    front_cap = spherical_shell(outer_diameter=outer_radius * 2,
                                wall_thickness=wall,
                                solid_angle=2 * math.pi,
                                center=front_cap_center)
    front_cap = rotatesolid(front_cap, 180, axis=point(0, 1, 0), cent=front_cap_center)

    back_cap_center = point(0, 0, z_start + outer_radius + cylinder_length)
    back_cap = spherical_shell(outer_diameter=outer_radius * 2,
                               wall_thickness=wall,
                               solid_angle=2 * math.pi,
                               center=back_cap_center)

    return merge_solids([body, front_cap, back_cap])


def make_engine(z_start):
    wall = 0.25
    chamber_d = 4.0
    chamber_l = 2.5
    throat_l = 1.5
    nozzle_l = 5.0

    nozzle = conic_tube(bottom_outer_diameter=6.0,
                        top_outer_diameter=2.2,
                        wall_thickness=wall,
                        length=nozzle_l,
                        base_point=point(0, 0, z_start),
                        include_caps=False)

    throat = tube(outer_diameter=2.2,
                  wall_thickness=wall,
                  length=throat_l,
                  base_point=point(0, 0, z_start + nozzle_l),
                  include_caps=False)

    chamber = conic_tube(bottom_outer_diameter=2.2,
                         top_outer_diameter=chamber_d,
                         wall_thickness=wall,
                         length=chamber_l,
                         base_point=point(0, 0, z_start + nozzle_l + throat_l),
                         include_caps=False)

    cap_center = z_start + nozzle_l + throat_l + chamber_l
    front_cap = spherical_shell(outer_diameter=chamber_d,
                                wall_thickness=wall,
                                solid_angle=2 * math.pi,
                                center=point(0, 0, cap_center))

    return merge_solids([nozzle, throat, chamber, front_cap])


def build_body_shell(cylinder_length, radius, wall):
    body = tube(outer_diameter=radius * 2,
                wall_thickness=wall,
                length=cylinder_length,
                center=point(0, 0, 0),
#                center=point(0, 0, cylinder_length / 2.0),
                include_caps=False)
    nose_center = point(0, 0, cylinder_length )
    nose = spherical_shell(outer_diameter=radius * 2,
                           wall_thickness=wall,
                           solid_angle=2 * math.pi,
                           center=nose_center)
    shell = merge_solids([body, nose])
    shell = cutaway_solid_x(shell, plane_x=0.0)
    bbox = solidbbox(shell)
    if bbox[0][2] > 0.0:
        shell = translatesolid(shell, vect(0, 0, -bbox[0][2], 0))
    return shell


def main():
    body_radius = 6.0
    wall = 0.3
    cylinder_length = 32.0
    body_shell = build_body_shell(cylinder_length, body_radius, wall)

    engine = make_engine(0.0)
    fuel_tank = make_tank(0.0, length=8.0, outer_radius=3.0, wall=0.25)
    oxidizer_tank = make_tank(0.0, length=9.0, outer_radius=3.8, wall=0.25)

    engine, oxidizer_tank, fuel_tank = stack_solids(
        [engine, 'space:2.0', oxidizer_tank, 'space:2.0', fuel_tank],
        axis='z',
        start=0.0,
        gap=0.0,
        align='center',
    )

    min_z = min(solidbbox(obj)[0][2] for obj in (engine, fuel_tank, oxidizer_tank))
    if abs(min_z) > epsilon:
        offset = -min_z
        engine = translatesolid(engine, vect(0, 0, offset, 0))
        fuel_tank = translatesolid(fuel_tank, vect(0, 0, offset, 0))
        oxidizer_tank = translatesolid(oxidizer_tank, vect(0, 0, offset, 0))

    top_target = cylinder_length - 2.0
    current_top = max(solidbbox(obj)[1][2] for obj in (fuel_tank, oxidizer_tank))
    adjust = top_target - current_top
    engine = translatesolid(engine, vect(0, 0, adjust, 0))
    fuel_tank = translatesolid(fuel_tank, vect(0, 0, adjust, 0))
    oxidizer_tank = translatesolid(oxidizer_tank, vect(0, 0, adjust, 0))

    min_z = min(solidbbox(obj)[0][2] for obj in (engine, fuel_tank, oxidizer_tank))
    if min_z < 0.0:
        engine = translatesolid(engine, vect(0, 0, -min_z, 0))
        fuel_tank = translatesolid(fuel_tank, vect(0, 0, -min_z, 0))
        oxidizer_tank = translatesolid(oxidizer_tank, vect(0, 0, -min_z, 0))

    engine = translatesolid(engine, vect(0, 0, -4.0, 0))

    components = [body_shell, engine, fuel_tank, oxidizer_tank]

    assembly = merge_solids(components)
    bbox = solidbbox(assembly)
    print('Bounding box:', bbox)

    stl_path = OUTPUT_DIR / 'rocket_internal_demo.stl'
    step_path = OUTPUT_DIR / 'rocket_internal_demo.step'
    write_stl(assembly, stl_path, binary=True, name='rocket_internal_demo')
    write_step(assembly, step_path, name='rocket_internal_demo')
    print('Exported STL to', stl_path)
    print('Exported STEP to', step_path)


if __name__ == '__main__':
    main()
