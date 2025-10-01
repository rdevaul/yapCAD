"""Advanced rocket demo featuring grid fins and exploded view."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple

from yapcad.geom import point, vect
from yapcad.geom3d import solid, surface, translatesolid, rotatesolid
from yapcad.geom3d_util import conic, extrude
from yapcad.io import write_stl, write_step


def merge_solids(parts: List[list]) -> list:
    surfaces: List = []
    for part in parts:
        surfaces.extend(part[1])
    return solid(surfaces)


def stack_conic(radius_base: float, radius_top: float, height: float, base_z: float) -> list:
    return conic(radius_base, radius_top, height, center=point(0, 0, base_z))


def create_engine(x: float, y: float, base_z: float, bell_radius: float, height: float) -> list:
    bell = conic(bell_radius, bell_radius * 0.4, height, center=point(x, y, base_z))
    nozzle = conic(bell_radius * 0.4, bell_radius * 0.4, height * 0.2, center=point(x, y, base_z + height))
    return merge_solids([bell, nozzle])


def create_box(width: float, depth: float, height: float, center: Tuple[float, float, float]) -> list:
    cx, cy, cz = center
    hx, hy, hz = width / 2.0, depth / 2.0, height / 2.0
    base_z = cz - hz
    poly = [
        point(cx - hx, cy - hy, base_z),
        point(cx + hx, cy - hy, base_z),
        point(cx + hx, cy + hy, base_z),
        point(cx - hx, cy + hy, base_z),
    ]
    normals = [[0, 0, 1, 0] for _ in poly]
    faces = [[0, 1, 2], [0, 2, 3]]
    boundary = [0, 1, 2, 3]
    base_surface = surface(poly, normals, faces, boundary, [])
    return extrude(base_surface, height, direction=vect(0, 0, 1, 0))


def create_grid_fin_panel(width: float, height: float, thickness: float, cells: int) -> list:
    bars: List[list] = []
    half_w, half_h = width / 2.0, height / 2.0

    frame = thickness

    # Outer frame top/bottom
    bars.append(create_box(thickness, width, thickness, (0.0, 0.0, half_h - thickness / 2.0)))
    bars.append(create_box(thickness, width, thickness, (0.0, 0.0, -half_h + thickness / 2.0)))

    # Outer frame sides
    bars.append(create_box(thickness, thickness, height, (0.0, -half_w + thickness / 2.0, 0.0)))
    bars.append(create_box(thickness, thickness, height, (0.0, half_w - thickness / 2.0, 0.0)))

    # Internal horizontal bars
    if cells > 1:
        step = height / cells
        for i in range(1, cells):
            z = -half_h + i * step
            bars.append(create_box(thickness, width - 2 * frame, thickness / 1.5, (0.0, 0.0, z)))

    # Internal vertical bars
    if cells > 1:
        step = width / cells
        for i in range(1, cells):
            y = -half_w + i * step
            bars.append(create_box(thickness, thickness / 1.5, height - 2 * frame, (0.0, y, 0.0)))

    return merge_solids(bars)


def build_components() -> Dict[str, list]:
    components: Dict[str, list] = {}

    booster_height = 18.0
    booster_radius = 2.6
    interstage_height = 3.0
    stage2_height = 8.0
    fairing_lower_height = 2.0
    fairing_upper_height = 2.5

    booster_body = stack_conic(booster_radius, booster_radius, booster_height, 0.0)

    # Engine cluster
    engines = [create_engine(0.0, 0.0, -2.5, 0.55, 2.2)]
    ring = 1.4
    for dx, dy in [(ring, 0.0), (-ring, 0.0), (0.0, ring), (0.0, -ring)]:
        engines.append(create_engine(dx, dy, -2.5, 0.45, 1.8))
    engine_cluster = merge_solids(engines)

    # Grid fins near top of booster
    fin_panel = create_grid_fin_panel(width=3.4, height=3.0, thickness=0.35, cells=4)
    fin_z = 14.0
    fin_offset = booster_radius + 0.25
    grid_fins: List[list] = []
    translated = translatesolid(fin_panel, point(fin_offset, 0.0, fin_z))
    grid_fins.append(translated)
    for angle in (90, 180, 270):
        grid_fins.append(rotatesolid(translated, angle, cent=point(0, 0, fin_z), axis=point(0, 0, 1)))
    grid_fin_cluster = merge_solids(grid_fins)

    booster = merge_solids([booster_body, engine_cluster, grid_fin_cluster])
    components['booster'] = booster

    interstage = stack_conic(booster_radius, 1.9, interstage_height, booster_height)
    components['interstage'] = interstage

    stage2 = stack_conic(1.9, 1.9, stage2_height, booster_height + interstage_height)
    upper_engine = create_engine(0.0, 0.0, booster_height + interstage_height - 1.5, 0.4, 1.6)
    payload_barrel = stack_conic(1.9, 1.9, 3.0, booster_height + interstage_height + stage2_height)
    fairing_lower = stack_conic(1.9, 1.1, fairing_lower_height, booster_height + interstage_height + stage2_height + 3.0)
    fairing_upper = stack_conic(1.1, 0.1, fairing_upper_height, booster_height + interstage_height + stage2_height + 3.0 + fairing_lower_height)

    upper_stage = merge_solids([stage2, payload_barrel, fairing_lower, fairing_upper, upper_engine])
    components['upper_stage'] = upper_stage

    return components


def build_assembly(components: Dict[str, list]) -> list:
    return merge_solids(list(components.values()))


def visualise_exploded(components: Dict[str, list]) -> None:
    try:
        from yapcad.pyglet_drawable import pygletDraw
    except Exception as exc:  # pragma: no cover - visual only
        print(f'Visualisation skipped: {exc}')
        return

    offsets = {
        'booster': 0.0,
        'interstage': 5.0,
        'upper_stage': 10.0,
    }

    palette = {
        'booster': 'silver',
        'interstage': 'olive',
        'upper_stage': 'navy',
    }

    dd = pygletDraw()
    for name, solid_obj in components.items():
        offset = offsets.get(name, 0.0)
        display_obj = translatesolid(solid_obj, point(0, 0, offset))
        dd.make_object(name, lighting=True, material='pearl')
        dd.linecolor = palette.get(name, 'gray')
        dd.draw(display_obj, name=name)
    dd.display()


def export_models(assembly: list, stem: str) -> None:
    stl_path = Path(f'{stem}.stl')
    step_path = Path(f'{stem}.step')
    write_stl(assembly, stl_path, binary=True, name=stem)
    write_step(assembly, step_path, name=stem)
    print(f'Exported STL: {stl_path.resolve()}')
    print(f'Exported STEP: {step_path.resolve()}')


def main() -> None:
    components = build_components()
    assembly = build_assembly(components)
    export_models(assembly, 'rocket_grid_demo')
    visualise_exploded(components)


if __name__ == '__main__':
    main()
