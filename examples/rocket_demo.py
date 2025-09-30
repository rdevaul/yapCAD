"""Simple generative rocket demo using yapCAD.

Builds a multi-stage rocket model, visualises it using pyglet, and exports
an STL file suitable for downstream slicers.
"""

from __future__ import annotations

from pathlib import Path

from yapcad.geom import point, vect
from yapcad.geom3d import solid, surface, translatesolid, rotatesolid
from yapcad.geom3d_util import conic, extrude
from yapcad.io.stl import write_stl


def stack_conic(baser: float, topr: float, height: float, z0: float) -> list:
    """Create a conic segment whose base sits at ``z0``."""
    return conic(baser, topr, height, center=point(0, 0, z0))


def create_engine(x: float, y: float) -> list:
    """Create a single engine bell located at ``(x, y)`` below the core."""
    bell = conic(0.5, 0.2, 2.0, center=point(x, y, -2.0))
    nozzle = conic(0.2, 0.2, 0.5, center=point(x, y, 0.0))
    return merge_solids([bell, nozzle])


def create_fin() -> list:
    """Create one guidance fin extending along +x and +z."""
    verts = [
        point(2.0, 0.0, 0.8),
        point(2.0, 0.0, 5.5),
        point(4.5, 0.0, 2.4),
    ]
    normals = [[0, 1, 0, 0] for _ in verts]
    faces = [[0, 1, 2]]
    boundary = [0, 1, 2]
    fin_face = surface(verts, normals, faces, boundary, [])
    fin = extrude(fin_face, 0.3, direction=vect(0, 1, 0, 0))
    return fin


def merge_solids(parts: list[list]) -> list:
    """Combine a list of solids by concatenating their surfaces."""
    surfaces = []
    for part in parts:
        surfaces.extend(part[1])
    return solid(surfaces)


def build_rocket() -> tuple[list[list], list[list]]:
    """Return (components, assembly) where components is a list of solids."""

    components: list[list] = []

    # Core stack
    stage1 = stack_conic(2.0, 2.0, 12.0, 0.0)
    interstage = stack_conic(2.0, 1.6, 2.0, 12.0)
    stage2 = stack_conic(1.6, 1.6, 8.0, 14.0)
    adapter = stack_conic(1.6, 1.2, 2.0, 22.0)
    payload_barrel = stack_conic(1.2, 1.2, 4.0, 24.0)
    fairing_lower = stack_conic(1.2, 0.6, 1.5, 28.0)
    fairing_upper = stack_conic(0.6, 0.05, 1.8, 29.5)

    components.extend([stage1, interstage, stage2, adapter, payload_barrel, fairing_lower, fairing_upper])

    # Engines: one central and four around it
    engines = []
    engines.append(create_engine(0.0, 0.0))
    ring_radius = 1.3
    for dx, dy in [(ring_radius, 0.0), (-ring_radius, 0.0), (0.0, ring_radius), (0.0, -ring_radius)]:
        engines.append(create_engine(dx, dy))
    components.extend(engines)

    # Guidance fins near the base
    base_fin = create_fin()
    fins = [base_fin]
    for angle in (90, 180, 270):
        fins.append(rotatesolid(base_fin, angle, cent=point(0, 0, 0), axis=point(0, 0, 1)))
    components.extend(fins)

    assembly = merge_solids(components)
    return components, assembly


def visualise(components: list[list]) -> None:
    """Render the rocket components with simple colouring."""
    from yapcad.pyglet_drawable import pygletDraw

    palette = [
        'silver', 'gray', 'maroon', 'olive', 'green', 'aqua', 'navy', 'purple', 'fuchsia', 'lime'
    ]
    dd = pygletDraw()
    for idx, comp in enumerate(components):
        color = palette[idx % len(palette)]
        dd.make_object(f'part{idx}', lighting=True, material='pearl')
        dd.linecolor = color
        dd.draw(comp, name=f'part{idx}')
    dd.display()


def export_stl(solid_obj: list, output_path: Path) -> Path:
    output_path = output_path.resolve()
    write_stl(solid_obj, output_path, binary=True, name='yapCAD_rocket')
    return output_path


def main():
    components, assembly = build_rocket()
    out_path = export_stl(assembly, Path('rocket_demo.stl'))
    print(f'STL written to {out_path}')

    try:
        visualise(components)
    except Exception as exc:
        print(f'Visualisation skipped: {exc}')


if __name__ == '__main__':
    main()
