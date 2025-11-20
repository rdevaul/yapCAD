"""STEP importer + boolean demo for yapCAD.

This example:
1. Loads a STEP file via the pythonocc-backed importer.
2. Scales the geometry so its largest bounding-box dimension becomes 10 units.
3. Recenters the geometry at the origin for easier viewing.
4. Shows the tessellated result in the pyglet viewer.
5. (Optional) Applies a boolean op with a radius-5 sphere centered on the bbox.

Requires the ``yapcad-brep`` conda environment so that pythonocc-core is
available; activate via ``conda activate yapcad-brep`` before running.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, List, Tuple

try:
    from yapcad.brep import occ_available
except ImportError:  # pragma: no cover - older installs
    def occ_available() -> bool:
        return False
from yapcad.geom import point
from yapcad.geom3d import solid, solid_boolean, solidbbox
from yapcad.geom3d_util import sphere
from yapcad.io.step_importer import import_step
from yapcad.boolean.trimesh_engine import engines_available


def _combined_bbox(geoms) -> Tuple[point, point]:
    bbox = None
    for geo in geoms:
        bb = geo.bbox
        if bbox is None:
            bbox = [point(bb[0][0], bb[0][1], bb[0][2]),
                    point(bb[1][0], bb[1][1], bb[1][2])]
        else:
            bbox[0] = point(min(bbox[0][0], bb[0][0]),
                            min(bbox[0][1], bb[0][1]),
                            min(bbox[0][2], bb[0][2]))
            bbox[1] = point(max(bbox[1][0], bb[1][0]),
                            max(bbox[1][1], bb[1][1]),
                            max(bbox[1][2], bb[1][2]))
    if bbox is None:
        raise ValueError("no geometry to compute bounding box")
    return bbox[0], bbox[1]


def _scale_and_center_geoms(geoms, target_size: float) -> Tuple[float, point]:
    bbox_min, bbox_max = _combined_bbox(geoms)
    dims = [bbox_max[i] - bbox_min[i] for i in range(3)]
    largest = max(dims)
    if largest <= 0:
        raise ValueError("unable to determine bounding box size for scaling")
    factor = target_size / largest
    bbox_center = point((bbox_min[0] + bbox_max[0]) / 2.0,
                        (bbox_min[1] + bbox_max[1]) / 2.0,
                        (bbox_min[2] + bbox_max[2]) / 2.0)
    for geo in geoms:
        geo.scale(factor, cent=bbox_center)
    # recenter around origin for easier viewing
    scaled_min, scaled_max = _combined_bbox(geoms)
    scaled_center = point((scaled_min[0] + scaled_max[0]) / 2.0,
                          (scaled_min[1] + scaled_max[1]) / 2.0,
                          (scaled_min[2] + scaled_max[2]) / 2.0)
    delta = point(-scaled_center[0], -scaled_center[1], -scaled_center[2])
    for geo in geoms:
        geo.translate(delta)
    return factor, point(0.0, 0.0, 0.0)


def _geometry_to_solid(geo):
    surface = geo.surface()
    return solid([surface])


def _merge_solids(solids: Iterable) -> object:
    solids = list(solids)
    if not solids:
        raise ValueError("no solids to merge")
    result = solids[0]
    for other in solids[1:]:
        result = solid_boolean(result, other, 'union')
    return result


def _render_scene(base_solid, boolean_result=None):
    from yapcad.pyglet_drawable import pygletDraw

    viewer = pygletDraw()
    bbox = solidbbox(base_solid)
    base_width = bbox[1][0] - bbox[0][0]
    spacing = max(base_width * 0.1, 5.0)

    viewer.make_object('scaled', lighting=True, material='pearl',
                       position=point(-(base_width / 2 + spacing / 2), 0.0, 0.0))
    viewer.draw_solid(base_solid, name='scaled')

    if boolean_result is not None:
        result_bbox = solidbbox(boolean_result)
        result_width = result_bbox[1][0] - result_bbox[0][0]
        offset = base_width / 2 + result_width / 2 + spacing
        viewer.make_object('boolean', lighting=True, material='gold',
                           position=point(offset, 0.0, 0.0))
        viewer.draw_solid(boolean_result, name='boolean')

    viewer.display()


def main():
    parser = argparse.ArgumentParser(description="STEP import + boolean viewer demo.")
    parser.add_argument('step_file', type=Path,
                        help='Path to the STEP file to import.')
    parser.add_argument('--target-size', type=float, default=10.0,
                        help='Scale imported geometry so the max bbox dimension equals this value.')
    parser.add_argument('--boolean', choices=['none', 'union', 'intersection', 'difference'],
                        default='none', help='Optional boolean op with a radius-5 sphere.')
    parser.add_argument('--engine', default='auto',
                        help='Boolean engine (auto picks first trimesh backend, otherwise native).')
    parser.add_argument('--sphere-radius', type=float, default=5.0,
                        help='Radius of the optional boolean sphere.')
    args = parser.parse_args()

    if not occ_available():  # pragma: no cover - requires OCC runtime
        raise SystemExit("pythonocc-core unavailable. Activate the yapcad-brep conda environment.")

    geoms = import_step(args.step_file)
    if not geoms:
        raise SystemExit("No solids found in STEP file.")

    scale_factor, center = _scale_and_center_geoms(geoms, args.target_size)
    print(f"Imported {len(geoms)} solids. Applied scale factor: {scale_factor:.3f}")

    solids = [_geometry_to_solid(geo) for geo in geoms]
    base_solid = _merge_solids(solids)

    boolean_result = None
    backend = args.engine
    if backend == 'auto':
        available = engines_available()
        backend = next(iter(available)) if available else None

    if args.boolean != 'none':
        sphere_solid = sphere(args.sphere_radius * 2.0, center=center)
        try:
            boolean_result = solid_boolean(base_solid, sphere_solid, args.boolean,
                                           engine=backend)
        except (ValueError, RuntimeError) as exc:
            raise SystemExit(
                "Boolean failed with the selected engine. Try setting "
                "--engine to a working trimesh backend (see trimesh.boolean.engines_available) "
                "or use --engine native."
            ) from exc

    _render_scene(base_solid, boolean_result)


if __name__ == '__main__':
    main()
