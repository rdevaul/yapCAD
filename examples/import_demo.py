"""Geometry importer + boolean demo for yapCAD.

This example:
1. Loads a STEP or STL file via the appropriate importer.
2. Scales the geometry so its largest bounding-box dimension becomes 10 units.
3. Recenters the geometry at the origin for easier viewing.
4. Shows the tessellated result in the pyglet viewer.
5. (Optional) Applies a boolean op with a radius-5 sphere centered on the bbox.

For STEP files, requires the ``yapcad-brep`` conda environment so that
pythonocc-core is available; activate via ``conda activate yapcad-brep``.
STL files work without the conda environment.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, List, Tuple

from yapcad.brep import BrepSolid, attach_brep_to_solid, has_brep_data, occ_available
from yapcad.geom import point
from yapcad.geom3d import solid, solid_boolean, solidbbox, issolid
from yapcad.geom3d_util import sphere
from yapcad.io.step_importer import import_step
from yapcad.io.stl import import_stl
from yapcad.boolean.trimesh_engine import engines_available


def _detect_format(path: Path) -> str:
    """Detect file format from extension."""
    suffix = path.suffix.lower()
    if suffix in ('.step', '.stp'):
        return 'step'
    elif suffix == '.stl':
        return 'stl'
    else:
        raise ValueError(f"Unsupported file format: {suffix}")


def _combined_bbox_from_geoms(geoms) -> Tuple[point, point]:
    """Compute combined bounding box from Geometry objects."""
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


def _combined_bbox_from_solid(sld) -> Tuple[point, point]:
    """Compute bounding box from a yapCAD solid."""
    bb = solidbbox(sld)
    return point(bb[0][0], bb[0][1], bb[0][2]), point(bb[1][0], bb[1][1], bb[1][2])


def _scale_and_center_geoms(geoms, target_size: float) -> Tuple[float, point]:
    """Scale and center Geometry objects (STEP import)."""
    bbox_min, bbox_max = _combined_bbox_from_geoms(geoms)
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
    scaled_min, scaled_max = _combined_bbox_from_geoms(geoms)
    scaled_center = point((scaled_min[0] + scaled_max[0]) / 2.0,
                          (scaled_min[1] + scaled_max[1]) / 2.0,
                          (scaled_min[2] + scaled_max[2]) / 2.0)
    delta = point(-scaled_center[0], -scaled_center[1], -scaled_center[2])
    for geo in geoms:
        geo.translate(delta)
    return factor, point(0.0, 0.0, 0.0)


def _scale_and_center_solid(sld, target_size: float) -> Tuple[object, float, point]:
    """Scale and center a yapCAD solid (STL import)."""
    from yapcad.geom3d import scale as scale_solid, translate as translate_solid

    bbox_min, bbox_max = _combined_bbox_from_solid(sld)
    dims = [bbox_max[i] - bbox_min[i] for i in range(3)]
    largest = max(dims)
    if largest <= 0:
        raise ValueError("unable to determine bounding box size for scaling")
    factor = target_size / largest
    bbox_center = point((bbox_min[0] + bbox_max[0]) / 2.0,
                        (bbox_min[1] + bbox_max[1]) / 2.0,
                        (bbox_min[2] + bbox_max[2]) / 2.0)

    # Scale around bbox center
    sld = scale_solid(sld, factor, cent=bbox_center)

    # Recenter around origin
    scaled_min, scaled_max = _combined_bbox_from_solid(sld)
    scaled_center = point((scaled_min[0] + scaled_max[0]) / 2.0,
                          (scaled_min[1] + scaled_max[1]) / 2.0,
                          (scaled_min[2] + scaled_max[2]) / 2.0)
    delta = point(-scaled_center[0], -scaled_center[1], -scaled_center[2])
    sld = translate_solid(sld, delta)

    return sld, factor, point(0.0, 0.0, 0.0)


def _geometry_to_solid(geo):
    """Convert a Geometry object to a yapCAD solid."""
    surface = geo.surface()
    sld = solid([surface])
    elem = geo.geom
    if isinstance(elem, BrepSolid):
        attach_brep_to_solid(sld, elem)
    return sld


def _merge_solids(solids: Iterable) -> object:
    """Merge multiple solids into one via union."""
    solids = list(solids)
    if not solids:
        raise ValueError("no solids to merge")
    result = solids[0]
    for other in solids[1:]:
        result = solid_boolean(result, other, 'union')
    return result


def _render_scene(base_solid, boolean_result=None):
    """Render the solid(s) in the pyglet viewer."""
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
    parser = argparse.ArgumentParser(
        description="Geometry import (STEP/STL) + boolean viewer demo.")
    parser.add_argument('input_file', type=Path,
                        help='Path to the STEP (.step, .stp) or STL (.stl) file to import.')
    parser.add_argument('--target-size', type=float, default=10.0,
                        help='Scale imported geometry so the max bbox dimension equals this value.')
    parser.add_argument('--boolean', choices=['none', 'union', 'intersection', 'difference'],
                        default='none', help='Optional boolean op with a sphere.')
    parser.add_argument('--engine', default='auto',
                        help='Boolean engine (auto picks best available, otherwise native).')
    parser.add_argument('--sphere-radius', type=float, default=5.0,
                        help='Radius of the optional boolean sphere.')
    args = parser.parse_args()

    # Detect file format
    file_format = _detect_format(args.input_file)
    print(f"Detected format: {file_format.upper()}")

    # Import based on format
    if file_format == 'step':
        geoms = import_step(args.input_file)
        if not geoms:
            raise SystemExit("No solids found in STEP file.")
        scale_factor, center = _scale_and_center_geoms(geoms, args.target_size)
        print(f"Imported {len(geoms)} solid(s). Applied scale factor: {scale_factor:.3f}")
        solids = [_geometry_to_solid(geo) for geo in geoms]
        base_solid = _merge_solids(solids)
    else:  # stl
        sld = import_stl(str(args.input_file))
        if not issolid(sld) or not sld[1]:
            raise SystemExit("No geometry found in STL file.")
        base_solid, scale_factor, center = _scale_and_center_solid(sld, args.target_size)
        # Count triangles from the surface
        surf = base_solid[1][0]
        n_faces = len(surf[3]) if len(surf) > 3 else 0
        n_verts = len(surf[1]) if len(surf) > 1 else 0
        print(f"Imported STL with {n_faces} triangles, {n_verts} vertices. "
              f"Applied scale factor: {scale_factor:.3f}")

    boolean_result = None
    backend = args.engine
    if backend == 'auto':
        if has_brep_data(base_solid) and occ_available():
            backend = 'occ'
        else:
            available = engines_available()
            backend = next(iter(available)) if available else 'native'

    if args.boolean != 'none':
        sphere_solid = sphere(args.sphere_radius * 2.0, center=center)
        print(f"Applying {args.boolean} with sphere (radius={args.sphere_radius}) "
              f"using engine: {backend}")
        try:
            boolean_result = solid_boolean(base_solid, sphere_solid, args.boolean,
                                           engine=backend)
        except (ValueError, RuntimeError) as exc:
            raise SystemExit(
                f"Boolean failed with engine '{backend}'. Try setting "
                "--engine to a working backend (see trimesh.boolean.engines_available) "
                "or use --engine native."
            ) from exc

    _render_scene(base_solid, boolean_result)


if __name__ == '__main__':
    main()
