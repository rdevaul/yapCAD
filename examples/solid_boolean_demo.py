"""Demonstration of yapCAD solid boolean operations.

This script can visualise the result via OpenGL or export the mesh to STL/STEP.
"""

import argparse
from pathlib import Path

from yapcad.geom import point
from yapcad.geom3d import solid_boolean, translatesolid
from yapcad.geom3d_util import prism, sphere, conic, tube
from yapcad.io.stl import write_stl
from yapcad.io.step import write_step


def _build_solids(kind: str):
    if kind == 'cube':
        a = prism(2, 2, 2)
        b = translatesolid(prism(2, 2, 2), point(0.75, 0.0, 0.0))
    elif kind == 'sphere':
        a = sphere(2.0)
        b = translatesolid(sphere(2.0), point(1.0, 0.0, 0.0))
    elif kind == 'box_hole':
        a = prism(2, 2, 2)
        b = conic(0.6, 0.6, 3.0, center=point(0.0, 0.0, -1.5))
    elif kind == 'tube_hole':
        a = tube(outer_diameter=3.0, wall_thickness=0.5, length=4.0,
                 base_point=point(0.0, 0.0, -2.0))
        b = conic(0.6, 0.6, 4.2, center=point(1.2, 0.0, -2.1))
    else:
        raise ValueError(f'unsupported shape kind: {kind!r}')
    return a, b


def _render_opengl(solids, result):
    from yapcad.pyglet_drawable import pygletDraw

    viewer = pygletDraw()
    viewer.make_object('obj1', lighting=True, linecolor='white',
                       material='pearl',
                       position=point(-2.5, 2.5, 0))

    viewer.make_object('obj2', lighting=True, linecolor='white',
                       material='obsidian',
                       position=point(2.5, 2.5, 0))

    viewer.make_object('rslt', lighting=True, linecolor='white',
                       material='gold',
                       position=point(0, 0, 0))

    viewer.draw_solid(solids[0], name='obj1')
    viewer.draw_solid(solids[1], name='obj2')
    viewer.draw_solid(result, name='rslt')
    viewer.display()


def _export_stl(result, output):
    path = Path(output)
    if path.suffix.lower() != '.stl':
        path = path.with_suffix('.stl')
    write_stl(result, path)
    print(f'Wrote STL to {path}')


def _export_step(result, output):
    path = Path(output)
    if path.suffix.lower() not in {'.step', '.stp'}:
        path = path.with_suffix('.step')
    write_step(result, path)
    print(f'Wrote STEP to {path}')


def main():
    parser = argparse.ArgumentParser(description='Solid boolean demonstration.')
    parser.add_argument('--mode', choices=['gl', 'stl', 'step'], default='gl',
                        help='visualise with OpenGL or export as STL/STEP')
    parser.add_argument('--operation', choices=['union', 'intersection', 'difference'],
                        default='union', help='boolean operation to apply')
    parser.add_argument('--shapes', choices=['cube', 'sphere', 'box_hole', 'tube_hole'],
                        default='cube',
                        help='primitive pair to use')
    parser.add_argument('--output', default='boolean_result',
                        help='output basename for STL/STEP export')
    parser.add_argument('--stitch', action='store_true',
                        help='experimentally stitch open edges in the result')
    parser.add_argument('--engine', default=None,
                        help='boolean engine to use (default: native)')
    args = parser.parse_args()

    solids = _build_solids(args.shapes)
    result = solid_boolean(solids[0], solids[1], args.operation,
                           stitch=args.stitch, engine=args.engine)

    if args.mode == 'gl':
        _render_opengl(solids, result)
    elif args.mode == 'stl':
        _export_stl(result, args.output)
    else:
        _export_step(result, args.output)


if __name__ == '__main__':
    main()
