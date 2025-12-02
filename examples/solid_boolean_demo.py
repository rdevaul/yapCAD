"""Demonstration of yapCAD solid boolean operations.

This script can visualise the result via OpenGL or export the mesh to STL/STEP.

STEP Export Modes:
- analytic: Uses native BREP with exact geometric definitions (PLANE, CYLINDER, etc.)
            Requires OCC engine and preserves analytic surface data.
- faceted:  Uses triangle mesh approximation. Works with any boolean engine.

The default STEP format can be set via the YAPCAD_STEP_FORMAT environment variable:
    export YAPCAD_STEP_FORMAT=analytic  # or 'faceted'

Examples:
    # Boolean with OCC engine, export analytic STEP
    python solid_boolean_demo.py --mode step --engine occ --step-format analytic

    # Boolean with native engine, export faceted STEP
    python solid_boolean_demo.py --mode step --step-format faceted

    # Use environment variable default
    export YAPCAD_STEP_FORMAT=analytic
    python solid_boolean_demo.py --mode step --engine occ
"""

import argparse
import os
from pathlib import Path

from yapcad.geom import point
from yapcad.geom3d import solid_boolean, translatesolid
from yapcad.geom3d_util import prism, sphere, conic, tube
from yapcad.brep import has_brep_data
from yapcad.io.stl import write_stl
from yapcad.io.step import write_step, write_step_analytic

# Environment variable for default STEP format
STEP_FORMAT_ENV_VAR = 'YAPCAD_STEP_FORMAT'
DEFAULT_STEP_FORMAT = 'faceted'  # Safe default that always works


def get_default_step_format():
    """Get the default STEP format from environment variable or fallback."""
    env_value = os.environ.get(STEP_FORMAT_ENV_VAR, '').lower().strip()
    if env_value in ('analytic', 'faceted'):
        return env_value
    return DEFAULT_STEP_FORMAT


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


def _export_step(result, output, step_format='faceted'):
    """Export result to STEP file.

    Parameters
    ----------
    result : solid
        The solid to export.
    output : str
        Output file path (basename or full path).
    step_format : str
        'analytic' for native BREP with exact geometry, or
        'faceted' for triangle mesh approximation.
    """
    path = Path(output)
    if path.suffix.lower() not in {'.step', '.stp'}:
        path = path.with_suffix('.step')

    if step_format == 'analytic':
        # Try analytic export, fall back to faceted if no BREP data
        success = write_step_analytic(result, str(path), fallback_to_faceted=True)
        if success:
            print(f'Wrote ANALYTIC STEP to {path}')
            _report_step_surfaces(path)
        else:
            print(f'Wrote FACETED STEP to {path} (analytic not available)')
            _report_step_surfaces(path)
    else:
        # Faceted export
        write_step(result, path)
        print(f'Wrote FACETED STEP to {path}')
        _report_step_surfaces(path)


def _report_step_surfaces(step_path):
    """Report surface types found in the STEP file."""
    try:
        with open(step_path, 'r') as f:
            content = f.read()

        surface_counts = {}
        for surface_type in ['PLANE', 'CYLINDRICAL_SURFACE', 'SPHERICAL_SURFACE',
                            'CONICAL_SURFACE', 'TOROIDAL_SURFACE', 'B_SPLINE_SURFACE']:
            count = content.count(f'{surface_type}(')
            if count > 0:
                surface_counts[surface_type] = count

        face_count = content.count('ADVANCED_FACE')
        if face_count > 0:
            surface_counts['FACES'] = face_count

        if surface_counts:
            print(f'  Surface types: {surface_counts}')
    except Exception:
        pass  # Don't fail if we can't read the file


def main():
    default_format = get_default_step_format()

    parser = argparse.ArgumentParser(
        description='Solid boolean demonstration.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f'''
STEP Format Options:
  analytic  - Native BREP with exact geometric definitions (PLANE, CYLINDER, etc.)
              Requires --engine occ to preserve analytic surface data.
  faceted   - Triangle mesh approximation. Works with any engine.

Environment Variable:
  Set {STEP_FORMAT_ENV_VAR} to 'analytic' or 'faceted' to change the default.
  Current default: {default_format}
'''
    )
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
                        help='boolean engine to use (occ for analytic, native otherwise)')
    parser.add_argument('--step-format', choices=['analytic', 'faceted'],
                        default=default_format,
                        help=f'STEP export format (default: {default_format}, '
                             f'set {STEP_FORMAT_ENV_VAR} env var to change)')
    args = parser.parse_args()

    # Warn if analytic format requested without OCC engine
    if args.mode == 'step' and args.step_format == 'analytic' and args.engine != 'occ':
        print(f"Note: --step-format analytic works best with --engine occ")
        print(f"      Without OCC engine, result may fall back to faceted export.")
        print()

    solids = _build_solids(args.shapes)
    engine = args.engine
    if engine == 'occ':
        if not (has_brep_data(solids[0]) and has_brep_data(solids[1])):
            print("OCC engine requested but one or both solids lack BREP metadata. "
                  "Falling back to native boolean.")
            engine = None

    result = solid_boolean(solids[0], solids[1], args.operation,
                           stitch=args.stitch, engine=engine)

    if args.mode == 'gl':
        _render_opengl(solids, result)
    elif args.mode == 'stl':
        _export_stl(result, args.output)
    else:
        _export_step(result, args.output, step_format=args.step_format)


if __name__ == '__main__':
    main()
