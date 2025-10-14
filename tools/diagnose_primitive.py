#!/usr/bin/env python3
"""Diagnose individual yapCAD primitives using trimesh helpers."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict

from mesh_diagnostics import diagnose_mesh, export_boundary

from yapcad.geom import point
from yapcad.geom3d_util import prism, sphere, conic, tube
from yapcad.boolean.trimesh_engine import _solid_to_mesh


def _build_primitive(name: str, params: Dict[str, Any]):
    if name == 'sphere':
        diameter = float(params.get('diameter', 2.0))
        depth = int(params.get('depth', 2))
        return sphere(diameter, depth=depth)
    if name == 'prism':
        length = float(params.get('length', 2.0))
        width = float(params.get('width', 2.0))
        height = float(params.get('height', 2.0))
        return prism(length, width, height)
    if name == 'conic':
        baser = float(params.get('base_radius', 1.0))
        topr = float(params.get('top_radius', 0.5))
        height = float(params.get('height', 2.0))
        angr = float(params.get('angular_resolution', 10.0))
        return conic(baser, topr, height, center=params.get('center', point(0, 0, 0)), angr=angr)
    if name == 'tube':
        outer_diameter = float(params.get('outer_diameter', 3.0))
        wall_thickness = float(params.get('wall_thickness', 0.5))
        length = float(params.get('length', 4.0))
        include_caps = bool(params.get('include_caps', True))
        base_point = params.get('base_point')
        if base_point is not None:
            base_point = point(*base_point)
        return tube(outer_diameter=outer_diameter, wall_thickness=wall_thickness,
                    length=length, base_point=base_point, include_caps=include_caps)
    raise ValueError(f'unsupported primitive {name!r}')


def main():
    parser = argparse.ArgumentParser(description='Diagnose individual yapCAD primitives.')
    parser.add_argument('--primitive', choices=['sphere', 'prism', 'conic', 'tube'], required=True)
    parser.add_argument('--params', help='JSON object with primitive parameters (e.g. {"diameter": 3.0})')
    parser.add_argument('--boundary-output', type=Path, help='optional path to export boundary edges (PLY)')
    parser.add_argument('--json', type=Path, help='write diagnostics JSON to path')
    parser.add_argument('--export-stl', type=Path, help='optional path to export the mesh as STL')
    args = parser.parse_args()

    params: Dict[str, Any] = {}
    if args.params:
        try:
            params = json.loads(args.params)
        except json.JSONDecodeError as exc:
            raise SystemExit(f'Failed to parse --params JSON: {exc}') from exc

    solid = _build_primitive(args.primitive, params)
    mesh = _solid_to_mesh(solid)
    diagnostics = diagnose_mesh(mesh)

    if args.boundary_output:
        export_boundary(mesh, args.boundary_output)

    if args.export_stl:
        args.export_stl.parent.mkdir(parents=True, exist_ok=True)
        mesh.export(str(args.export_stl))

    print(json.dumps(diagnostics.to_dict(), indent=2))

    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(diagnostics.to_dict(), indent=2))


if __name__ == '__main__':
    main()
