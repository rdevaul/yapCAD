#!/usr/bin/env python3
"""Diagnose boolean meshes produced by yapCAD using trimesh helpers."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from mesh_diagnostics import MeshDiagnostics, diagnose_mesh, export_boundary

from yapcad.geom import point
from yapcad.geom3d import solid_boolean, translatesolid
from yapcad.geom3d_util import prism, sphere, conic, tube

from yapcad.boolean.trimesh_engine import _solid_to_mesh


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


def main():
    parser = argparse.ArgumentParser(description='Diagnose boolean mesh quality.')
    parser.add_argument('--shapes', choices=['cube', 'sphere', 'box_hole', 'tube_hole'], required=True)
    parser.add_argument('--operation', choices=['union', 'intersection', 'difference'], default='union')
    parser.add_argument('--engine', default=None, help='boolean engine (e.g., native, trimesh:manifold)')
    parser.add_argument('--stitch', action='store_true', help='enable stitch flag for boolean call')
    parser.add_argument('--boundary-output', type=Path, help='optional path to export boundary edges (PLY)')
    parser.add_argument('--json', type=Path, help='write diagnostics JSON to path')
    args = parser.parse_args()

    solids = _build_solids(args.shapes)
    result = solid_boolean(solids[0], solids[1], args.operation,
                           stitch=args.stitch, engine=args.engine)

    mesh = _solid_to_mesh(result)
    diagnostics = diagnose_mesh(mesh)

    if args.boundary_output:
        export_boundary(mesh, args.boundary_output)

    print(json.dumps(diagnostics.to_dict(), indent=2))

    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(diagnostics.to_dict(), indent=2))


if __name__ == '__main__':
    main()
