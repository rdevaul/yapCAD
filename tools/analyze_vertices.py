#!/usr/bin/env python3
"""Analyze yapCAD solids for orphaned and duplicated vertices.

This tool operates directly on yapCAD's internal representation,
analyzing vertex usage and identifying potential geometry issues with
configurable epsilon tolerance for determining vertex proximity.

Unlike trimesh-based diagnostics, this operates at the yapCAD level
and can detect:
- Truly orphaned vertices (not referenced by any face)
- Near-duplicate vertices (within epsilon distance)
- Vertex usage statistics per surface

Example:
    PYTHONPATH=./src python tools/analyze_vertices.py \\
        --primitive tube \\
        --params '{"outer_diameter": 3.0, "wall_thickness": 0.5, "length": 4.0}' \\
        --epsilon 1e-9
"""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, List, Set, Tuple

from yapcad.geom import point
from yapcad.geom3d import issolid, issurface
from yapcad.geom3d_util import prism, sphere, conic, tube


@dataclass
class VertexAnalysis:
    """Results of vertex analysis for a yapCAD solid."""
    total_vertices: int
    unique_vertices: int  # After de-duplication with epsilon
    orphaned_vertices: int  # Not referenced by any face
    duplicate_groups: int  # Number of groups of duplicate vertices
    max_duplicates_in_group: int  # Largest duplicate group
    surfaces_analyzed: int
    epsilon: float

    # Per-surface statistics
    surface_stats: List[Dict[str, Any]]

    def to_dict(self) -> dict:
        return asdict(self)


def _vertex_key(v: list, epsilon: float) -> tuple:
    """Create a hashable key for a vertex with epsilon precision."""
    if epsilon <= 0:
        epsilon = 1e-10
    scale = 1.0 / epsilon
    return (
        round(v[0] * scale),
        round(v[1] * scale),
        round(v[2] * scale),
    )


def _vertex_distance(v1: list, v2: list) -> float:
    """Calculate Euclidean distance between two vertices."""
    dx = v1[0] - v2[0]
    dy = v1[1] - v2[1]
    dz = v1[2] - v2[2]
    return math.sqrt(dx*dx + dy*dy + dz*dz)


def analyze_solid_vertices(solid, epsilon: float = 1e-9) -> VertexAnalysis:
    """Analyze vertices in a yapCAD solid.

    Args:
        solid: A yapCAD solid object
        epsilon: Distance threshold for considering vertices as duplicates

    Returns:
        VertexAnalysis object with detailed vertex statistics
    """
    if not issolid(solid):
        raise ValueError('Input must be a yapCAD solid')

    surfaces = solid[1]  # Extract surfaces from solid

    # Global vertex tracking
    all_vertices: List[list] = []
    vertex_to_indices: Dict[tuple, List[int]] = defaultdict(list)
    referenced_indices: Set[int] = set()

    # Per-surface statistics
    surface_stats = []

    for surf_idx, surf in enumerate(surfaces):
        if not issurface(surf):
            continue

        surf_vertices = surf[1]
        surf_faces = surf[3]

        # Track which local vertices are referenced by faces
        local_referenced = set()
        for face in surf_faces:
            for vertex_idx in face:
                local_referenced.add(vertex_idx)

        local_orphaned = len(surf_vertices) - len(local_referenced)

        # Add vertices to global list
        start_idx = len(all_vertices)
        for local_idx, v in enumerate(surf_vertices):
            global_idx = start_idx + local_idx
            all_vertices.append(v)

            # Track by epsilon-rounded key
            key = _vertex_key(v, epsilon)
            vertex_to_indices[key].append(global_idx)

            # Track if this vertex is referenced
            if local_idx in local_referenced:
                referenced_indices.add(global_idx)

        surface_stats.append({
            'surface_index': surf_idx,
            'vertices': len(surf_vertices),
            'faces': len(surf_faces),
            'orphaned_vertices': local_orphaned,
            'referenced_vertices': len(local_referenced),
        })

    # Analyze duplicates
    duplicate_groups = 0
    max_duplicates = 0

    for key, indices in vertex_to_indices.items():
        if len(indices) > 1:
            # Verify they're actually close (not just hash collisions)
            verified_group = [indices[0]]
            for idx in indices[1:]:
                if any(_vertex_distance(all_vertices[idx], all_vertices[ref_idx]) < epsilon
                       for ref_idx in verified_group):
                    verified_group.append(idx)

            if len(verified_group) > 1:
                duplicate_groups += 1
                max_duplicates = max(max_duplicates, len(verified_group))

    total_vertices = len(all_vertices)
    orphaned_vertices = total_vertices - len(referenced_indices)
    unique_vertices = total_vertices - (sum(len(indices) - 1 for indices in vertex_to_indices.values()
                                             if len(indices) > 1))

    return VertexAnalysis(
        total_vertices=total_vertices,
        unique_vertices=unique_vertices,
        orphaned_vertices=orphaned_vertices,
        duplicate_groups=duplicate_groups,
        max_duplicates_in_group=max_duplicates,
        surfaces_analyzed=len(surfaces),
        epsilon=epsilon,
        surface_stats=surface_stats,
    )


def _build_primitive(name: str, params: Dict[str, Any]):
    """Build a yapCAD primitive from name and parameters."""
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
        return conic(baser, topr, height, center=point(0, 0, 0), angr=angr)
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
    parser = argparse.ArgumentParser(
        description='Analyze yapCAD solid vertices for orphans and duplicates.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('--primitive', choices=['sphere', 'prism', 'conic', 'tube'], required=True)
    parser.add_argument('--params', help='JSON object with primitive parameters')
    parser.add_argument('--epsilon', type=float, default=1e-9,
                       help='Distance threshold for duplicate detection (default: 1e-9)')
    parser.add_argument('--json-output', type=Path, help='Write JSON results to file')

    args = parser.parse_args()

    params: Dict[str, Any] = {}
    if args.params:
        try:
            params = json.loads(args.params)
        except json.JSONDecodeError as exc:
            raise SystemExit(f'Failed to parse --params JSON: {exc}') from exc

    solid = _build_primitive(args.primitive, params)
    analysis = analyze_solid_vertices(solid, epsilon=args.epsilon)

    if args.json_output:
        args.json_output.parent.mkdir(parents=True, exist_ok=True)
        args.json_output.write_text(json.dumps(analysis.to_dict(), indent=2))
        print(f'Results written to {args.json_output}')
    else:
        print(json.dumps(analysis.to_dict(), indent=2))


if __name__ == '__main__':
    main()
