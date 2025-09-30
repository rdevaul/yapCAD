"""Validation helpers for yapCAD geometry."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Iterable, List, Sequence

from yapcad.geom import dist, epsilon, point
from yapcad.geometry_utils import triangle_normal


def is_closed_polygon(points: Sequence[Sequence[float]], tol: float = epsilon) -> bool:
    """Return ``True`` if a polyline is closed within ``tol``."""

    if not points:
        return False
    first = points[0]
    last = points[-1]
    if len(first) < 3 or len(last) < 3:
        return False
    return dist(first, last) <= tol


def faces_oriented(surface: Sequence) -> "CheckResult":
    from yapcad.geom3d import issurface

    if not issurface(surface):
        raise ValueError('faces_oriented expects a surface')

    verts = surface[1]
    faces = surface[3]

    reference = None
    inconsistent = []

    for idx, face in enumerate(faces):
        if len(face) != 3:
            continue
        v0 = verts[face[0]]
        v1 = verts[face[1]]
        v2 = verts[face[2]]
        normal = triangle_normal((v0[0], v0[1], v0[2]),
                                 (v1[0], v1[1], v1[2]),
                                 (v2[0], v2[1], v2[2]))
        if normal is None:
            continue
        if reference is None:
            reference = normal
            continue
        dot = reference[0] * normal[0] + reference[1] * normal[1] + reference[2] * normal[2]
        if dot < -epsilon:
            inconsistent.append(idx)

    if reference is None:
        return CheckResult(True, ['no non-degenerate faces found'])
    if inconsistent:
        return CheckResult(False, [f'inconsistent face orientation indices: {inconsistent}'])
    return CheckResult(True, [])


def surface_watertight(surface: Sequence) -> "CheckResult":
    from yapcad.geom3d import issurface

    if not issurface(surface):
        raise ValueError('surface_watertight expects a surface')

    faces = surface[3]
    edges = Counter()

    for face in faces:
        if len(face) != 3:
            continue
        a, b, c = face
        edges[_edge_key(a, b)] += 1
        edges[_edge_key(b, c)] += 1
        edges[_edge_key(c, a)] += 1

    boundary = [edge for edge, count in edges.items() if count == 1]
    invalid = [edge for edge, count in edges.items() if count > 2]

    warnings: List[str] = []
    ok = True
    if boundary:
        ok = False
        warnings.append(f'{len(boundary)} boundary edges detected')
    if invalid:
        ok = False
        warnings.append(f'edges with multiplicity >2: {invalid}')

    return CheckResult(ok, warnings)


def _edge_key(a: int, b: int) -> tuple[int, int]:
    return (a, b) if a < b else (b, a)


@dataclass
class CheckResult:
    ok: bool
    warnings: List[str]

    def __bool__(self) -> bool:
        return self.ok


__all__ = [
    'CheckResult',
    'is_closed_polygon',
    'faces_oriented',
    'surface_watertight',
]
