"""Triangulation helpers for yapCAD surfaces.

We delegate to ``mapbox-earcut`` (the fast ear clipping implementation
used by Mapbox GL) to keep the logic compact and reliable.  The helper
routine in this file simply normalises yapCAD polygon inputs into the
format expected by earcut and converts the resulting indices back into
triangle vertex tuples.
"""

from __future__ import annotations

from typing import Iterable, List, Sequence, Tuple

import numpy as np

try:
    import mapbox_earcut as _earcut
except ImportError as exc:  # pragma: no cover - import guard
    raise ImportError(
        "mapbox-earcut must be installed to triangulate polygons with holes"
    ) from exc

from yapcad.geom import epsilon

Point2D = Tuple[float, float]


def triangulate_polygon(outer: Sequence[Sequence[float]],
                        holes: Iterable[Sequence[Sequence[float]]] | None = None
                        ) -> List[List[Point2D]]:
    """Return triangles covering ``outer`` minus any ``holes``.

    ``outer`` and each entry in ``holes`` is expected to be a sequence of
    XY-like points.  Degenerate loops (fewer than three distinct points)
    are ignored.  The returned triangles are lists of three ``(x, y)``
    pairs.  Downstream code is responsible for lifting the coordinates
    into 3D and enforcing winding to match the target surface normal.
    """

    if holes is None:
        holes = []

    outer_loop = _prepare_loop(outer, want_ccw=True)
    if len(outer_loop) < 3:
        return []

    point_map: List[Point2D] = []

    ring_ends: List[int] = []

    def _append(loop: Sequence[Point2D]) -> None:
        for x, y in loop:
            point_map.append((x, y))
        ring_ends.append(len(point_map))

    _append(outer_loop)

    for hole in holes:
        loop = _prepare_loop(hole, want_ccw=False)
        if len(loop) < 3:
            continue
        _append(loop)

    vertices = np.asarray(point_map, dtype=np.float32)
    ring_array = np.asarray(ring_ends, dtype=np.uint32)
    indices = _earcut.triangulate_float32(vertices, ring_array)
    triangles: List[List[Point2D]] = []
    for i in range(0, len(indices), 3):
        triangles.append([point_map[indices[i]],
                          point_map[indices[i + 1]],
                          point_map[indices[i + 2]]])
    return triangles


def _prepare_loop(points: Sequence[Sequence[float]], *, want_ccw: bool) -> List[Point2D]:
    loop: List[Point2D] = []
    for pt in points:
        x, y = float(pt[0]), float(pt[1])
        if loop and _near(loop[-1], (x, y)):
            continue
        loop.append((x, y))
    if loop and _near(loop[0], loop[-1]):
        loop.pop()
    if len(loop) < 3:
        return loop
    area = _signed_area(loop)
    if want_ccw and area < 0:
        loop.reverse()
    elif not want_ccw and area > 0:
        loop.reverse()
    return loop


def _near(p1: Point2D, p2: Point2D) -> bool:
    return abs(p1[0] - p2[0]) <= epsilon and abs(p1[1] - p2[1]) <= epsilon


def _signed_area(loop: Sequence[Point2D]) -> float:
    total = 0.0
    for i, (x0, y0) in enumerate(loop):
        x1, y1 = loop[(i + 1) % len(loop)]
        total += x0 * y1 - x1 * y0
    return total / 2.0
