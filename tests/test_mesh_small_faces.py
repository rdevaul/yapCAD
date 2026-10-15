"""Mesh transport must not punch holes by dropping small valid triangles."""

from collections import Counter
import io
import math
import struct

import pytest
from yapcad.geom import point
from yapcad.geom3d import solid, surface
from yapcad.geometry_utils import triangle_normal
from yapcad.io.stl import write_stl
from yapcad.mesh import mesh_view


@pytest.mark.parametrize("size", [1.0, 1e-3, 1e-4])
def test_small_closed_mesh_survives_view_and_stl(size):
    vertices = [
        point(0, 0, 0),
        point(size, 0, 0),
        point(0, size, 0),
        point(0, 0, size),
    ]
    faces = [[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]]
    obj = solid(
        [surface(vertices, [[0, 0, 1, 0] for _ in vertices], faces)], [], []
    )
    triangles = list(mesh_view(obj))
    assert len(triangles) == 4
    edges = Counter()
    for normal, *corners in triangles:
        assert math.hypot(*normal) == pytest.approx(1.0)
        for i in range(3):
            edges[tuple(sorted((corners[i], corners[(i + 1) % 3])))] += 1
    assert set(edges.values()) == {2}
    output = io.BytesIO()
    write_stl(obj, output)
    assert struct.unpack_from("<I", output.getvalue(), 80)[0] == 4
    assert len(output.getvalue()) == 84 + 4 * 50


def test_exactly_degenerate_face_is_still_skipped():
    vertices = [point(0, 0, 0), point(1e-4, 0, 0), point(2e-4, 0, 0)]
    obj = solid(
        [surface(vertices, [[0, 0, 1, 0] for _ in vertices], [[0, 1, 2]])],
        [],
        [],
    )
    assert list(mesh_view(obj)) == []


def test_normal_default_tolerance_is_unchanged():
    points = [(0, 0, 0), (1e-4, 0, 0), (0, 1e-4, 0)]
    assert triangle_normal(*points) is None
    assert triangle_normal(*points, tol=0.0) == (0.0, 0.0, 1.0)
