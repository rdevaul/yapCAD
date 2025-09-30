import math

from yapcad.geom import point
from yapcad.geom3d import surface
from yapcad.geometry_checks import (
    CheckResult,
    faces_oriented,
    is_closed_polygon,
    surface_watertight,
)
from yapcad.geometry_utils import triangle_normal


def test_is_closed_polygon_true():
    pts = [
        point(0, 0),
        point(1, 0),
        point(1, 1),
        point(0, 1),
        point(0, 0),
    ]
    assert is_closed_polygon(pts)


def test_is_closed_polygon_false():
    pts = [
        point(0, 0),
        point(1, 0),
        point(1, 1),
    ]
    assert not is_closed_polygon(pts)


def _make_square_surface():
    verts = [
        point(0, 0, 0),
        point(1, 0, 0),
        point(1, 1, 0),
        point(0, 1, 0),
    ]
    normals = [[0, 0, 1, 0] for _ in verts]
    faces = [
        [0, 1, 2],
        [0, 2, 3],
    ]
    return surface(verts, normals, faces, [], [])


def _make_tetra_surface():
    verts = [
        point(0, 0, 0),
        point(1, 0, 0),
        point(0, 1, 0),
        point(0, 0, 1),
    ]
    normals = [[0, 0, 1, 0] for _ in verts]
    faces = [
        [0, 1, 2],
        [0, 1, 3],
        [1, 2, 3],
        [2, 0, 3],
    ]
    return surface(verts, normals, faces, [], [])


def test_faces_oriented_ok():
    surf = _make_square_surface()
    result = faces_oriented(surf)
    assert isinstance(result, CheckResult)
    assert result.ok


def test_faces_oriented_detects_flip():
    surf = _make_square_surface()
    surf[3][1] = [surf[3][1][0], surf[3][1][2], surf[3][1][1]]
    result = faces_oriented(surf)
    assert not result.ok
    assert result.warnings


def test_surface_watertight_detects_boundary():
    surf = _make_square_surface()
    result = surface_watertight(surf)
    assert not result.ok
    assert any('boundary' in msg for msg in result.warnings)


def test_surface_watertight_closed_mesh():
    surf = _make_tetra_surface()
    result = surface_watertight(surf)
    assert result.ok
