import math

from yapcad.geom import point
from yapcad.geom3d import poly2surfaceXY
from yapcad.geometry_utils import (
    Triangle,
    orient_triangle,
    to_point,
    to_vec3,
    triangle_area,
    triangle_centroid,
    triangle_is_degenerate,
    triangle_normal,
    triangles_from_mesh,
)
from yapcad.mesh import mesh_view


def test_to_vec3_and_to_point_roundtrip():
    pt = point(1.2, -3.4, 5.6)
    vec = to_vec3(pt)
    assert vec == (1.2, -3.4, 5.6)
    lifted = to_point(vec)
    assert lifted == (1.2, -3.4, 5.6, 1.0)


def test_triangle_normal_and_area():
    v0 = (0.0, 0.0, 0.0)
    v1 = (1.0, 0.0, 0.0)
    v2 = (0.0, 1.0, 0.0)

    n = triangle_normal(v0, v1, v2)
    assert n == (0.0, 0.0, 1.0)
    assert math.isclose(triangle_area(v0, v1, v2), 0.5)


def test_triangle_is_degenerate_when_colinear():
    v0 = (0.0, 0.0, 0.0)
    v1 = (1.0, 1.0, 1.0)
    v2 = (2.0, 2.0, 2.0)
    assert triangle_is_degenerate(v0, v1, v2)
    assert triangle_normal(v0, v1, v2) is None


def test_orient_triangle_flips_for_reference():
    v0 = (0.0, 0.0, 0.0)
    v1 = (1.0, 0.0, 0.0)
    v2 = (0.0, 1.0, 0.0)

    oriented = orient_triangle(v0, v1, v2, preferred_normal=(0.0, 0.0, -1.0))
    assert oriented[1] == v2 and oriented[2] == v1


def test_triangle_centroid():
    v0 = (0.0, 0.0, 0.0)
    v1 = (3.0, 0.0, 0.0)
    v2 = (0.0, 3.0, 0.0)
    cx, cy, cz = triangle_centroid(v0, v1, v2)
    assert math.isclose(cx, 1.0)
    assert math.isclose(cy, 1.0)
    assert math.isclose(cz, 0.0)


def test_triangles_from_mesh_produces_dataclasses():
    poly = [
        point(0, 0),
        point(1, 0),
        point(0, 1),
        point(0, 0),
    ]
    surf, _ = poly2surfaceXY(poly)
    tris = list(triangles_from_mesh(mesh_view(surf)))
    assert len(tris) == 1
    assert isinstance(tris[0], Triangle)
    assert tris[0].normal == (0.0, 0.0, 1.0)
