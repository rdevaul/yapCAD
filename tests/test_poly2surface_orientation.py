import pytest

from yapcad.geom import point, epsilon
from yapcad.geom3d import poly2surfaceXY, poly2surface
from yapcad.poly import Rect
from yapcad.combine import Boolean
from yapcad.geom_util import geomlist2poly_with_holes


def test_poly2surfacexy_generates_ccw_triangles():
    poly = [
        point(0, 0),
        point(4, 0),
        point(4, 4),
        point(2, 4),
        point(2, 2),
        point(1, 2),
        point(1, 4),
        point(0, 4),
        point(0, 0),
    ]

    surf, _ = poly2surfaceXY(poly)

    verts = surf[1]
    for face in surf[3]:
        p0, p1, p2 = [verts[i] for i in face]
        area = ((p1[0] - p0[0]) * (p2[1] - p0[1]) -
                (p2[0] - p0[0]) * (p1[1] - p0[1])) / 2.0
        assert area >= -epsilon


def test_poly2surfacexy_with_explicit_hole():
    outer = [
        point(0, 0),
        point(6, 0),
        point(6, 6),
        point(0, 6),
        point(0, 0),
    ]

    hole = [
        point(2, 2),
        point(2, 4),
        point(4, 4),
        point(4, 2),
        point(2, 2),
    ]

    surf, _ = poly2surfaceXY(outer, holepolys=[hole])

    expected_area = 36 - 4
    verts = surf[1]
    faces = surf[3]
    signed_area = 0.0
    for face in faces:
        p0, p1, p2 = [verts[i] for i in face]
        tri_area = ((p1[0] - p0[0]) * (p2[1] - p0[1]) -
                    (p2[0] - p0[0]) * (p1[1] - p0[1])) / 2.0
        assert tri_area >= -epsilon
        signed_area += tri_area
    assert signed_area == pytest.approx(expected_area, rel=1e-6)


def test_boolean_difference_surface_has_hole():
    outer = Rect(10, 10)
    inner = Rect(4, 4)

    shape = Boolean('difference', [outer, inner])
    outer_poly, hole_polys = geomlist2poly_with_holes(shape.geom)
    surf = poly2surface(outer_poly, holepolys=hole_polys)[0]

    verts = surf[1]
    faces = surf[3]
    signed_area = 0.0
    for face in faces:
        p0, p1, p2 = [verts[i] for i in face]
        tri_area = ((p1[0] - p0[0]) * (p2[1] - p0[1]) -
                    (p2[0] - p0[0]) * (p1[1] - p0[1])) / 2.0
        assert tri_area >= -epsilon
        signed_area += tri_area
    expected_area = 100 - 16
    assert signed_area == pytest.approx(expected_area, rel=1e-6)
