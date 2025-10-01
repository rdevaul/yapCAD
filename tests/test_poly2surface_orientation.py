import pytest

from yapcad.geom import point, epsilon
from yapcad.geom3d import poly2surfaceXY, poly2surface
from yapcad.poly import Rect
from yapcad.combine import Boolean
from yapcad.geom_util import geomlist2poly_with_holes, geomlist2poly_components


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


def _surface_signed_area(surface):
    verts = surface[1]
    faces = surface[3]
    total = 0.0
    for face in faces:
        p0, p1, p2 = [verts[i] for i in face]
        tri_area = ((p1[0] - p0[0]) * (p2[1] - p0[1]) -
                    (p2[0] - p0[0]) * (p1[1] - p0[1])) / 2.0
        total += tri_area
    return total


def _polygon_signed_area(loop):
    total = 0.0
    for idx in range(1, len(loop)):
        x0, y0 = loop[idx - 1][0], loop[idx - 1][1]
        x1, y1 = loop[idx][0], loop[idx][1]
        total += x0 * y1 - x1 * y0
    return total / 2.0


def test_poly2surfacexy_with_multiple_holes():
    outer = [
        point(0, 0),
        point(10, 0),
        point(10, 6),
        point(0, 6),
        point(0, 0),
    ]

    hole_a = [
        point(2, 2),
        point(4, 2),
        point(4, 4),
        point(2, 4),
        point(2, 2),
    ]

    hole_b = [
        point(6, 1),
        point(8, 1),
        point(8, 3),
        point(6, 3),
        point(6, 1),
    ]

    surf, _ = poly2surfaceXY(outer, holepolys=[hole_a, hole_b])
    signed_area = _surface_signed_area(surf)
    expected = 60 - 4 - 4
    assert signed_area == pytest.approx(expected, rel=1e-6)


def test_geomlist2poly_components_disjoint_union():
    outer1 = Rect(4, 4)
    outer2 = Rect(2, 2, center=point(6, 0))

    shape = Boolean('union', [outer1, outer2])
    components = geomlist2poly_components(shape.geom)

    assert len(components) == 2

    areas = []
    for outer_poly, holes in components:
        assert holes == []
        surf, _ = poly2surfaceXY(outer_poly, holepolys=holes)
        areas.append(_surface_signed_area(surf))

    assert pytest.approx(sum(areas), rel=1e-6) == 16 + 4


def test_geomlist2poly_components_nested_island():
    frame = Boolean('difference', [Rect(10, 10), Rect(4, 4)])
    island = Rect(3, 3)
    combined = Boolean('union', [frame, island])

    components = geomlist2poly_components(combined.geom)
    assert len(components) == 2

    components.sort(key=lambda comp: abs(_polygon_signed_area(comp[0])), reverse=True)

    outer_poly, holes = components[0]
    assert len(holes) == 1
    outer_surface, _ = poly2surfaceXY(outer_poly, holepolys=holes)
    outer_area = _surface_signed_area(outer_surface)
    expected_outer = 100 - 16
    assert outer_area == pytest.approx(expected_outer, rel=1e-6)

    island_poly, island_holes = components[1]
    assert island_holes == []
    island_surface, _ = poly2surfaceXY(island_poly)
    island_area = _surface_signed_area(island_surface)
    assert island_area == pytest.approx(9, rel=1e-6)
