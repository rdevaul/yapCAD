import os
import pytest

from yapcad.combine import Boolean
from yapcad.geom import arc, line, point, translate
from yapcad.geom3d import poly2surfaceXY, solid, translatesurface, surf2lines
from yapcad.geom3d_util import conic, sphere
from yapcad.geom_util import geomlist2poly, geomlist2poly_components
from yapcad.mesh import mesh_view
from yapcad.poly import Polygon, Rect

VISUALTEST = os.environ.get('VISUALTEST', 'false').lower() in ('true', '1', 'yes')
if VISUALTEST:
    from yapcad.pyglet_drawable import pygletDraw

def _strip_w(pt):
    return tuple(round(coord, 6) for coord in pt[:3])


def test_mesh_view_surface_basic():
    poly = [
        point(0, 0),
        point(1, 0),
        point(0, 1),
        point(0, 0),
    ]

    surf, _ = poly2surfaceXY(poly)

    tris = list(mesh_view(surf))
    assert len(tris) == 1

    normal, v0, v1, v2 = tris[0]
    assert normal == (0.0, 0.0, 1.0)
    assert _strip_w(surf[1][0]) in (v0, v1, v2)


def test_mesh_view_reorients_against_stored_normal():
    poly = [
        point(0, 0),
        point(1, 0),
        point(0, 1),
        point(0, 0),
    ]

    surf, _ = poly2surfaceXY(poly)
    surf[2] = [[0, 0, -1, 0] for _ in surf[2]]

    tri = next(mesh_view(surf))
    normal = tri[0]
    assert normal == (0.0, 0.0, 1.0)


def test_mesh_view_on_solid():
    poly = [
        point(0, 0),
        point(1, 0),
        point(0, 1),
        point(0, 0),
    ]

    surf, _ = poly2surfaceXY(poly)
    sld = solid([surf], [], [])

    tris = list(mesh_view(sld))
    assert len(tris) == 1


@pytest.mark.visual
def test_mesh_view_visual_normals():
    if not VISUALTEST:
        pytest.skip("Visual tests disabled (set VISUALTEST=true to enable)")

    outer = [
        point(-1, -1),
        point(1, -1),
        point(1, 1),
        point(-1, 1),
        point(-1, -1),
    ]

    hole = [
        point(-0.4, -0.4),
        point(0.4, -0.4),
        point(0.4, 0.4),
        point(-0.4, 0.4),
        point(-0.4, -0.4),
    ]

    surf, _ = poly2surfaceXY(outer, holepolys=[hole])
    sld = solid([surf], [], [])

    dd = pygletDraw()

    def draw_normals(obj, color):
        dd.linecolor = color
        for normal, v0, v1, v2 in mesh_view(obj):
            cx = (v0[0] + v1[0] + v2[0]) / 3.0
            cy = (v0[1] + v1[1] + v2[1]) / 3.0
            cz = (v0[2] + v1[2] + v2[2]) / 3.0
            start = point(cx, cy, cz + 0.01)
            end = point(
                cx + normal[0] * 0.5,
                cy + normal[1] * 0.5,
                cz + normal[2] * 0.5 + 0.01,
            )
            dd.draw(line(start, end))

    dd.linecolor = 'silver'
    dd.polystyle = 'both'
    dd.draw_surface(surf)
    draw_normals(sld, 'yellow')

    shp = sphere(1.5, center=point(3, 0, 0))
    dd.linecolor = 'white'
    dd.draw(shp)
    draw_normals(shp, 'aqua')

    cone = conic(1.0, 0.5, 2.5, center=point(-3, 0, 0))
    dd.linecolor = 'green'
    dd.draw(cone)
    draw_normals(cone, 'blue')

    dd.display()


@pytest.mark.visual
def test_mesh_view_visual_multi_component_tessellation():
    if not VISUALTEST:
        pytest.skip("Visual tests disabled (set VISUALTEST=true to enable)")

    dd = pygletDraw()
    dd.polystyle = 'both'

    multi_outer = [
        point(-4, -2),
        point(0, -2),
        point(2, 0),
        point(0, 2),
        point(-4, 2),
        point(-4, -2),
    ]

    hole_one = [
        point(-3.0, -1.0),
        point(-2.0, -1.0),
        point(-2.0, 0.0),
        point(-3.0, 0.0),
        point(-3.0, -1.0),
    ]

    hole_two = [
        point(-1.5, 0.5),
        point(-0.5, 0.5),
        point(-0.5, 1.5),
        point(-1.5, 1.5),
        point(-1.5, 0.5),
    ]

    multi_surface, _ = poly2surfaceXY(multi_outer, holepolys=[hole_one, hole_two])
    dd.linecolor = 'orange'
    dd.draw_surface(multi_surface)

    outer_a = Rect(4, 4)
    outer_b = Rect(2, 2, center=point(5.5, 0))
    union_shape = Boolean('union', [outer_a, outer_b])
    components = geomlist2poly_components(union_shape.geom)

    palette = ['cyan', 'magenta', 'yellow', 'lime']

    for idx, (outer_poly, holes) in enumerate(components):
        surf, _ = poly2surfaceXY(outer_poly, holepolys=holes)
        shifted = translatesurface(surf, point(6, idx * 3.5, 0))
        dd.linecolor = palette[idx % len(palette)]
        dd.draw_surface(shifted)

    def _circle_polygon(cx, cy, radius, segments=48):
        arc_geom = [arc(point(cx, cy), radius)]
        sampled = geomlist2poly(arc_geom, minang=360.0 / segments)
        return Polygon(sampled)

    base_square = Rect(6, 6)
    circle_a = _circle_polygon(-1.6, -0.8, 1.2)
    circle_b = _circle_polygon(1.6, 0.8, 1.0)
    carved_first = Boolean('difference', [base_square, circle_a])
    carved = Boolean('difference', [carved_first, circle_b])
    carved_components = geomlist2poly_components(carved.geom)

    if carved_components:
        carved_outer, carved_holes = carved_components[0]
        carved_surface, _ = poly2surfaceXY(carved_outer, holepolys=carved_holes)
        carved_surface = translatesurface(carved_surface, point(-6, -4, 0))
        dd.linecolor = 'lightsteelblue'
        dd.draw_surface(carved_surface)

        edge_lines = surf2lines(carved_surface)
        edge_lines = translate(edge_lines, point(0, 0, 0.15))
        dd.linecolor = 'yellow'
        dd.draw(edge_lines)

    dd.display()
