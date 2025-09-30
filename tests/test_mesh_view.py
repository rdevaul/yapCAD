import os
import pytest

from yapcad.geom import line, point
from yapcad.geom3d import poly2surfaceXY, solid
from yapcad.geom3d_util import conic, sphere
from yapcad.mesh import mesh_view
from yapcad.pyglet_drawable import pygletDraw


VISUALTEST = os.environ.get('VISUALTEST', 'false').lower() in ('true', '1', 'yes')


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
