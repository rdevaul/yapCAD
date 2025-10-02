import io

from yapcad.geom import point, vect
from yapcad.geom3d import poly2surfaceXY, surface, solid, translatesolid
from yapcad.geom3d_util import extrude
from yapcad.io.step import write_step


def _make_surface():
    poly = [
        point(0, 0),
        point(1, 0),
        point(0, 1),
        point(0, 0),
    ]
    surf, _ = poly2surfaceXY(poly)
    return surf


def test_write_step_basic():
    surf = _make_surface()
    buf = io.StringIO()
    write_step(surf, buf, name='test_step')
    data = buf.getvalue()
    assert data.startswith('ISO-10303-21;')
    assert 'OPEN_SHELL' in data
    assert 'SHELL_BASED_SURFACE_MODEL' in data
    assert 'MANIFOLD_SOLID_BREP' not in data
    assert 'ADVANCED_FACE' in data
    assert 'test_step' in data


def _make_box(origin=None):
    if origin is None:
        origin = (0.0, 0.0, 0.0)
    x0, y0, z0 = origin
    base = [
        point(x0, y0, z0),
        point(x0 + 1.0, y0, z0),
        point(x0 + 1.0, y0 + 1.0, z0),
        point(x0, y0 + 1.0, z0),
    ]
    normals = [vect(0, 0, 1, 0) for _ in base]
    faces = [[0, 1, 2], [0, 2, 3]]
    boundary = [0, 1, 2, 3]
    base_surface = surface(base, normals, faces, boundary, [])
    return extrude(base_surface, distance=1.0)


def test_write_step_solid_component():
    box = _make_box()
    buf = io.StringIO()
    write_step(box, buf, name='solid_box')
    data = buf.getvalue()
    assert 'MANIFOLD_SOLID_BREP' in data
    assert data.count('MANIFOLD_SOLID_BREP') == 1
    assert 'CLOSED_SHELL' in data
    assert 'OPEN_SHELL' not in data
    assert max(len(line) for line in data.splitlines()) <= 256


def test_write_step_multiple_components():
    box1 = _make_box()
    box2 = translatesolid(_make_box(), point(2.5, 0.0, 0.0))
    assembly = solid(box1[1] + box2[1])
    buf = io.StringIO()
    write_step(assembly, buf, name='assembly')
    data = buf.getvalue()
    assert data.count('MANIFOLD_SOLID_BREP') == 2
    assert 'SHELL_BASED_SURFACE_MODEL' not in data
    assert max(len(line) for line in data.splitlines()) <= 256
