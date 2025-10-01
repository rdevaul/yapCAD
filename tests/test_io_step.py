import io

from yapcad.geom import point
from yapcad.geom3d import poly2surfaceXY
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
    assert 'MANIFOLD_SOLID_BREP' in data
    assert 'CLOSED_SHELL' in data
    assert 'ADVANCED_FACE' in data
    assert 'test_step' in data

