import io
import struct

from yapcad.geom import point
from yapcad.geom3d import poly2surfaceXY
from yapcad.io.stl import write_stl


def _make_surface():
    poly = [
        point(0, 0),
        point(1, 0),
        point(0, 1),
        point(0, 0),
    ]
    surf, _ = poly2surfaceXY(poly)
    return surf


def test_write_stl_binary(tmp_path):
    surf = _make_surface()
    path = tmp_path / 'tri.stl'
    write_stl(surf, path, binary=True, name='test')

    data = path.read_bytes()
    assert len(data) == 80 + 4 + 50  # header + count + one triangle
    assert data[0:4] == b'test'
    count = struct.unpack('<I', data[80:84])[0]
    assert count == 1


def test_write_stl_ascii():
    surf = _make_surface()
    buf = io.StringIO()
    write_stl(surf, buf, binary=False, name='ascii_test')

    text = buf.getvalue()
    assert 'solid ascii_test' in text
    assert 'facet normal' in text
    assert 'vertex' in text
    assert text.strip().endswith('endsolid ascii_test')
