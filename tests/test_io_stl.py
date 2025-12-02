import io
import struct

from yapcad.geom import point
from yapcad.geom3d import poly2surfaceXY, issolid, issurface
from yapcad.geom3d_util import prism
from yapcad.io.stl import write_stl, read_stl, import_stl


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


# ---------------------------------------------------------------------------
# STL Import Tests
# ---------------------------------------------------------------------------


def test_read_stl_binary_roundtrip(tmp_path):
    """Test binary STL round-trip: export then import."""
    box = prism(2, 2, 2)
    path = tmp_path / 'box.stl'
    write_stl(box, path, binary=True)

    imported = read_stl(path)
    assert issolid(imported)
    assert len(imported[1]) == 1  # One surface
    surf = imported[1][0]
    assert issurface(surf)
    # A 2x2x2 box has 8 vertices and 12 triangles (2 per face * 6 faces)
    assert len(surf[1]) == 8  # Deduplicated vertices
    assert len(surf[3]) == 12  # Triangular faces


def test_read_stl_ascii_roundtrip(tmp_path):
    """Test ASCII STL round-trip: export then import."""
    box = prism(2, 2, 2)
    path = tmp_path / 'box_ascii.stl'
    write_stl(box, path, binary=False)

    imported = read_stl(path)
    assert issolid(imported)
    assert len(imported[1]) == 1
    surf = imported[1][0]
    assert len(surf[1]) == 8
    assert len(surf[3]) == 12


def test_read_stl_no_deduplicate(tmp_path):
    """Test import without vertex deduplication."""
    box = prism(2, 2, 2)
    path = tmp_path / 'box.stl'
    write_stl(box, path, binary=True)

    imported = read_stl(path, deduplicate=False)
    assert issolid(imported)
    surf = imported[1][0]
    # Without deduplication: 12 triangles * 3 vertices = 36 vertices
    assert len(surf[1]) == 36
    assert len(surf[3]) == 12


def test_import_stl_alias(tmp_path):
    """Test import_stl function alias."""
    box = prism(1, 1, 1)
    path = tmp_path / 'cube.stl'
    write_stl(box, path)

    imported = import_stl(str(path))
    assert issolid(imported)


def test_read_stl_surface(tmp_path):
    """Test importing a single triangle surface."""
    surf = _make_surface()
    path = tmp_path / 'tri.stl'
    write_stl(surf, path, binary=True)

    imported = read_stl(path)
    assert issolid(imported)
    imported_surf = imported[1][0]
    # Single triangle: 3 vertices, 1 face
    assert len(imported_surf[1]) == 3
    assert len(imported_surf[3]) == 1


def test_read_stl_from_file_object(tmp_path):
    """Test reading from an open file object."""
    box = prism(1, 1, 1)
    path = tmp_path / 'box.stl'
    write_stl(box, path, binary=True)

    with open(path, 'rb') as f:
        imported = read_stl(f)
    assert issolid(imported)


def test_read_stl_empty(tmp_path):
    """Test that empty STL returns empty solid."""
    # Create minimal valid binary STL with 0 triangles
    path = tmp_path / 'empty.stl'
    with open(path, 'wb') as f:
        f.write(b' ' * 80)  # Header
        f.write(struct.pack('<I', 0))  # 0 triangles

    imported = read_stl(path)
    assert issolid(imported)
    # Empty solid has empty surfaces list
    assert len(imported[1]) == 0
