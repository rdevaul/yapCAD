import io
import os
import tempfile

import pytest

from yapcad.geom import point, vect
from yapcad.geom3d import poly2surfaceXY, surface, solid, translatesolid
from yapcad.geom3d_util import extrude, prism
from yapcad.io.step import write_step, write_step_analytic

# Check if OCC is available for analytic STEP tests
try:
    from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder
    _OCC_AVAILABLE = True
except ImportError:
    _OCC_AVAILABLE = False


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


def _make_mesh_only_box():
    """Create a minimal mesh-only solid without any BREP data.

    This creates a solid directly using low-level surface/solid functions,
    bypassing prism() and extrude() which now attach BREP data.
    Used for testing fallback behavior when BREP is unavailable.
    """
    # Create bottom face (two triangles)
    verts_bottom = [point(0, 0, 0), point(1, 0, 0), point(1, 1, 0), point(0, 1, 0)]
    norms_bottom = [vect(0, 0, -1, 0)] * 4
    faces_bottom = [[0, 2, 1], [0, 3, 2]]
    bottom = surface(verts_bottom, norms_bottom, faces_bottom, [0, 1, 2, 3], [])

    # Create top face (two triangles)
    verts_top = [point(0, 0, 1), point(1, 0, 1), point(1, 1, 1), point(0, 1, 1)]
    norms_top = [vect(0, 0, 1, 0)] * 4
    faces_top = [[0, 1, 2], [0, 2, 3]]
    top = surface(verts_top, norms_top, faces_top, [0, 1, 2, 3], [])

    return solid([bottom, top])


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


# =========================================================================
# Analytic STEP Export Tests
# =========================================================================

@pytest.mark.skipif(not _OCC_AVAILABLE, reason="OCC not available")
def test_write_step_analytic_from_brep_solid():
    """Test analytic STEP export from imported BrepSolid."""
    from yapcad.io.step_importer import import_step
    from OCC.Core.STEPControl import STEPControl_Writer, STEPControl_AsIs

    # Create OCC box and write to temp STEP file
    occ_box = BRepPrimAPI_MakeBox(10, 10, 10).Shape()
    with tempfile.NamedTemporaryFile(suffix='.step', delete=False) as f:
        temp_step = f.name

    writer = STEPControl_Writer()
    writer.Transfer(occ_box, STEPControl_AsIs)
    writer.Write(temp_step)

    try:
        # Import (returns Geometry(BrepSolid))
        parts = import_step(temp_step)
        assert len(parts) > 0

        # Export using analytic
        with tempfile.NamedTemporaryFile(suffix='.step', delete=False) as f:
            out_step = f.name

        result = write_step_analytic(parts[0], out_step, name='TestBox')
        assert result is True

        # Verify output
        with open(out_step, 'r') as f:
            content = f.read()

        assert 'PLANE(' in content
        assert 'ADVANCED_BREP' in content
        assert 'MANIFOLD_SOLID_BREP' in content

        os.unlink(out_step)
    finally:
        os.unlink(temp_step)


@pytest.mark.skipif(not _OCC_AVAILABLE, reason="OCC not available")
def test_write_step_analytic_from_native_brep():
    """Test analytic STEP export from yapCAD solid with native BREP."""
    from yapcad.native_brep import attach_native_brep_to_solid, has_native_brep
    from yapcad.occ_native_convert import occ_solid_to_native_brep

    # Create a yapCAD prism solid
    solid = prism(10, 10, 10)
    assert not has_native_brep(solid)

    # Create OCC box and convert to native BREP
    occ_box = BRepPrimAPI_MakeBox(10, 10, 10).Solid()
    _, graph = occ_solid_to_native_brep(occ_box)

    # Attach native BREP
    attach_native_brep_to_solid(solid, graph)
    assert has_native_brep(solid)

    # Export using analytic
    with tempfile.NamedTemporaryFile(suffix='.step', delete=False) as f:
        out_step = f.name

    result = write_step_analytic(solid, out_step, name='NativeBox')
    assert result is True

    # Verify output
    with open(out_step, 'r') as f:
        content = f.read()

    assert 'PLANE(' in content
    assert 'ADVANCED_BREP' in content

    os.unlink(out_step)


@pytest.mark.skipif(not _OCC_AVAILABLE, reason="OCC not available")
def test_write_step_analytic_cylinder():
    """Test analytic STEP export with curved surfaces."""
    from yapcad.io.step_importer import import_step
    from OCC.Core.STEPControl import STEPControl_Writer, STEPControl_AsIs

    # Create OCC cylinder
    occ_cyl = BRepPrimAPI_MakeCylinder(5.0, 20.0).Shape()
    with tempfile.NamedTemporaryFile(suffix='.step', delete=False) as f:
        temp_step = f.name

    writer = STEPControl_Writer()
    writer.Transfer(occ_cyl, STEPControl_AsIs)
    writer.Write(temp_step)

    try:
        # Import
        parts = import_step(temp_step)

        # Export using analytic
        with tempfile.NamedTemporaryFile(suffix='.step', delete=False) as f:
            out_step = f.name

        result = write_step_analytic(parts[0], out_step, name='Cylinder')
        assert result is True

        # Verify output contains cylindrical surface
        with open(out_step, 'r') as f:
            content = f.read()

        assert 'CYLINDRICAL_SURFACE(' in content
        assert 'PLANE(' in content  # End caps
        assert 'ADVANCED_BREP' in content

        os.unlink(out_step)
    finally:
        os.unlink(temp_step)


@pytest.mark.skipif(not _OCC_AVAILABLE, reason="OCC not available")
def test_write_step_analytic_fallback():
    """Test that analytic export falls back to faceted when no BREP data."""
    # Create a mesh-only solid (no embedded BREP)
    # Note: prism() and extrude() now attach BREP data, so we use
    # _make_mesh_only_box() which creates a solid via low-level functions
    mesh_only_solid = _make_mesh_only_box()

    # Export using analytic with fallback enabled
    with tempfile.NamedTemporaryFile(suffix='.step', delete=False) as f:
        out_step = f.name

    result = write_step_analytic(mesh_only_solid, out_step, fallback_to_faceted=True)
    assert result is False  # Fell back to faceted

    # Verify file was created
    assert os.path.exists(out_step)
    os.unlink(out_step)


@pytest.mark.skipif(not _OCC_AVAILABLE, reason="OCC not available")
def test_write_step_analytic_no_fallback_error():
    """Test that analytic export raises error when no BREP and fallback disabled."""
    # Create a mesh-only solid (no embedded BREP)
    # Note: prism() and extrude() now attach BREP data, so we use
    # _make_mesh_only_box() which creates a solid via low-level functions
    mesh_only_solid = _make_mesh_only_box()

    with tempfile.NamedTemporaryFile(suffix='.step', delete=False) as f:
        out_step = f.name

    with pytest.raises(ValueError, match='no native BREP data'):
        write_step_analytic(mesh_only_solid, out_step, fallback_to_faceted=False)

    # Clean up if file was created
    if os.path.exists(out_step):
        os.unlink(out_step)
