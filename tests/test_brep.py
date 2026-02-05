import math

import pytest

from yapcad.geom import point
from yapcad.geometry import Geometry
from yapcad.geom3d import solid, translatesolid, rotatesolid, solid_boolean
from yapcad.geom3d_util import (
    sphere,
    extrude,
    poly2surfaceXY,
    rectangularPlane,
    prism,
    conic,
    tube,
    conic_tube,
    makeRevolutionSolid,
    makeRevolutionThetaSamplingSurface,
    spherical_shell,
    makeLoftSolid,
)
from yapcad.io.geometry_json import geometry_from_json, geometry_to_json
from yapcad.brep import (
    BrepSolid,
    attach_brep_to_solid,
    brep_from_solid,
    has_brep_data,
    occ_available,
    fillet_all_edges,
    chamfer_all_edges,
)

try:  # pragma: no cover - exercised when pythonocc-core is installed
    from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
except ImportError:  # pragma: no cover
    BRepPrimAPI_MakeBox = None

pytestmark = pytest.mark.skipif(not occ_available(), reason="pythonocc-core not available")


def _make_brep_box():
    if BRepPrimAPI_MakeBox is None:
        pytest.skip("pythonocc-core not available")
    shape = BRepPrimAPI_MakeBox(10, 10, 10).Shape()
    return Geometry(BrepSolid(shape))


def _brep_box_solid():
    if BRepPrimAPI_MakeBox is None:
        pytest.skip("pythonocc-core not available")
    shape = BRepPrimAPI_MakeBox(5, 5, 5).Shape()
    brep = BrepSolid(shape)
    surface = brep.tessellate()
    sld = solid([surface])
    attach_brep_to_solid(sld, brep)
    return sld


def test_brep_tessellate():
    """BrepSolid objects should tessellate into a yapCAD surface."""
    geo = _make_brep_box()
    surface = geo.surface()
    assert isinstance(surface, list)
    assert len(surface) == 6
    assert surface[0] == 'surface'
    assert isinstance(surface[1], list)
    assert isinstance(surface[2], list)
    assert isinstance(surface[3], list)


def test_brep_translate_updates_center():
    geo = _make_brep_box()
    original_center = geo.center
    geo.translate(point(5, 0, 0))
    cx = geo.center[0]
    assert math.isclose(cx, original_center[0] + 5.0, rel_tol=1e-9)


def test_brep_uniform_scale_updates_bbox():
    geo = _make_brep_box()
    geo.scale(2.0)
    bbox = geo.bbox
    assert math.isclose(bbox[1][0], 20.0, rel_tol=1e-9, abs_tol=1e-6)
    assert math.isclose(bbox[1][1], 20.0, rel_tol=1e-9, abs_tol=1e-6)
    assert math.isclose(bbox[1][2], 20.0, rel_tol=1e-9, abs_tol=1e-6)


def test_brep_anisotropic_scale_not_supported():
    geo = _make_brep_box()
    with pytest.raises(NotImplementedError):
        geo.scale(2.0, sy=1.0)


def test_brep_rotate_about_z_changes_bbox():
    geo = _make_brep_box()
    geo.rotate(90.0, cent=point(0, 0, 0), axis=point(0, 0, 1))
    bbox = geo.bbox
    assert math.isclose(bbox[0][0], -10.0, rel_tol=1e-9, abs_tol=1e-6)
    assert math.isclose(bbox[1][0], 0.0, rel_tol=1e-9, abs_tol=1e-6)


def test_brep_mirror_flips_bbox():
    geo = _make_brep_box()
    geo.translate(point(10, 0, 0))
    geo.mirror('yz')
    bbox = geo.bbox
    assert bbox[1][0] <= 0.0
    assert bbox[0][0] < 0.0


def test_brep_transform_not_supported():
    geo = _make_brep_box()
    with pytest.raises(NotImplementedError):
        geo.transform([[1, 0, 0], [0, 1, 0], [0, 0, 1]])


def test_brep_metadata_roundtrip():
    solid_obj = _brep_box_solid()
    restored = brep_from_solid(solid_obj)
    assert restored is not None
    assert has_brep_data(solid_obj)


def test_brep_serialization_roundtrip():
    solid_obj = _brep_box_solid()
    doc = geometry_to_json([solid_obj])
    reconstructed = geometry_from_json(doc)
    assert reconstructed
    restored = brep_from_solid(reconstructed[0])
    assert restored is not None


def test_occ_boolean_union():
    solid_a = _brep_box_solid()
    solid_b = translatesolid(_brep_box_solid(), point(2.0, 0.0, 0.0))
    result = solid_boolean(solid_a, solid_b, 'union', engine='occ')
    assert has_brep_data(result)


def test_sphere_attaches_brep():
    if not occ_available():
        pytest.skip("pythonocc-core not available")
    sld = sphere(10.0)
    assert has_brep_data(sld)


def test_transforms_preserve_brep():
    if not occ_available():
        pytest.skip("pythonocc-core not available")
    sld = sphere(6.0)
    moved = translatesolid(sld, point(5.0, 0.0, 0.0))
    assert has_brep_data(moved)
    rotated = rotatesolid(sld, 45.0, axis=point(0.0, 0.0, 1.0))
    assert has_brep_data(rotated)


def test_extrude_attaches_brep():
    if not occ_available():
        pytest.skip("pythonocc-core not available")
    surf = rectangularPlane(1.0, 1.0)
    sld = extrude(surf, 2.0)
    assert has_brep_data(sld)


def test_prism_attaches_brep():
    if not occ_available():
        pytest.skip("pythonocc-core not available")
    sld = prism(2.0, 3.0, 4.0)
    assert has_brep_data(sld)


def test_conic_attaches_brep():
    if not occ_available():
        pytest.skip("pythonocc-core not available")
    sld = conic(1.0, 0.5, 3.0)
    assert has_brep_data(sld)


def test_tube_attaches_brep():
    if not occ_available():
        pytest.skip("pythonocc-core not available")
    sld = tube(2.0, 0.5, 4.0)
    assert has_brep_data(sld)


def test_conic_tube_attaches_brep():
    if not occ_available():
        pytest.skip("pythonocc-core not available")
    sld = conic_tube(3.0, 2.0, 0.25, 5.0)
    assert has_brep_data(sld)


def test_revolution_solid_attaches_brep():
    if not occ_available():
        pytest.skip("pythonocc-core not available")

    def contour(z):
        return 1.0 + 0.2 * z

    sld = makeRevolutionSolid(contour, 0.0, 2.0, steps=12, arcSamples=48)
    assert has_brep_data(sld)

    moved = translatesolid(sld, point(2.0, 0.0, 0.0))
    result = solid_boolean(sld, moved, 'union', engine='occ')
    assert has_brep_data(result)


def test_revolution_occ_boolean_union():
    if not occ_available():
        pytest.skip("pythonocc-core not available")

    def contour(z):
        return 1.0  # unit cylinder profile

    rev_solid = makeRevolutionSolid(contour, 0.0, 1.0, steps=12, arcSamples=48)
    box = prism(0.6, 0.6, 1.0)
    box = translatesolid(box, point(0.0, 0.0, 0.0))

    result = solid_boolean(rev_solid, box, 'union', engine='occ')
    assert has_brep_data(result)


def test_theta_revolution_axisymmetric_occ_boolean():
    if not occ_available():
        pytest.skip("pythonocc-core not available")

    def contour(z0, z1, theta):
        # Axisymmetric profile independent of theta
        return [(z0, 1.0), (z1, 1.0)]

    surf, brep_shape = makeRevolutionThetaSamplingSurface(
        contour, 0.0, 2.0, arcSamples=32, return_brep=True
    )
    sld = solid([surf])
    if brep_shape is not None:
        attach_brep_to_solid(sld, BrepSolid(brep_shape))

    box = prism(1.0, 1.0, 0.5)
    box = translatesolid(box, point(0.0, 0.0, 1.0))

    result = solid_boolean(sld, box, 'union', engine='occ')
    assert has_brep_data(result)


def test_extrude_with_hole_attaches_brep_and_occ_boolean():
    if not occ_available():
        pytest.skip("pythonocc-core not available")
    # Outer square with inner hole
    outer = [
        point(-1, -1, 0),
        point(1, -1, 0),
        point(1, 1, 0),
        point(-1, 1, 0),
        point(-1, -1, 0),
    ]
    inner = [
        point(-0.3, -0.3, 0),
        point(0.3, -0.3, 0),
        point(0.3, 0.3, 0),
        point(-0.3, 0.3, 0),
        point(-0.3, -0.3, 0),
    ]
    surf, _ = poly2surfaceXY(outer, holepolys=[inner])
    sld = extrude(surf, 1.0)
    assert has_brep_data(sld)

    blocker = prism(0.4, 0.4, 1.0)
    blocker = translatesolid(blocker, point(0.0, 0.0, 0.5))
    result = solid_boolean(sld, blocker, 'union', engine='occ')
    assert has_brep_data(result)


def test_spherical_shell_brep_and_occ_boolean():
    if not occ_available():
        pytest.skip("pythonocc-core not available")
    shell = spherical_shell(outer_diameter=2.0, wall_thickness=0.2, solid_angle=math.pi)
    assert has_brep_data(shell)
    cutter = conic(0.3, 0.3, 1.0, center=point(0.0, 0.0, -0.5))
    result = solid_boolean(shell, cutter, 'difference', engine='occ')
    assert has_brep_data(result)


def test_loft_solid_brep_and_occ_boolean():
    if not occ_available():
        pytest.skip("pythonocc-core not available")
    lower = [
        point(-1, -1, 0),
        point(1, -1, 0),
        point(1, 1, 0),
        point(-1, 1, 0),
    ]
    upper = [
        point(-0.5, -0.5, 1.0),
        point(1.5, -0.5, 1.0),
        point(1.5, 1.5, 1.0),
        point(-0.5, 1.5, 1.0),
    ]
    sld = makeLoftSolid(lower, upper)
    assert has_brep_data(sld)

    cutter = prism(0.5, 0.5, 0.5)
    cutter = translatesolid(cutter, point(0.0, 0.0, 0.25))
    result = solid_boolean(sld, cutter, 'union', engine='occ')
    assert has_brep_data(result)


def test_fillet_all_edges_on_box():
    """Test fillet_all_edges on a simple box."""
    if not occ_available():
        pytest.skip("pythonocc-core not available")
    if BRepPrimAPI_MakeBox is None:
        pytest.skip("pythonocc-core not available")

    # Create a 10x10x10 box
    shape = BRepPrimAPI_MakeBox(10.0, 10.0, 10.0).Shape()
    brep = BrepSolid(shape)

    # Apply fillet with radius 1.0
    filleted = fillet_all_edges(brep, 1.0)

    # Verify we got a valid shape back
    assert filleted is not None
    assert isinstance(filleted, BrepSolid)

    # Verify tessellation works
    surface = filleted.tessellate()
    assert surface[0] == 'surface'
    assert len(surface[1]) > 0  # Has vertices
    assert len(surface[3]) > 0  # Has triangles


def test_chamfer_all_edges_on_box():
    """Test chamfer_all_edges on a simple box."""
    if not occ_available():
        pytest.skip("pythonocc-core not available")
    if BRepPrimAPI_MakeBox is None:
        pytest.skip("pythonocc-core not available")

    # Create a 10x10x10 box
    shape = BRepPrimAPI_MakeBox(10.0, 10.0, 10.0).Shape()
    brep = BrepSolid(shape)

    # Apply chamfer with distance 0.5
    chamfered = chamfer_all_edges(brep, 0.5)

    # Verify we got a valid shape back
    assert chamfered is not None
    assert isinstance(chamfered, BrepSolid)

    # Verify tessellation works
    surface = chamfered.tessellate()
    assert surface[0] == 'surface'
    assert len(surface[1]) > 0  # Has vertices
    assert len(surface[3]) > 0  # Has triangles


def test_fillet_on_prism():
    """Test fillet on a prism created via yapCAD's prism function."""
    if not occ_available():
        pytest.skip("pythonocc-core not available")

    # Create a prism which already has BREP attached
    sld = prism(20.0, 10.0, 5.0)
    assert has_brep_data(sld)

    # Get the BREP and fillet it
    brep = brep_from_solid(sld)
    filleted = fillet_all_edges(brep, 0.5)

    # Verify result
    assert filleted is not None
    surface = filleted.tessellate()
    assert surface[0] == 'surface'


def test_chamfer_on_prism():
    """Test chamfer on a prism created via yapCAD's prism function."""
    if not occ_available():
        pytest.skip("pythonocc-core not available")

    # Create a prism which already has BREP attached
    sld = prism(20.0, 10.0, 5.0)
    assert has_brep_data(sld)

    # Get the BREP and chamfer it
    brep = brep_from_solid(sld)
    chamfered = chamfer_all_edges(brep, 0.3)

    # Verify result
    assert chamfered is not None
    surface = chamfered.tessellate()
    assert surface[0] == 'surface'


def test_fillet_preserves_solid_creation():
    """Test that filleted BREP can be attached to a new yapCAD solid."""
    if not occ_available():
        pytest.skip("pythonocc-core not available")
    if BRepPrimAPI_MakeBox is None:
        pytest.skip("pythonocc-core not available")

    # Create and fillet a box
    shape = BRepPrimAPI_MakeBox(10.0, 10.0, 10.0).Shape()
    brep = BrepSolid(shape)
    filleted = fillet_all_edges(brep, 1.0)

    # Create a yapCAD solid from the filleted BREP
    surface = filleted.tessellate()
    sld = solid([surface])
    attach_brep_to_solid(sld, filleted)

    # Verify BREP data is preserved
    assert has_brep_data(sld)

    # Verify we can get BREP back
    restored = brep_from_solid(sld)
    assert restored is not None
