import math

import pytest

from yapcad.geom import point
from yapcad.geometry import Geometry
from yapcad.brep import BrepSolid, occ_available

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
