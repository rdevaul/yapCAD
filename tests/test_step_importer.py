import math
from pathlib import Path

import pytest

from yapcad.brep import occ_available

pytestmark = pytest.mark.skipif(not occ_available(), reason="pythonocc-core not available")

from yapcad.io.step_importer import import_step

try:  # pragma: no cover
    from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
    from OCC.Core.STEPControl import STEPControl_Writer, STEPControl_AsIs
    from OCC.Core.IFSelect import IFSelect_RetDone
except ImportError:  # pragma: no cover
    BRepPrimAPI_MakeBox = None
    STEPControl_Writer = None
    STEPControl_AsIs = None
    IFSelect_RetDone = None


def _write_step_box(path: Path):
    if BRepPrimAPI_MakeBox is None:
        pytest.skip("pythonocc-core not available")
    shape = BRepPrimAPI_MakeBox(5, 5, 5).Shape()
    writer = STEPControl_Writer()
    assert writer.Transfer(shape, STEPControl_AsIs) == IFSelect_RetDone
    assert writer.Write(str(path)) == IFSelect_RetDone


def test_import_step_returns_geometry(tmp_path):
    step_path = tmp_path / "box.step"
    _write_step_box(step_path)
    geoms = import_step(step_path)
    assert len(geoms) == 1
    bbox = geoms[0].bbox
    assert math.isclose(bbox[1][0], 5.0, abs_tol=1e-6)
    assert math.isclose(bbox[1][1], 5.0, abs_tol=1e-6)
    assert math.isclose(bbox[1][2], 5.0, abs_tol=1e-6)
