"""STEP import utilities backed by pythonocc-core."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, List

try:  # pragma: no cover - runtime guarded
    from OCC.Core.STEPControl import STEPControl_Reader, STEPControl_AsIs
    from OCC.Core.IFSelect import IFSelect_RetDone
    from OCC.Core.TopAbs import TopAbs_SOLID
    from OCC.Core.TopExp import TopExp_Explorer
    from OCC.Core.TopoDS import topods
except ImportError:  # pragma: no cover
    STEPControl_Reader = STEPControl_AsIs = IFSelect_RetDone = None
    TopAbs_SOLID = None
    TopExp_Explorer = None
    topods = None

from yapcad.brep import BrepSolid, require_occ
from yapcad.geometry import Geometry


def _iter_solids(shape) -> Iterable:
    """Yield TopoDS_Solid objects contained in ``shape``."""
    if shape.ShapeType() == TopAbs_SOLID:
        yield topods.Solid(shape)
        return
    explorer = TopExp_Explorer(shape, TopAbs_SOLID)
    while explorer.More():
        yield topods.Solid(explorer.Current())
        explorer.Next()


def import_step(path) -> List[Geometry]:
    """
    Read a STEP file and return a list of yapCAD Geometry objects.

    :param path: filesystem path or string to the STEP file.
    """
    require_occ()
    reader = STEPControl_Reader()
    status = reader.ReadFile(str(Path(path)))
    if status != IFSelect_RetDone:
        raise RuntimeError(f"STEP import failed with status {status}")
    if reader.TransferRoots() != IFSelect_RetDone:
        raise RuntimeError("STEP transfer failed")
    geometries: List[Geometry] = []
    for idx in range(1, reader.NbShapes() + 1):
        shape = reader.Shape(idx)
        for solid in _iter_solids(shape):
            geometries.append(Geometry(BrepSolid(solid)))
    return geometries


__all__ = ["import_step"]
