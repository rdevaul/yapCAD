"""OCC-backed boolean engine operating on BrepSolid metadata."""

from __future__ import annotations

from typing import Optional

try:  # pragma: no cover - optional dependency
    from OCC.Core.BRepAlgoAPI import (
        BRepAlgoAPI_Common,
        BRepAlgoAPI_Cut,
        BRepAlgoAPI_Fuse,
    )
    from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeSolid
    from OCC.Core.TopExp import TopExp_Explorer
    from OCC.Core.TopAbs import TopAbs_SOLID, TopAbs_SHELL
    from OCC.Core.TopoDS import topods, TopoDS_Compound, TopoDS_Builder
    from OCC.Core.BRep import BRep_Builder
except ImportError:  # pragma: no cover
    BRepAlgoAPI_Common = BRepAlgoAPI_Cut = BRepAlgoAPI_Fuse = None
    BRepBuilderAPI_MakeSolid = None
    TopExp_Explorer = None
    TopAbs_SOLID = TopAbs_SHELL = None
    topods = TopoDS_Compound = TopoDS_Builder = BRep_Builder = None

from yapcad.brep import (
    BrepSolid,
    attach_brep_to_solid,
    brep_from_solid,
    occ_available,
    require_occ,
)


def is_available() -> bool:
    return occ_available() and BRepAlgoAPI_Fuse is not None


def _make_builder(op: str, a: BrepSolid, b: BrepSolid):
    require_occ()
    if BRepAlgoAPI_Fuse is None:
        raise RuntimeError("pythonocc-core boolean operators unavailable")
    if op == 'union':
        return BRepAlgoAPI_Fuse(a.shape, b.shape)
    if op == 'intersection':
        return BRepAlgoAPI_Common(a.shape, b.shape)
    if op == 'difference':
        return BRepAlgoAPI_Cut(a.shape, b.shape)
    raise ValueError(f"unsupported OCC boolean operation '{op}'")


def solid_boolean(a, b, operation: str):
    """Perform an OCC boolean; requires both solids to carry BREP metadata.

    Handles both single-solid results and compound results (e.g., union of
    disconnected solids produces a compound containing multiple solids).
    """
    if not is_available():
        raise RuntimeError("OCC booleans unavailable; activate yapcad-brep environment")
    brep_a = brep_from_solid(a)
    brep_b = brep_from_solid(b)
    if brep_a is None or brep_b is None:
        raise ValueError("OCC boolean requires solids with BREP metadata")
    builder = _make_builder(operation.lower(), brep_a, brep_b)
    builder.Build()
    if not builder.IsDone():
        raise RuntimeError("OCC boolean operation failed")
    result_shape = builder.Shape()

    # Collect all solids from the result (may be compound for disconnected unions)
    solids = []
    exp = TopExp_Explorer(result_shape, TopAbs_SOLID)
    while exp.More():
        solid = exp.Current()
        if topods is not None:
            solid = topods.Solid(solid)
        solids.append(solid)
        exp.Next()

    if not solids:
        raise RuntimeError("OCC boolean produced no solids")

    # If single solid, use it directly
    if len(solids) == 1:
        result_solid = solids[0]
    else:
        # Multiple solids: create a compound to hold them all
        compound = TopoDS_Compound()
        builder_compound = BRep_Builder()
        builder_compound.MakeCompound(compound)
        for solid in solids:
            builder_compound.Add(compound, solid)
        result_solid = compound

    result_brep = BrepSolid(result_solid)
    surface = result_brep.tessellate()
    from yapcad.geom3d import solid as geom3d_solid  # circular-safe import

    result = geom3d_solid([surface])
    attach_brep_to_solid(result, result_brep)
    return result


__all__ = ["solid_boolean", "is_available"]
