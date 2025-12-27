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
    """Perform an OCC boolean, falling back to mesh engine if BREP unavailable.

    Handles both single-solid results and compound results (e.g., union of
    disconnected solids produces a compound containing multiple solids).

    For union operations, includes sanity check: if the result volume is
    significantly less than expected (more than 10% loss), falls back to
    creating a simple compound containing both shapes. This works around
    OCC boolean issues with certain lofted/swept shapes.

    If either solid lacks BREP metadata, automatically falls back to the
    native mesh boolean engine rather than raising an error.
    """
    import os
    from OCC.Core.GProp import GProp_GProps
    from OCC.Core.BRepGProp import brepgprop

    _debug = os.environ.get('DEBUG_OCC_BOOLEAN', '').lower() in ('1', 'true', 'yes')

    if not is_available():
        raise RuntimeError("OCC booleans unavailable; activate yapcad-brep environment")
    brep_a = brep_from_solid(a)
    brep_b = brep_from_solid(b)
    if brep_a is None or brep_b is None:
        # Fall back to native mesh engine when BREP metadata is missing
        if _debug:
            print(f'[OCC_BOOL] BREP missing (a={brep_a is not None}, b={brep_b is not None}), '
                  f'falling back to native mesh engine')
        from yapcad.boolean import native as native_engine
        return native_engine.solid_boolean(a, b, operation)

    # Get input volumes for sanity check
    props_a = GProp_GProps()
    brepgprop.VolumeProperties(brep_a.shape, props_a)
    vol_a = abs(props_a.Mass())

    props_b = GProp_GProps()
    brepgprop.VolumeProperties(brep_b.shape, props_b)
    vol_b = abs(props_b.Mass())

    if _debug:
        print(f'[OCC_BOOL] Input A: type={brep_a.shape.ShapeType()}, vol={vol_a:.2f}')
        print(f'[OCC_BOOL] Input B: type={brep_b.shape.ShapeType()}, vol={vol_b:.2f}')

    builder = _make_builder(operation.lower(), brep_a, brep_b)
    builder.Build()
    if not builder.IsDone():
        raise RuntimeError("OCC boolean operation failed")
    result_shape = builder.Shape()

    props_r = GProp_GProps()
    brepgprop.VolumeProperties(result_shape, props_r)
    vol_r = abs(props_r.Mass())

    if _debug:
        print(f'[OCC_BOOL] Result: type={result_shape.ShapeType()}, vol={vol_r:.2f}')

    # For union operations, sanity check the result volume
    # If we lost more than 10% of the expected volume, fall back to compound
    if operation.lower() == 'union':
        expected_vol = vol_a + vol_b  # Upper bound (no overlap)
        min_expected_vol = max(vol_a, vol_b)  # Lower bound (complete overlap)
        if vol_r < min_expected_vol * 0.9:
            if _debug:
                print(f'[OCC_BOOL] Volume sanity check failed! '
                      f'Result {vol_r:.2f} < min expected {min_expected_vol * 0.9:.2f}')
                print(f'[OCC_BOOL] Falling back to simple compound')
            # OCC fuse failed silently - create compound instead
            compound = TopoDS_Compound()
            builder_compound = BRep_Builder()
            builder_compound.MakeCompound(compound)
            # Add all solids from both inputs
            exp_a = TopExp_Explorer(brep_a.shape, TopAbs_SOLID)
            while exp_a.More():
                builder_compound.Add(compound, exp_a.Current())
                exp_a.Next()
            exp_b = TopExp_Explorer(brep_b.shape, TopAbs_SOLID)
            while exp_b.More():
                builder_compound.Add(compound, exp_b.Current())
                exp_b.Next()
            result_shape = compound
            if _debug:
                props_c = GProp_GProps()
                brepgprop.VolumeProperties(result_shape, props_c)
                print(f'[OCC_BOOL] Compound vol: {abs(props_c.Mass()):.2f}')

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
