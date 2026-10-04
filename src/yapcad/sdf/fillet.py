"""Round every edge of an SDF-authored part.

``fillet(part, r)`` means what OCC's fillet-all-edges means: every sharp
edge of the part becomes a radius-``r`` round.  For the analytic primitives
that is exact and closed-form -- a box becomes a :func:`rounded_box`, a
cylinder a :func:`rounded_cylinder` -- and it stays exact through similarity
transforms, whose uniform scale divides the radius.  Both rounded
primitives replay through OCC, so a filleted primitive still exports as
analytic STEP.  Bodies of a :func:`compound` are filleted independently.

Combined fields are refused rather than approximated.  The tempting
rewrite -- round each operand, blend each boolean -- rounds edges the
boolean removes: two boxes sharing a face come out with a groove along the
seam, where OCC's fuse-then-fillet leaves a clean box.  A correct general
fillet needs the distance to the combined part recomputed, which is what
the planned sampled-grid node (re-distancing on a grid) will provide; until
then, fillet primitives before combining them, as most designs do, or use
:func:`smooth_union` and friends for blends.

Mesh-only solids are not handled here: see ``docs/SDF-DESIGN.md`` for the
plan.
"""

import math

import numpy as np

from yapcad.construction import get_construction
from yapcad.sdf.node import SdfError, from_construction, make_node
from yapcad.sdf.ops import compound
from yapcad.sdf.primitives import rounded_box, rounded_cylinder

__all__ = ["fillet", "fillet_solid"]

#: Primitives with no edges: a fillet leaves them unchanged.
_SMOOTH = {"sphere", "capsule", "torus", "rounded_box", "rounded_cylinder"}


def _uniform_scale(matrix, kind):
    """The uniform scale of a similarity transform, or raise."""
    linear = np.asarray([row[:3] for row in matrix[:3]], dtype=float)
    sv = np.linalg.svd(linear, compute_uv=False)
    if float(sv.max() - sv.min()) > 1e-9 * max(1.0, float(sv.max())):
        raise SdfError(
            f"fillet: cannot round edges through a non-uniform scale (in a "
            f"{kind}); the rounds would not be circular"
        )
    return float(sv.max())


def _round(node, radius):
    kind = node.kind
    p = node.p
    if kind == "transform":
        s = _uniform_scale(p["matrix"], kind)
        return make_node("transform", p,
                         (_round(node.children[0], radius / s),))
    if kind == "compound":
        return compound(*(_round(c, radius) for c in node.children))
    if kind in _SMOOTH:
        return node
    if kind == "box":
        if radius > 0.5 * min(p["size"]):
            raise SdfError(
                f"fillet: radius {radius:g} exceeds half the smallest edge "
                f"{0.5 * min(p['size']):g} of a box"
            )
        return rounded_box(p["size"], radius)
    if kind == "cylinder":
        if radius > p["radius"] or 2.0 * radius > p["height"]:
            raise SdfError(
                f"fillet: radius {radius:g} exceeds the radius "
                f"{p['radius']:g} or half the height {0.5 * p['height']:g} "
                f"of a cylinder"
            )
        return rounded_cylinder(p["radius"], p["height"], radius)
    if kind in ("union", "intersect", "subtract"):
        raise SdfError(
            f"fillet: rounding the edges of a combined field ({kind}) is not "
            f"supported yet. Rounding each operand would also round edges "
            f"the {kind} removes; a correct fillet needs the planned "
            f"sampled-grid re-distancing. Fillet the primitives before "
            f"combining them, or use smooth_{kind} for a blend."
        )
    raise SdfError(f"fillet: cannot round the edges of a {kind} node yet")


def fillet(node, radius):
    """``node`` with every edge rounded to ``radius``; see the module notes.

    Raises :class:`~yapcad.sdf.node.SdfError` for the cases it cannot do
    exactly -- combined fields, non-uniform scales, cones -- rather than
    returning something that is not a fillet.
    """
    try:
        r = float(radius)
    except (TypeError, ValueError):
        raise SdfError("fillet: radius must be a number") from None
    if not math.isfinite(r) or r <= 0.0:
        raise SdfError("fillet: radius must be positive and finite")
    return _round(node, r)


def fillet_solid(solid, radius):
    """Fillet an SDF-authored solid and mesh the result.

    The result is meshed no coarser than ``solid`` was, the rule the SDF
    booleans follow, and keeps a derived BREP when ``solid`` had one and
    the filleted tree still replays.
    """
    from yapcad.sdf.booleans import (
        DEFAULT_RESOLUTION,
        MAX_RESOLUTION,
        _cell_size,
        _derived_from,
        is_sdf_solid,
    )
    from yapcad.sdf.convert import to_solid

    if not is_sdf_solid(solid):
        raise SdfError("fillet_solid: the solid is not SDF-authored")
    tree = from_construction(get_construction(solid))
    rounded = fillet(tree, radius)
    cell = _cell_size(solid, tree)
    lo, hi = rounded.bounds
    longest = max(hi[i] - lo[i] for i in range(3))
    if cell and cell > 0.0 and longest > 0.0:
        resolution = math.ceil(longest / cell)
    else:
        resolution = DEFAULT_RESOLUTION
    resolution = int(min(max(resolution, 2), MAX_RESOLUTION))
    brep = "auto" if _derived_from(solid, tree) else False
    return to_solid(rounded, resolution=resolution, brep=brep)
