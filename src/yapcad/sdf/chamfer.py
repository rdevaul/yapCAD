"""Bevel every edge of an SDF-authored part.

``chamfer(part, d)`` means what OCC's chamfer-all-edges means: every sharp
edge of the part becomes a 45 degree bevel, cut back ``d`` along each face
the edge joins.  For the analytic primitives that is exact and closed form
-- a box becomes a :func:`chamfered_box`, a cylinder a
:func:`chamfered_cylinder` -- and it stays exact through similarity
transforms, whose uniform scale divides the distance.  Both chamfered
primitives replay through OCC, so a chamfered primitive still exports as
analytic STEP.  Bodies of a :func:`compound` are chamfered independently.

This is :mod:`yapcad.sdf.fillet` with a bevel in place of a round, and it
stops at the same wall for the same reason.  Combined fields are refused
rather than approximated: bevelling each operand would bevel edges the
boolean removes, leaving a groove along a seam where OCC's fuse-then-chamfer
leaves a clean face.  A correct general chamfer needs the distance to the
combined part recomputed, which is what the planned sampled-grid node
(``docs/SDF-DESIGN.md`` section 4.2) will provide; until then, chamfer the
primitives before combining them.

Mesh-only solids are not handled here either; see ``docs/SDF-DESIGN.md``.
"""

import math

from yapcad.construction import get_construction
from yapcad.sdf.fillet import _uniform_scale
from yapcad.sdf.node import SdfError, from_construction, make_node
from yapcad.sdf.ops import compound
from yapcad.sdf.primitives import chamfered_box, chamfered_cylinder

__all__ = ["chamfer", "chamfer_solid"]

#: Primitives with no sharp edge for a bevel to sit on: a chamfer leaves
#: them unchanged.  ``rounded_box`` and ``rounded_cylinder`` belong here --
#: their edges are already rounds -- which is also what OCC does with them.
_SMOOTH = {"sphere", "capsule", "torus", "rounded_box", "rounded_cylinder"}


def _bevel(node, distance):
    kind = node.kind
    p = node.p
    if kind == "transform":
        s = _uniform_scale(p["matrix"], kind, op="chamfer", verb="bevel",
                           consequence="bevels would not be at 45 degrees")
        return make_node("transform", p,
                         (_bevel(node.children[0], distance / s),))
    if kind == "compound":
        return compound(*(_bevel(c, distance) for c in node.children))
    if kind in _SMOOTH:
        return node
    if kind == "box":
        if distance > 0.5 * min(p["size"]):
            raise SdfError(
                f"chamfer: distance {distance:g} exceeds half the smallest "
                f"edge {0.5 * min(p['size']):g} of a box"
            )
        return chamfered_box(p["size"], distance)
    if kind == "cylinder":
        if distance > p["radius"] or 2.0 * distance > p["height"]:
            raise SdfError(
                f"chamfer: distance {distance:g} exceeds the radius "
                f"{p['radius']:g} or half the height "
                f"{0.5 * p['height']:g} of a cylinder"
            )
        return chamfered_cylinder(p["radius"], p["height"], distance)
    if kind in ("chamfered_box", "chamfered_cylinder"):
        raise SdfError(
            f"chamfer: a {kind} is already chamfered; bevelling the edges "
            f"its first chamfer created is not supported. Chamfer once, "
            f"with the distance you want."
        )
    if kind in ("union", "intersect", "subtract"):
        raise SdfError(
            f"chamfer: bevelling the edges of a combined field ({kind}) is "
            f"not supported yet. Bevelling each operand would also bevel "
            f"edges the {kind} removes; a correct chamfer needs the planned "
            f"sampled-grid re-distancing. Chamfer the primitives before "
            f"combining them."
        )
    raise SdfError(f"chamfer: cannot bevel the edges of a {kind} node yet")


def chamfer(node, distance):
    """``node`` with every edge bevelled by ``distance``; see the module
    notes.

    Raises :class:`~yapcad.sdf.node.SdfError` for the cases it cannot do
    exactly -- combined fields, non-uniform scales, cones -- rather than
    returning something that is not a chamfer.
    """
    try:
        d = float(distance)
    except (TypeError, ValueError):
        raise SdfError("chamfer: distance must be a number") from None
    if not math.isfinite(d) or d <= 0.0:
        raise SdfError("chamfer: distance must be positive and finite")
    return _bevel(node, d)


def chamfer_solid(solid, distance):
    """Chamfer an SDF-authored solid and mesh the result.

    The result is meshed no coarser than ``solid`` was, the rule the SDF
    booleans follow, and keeps a derived BREP when ``solid`` had one and the
    chamfered tree still replays.
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
        raise SdfError("chamfer_solid: the solid is not SDF-authored")
    tree = from_construction(get_construction(solid))
    bevelled = chamfer(tree, distance)
    cell = _cell_size(solid, tree)
    lo, hi = bevelled.bounds
    longest = max(hi[i] - lo[i] for i in range(3))
    if cell and cell > 0.0 and longest > 0.0:
        resolution = math.ceil(longest / cell)
    else:
        resolution = DEFAULT_RESOLUTION
    resolution = int(min(max(resolution, 2), MAX_RESOLUTION))
    brep = "auto" if _derived_from(solid, tree) else False
    return to_solid(bevelled, resolution=resolution, brep=brep)
