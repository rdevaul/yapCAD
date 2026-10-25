"""Booleans between SDF-authored solids, done on the fields.

Design document §9 calls for booleans to be polymorphic: when both operands
were authored as fields, they should stay fields.  The combination is then
arithmetic on scalars -- ``min``, ``max`` -- which cannot fail, has no
coplanar-face or sliver cases, and keeps the result SDF-authoritative so it
can still be re-meshed, replayed through OCC, or combined again.

This matters in practice, not just in principle.  Handed two dual-contoured
meshes, the native mesh boolean engine returns an open, unusable mesh even
for two boxes sharing a face, takes tens of seconds at modest resolution, and
at higher resolution times out or raises on a degenerate face.  Every yapCAD
boolean goes through :func:`yapcad.geom3d.solid_boolean`, the DSL's included,
so routing SDF pairs here fixes all of them at once.

The one real decision is the resolution the result is meshed at.  The rule
is: *never coarser than either operand*.  Each operand's cell size is
recovered from the meshing parameters it records, or estimated from its mesh
when it has none (a transformed solid drops them), and the smaller cell is
carried over the result's extent.  :data:`MAX_RESOLUTION`, which comes from
:mod:`yapcad.sdf.convert`, caps it, because unioning a finely meshed small
part with a large one would otherwise ask for an unbounded lattice.

That rule reads the operands, and a boolean can make a feature thinner than
anything either operand had -- two surfaces brought nearly together, or a
cut that leaves a sliver.  Nothing here can know that in advance, so the
cell is a starting point: :func:`yapcad.sdf.to_solid` refines it when the
surface does not show up at it, and says so.  ``cell_mm`` is the way to
skip that extra pass once the right cell is known.

Mixed operands -- one field, one BREP or mesh -- are not handled here yet;
they still go to the existing engines.  That is the rest of Phase 4.
"""

import math

import numpy as np

from yapcad.construction import get_construction
from yapcad.sdf.convert import MAX_RESOLUTION, to_solid
from yapcad.sdf.node import (
    INF,
    SdfError,
    from_construction,
    is_sdf_construction,
    meshing_from_construction,
    to_construction,
)
from yapcad.sdf.ops import intersect, subtract, union

__all__ = [
    "DEFAULT_RESOLUTION",
    "MAX_RESOLUTION",
    "OPERATIONS",
    "combine_all",
    "combine_solids",
    "is_sdf_solid",
]

#: yapCAD's boolean vocabulary, mapped onto field operations.
OPERATIONS = {
    "union": union,
    "intersection": intersect,
    "difference": subtract,
}

#: Used only when neither operand's cell size can be recovered.
DEFAULT_RESOLUTION = 64


def is_sdf_solid(solid):
    """True if ``solid`` was authored as a field (carries an SDF tree)."""
    try:
        record = get_construction(solid)
    except ValueError:
        return False
    return is_sdf_construction(record)


def _tree(solid, role):
    if not is_sdf_solid(solid):
        raise SdfError(f"the {role} operand is not an SDF-authored solid")
    return from_construction(get_construction(solid))


def _finite(bounds):
    lo, hi = bounds
    return all(v not in (INF, -INF) for v in (*lo, *hi))


def _longest(bounds):
    lo, hi = bounds
    return max(hi[i] - lo[i] for i in range(3))


def _cell_size(solid, tree):
    """The cell size ``solid`` was meshed at, or ``None`` if unknowable.

    Exact when the solid records its meshing parameters, since dual
    contouring's spacing is the region's longest extent over the resolution.
    Otherwise estimated as the median mesh edge length, which for a
    dual-contoured mesh is close to one cell.
    """
    meshing = meshing_from_construction(get_construction(solid))
    if meshing and meshing.get("resolution"):
        if meshing.get("bounds"):
            flat = meshing["bounds"]
            region = (tuple(flat[:3]), tuple(flat[3:]))
        else:
            region = tree.bounds
        if _finite(region) and _longest(region) > 0.0:
            return _longest(region) / float(meshing["resolution"])

    edges = []
    for surf in solid[1]:
        verts = np.asarray([v[:3] for v in surf[1]], dtype=float)
        faces = np.asarray(surf[3], dtype=np.int64)
        if len(faces) == 0:
            continue
        for i, j in ((0, 1), (1, 2), (2, 0)):
            delta = verts[faces[:, i]] - verts[faces[:, j]]
            edges.append(np.linalg.norm(delta, axis=1))
    if not edges:
        return None
    lengths = np.concatenate(edges)
    lengths = lengths[lengths > 0.0]
    return float(np.median(lengths)) if lengths.size else None


def _mesh_bounds(solid):
    """Axis-aligned bounds of a solid's mesh vertices."""
    pts = np.asarray([v[:3] for surf in solid[1] for v in surf[1]],
                     dtype=float)
    return tuple(pts.min(axis=0)), tuple(pts.max(axis=0))


def _region(solids, node, operation):
    """Where the result can lie: the tree's bounds when finite.

    A tree can be unbounded even though its solid is finite -- a field that
    was meshed with explicit bounds.  Then the operands' own meshes bound the
    result, combined the same way the operation combines the solids.
    """
    if _finite(node.bounds):
        return node.bounds, False
    boxes = [_mesh_bounds(s) for s in solids]
    if operation == "union":
        lo = tuple(min(b[0][i] for b in boxes) for i in range(3))
        hi = tuple(max(b[1][i] for b in boxes) for i in range(3))
    elif operation == "intersection":
        lo = tuple(max(b[0][i] for b in boxes) for i in range(3))
        hi = tuple(min(b[1][i] for b in boxes) for i in range(3))
    else:
        lo, hi = boxes[0]
    return (lo, hi), True


def _derived_from(solid, tree):
    """True if ``solid`` carries a BREP replayed from ``tree`` itself."""
    from yapcad.metadata import get_solid_metadata
    meta = get_solid_metadata(solid) or {}
    marker = ((meta.get("brep") or {}).get("derivedFrom") or {}).get("sdf")
    return marker == tree.digest


def combine_all(solids, operation, cell_mm=None):
    """Combine any number of SDF-authored solids in one field operation.

    :param solids: two or more SDF-authored solids.
    :param operation: ``"union"``, ``"intersection"`` or ``"difference"``
        (the first solid minus all the rest), matching
        :func:`yapcad.geom3d.solid_boolean`.
    :param cell_mm: mesh the result at this cell size instead of the
        operand-derived one -- what a previous run's refinement warning
        reported, when the result is known to be finer than its operands.
    :returns: a new SDF-authoritative solid whose tree is the combination.

    Fields combine n-ary for free, so this meshes once, where folding
    :func:`combine_solids` pairwise would mesh -- and discard -- every
    intermediate result.  If every operand carries a derived BREP, the
    result does too whenever the combined tree is still replayable, so exact
    geometry stays exact.
    """
    if operation not in OPERATIONS:
        raise ValueError(f"unsupported boolean operation {operation!r}")
    solids = list(solids)
    if len(solids) < 2:
        raise SdfError("a boolean needs at least two operands")
    roles = ["first", "second"] + [f"operand {i + 1}"
                                   for i in range(2, len(solids))]
    trees = [_tree(s, role) for s, role in zip(solids, roles)]
    node = OPERATIONS[operation](*trees)

    region, explicit = _region(solids, node, operation)
    lo, hi = region
    if any(lo[i] > hi[i] for i in range(3)):
        # Disjoint bounds for an intersection: the result is empty.  yapCAD
        # allows an empty solid precisely for empty CSG results; it still
        # carries the tree, which remains the truth about the part.
        from yapcad.geom3d import solid
        return solid([], [], to_construction(node))

    exact = all(_derived_from(s, t) for s, t in zip(solids, trees))
    brep = "auto" if exact else False
    bounds = region if explicit else None
    if cell_mm is not None:
        # An explicit cell is the caller's call, and is not capped: the cap
        # is there to stop a *derived* resolution running away.
        return to_solid(node, cell_mm=cell_mm, bounds=bounds, brep=brep)

    cells = [c for c in (_cell_size(s, t) for s, t in zip(solids, trees))
             if c is not None and c > 0.0]
    longest = _longest(region)
    if cells and longest > 0.0:
        resolution = math.ceil(longest / min(cells))
    else:
        resolution = DEFAULT_RESOLUTION
    resolution = int(min(max(resolution, 2), MAX_RESOLUTION))
    return to_solid(node, resolution=resolution, bounds=bounds, brep=brep)


def combine_solids(a, b, operation, cell_mm=None):
    """Combine two SDF-authored solids; see :func:`combine_all`."""
    return combine_all([a, b], operation, cell_mm=cell_mm)
