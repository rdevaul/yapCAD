"""Replay CSG-exact SDF trees through OpenCASCADE for an exact BREP.

Design document §5.2: when a tree holds only analytic primitives, hard
booleans and similarity transforms, it is isomorphic to a CSG tree, and the
right move is not to convert the field at all but to *replay the same tree*
through OCC.  ``sphere`` becomes ``BRepPrimAPI_MakeSphere``, ``union`` becomes
``BRepAlgoAPI_Fuse``, and the result is an exact analytic BREP -- planes,
cylinders, spheres, cones, tori -- suitable for STEP export and for
downstream fillets.

Each node kind contributes an ``"occ"`` backend to its existing
:class:`~yapcad.sdf.node.NodeSpec` through
:func:`~yapcad.sdf.node.register_backend`, so this is the same seam the
numpy evaluator uses rather than a parallel registry.  Importing this module
requires pythonocc-core; importing :mod:`yapcad.sdf` does not.

Two parts of the replay are not a one-to-one call and are worth knowing:

``rounded_box``
    Built as a box with every edge filleted to the same radius, which is
    exactly the Minkowski rounding the field describes.  The classifier
    declines a radius equal to half the smallest edge, where the fillet has
    no face left to run along.

``half_space``
    OCC booleans want finite operands, so a half space is replayed as a box
    lying on the correct side of its plane, sized to cover the region where
    the result can exist.  That region is the root's bounds, carried into
    each node's local frame through any transforms above it.  Replacing a
    leaf with anything that agrees with it inside that region cannot change
    the result inside it, and the result lies inside it by construction;
    :func:`to_brep` also checks the replayed shape against the root's bounds
    after the fact, so a violated assumption is loud rather than silent.
"""

import math
from dataclasses import replace

import numpy as np

from yapcad.sdf.node import (
    INF,
    SdfError,
    analyze,
    bounds_transform,
    digest,
    get_spec,
    register_backend,
    walk,
)

try:  # pragma: no cover - exercised only where OCC is installed
    from OCC.Core.BRepAlgoAPI import (
        BRepAlgoAPI_Common,
        BRepAlgoAPI_Cut,
        BRepAlgoAPI_Fuse,
    )
    from OCC.Core.BRepBndLib import brepbndlib
    from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_Transform
    from OCC.Core.BRep import BRep_Builder
    from OCC.Core.BRepAdaptor import BRepAdaptor_Curve
    from OCC.Core.BRepCheck import BRepCheck_Analyzer
    from OCC.Core.BRepFilletAPI import BRepFilletAPI_MakeFillet
    from OCC.Core.BRepGProp import brepgprop
    from OCC.Core.BRepPrimAPI import (
        BRepPrimAPI_MakeBox,
        BRepPrimAPI_MakeCone,
        BRepPrimAPI_MakeCylinder,
        BRepPrimAPI_MakeSphere,
        BRepPrimAPI_MakeTorus,
    )
    from OCC.Core.Bnd import Bnd_Box
    from OCC.Core.GProp import GProp_GProps
    from OCC.Core.TopAbs import TopAbs_EDGE
    from OCC.Core.TopExp import TopExp_Explorer, topexp
    from OCC.Core.TopTools import TopTools_IndexedMapOfShape
    from OCC.Core.GeomAbs import GeomAbs_Circle
    from OCC.Core.TopoDS import TopoDS_Compound, topods
    from OCC.Core.gp import gp_Ax2, gp_Dir, gp_Pnt, gp_Trsf
    _OCC_ERROR = None
except ImportError as exc:  # pragma: no cover
    _OCC_ERROR = exc


def occ_available():
    """True when pythonocc-core can be imported."""
    return _OCC_ERROR is None


def _require_occ():
    if _OCC_ERROR is not None:
        raise RuntimeError(
            "replaying an SDF tree as BREP requires pythonocc-core "
            f"({_OCC_ERROR})"
        )


# ---------------------------------------------------------------------------
# Classification helpers (no OCC needed)
# ---------------------------------------------------------------------------


def _blocks_on_its_own(node):
    """True if ``node`` would be non-CSG even with every child CSG-exact.

    Asked of the node's own analysis with the children's ``csg_exact``
    forced true, so a ``smooth_union`` is reported even when one of its
    operands is also non-CSG -- each is a separate reason, and fixing only
    one would leave the part still unreplayable.
    """
    if node.csg_exact:
        return False
    spec = get_spec(node.kind)
    hypothetical = tuple(
        replace(analyze(child), csg_exact=True) for child in node.children
    )
    return not spec.analyze(node.p, hypothetical).csg_exact


def csg_blockers(node):
    """The nodes that stop ``node`` from being replayed as exact CSG.

    A blocker is a node that is non-CSG in its own right -- a smooth blend,
    a lattice, an offset, a non-similarity transform, a spindle torus --
    rather than an ancestor that merely inherits the loss.  Returned
    children-first, so the list reads in construction order.  Empty exactly
    when ``node.csg_exact``.

    This is what makes the §5.3 promise concrete: authoring can say "this
    part uses ``smooth_union`` at such-and-such node; it will not round-trip
    to STEP" before anyone attempts an export.
    """
    return [n for n in walk(node) if _blocks_on_its_own(n)]


def describe_blockers(node):
    """A human-readable account of :func:`csg_blockers`, or ``""``."""
    lines = []
    for blocker in csg_blockers(node):
        shown = ", ".join(f"{k}={v!r}" for k, v in blocker.params)
        lines.append(f"  - {blocker.kind}({shown})  [{blocker.digest}]")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Small OCC utilities
# ---------------------------------------------------------------------------


def _boolean(algorithm, a, b, kind):
    op = algorithm(a, b)
    if not op.IsDone() or op.HasErrors():
        raise SdfError(f"OCC {kind} failed while replaying the SDF tree")
    return op.Shape()


def _trsf(matrix):
    """A ``gp_Trsf`` from a row-major affine 4x4 similarity matrix."""
    t = gp_Trsf()
    t.SetValues(*[float(matrix[r][c]) for r in range(3) for c in range(4)])
    return t


def _transform(shape, matrix):
    op = BRepBuilderAPI_Transform(shape, _trsf(matrix), True)
    if not op.IsDone():
        raise SdfError("OCC transform failed while replaying the SDF tree")
    return op.Shape()


def _frame(normal):
    """Right-handed orthonormal ``(u, v, n)`` with ``n`` along ``normal``."""
    n = np.asarray(normal, dtype=float)
    n = n / np.linalg.norm(n)
    helper = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.9 else \
        np.array([0.0, 1.0, 0.0])
    u = np.cross(helper, n)
    u = u / np.linalg.norm(u)
    v = np.cross(n, u)
    return u, v, n


def volume(shape):
    """Enclosed volume of an OCC shape."""
    _require_occ()
    props = GProp_GProps()
    brepgprop.VolumeProperties(shape, props)
    return float(props.Mass())


def bounding_box(shape):
    """Tight ``((x0, y0, z0), (x1, y1, z1))`` of an OCC shape.

    ``AddOptimal`` rather than ``Add``: the fast path bounds a trimmed face
    by its whole underlying surface, so a sphere cut down by a box reports
    the sphere's full extent.
    """
    _require_occ()
    box = Bnd_Box()
    brepbndlib.AddOptimal(shape, box, False, False)
    x0, y0, z0, x1, y1, z1 = box.Get()
    return (x0, y0, z0), (x1, y1, z1)


# ---------------------------------------------------------------------------
# Per-kind replay backends
# ---------------------------------------------------------------------------
#
# Signature: (params, children, build, region) -> TopoDS_Shape.  ``build`` is
# ``build(child, region)``; ``region`` is the local-frame box the result must
# be correct within, used only by half_space.


def _occ_sphere(params, _children, _build, _region):
    return BRepPrimAPI_MakeSphere(gp_Pnt(0.0, 0.0, 0.0),
                                  params["radius"]).Shape()


def _occ_box(params, _children, _build, _region):
    hx, hy, hz = (v * 0.5 for v in params["size"])
    return BRepPrimAPI_MakeBox(gp_Pnt(-hx, -hy, -hz),
                               gp_Pnt(hx, hy, hz)).Shape()


def _occ_rounded_box(params, children, build, region):
    shape = _occ_box(params, children, build, region)
    radius = params["radius"]
    if radius == 0.0:
        return shape
    fillet = BRepFilletAPI_MakeFillet(shape)
    explorer = TopExp_Explorer(shape, TopAbs_EDGE)
    while explorer.More():
        # Every edge, and no silent skipping: a box with one edge left
        # sharp is a different part, not an approximation of this one.
        fillet.Add(radius, topods.Edge(explorer.Current()))
        explorer.Next()
    fillet.Build()
    if not fillet.IsDone():
        raise SdfError(
            f"OCC could not fillet a rounded_box of radius {radius}"
        )
    return fillet.Shape()


def _occ_cylinder(params, _children, _build, _region):
    h = params["height"]
    axis = gp_Ax2(gp_Pnt(0.0, 0.0, -0.5 * h), gp_Dir(0.0, 0.0, 1.0))
    return BRepPrimAPI_MakeCylinder(axis, params["radius"], h).Shape()


def _occ_rounded_cylinder(params, children, build, region):
    shape = _occ_cylinder(params, children, build, region)
    radius = params["edge_radius"]
    fillet = BRepFilletAPI_MakeFillet(shape)
    edges = TopTools_IndexedMapOfShape()
    topexp.MapShapes(shape, TopAbs_EDGE, edges)
    count = 0
    for k in range(1, edges.Size() + 1):
        edge = topods.Edge(edges.FindKey(k))
        # The two cap circles; the seam line is not an edge of the part.
        if BRepAdaptor_Curve(edge).GetType() == GeomAbs_Circle:
            fillet.Add(radius, edge)
            count += 1
    fillet.Build()
    if count != 2 or not fillet.IsDone():
        raise SdfError(
            f"OCC could not fillet a rounded_cylinder of edge radius {radius}"
        )
    return fillet.Shape()


def _occ_compound(_params, children, build, region):
    shape = TopoDS_Compound()
    builder = BRep_Builder()
    builder.MakeCompound(shape)
    for child in children:
        builder.Add(shape, build(child, region))
    return shape


def _occ_cone(params, _children, _build, _region):
    h = params["height"]
    r1 = params["radius_bottom"]
    r2 = params["radius_top"]
    axis = gp_Ax2(gp_Pnt(0.0, 0.0, -0.5 * h), gp_Dir(0.0, 0.0, 1.0))
    if r1 == r2:
        # OCC refuses a cone with two identical radii; it is a cylinder.
        return BRepPrimAPI_MakeCylinder(axis, r1, h).Shape()
    return BRepPrimAPI_MakeCone(axis, r1, r2, h).Shape()


def _occ_torus(params, _children, _build, _region):
    axis = gp_Ax2(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(0.0, 0.0, 1.0))
    return BRepPrimAPI_MakeTorus(axis, params["major_radius"],
                                 params["minor_radius"]).Shape()


def _occ_capsule(params, _children, _build, _region):
    a = np.asarray(params["start"], dtype=float)
    b = np.asarray(params["end"], dtype=float)
    r = params["radius"]
    length = float(np.linalg.norm(b - a))
    first = BRepPrimAPI_MakeSphere(gp_Pnt(*a), r).Shape()
    if length == 0.0:
        return first
    axis = gp_Ax2(gp_Pnt(*a), gp_Dir(*(b - a)))
    body = BRepPrimAPI_MakeCylinder(axis, r, length).Shape()
    body = _boolean(BRepAlgoAPI_Fuse, body, first, "capsule fuse")
    last = BRepPrimAPI_MakeSphere(gp_Pnt(*b), r).Shape()
    return _boolean(BRepAlgoAPI_Fuse, body, last, "capsule fuse")


def _occ_half_space(params, _children, _build, region):
    lo, hi = region
    if any(v in (INF, -INF) for v in (*lo, *hi)):
        raise SdfError(
            "a half_space can only be replayed inside a bounded tree; "
            "intersect it with a solid"
        )
    u, v, n = _frame(params["normal"])
    offset = params["offset"]
    corners = np.array([[x, y, z] for x in (lo[0], hi[0])
                        for y in (lo[1], hi[1]) for z in (lo[2], hi[2])])
    pu, pv, pn = corners @ u, corners @ v, corners @ n
    margin = 0.1 * float(np.linalg.norm(np.asarray(hi) - np.asarray(lo)))
    margin = max(margin, 1e-3)

    u0 = pu.min() - margin
    v0 = pv.min() - margin
    # The box ends exactly on the plane, and starts below the region. If the
    # plane misses the region entirely it becomes a slab outside it, which
    # is still correct: it agrees with the half space inside the region.
    n0 = min(pn.min(), offset) - margin
    origin = u * u0 + v * v0 + n * n0
    frame = gp_Ax2(gp_Pnt(*origin), gp_Dir(*n), gp_Dir(*u))
    return BRepPrimAPI_MakeBox(
        frame,
        float(pu.max() - pu.min() + 2.0 * margin),
        float(pv.max() - pv.min() + 2.0 * margin),
        float(offset - n0),
    ).Shape()


def _occ_union(_params, children, build, region):
    shape = build(children[0], region)
    for child in children[1:]:
        shape = _boolean(BRepAlgoAPI_Fuse, shape, build(child, region),
                         "union")
    return shape


def _occ_intersect(_params, children, build, region):
    shape = build(children[0], region)
    for child in children[1:]:
        shape = _boolean(BRepAlgoAPI_Common, shape, build(child, region),
                         "intersect")
    return shape


def _occ_subtract(_params, children, build, region):
    shape = build(children[0], region)
    for child in children[1:]:
        shape = _boolean(BRepAlgoAPI_Cut, shape, build(child, region),
                         "subtract")
    return shape


def _occ_transform(params, children, build, region):
    matrix = params["matrix"]
    inverse = np.linalg.inv(np.asarray(matrix, dtype=float))
    local = bounds_transform(region, tuple(tuple(r) for r in inverse))
    return _transform(build(children[0], local), matrix)


_BACKENDS = {
    "sphere": _occ_sphere,
    "box": _occ_box,
    "rounded_box": _occ_rounded_box,
    "rounded_cylinder": _occ_rounded_cylinder,
    "compound": _occ_compound,
    "cylinder": _occ_cylinder,
    "cone": _occ_cone,
    "torus": _occ_torus,
    "capsule": _occ_capsule,
    "half_space": _occ_half_space,
    "union": _occ_union,
    "intersect": _occ_intersect,
    "subtract": _occ_subtract,
    "transform": _occ_transform,
}

for _kind, _fn in _BACKENDS.items():
    register_backend(_kind, "occ", _fn)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def to_occ_shape(node):
    """Replay ``node`` and return the raw ``TopoDS_Shape``.

    Raises :class:`~yapcad.sdf.node.SdfError` for a tree that is not
    ``csg_exact`` -- naming the nodes responsible -- or that is unbounded,
    and when the replayed shape is invalid or escapes the tree's bounds.
    """
    _require_occ()
    if not node.csg_exact:
        raise SdfError(
            "this SDF tree cannot be replayed as exact CSG, because of:\n"
            + describe_blockers(node)
            + "\nMesh it with to_solid() instead; STEP export will then be "
            "faceted rather than analytic."
        )
    region = node.bounds
    lo, hi = region
    if any(v in (INF, -INF) for v in (*lo, *hi)):
        raise SdfError(
            "an unbounded SDF tree has no finite BREP; intersect it with a "
            "solid first"
        )

    cache = {}

    def build(child, child_region):
        key = (digest(child), child_region)
        if key not in cache:
            spec_fn = get_spec(child.kind).backends.get("occ")
            if spec_fn is None:
                raise SdfError(f"no OCC replay for node kind {child.kind!r}")
            cache[key] = spec_fn(child.p, child.children, build,
                                 child_region)
        return cache[key]

    shape = build(node, region)

    if not BRepCheck_Analyzer(shape).IsValid():
        raise SdfError("OCC produced an invalid shape replaying the SDF tree")
    _check_containment(shape, region)
    return shape


def _check_containment(shape, region):
    """The replayed shape must lie within the tree's own bounds.

    The field's bounds are conservative, so a correct replay always fits.
    A shape that does not fit means a half-space stand-in leaked past the
    region it was sized for -- the one assumption the replay rests on.
    """
    lo, hi = region
    diag = math.dist(lo, hi)
    tolerance = 1e-6 * diag + 1e-6
    (x0, y0, z0), (x1, y1, z1) = bounding_box(shape)
    # Bnd_Box enlarges by the shape tolerance; allow for it generously.
    slack = tolerance + 1e-3 * diag
    if (x0 < lo[0] - slack or y0 < lo[1] - slack or z0 < lo[2] - slack
            or x1 > hi[0] + slack or y1 > hi[1] + slack
            or z1 > hi[2] + slack):
        raise SdfError(
            "the replayed BREP extends beyond the SDF tree's bounds; this is "
            "a replay bug, please report it with the tree"
        )


def to_brep(node):
    """Replay ``node`` as a :class:`yapcad.brep.BrepSolid`."""
    from yapcad.brep import BrepSolid
    return BrepSolid(to_occ_shape(node))


__all__ = [
    "occ_available",
    "csg_blockers",
    "describe_blockers",
    "to_occ_shape",
    "to_brep",
    "volume",
    "bounding_box",
]
