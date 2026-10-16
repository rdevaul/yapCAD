"""SDF combinators and domain operators.

Hard CSG (:func:`union`, :func:`intersect`, :func:`subtract`) is the Phase 1
core; the smooth blends and :func:`transform` are here too because they are
what exercises the Lipschitz and exactness propagation in anger.

A note on ``exact`` for the boolean operators.  ``min(a, b)`` is *not* the
true signed distance to the union: inside the union it reports the distance
to the nearer operand's surface, which may lie in the interior of the other
operand and therefore not on the union's boundary at all.  It remains a valid
*bound*, though, and that is what the field machinery guarantees: ``min`` of
two 1-Lipschitz functions is 1-Lipschitz, and its zero set is exactly the
union's boundary, so ``|min(a,b)(p)| <= d(p)`` holds everywhere.  Hard CSG is
therefore ``exact=False, lipschitz=1``, which is precisely the distinction
the design document asks the schema to record.
"""

import math

import numpy as np

from yapcad.sdf.node import (
    EMPTY_BOUNDS,
    FieldProps,
    NodeSpec,
    SdfError,
    analyze,
    bounds_expand,
    bounds_intersect,
    bounds_is_empty,
    bounds_transform,
    bounds_union,
    make_node,
    register,
)


def _child_nodes(kind, children):
    """Validate and return a non-empty tuple of child nodes."""
    kids = tuple(children)
    if not kids:
        raise SdfError(f"{kind}: needs at least one operand")
    return kids


def _blend_radius(kind, k):
    """Validate a smooth-blend radius."""
    try:
        v = float(k)
    except (TypeError, ValueError):
        raise SdfError(f"{kind}: blend radius must be a number") from None
    if not math.isfinite(v) or v <= 0.0:
        raise SdfError(f"{kind}: blend radius must be positive and finite")
    return v


def _max_lipschitz(child_props):
    return max(c.lipschitz for c in child_props)


# ---------------------------------------------------------------------------
# Hard CSG
# ---------------------------------------------------------------------------


def _analyze_union(_params, child_props):
    bounds = EMPTY_BOUNDS
    for c in child_props:
        bounds = bounds_union(bounds, c.bounds)
    return FieldProps(
        exact=False,
        lipschitz=_max_lipschitz(child_props),
        bounds=bounds,
        csg_exact=all(c.csg_exact for c in child_props),
    )


#: Children above which a union culls by bounds rather than evaluating all of
#: them.  Below it the scan is cheaper than the bookkeeping: a union of a few
#: operands pays more to decide what to skip than to evaluate it.
_UNION_CULL_MIN_CHILDREN = 8

#: Points below which a union stops splitting the batch.  Every level re-tests
#: the children's bounds, so splitting too far costs more in bookkeeping than
#: it saves in skipped children: on an engraved line, simplification takes
#: 6.4 s at this value against 20 s unsplit, but 62 s at a threshold of 16.
_UNION_SPLIT_MIN_POINTS = 1024


def _batch_bounds(p):
    """Axis-aligned bounds of a point batch, or None if it is empty."""
    if p.size == 0:
        return None
    flat = p.reshape(-1, p.shape[-1])[:, :3]
    return flat.min(axis=0), flat.max(axis=0)


def _bounds_gap(batch, child_bounds):
    """Least possible distance from any point of ``batch`` to ``child_bounds``.

    Zero when the boxes touch or overlap, which is the case where nothing may
    be skipped: a point inside a child reads a negative value that no bound on
    distance constrains.
    """
    (blo, bhi), (clo, chi) = batch, child_bounds
    gap = 0.0
    for axis in range(3):
        low = clo[axis] - bhi[axis]
        high = blo[axis] - chi[axis]
        d = low if low > high else high
        if d > 0.0:
            gap += d * d
    return math.sqrt(gap)


def _eval_union(_params, p, children, ev):
    """``min`` of the children, skipping those that cannot lower it.

    Evaluating every child at every point is wasteful on unions with many
    operands -- engraved lettering reaches a couple of hundred -- when any one
    point is near only a few of them.  Two things make the saving available:
    a child whose bounds lie farther from the points than the smallest value
    found so far cannot win the ``min``, and a batch that is spread across the
    whole part can be split until its parts are local enough for that test to
    bite.  Callers do spread their points: the simplifier scores a round of
    deliberately non-adjacent collapses in one call.

    The result is unchanged bit for bit.  The children are 1-Lipschitz and
    bounded, so outside a child's bounds its field is at least the distance to
    them, and ``min`` is exact and order-independent in floating point.
    """
    if len(children) < _UNION_CULL_MIN_CHILDREN:
        out = ev(children[0], p)
        for child in children[1:]:
            out = np.minimum(out, ev(child, p))
        return out

    flat = p.reshape(-1, p.shape[-1])
    bounds = [analyze(child).bounds for child in children]
    out = np.empty(len(flat), dtype=float)
    _union_recurse(flat, np.arange(len(flat)), children, bounds, ev, out)
    return out.reshape(p.shape[:-1]) if p.ndim > 1 else out


def _union_recurse(flat, idx, children, bounds, ev, out):
    """Evaluate ``idx`` of ``flat``, culling and splitting as it pays to."""
    points = flat[idx]
    batch = _batch_bounds(points)

    keep, gaps = [], []
    if batch is not None:
        for child, cb in zip(children, bounds):
            gap = 0.0 if bounds_is_empty(cb) else _bounds_gap(batch, cb)
            keep.append(child)
            gaps.append(gap)
        order = np.argsort(np.asarray(gaps), kind="stable")
        keep = [keep[i] for i in order]
        gaps = [gaps[i] for i in order]
    else:
        keep, gaps = list(children), [0.0] * len(children)

    # Splitting only helps while the points are spread over many children and
    # there are enough of them to amortise the bookkeeping.
    if (len(idx) > _UNION_SPLIT_MIN_POINTS
            and len(keep) > _UNION_CULL_MIN_CHILDREN
            and batch is not None):
        lo, hi = batch
        axis = int(np.argmax(hi - lo))
        if hi[axis] - lo[axis] > 0.0:
            coord = flat[idx, axis]
            mid = np.median(coord)
            left = idx[coord <= mid]
            right = idx[coord > mid]
            if len(left) and len(right):
                _union_recurse(flat, left, children, bounds, ev, out)
                _union_recurse(flat, right, children, bounds, ev, out)
                return

    value = ev(keep[0], points)
    for gap, child in zip(gaps[1:], keep[1:]):
        if gap > 0.0 and gap >= value.max():
            break
        value = np.minimum(value, ev(child, points))
    out[idx] = value


register(NodeSpec(
    kind="union",
    min_children=1,
    max_children=None,
    analyze=_analyze_union,
    backends={"numpy": _eval_union},
))


def union(*children):
    """The union of every operand: ``min`` of the fields."""
    kids = _child_nodes("union", children)
    if len(kids) == 1:
        return kids[0]
    return make_node("union", {}, kids)


# ---------------------------------------------------------------------------
# Compound: separate bodies
# ---------------------------------------------------------------------------


register(NodeSpec(
    kind="compound",
    min_children=1,
    max_children=None,
    analyze=_analyze_union,
    backends={"numpy": _eval_union},
))


def compound(*children):
    """Several bodies kept separate, as the DSL's ``compound`` keeps them.

    As a field it *is* a union -- the region the bodies occupy -- so it
    combines in booleans exactly like one.  What it adds is the intent:
    :func:`yapcad.sdf.to_solid` meshes each body on its own, so touching or
    overlapping bodies stay distinct, closed meshes instead of being fused,
    and the OCC replay builds a compound rather than a fused solid.
    """
    kids = _child_nodes("compound", children)
    if len(kids) == 1:
        return kids[0]
    return make_node("compound", {}, kids)


def _analyze_intersect(_params, child_props):
    bounds = child_props[0].bounds
    for c in child_props[1:]:
        bounds = bounds_intersect(bounds, c.bounds)
    return FieldProps(
        exact=False,
        lipschitz=_max_lipschitz(child_props),
        bounds=bounds,
        csg_exact=all(c.csg_exact for c in child_props),
    )


def _eval_intersect(_params, p, children, ev):
    out = ev(children[0], p)
    for child in children[1:]:
        out = np.maximum(out, ev(child, p))
    return out


register(NodeSpec(
    kind="intersect",
    min_children=1,
    max_children=None,
    analyze=_analyze_intersect,
    backends={"numpy": _eval_intersect},
))


def intersect(*children):
    """The intersection of every operand: ``max`` of the fields."""
    kids = _child_nodes("intersect", children)
    if len(kids) == 1:
        return kids[0]
    return make_node("intersect", {}, kids)


def _analyze_subtract(_params, child_props):
    # Removing material can only shrink the result, so the first operand's
    # bound still holds.
    return FieldProps(
        exact=False,
        lipschitz=_max_lipschitz(child_props),
        bounds=child_props[0].bounds,
        csg_exact=all(c.csg_exact for c in child_props),
    )


def _eval_subtract(_params, p, children, ev):
    out = ev(children[0], p)
    for child in children[1:]:
        out = np.maximum(out, -ev(child, p))
    return out


register(NodeSpec(
    kind="subtract",
    min_children=2,
    max_children=None,
    analyze=_analyze_subtract,
    backends={"numpy": _eval_subtract},
))


def subtract(target, *tools):
    """Remove every ``tool`` from ``target``: ``max(target, -tool...)``."""
    if not tools:
        raise SdfError("subtract: needs at least one tool operand")
    return make_node("subtract", {}, (target,) + tuple(tools))


# ---------------------------------------------------------------------------
# Smooth blends
# ---------------------------------------------------------------------------
#
# All three are built on the polynomial smooth minimum
#
#     h = clamp(0.5 + 0.5*(b - a)/k, 0, 1)
#     smin = lerp(b, a, h) - k*h*(1 - h)
#
# whose partial derivatives work out to exactly `h` and `1 - h`.  The gradient
# is therefore a convex combination `h*grad(a) + (1-h)*grad(b)`, so
# `|grad| <= max(La, Lb)`: the polynomial blend preserves the Lipschitz bound
# even though it destroys exactness.  That is why these carry
# `lipschitz = max(children)` rather than something inflated.


def _smin(a, b, k):
    """Polynomial smooth minimum of two value arrays."""
    h = np.clip(0.5 + 0.5 * (b - a) / k, 0.0, 1.0)
    return b + (a - b) * h - k * h * (1.0 - h)


def _analyze_smooth_union(params, child_props):
    k = params["radius"]
    bounds = EMPTY_BOUNDS
    for c in child_props:
        bounds = bounds_union(bounds, c.bounds)
    return FieldProps(
        exact=False,
        lipschitz=_max_lipschitz(child_props),
        bounds=bounds_expand(bounds, k),
        csg_exact=False,
    )


def _eval_smooth_union(params, p, children, ev):
    return _smin(ev(children[0], p), ev(children[1], p), params["radius"])


register(NodeSpec(
    kind="smooth_union",
    min_children=2,
    max_children=2,
    analyze=_analyze_smooth_union,
    backends={"numpy": _eval_smooth_union},
))


def smooth_union(a, b, radius):
    """Union of ``a`` and ``b`` with a fillet of scale ``radius``."""
    return make_node(
        "smooth_union",
        {"radius": _blend_radius("smooth_union", radius)},
        (a, b),
    )


def _analyze_smooth_intersect(params, child_props):
    k = params["radius"]
    bounds = child_props[0].bounds
    for c in child_props[1:]:
        bounds = bounds_intersect(bounds, c.bounds)
    return FieldProps(
        exact=False,
        lipschitz=_max_lipschitz(child_props),
        bounds=bounds_expand(bounds, k),
        csg_exact=False,
    )


def _eval_smooth_intersect(params, p, children, ev):
    a = ev(children[0], p)
    b = ev(children[1], p)
    # smax(a, b) == -smin(-a, -b), which keeps one blend kernel authoritative.
    return -_smin(-a, -b, params["radius"])


register(NodeSpec(
    kind="smooth_intersect",
    min_children=2,
    max_children=2,
    analyze=_analyze_smooth_intersect,
    backends={"numpy": _eval_smooth_intersect},
))


def smooth_intersect(a, b, radius):
    """Intersection of ``a`` and ``b`` with a fillet of scale ``radius``."""
    return make_node(
        "smooth_intersect",
        {"radius": _blend_radius("smooth_intersect", radius)},
        (a, b),
    )


def _analyze_smooth_subtract(params, child_props):
    return FieldProps(
        exact=False,
        lipschitz=_max_lipschitz(child_props),
        bounds=bounds_expand(child_props[0].bounds, params["radius"]),
        csg_exact=False,
    )


def _eval_smooth_subtract(params, p, children, ev):
    a = ev(children[0], p)
    b = ev(children[1], p)
    return -_smin(-a, b, params["radius"])


register(NodeSpec(
    kind="smooth_subtract",
    min_children=2,
    max_children=2,
    analyze=_analyze_smooth_subtract,
    backends={"numpy": _eval_smooth_subtract},
))


def smooth_subtract(target, tool, radius):
    """Remove ``tool`` from ``target`` with a fillet of scale ``radius``."""
    return make_node(
        "smooth_subtract",
        {"radius": _blend_radius("smooth_subtract", radius)},
        (target, tool),
    )


# ---------------------------------------------------------------------------
# Offset and shell
# ---------------------------------------------------------------------------


def _analyze_offset(params, child_props):
    (c,) = child_props
    d = params["distance"]
    return FieldProps(
        exact=c.exact,
        lipschitz=c.lipschitz,
        bounds=bounds_expand(c.bounds, max(d, 0.0)),
        # OCC can offset a shape, but the design document restricts exact CSG
        # replay to primitives, hard booleans and similarity transforms.
        csg_exact=False,
    )


def _eval_offset(params, p, children, ev):
    return ev(children[0], p) - params["distance"]


register(NodeSpec(
    kind="offset",
    min_children=1,
    max_children=1,
    analyze=_analyze_offset,
    backends={"numpy": _eval_offset},
))


def offset(child, distance):
    """Grow ``child`` outward by ``distance`` (negative shrinks it).

    Exactness is inherited from ``child``, with the usual caveat that an
    inward offset stops being a true distance once it passes the solid's
    medial axis and the level set self-intersects.
    """
    try:
        d = float(distance)
    except (TypeError, ValueError):
        raise SdfError("offset: distance must be a number") from None
    if not math.isfinite(d):
        raise SdfError("offset: distance must be finite")
    return make_node("offset", {"distance": d}, (child,))


def _analyze_shell(params, child_props):
    (c,) = child_props
    return FieldProps(
        exact=c.exact,
        lipschitz=c.lipschitz,
        bounds=bounds_expand(c.bounds, 0.5 * params["thickness"]),
        csg_exact=False,
    )


def _eval_shell(params, p, children, ev):
    return np.abs(ev(children[0], p)) - 0.5 * params["thickness"]


register(NodeSpec(
    kind="shell",
    min_children=1,
    max_children=1,
    analyze=_analyze_shell,
    backends={"numpy": _eval_shell},
))


def shell(child, thickness):
    """A wall of ``thickness`` straddling ``child``'s surface."""
    t = float(thickness)
    if not math.isfinite(t) or t <= 0.0:
        raise SdfError("shell: thickness must be positive and finite")
    return make_node("shell", {"thickness": t}, (child,))


# ---------------------------------------------------------------------------
# Transform
# ---------------------------------------------------------------------------


def _as_matrix4(value):
    """Coerce a yapCAD ``xform.Matrix`` or a 4x4 sequence to nested tuples."""
    if hasattr(value, "get") and hasattr(value, "trans"):
        # An xform.Matrix: read through get(), which honours its transpose
        # flag.  Its raw ``.m`` does not, and would silently use M^T.
        rows = [[value.get(i, j) for j in range(4)] for i in range(4)]
    else:
        rows = getattr(value, "m", value)
    try:
        rows = [list(r) for r in rows]
    except TypeError:
        raise SdfError("transform: matrix must be a 4x4 sequence") from None
    if len(rows) != 4 or any(len(r) != 4 for r in rows):
        raise SdfError("transform: matrix must be 4x4")
    out = []
    for r in rows:
        vals = []
        for v in r:
            f = float(v)
            if not math.isfinite(f):
                raise SdfError("transform: matrix entries must be finite")
            vals.append(f)
        out.append(tuple(vals))
    bottom = out[3]
    if bottom != (0.0, 0.0, 0.0, 1.0):
        raise SdfError(
            "transform: matrix must be affine with bottom row [0,0,0,1], "
            f"got {list(bottom)}"
        )
    return tuple(out)


def _singular_range(matrix):
    """Return ``(sigma_min, sigma_max)`` of the 3x3 linear part."""
    linear = np.asarray([row[:3] for row in matrix[:3]], dtype=float)
    sv = np.linalg.svd(linear, compute_uv=False)
    return float(sv.min()), float(sv.max())


def _analyze_transform(params, child_props):
    (c,) = child_props
    matrix = params["matrix"]
    s_min, s_max = _singular_range(matrix)
    if s_min <= 0.0:
        raise SdfError("transform: matrix is singular and cannot be inverted")
    # The node evaluates g(p) = s_min * f(M^-1 p).  Then
    #   |grad g| = s_min * |(M^-1)^T grad f| <= s_min * (1/s_min) * L_f = L_f,
    # so a non-uniform scale costs accuracy (the field underestimates by up
    # to s_max/s_min) but never invalidates the Lipschitz bound.
    similarity = (s_max - s_min) <= 1e-12 * max(1.0, s_max)
    return FieldProps(
        exact=c.exact and similarity,
        lipschitz=c.lipschitz,
        bounds=bounds_transform(c.bounds, matrix),
        csg_exact=c.csg_exact and similarity,
    )


def _eval_transform(params, p, children, ev):
    matrix = np.asarray(params["matrix"], dtype=float)
    inverse = np.linalg.inv(matrix)
    local = p @ inverse[:3, :3].T + inverse[:3, 3]
    s_min, _ = _singular_range(params["matrix"])
    return ev(children[0], local) * s_min


register(NodeSpec(
    kind="transform",
    min_children=1,
    max_children=1,
    analyze=_analyze_transform,
    backends={"numpy": _eval_transform},
))


def transform(child, matrix):
    """Place ``child`` by the affine 4x4 ``matrix``.

    Accepts a :class:`yapcad.xform.Matrix` or any row-major 4x4 sequence.
    Rigid and uniform-scale transforms preserve exactness; a non-uniform
    scale does not, and is marked accordingly.
    """
    return make_node("transform", {"matrix": _as_matrix4(matrix)}, (child,))


def translate(child, delta):
    """Translate ``child`` by the 3-vector ``delta``."""
    seq = list(delta)
    if len(seq) < 3:
        raise SdfError("translate: delta needs at least 3 components")
    dx, dy, dz = (float(v) for v in seq[:3])
    return transform(child, (
        (1.0, 0.0, 0.0, dx),
        (0.0, 1.0, 0.0, dy),
        (0.0, 0.0, 1.0, dz),
        (0.0, 0.0, 0.0, 1.0),
    ))


def rotate(child, axis, angle):
    """Rotate ``child`` by ``angle`` degrees about ``axis``.

    Degrees, to match :func:`yapcad.xform.Rotation` and the rest of yapCAD.
    """
    seq = list(axis)
    if len(seq) < 3:
        raise SdfError("rotate: axis needs at least 3 components")
    ux, uy, uz = (float(v) for v in seq[:3])
    length = math.sqrt(ux * ux + uy * uy + uz * uz)
    if length == 0.0:
        raise SdfError("rotate: axis must be non-zero")
    ux, uy, uz = ux / length, uy / length, uz / length
    rad = math.radians(float(angle))
    c = math.cos(rad)
    s = math.sin(rad)
    t = 1.0 - c
    return transform(child, (
        (c + ux * ux * t, ux * uy * t - uz * s, ux * uz * t + uy * s, 0.0),
        (uy * ux * t + uz * s, c + uy * uy * t, uy * uz * t - ux * s, 0.0),
        (uz * ux * t - uy * s, uz * uy * t + ux * s, c + uz * uz * t, 0.0),
        (0.0, 0.0, 0.0, 1.0),
    ))


def scale(child, factor):
    """Scale ``child`` by a scalar or a 3-sequence of per-axis factors."""
    if isinstance(factor, (int, float)):
        sx = sy = sz = float(factor)
    else:
        seq = list(factor)
        if len(seq) < 3:
            raise SdfError("scale: factor needs at least 3 components")
        sx, sy, sz = (float(v) for v in seq[:3])
    if 0.0 in (sx, sy, sz):
        raise SdfError("scale: factors must be non-zero")
    return transform(child, (
        (sx, 0.0, 0.0, 0.0),
        (0.0, sy, 0.0, 0.0),
        (0.0, 0.0, sz, 0.0),
        (0.0, 0.0, 0.0, 1.0),
    ))
