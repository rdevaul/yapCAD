"""Analytic SDF primitives.

Every primitive here is defined in a canonical pose — centred on the origin
and aligned to the axes — and is placed with :func:`yapcad.sdf.ops.transform`
or its ``translate``/``rotate``/``scale`` helpers.  Baking a position into
each primitive would give two ways to express the same solid, and therefore
two different content digests for it, which defeats the DAG deduplication
that :func:`yapcad.sdf.node.tree_to_json` depends on.

The exceptions are :func:`capsule`, which is naturally defined by its two end
points, and :func:`half_space`, which is defined by its plane.

All of these are exact distance fields with a Lipschitz constant of 1, except
:func:`gyroid` and :func:`schwarz_p`, whose constants are derived at their
definitions.  Each evaluator takes points of shape ``(N, 3)`` and returns
values of shape ``(N,)``, per the vectorisation requirement in design
document §3.1.
"""

import math

import numpy as np

from yapcad.sdf.node import (
    UNBOUNDED,
    FieldProps,
    NodeSpec,
    SdfError,
    make_node,
    register,
)


def _positive(kind, name, value):
    """Validate that a dimension is a positive, finite number."""
    try:
        v = float(value)
    except (TypeError, ValueError):
        raise SdfError(
            f"{kind}: {name} must be a number, got {value!r}"
        ) from None
    if not math.isfinite(v) or v <= 0.0:
        raise SdfError(f"{kind}: {name} must be positive and finite, got {v}")
    return v


def _nonnegative(kind, name, value):
    """Validate that a dimension is a non-negative, finite number."""
    try:
        v = float(value)
    except (TypeError, ValueError):
        raise SdfError(
            f"{kind}: {name} must be a number, got {value!r}"
        ) from None
    if not math.isfinite(v) or v < 0.0:
        raise SdfError(f"{kind}: {name} must be non-negative, got {v}")
    return v


def _point3(kind, name, value):
    """Coerce a yapCAD point or any 3-sequence to a plain 3-tuple of floats."""
    try:
        seq = list(value)
    except TypeError:
        raise SdfError(f"{kind}: {name} must be a sequence") from None
    if len(seq) < 3:
        raise SdfError(f"{kind}: {name} needs at least 3 components")
    # yapCAD points carry a homogeneous w at index 3; ignore it.
    out = []
    for v in seq[:3]:
        f = float(v)
        if not math.isfinite(f):
            raise SdfError(f"{kind}: {name} components must be finite")
        out.append(f)
    return tuple(out)


def _exact(bounds):
    """FieldProps for an exact analytic primitive."""
    return FieldProps(
        exact=True, lipschitz=1.0, bounds=bounds, csg_exact=True
    )


# ---------------------------------------------------------------------------
# sphere
# ---------------------------------------------------------------------------


def _analyze_sphere(params, _children):
    r = params["radius"]
    return _exact(((-r, -r, -r), (r, r, r)))


def _eval_sphere(params, p, _children, _ev):
    return np.linalg.norm(p, axis=1) - params["radius"]


register(NodeSpec(
    kind="sphere",
    min_children=0,
    max_children=0,
    analyze=_analyze_sphere,
    backends={"numpy": _eval_sphere},
))


def sphere(radius):
    """A sphere of ``radius`` centred on the origin."""
    r = _positive("sphere", "radius", radius)
    return make_node("sphere", {"radius": r})


# ---------------------------------------------------------------------------
# box and rounded box
# ---------------------------------------------------------------------------


def _size3(kind, size):
    """Accept a scalar cube edge or a 3-sequence of edge lengths."""
    if isinstance(size, (int, float)):
        s = _positive(kind, "size", size)
        return (s, s, s)
    vals = _point3(kind, "size", size)
    for v in vals:
        if v <= 0.0:
            raise SdfError(f"{kind}: every size component must be positive")
    return vals


def _analyze_box(params, _children):
    hx, hy, hz = (v * 0.5 for v in params["size"])
    return _exact(((-hx, -hy, -hz), (hx, hy, hz)))


def _box_distance(p, half):
    """Exact distance to an origin-centred box of the given half-extents."""
    q = np.abs(p) - np.asarray(half, dtype=float)
    outside = np.linalg.norm(np.maximum(q, 0.0), axis=1)
    inside = np.minimum(np.max(q, axis=1), 0.0)
    return outside + inside


def _eval_box(params, p, _children, _ev):
    return _box_distance(p, [v * 0.5 for v in params["size"]])


register(NodeSpec(
    kind="box",
    min_children=0,
    max_children=0,
    analyze=_analyze_box,
    backends={"numpy": _eval_box},
))


def box(size):
    """An axis-aligned box centred on the origin.

    :param size: a scalar edge length, or a 3-sequence of edge lengths.
    """
    return make_node("box", {"size": _size3("box", size)})


def _analyze_rounded_box(params, _children):
    hx, hy, hz = (v * 0.5 for v in params["size"])
    props = _exact(((-hx, -hy, -hz), (hx, hy, hz)))
    # Replayed as a box with every edge filleted, which is exact -- except
    # when the radius consumes the whole of the smallest edge and the
    # fillet has no face left to run along.  OCC refuses that, so the
    # classifier must not promise it.
    if params["radius"] >= 0.5 * min(params["size"]):
        props = FieldProps(exact=props.exact, lipschitz=props.lipschitz,
                           bounds=props.bounds, csg_exact=False)
    return props


def _eval_rounded_box(params, p, _children, _ev):
    r = params["radius"]
    half = [v * 0.5 - r for v in params["size"]]
    return _box_distance(p, half) - r


register(NodeSpec(
    kind="rounded_box",
    min_children=0,
    max_children=0,
    analyze=_analyze_rounded_box,
    backends={"numpy": _eval_rounded_box},
))


def rounded_box(size, radius):
    """A box with edges and corners rounded to ``radius``.

    The overall extent is ``size``: the rounding is taken out of the box
    rather than added to it, so a rounded box and a box of the same size
    occupy the same bounding volume.
    """
    dims = _size3("rounded_box", size)
    r = _nonnegative("rounded_box", "radius", radius)
    if r > 0.5 * min(dims):
        raise SdfError(
            f"rounded_box: radius {r} exceeds half the smallest edge "
            f"{0.5 * min(dims)}"
        )
    return make_node("rounded_box", {"size": dims, "radius": r})


# ---------------------------------------------------------------------------
# cylinder
# ---------------------------------------------------------------------------


def _analyze_cylinder(params, _children):
    r = params["radius"]
    h = params["height"] * 0.5
    return _exact(((-r, -r, -h), (r, r, h)))


def _eval_cylinder(params, p, _children, _ev):
    d_radial = np.hypot(p[:, 0], p[:, 1]) - params["radius"]
    d_axial = np.abs(p[:, 2]) - 0.5 * params["height"]
    outside = np.hypot(
        np.maximum(d_radial, 0.0), np.maximum(d_axial, 0.0)
    )
    inside = np.minimum(np.maximum(d_radial, d_axial), 0.0)
    return outside + inside


register(NodeSpec(
    kind="cylinder",
    min_children=0,
    max_children=0,
    analyze=_analyze_cylinder,
    backends={"numpy": _eval_cylinder},
))


def cylinder(radius, height):
    """A capped cylinder on the z axis, centred on the origin."""
    return make_node("cylinder", {
        "radius": _positive("cylinder", "radius", radius),
        "height": _positive("cylinder", "height", height),
    })


# ---------------------------------------------------------------------------
# capsule
# ---------------------------------------------------------------------------


def _analyze_capsule(params, _children):
    a = params["start"]
    b = params["end"]
    r = params["radius"]
    lo = tuple(min(a[i], b[i]) - r for i in range(3))
    hi = tuple(max(a[i], b[i]) + r for i in range(3))
    return _exact((lo, hi))


def _eval_capsule(params, p, _children, _ev):
    a = np.asarray(params["start"], dtype=float)
    b = np.asarray(params["end"], dtype=float)
    r = params["radius"]
    ba = b - a
    denom = float(np.dot(ba, ba))
    pa = p - a
    if denom == 0.0:
        # Degenerate capsule: the two end points coincide, so it is a sphere.
        return np.linalg.norm(pa, axis=1) - r
    t = np.clip((pa @ ba) / denom, 0.0, 1.0)
    return np.linalg.norm(pa - t[:, None] * ba, axis=1) - r


register(NodeSpec(
    kind="capsule",
    min_children=0,
    max_children=0,
    analyze=_analyze_capsule,
    backends={"numpy": _eval_capsule},
))


def capsule(start, end, radius):
    """A capsule: the set of points within ``radius`` of segment start-end."""
    return make_node("capsule", {
        "start": _point3("capsule", "start", start),
        "end": _point3("capsule", "end", end),
        "radius": _positive("capsule", "radius", radius),
    })


# ---------------------------------------------------------------------------
# torus
# ---------------------------------------------------------------------------


def _analyze_torus(params, _children):
    outer = params["major_radius"] + params["minor_radius"]
    r = params["minor_radius"]
    props = _exact(((-outer, -outer, -r), (outer, outer, r)))
    # A spindle or horn torus (minor >= major) is a fine field -- the union
    # of the swept tube -- but its BREP surface self-intersects: OCC builds
    # it, reports it valid, and counts the overlap twice in the volume.  So
    # it is not replayable, whatever the field says.
    if params["minor_radius"] >= params["major_radius"]:
        props = FieldProps(exact=props.exact, lipschitz=props.lipschitz,
                           bounds=props.bounds, csg_exact=False)
    return props


def _eval_torus(params, p, _children, _ev):
    radial = np.hypot(p[:, 0], p[:, 1]) - params["major_radius"]
    return np.hypot(radial, p[:, 2]) - params["minor_radius"]


register(NodeSpec(
    kind="torus",
    min_children=0,
    max_children=0,
    analyze=_analyze_torus,
    backends={"numpy": _eval_torus},
))


def torus(major_radius, minor_radius):
    """A torus in the xy plane, centred on the origin.

    :param major_radius: distance from the origin to the tube centre.
    :param minor_radius: tube radius.
    """
    return make_node("torus", {
        "major_radius": _positive("torus", "major_radius", major_radius),
        "minor_radius": _positive("torus", "minor_radius", minor_radius),
    })


# ---------------------------------------------------------------------------
# cone (capped, i.e. a frustum)
# ---------------------------------------------------------------------------


def _analyze_cone(params, _children):
    r = max(params["radius_bottom"], params["radius_top"])
    h = params["height"] * 0.5
    return _exact(((-r, -r, -h), (r, r, h)))


def _eval_cone(params, p, _children, _ev):
    # Exact capped-cone field, reduced to the 2D (radial, axial) half-plane:
    # the distance is to whichever of the two caps or the slanted side is
    # nearest, with the sign taken from being inside both.
    h = 0.5 * params["height"]
    r1 = params["radius_bottom"]
    r2 = params["radius_top"]
    qx = np.hypot(p[:, 0], p[:, 1])
    qy = p[:, 2]

    # ca: offset to the nearest cap disc.
    cap_radius = np.where(qy < 0.0, r1, r2)
    ca_x = qx - np.minimum(qx, cap_radius)
    ca_y = np.abs(qy) - h

    # cb: offset to the slanted side, as a clamped projection onto it.
    k1x, k1y = r2, h
    k2x, k2y = r2 - r1, 2.0 * h
    k2_sq = k2x * k2x + k2y * k2y
    t = np.clip(((k1x - qx) * k2x + (k1y - qy) * k2y) / k2_sq, 0.0, 1.0)
    cb_x = qx - k1x + k2x * t
    cb_y = qy - k1y + k2y * t

    inside = (cb_x < 0.0) & (ca_y < 0.0)
    sign = np.where(inside, -1.0, 1.0)
    nearest_sq = np.minimum(
        ca_x * ca_x + ca_y * ca_y, cb_x * cb_x + cb_y * cb_y
    )
    return sign * np.sqrt(nearest_sq)


register(NodeSpec(
    kind="cone",
    min_children=0,
    max_children=0,
    analyze=_analyze_cone,
    backends={"numpy": _eval_cone},
))


def cone(radius_bottom, radius_top, height):
    """A capped cone (frustum) on the z axis, centred on the origin.

    ``radius_bottom`` applies at ``z = -height/2`` and ``radius_top`` at
    ``z = +height/2``.  A true point-tipped cone is the ``radius_top = 0``
    case.
    """
    rb = _nonnegative("cone", "radius_bottom", radius_bottom)
    rt = _nonnegative("cone", "radius_top", radius_top)
    if rb == 0.0 and rt == 0.0:
        raise SdfError("cone: at least one radius must be non-zero")
    return make_node("cone", {
        "radius_bottom": rb,
        "radius_top": rt,
        "height": _positive("cone", "height", height),
    })


# ---------------------------------------------------------------------------
# half space
# ---------------------------------------------------------------------------


def _analyze_half_space(_params, _children):
    return FieldProps(
        exact=True, lipschitz=1.0, bounds=UNBOUNDED, csg_exact=True
    )


def _eval_half_space(params, p, _children, _ev):
    n = np.asarray(params["normal"], dtype=float)
    return p @ n - params["offset"]


register(NodeSpec(
    kind="half_space",
    min_children=0,
    max_children=0,
    analyze=_analyze_half_space,
    backends={"numpy": _eval_half_space},
))


def half_space(normal, offset=0.0):
    """The half space on the negative-``normal`` side of a plane.

    The solid is ``dot(p, normal) <= offset``; the normal points *out* of the
    material, matching the sign convention of the rest of the module.
    ``normal`` is normalised, which is what keeps the field exact.
    """
    n = _point3("half_space", "normal", normal)
    length = math.sqrt(sum(v * v for v in n))
    if length == 0.0:
        raise SdfError("half_space: normal must be non-zero")
    unit = tuple(v / length for v in n)
    try:
        d = float(offset)
    except (TypeError, ValueError):
        raise SdfError("half_space: offset must be a number") from None
    if not math.isfinite(d):
        raise SdfError("half_space: offset must be finite")
    return make_node("half_space", {"normal": unit, "offset": d})


# ---------------------------------------------------------------------------
# Triply periodic minimal surfaces
# ---------------------------------------------------------------------------
#
# These are the first nodes whose field is not a distance, and they are
# included in Phase 1 deliberately: they are the test case that proves the
# Lipschitz machinery is doing real work rather than propagating 1.0
# everywhere.  Their `thickness` is in *field* units, not millimetres; divide
# by the node's Lipschitz constant for an approximate wall thickness.


def _analyze_gyroid(params, _children):
    w = 2.0 * math.pi / params["period"]
    # f = sin(wx)cos(wy) + sin(wy)cos(wz) + sin(wz)cos(wx)
    # df/dx = w[cos(wx)cos(wy) - sin(wz)sin(wx)], so |df/dx| <= 2w, and the
    # same bound holds on each axis; hence |grad f| <= 2*sqrt(3)*w.
    return FieldProps(
        exact=False,
        lipschitz=2.0 * math.sqrt(3.0) * w,
        bounds=UNBOUNDED,
        csg_exact=False,
    )


def _eval_gyroid(params, p, _children, _ev):
    w = 2.0 * math.pi / params["period"]
    x, y, z = p[:, 0] * w, p[:, 1] * w, p[:, 2] * w
    f = np.sin(x) * np.cos(y) + np.sin(y) * np.cos(z) + np.sin(z) * np.cos(x)
    return np.abs(f) - params["thickness"]


register(NodeSpec(
    kind="gyroid",
    min_children=0,
    max_children=0,
    analyze=_analyze_gyroid,
    backends={"numpy": _eval_gyroid},
))


def gyroid(period, thickness):
    """A gyroid sheet of the given period, infinite in every direction.

    Intersect it with a bounding solid to get a finite lattice.  The field is
    not metric: see the note on ``thickness`` above.
    """
    return make_node("gyroid", {
        "period": _positive("gyroid", "period", period),
        "thickness": _positive("gyroid", "thickness", thickness),
    })


def _analyze_schwarz_p(params, _children):
    w = 2.0 * math.pi / params["period"]
    # f = cos(wx) + cos(wy) + cos(wz); df/dx = -w sin(wx), so |df/dx| <= w
    # on each axis and |grad f| <= sqrt(3)*w.
    return FieldProps(
        exact=False,
        lipschitz=math.sqrt(3.0) * w,
        bounds=UNBOUNDED,
        csg_exact=False,
    )


def _eval_schwarz_p(params, p, _children, _ev):
    w = 2.0 * math.pi / params["period"]
    f = np.cos(p[:, 0] * w) + np.cos(p[:, 1] * w) + np.cos(p[:, 2] * w)
    return np.abs(f) - params["thickness"]


register(NodeSpec(
    kind="schwarz_p",
    min_children=0,
    max_children=0,
    analyze=_analyze_schwarz_p,
    backends={"numpy": _eval_schwarz_p},
))


def schwarz_p(period, thickness):
    """A Schwarz P sheet of the given period, infinite in every direction."""
    return make_node("schwarz_p", {
        "period": _positive("schwarz_p", "period", period),
        "thickness": _positive("schwarz_p", "thickness", thickness),
    })
