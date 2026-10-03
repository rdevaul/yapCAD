"""Planar profiles and the solids swept from them.

Gears and fasteners are, at heart, a 2D outline carried through space: a
spur gear's tooth profile extruded, a bevel gear's profile scaled toward an
apex, a hex nut's hexagon extruded and chamfered.  These nodes give that
pattern exact fields:

``polygon``
    The exact signed distance to a closed 2D polygon, constant along z.
    An ``n``-fold rotationally symmetric polygon -- a gear -- evaluates only
    the segments near one sector, after rotating each point into it, which
    is exact (each point's nearest segment lies in its own sector or a
    neighbour) and ``n/3`` times cheaper than testing every segment.
``extrude``
    A planar child extruded symmetrically about z = 0.  Exact when the
    child is.
``apex_extrude``
    The cone through the origin over a planar profile given in the plane
    ``z = z_ref``, clipped to ``z_lo <= z <= z_hi``: every section is the
    profile scaled by ``z / z_ref``.  This is exactly how a straight bevel
    gear's teeth are built (see :mod:`yapcad.gears.bevel`).

None of them replays through OCC directly; semantic nodes built on them,
like ``straight_bevel_gear``, replay through their own generators.
"""

import math

import numpy as np

from yapcad.sdf.node import INF, FieldProps, NodeSpec, SdfError, make_node, \
    register

__all__ = ["polygon", "extrude", "apex_extrude"]

#: Points evaluated per batch against all segments, bounding memory at
#: about ``BATCH x segments`` doubles per temporary.
BATCH = 4096


# ---------------------------------------------------------------------------
# polygon
# ---------------------------------------------------------------------------


def _segment_field(px, py, ax, ay, bx, by, parity=True):
    """Squared distance from each point to the nearest segment, and (with
    ``parity``) whether a ray from the point crosses an odd number of them.

    The ray runs radially outward from each point -- along ``p / |p|``, or
    along +x from the origin -- so after the symmetry fold it stays inside
    the point's own sector.  Each segment is tested with a half-open
    straddle rule on the side of the ray its endpoints fall, so a ray
    through a shared vertex is counted exactly once.
    """
    n = px.shape[0]
    d2 = np.empty(n)
    odd = np.zeros(n, dtype=bool)
    ex, ey = bx - ax, by - ay
    ee = ex * ex + ey * ey
    ee = np.where(ee > 0.0, ee, 1e-300)
    for s in range(0, n, BATCH):
        X = px[s:s + BATCH, None]
        Y = py[s:s + BATCH, None]
        wx, wy = X - ax, Y - ay
        t = np.clip((wx * ex + wy * ey) / ee, 0.0, 1.0)
        dx, dy = wx - t * ex, wy - t * ey
        d2[s:s + BATCH] = np.min(dx * dx + dy * dy, axis=1)
        if not parity:
            continue
        length = np.hypot(X, Y)
        at_origin = length == 0.0
        ux = np.where(at_origin, 1.0, X / np.where(at_origin, 1.0, length))
        uy = np.where(at_origin, 0.0, Y / np.where(at_origin, 1.0, length))
        # Side of the ray each endpoint lies on, and where the segment
        # meets the ray's line, as a distance along the ray.
        ca = ux * (ay - Y) - uy * (ax - X)
        cb = ux * (by - Y) - uy * (bx - X)
        straddle = (ca >= 0.0) != (cb >= 0.0)
        denom = ux * ey - uy * ex
        with np.errstate(divide="ignore", invalid="ignore"):
            along = (wy * ex - wx * ey) / denom     # cross(a - p, e) / denom
        hit = straddle & (along > 0.0)
        odd[s:s + BATCH] = (np.count_nonzero(hit, axis=1) % 2) == 1
    return d2, odd


_SEGMENTS = {}


def _split(a, b):
    return (a[:, 0].copy(), a[:, 1].copy(), b[:, 0].copy(), b[:, 1].copy())


def _segments(points, symmetry):
    """The polygon's segments, organised for evaluation.

    Returns ``(core, rest, boxes)``.  For an asymmetric polygon ``core`` is
    every segment and ``rest`` is empty.  For an ``n``-fold symmetric one,
    a folded point lies within half a sector of angle 0, so

    * ``core`` is every segment whose angular extent reaches into that
      sector.  A point's radial ray stays at its own angle, so it can only
      cross these: they alone decide the sign, and give a first distance.
    * ``rest`` is every other segment still within reach -- one more sector
      either side -- with ``boxes``, the bounding box of each side's share.
      A point needs them only if a box is nearer than its core distance.
    """
    key = (points, symmetry)
    hit = _SEGMENTS.get(key)
    if hit is not None:
        return hit
    a = np.asarray(points, dtype=float)
    b = np.roll(a, -1, axis=0)
    if symmetry == 1:
        hit = (_split(a, b), None, ())
    else:
        half = math.pi / symmetry
        aa = np.arctan2(a[:, 1], a[:, 0])
        ab = np.arctan2(b[:, 1], b[:, 0])
        sweep = (ab - aa + math.pi) % (2 * math.pi) - math.pi
        lo = np.minimum(aa, aa + sweep)
        hi = np.maximum(aa, aa + sweep)
        # Unwrap each extent to lie around angle 0 where possible.
        shift = np.round(0.5 * (lo + hi) / (2 * math.pi)) * 2 * math.pi
        lo, hi = lo - shift, hi - shift
        eps = 1e-9
        core = (hi >= -half - eps) & (lo <= half + eps)
        span = float(np.max(hi - lo))
        reach = (hi >= -3 * half - span) & (lo <= 3 * half + span)
        rest = reach & ~core
        boxes = []
        mid = 0.5 * (lo + hi)
        for side in (mid < 0.0, mid >= 0.0):
            pick = rest & side
            if pick.any():
                pts = np.vstack([a[pick], b[pick]])
                boxes.append((pts.min(axis=0), pts.max(axis=0)))
        hit = (_split(a[core], b[core]),
               _split(a[rest], b[rest]) if rest.any() else None,
               tuple(boxes))
    if len(_SEGMENTS) > 64:
        _SEGMENTS.clear()
    _SEGMENTS[key] = hit
    return hit


def _box_distance2(px, py, box):
    (x0, y0), (x1, y1) = box
    dx = np.maximum(np.maximum(x0 - px, px - x1), 0.0)
    dy = np.maximum(np.maximum(y0 - py, py - y1), 0.0)
    return dx * dx + dy * dy


def _analyze_polygon(params, _children):
    pts = np.asarray(params["points"], dtype=float)
    if params["symmetry"] > 1:
        r = float(np.hypot(pts[:, 0], pts[:, 1]).max())
        lo, hi = (-r, -r), (r, r)
    else:
        lo, hi = pts.min(axis=0), pts.max(axis=0)
    return FieldProps(
        exact=True,
        lipschitz=1.0,
        bounds=((float(lo[0]), float(lo[1]), -INF),
                (float(hi[0]), float(hi[1]), INF)),
        csg_exact=False,
    )


def _eval_polygon(params, p, _children, _ev):
    px, py = p[:, 0], p[:, 1]
    n = params["symmetry"]
    if n > 1:
        sector = 2.0 * math.pi / n
        k = np.round(np.arctan2(py, px) / sector)
        c, s = np.cos(-k * sector), np.sin(-k * sector)
        px, py = c * px - s * py, s * px + c * py
    core, rest, boxes = _segments(params["points"], n)
    d2, odd = _segment_field(px, py, *core)
    if rest is not None:
        near = np.zeros(len(px), dtype=bool)
        for box in boxes:
            near |= _box_distance2(px, py, box) < d2
        if near.any():
            d2r, _ = _segment_field(px[near], py[near], *rest, parity=False)
            d2[near] = np.minimum(d2[near], d2r)
    return np.where(odd, -1.0, 1.0) * np.sqrt(d2)


register(NodeSpec(
    kind="polygon",
    min_children=0,
    max_children=0,
    analyze=_analyze_polygon,
    backends={"numpy": _eval_polygon},
))


def _signed_area(pts):
    x = np.asarray([q[0] for q in pts])
    y = np.asarray([q[1] for q in pts])
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def polygon(points, symmetry=1):
    """A closed 2D polygon in the xy plane, as a field constant along z.

    :param points: the vertices, in either winding; a repeated closing
        vertex is dropped.  The polygon must not self-intersect.
    :param symmetry: ``n`` when the polygon is unchanged by rotation through
        ``2 pi / n`` about the origin -- a gear's tooth count -- which lets
        evaluation skip all but the nearby segments.  Not checked beyond a
        sample, so pass it only when it is true.
    """
    try:
        pts = [(float(q[0]), float(q[1])) for q in points]
    except (TypeError, ValueError, IndexError):
        raise SdfError("polygon: points must be (x, y) pairs") from None
    if len(pts) > 1 and pts[0] == pts[-1]:
        pts = pts[:-1]
    if len(pts) < 3:
        raise SdfError("polygon: needs at least three vertices")
    if not all(math.isfinite(v) for q in pts for v in q):
        raise SdfError("polygon: vertices must be finite")
    if abs(_signed_area(pts)) == 0.0:
        raise SdfError("polygon: the vertices enclose no area")
    n = int(symmetry)
    if n < 1:
        raise SdfError("polygon: symmetry must be a positive integer")
    if n > 1:
        if len(pts) % n:
            raise SdfError(
                f"polygon: {len(pts)} vertices cannot be {n}-fold symmetric"
            )
        a = 2.0 * math.pi / n
        c, s = math.cos(a), math.sin(a)
        x, y = pts[0]
        rotated = (c * x - s * y, s * x + c * y)
        scale = max(math.hypot(*q) for q in pts)
        if min(math.dist(rotated, q) for q in pts) > 1e-9 * scale:
            raise SdfError(f"polygon: the vertices are not {n}-fold symmetric")
    return make_node("polygon", {"points": pts, "symmetry": n})


# ---------------------------------------------------------------------------
# extrude
# ---------------------------------------------------------------------------


def _analyze_extrude(params, child_props):
    (child,) = child_props
    h = 0.5 * params["height"]
    (x0, y0, _), (x1, y1, _) = child.bounds
    return FieldProps(
        exact=child.exact,
        lipschitz=child.lipschitz,
        bounds=((x0, y0, -h), (x1, y1, h)),
        csg_exact=False,
    )


def _eval_extrude(params, p, children, ev):
    flat = p.copy()
    flat[:, 2] = 0.0
    d = ev(children[0], flat)
    w = np.abs(p[:, 2]) - 0.5 * params["height"]
    return (np.minimum(np.maximum(d, w), 0.0)
            + np.hypot(np.maximum(d, 0.0), np.maximum(w, 0.0)))


register(NodeSpec(
    kind="extrude",
    min_children=1,
    max_children=1,
    analyze=_analyze_extrude,
    backends={"numpy": _eval_extrude},
))


def extrude(child, height):
    """Extrude planar ``child`` -- evaluated in the plane z = 0 -- to
    ``height``, centred on z = 0."""
    try:
        h = float(height)
    except (TypeError, ValueError):
        raise SdfError("extrude: height must be a number") from None
    if not math.isfinite(h) or h <= 0.0:
        raise SdfError("extrude: height must be positive and finite")
    return make_node("extrude", {"height": h}, (child,))


# ---------------------------------------------------------------------------
# apex_extrude
# ---------------------------------------------------------------------------


def _apex_extent(params, child):
    """Half-width of the result's xy footprint, and of the profile."""
    (x0, y0, _), (x1, y1, _) = child.bounds
    profile = max(abs(x0), abs(x1), abs(y0), abs(y1))
    return profile * params["z_hi"] / params["z_ref"], profile


def _analyze_apex(params, child_props):
    (child,) = child_props
    zr, zl, zh = params["z_ref"], params["z_lo"], params["z_hi"]
    footprint, profile = _apex_extent(params, child)
    # g = (z/zr) s(u) with u = zr (x, y)/z; |grad g|^2 = |grad s|^2 +
    # ((s - grad s . u)/zr)^2.  Within the region meshed -- the bounds, plus
    # padding -- z >= zl/2 after the clamp below, so |u| <= 2 R zr / zl,
    # and |s| <= |u| + profile for a polygon around the origin.
    u_max = 2.0 * footprint * zr / zl
    lipschitz = child.lipschitz * math.sqrt(
        1.0 + ((2.0 * u_max + profile) / zr) ** 2)
    return FieldProps(
        exact=False,
        lipschitz=lipschitz,
        bounds=((-footprint, -footprint, zl), (footprint, footprint, zh)),
        csg_exact=False,
    )


def _eval_apex(params, p, children, ev):
    zr, zl, zh = params["z_ref"], params["z_lo"], params["z_hi"]
    z = np.maximum(p[:, 2], 0.5 * zl)      # keep the projection finite
    scale = zr / z
    projected = np.column_stack((p[:, 0] * scale, p[:, 1] * scale,
                                 np.zeros(len(p))))
    cone = ev(children[0], projected) / scale
    slab = np.abs(p[:, 2] - 0.5 * (zl + zh)) - 0.5 * (zh - zl)
    return np.maximum(cone, slab)


register(NodeSpec(
    kind="apex_extrude",
    min_children=1,
    max_children=1,
    analyze=_analyze_apex,
    backends={"numpy": _eval_apex},
))


def apex_extrude(child, z_ref, z_lo, z_hi):
    """The cone through the origin over planar ``child``, which gives the
    section at ``z = z_ref``, clipped to ``z_lo <= z <= z_hi``."""
    try:
        zr, zl, zh = float(z_ref), float(z_lo), float(z_hi)
    except (TypeError, ValueError):
        raise SdfError("apex_extrude: planes must be numbers") from None
    if not all(math.isfinite(v) for v in (zr, zl, zh)):
        raise SdfError("apex_extrude: planes must be finite")
    if not 0.0 < zl < zh or zr <= 0.0:
        raise SdfError("apex_extrude: needs 0 < z_lo < z_hi and z_ref > 0")
    return make_node("apex_extrude",
                     {"z_ref": zr, "z_lo": zl, "z_hi": zh}, (child,))
