"""Repeating a field on a rectilinear lattice.

``repeat(part, period)`` is every copy of ``part`` translated by a whole
number of periods along each repeated axis -- a row, a sheet or a block of
them, endless or counted.  It is a domain operator: a query point is folded
back into the cell it falls in and the child is evaluated there, so a
pattern of ten thousand copies costs what a handful cost and a tree holding
one stays a few hundred bytes.

Folding alone is only right when each copy stays inside its own cell.  A
copy that reaches into its neighbours' cells -- interlocking stitches,
overlapping scales, a lattice of touching struts -- can be the nearest
surface to a point that folds into the cell next door.  So the child's
bounds decide how many neighbouring copies each point is tested against:
none for a copy that fits its cell, one ring for a copy up to a period and a
half wide, and so on.  Copies beyond that ring are not evaluated, and the
field is capped by a lower bound on the distance to their bounding boxes.
The cap keeps the field a lower bound on distance, which is what the mesher
and the ray tracer step by, and it is what keeps the field continuous: a
copy that starts or stops being tested as a point crosses from one cell to
the next does so underneath it.
"""

import itertools
import math

import numpy as np

from yapcad.sdf.node import (
    INF,
    FieldProps,
    NodeSpec,
    SdfError,
    make_node,
    register,
)

__all__ = ["repeat", "repeat_polar"]

#: The least gap, in periods, kept between a point and the nearest copy that
#: is not evaluated.  Smaller makes the cap bite close to the surface.
_MARGIN = 0.25


def _axes(params):
    """(axis, period, count) for each repeated axis; count 0 is endless."""
    return [(i, params["period"][i], params["count"][i])
            for i in range(3) if params["period"][i] > 0.0]


def _reach(period, half_width):
    """How many neighbouring copies each side a point must be tested
    against, for a copy ``half_width`` either side of its cell's centre."""
    reach = max(0, math.ceil(half_width / period - 0.5))
    if (reach + 0.5) * period - half_width < _MARGIN * period:
        reach += 1
    return reach


def _analyze_repeat(params, child_props):
    (child,) = child_props
    lo, hi = list(child.bounds[0]), list(child.bounds[1])
    for axis, period, count in _axes(params):
        if lo[axis] == -INF or hi[axis] == INF:
            raise SdfError(
                "repeat: the part is unbounded along "
                f"{'xyz'[axis]}, the axis it is repeated on; bound it first"
            )
        if count:
            hi[axis] += (count - 1) * period
        else:
            lo[axis], hi[axis] = -INF, INF
    return FieldProps(
        exact=False,
        lipschitz=child.lipschitz,
        bounds=(tuple(lo), tuple(hi)),
        csg_exact=False,
    )


def _eval_repeat(params, p, children, ev):
    child = children[0]
    (lo, hi) = child.bounds
    axes = _axes(params)
    if not axes:
        return ev(child, p)

    # For each repeated axis: the cell each point falls in, the copies either
    # side of it worth testing, and the cap -- a lower bound on the distance
    # to every copy that is not tested.
    home, offsets, cap = {}, [], np.full(len(p), INF)
    for axis, period, count in axes:
        centre = 0.5 * (lo[axis] + hi[axis])
        half = 0.5 * (hi[axis] - lo[axis])
        reach = _reach(period, half)
        s = p[:, axis] - centre
        cell = np.rint(s / period)
        home[axis] = np.clip(cell, 0, count - 1) if count else cell
        offsets.append(range(-reach, reach + 1))
        # The nearest copy not tested has its box (reach + 1) periods from
        # this cell's centre, less its half width.  Measured from the cell
        # the point is in, that distance is least on the cell's boundary,
        # where it is `edge` -- and there it is the same from both sides,
        # which is what keeps the field continuous as the set of copies
        # tested changes from one cell to the next.
        edge = (reach + 0.5) * period - half
        bound = edge + 0.5 * period - np.abs(s - cell * period)
        if count:
            # Past either end of a counted pattern no copy is nearer than
            # the end cell's boundary, and the bound grows with distance.
            before = -0.5 * period - s
            after = s - (count - 0.5) * period
            past = np.maximum(before, after)
            bound = np.where(past > 0.0, edge + past, bound)
        cap = np.minimum(cap, bound)

    out = np.full(len(p), INF)
    for combo in itertools.product(*offsets):
        q = p.copy()
        live = np.ones(len(p), dtype=bool)
        for (axis, period, count), step in zip(axes, combo):
            k = home[axis] + step
            if count:
                live &= (k >= 0) & (k <= count - 1)
            q[:, axis] -= k * period
        if live.all():
            out = np.minimum(out, ev(child, q))
        elif live.any():
            out[live] = np.minimum(out[live], ev(child, q[live]))
    return np.minimum(out, cap)


register(NodeSpec(
    kind="repeat",
    min_children=1,
    max_children=1,
    analyze=_analyze_repeat,
    backends={"numpy": _eval_repeat},
))


def _triple(name, value, cast):
    if isinstance(value, (int, float)):
        value = (value, value, value)
    try:
        out = tuple(cast(v) for v in value)
    except (TypeError, ValueError):
        raise SdfError(f"repeat: {name} must be a number or three") from None
    if len(out) != 3:
        raise SdfError(f"repeat: {name} must be a number or three")
    return out


def repeat(child, period, count=None):
    """Copies of ``child`` on a rectilinear lattice.

    :param period: the spacing of the copies: a number for the same spacing
        on all three axes, or ``(px, py, pz)``.  A period of 0 leaves that
        axis alone, so ``(4, 0, 0)`` is a row, ``(4, 3, 0)`` a sheet and
        ``(4, 3, 5)`` a block.
    :param count: how many copies along each repeated axis: a number, or
        ``(nx, ny, nz)``.  The first copy is ``child`` where it stands and
        the rest follow along +x, +y and +z.  ``None`` or 0 on an axis is
        endless in both directions; such a field is unbounded there and has
        to be intersected with a bounded part before it is meshed.

    The copies may overlap; see the module notes for what that costs.
    """
    px = _triple("period", period, float)
    if not all(math.isfinite(v) and v >= 0.0 for v in px):
        raise SdfError("repeat: periods must be finite and not negative")
    if not any(v > 0.0 for v in px):
        raise SdfError("repeat: at least one period must be positive")
    nx = (0, 0, 0) if count is None else _triple("count", count, int)
    if any(v < 0 for v in nx):
        raise SdfError("repeat: counts must not be negative")
    nx = tuple(n if p > 0.0 else 1 for p, n in zip(px, nx))
    return make_node("repeat", {"period": list(px), "count": list(nx)},
                     (child,))


# ---------------------------------------------------------------------------
# repeat_polar
# ---------------------------------------------------------------------------


def _polar_layout(child, count):
    """(centre angle, angular half width, reach) of the child's copies.

    The child is read through its bounding box.  A box that reaches the
    axis subtends every direction, and then every copy has to be tested.
    """
    step = 2.0 * math.pi / count
    (x0, y0, _), (x1, y1, _) = child.bounds
    finite = all(v not in (INF, -INF) for v in (x0, y0, x1, y1))
    if not finite or (x0 <= 0.0 <= x1 and y0 <= 0.0 <= y1):
        return 0.0, math.pi, count
    middle = math.atan2(0.5 * (y0 + y1), 0.5 * (x0 + x1))
    half = max(abs(math.remainder(math.atan2(y, x) - middle, 2.0 * math.pi))
               for x in (x0, x1) for y in (y0, y1))
    reach = max(0, math.ceil(half / step - 0.5))
    if (reach + 0.5) * step - half < _MARGIN * step:
        reach += 1
    return middle, half, reach


def _analyze_repeat_polar(_params, child_props):
    (child,) = child_props
    (x0, y0, z0), (x1, y1, z1) = child.bounds
    if INF in (x1, y1) or -INF in (x0, y0):
        raise SdfError(
            "repeat_polar: the part is unbounded across the axis it is "
            "repeated round; bound it first"
        )
    r = max(math.hypot(x, y) for x in (x0, x1) for y in (y0, y1))
    return FieldProps(
        exact=False,
        lipschitz=child.lipschitz,
        bounds=((-r, -r, z0), (r, r, z1)),
        csg_exact=False,
    )


def _eval_repeat_polar(params, p, children, ev):
    child = children[0]
    count = params["count"]
    step = 2.0 * math.pi / count
    middle, half, reach = _polar_layout(child, count)

    def turned(k):
        """The points seen from copy ``k``'s own position."""
        c, s = np.cos(-k * step), np.sin(-k * step)
        q = p.copy()
        q[:, 0] = c * p[:, 0] - s * p[:, 1]
        q[:, 1] = s * p[:, 0] + c * p[:, 1]
        return q

    if 2 * reach + 1 >= count:
        out = ev(child, p)
        for k in range(1, count):
            out = np.minimum(out, ev(child, turned(k)))
        return out

    rho = np.hypot(p[:, 0], p[:, 1])
    angle = np.arctan2(p[:, 1], p[:, 0]) - middle
    sector = np.rint(angle / step)
    u = angle - sector * step
    out = np.full(len(p), INF)
    for offset in range(-reach, reach + 1):
        out = np.minimum(out, ev(child, turned(sector + offset)))
    # The nearest copy not tested lies in a wedge that starts this far round
    # from the point; the distance to that wedge bounds the distance to the
    # copy, and is the same seen from either side of a sector boundary.
    gap = (reach + 1) * step - half - np.abs(u)
    cap = np.where(gap < 0.5 * math.pi, rho * np.sin(gap), rho)
    # Close to the axis every wedge is close, and that bound falls to
    # nothing.  But no copy comes nearer the axis than the part's own
    # bounds do, which bounds the distance from the other side.
    (x0, y0, _), (x1, y1, _) = child.bounds
    inner = math.hypot(max(x0, -x1, 0.0), max(y0, -y1, 0.0))
    return np.minimum(out, np.maximum(cap, inner - rho))


register(NodeSpec(
    kind="repeat_polar",
    min_children=1,
    max_children=1,
    analyze=_analyze_repeat_polar,
    backends={"numpy": _eval_repeat_polar},
))


def repeat_polar(child, count):
    """``count`` copies of ``child`` evenly spaced round the z axis.

    The first copy is ``child`` where it stands; the rest are turned from
    it by whole multiples of ``360 / count`` degrees.  Fins, spokes, bolt
    circles, portholes: draw one, away from the axis, and repeat it.

    Like :func:`repeat` this folds a query point into its sector rather
    than evaluating every copy, tests as many neighbouring copies as the
    part's angular width calls for, and caps the field by a lower bound on
    the distance to the rest.  A part whose bounds reach the axis has no
    sector of its own, and every copy is tested.
    """
    try:
        n = int(count)
    except (TypeError, ValueError):
        raise SdfError("repeat_polar: count must be a whole number") from None
    if n != count or n < 1:
        raise SdfError("repeat_polar: count must be a whole number, 1 or more")
    if n == 1:
        return child
    return make_node("repeat_polar", {"count": n}, (child,))
