"""Helical threads and hex nuts as fields.

A screw thread is a surface ``r = R(u)`` in the coordinates ``u = z - hand *
lead * theta / 2 pi``, folded modulo the pitch, where ``R`` is the thread
profile over one pitch: crest flat, flank, root flat, flank.  The ``thread``
node's field is the distance from ``(u, r)`` to that profile in the
``(u, r)`` half-plane, negative where ``r < R(u)`` -- inside the thread.

That is not a true 3D distance, because the helical map stretches space:
its largest singular value is ``sqrt(1 + (lead / 2 pi r)^2)``, which is
1.002 at an M8 thread's minor radius, so the field is very nearly exact on
the flanks.  Toward the axis the stretch grows without bound, so inside a
guard radius -- half the minor radius, far from any surface -- the angular
dependence is faded out, which keeps the Lipschitz bound finite.

``hex_nut`` and ``hex_bolt`` are semantic nodes storing a fastener's
dimensions.  A nut is a hexagonal prism minus an internal thread; a bolt is
an external thread from its tip at ``z = 0`` to the thread length, then a
plain shank, a washer face and a hex head -- both placed as
:mod:`yapcad.fasteners_legacy` places them, axis +z, flats parallel to the
x axis.  The legacy meshes taper their threads out over the last fifth of
a lead at each end; the fields end them square.
"""

import math
from functools import lru_cache

import numpy as np

from yapcad.sdf.node import INF, FieldProps, NodeSpec, SdfError, analyze, \
    make_node, register
from yapcad.sdf.ops import intersect, subtract, translate, union
from yapcad.sdf.planar import _segment_field, extrude, polygon
from yapcad.sdf.primitives import box, cylinder

__all__ = ["thread", "thread_profile", "hex_nut", "hex_bolt"]


# ---------------------------------------------------------------------------
# thread
# ---------------------------------------------------------------------------


def _guard_radius(profile):
    return 0.5 * min(r for _, r in profile)


def _analyze_thread(params, _children):
    profile = params["profile"]
    r_max = max(r for _, r in profile)
    r_min = min(r for _, r in profile)
    guard = _guard_radius(profile)
    lead = params["pitch"] * params["starts"]
    # Outside the guard: the helical map's stretch at the guard radius.
    # Inside it: the radial term picks up the profile's depth over the
    # guard radius, and the faded angular term is bounded by its value at
    # the guard radius.
    radial = 1.0 + (r_max - r_min) / guard
    angular = lead / (2.0 * math.pi * guard)
    return FieldProps(
        exact=False,
        lipschitz=math.sqrt(radial ** 2 + angular ** 2),
        bounds=((-r_max, -r_max, -INF), (r_max, r_max, INF)),
        csg_exact=False,
    )


@lru_cache(maxsize=32)
def _period_segments(profile, pitch):
    """The profile's segments over three periods, so every folded ``u`` in
    [0, pitch) sees its nearest segment."""
    prof = np.asarray(profile, dtype=float)
    a, b = [], []
    for k in (-1, 0, 1):
        q = prof + np.array([k * pitch, 0.0])
        a.append(q[:-1])
        b.append(q[1:])
    a, b = np.vstack(a), np.vstack(b)
    return a[:, 0], a[:, 1], b[:, 0], b[:, 1]


def _eval_thread(params, p, _children, _ev):
    pitch = params["pitch"]
    lead = pitch * params["starts"]
    hand = 1.0 if params["hand"] == "right" else -1.0
    profile = params["profile"]
    prof = np.asarray(profile, dtype=float)
    guard = _guard_radius(profile)

    r = np.hypot(p[:, 0], p[:, 1])
    theta = np.arctan2(p[:, 1], p[:, 0])
    u = np.mod(p[:, 2] - hand * lead * theta / (2.0 * math.pi), pitch)
    rr = np.maximum(r, guard)
    # Unsigned (u, r) distance to the profile; the sign is which side of
    # the graph r = R(u) the point lies.
    d2, _ = _segment_field(u, rr, *_period_segments(profile, pitch))
    radius_at_u = np.interp(u, prof[:, 0], prof[:, 1])
    d = np.where(rr < radius_at_u, -1.0, 1.0) * np.sqrt(d2)

    # Inside the guard radius, blend toward a value with no angular
    # dependence, reaching it on the axis: continuous at the guard radius,
    # and the 1/r growth of the angular derivative is cancelled by w ~ r.
    # At the guard radius d is about guard - R(u), so the blend target is
    # its mean, and the axis value comes out near -mean(R): well inside.
    w = np.clip(r / guard, 0.0, 1.0)
    target = guard - float(prof[:, 1].mean())
    inner = (r - guard) + w * d + (1.0 - w) * target
    return np.where(r >= guard, d, inner)


register(NodeSpec(
    kind="thread",
    min_children=0,
    max_children=0,
    analyze=_analyze_thread,
    backends={"numpy": _eval_thread},
))


def thread(profile, pitch, starts=1, hand="right"):
    """An infinite helical thread about the z axis.

    :param profile: ``(u, r)`` points over one pitch, ``u`` from 0 to
        ``pitch`` and ``R(0) == R(pitch)``; see :func:`thread_profile`.
    :param starts: thread starts; the lead is ``starts * pitch``.
    :param hand: ``"right"`` or ``"left"``.
    """
    try:
        prof = [(float(u), float(r)) for u, r in profile]
        pitch = float(pitch)
        starts = int(starts)
    except (TypeError, ValueError):
        raise SdfError("thread: profile must be (u, r) pairs") from None
    if not math.isfinite(pitch) or pitch <= 0.0:
        raise SdfError("thread: pitch must be positive and finite")
    if starts < 1:
        raise SdfError("thread: starts must be at least 1")
    if hand not in ("right", "left"):
        raise SdfError("thread: hand must be 'right' or 'left'")
    if len(prof) < 3 or prof[0][0] != 0.0 or \
            abs(prof[-1][0] - pitch) > 1e-12 * pitch:
        raise SdfError("thread: profile must run from u = 0 to u = pitch")
    if any(b[0] < a[0] for a, b in zip(prof, prof[1:])):
        raise SdfError("thread: profile u must not decrease")
    if prof[0][1] != prof[-1][1] or min(r for _, r in prof) <= 0.0:
        raise SdfError("thread: profile must be periodic and off the axis")
    return make_node("thread", {"profile": prof, "pitch": pitch,
                                "starts": starts, "hand": hand})


def thread_profile(thread_spec):
    """The ``(u, r)`` profile of a :class:`yapcad.threadgen.ThreadProfile`,
    exactly as :func:`yapcad.threadgen._radius_at` defines it."""
    if thread_spec.taper_ratio:
        raise SdfError("thread_profile: tapered threads are not supported")
    pitch = thread_spec.P_pitch
    crest = thread_spec.D_nominal / 2.0
    root = max(crest - thread_spec.thread_depth_ratio * pitch, 0.0)
    crest_half = max(0.0, thread_spec.crest_flat_ratio * pitch / 2.0)
    root_half = max(0.0, thread_spec.root_flat_ratio * pitch / 2.0)
    flank = max(0.0, (pitch - 2 * crest_half - 2 * root_half) / 2.0)
    x1 = crest_half
    x2 = x1 + flank
    x3 = x2 + 2 * root_half
    x4 = pitch - crest_half
    return [(0.0, crest), (x1, crest), (x2, root), (x3, root), (x4, crest),
            (pitch, crest)]


# ---------------------------------------------------------------------------
# hex nut
# ---------------------------------------------------------------------------


_NUT_FIELDS = ("diameter", "pitch", "width_flat", "thickness", "starts",
               "hand", "crest_flat_ratio", "root_flat_ratio",
               "thread_depth_ratio")


@lru_cache(maxsize=32)
def _nut_tree(items):
    from yapcad.threadgen import ThreadProfile
    p = dict(items)
    circum = p["width_flat"] / 2.0 / math.cos(math.pi / 6.0)
    hexagon = polygon([(circum * math.cos(math.pi / 6 + k * math.pi / 3),
                        circum * math.sin(math.pi / 6 + k * math.pi / 3))
                       for k in range(6)], symmetry=6)
    body = translate(extrude(hexagon, p["thickness"]),
                     (0.0, 0.0, p["thickness"] / 2.0))
    spec = ThreadProfile(D_nominal=p["diameter"], P_pitch=p["pitch"],
                         crest_flat_ratio=p["crest_flat_ratio"],
                         root_flat_ratio=p["root_flat_ratio"],
                         thread_depth_ratio=p["thread_depth_ratio"],
                         internal=True)
    hole = thread(thread_profile(spec), p["pitch"], starts=p["starts"],
                  hand=p["hand"])
    return subtract(body, hole)


def _analyze_nut(params, _children):
    return analyze(_nut_tree(tuple(sorted(params.items()))))


def _eval_nut(params, p, _children, ev):
    return ev(_nut_tree(tuple(sorted(params.items()))), p)


register(NodeSpec(
    kind="hex_nut",
    min_children=0,
    max_children=0,
    analyze=_analyze_nut,
    backends={"numpy": _eval_nut},
))


def hex_nut(thread_spec, width_flat, thickness):
    """A hex nut with an internal thread from a
    :class:`yapcad.threadgen.ThreadProfile`."""
    try:
        wf, t = float(width_flat), float(thickness)
    except (TypeError, ValueError):
        raise SdfError("hex_nut: width_flat and thickness must be numbers") \
            from None
    if not (math.isfinite(wf) and math.isfinite(t)) or wf <= 0 or t <= 0:
        raise SdfError("hex_nut: width_flat and thickness must be positive")
    if thread_spec.D_nominal >= wf:
        raise SdfError("hex_nut: the thread is wider than the flats")
    if thread_spec.handedness not in ("right", "left"):
        raise SdfError("hex_nut: handedness must be 'right' or 'left'")
    params = {
        "diameter": float(thread_spec.D_nominal),
        "pitch": float(thread_spec.P_pitch),
        "width_flat": wf,
        "thickness": t,
        "starts": int(max(1, thread_spec.starts)),
        "hand": thread_spec.handedness,
        "crest_flat_ratio": float(thread_spec.crest_flat_ratio),
        "root_flat_ratio": float(thread_spec.root_flat_ratio),
        "thread_depth_ratio": float(thread_spec.thread_depth_ratio),
    }
    return make_node("hex_nut", params)


# ---------------------------------------------------------------------------
# hex bolt
# ---------------------------------------------------------------------------


def _hexagon(width_flat):
    circum = width_flat / 2.0 / math.cos(math.pi / 6.0)
    return polygon([(circum * math.cos(math.pi / 6 + k * math.pi / 3),
                     circum * math.sin(math.pi / 6 + k * math.pi / 3))
                    for k in range(6)], symmetry=6)


def _slab(r, z0, z1):
    """A box of half-width ``r`` spanning ``z0 <= z <= z1``."""
    return translate(box((2 * r, 2 * r, z1 - z0)), (0.0, 0.0, 0.5 * (z0 + z1)))


def _rod(r, z0, z1):
    return translate(cylinder(r, z1 - z0), (0.0, 0.0, 0.5 * (z0 + z1)))


@lru_cache(maxsize=32)
def _bolt_tree(items):
    from yapcad.threadgen import ThreadProfile
    p = dict(items)
    spec = ThreadProfile(D_nominal=p["diameter"], P_pitch=p["pitch"],
                         crest_flat_ratio=p["crest_flat_ratio"],
                         root_flat_ratio=p["root_flat_ratio"],
                         thread_depth_ratio=p["thread_depth_ratio"])
    radius = p["diameter"] / 2.0
    lt, ls = p["thread_length"], p["shank_length"]
    parts = [intersect(thread(thread_profile(spec), p["pitch"],
                              starts=p["starts"], hand=p["hand"]),
                       _slab(radius, 0.0, lt))]
    if ls > lt:
        parts.append(_rod(p["shank_diameter"] / 2.0, lt, ls))
    top = ls
    if p["washer_thickness"] > 0.0:
        parts.append(_rod(p["washer_diameter"] / 2.0, ls,
                          ls + p["washer_thickness"]))
        top += p["washer_thickness"]
    parts.append(translate(extrude(_hexagon(p["head_flat"]),
                                   p["head_height"]),
                           (0.0, 0.0, top + 0.5 * p["head_height"])))
    return union(*parts)


def _analyze_bolt(params, _children):
    return analyze(_bolt_tree(tuple(sorted(params.items()))))


def _eval_bolt(params, p, _children, ev):
    return ev(_bolt_tree(tuple(sorted(params.items()))), p)


register(NodeSpec(
    kind="hex_bolt",
    min_children=0,
    max_children=0,
    analyze=_analyze_bolt,
    backends={"numpy": _eval_bolt},
))


def hex_bolt(thread_spec, thread_length, shank_length, head_height,
             head_flat, washer_thickness=0.5, washer_diameter=None,
             shank_diameter=None):
    """A hex bolt with an external thread from a
    :class:`yapcad.threadgen.ThreadProfile`; dimensions as
    :class:`yapcad.fasteners_legacy.HexCapScrewSpec` names them."""
    try:
        values = [float(v) for v in (thread_length, shank_length,
                                     head_height, head_flat,
                                     washer_thickness)]
    except (TypeError, ValueError):
        raise SdfError("hex_bolt: dimensions must be numbers") from None
    lt, ls, hh, hf, wt = values
    if not all(math.isfinite(v) for v in values) or min(lt, ls, hh, hf) <= 0 \
            or wt < 0:
        raise SdfError("hex_bolt: dimensions must be positive")
    if lt > ls:
        raise SdfError("hex_bolt: the thread is longer than the shank")
    if thread_spec.handedness not in ("right", "left"):
        raise SdfError("hex_bolt: handedness must be 'right' or 'left'")
    d = float(thread_spec.D_nominal)
    params = {
        "diameter": d,
        "pitch": float(thread_spec.P_pitch),
        "thread_length": lt,
        "shank_length": ls,
        "head_height": hh,
        "head_flat": hf,
        "washer_thickness": wt,
        "washer_diameter": float(washer_diameter or hf),
        "shank_diameter": float(shank_diameter or d),
        "starts": int(max(1, thread_spec.starts)),
        "hand": thread_spec.handedness,
        "crest_flat_ratio": float(thread_spec.crest_flat_ratio),
        "root_flat_ratio": float(thread_spec.root_flat_ratio),
        "thread_depth_ratio": float(thread_spec.thread_depth_ratio),
    }
    return make_node("hex_bolt", params)
