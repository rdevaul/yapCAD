"""Spline helpers for yapCAD.

Provides evaluation and sampling routines for spline primitives defined in
:mod:`yapcad.geom`, including Catmull-Rom and NURBS curves.
"""

from __future__ import annotations

from math import pow
from typing import Iterable, List, Sequence, Tuple

from yapcad.geom import point
from yapcad.geometry_utils import to_vec3

Vec3 = Tuple[float, float, float]


def is_catmullrom(curve) -> bool:
    """Return ``True`` if *curve* is a Catmull-Rom spline definition."""

    return isinstance(curve, list) and len(curve) == 3 and curve[0] == 'catmullrom'




def evaluate_catmullrom(curve, u: float) -> list:
    """Evaluate a Catmull-Rom spline at parameter ``u`` in ``[0, 1]``."""

    if not is_catmullrom(curve):
        raise ValueError('curve is not a Catmull-Rom spline')
    _, pts, meta = curve
    ctrl = [point(p) for p in pts]
    count = len(ctrl)
    if count == 0:
        raise ValueError('Catmull-Rom spline has no control points')
    if count == 1:
        return point(ctrl[0])

    closed = bool(meta.get('closed', False))
    segment_count = count if closed else count - 1
    if segment_count <= 0:
        return point(ctrl[-1])

    u_clamped = max(0.0, min(1.0, float(u)))
    span = u_clamped * segment_count
    idx = int(span)
    tau = span - idx
    if idx >= segment_count:
        idx = segment_count - 1
        tau = 1.0

    p0 = ctrl[(idx - 1) % count] if closed else ctrl[max(idx - 1, 0)]
    p1 = ctrl[idx % count]
    p2 = ctrl[(idx + 1) % count] if closed else ctrl[min(idx + 1, count - 1)]
    p3 = ctrl[(idx + 2) % count] if closed else ctrl[min(idx + 2, count - 1)]

    alpha = float(meta.get('alpha', 0.5))
    return _catmullrom_point(p0, p1, p2, p3, alpha, tau)


def evaluate_nurbs(curve, u: float) -> list:
    """Evaluate a NURBS curve at parameter ``u`` in ``[0, 1]``."""

    if not is_nurbs(curve):
        raise ValueError('curve is not a NURBS definition')
    _, ctrl, meta = curve
    degree = int(meta['degree'])
    knots = meta['knots']
    weights = meta['weights']
    u_start = knots[degree]
    u_end = knots[-degree - 1]
    u_clamped = max(0.0, min(1.0, float(u)))
    real_u = u_start + (u_end - u_start) * u_clamped
    return _nurbs_point([point(p) for p in ctrl], weights, knots, degree, real_u)

def sample_catmullrom(curve, *, segments_per_span: int = 12) -> List[list]:
    """Sample a Catmull-Rom spline into a list of :func:`point` values."""

    if segments_per_span < 1:
        raise ValueError('segments_per_span must be >= 1')
    if not is_catmullrom(curve):
        raise ValueError('curve is not a Catmull-Rom spline')

    _, pts, meta = curve
    alpha = float(meta.get('alpha', 0.5))
    closed = bool(meta.get('closed', False))
    ctrl = [point(p) for p in pts]
    count = len(ctrl)
    if count < 2:
        raise ValueError('Catmull-Rom spline needs at least 2 control points')

    samples: List[list] = []
    segment_count = count if closed else count - 1
    for i in range(segment_count):
        p0 = ctrl[(i - 1) % count] if closed else ctrl[max(i - 1, 0)]
        p1 = ctrl[i % count]
        p2 = ctrl[(i + 1) % count] if closed else ctrl[min(i + 1, count - 1)]
        p3 = ctrl[(i + 2) % count] if closed else ctrl[min(i + 2, count - 1)]

        if not samples:
            samples.append(point(p1))

        for step in range(1, segments_per_span + 1):
            tau = step / segments_per_span
            samples.append(_catmullrom_point(p0, p1, p2, p3, alpha, tau))

    return samples


def _catmullrom_point(p0, p1, p2, p3, alpha: float, tau: float) -> list:
    v0 = to_vec3(p0)
    v1 = to_vec3(p1)
    v2 = to_vec3(p2)
    v3 = to_vec3(p3)

    def tj(ti: float, pa: Vec3, pb: Vec3) -> float:
        delta = ((pb[0] - pa[0]) ** 2 + (pb[1] - pa[1]) ** 2 + (pb[2] - pa[2]) ** 2) ** 0.5
        return ti + pow(delta, alpha)

    t0 = 0.0
    t1 = tj(t0, v0, v1)
    t2 = tj(t1, v1, v2)
    t3 = tj(t2, v2, v3)

    if t2 - t1 < 1e-12:
        return point(*v2, 1.0)

    t = t1 + (t2 - t1) * tau

    A1 = _catmull_blend(v0, v1, t0, t1, t)
    A2 = _catmull_blend(v1, v2, t1, t2, t)
    A3 = _catmull_blend(v2, v3, t2, t3, t)

    B1 = _catmull_blend(A1, A2, t0, t2, t)
    B2 = _catmull_blend(A2, A3, t1, t3, t)

    C = _catmull_blend(B1, B2, t1, t2, t)
    return point(C[0], C[1], C[2])


def _catmull_blend(a: Vec3, b: Vec3, t0: float, t1: float, t: float) -> Vec3:
    denom = t1 - t0
    if abs(denom) < 1e-12:
        return b
    w0 = (t1 - t) / denom
    w1 = (t - t0) / denom
    return (
        a[0] * w0 + b[0] * w1,
        a[1] * w0 + b[1] * w1,
        a[2] * w0 + b[2] * w1,
    )


def is_nurbs(curve) -> bool:
    """Return ``True`` if *curve* is a NURBS definition."""

    return isinstance(curve, list) and len(curve) == 3 and curve[0] == 'nurbs'


def sample_nurbs(curve, *, samples: int = 64) -> List[list]:
    """Sample a NURBS curve into :func:`point` values."""

    if samples < 2:
        raise ValueError('samples must be >= 2')
    if not is_nurbs(curve):
        raise ValueError('curve is not a NURBS definition')

    _, ctrl_points, meta = curve
    degree = int(meta['degree'])
    knots: Sequence[float] = meta['knots']
    weights: Sequence[float] = meta['weights']
    ctrl = [point(p) for p in ctrl_points]

    u_start = knots[degree]
    u_end = knots[-degree - 1]

    samples_out: List[list] = []
    for i in range(samples):
        if i == samples - 1:
            u = u_end
        else:
            u = u_start + (u_end - u_start) * (i / (samples - 1))
        samples_out.append(_nurbs_point(ctrl, weights, knots, degree, u))
    return samples_out


def _nurbs_point(ctrl: Sequence[list], weights: Sequence[float], knots: Sequence[float], degree: int, u: float) -> list:
    n = len(ctrl) - 1
    numerator = [0.0, 0.0, 0.0]
    denominator = 0.0
    for i in range(n + 1):
        basis = _nip(i, degree, u, knots)
        if basis == 0.0:
            continue
        w = weights[i] * basis
        v = to_vec3(ctrl[i])
        numerator[0] += w * v[0]
        numerator[1] += w * v[1]
        numerator[2] += w * v[2]
        denominator += w
    if denominator == 0.0:
        return point(ctrl[0])
    return point(numerator[0] / denominator, numerator[1] / denominator, numerator[2] / denominator)


def _nip(i: int, p: int, u: float, knots: Sequence[float]) -> float:
    if p == 0:
        if knots[i] <= u < knots[i + 1] or (u == knots[-1] and knots[i] < knots[i + 1]):
            return 1.0
        return 0.0

    left = 0.0
    denom = knots[i + p] - knots[i]
    if denom != 0.0:
        left = (u - knots[i]) / denom * _nip(i, p - 1, u, knots)

    right = 0.0
    denom = knots[i + p + 1] - knots[i + 1]
    if denom != 0.0:
        right = (knots[i + p + 1] - u) / denom * _nip(i + 1, p - 1, u, knots)

    return left + right


__all__ = [
    'is_catmullrom',
    'sample_catmullrom',
    'evaluate_catmullrom',
    'is_nurbs',
    'sample_nurbs',
    'evaluate_nurbs',
]
