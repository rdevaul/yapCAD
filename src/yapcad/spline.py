"""Spline helpers for yapCAD.

Provides evaluation and sampling routines for spline primitives defined in
:mod:`yapcad.geom`, including Bezier, B-spline, Catmull-Rom and NURBS curves.
"""

from __future__ import annotations

from math import pow, factorial
from typing import Iterable, List, Sequence, Tuple

from yapcad.geom import point
from yapcad.geometry_utils import to_vec3

Vec3 = Tuple[float, float, float]


# =============================================================================
# Bezier Curve Support
# =============================================================================

def is_bezier(curve) -> bool:
    """Return ``True`` if *curve* is a Bezier curve definition."""

    return isinstance(curve, list) and len(curve) == 3 and curve[0] == 'bezier'


def _binomial(n: int, k: int) -> int:
    """Compute binomial coefficient C(n, k) = n! / (k! * (n-k)!)."""
    if k < 0 or k > n:
        return 0
    return factorial(n) // (factorial(k) * factorial(n - k))


def _bernstein(n: int, i: int, t: float) -> float:
    """Compute Bernstein basis polynomial B_{i,n}(t).

    The Bernstein basis polynomials are:
        B_{i,n}(t) = C(n,i) * t^i * (1-t)^(n-i)

    Parameters
    ----------
    n : int
        Degree of the polynomial
    i : int
        Index (0 <= i <= n)
    t : float
        Parameter value (typically in [0, 1])

    Returns
    -------
    float
        Value of the Bernstein polynomial at t
    """
    return _binomial(n, i) * (t ** i) * ((1.0 - t) ** (n - i))


def bezier_point(control_points, t: float) -> list:
    """Evaluate a Bezier curve at parameter t using De Casteljau's algorithm.

    Parameters
    ----------
    control_points : list
        List of control points (yapCAD points or coordinate tuples).
    t : float
        Parameter value in [0, 1]. t=0 returns first control point,
        t=1 returns last control point.

    Returns
    -------
    list
        Point on the curve at parameter t as a yapCAD point.

    Examples
    --------
    >>> cp = [point(0, 0), point(10, 20), point(30, 20), point(40, 0)]
    >>> p = bezier_point(cp, 0.5)  # Point at middle of curve
    """
    pts = [to_vec3(point(p)) for p in control_points]
    n = len(pts)
    if n == 0:
        raise ValueError('bezier_point requires at least one control point')
    if n == 1:
        return point(pts[0][0], pts[0][1], pts[0][2])

    t = max(0.0, min(1.0, float(t)))

    # De Casteljau's algorithm - numerically stable recursive subdivision
    work = list(pts)
    for level in range(n - 1, 0, -1):
        new_work = []
        for i in range(level):
            x = (1.0 - t) * work[i][0] + t * work[i + 1][0]
            y = (1.0 - t) * work[i][1] + t * work[i + 1][1]
            z = (1.0 - t) * work[i][2] + t * work[i + 1][2]
            new_work.append((x, y, z))
        work = new_work

    return point(work[0][0], work[0][1], work[0][2])


def evaluate_bezier(curve, t: float) -> list:
    """Evaluate a Bezier curve definition at parameter ``t`` in [0, 1].

    Parameters
    ----------
    curve : list
        A Bezier curve definition from :func:`yapcad.geom.bezier`.
    t : float
        Parameter value in [0, 1].

    Returns
    -------
    list
        Point on the curve at parameter t.
    """
    if not is_bezier(curve):
        raise ValueError('curve is not a Bezier definition')
    _, pts, _ = curve
    return bezier_point(pts, t)


def bezier_tangent(control_points, t: float) -> list:
    """Compute the tangent vector of a Bezier curve at parameter t.

    The tangent is the first derivative of the curve with respect to t.
    For a Bezier curve of degree n, the derivative is a Bezier curve of
    degree n-1 with control points: n * (P_{i+1} - P_i).

    Parameters
    ----------
    control_points : list
        List of control points (yapCAD points or coordinate tuples).
    t : float
        Parameter value in [0, 1].

    Returns
    -------
    list
        Tangent vector at parameter t as a yapCAD direction vector [x, y, z, 0].

    Examples
    --------
    >>> cp = [point(0, 0), point(10, 20), point(30, 20), point(40, 0)]
    >>> tangent = bezier_tangent(cp, 0.5)
    """
    pts = [to_vec3(point(p)) for p in control_points]
    n = len(pts)
    if n < 2:
        # Constant point has zero derivative
        return [0.0, 0.0, 0.0, 0.0]

    # Compute derivative control points: n * (P_{i+1} - P_i)
    degree = n - 1
    deriv_pts = []
    for i in range(degree):
        dx = degree * (pts[i + 1][0] - pts[i][0])
        dy = degree * (pts[i + 1][1] - pts[i][1])
        dz = degree * (pts[i + 1][2] - pts[i][2])
        deriv_pts.append(point(dx, dy, dz))

    # Evaluate the derivative curve at t
    tang = bezier_point(deriv_pts, t)
    # Return as direction vector (w=0)
    return [tang[0], tang[1], tang[2], 0.0]


def bezier_curve(control_points, segments: int = 32) -> List[list]:
    """Sample a Bezier curve as a polyline.

    Parameters
    ----------
    control_points : list
        List of control points (yapCAD points or coordinate tuples).
    segments : int, optional
        Number of line segments in the output polyline (default 32).
        More segments = smoother curve approximation.

    Returns
    -------
    list
        List of points suitable for use as a yapCAD polyline.

    Examples
    --------
    >>> cp = [point(0, 0), point(10, 20), point(30, 20), point(40, 0)]
    >>> polyline = bezier_curve(cp, segments=32)
    >>> # polyline can be used anywhere a yapCAD polyline is expected
    """
    if segments < 1:
        raise ValueError('segments must be >= 1')

    samples = []
    for i in range(segments + 1):
        t = i / segments
        samples.append(bezier_point(control_points, t))
    return samples


def sample_bezier(curve, *, segments: int = 32) -> List[list]:
    """Sample a Bezier curve definition into a list of points.

    Parameters
    ----------
    curve : list
        A Bezier curve definition from :func:`yapcad.geom.bezier`.
    segments : int, optional
        Number of line segments in output (default 32).

    Returns
    -------
    list
        List of points suitable for use as a yapCAD polyline.
    """
    if not is_bezier(curve):
        raise ValueError('curve is not a Bezier definition')
    _, pts, _ = curve
    return bezier_curve(pts, segments)


# =============================================================================
# B-spline Curve Support
# =============================================================================

def is_bspline(curve) -> bool:
    """Return ``True`` if *curve* is a B-spline curve definition."""

    return isinstance(curve, list) and len(curve) == 3 and curve[0] == 'bspline'


def _uniform_knot_vector(n: int, degree: int, closed: bool = False) -> List[float]:
    """Generate a uniform knot vector for B-spline.

    Parameters
    ----------
    n : int
        Number of control points
    degree : int
        Degree of the B-spline
    closed : bool
        If True, generate periodic knot vector for closed curve

    Returns
    -------
    list
        Uniform knot vector
    """
    if closed:
        # Periodic knot vector for closed curves
        m = n + degree + 1
        return [float(i - degree) / n for i in range(m)]
    else:
        # Clamped uniform knot vector (curve passes through endpoints)
        # First degree+1 knots = 0, last degree+1 knots = 1
        m = n + degree + 1
        knots = []
        for i in range(m):
            if i < degree + 1:
                knots.append(0.0)
            elif i >= m - degree - 1:
                knots.append(1.0)
            else:
                knots.append((i - degree) / (n - degree))
        return knots


def _bspline_basis(i: int, p: int, u: float, knots: Sequence[float]) -> float:
    """Compute B-spline basis function N_{i,p}(u) using Cox-de Boor recursion.

    Parameters
    ----------
    i : int
        Basis function index
    p : int
        Degree
    u : float
        Parameter value
    knots : sequence
        Knot vector

    Returns
    -------
    float
        Value of basis function N_{i,p}(u)
    """
    if p == 0:
        # Base case: N_{i,0}(u) = 1 if knots[i] <= u < knots[i+1], else 0
        if i + 1 >= len(knots):
            return 0.0
        if knots[i] <= u < knots[i + 1]:
            return 1.0
        # Handle the special case where u equals the last knot
        if u == knots[-1] and knots[i] <= u <= knots[i + 1]:
            return 1.0
        return 0.0

    # Recursive case: Cox-de Boor formula
    left = 0.0
    denom1 = knots[i + p] - knots[i]
    if denom1 > 1e-12:
        left = (u - knots[i]) / denom1 * _bspline_basis(i, p - 1, u, knots)

    right = 0.0
    if i + p + 1 < len(knots):
        denom2 = knots[i + p + 1] - knots[i + 1]
        if denom2 > 1e-12:
            right = (knots[i + p + 1] - u) / denom2 * _bspline_basis(i + 1, p - 1, u, knots)

    return left + right


def bspline_point(control_points, t: float, degree: int = 3, closed: bool = False) -> list:
    """Evaluate a B-spline curve at parameter t.

    Parameters
    ----------
    control_points : list
        List of control points.
    t : float
        Parameter value in [0, 1].
    degree : int, optional
        Degree of B-spline (default 3 for cubic).
    closed : bool, optional
        If True, treat as closed (periodic) B-spline.

    Returns
    -------
    list
        Point on the curve at parameter t.
    """
    pts = [to_vec3(point(p)) for p in control_points]
    n = len(pts)
    if n < degree + 1:
        raise ValueError(f'B-spline needs at least {degree + 1} control points for degree {degree}')

    t = max(0.0, min(1.0, float(t)))

    if closed:
        # For closed B-spline, extend control points cyclically
        extended_pts = pts + pts[:degree]
        knots = _uniform_knot_vector(n, degree, closed=True)
        # Map t to the active knot domain
        u_min = knots[degree]
        u_max = knots[n]
        u = u_min + t * (u_max - u_min)

        x = y = z = 0.0
        for i in range(len(extended_pts)):
            basis = _bspline_basis(i, degree, u, knots)
            x += basis * extended_pts[i][0]
            y += basis * extended_pts[i][1]
            z += basis * extended_pts[i][2]
    else:
        # Clamped B-spline
        knots = _uniform_knot_vector(n, degree, closed=False)
        # Map t to knot domain [0, 1] (clamped knots use this range)
        u = t

        x = y = z = 0.0
        for i in range(n):
            basis = _bspline_basis(i, degree, u, knots)
            x += basis * pts[i][0]
            y += basis * pts[i][1]
            z += basis * pts[i][2]

    return point(x, y, z)


def evaluate_bspline(curve, t: float) -> list:
    """Evaluate a B-spline curve definition at parameter ``t`` in [0, 1].

    Parameters
    ----------
    curve : list
        A B-spline curve definition from :func:`yapcad.geom.bspline`.
    t : float
        Parameter value in [0, 1].

    Returns
    -------
    list
        Point on the curve at parameter t.
    """
    if not is_bspline(curve):
        raise ValueError('curve is not a B-spline definition')
    _, pts, meta = curve
    degree = meta.get('degree', 3)
    closed = meta.get('closed', False)
    return bspline_point(pts, t, degree=degree, closed=closed)


def bspline_curve(control_points, degree: int = 3, closed: bool = False,
                  segments: int = 64) -> List[list]:
    """Sample a B-spline curve as a polyline.

    Parameters
    ----------
    control_points : list
        List of control points.
    degree : int, optional
        Degree of B-spline (default 3).
    closed : bool, optional
        If True, create closed B-spline.
    segments : int, optional
        Number of line segments in output (default 64).

    Returns
    -------
    list
        List of points suitable for use as a yapCAD polyline.
    """
    if segments < 1:
        raise ValueError('segments must be >= 1')

    samples = []
    for i in range(segments + 1):
        t = i / segments
        samples.append(bspline_point(control_points, t, degree, closed))

    if closed:
        # Ensure the curve actually closes
        samples[-1] = samples[0]

    return samples


def sample_bspline(curve, *, segments: int = 64) -> List[list]:
    """Sample a B-spline curve definition into a list of points.

    Parameters
    ----------
    curve : list
        A B-spline curve definition from :func:`yapcad.geom.bspline`.
    segments : int, optional
        Number of line segments in output (default 64).

    Returns
    -------
    list
        List of points suitable for use as a yapCAD polyline.
    """
    if not is_bspline(curve):
        raise ValueError('curve is not a B-spline definition')
    _, pts, meta = curve
    degree = meta.get('degree', 3)
    closed = meta.get('closed', False)
    return bspline_curve(pts, degree=degree, closed=closed, segments=segments)


def bspline_tangent(control_points, t: float, degree: int = 3, closed: bool = False) -> list:
    """Compute the tangent vector of a B-spline curve at parameter t.

    Uses numerical differentiation for robustness.

    Parameters
    ----------
    control_points : list
        List of control points.
    t : float
        Parameter value in [0, 1].
    degree : int, optional
        Degree of B-spline (default 3).
    closed : bool, optional
        If True, treat as closed B-spline.

    Returns
    -------
    list
        Tangent vector at parameter t as [x, y, z, 0].
    """
    h = 1e-6
    t1 = max(0.0, t - h)
    t2 = min(1.0, t + h)

    if t1 == t2:
        return [0.0, 0.0, 0.0, 0.0]

    p1 = bspline_point(control_points, t1, degree, closed)
    p2 = bspline_point(control_points, t2, degree, closed)

    dx = (p2[0] - p1[0]) / (t2 - t1)
    dy = (p2[1] - p1[1]) / (t2 - t1)
    dz = (p2[2] - p1[2]) / (t2 - t1)

    return [dx, dy, dz, 0.0]


# =============================================================================
# Catmull-Rom Spline Support
# =============================================================================

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
    # Bezier curves
    'is_bezier',
    'bezier_point',
    'bezier_curve',
    'bezier_tangent',
    'sample_bezier',
    'evaluate_bezier',
    # B-splines
    'is_bspline',
    'bspline_point',
    'bspline_curve',
    'bspline_tangent',
    'sample_bspline',
    'evaluate_bspline',
    # Catmull-Rom splines
    'is_catmullrom',
    'sample_catmullrom',
    'evaluate_catmullrom',
    # NURBS
    'is_nurbs',
    'sample_nurbs',
    'evaluate_nurbs',
]
