"""The numpy evaluation backend for the SDF node DAG.

This is the reference backend and, for now, the only one.  It deliberately
drives evaluation through the ``backends`` mapping on each
:class:`~yapcad.sdf.node.NodeSpec` rather than calling methods on nodes, so
that the GLSL/WGSL emitter planned for Phase 5 registers alongside it on the
same specs instead of growing into a parallel implementation that can drift.

Everything is array-at-a-time: a call evaluates ``(N, 3)`` points and returns
``(N,)`` values.  The scalar, point-at-a-time idiom used elsewhere in yapCAD
is not available here by design — a preview samples the field at 10^6 points
or more, and per-point Python would make that hopeless.
"""

import numpy as np

from yapcad.sdf.node import (
    INF,
    SdfError,
    bounds_expand,
    bounds_is_empty,
    get_spec,
)

DEFAULT_BACKEND = "numpy"


def _as_points(points):
    """Normalise input to a contiguous ``(N, 3)`` float array.

    Accepts ``(N, 3)`` or ``(N, 4)`` arrays — the latter being yapCAD's
    homogeneous point form, whose ``w`` is dropped — and a single bare point
    of length 3 or 4, which is reported back so the caller can return a
    scalar.
    """
    arr = np.asarray(points, dtype=float)
    if arr.ndim == 1:
        if arr.shape[0] not in (3, 4):
            raise SdfError(
                f"a single point must have 3 or 4 components, "
                f"got {arr.shape[0]}"
            )
        return np.ascontiguousarray(arr[None, :3]), True
    if arr.ndim == 2:
        if arr.shape[1] not in (3, 4):
            raise SdfError(
                f"points must have shape (N, 3) or (N, 4), got {arr.shape}"
            )
        return np.ascontiguousarray(arr[:, :3]), False
    raise SdfError(f"points must be 1- or 2-dimensional, got {arr.ndim} dims")


def _eval(node, pts, cache, backend):
    """Evaluate ``node`` at ``pts``, memoising on the DAG."""
    # Keyed by node identity *and* point-set identity, because a domain
    # operator evaluates its child at resampled points.  The point array is
    # kept alive in the cache value so that its id() cannot be recycled for a
    # different array while the entry is live.
    key = (node, id(pts))
    hit = cache.get(key)
    if hit is not None:
        return hit[1]

    spec = get_spec(node.kind)
    try:
        fn = spec.backends[backend]
    except KeyError:
        raise SdfError(
            f"node kind {node.kind!r} has no {backend!r} backend"
        ) from None

    def ev(child, query):
        return _eval(child, np.ascontiguousarray(query), cache, backend)

    values = np.asarray(fn(node.p, pts, node.children, ev), dtype=float)
    if values.shape != (pts.shape[0],):
        raise SdfError(
            f"node kind {node.kind!r} returned shape {values.shape}, "
            f"expected {(pts.shape[0],)}"
        )
    cache[key] = (pts, values)
    return values


def evaluate(node, points, backend=DEFAULT_BACKEND):
    """Evaluate the field of ``node`` at ``points``.

    :param points: ``(N, 3)`` or ``(N, 4)`` array, or a single point.
    :returns: an ``(N,)`` array, or a float for a single point.

    Negative values are inside the solid, positive outside, zero on the
    surface.  The magnitude is the true distance only when ``node.exact``;
    otherwise it is a lower bound, scaled by ``node.lipschitz`` — see
    :func:`safe_step`.
    """
    pts, single = _as_points(points)
    values = _eval(node, pts, {}, backend)
    return float(values[0]) if single else values


def safe_step(node, values):
    """The distance it is always safe to advance a ray, given ``values``.

    This is the operational meaning of the Lipschitz constant: ``|f| / L`` is
    a lower bound on the true distance to the surface, so a sphere tracer
    that steps by this amount can never tunnel through geometry, whether the
    field is exact or not.
    """
    return np.abs(np.asarray(values, dtype=float)) / node.lipschitz


def _default_epsilon(node):
    """A finite-difference step scaled to the node's extent."""
    lo, hi = node.bounds
    spans = [hi[i] - lo[i] for i in range(3)]
    finite = [s for s in spans if s not in (INF, -INF) and s == s]
    if finite and max(finite) > 0.0:
        return 1e-6 * max(finite)
    return 1e-6


def gradient(node, points, eps=None, backend=DEFAULT_BACKEND):
    """Numerical gradient of the field, by central differences.

    :returns: an ``(N, 3)`` array, or a ``(3,)`` array for a single point.

    Six evaluations per point rather than the cheaper four-point tetrahedron
    scheme, because dual contouring in Phase 2 places vertices by solving
    against these normals and the extra accuracy is worth the samples.
    """
    pts, single = _as_points(points)
    h = _default_epsilon(node) if eps is None else float(eps)
    if h <= 0.0:
        raise SdfError("gradient: eps must be positive")
    cache = {}
    out = np.empty_like(pts)
    for axis in range(3):
        step = np.zeros(3)
        step[axis] = h
        ahead = _eval(node, np.ascontiguousarray(pts + step), cache, backend)
        behind = _eval(node, np.ascontiguousarray(pts - step), cache, backend)
        out[:, axis] = (ahead - behind) / (2.0 * h)
    return out[0] if single else out


def normal(node, points, eps=None, backend=DEFAULT_BACKEND):
    """Unit surface normal, pointing out of the solid.

    Degenerate points — the centre of a sphere, a point on the medial axis —
    have no defined gradient; those rows come back as zero rather than as
    ``nan``, so downstream meshing can detect and skip them.
    """
    grad = gradient(node, points, eps=eps, backend=backend)
    single = grad.ndim == 1
    g = grad[None, :] if single else grad
    lengths = np.linalg.norm(g, axis=1)
    safe = lengths > 0.0
    unit = np.zeros_like(g)
    unit[safe] = g[safe] / lengths[safe, None]
    return unit[0] if single else unit


def sample_grid(node, bounds=None, resolution=32, padding=0.0,
                backend=DEFAULT_BACKEND):
    """Sample the field on a regular grid.

    :param bounds: ``((x0,y0,z0), (x1,y1,z1))``; defaults to the node's own
        bounds, which must then be finite.
    :param resolution: sample count, as a scalar or a per-axis 3-sequence.
    :param padding: grow the sampled region by this much on every side, so
        that a surface touching the bound is not clipped by it.
    :returns: ``(values, (xs, ys, zs))`` with ``values`` of shape
        ``(nx, ny, nz)`` indexed ``[i, j, k]``.

    This is the Tier 0 path of design document §6: sample, then hand the grid
    to a mesher.  Phase 2 consumes it.
    """
    region = node.bounds if bounds is None else bounds
    if bounds_is_empty(region):
        raise SdfError("sample_grid: the region to sample is empty")
    if padding:
        region = bounds_expand(region, float(padding))
    lo, hi = region
    for i in range(3):
        if lo[i] in (INF, -INF) or hi[i] in (INF, -INF):
            raise SdfError(
                "sample_grid: this field is unbounded; supply explicit bounds"
            )

    if isinstance(resolution, int):
        counts = (resolution, resolution, resolution)
    else:
        counts = tuple(int(v) for v in resolution)
        if len(counts) != 3:
            raise SdfError("sample_grid: resolution needs 3 components")
    if any(c < 2 for c in counts):
        raise SdfError("sample_grid: resolution must be at least 2 per axis")

    axes = [np.linspace(lo[i], hi[i], counts[i]) for i in range(3)]
    mesh = np.meshgrid(*axes, indexing="ij")
    pts = np.stack([m.ravel() for m in mesh], axis=1)
    values = evaluate(node, pts, backend=backend)
    return values.reshape(counts), tuple(axes)
