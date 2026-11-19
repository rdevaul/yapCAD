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

import os
import threading
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from yapcad.sdf.node import (
    INF,
    SdfError,
    bounds_expand,
    bounds_is_empty,
    get_spec,
)

DEFAULT_BACKEND = "numpy"

#: The compiled backend: the tree flattened to a tape and run by one numba
#: kernel (:mod:`yapcad.sdf.tape`).  Much faster on a large tree, and equal
#: to numpy to rounding rather than bit for bit.
NUMBA_BACKEND = "numba"

#: The tape on a GPU through WebGPU (:mod:`yapcad.sdf.tape_wgpu`).  Single
#: precision and experimental: for previews, by name only.
WGPU_BACKEND = "wgpu"

#: Use the compiled backend where it applies and numpy where it does not:
#: numba is not installed, or the tree has a node kind the tape cannot
#: express yet.
AUTO_BACKEND = "auto"


def _compiled(node, backend):
    """The tape evaluator for ``node`` under ``backend``, or ``None`` to
    walk the tree with numpy.

    ``"auto"`` falls back quietly.  ``"numba"`` was asked for by name, so
    it says why it cannot be had.
    """
    if backend not in (NUMBA_BACKEND, AUTO_BACKEND, WGPU_BACKEND):
        return None
    from yapcad.sdf import tape
    if backend == WGPU_BACKEND:
        from yapcad.sdf import tape_wgpu
        flat = tape.compiled(node)
        return lambda pts: tape_wgpu.run(flat, pts)
    if not tape.numba_available():
        if backend == AUTO_BACKEND:
            return None
        raise SdfError(
            "the numba backend needs numba, which is not installed; "
            "pip install \"yapcad[fast]\", or use backend=\"auto\" to "
            "fall back to numpy")
    try:
        flat = tape.compiled(node)
    except tape.TapeUnsupported:
        if backend == AUTO_BACKEND:
            return None
        raise
    return lambda pts: tape.run(flat, pts)


def mesh_backend(node, backend=None):
    """The backend that meshes and simplifies ``node``: ``"numba"`` or
    ``"numpy"``.

    ``None`` follows ``YAPCAD_SDF_BACKEND``, which is ``auto`` unless set:
    the compiled backend when numba is installed and the tree can be
    flattened, numpy otherwise.  Set it to ``numpy`` to mesh with the
    reference backend whatever is installed -- the choice for a mesh whose
    bytes have to match one made on another machine, since the compiled
    backend agrees with numpy to rounding and not bit for bit.  The backend
    that made a mesh is written into its meshing record.

    The GPU backend is refused: it computes in single precision, which is
    enough for an image and not for a mesh.
    """
    if backend is None:
        backend = (os.environ.get("YAPCAD_SDF_BACKEND", "").strip().lower()
                   or AUTO_BACKEND)
    if backend == WGPU_BACKEND:
        raise SdfError(
            "the wgpu backend computes in single precision and is for "
            "previews; mesh with numba or numpy")
    if backend not in (AUTO_BACKEND, NUMBA_BACKEND, DEFAULT_BACKEND):
        raise SdfError(
            f"unknown meshing backend {backend!r}; expected auto, numba or "
            f"numpy")
    return resolve_backend(node, backend)


def resolve_backend(node, backend=AUTO_BACKEND):
    """The backend that will actually evaluate ``node``: ``"numba"`` or
    ``"numpy"`` for ``"auto"``, and the name itself otherwise."""
    if backend == AUTO_BACKEND:
        return NUMBA_BACKEND if _compiled(node, backend) else DEFAULT_BACKEND
    return backend


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


#: Points per piece when a batch is split across threads.  Large enough
#: that numpy's per-call overhead stays small against the arithmetic, small
#: enough that a piece's intermediates stay in cache.
THREAD_CHUNK = 25000

_pool = None
_pool_key = None
_pool_lock = threading.Lock()


def thread_count():
    """Threads a field evaluation is spread over.

    ``YAPCAD_SDF_THREADS`` sets it; ``1`` turns threading off.  By default
    it is the number of CPUs this process may use, up to eight -- past
    that the pieces contend for memory bandwidth more than they gain.
    """
    value = os.environ.get("YAPCAD_SDF_THREADS", "").strip()
    if value:
        try:
            return max(1, int(value))
        except ValueError:
            raise SdfError(
                f"YAPCAD_SDF_THREADS must be a whole number, got {value!r}"
            ) from None
    try:
        available = len(os.sched_getaffinity(0))
    except AttributeError:      # not every platform has it
        available = os.cpu_count() or 1
    return max(1, min(available, 8))


def _executor(workers):
    """The shared pool, rebuilt if the worker count or the process changed.

    A pool does not survive a fork: the child would inherit its bookkeeping
    without its threads and wait on them for ever, and parts are meshed in
    worker processes often enough for that to matter.
    """
    global _pool, _pool_key
    key = (os.getpid(), workers)
    with _pool_lock:
        if _pool_key != key:
            _pool = ThreadPoolExecutor(
                max_workers=workers, thread_name_prefix="yapcad-sdf")
            _pool_key = key
        return _pool


def evaluate_batches(node, batches, backend=DEFAULT_BACKEND):
    """Evaluate ``node`` at each ``(N, 3)`` array of ``batches``.

    The batches are spread over threads -- numpy releases the interpreter
    lock inside its array arithmetic, which is where a field evaluation
    spends its time.  The result is a list in the order of ``batches``, and
    each entry is exactly what evaluating that batch alone would give: the
    threads share nothing but the node, which is immutable.
    """
    batches = list(batches)
    fast = _compiled(node, backend)
    if fast is not None:
        # The compiled kernel already uses every core, and takes one
        # caller at a time.
        return [fast(pts) for pts in batches]
    backend = _walked(backend)
    workers = min(thread_count(), len(batches))
    if workers < 2 or threading.current_thread().name.startswith(
            "yapcad-sdf"):
        return [_eval(node, pts, {}, backend) for pts in batches]
    return list(_executor(workers).map(
        lambda pts: _eval(node, pts, {}, backend), batches))


def _eval_split(node, pts, backend):
    """Evaluate one array of points, in pieces across threads if it is
    large enough to be worth it.

    A point's value does not depend on what else is in its batch, so the
    pieces put back together are the values the whole batch would give.
    """
    if len(pts) < 2 * THREAD_CHUNK or thread_count() < 2:
        return _eval(node, pts, {}, backend)
    pieces = [np.ascontiguousarray(pts[start:start + THREAD_CHUNK])
              for start in range(0, len(pts), THREAD_CHUNK)]
    return np.concatenate(evaluate_batches(node, pieces, backend))


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
    fast = _compiled(node, backend)
    if fast is not None:
        values = fast(pts)
    else:
        values = _eval_split(node, pts, _walked(backend))
    return float(values[0]) if single else values


def _walked(backend):
    """The per-node backend a tree walk uses: ``"auto"`` that could not be
    compiled walks with numpy."""
    return DEFAULT_BACKEND if backend == AUTO_BACKEND else backend


def safe_step(node, values):
    """The distance it is always safe to advance a ray, given ``values``.

    This is the operational meaning of the Lipschitz constant: ``|f| / L`` is
    a lower bound on the true distance to the surface, so a sphere tracer
    that steps by this amount can never tunnel through geometry, whether the
    field is exact or not.
    """
    return np.abs(np.asarray(values, dtype=float)) / node.lipschitz


#: How much wider the finite-difference step is in single precision, as a
#: multiple of the default: 1e-3 of the extent instead of 1e-6.  Wide enough
#: for the difference to carry four digits, and still well under a pixel of
#: any frame an image is drawn at.
_SINGLE_PRECISION_STEP = 1000.0


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
    if backend == WGPU_BACKEND:
        # Single precision resolves about 1e-7 of the extent, so a step
        # of 1e-6 of it differences two values that agree to a digit.
        h = max(h, _SINGLE_PRECISION_STEP * _default_epsilon(node))
    fast = _compiled(node, backend)
    if fast is None:
        walked = _walked(backend)

        def fast(query):
            return _eval_split(node, np.ascontiguousarray(query), walked)
    out = np.empty_like(pts)
    for axis in range(3):
        step = np.zeros(3)
        step[axis] = h
        out[:, axis] = (fast(pts + step) - fast(pts - step)) / (2.0 * h)
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
