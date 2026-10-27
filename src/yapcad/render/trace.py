"""Sphere tracing a field, in numpy, headlessly.

This is the same step rule the rest of :mod:`yapcad.sdf` is built on.  Because
hard CSG produces a conservative *lower* bound on the true distance rather
than an exact one, a march that advances by ``|f| / L`` can never step through
a surface: ``exact=False`` costs iterations, never correctness.  That is why a
tree of boxes and subtractions can be ray-marched at all.

Two things carry the performance, and both are here:

* **Ray/box clipping.** Rays start at the field's bounding box rather than at
  the eye, and a ray that misses the box is never marched.
* **Live-ray compaction.** Each iteration evaluates only the rays still in
  flight.  Rays finish at wildly different times -- a silhouette ray takes far
  longer than one that lands on a flat face -- so without compaction every
  iteration pays for the slowest ray in the frame.

On a 105-node part at 512x512 the two together were worth 43x -- 1.8 s
against 79 s -- of which compaction is 37x and clipping 2.9x.  What they do
*not* fix is the cost of walking a large tree per evaluation, which is why
assemblies are rendered one part at a time and composited: see
:mod:`yapcad.render.buffers`.
"""

import math
import warnings

import numpy as np

from yapcad.sdf.evaluate import DEFAULT_BACKEND, evaluate, normal

from yapcad.render.buffers import PartBuffer

#: Rays evaluated per call.  Bounds peak memory on a large frame without
#: measurably slowing a small one.
DEFAULT_CHUNK = 1 << 18

#: Hit tolerance as a fraction of a pixel.  Below about a quarter of a pixel
#: the extra iterations buy nothing a viewer can see.
TOLERANCE_FRACTION = 0.45


def slab(origins, dirs, bounds):
    """Ray/box intersection parameters, one pair per ray.

    :return: ``(enter, exit)``; a ray misses the box where
        ``enter >= exit``.  Handles axis-parallel rays and infinite bounds,
        both of which a naive division turns into ``nan``.
    """
    lo = np.asarray(bounds[0], dtype=float)
    hi = np.asarray(bounds[1], dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        inv = 1.0 / dirs
        ta = (lo - origins) * inv
        tb = (hi - origins) * inv
    near = np.minimum(ta, tb)
    far = np.maximum(ta, tb)
    flat = dirs == 0.0
    if flat.any():
        # A ray with no component along an axis is either inside that slab
        # for its whole length or outside it for all of it.  The override has
        # to land after the min/max, or an empty slab comes back inverted and
        # the box swallows every ray.
        inside = (origins >= lo) & (origins <= hi)
        near = np.where(flat, np.where(inside, -np.inf, np.inf), near)
        far = np.where(flat, np.where(inside, np.inf, -np.inf), far)
    return np.maximum(near.max(axis=1), 0.0), far.min(axis=1)


def _chunked(fn, points, chunk):
    if points.shape[0] <= chunk:
        return fn(points)
    return np.concatenate([fn(points[i:i + chunk])
                           for i in range(0, points.shape[0], chunk)])


def march(node, origins, dirs, bounds, far, tol, max_steps=256,
          backend=DEFAULT_BACKEND, chunk=DEFAULT_CHUNK):
    """March every ray and return ``(depth, stats)``.

    ``depth`` is ``inf`` where the ray missed.  ``stats`` reports the
    iterations used, the rays that ran out of iterations while still inside
    the box, and the evaluations spent -- the last being the number to watch
    when a render feels slow.
    """
    lipschitz = float(node.lipschitz)
    if not math.isfinite(lipschitz) or lipschitz <= 0.0:
        warnings.warn(
            f"field reports a Lipschitz constant of {node.lipschitz}; "
            "marching with L = 1, which may step through thin features",
            RuntimeWarning, stacklevel=3)
        lipschitz = 1.0

    enter, leave = slab(origins, dirs, bounds)
    leave = np.minimum(leave, float(far))
    depth = np.full(origins.shape[0], np.inf)
    t = enter.copy()
    live = np.nonzero(enter < leave)[0]
    min_step = 0.5 * tol
    evaluations = 0
    used = 0
    while used < max_steps and live.size:
        used += 1
        points = origins[live] + t[live, None] * dirs[live]
        values = _chunked(lambda p: evaluate(node, p, backend=backend),
                          points, chunk)
        evaluations += values.size
        hit = np.abs(values) < tol
        if hit.any():
            depth[live[hit]] = t[live[hit]]
        t[live] += np.maximum(np.abs(values) / lipschitz, min_step)
        live = live[~hit & (t[live] < leave[live])]
    stats = {"rays": int(origins.shape[0]),
             "clipped": int(origins.shape[0] - int((enter < leave).sum())),
             "steps": int(used),
             "evaluations": int(evaluations),
             "unconverged": int(live.size),
             "tolerance": float(tol)}
    return depth, stats


def render_part(node, camera, name="", group="", tol=None, max_steps=256,
                backend=DEFAULT_BACKEND, chunk=DEFAULT_CHUNK, bounds=None,
                epsilon=None):
    """Fill a :class:`~yapcad.render.buffers.PartBuffer` from a field.

    This is the backend side of the contract in
    :mod:`yapcad.render.buffers`: one part in, one depth-and-normal buffer
    out, no knowledge of what else is in the scene.

    :param node: the field to trace.
    :param camera: the shared :class:`~yapcad.render.camera.Camera`.
    :param bounds: clip box; defaults to the field's own bounds, falling back
        to the camera's reach for an unbounded field.
    :param tol: hit tolerance in millimetres; defaults to a fraction of a
        pixel.
    :param epsilon: finite-difference step for the normals; defaults to
        :mod:`yapcad.sdf`'s own choice for the tree.
    """
    origins, dirs = camera.rays()
    box = node.bounds if bounds is None else bounds
    lo = np.asarray(box[0], dtype=float)
    hi = np.asarray(box[1], dtype=float)
    if not (np.all(np.isfinite(lo)) and np.all(np.isfinite(hi))):
        # An unbounded field still has to stop somewhere: fall back to the
        # region the camera can see at all.
        reach = camera.far
        centre = np.asarray(camera.target, dtype=float)
        lo = np.where(np.isfinite(lo), lo, centre - reach)
        hi = np.where(np.isfinite(hi), hi, centre + reach)
    tolerance = TOLERANCE_FRACTION * camera.mm_per_px if tol is None \
        else float(tol)
    depth, stats = march(node, origins, dirs, (lo, hi), camera.far,
                         tolerance, max_steps=max_steps, backend=backend,
                         chunk=chunk)

    normals = np.zeros((depth.size, 3))
    index = np.nonzero(np.isfinite(depth))[0]
    if index.size:
        points = origins[index] + depth[index, None] * dirs[index]
        normals[index] = _chunked(
            lambda p: normal(node, p, eps=epsilon, backend=backend),
            points, chunk)
    stats["hits"] = int(index.size)

    if stats["unconverged"] > 0.01 * stats["rays"]:
        warnings.warn(
            f"{stats['unconverged']} of {stats['rays']} rays did not converge "
            f"in {max_steps} steps while rendering {name or 'the field'}; the "
            "image may have holes along grazing surfaces. Raise max_steps, or "
            "tol if the field is not Lipschitz-tight.",
            RuntimeWarning, stacklevel=2)

    height, width = camera.height, camera.width
    return PartBuffer(depth=depth.reshape(height, width),
                      normal=normals.reshape(height, width, 3),
                      name=name, group=group, stats=stats)
