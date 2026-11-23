"""Octree sampling of the dual-contouring corner lattice.

This is stage A of ``docs/sdf-adaptive-contouring.md``: a pointer-free
octree over the corner lattice that :func:`yapcad.sdf.contour.dual_contour`
already builds, used to decide which corners need a field sample.  The
output mesh is unchanged -- byte for byte -- but only corners near the
surface are evaluated, where chunk pruning evaluates every corner of any
32-cell chunk the surface comes near.

**Blocks and keys.**  A block at size ``s`` (a power of two, in cells) and
key ``(i, j, k)`` covers the corners ``i*s .. (i+1)*s`` along x, and so on,
clipped to the lattice.  The root block covers everything.  Keys are
integers on the block's own level, so the tree needs no pointers and is
walked a level at a time in key order, which is what keeps it
deterministic.

**The test (R1 of the design note).**  For a block with centre ``c`` and
half-diagonal ``r``, a field value ``|f(c)| > L * (r + h)``, with ``L`` the
node's Lipschitz bound and ``h`` the cell, proves there is no surface
within one cell of any corner of the block.  Then no lattice edge from
those corners crosses the surface, and every corner has the sign of
``f(c)``, so it is filled with that value instead of sampled.  Blocks the
test cannot prove empty are split, down to ``leaf`` cells, whose corners
are then sampled exactly.

**Why the output is identical.**  Dual contouring reads the corner values
through their signs, except at the two ends of an edge the surface
crosses, where it interpolates.  Every corner here has its true sign, and
both ends of every crossing edge lie within a cell of the surface, so they
are never in a proved block and are always sampled exactly.  A point's
value does not depend on the batch it is evaluated in
(:func:`yapcad.sdf.evaluate.evaluate_batches`), so those samples are the
values the dense lattice would hold.
"""

import math
from collections import namedtuple

import numpy as np

from yapcad.sdf.evaluate import DEFAULT_BACKEND, evaluate, evaluate_batches
from yapcad.sdf.node import SdfError

#: What :func:`sample_corners_octree` returns: the corner lattice, and how
#: many block centres and corners it evaluated to fill it.
OctreeSample = namedtuple(
    "OctreeSample", ("values", "centre_evaluations", "corner_evaluations"))

#: Corners evaluated per batch.
_BATCH = 1 << 16

#: Proved blocks at least this many cells across are filled one slice at a
#: time; smaller ones, which are far more numerous, in one vectorised write.
_SLICE_FILL = 8

#: The 8 children of a block, in key order.
_CHILDREN = np.array(
    [(x, y, z) for x in (0, 1) for y in (0, 1) for z in (0, 1)],
    dtype=np.int64)


def root_size(counts):
    """Cells across the root block: the least power of two covering the
    lattice."""
    longest = max(int(c) for c in counts)
    return 1 << max(0, math.ceil(math.log2(longest)))


def sample_corners_octree(node, origin, spacing, counts, leaf=2,
                          backend=DEFAULT_BACKEND):
    """Fill the ``(nx+1, ny+1, nz+1)`` corner lattice, sampling only the
    corners the Lipschitz bound cannot prove far from the surface.

    :param leaf: the block size, in cells, below which blocks are not
        split; their corners are all sampled.  A power of two.
    :returns: an :class:`OctreeSample`.

    A node with no usable Lipschitz bound cannot prove anything, so every
    corner is sampled.
    """
    leaf = int(leaf)
    if leaf < 1 or leaf & (leaf - 1):
        raise SdfError(f"octree leaf must be a power of two, got {leaf}")
    origin = np.asarray(origin, dtype=float)
    spacing = float(spacing)
    counts = np.asarray(counts, dtype=np.int64)
    shape = tuple(int(c) + 1 for c in counts)
    lipschitz = node.lipschitz
    if not (math.isfinite(lipschitz) and lipschitz > 0.0):
        from yapcad.sdf.contour import sample_corners
        dense = sample_corners(node, origin, spacing, tuple(counts),
                               prune=False, backend=backend)
        return OctreeSample(dense, 0, dense.size)

    values = np.full(shape, np.nan, dtype=float)
    size = root_size(counts)
    keys = np.zeros((1, 3), dtype=np.int64)
    centre_evaluations = 0

    while size > leaf:
        lo = keys * size
        hi = np.minimum(lo + size, counts)
        centres = origin + 0.5 * (lo + hi).astype(float) * spacing
        radii = 0.5 * spacing * np.linalg.norm((hi - lo).astype(float),
                                               axis=1)
        found = np.asarray(evaluate(node, centres, backend=backend),
                           dtype=float).reshape(-1)
        centre_evaluations += len(keys)
        far = np.abs(found) > lipschitz * (radii + spacing)
        _fill(values, lo[far], hi[far], found[far], size)

        near = keys[~far]
        size //= 2
        keys = (2 * near[:, None, :] + _CHILDREN[None, :, :]).reshape(-1, 3)
        # Children that start past the last corner hold nothing: the corner
        # on the lattice's far face belongs to the block below it.
        keys = keys[np.all(keys * size < counts, axis=1)]

    corner_evaluations = _sample_leaves(node, values, origin, spacing,
                                        counts, keys, size, backend)
    if np.isnan(values).any():      # pragma: no cover - an internal error
        raise SdfError("octree sampling left corners unfilled")
    return OctreeSample(values, centre_evaluations, corner_evaluations)


def _fill(values, lo, hi, found, size):
    """Write each proved block's centre value over its corners."""
    if len(lo) == 0:
        return
    if size >= _SLICE_FILL:
        for (x0, y0, z0), (x1, y1, z1), v in zip(lo, hi, found):
            values[x0:x1 + 1, y0:y1 + 1, z0:z1 + 1] = v
        return
    for idx, inside in _block_corners(lo, hi, size):
        values[idx[inside, 0], idx[inside, 1], idx[inside, 2]] = \
            found[inside]


def _block_corners(lo, hi, size):
    """Yield, for each corner offset in a block, the corner index of every
    block at that offset and whether it lies inside the (clipped) block.

    One offset at a time keeps the temporaries at the size of the block
    list rather than ``(size + 1) ** 3`` times it.
    """
    for dx in range(size + 1):
        for dy in range(size + 1):
            for dz in range(size + 1):
                idx = lo + np.array((dx, dy, dz), dtype=np.int64)
                yield idx, np.all(idx <= hi, axis=1)


def _sample_leaves(node, values, origin, spacing, counts, keys, size,
                   backend):
    """Sample every corner of the unproved leaf blocks, each once."""
    if len(keys) == 0:
        return 0
    lo = keys * size
    hi = np.minimum(lo + size, counts)
    wanted = np.zeros(values.shape, dtype=bool)
    for idx, inside in _block_corners(lo, hi, size):
        wanted[idx[inside, 0], idx[inside, 1], idx[inside, 2]] = True
    # C order, so each corner once and the batches in a fixed order.
    flat = np.flatnonzero(wanted)
    del wanted
    ijk = np.unravel_index(flat, values.shape)
    # The same arithmetic as the dense lattice, so the same coordinates.
    points = np.empty((len(flat), 3), dtype=float)
    for a in range(3):
        points[:, a] = origin[a] + ijk[a].astype(float) * spacing
    batches = [points[start:start + _BATCH]
               for start in range(0, len(points), _BATCH)]
    sampled = np.concatenate(
        evaluate_batches(node, batches, backend=backend))
    values.reshape(-1)[flat] = sampled
    return len(flat)
