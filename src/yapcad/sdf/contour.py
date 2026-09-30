"""Dual contouring: SDF node DAG to an indexed triangle mesh.

Design document §5.1 calls for dual contouring rather than marching cubes,
for one decisive reason: dual contouring places one vertex per surface cell
by solving for the point that best fits the surface planes crossing that
cell, so a box meshes with actual edges and corners instead of the stair-step
chamfer marching cubes produces.  It needs surface normals to do that, which
an SDF supplies for free from its gradient.

Three properties of this implementation are load-bearing.

**It is deterministic.**  yapCAD signs packages, and a preview that
regenerates differently run to run would churn package hashes.  Vertex order
is the C-order of the active-cell mask and triangle order follows the edge
axes in turn, so identical parameters give byte-identical output.  Nothing
here iterates a dict or a set.

**It is crack-free.**  Every surface-crossing grid edge emits exactly one
quad joining the vertices of the four cells around it, so the result is a
closed manifold whenever the sampled region encloses the surface.

**It bounds its own memory.**  The corner lattice is filled in chunks rather
than as one allocation, and chunks that provably cannot contain the surface
are skipped using the node's Lipschitz constant — the first place the Phase 1
bookkeeping pays for itself. A chunk whose centre satisfies
``|f| > L * (circumradius + margin)`` has no surface within ``margin`` of it;
with ``margin`` at least one cell diagonal, no edge touching that chunk can
cross, so filling it with the centre value is sign-correct and safe.

What this is not, yet: the octree here is uniform in *output* resolution.
Chunk pruning makes the work adaptive — empty space costs almost nothing —
but every emitted cell is the same size.  Genuinely adaptive output, with
variable leaf depth and QEF-error-driven collapse, needs the crack-patching
traversal of Ju et al. and is deferred; see the design document.
"""

import math
from collections import namedtuple

import numpy as np

from yapcad.geom import epsilon
from yapcad.sdf.evaluate import DEFAULT_BACKEND, evaluate, gradient
from yapcad.sdf.node import INF, SdfError, bounds_is_empty

#: The result of :func:`dual_contour`.
#:
#: ``vertices`` and ``normals`` are ``(V, 3)``; ``triangles`` is ``(T, 3)``
#: of vertex indices, wound counter-clockwise seen from outside the solid.
Mesh = namedtuple("Mesh", ("vertices", "normals", "triangles"))

#: Eigenvalues below this fraction of a cell's largest are treated as
#: unconstrained directions in the QEF solve (Ju et al.'s truncation).
QEF_TRUNCATION = 1e-2


# For an edge along axis ``a``, the four cells sharing it, in counter-
# clockwise order seen from ``+a``.  The perpendicular axes are taken as
# ``(a+1) % 3`` and ``(a+2) % 3`` to preserve handedness, and the cyclic
# order of the four quadrants around the edge is then always
# ``(0,-1), (0,0), (-1,0), (-1,-1)`` in those two axes.
_QUAD_CELL_OFFSETS = (
    ((0, 0, -1), (0, 0, 0), (0, -1, 0), (0, -1, -1)),   # x edges
    ((-1, 0, 0), (0, 0, 0), (0, 0, -1), (-1, 0, -1)),   # y edges
    ((0, -1, 0), (0, 0, 0), (-1, 0, 0), (-1, -1, 0)),   # z edges
)


def _resolve_region(node, bounds, resolution, padding):
    """Return ``(origin, spacing, counts)`` for a lattice of cubic cells.

    Cells are kept cubic because the QEF solve weights all three axes
    equally; anisotropic cells would bias vertex placement.
    """
    region = node.bounds if bounds is None else bounds
    if bounds_is_empty(region):
        raise SdfError("dual_contour: the region to mesh is empty")
    lo, hi = region
    for i in range(3):
        if lo[i] in (INF, -INF) or hi[i] in (INF, -INF):
            raise SdfError(
                "dual_contour: this field is unbounded; supply explicit "
                "bounds (intersect it with a solid, or pass bounds=)"
            )
    lo = np.asarray(lo, dtype=float)
    hi = np.asarray(hi, dtype=float)
    extent = hi - lo
    if not np.all(extent >= 0.0):
        raise SdfError("dual_contour: bounds are inverted")

    resolution = int(resolution)
    if resolution < 2:
        raise SdfError("dual_contour: resolution must be at least 2")

    longest = float(np.max(extent))
    if longest <= 0.0:
        raise SdfError("dual_contour: the region to mesh has no extent")
    spacing = longest / resolution

    # Pad outward so the surface never sits exactly on the lattice boundary,
    # where a cell would be clipped and the mesh left open.
    pad = 2.0 * spacing if padding is None else float(padding)
    if pad < 0.0:
        raise SdfError("dual_contour: padding must be non-negative")

    counts = np.maximum(
        np.ceil((extent + 2.0 * pad) / spacing).astype(int), 2
    )
    # Recentre so the padded lattice straddles the region symmetrically.
    origin = 0.5 * (lo + hi) - 0.5 * counts * spacing
    return origin, spacing, tuple(int(c) for c in counts)


def _chunk_spans(n, size):
    """Split ``range(n)`` into contiguous spans of at most ``size``."""
    return [(start, min(start + size, n)) for start in range(0, n, size)]


def sample_corners(node, origin, spacing, counts, chunk=32, prune=True,
                   backend=DEFAULT_BACKEND):
    """Fill the ``(nx+1, ny+1, nz+1)`` corner lattice with field values.

    Evaluates in chunks so that peak memory stays proportional to the chunk
    rather than to the whole lattice, and skips chunks the Lipschitz bound
    proves the surface cannot reach.
    """
    shape = tuple(c + 1 for c in counts)
    values = np.empty(shape, dtype=float)

    spans = [_chunk_spans(shape[a], max(int(chunk), 2)) for a in range(3)]
    blocks = [
        (sx, sy, sz) for sx in spans[0] for sy in spans[1] for sz in spans[2]
    ]

    # A chunk is safe to skip only if no edge *touching* it can cross, so the
    # margin has to cover one cell beyond the chunk in every direction.
    margin = spacing * math.sqrt(3.0) * 2.0
    lipschitz = node.lipschitz

    if prune:
        centres = np.empty((len(blocks), 3), dtype=float)
        radii = np.empty(len(blocks), dtype=float)
        for n, (sx, sy, sz) in enumerate(blocks):
            starts = np.array([sx[0], sy[0], sz[0]], dtype=float)
            stops = np.array([sx[1] - 1, sy[1] - 1, sz[1] - 1], dtype=float)
            low = origin + starts * spacing
            high = origin + stops * spacing
            centres[n] = 0.5 * (low + high)
            radii[n] = 0.5 * float(np.linalg.norm(high - low))
        centre_values = evaluate(node, centres, backend=backend)
        far = np.abs(centre_values) > lipschitz * (radii + margin)
    else:
        centre_values = None
        far = np.zeros(len(blocks), dtype=bool)

    for n, (sx, sy, sz) in enumerate(blocks):
        region = (slice(*sx), slice(*sy), slice(*sz))
        if far[n]:
            # Sign-correct by construction: no surface within `margin`, so
            # the field cannot change sign anywhere in this chunk.
            values[region] = centre_values[n]
            continue
        axes = [
            origin[a] + np.arange(span[0], span[1], dtype=float) * spacing
            for a, span in enumerate((sx, sy, sz))
        ]
        grid = np.meshgrid(*axes, indexing="ij")
        pts = np.stack([g.ravel() for g in grid], axis=1)
        values[region] = evaluate(node, pts, backend=backend).reshape(
            tuple(s[1] - s[0] for s in (sx, sy, sz))
        )
    return values


def _edge_crossings(values, axis):
    """Boolean array of sign-changing edges along ``axis``.

    ``f < 0`` is inside; a sample of exactly zero counts as outside, which
    keeps the inside test consistent and avoids double-counting a corner the
    surface passes exactly through.
    """
    inside = values < 0.0
    lower = [slice(None)] * 3
    upper = [slice(None)] * 3
    lower[axis] = slice(0, -1)
    upper[axis] = slice(1, None)
    return inside[tuple(lower)] != inside[tuple(upper)]


def _crossing_points(values, axis, idx, origin, spacing):
    """Linear-interpolated surface crossings on the given edges.

    Along one cell of an exact distance field the function is very nearly
    linear, so plain interpolation is accurate and — unlike a bisection loop
    with a tolerance — exactly reproducible.
    """
    lower = [slice(None)] * 3
    upper = [slice(None)] * 3
    lower[axis] = slice(0, -1)
    upper[axis] = slice(1, None)
    f0 = values[tuple(lower)][idx]
    f1 = values[tuple(upper)][idx]

    denom = f0 - f1
    # A sign change guarantees f0 != f1, but guard the division anyway so a
    # denormal cannot produce a nan that silently corrupts the mesh.
    t = np.where(np.abs(denom) > 0.0, f0 / np.where(denom == 0.0, 1.0, denom),
                 0.5)
    t = np.clip(t, 0.0, 1.0)

    base = np.stack(idx, axis=1).astype(float)
    points = origin + base * spacing
    points[:, axis] += t * spacing
    return points


def _incident_cells(idx, offset, counts):
    """Cell indices adjacent to edges ``idx`` at ``offset``, and a validity
    mask for those that fall inside the lattice."""
    ci = idx[0] + offset[0]
    cj = idx[1] + offset[1]
    ck = idx[2] + offset[2]
    ok = (
        (ci >= 0) & (ci < counts[0])
        & (cj >= 0) & (cj < counts[1])
        & (ck >= 0) & (ck < counts[2])
    )
    return ci, cj, ck, ok


def _edge_normals(node, points, axis, eps, backend):
    """Unit surface normals at the crossings, from the field gradient.

    Where the gradient vanishes — a crossing that lands on the medial axis —
    fall back to the edge direction, which is at least a plane the surface
    genuinely crosses and keeps the QEF matrix non-singular.
    """
    grad = gradient(node, points, eps=eps, backend=backend)
    lengths = np.linalg.norm(grad, axis=1)
    usable = lengths > 0.0
    normals = np.zeros_like(grad)
    normals[usable] = grad[usable] / lengths[usable, None]
    if not np.all(usable):
        fallback = np.zeros(3)
        fallback[axis] = 1.0
        normals[~usable] = fallback
    return normals


def dual_contour(node, resolution=64, bounds=None, padding=None,
                 chunk=32, prune=True, backend=DEFAULT_BACKEND):
    """Mesh the surface of ``node`` by dual contouring.

    :param resolution: cells across the longest axis of the region.
    :param bounds: region to mesh; defaults to the node's own bounds, which
        must then be finite.
    :param padding: outward margin; defaults to two cells, so the surface
        never lands on the lattice boundary.
    :param chunk: corner-lattice chunk edge, trading memory for call count.
    :param prune: use the Lipschitz bound to skip empty chunks.
    :returns: a :class:`Mesh`.

    The result is deterministic: identical arguments give identical arrays.
    """
    origin, spacing, counts = _resolve_region(node, bounds, resolution,
                                              padding)
    values = sample_corners(node, origin, spacing, counts, chunk=chunk,
                            prune=prune, backend=backend)

    crossings = [_edge_crossings(values, axis) for axis in range(3)]
    edge_idx = [np.nonzero(crossings[axis]) for axis in range(3)]

    # Which cells the surface passes through.
    active = np.zeros(counts, dtype=bool)
    for axis in range(3):
        for offset in _QUAD_CELL_OFFSETS[axis]:
            ci, cj, ck, ok = _incident_cells(edge_idx[axis], offset, counts)
            active[ci[ok], cj[ok], ck[ok]] = True

    active_idx = np.nonzero(active)
    n_active = int(active_idx[0].size)
    if n_active == 0:
        empty = np.zeros((0, 3), dtype=float)
        return Mesh(empty, empty, np.zeros((0, 3), dtype=np.int64))

    # C-order of the mask, so vertex numbering is reproducible.
    cell_index = np.full(counts, -1, dtype=np.int64)
    cell_index[active_idx] = np.arange(n_active, dtype=np.int64)

    ata = np.zeros((n_active, 9), dtype=float)
    atb = np.zeros((n_active, 3), dtype=float)
    mass = np.zeros((n_active, 3), dtype=float)
    tally = np.zeros(n_active, dtype=float)
    normal_sum = np.zeros((n_active, 3), dtype=float)

    eps = 1e-4 * spacing
    crossing_points = []
    for axis in range(3):
        idx = edge_idx[axis]
        if idx[0].size == 0:
            crossing_points.append(None)
            continue
        pts = _crossing_points(values, axis, idx, origin, spacing)
        nrm = _edge_normals(node, pts, axis, eps, backend)
        crossing_points.append(pts)

        outer = (nrm[:, :, None] * nrm[:, None, :]).reshape(-1, 9)
        plane_d = np.sum(nrm * pts, axis=1)
        nb = nrm * plane_d[:, None]

        for offset in _QUAD_CELL_OFFSETS[axis]:
            ci, cj, ck, ok = _incident_cells(idx, offset, counts)
            cells = cell_index[ci[ok], cj[ok], ck[ok]]
            for comp in range(9):
                ata[:, comp] += np.bincount(
                    cells, weights=outer[ok, comp], minlength=n_active
                )
            for comp in range(3):
                atb[:, comp] += np.bincount(
                    cells, weights=nb[ok, comp], minlength=n_active
                )
                mass[:, comp] += np.bincount(
                    cells, weights=pts[ok, comp], minlength=n_active
                )
                normal_sum[:, comp] += np.bincount(
                    cells, weights=nrm[ok, comp], minlength=n_active
                )
            tally += np.bincount(cells, minlength=n_active)

    # Solve each cell's quadratic error function about its mass point, using
    # a pseudo-inverse with singular-value truncation rather than uniform
    # damping.  ATA is symmetric positive semi-definite, and its rank is the
    # number of directions the surface actually constrains: one inside a
    # flat cell, two along an edge, three at a corner.  Keeping only the
    # well-conditioned directions solves exactly for those and leaves the
    # vertex at the mass point along the rest.  Uniform Tikhonov damping
    # instead lets an ill-conditioned flat cell slide its vertex far along
    # the surface, which is what drives vertices out of their cells and
    # produces coincident points where two neighbours both clamp to the
    # corner they share.
    matrices = ata.reshape(n_active, 3, 3)
    centre = mass / tally[:, None]
    residual = atb - np.einsum("nij,nj->ni", matrices, centre)

    evals, evecs = np.linalg.eigh(matrices)
    largest = evals[:, 2:3]
    keep = evals > QEF_TRUNCATION * largest
    usable = keep & (evals > 0.0)
    inverted = np.where(
        usable, 1.0 / np.where(evals > 0.0, evals, 1.0), 0.0
    )
    projected = np.einsum("nik,ni->nk", evecs, residual) * inverted
    vertices = centre + np.einsum("nik,nk->ni", evecs, projected)
    # A cell with no usable direction keeps its mass point.
    vertices = np.where(largest > 0.0, vertices, centre)

    # Keep each vertex in its own cell; an unclamped QEF can place a vertex
    # far outside on a nearly-degenerate configuration and tangle the mesh.
    cell_lo = origin + np.stack(active_idx, axis=1).astype(float) * spacing
    vertices = np.clip(vertices, cell_lo, cell_lo + spacing)
    vertices = _separate_coincident(vertices, cell_lo, spacing)

    normals = _vertex_normals(node, vertices, normal_sum, eps, backend)
    triangles = _emit_triangles(values, edge_idx, cell_index, counts)
    return Mesh(vertices, normals, triangles)


def _separate_coincident(vertices, cell_lo, spacing):
    """Pull any vertices sharing a position strictly inside their own cells.

    Clamping to the closed cell lets two neighbours land on the face they
    share.  That matters more than it looks: yapCAD's
    :func:`yapcad.geom3d.issolidclosed` keys edges by vertex *position*
    rather than index, so a coincident pair merges two distinct edges, the
    solid reads as open, and ``volumeof`` refuses it — even though the index
    topology is perfectly sound.

    Only the colliding vertices are moved, and only as far as one position
    quantum needs.  Everything else keeps the position the QEF chose, so a
    box whose faces land exactly on lattice planes still meshes to its exact
    volume.  After the move each colliding vertex sits at least ``inset``
    inside its own cell, so it is at least ``inset`` from anything in a
    neighbouring cell and at least ``2 * inset`` from another moved vertex —
    both comfortably more than the quantum.
    """
    keys = np.round(np.asarray(vertices) / epsilon).astype(np.int64)
    _, inverse, counts = np.unique(
        keys, axis=0, return_inverse=True, return_counts=True
    )
    colliding = counts[inverse.ravel()] > 1
    if not np.any(colliding):
        return vertices

    inset = min(4.0 * epsilon, 0.25 * spacing)
    moved = vertices.copy()
    moved[colliding] = np.clip(
        vertices[colliding],
        cell_lo[colliding] + inset,
        cell_lo[colliding] + spacing - inset,
    )
    return moved


def _vertex_normals(node, vertices, normal_sum, eps, backend):
    """Field gradient at each vertex, falling back to the accumulated
    edge normals where the gradient vanishes."""
    grad = gradient(node, vertices, eps=eps, backend=backend)
    lengths = np.linalg.norm(grad, axis=1)
    usable = lengths > 0.0
    normals = np.zeros_like(grad)
    normals[usable] = grad[usable] / lengths[usable, None]
    if not np.all(usable):
        spare = normal_sum[~usable]
        spare_len = np.linalg.norm(spare, axis=1)
        good = spare_len > 0.0
        replacement = np.tile(np.array([0.0, 0.0, 1.0]), (spare.shape[0], 1))
        replacement[good] = spare[good] / spare_len[good, None]
        normals[~usable] = replacement
    return normals


def _emit_triangles(values, edge_idx, cell_index, counts):
    """One quad per interior crossing edge, split into two triangles.

    Winding is taken from the direction the field increases along the edge,
    so every triangle comes out counter-clockwise seen from outside the
    solid — yapCAD's convention, and what makes ``volumeof`` positive.
    """
    pieces = []
    for axis in range(3):
        idx = edge_idx[axis]
        if idx[0].size == 0:
            continue

        corners = []
        valid = np.ones(idx[0].size, dtype=bool)
        for offset in _QUAD_CELL_OFFSETS[axis]:
            ci, cj, ck, ok = _incident_cells(idx, offset, counts)
            valid &= ok
            corners.append((ci, cj, ck))
        if not np.any(valid):
            continue

        quad = np.stack(
            [cell_index[c[0][valid], c[1][valid], c[2][valid]]
             for c in corners],
            axis=1,
        )

        lower = [slice(None)] * 3
        lower[axis] = slice(0, -1)
        inside_first = (values[tuple(lower)][idx] < 0.0)[valid]
        flip = ~inside_first
        quad[flip] = quad[flip][:, ::-1]

        pieces.append(quad[:, [0, 1, 2]])
        pieces.append(quad[:, [0, 2, 3]])

    if not pieces:
        return np.zeros((0, 3), dtype=np.int64)
    return np.concatenate(pieces, axis=0)


def manifold_defects(mesh):
    """Count the ways ``mesh`` fails to be a closed manifold.

    :returns: ``{"boundary": n, "nonmanifold": n, "coincident": n}`` — edges
        used by exactly one triangle, edges used by three or more, and
        vertices sharing a position with another vertex.

    Dual contouring places **one** vertex per cell, so a cell that the
    surface passes through twice cannot represent both sheets and the result
    is non-manifold there.  In practice that means a wall thinner than a
    cell, or two lattice sheets approaching within a cell of each other.
    The cure is more resolution, or Manifold Dual Contouring (Schaefer and
    Ju), which splits such a cell into one vertex per surface component and
    is not implemented here.

    Coincident vertices are counted separately because yapCAD's own
    :func:`yapcad.geom3d.issolidclosed` keys edges by vertex *position*
    rather than index, so two vertices at the same point merge two distinct
    edges and the solid reads as open even when the index topology is sound.
    ``volumeof`` then refuses it. Index manifoldness alone is therefore not
    enough to promise the downstream a usable solid.

    Checking is cheap next to meshing, and design document §5.3 argues that
    a tool which states its limits beats one that degrades silently — so
    :func:`~yapcad.sdf.convert.to_solid` runs this by default.
    """
    triangles = np.asarray(mesh.triangles)
    if triangles.size == 0:
        return {"boundary": 0, "nonmanifold": 0, "coincident": 0}
    edges = np.concatenate([
        triangles[:, [0, 1]], triangles[:, [1, 2]], triangles[:, [2, 0]]
    ])
    _, counts = np.unique(np.sort(edges, axis=1), axis=0, return_counts=True)

    # Quantised the same way geom3d._point_to_key quantises, so this agrees
    # with what issolidclosed will conclude.
    keys = np.round(np.asarray(mesh.vertices) / epsilon).astype(np.int64)
    distinct = len(np.unique(keys, axis=0))
    return {
        "boundary": int(np.sum(counts == 1)),
        "nonmanifold": int(np.sum(counts > 2)),
        "coincident": int(len(mesh.vertices) - distinct),
    }


def is_manifold(mesh):
    """True when ``mesh`` is a closed manifold with no coincident vertices."""
    return not any(manifold_defects(mesh).values())
