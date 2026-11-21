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
from yapcad.geometry_utils import DEGENERATE_ASPECT, degenerate_faces
from yapcad.sdf.evaluate import (
    DEFAULT_BACKEND,
    evaluate,
    evaluate_batches,
    gradient,
    thread_count,
)
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


def _resolve_region(node, bounds, resolution, padding, cell_mm=None):
    """Return ``(origin, spacing, counts)`` for a lattice of cubic cells.

    The cell comes either from ``resolution`` -- cells across the longest
    axis -- or straight from ``cell_mm``, which takes precedence.  Cells are
    kept cubic because the QEF solve weights all three axes equally;
    anisotropic cells would bias vertex placement.
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

    longest = float(np.max(extent))
    if longest <= 0.0:
        raise SdfError("dual_contour: the region to mesh has no extent")
    if cell_mm is None:
        resolution = int(resolution)
        if resolution < 2:
            raise SdfError("dual_contour: resolution must be at least 2")
        spacing = longest / resolution
    else:
        spacing = float(cell_mm)
        if not math.isfinite(spacing) or spacing <= 0.0:
            raise SdfError(
                "dual_contour: cell_mm must be positive and finite"
            )

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

    def corners(block):
        axes = [
            origin[a] + np.arange(span[0], span[1], dtype=float) * spacing
            for a, span in enumerate(block)
        ]
        grid = np.meshgrid(*axes, indexing="ij")
        return np.stack([g.ravel() for g in grid], axis=1)

    near = []
    for n, (sx, sy, sz) in enumerate(blocks):
        if far[n]:
            # Sign-correct by construction: no surface within `margin`, so
            # the field cannot change sign anywhere in this chunk.
            values[slice(*sx), slice(*sy), slice(*sz)] = centre_values[n]
        else:
            near.append(n)

    # The chunks that may hold surface are independent of one another, so
    # they are evaluated across threads -- a wave of them at a time, which
    # keeps peak memory proportional to a few chunks rather than to all.
    wave = max(1, 2 * thread_count())
    for start in range(0, len(near), wave):
        group = [blocks[n] for n in near[start:start + wave]]
        sampled = evaluate_batches(node, [corners(block) for block in group],
                                   backend=backend)
        for (sx, sy, sz), block_values in zip(group, sampled):
            values[slice(*sx), slice(*sy), slice(*sz)] = block_values.reshape(
                tuple(s[1] - s[0] for s in (sx, sy, sz)))
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


# ---------------------------------------------------------------------------
# Manifold dual contouring: one vertex per surface component of a cell
# ---------------------------------------------------------------------------
#
# Plain dual contouring places one vertex per cell.  Where the surface
# passes through a cell twice -- a thin wall, or the thin wedge of material
# near any acute edge, such as a helical gear's tooth tip -- that one vertex
# serves two sheets and the mesh comes out non-manifold, at any resolution.
# Following Schaefer, Ju and Warren's Manifold Dual Contouring, each cell
# instead gets a vertex per connected component of the surface within it.
#
# A cell's sign configuration fixes those components.  Two crossing edges of
# the cell belong to one component when they bound one surface segment on a
# face the cell shares: a face with two crossings has one segment joining
# them; a face with four -- inside corners on a diagonal -- has two, each
# cutting off one inside corner.  That rule depends only on the face, so the
# two cells sharing a face always agree, which is what keeps the result
# manifold.  A cell with a single component gets exactly the vertex plain
# dual contouring would give it.

#: Local cell edges as (axis, u, v): u and v are the edge's offsets along the
#: axes (axis + 1) % 3 and (axis + 2) % 3.  Index = axis * 4 + u * 2 + v.
_CELL_EDGES = [(a, u, v) for a in range(3) for u in (0, 1) for v in (0, 1)]


def _edge_corners(edge):
    """The two corners of a local edge, as (x, y, z) offsets."""
    a, u, v = edge
    ends = []
    for t in (0, 1):
        c = [0, 0, 0]
        c[a] = t
        c[(a + 1) % 3] = u
        c[(a + 2) % 3] = v
        ends.append(tuple(c))
    return ends


def _corner_bit(c):
    return c[0] | (c[1] << 1) | (c[2] << 2)


def _component_table():
    """For each of the 256 sign configurations, the component of each of
    the 12 cell edges (-1 where the edge has no crossing), and the count."""
    edge_ends = [_edge_corners(e) for e in _CELL_EDGES]
    table = np.full((256, 12), -1, dtype=np.int64)
    count = np.zeros(256, dtype=np.int64)
    for config in range(256):
        def inside(c):
            return bool((config >> _corner_bit(c)) & 1)
        crossing = [inside(a) != inside(b) for a, b in edge_ends]
        parent = list(range(12))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for normal in range(3):
            for side in (0, 1):
                on_face = [e for e, ends in enumerate(edge_ends)
                           if _CELL_EDGES[e][0] != normal
                           and ends[0][normal] == side and crossing[e]]
                if len(on_face) == 2:
                    a, b = on_face
                    parent[find(a)] = find(b)
                elif len(on_face) == 4:
                    # Pair the two edges around each inside corner.
                    corners = {c for e in on_face for c in edge_ends[e]}
                    for c in sorted(corners):
                        if not inside(c):
                            continue
                        pair = [e for e in on_face if c in edge_ends[e]]
                        parent[find(pair[0])] = find(pair[1])
        roots = {}
        for e in range(12):
            if crossing[e]:
                root = find(e)
                table[config, e] = roots.setdefault(root, len(roots))
        count[config] = len(roots)
    return table, count


_COMPONENT, _COMPONENT_COUNT = _component_table()


def _cell_configs(values, cells):
    """Sign configuration (bit per corner, set when inside) of each cell."""
    ci, cj, ck = cells
    config = np.zeros(ci.size, dtype=np.int64)
    for c in range(8):
        dx, dy, dz = c & 1, (c >> 1) & 1, (c >> 2) & 1
        config |= (values[ci + dx, cj + dy, ck + dz] < 0.0).astype(
            np.int64) << c
    return config


def _local_edge(axis, offset):
    """Local index, within an incident cell, of a grid edge along ``axis``
    reached from that cell at ``offset`` (each component 0 or -1)."""
    u = -offset[(axis + 1) % 3]
    v = -offset[(axis + 2) % 3]
    return axis * 4 + u * 2 + v


def proves_empty(node, resolution=64, bounds=None, padding=None, chunk=32,
                 backend=DEFAULT_BACKEND, cell_mm=None):
    """True when this lattice *proves* that nothing in the region is inside
    ``node``, so an empty mesh of it is the right answer and not a miss.

    It is the chunk-pruning argument applied to the whole lattice: a field
    that never changes faster than its Lipschitz constant ``L`` cannot fall
    to zero within ``r`` of a sample whose value exceeds ``L * r``, so if
    every corner sample does, with ``r`` the half-diagonal of a cell, there
    is no interior anywhere between the samples either.

    ``False`` means *not proved*, not *occupied*: the surface may be there
    and too thin to show on this lattice, which is exactly the case
    :func:`yapcad.sdf.to_solid` responds to by meshing finer.

    Pruning is off here on purpose.  A pruned chunk is filled with its
    centre value, which can read higher than the values it stands for, and
    this proof needs a true minimum.
    """
    lipschitz = node.lipschitz
    if not math.isfinite(lipschitz) or lipschitz <= 0.0:
        return False
    origin, spacing, counts = _resolve_region(node, bounds, resolution,
                                              padding, cell_mm)
    values = sample_corners(node, origin, spacing, counts, chunk=chunk,
                            prune=False, backend=backend)
    return float(values.min()) > lipschitz * 0.5 * math.sqrt(3.0) * spacing


def dual_contour(node, resolution=64, bounds=None, padding=None,
                 chunk=32, prune=True, backend=DEFAULT_BACKEND,
                 cell_mm=None):
    """Mesh the surface of ``node`` by dual contouring.

    :param resolution: cells across the longest axis of the region.
    :param bounds: region to mesh; defaults to the node's own bounds, which
        must then be finite.
    :param padding: outward margin; defaults to two cells, so the surface
        never lands on the lattice boundary.
    :param chunk: corner-lattice chunk edge, trading memory for call count.
    :param prune: use the Lipschitz bound to skip empty chunks.
    :param cell_mm: the cell size to mesh at, in model units, instead of
        deriving one from ``resolution``.
    :returns: a :class:`Mesh`.

    An empty :class:`Mesh` means the lattice found no surface, which is
    either empty space or a part thinner than a cell; :func:`proves_empty`
    is the test that tells those apart, and
    :func:`yapcad.sdf.to_solid` uses it to decide whether to mesh finer.

    The result is deterministic: identical arguments give identical arrays.
    """
    origin, spacing, counts = _resolve_region(node, bounds, resolution,
                                              padding, cell_mm)
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

    # C-order of the mask, then component order within a cell, so vertex
    # numbering is reproducible.  A cell holds one vertex per surface
    # component; almost every cell has exactly one.
    config = _cell_configs(values, active_idx)
    components = _COMPONENT_COUNT[config]
    first_vertex = np.concatenate([[0], np.cumsum(components)[:-1]])
    n_vertices = int(components.sum())
    cell_of_vertex = np.repeat(np.arange(n_active), components)
    cell_index = np.full(counts, -1, dtype=np.int64)
    cell_index[active_idx] = np.arange(n_active, dtype=np.int64)
    cell_config = np.zeros(counts, dtype=np.int64)
    cell_config[active_idx] = config
    cell_first = np.zeros(counts, dtype=np.int64)
    cell_first[active_idx] = first_vertex

    def vertex_of(ci, cj, ck, axis, offset):
        """The vertex, in cells (ci, cj, ck), of the component containing
        the grid edge they reach at ``offset``."""
        local = _local_edge(axis, offset)
        comp = _COMPONENT[cell_config[ci, cj, ck], local]
        return cell_first[ci, cj, ck] + comp

    ata = np.zeros((n_vertices, 9), dtype=float)
    atb = np.zeros((n_vertices, 3), dtype=float)
    mass = np.zeros((n_vertices, 3), dtype=float)
    tally = np.zeros(n_vertices, dtype=float)
    normal_sum = np.zeros((n_vertices, 3), dtype=float)

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
            cells = vertex_of(ci[ok], cj[ok], ck[ok], axis, offset)
            for comp in range(9):
                ata[:, comp] += np.bincount(
                    cells, weights=outer[ok, comp], minlength=n_vertices
                )
            for comp in range(3):
                atb[:, comp] += np.bincount(
                    cells, weights=nb[ok, comp], minlength=n_vertices
                )
                mass[:, comp] += np.bincount(
                    cells, weights=pts[ok, comp], minlength=n_vertices
                )
                normal_sum[:, comp] += np.bincount(
                    cells, weights=nrm[ok, comp], minlength=n_vertices
                )
            tally += np.bincount(cells, minlength=n_vertices)

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
    matrices = ata.reshape(n_vertices, 3, 3)
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
    cell_lo = origin + np.stack(active_idx, axis=1).astype(float)[
        cell_of_vertex] * spacing
    vertices = np.clip(vertices, cell_lo, cell_lo + spacing)

    triangles = _emit_triangles(values, edge_idx, vertex_of, counts, vertices)
    vertices, triangles, kept = _weld_coincident(
        vertices, triangles, cell_lo, centre, spacing)
    triangles = _flip_degenerate(vertices, triangles)
    normals = _vertex_normals(node, vertices, normal_sum[kept], eps, backend)
    return Mesh(vertices, normals, triangles)


def _position_keys(vertices):
    """Vertex positions quantised the way :func:`yapcad.geom3d._point_to_key`
    quantises them, so a collision here is one ``issolidclosed`` would see."""
    return np.round(np.asarray(vertices) / epsilon).astype(np.int64)


def _weld_coincident(vertices, triangles, cell_lo, centre, spacing):
    """Merge vertices that share a position, or move them apart.

    Clamping each vertex to its closed cell lets two neighbouring cells put
    their vertices at the same point on the face, edge or corner they share.
    Two vertices at one position are a defect for every consumer that keys
    the mesh by position -- :func:`yapcad.geom3d.issolidclosed`, and any STL
    reader, which has nothing else to go on.

    Where the two vertices are joined by a mesh edge and collapsing that
    edge keeps the surface a manifold (the *link condition*: the only
    vertices they both neighbour are the two opposite that edge), they are
    merged and the two triangles on the edge, which have no area, are
    dropped.  Nothing moves, so the surface is exactly the one the QEF
    chose.

    What cannot be merged -- two surface components of one cell that meet at
    a point, or two sheets that touch without sharing an edge -- is moved
    apart by one position quantum instead, as it always was.  Merging is
    what fixes #121: the 40 nm edges that nudge used to leave between
    neighbouring cells collapsed when the mesh was stored as float32.

    :returns: ``(vertices, triangles, kept)``, where ``kept`` indexes the
        surviving vertices in the input, so per-vertex arrays can follow.
    """
    n = len(vertices)
    keys = _position_keys(vertices)
    _, inverse, counts = np.unique(
        keys, axis=0, return_inverse=True, return_counts=True)
    inverse = inverse.ravel()
    colliding = np.nonzero(counts[inverse] > 1)[0]
    if colliding.size == 0 or triangles.size == 0:
        return vertices, triangles, np.arange(n)

    tris = triangles.copy()
    alive = np.ones(len(tris), dtype=bool)
    is_colliding = np.zeros(n, dtype=bool)
    is_colliding[colliding] = True
    incident = {int(v): set() for v in colliding}
    for t in np.nonzero(is_colliding[tris].any(axis=1))[0]:
        for v in tris[t]:
            if is_colliding[v]:
                incident[int(v)].add(int(t))

    def neighbours(v):
        out = set()
        for t in incident[v]:
            out.update(int(x) for x in tris[t])
        out.discard(v)
        return out

    def collapse(u, v):
        """Merge ``v`` into ``u`` if the edge between them allows it."""
        shared = incident[u] & incident[v]
        if len(shared) != 2:
            return False
        opposite = set()
        for t in shared:
            opposite.update(int(x) for x in tris[t])
        opposite -= {u, v}
        if len(opposite) != 2 or neighbours(u) & neighbours(v) != opposite:
            return False
        for t in shared:
            alive[t] = False
            for w in tris[t]:
                if is_colliding[w] and int(w) not in (u, v):
                    incident[int(w)].discard(t)
        for t in incident[v] - shared:
            tris[t][tris[t] == v] = u
        incident[u] = (incident[u] | incident[v]) - shared
        incident[v] = set()
        return True

    removed = np.zeros(n, dtype=bool)
    stranded = []
    order = colliding[np.argsort(inverse[colliding], kind="stable")]
    groups = np.split(order, np.nonzero(np.diff(inverse[order]))[0] + 1)
    for group in groups:
        # Several cells can put a vertex on one lattice point -- a face that
        # lands on a lattice plane does it at every corner -- and a vertex
        # may be joined by an edge to only some of the others.  Merge pairs
        # until no collapse is allowed, in a fixed order for determinism.
        survivors = [int(x) for x in group]
        merged = True
        while merged and len(survivors) > 1:
            merged = False
            for a in range(len(survivors)):
                for b in range(a + 1, len(survivors)):
                    u, v = survivors[a], survivors[b]
                    if collapse(u, v):
                        removed[v] = True
                        survivors.pop(b)
                        merged = True
                        break
                if merged:
                    break
        if len(survivors) > 1:
            stranded.extend(survivors)

    vertices = vertices.copy()
    if stranded:
        # Moved only as far as one position quantum needs, so a face that
        # lands on a lattice plane keeps its exact position.  This distance
        # does not survive float32 far from the origin, which is why every
        # vertex that can be merged is merged instead: what reaches here is
        # rare -- two sheets that touch without sharing an edge, such as the
        # zero-thickness wall where two unioned parts share a face.
        s = np.asarray(stranded)
        inset = min(4.0 * epsilon, 0.25 * spacing)
        toward = centre[s] - vertices[s]
        length = np.linalg.norm(toward, axis=1, keepdims=True)
        step = np.where(length > 0.0,
                        toward * (inset / np.where(length > 0.0, length, 1.0)),
                        0.0)
        # The step towards each vertex's own mass point separates two
        # components of one cell; the inset separates neighbouring cells.
        vertices[s] = np.clip(vertices[s] + step,
                              cell_lo[s] + inset, cell_lo[s] + spacing - inset)

    kept = np.nonzero(~removed)[0]
    remap = np.full(n, -1, dtype=np.int64)
    remap[kept] = np.arange(kept.size)
    return vertices[kept], remap[tris[alive]], kept


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


def _shape_quality(a, b, c):
    """Twice a triangle's area over its longest edge squared, per row:
    dimensionless, 0 for a degenerate triangle and about 0.87 at best."""
    cross = np.linalg.norm(np.cross(b - a, c - a), axis=1)
    longest = np.maximum.reduce([
        np.einsum("ij,ij->i", b - a, b - a),
        np.einsum("ij,ij->i", c - b, c - b),
        np.einsum("ij,ij->i", a - c, a - c),
    ])
    return np.where(longest > 0.0, cross / np.where(longest > 0.0, longest, 1.0), 0.0)


#: Shape quality below which a triangle has no meaningful normal.
_DEGENERATE = DEGENERATE_ASPECT


def _flip_degenerate(vertices, triangles, rounds=8):
    """Remove zero-area triangles by flipping their longest edge.

    A triangle whose vertices lie on one line -- four vertices of a quad
    clamped onto a lattice line along a sharp edge, where neither diagonal
    helps -- has a third vertex on its longest edge.  Flipping that edge
    with the triangle across it splits the neighbour at that vertex: two
    triangles covering exactly the neighbour's area, so the surface does not
    move and nothing folds.  A flip is skipped when the new edge already
    exists, and the pass repeats for chains of such triangles.
    """
    tris = np.asarray(triangles).copy()
    if tris.size == 0:
        return tris
    for _ in range(rounds):
        P = vertices[tris]
        quality = _shape_quality(P[:, 0], P[:, 1], P[:, 2])
        bad = np.nonzero(quality < _DEGENERATE)[0]
        if bad.size == 0:
            break
        edges = {}
        for t, (a, b, c) in enumerate(tris.tolist()):
            edges[(a, b)] = t
            edges[(b, c)] = t
            edges[(c, a)] = t
        done = set()
        flipped = 0
        for t in bad.tolist():
            if t in done:
                continue
            tri = tris[t].tolist()
            lengths = [np.sum((vertices[tri[(k + 1) % 3]] - vertices[tri[k]]) ** 2)
                       for k in range(3)]
            k = int(np.argmax(lengths))
            a, b, c = tri[k], tri[(k + 1) % 3], tri[(k + 2) % 3]
            t2 = edges.get((b, a))
            if t2 is None or t2 in done:
                continue
            other = tris[t2].tolist()
            d = other[(other.index(b) + 2) % 3]
            if d == c or (c, d) in edges or (d, c) in edges:
                continue
            new1, new2 = [a, d, c], [d, b, c]
            q = _shape_quality(vertices[[a, d]], vertices[[d, b]], vertices[[c, c]])
            if q.min() < _DEGENERATE:
                continue
            tris[t], tris[t2] = new1, new2
            for key in ((a, b), (b, c), (c, a), (b, a), (a, d), (d, b)):
                edges.pop(key, None)
            for e in ((a, d), (d, c), (c, a), (d, b), (b, c), (c, d)):
                edges[e] = t if e in ((a, d), (d, c), (c, a)) else t2
            done.update((t, t2))
            flipped += 1
        if not flipped:
            break
    return tris


def _emit_triangles(values, edge_idx, vertex_of, counts, vertices=None):
    """One quad per interior crossing edge, split into two triangles.

    Winding is taken from the direction the field increases along the edge,
    so every triangle comes out counter-clockwise seen from outside the
    solid — yapCAD's convention, and what makes ``volumeof`` positive.

    Given ``vertices``, each quad is split along whichever diagonal leaves
    its worse triangle the better shaped.  A fixed diagonal makes a
    zero-area triangle whenever three of the quad's vertices are clamped
    onto one line of the lattice, and that triangle is load-bearing: it
    closes the mesh by index while having no normal for anything that
    reads it (#85).
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
            corners.append((ci, cj, ck, offset))
        if not np.any(valid):
            continue

        # Each quad corner is the vertex of the component, in that cell,
        # that this edge's crossing belongs to.
        quad = np.stack(
            [vertex_of(c[0][valid], c[1][valid], c[2][valid], axis, c[3])
             for c in corners],
            axis=1,
        )

        lower = [slice(None)] * 3
        lower[axis] = slice(0, -1)
        inside_first = (values[tuple(lower)][idx] < 0.0)[valid]
        flip = ~inside_first
        quad[flip] = quad[flip][:, ::-1]

        if vertices is None:
            other = np.zeros(len(quad), dtype=bool)
        else:
            p = [vertices[quad[:, k]] for k in range(4)]
            first = np.minimum(_shape_quality(p[0], p[1], p[2]),
                               _shape_quality(p[0], p[2], p[3]))
            second = np.minimum(_shape_quality(p[0], p[1], p[3]),
                                _shape_quality(p[1], p[2], p[3]))
            other = second > first
        keep = quad[~other]
        swap = quad[other]
        pieces.extend([keep[:, [0, 1, 2]], keep[:, [0, 2, 3]],
                       swap[:, [0, 1, 3]], swap[:, [1, 2, 3]]])

    if not pieces:
        return np.zeros((0, 3), dtype=np.int64)
    return np.concatenate(pieces, axis=0)


def manifold_defects(mesh):
    """Count the ways ``mesh`` fails to be a closed manifold.

    :returns: ``{"boundary": n, "nonmanifold": n, "coincident": n,
        "degenerate": n}`` — edges used by exactly one triangle, edges used
        by three or more, vertices sharing a position with another vertex,
        and triangles with no meaningful normal
        (:func:`yapcad.geometry_utils.triangle_normal`).

    Dual contouring places **one** vertex per cell, so a cell that the
    surface passes through twice cannot represent both sheets and the result
    is non-manifold there.  In practice that means a wall thinner than a
    cell, or two lattice sheets approaching within a cell of each other.
    :func:`dual_contour` already places one vertex per surface *component*
    of a cell (the vertex-splitting half of Manifold Dual Contouring), so
    what remains is a cell whose components the corner signs cannot tell
    apart; the cure is more resolution, or an adaptive grid.

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
        return {"boundary": 0, "nonmanifold": 0, "coincident": 0,
                "degenerate": 0}
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
        "degenerate": int(np.sum(
            degenerate_faces(np.asarray(mesh.vertices), triangles))),
    }


def is_manifold(mesh):
    """True when ``mesh`` is a closed manifold with no coincident vertices."""
    return not any(manifold_defects(mesh).values())
