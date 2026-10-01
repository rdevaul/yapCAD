"""Mesh booleans by split, classify, select and repair.

This is the core of the native mesh boolean engine.  The previous core
clipped each triangle against the *intersection of half-spaces* of nearby
target faces.  That intersection is always convex, so the result was wrong
wherever the target is concave, and because each triangle was cut against
its own plane set, neighbours' cut points did not coincide and results came
back with true gaps.  It never classified coincident faces at all.

The structure here is the standard one:

1. **Split.**  Each triangle is cut by the planes of the other mesh's
   triangles that actually cross it, and -- where a triangle of the other
   mesh is coplanar with it -- by the lines of that triangle's edges.  Every
   resulting convex fragment then lies entirely inside, outside, or on the
   other surface.  Cutting by an infinite plane over-fragments a little; it
   never misclassifies, because classification does not use the planes.

2. **Classify.**  A fragment lying on a coplanar triangle of the other mesh
   is SAME or OPPOSITE by the sign of the normals' dot product.  Any other
   fragment is IN or OUT by the *generalised winding number* of its interior
   point -- the solid angle the other mesh subtends, over 4 pi.  That is far
   more robust than ray casting: no ray can graze an edge, and it degrades
   gracefully on slightly open input.  Triangles that no triangle of the
   other mesh comes near are classified a connected region at a time.

3. **Select** by the usual table::

       union         A: OUT, SAME      B: OUT
       intersection  A: IN,  SAME      B: IN
       difference    A: OUT, OPPOSITE  B: IN (reversed)

4. **Repair.**  Both meshes' cut points lie on the same intersection lines,
   so after welding, every remaining crack is a vertex lying exactly on
   another triangle's edge: a T-junction, which is fixed exactly by
   splitting that edge.  Nothing is snapped and nothing is dropped.

Fragments are triangulated by clipping only *true* corners, so the vertices
that sit on straight fragment edges -- precisely the ones neighbours share
-- are kept and no zero-area triangle is ever emitted.
"""

import math

import numpy as np

from yapcad.geom import epsilon

IN, OUT, SAME, OPPOSITE = "in", "out", "same", "opposite"

KEEP = {
    "union": ({OUT, SAME}, {OUT}),
    "intersection": ({IN, SAME}, {IN}),
    "difference": ({OUT, OPPOSITE}, {IN}),
}


# ---------------------------------------------------------------------------
# Input
# ---------------------------------------------------------------------------


def _raw_triangles(solid):
    """Every triangle of a solid, unfiltered (see manifold_engine)."""
    tris = []
    for surf in solid[1]:
        verts = surf[1]
        for face in surf[3]:
            if len(face) == 3:
                tris.append((verts[face[0]], verts[face[1]], verts[face[2]]))
    if tris or solid[1]:
        return tris
    from yapcad.mesh import mesh_view
    return [(v0, v1, v2) for _, v0, v1, v2 in mesh_view(solid)]


def weld(points, quantum=epsilon):
    """Merge points sharing a quantised position key; first occurrence wins.

    The default quantum is yapCAD's own, the one ``geom3d.issolidclosed``
    keys vertices by, so a mesh welded with it means by "closed" exactly
    what yapCAD will check.  The algorithm itself welds finer, at its
    geometric tolerance, and only re-welds at this quantum at the end.

    Returns ``(vertices, index)`` with ``vertices[index[i]]`` the welded
    position of ``points[i]``.  Deterministic: vertex order is first-seen.
    """
    pts = np.asarray(points, dtype=np.float64).reshape(-1, 3)
    if len(pts) == 0:
        return np.zeros((0, 3)), np.zeros(0, dtype=np.int64)
    keys = np.round(pts / quantum).astype(np.int64)
    _, first, inverse = np.unique(keys, axis=0, return_index=True,
                                  return_inverse=True)
    inverse = inverse.ravel()
    order = np.argsort(first, kind="stable")
    remap = np.empty_like(order)
    remap[order] = np.arange(len(order))
    return pts[first[order]], remap[inverse]


def solid_to_mesh(solid):
    """An indexed triangle mesh ``(V, F)`` for a yapCAD solid."""
    tris = _raw_triangles(solid)
    if not tris:
        return np.zeros((0, 3)), np.zeros((0, 3), dtype=np.int64)
    V, index = weld([[p[0], p[1], p[2]] for tri in tris for p in tri])
    F = index.reshape(-1, 3)
    keep = (F[:, 0] != F[:, 1]) & (F[:, 1] != F[:, 2]) & (F[:, 0] != F[:, 2])
    return V, F[keep]


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------


class _Mesh:
    """An indexed mesh with the per-triangle data the algorithm needs."""

    def __init__(self, V, F):
        self.V = V
        self.F = F
        self.P = V[F] if len(F) else np.zeros((0, 3, 3))   # (m, 3, 3)
        if len(F):
            n = np.cross(self.P[:, 1] - self.P[:, 0],
                         self.P[:, 2] - self.P[:, 0])
            length = np.linalg.norm(n, axis=1)
            self.area2 = length
            self.N = np.divide(n, length[:, None], out=np.zeros_like(n),
                               where=length[:, None] > 0)
            self.D = np.einsum("ij,ij->i", self.N, self.P[:, 0])
            self.lo = self.P.min(axis=1)
            self.hi = self.P.max(axis=1)
        else:
            self.area2 = np.zeros(0)
            self.N = np.zeros((0, 3))
            self.D = np.zeros(0)
            self.lo = np.zeros((0, 3))
            self.hi = np.zeros((0, 3))

    def __len__(self):
        return len(self.F)


def _scale(*meshes):
    pts = [m.V for m in meshes if len(m.V)]
    if not pts:
        return 1.0
    allp = np.concatenate(pts)
    return max(float(np.max(allp.max(axis=0) - allp.min(axis=0))), 1e-12)


def candidate_pairs(a, b, pad):
    """For each triangle of ``a``, the triangles of ``b`` whose boxes touch.

    Chunked numpy rather than a tree: deterministic, and fast enough for the
    mesh sizes a pure-Python engine is used at.
    """
    result = [[] for _ in range(len(a))]
    if len(a) == 0 or len(b) == 0:
        return result
    rows = max(1, 4_000_000 // max(len(b), 1))
    for start in range(0, len(a), rows):
        stop = min(start + rows, len(a))
        lo = a.lo[start:stop, None, :] - pad
        hi = a.hi[start:stop, None, :] + pad
        hit = np.all((lo <= b.hi[None, :, :]) & (b.lo[None, :, :] <= hi),
                     axis=2)
        for k, row in enumerate(hit):
            result[start + k] = np.nonzero(row)[0].tolist()
    return result


def winding_numbers(points, mesh, chunk=2_000_000):
    """Generalised winding number of each point with respect to ``mesh``.

    About 1 inside a closed, outward-oriented mesh, 0 outside.  Uses the
    Van Oosterom-Strackee solid-angle formula per triangle.
    """
    points = np.asarray(points, dtype=np.float64).reshape(-1, 3)
    out = np.zeros(len(points))
    if len(points) == 0 or len(mesh) == 0:
        return out
    rows = max(1, chunk // len(mesh))
    A0, A1, A2 = mesh.P[:, 0], mesh.P[:, 1], mesh.P[:, 2]
    for start in range(0, len(points), rows):
        q = points[start:start + rows, None, :]
        a, b, c = A0[None] - q, A1[None] - q, A2[None] - q
        la = np.linalg.norm(a, axis=2)
        lb = np.linalg.norm(b, axis=2)
        lc = np.linalg.norm(c, axis=2)
        det = np.einsum("ijk,ijk->ij", a, np.cross(b, c))
        den = (la * lb * lc + np.einsum("ijk,ijk->ij", a, b) * lc
               + np.einsum("ijk,ijk->ij", b, c) * la
               + np.einsum("ijk,ijk->ij", c, a) * lb)
        out[start:start + rows] = (2.0 * np.arctan2(det, den)).sum(axis=1)
    return out / (4.0 * math.pi)


# ---------------------------------------------------------------------------
# Polygons
# ---------------------------------------------------------------------------


def split_polygon(poly, normal, offset, tol):
    """Split a convex polygon by a plane; returns ``(front, back)``.

    Vertices within ``tol`` of the plane go to both sides.  A side with
    fewer than three vertices is returned as ``None``.
    """
    d = [float(np.dot(normal, p)) - offset for p in poly]
    if all(x >= -tol for x in d) or all(x <= tol for x in d):
        return (poly, None) if any(x > tol for x in d) or \
            all(abs(x) <= tol for x in d) else (None, poly)
    front, back = [], []
    n = len(poly)
    for i in range(n):
        p, q = poly[i], poly[(i + 1) % n]
        dp, dq = d[i], d[(i + 1) % n]
        if dp >= -tol:
            front.append(p)
        if dp <= tol:
            back.append(p)
        if (dp > tol and dq < -tol) or (dp < -tol and dq > tol):
            t = dp / (dp - dq)
            x = p + t * (q - p)
            front.append(x)
            back.append(x)
    return (front if len(front) >= 3 else None,
            back if len(back) >= 3 else None)


def triangulate_convex(poly, normal, tol):
    """Triangulate a convex polygon without losing any vertex.

    :param tol: a *length*: a vertex closer than this to the chord of its
        neighbours is straight, and a triangle whose height above its
        longest edge is below it is degenerate.  A height test, unlike an
        area test, means the same thing at every size.
    :returns: ``(points, triangles)``: the polygon's points -- with one added
        at the end only in the fallback case -- and index triples into them.

    Fragments carry vertices on their straight edges -- exactly the points
    their neighbours share -- and those must survive.  A fan from a polygon
    vertex cannot guarantee it: a side on the apex's own line yields
    collinear, zero-area triangles.  Ear clipping can, given one rule: clip
    a corner only if what remains still has three true corners.  Without
    that rule, clipping a corner can leave a chord lying along a straight
    chain, and the collinear remainder -- with its vertices -- is lost.

    No extra point is added in the normal course.  An earlier version fanned
    around the centroid instead, but on a sliver the centroid lands within
    tolerance of the long edge, its fan triangles fall below the height
    test, and the hole they leave cannot be repaired: the centroid is not
    one of the vertices the T-junction repair searches.
    """
    n = len(poly)
    if n < 3:
        return list(poly), []
    length = float(np.linalg.norm(normal))
    if length == 0.0:
        return list(poly), []
    unit = normal / length

    def height(i, j, k):
        """Signed height of the triangle ``i, j, k`` over its longest edge."""
        twice_area = float(np.dot(np.cross(poly[j] - poly[i],
                                           poly[k] - poly[i]), unit))
        base = max(float(np.linalg.norm(poly[k] - poly[i])),
                   float(np.linalg.norm(poly[j] - poly[i])),
                   float(np.linalg.norm(poly[k] - poly[j])))
        return twice_area / base if base > 0.0 else 0.0

    def corners(ring):
        m = len(ring)
        return sum(height(ring[r - 1], ring[r], ring[(r + 1) % m]) > tol
                   for r in range(m))

    def strict_ear(k, ring):
        i0, i1, i2 = ring[k - 1], ring[k], ring[(k + 1) % len(ring)]
        if height(i0, i1, i2) <= tol:
            return False                          # straight, or a sliver
        rest = ring[:k] + ring[k + 1:]
        if len(rest) > 3 and corners(rest) < 3:
            return False
        return len(rest) > 3 or height(*rest) > tol

    def loose_ear(k, ring):
        i0, i1, i2 = ring[k - 1], ring[k], ring[(k + 1) % len(ring)]
        return height(i0, i1, i2) > 0.0

    # Degenerate means the whole polygon has no height, not that it has
    # fewer than three corners: two genuine corners a hair apart -- the
    # same intersection point computed twice -- read as straight by the
    # corner test, yet the polygon between them has real area.
    centroid = np.mean(poly, axis=0)
    span = max(float(np.linalg.norm(p - centroid)) for p in poly)
    area2 = float(np.dot(sum(np.cross(poly[r - 1] - centroid,
                                      poly[r] - centroid)
                             for r in range(n)), unit))
    if span == 0.0 or area2 / (2.0 * span) <= tol:
        return list(poly), []                     # degenerate: no area
    ring = list(range(n))
    tris = []
    while len(ring) > 3:
        # Prefer ears that leave no sliver; failing that, any ear with
        # positive area.  A near-zero triangle between near-coincident
        # points is harmless -- the final weld collapses it -- whereas
        # leaving the polygon untriangulated would open a hole.
        for ear in (strict_ear, loose_ear):
            k = next((k for k in range(len(ring)) if ear(k, ring)), None)
            if k is not None:
                tris.append((ring[k - 1], ring[k], ring[(k + 1) % len(ring)]))
                ring = ring[:k] + ring[k + 1:]
                break
        else:
            break
    if len(ring) == 3:
        if height(*ring) > 0.0:
            tris.append(tuple(ring))
        return list(poly), tris

    # No valid ear: fall back to a centroid fan over what remains.  This has
    # not been seen on real input; it is here so that a polygon is never
    # silently dropped.
    pts = list(poly) + [np.mean([poly[r] for r in ring], axis=0)]
    c = len(pts) - 1
    for r in range(len(ring)):
        t = (c, ring[r], ring[(r + 1) % len(ring)])
        a, b, d = pts[t[0]], pts[t[1]], pts[t[2]]
        if float(np.dot(np.cross(b - a, d - a), unit)) > 0.0:
            tris.append(t)
    return pts, tris


# ---------------------------------------------------------------------------
# Split and classify one mesh against the other
# ---------------------------------------------------------------------------


def _plane_key(normal, offset, quantum):
    """A sign-normalised key, so a plane and its negation split once."""
    n = normal
    k = int(np.argmax(np.abs(n) > 1e-12))
    if n[k] < 0:
        n, offset = -n, -offset
    return tuple(np.round(np.append(n, offset) / quantum).astype(np.int64))


def _point_in_triangle(p, tri, normal, tol):
    """Is ``p`` (already in the triangle's plane) inside it, edges included?"""
    for i in range(3):
        a, b = tri[i], tri[(i + 1) % 3]
        if float(np.dot(np.cross(b - a, p - a), normal)) < -tol:
            return False
    return True


def fragments(src, other, pairs, tol, cos_tol):
    """Split every triangle of ``src`` against ``other``.

    Returns ``(polys, owner, coplanar, isolated)``: the convex fragments,
    the source triangle each came from, per source triangle the indices of
    ``other`` triangles coplanar with it, and per source triangle whether
    nothing in ``other`` came near it at all.
    """
    polys, owner = [], []
    coplanar = [[] for _ in range(len(src))]
    isolated = np.zeros(len(src), dtype=bool)
    quantum = max(tol, 1e-12)
    for i in range(len(src)):
        tri = [src.P[i, 0], src.P[i, 1], src.P[i, 2]]
        cands = pairs[i]
        if not cands:
            isolated[i] = True
            polys.append(tri)
            owner.append(i)
            continue
        n_i = src.N[i]
        cutters = {}
        for j in cands:
            if other.area2[j] <= 0.0:
                continue
            n_j, d_j = other.N[j], other.D[j]
            dist_i = src.P[i] @ n_j - d_j
            aligned = abs(float(np.dot(n_i, n_j))) >= cos_tol
            if aligned and np.all(np.abs(dist_i) <= tol):
                coplanar[i].append(j)
                # Split along the lines of the coplanar triangle's edges.
                for e in range(3):
                    p, q = other.P[j, e], other.P[j, (e + 1) % 3]
                    m = np.cross(q - p, n_j)
                    length = np.linalg.norm(m)
                    if length > 0:
                        m = m / length
                        cutters.setdefault(_plane_key(m, float(m @ p),
                                                      quantum),
                                           (m, float(m @ p)))
                continue
            # Only planes of triangles that can actually cross this one.
            if np.all(dist_i > tol) or np.all(dist_i < -tol):
                continue
            dist_j = other.P[j] @ n_i - src.D[i]
            if np.all(dist_j > tol) or np.all(dist_j < -tol):
                continue
            cutters.setdefault(_plane_key(n_j, d_j, quantum), (n_j, d_j))
        pieces = [tri]
        for key in sorted(cutters):
            normal, offset = cutters[key]
            next_pieces = []
            for piece in pieces:
                front, back = split_polygon(piece, normal, offset, tol)
                if front is not None:
                    next_pieces.append(front)
                if back is not None:
                    next_pieces.append(back)
            pieces = next_pieces
        for piece in pieces:
            polys.append(piece)
            owner.append(i)
    return polys, owner, coplanar, isolated


def _components(F, members):
    """Connected components, through shared edges, of the triangles in
    ``members`` (a boolean mask).  Returns a label per triangle, -1 if not
    a member."""
    label = np.full(len(F), -1, dtype=np.int64)
    parent = {}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    edge_owner = {}
    for t in np.nonzero(members)[0].tolist():
        parent[t] = t
        a, b, c = F[t].tolist()
        for e in ((a, b), (b, c), (c, a)):
            key = (min(e), max(e))
            if key in edge_owner:
                ra, rb = find(t), find(edge_owner[key])
                if ra != rb:
                    parent[max(ra, rb)] = min(ra, rb)
            else:
                edge_owner[key] = t
    for t in parent:
        label[t] = find(t)
    return label


def classify(src, other, polys, owner, coplanar, isolated, tol):
    """IN / OUT / SAME / OPPOSITE for every fragment of ``src``."""
    classes = [None] * len(polys)
    queries, query_of = [], []

    # Isolated triangles: one winding test per connected region.
    region = _components(src.F, isolated)
    region_query = {}
    for k, i in enumerate(owner):
        if isolated[i]:
            r = int(region[i])
            if r not in region_query:
                region_query[r] = len(queries)
                queries.append(src.P[i].mean(axis=0))
            query_of.append(region_query[r])
            continue
        centre = np.mean(polys[k], axis=0)
        found = None
        for j in coplanar[i]:
            if _point_in_triangle(centre, other.P[j], other.N[j], tol):
                found = SAME if float(np.dot(src.N[i], other.N[j])) > 0 \
                    else OPPOSITE
                break
        if found is not None:
            classes[k] = found
            query_of.append(-1)
        else:
            query_of.append(len(queries))
            queries.append(centre)

    w = winding_numbers(np.array(queries) if queries else np.zeros((0, 3)),
                        other)
    for k, q in enumerate(query_of):
        if q >= 0:
            classes[k] = IN if w[q] > 0.5 else OUT
    return classes


# ---------------------------------------------------------------------------
# Assemble and repair
# ---------------------------------------------------------------------------


def _polygon_area2(loop2d):
    """Twice the signed area of a 2D loop (positive counter-clockwise)."""
    x, y = loop2d[:, 0], loop2d[:, 1]
    return float(np.dot(x, np.roll(y, -1)) - np.dot(np.roll(x, -1), y))


def _inside_2d(pt, loop2d):
    """Even-odd point-in-polygon test in 2D."""
    x, y = pt
    inside = False
    n = len(loop2d)
    for i in range(n):
        x0, y0 = loop2d[i]
        x1, y1 = loop2d[(i + 1) % n]
        if (y0 > y) != (y1 > y):
            if x < x0 + (y - y0) * (x1 - x0) / (y1 - y0):
                inside = not inside
    return inside


def merge_fragments(polys, origin, normal, tol, flat_tol):
    """Re-triangulate the region a triangle's kept fragments cover.

    Cutting by infinite planes slices a triangle where its cutter never was,
    so a large face crossed by a finely tessellated part can shatter into
    many slivers even though most of its fragments end up in the result.
    The kept fragments of one triangle share their split points exactly --
    every cutter splits every fragment -- so their interior edges cancel,
    leaving the region's boundary loops, which earcut triangulates afresh.
    Vertices that only existed on phantom cuts disappear; the ones shared
    with neighbouring triangles come back in the T-junction repair.

    Returns 3D triangles, or ``None`` to mean "emit the fragments as they
    are" -- for a pinched boundary, an ambiguous hole, or a result whose
    area does not match the fragments'.  Merging is an optimisation, so it
    only ever happens when it is provably area-preserving.
    """
    import mapbox_earcut

    n = normal / np.linalg.norm(normal)
    helper = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.9 else \
        np.array([0.0, 1.0, 0.0])
    u = np.cross(helper, n)
    u /= np.linalg.norm(u)
    v = np.cross(n, u)

    pts, index = weld([p for poly in polys for p in poly], quantum=tol)
    ids, k = [], 0
    for poly in polys:
        ids.append([int(i) for i in index[k:k + len(poly)]])
        k += len(poly)

    edges = set()
    for ring in ids:
        ring = [r for j, r in enumerate(ring) if r != ring[j - 1]]
        for j in range(len(ring)):
            edges.add((ring[j], ring[(j + 1) % len(ring)]))
    boundary = {(a, b) for a, b in edges if (b, a) not in edges and a != b}
    if not boundary:
        return None
    nxt = {}
    for a, b in boundary:
        if a in nxt:
            return None                       # pinch point: do not guess
        nxt[a] = b

    loops, seen = [], set()
    for start in sorted(nxt):
        if start in seen:
            continue
        loop, cur = [], start
        while cur not in seen:
            seen.add(cur)
            loop.append(cur)
            cur = nxt.get(cur)
            if cur is None:
                return None
        if cur != start or len(loop) < 3:
            return None
        loops.append(loop)

    flat = (pts - origin) @ np.stack([u, v], axis=1)
    outers, holes = [], []
    for loop in loops:
        area = _polygon_area2(flat[loop])
        (outers if area > 0 else holes).append(loop)
    if not outers:
        return None
    groups = {i: [] for i in range(len(outers))}
    for hole in holes:
        owners = [i for i, outer in enumerate(outers)
                  if _inside_2d(flat[hole[0]], flat[outer])]
        if len(owners) != 1:
            return None
        groups[owners[0]].append(hole)

    tris = []
    for i, outer in enumerate(outers):
        rings = [outer] + groups[i]
        order = [vid for ring in rings for vid in ring]
        ends = np.cumsum([len(ring) for ring in rings]).astype(np.uint32)
        idx = mapbox_earcut.triangulate_float64(
            np.ascontiguousarray(flat[order], dtype=np.float64), ends)
        for t in range(0, len(idx), 3):
            a, b, c = (order[int(idx[t + j])] for j in range(3))
            pa, pb, pc = pts[a], pts[b], pts[c]
            twice = float(np.dot(np.cross(pb - pa, pc - pa), n))
            base = max(float(np.linalg.norm(pb - pa)),
                       float(np.linalg.norm(pc - pa)),
                       float(np.linalg.norm(pc - pb)))
            if base == 0.0 or abs(twice) / base <= flat_tol:
                # earcut can emit a zero-area triangle along a straight run
                # of boundary points; it covers nothing, and the point it
                # sits on comes back in the T-junction repair.
                continue
            if twice < 0:
                pb, pc = pc, pb
            tris.append((pa, pb, pc))

    def area(tri_list):
        return sum(0.5 * float(np.linalg.norm(np.cross(t[1] - t[0],
                                                       t[2] - t[0])))
                   for t in tri_list)

    target = sum(0.5 * abs(_polygon_area2(((np.asarray(poly) - origin)
                                           @ np.stack([u, v], axis=1))))
                 for poly in polys)
    if abs(area(tris) - target) > 1e-9 * max(target, tol * tol) + tol * tol:
        return None
    return tris


def _repair_t_junctions(V, F, tol, flat_tol):
    """Split every edge that has a welded vertex lying strictly inside it.

    After welding, a crack can only be a vertex on another triangle's edge;
    inserting it there makes the mesh conforming without moving anything.
    """
    if len(F) == 0:
        return V, F
    cell = max(float(np.median(np.linalg.norm(
        V[F[:, 1]] - V[F[:, 0]], axis=1))), tol * 10.0)
    grid = {}
    keys = np.floor(V / cell).astype(np.int64)
    for v, key in enumerate(map(tuple, keys)):
        grid.setdefault(key, []).append(v)

    def on_edge(a, b):
        p, q = V[a], V[b]
        lo = np.floor((np.minimum(p, q) - tol) / cell).astype(np.int64)
        hi = np.floor((np.maximum(p, q) + tol) / cell).astype(np.int64)
        d = q - p
        dd = float(d @ d)
        if dd == 0.0:
            return []
        found = []
        for x in range(lo[0], hi[0] + 1):
            for y in range(lo[1], hi[1] + 1):
                for z in range(lo[2], hi[2] + 1):
                    for v in grid.get((x, y, z), ()):
                        if v == a or v == b:
                            continue
                        t = float((V[v] - p) @ d) / dd
                        # Within tol of an endpoint is the endpoint itself
                        # (the cluster weld has merged such points); only
                        # genuinely interior points split the edge.
                        span = math.sqrt(dd)
                        if t * span <= tol or (1.0 - t) * span <= tol:
                            continue
                        if np.linalg.norm(V[v] - (p + t * d)) <= tol:
                            found.append((t, v))
        found.sort()
        return [v for _, v in found]

    cache = {}
    out = []
    extra = []
    for a, b, c in F.tolist():
        ring, lines = [a], [frozenset((0, 2))]
        corner_lines = {b: frozenset((0, 1)), c: frozenset((1, 2))}
        for edge, (s, e) in enumerate(((a, b), (b, c), (c, a))):
            key = (s, e)
            if key not in cache:
                rev = cache.get((e, s))
                cache[key] = rev[::-1] if rev is not None else on_edge(s, e)
            ring.extend(cache[key])
            lines.extend([frozenset((edge,))] * len(cache[key]))
            if e != a:
                ring.append(e)
                lines.append(corner_lines[e])
        if len(ring) == 3:
            out.append((a, b, c))
            continue
        tris = _triangulate_labelled(lines)
        if tris is not None:
            out.extend((ring[i0], ring[i1], ring[i2]) for i0, i1, i2 in tris)
            continue
        # Not seen on real input: fan round the centroid, keeping every
        # triangle -- this is a repair, and dropping one would open a hole.
        centre = len(V) + len(extra)
        extra.append(np.mean(V[ring], axis=0))
        out.extend((centre, ring[r], ring[(r + 1) % len(ring)])
                   for r in range(len(ring)))
    F2 = np.asarray(out, dtype=np.int64).reshape(-1, 3)
    if extra:
        V = np.vstack([V, np.asarray(extra)])
    return V, F2


def solid_boolean(a, b, operation):
    """``(V, F)`` of the boolean of two yapCAD solids."""
    if operation not in KEEP:
        raise ValueError(f"unsupported solid boolean operation {operation!r}")
    A = _Mesh(*solid_to_mesh(a))
    B = _Mesh(*solid_to_mesh(b))
    scale = _scale(A, B)
    tol = 1e-9 * scale
    cos_tol = 1.0 - 1e-9
    pad = 10.0 * tol
    # Straightness and degeneracy are judged by height, in length units.
    flat_tol = 10.0 * tol

    keep_a, keep_b = KEEP[operation]
    triangles = []
    for src, other, keep, flip in ((A, B, keep_a, False),
                                   (B, A, keep_b, operation == "difference")):
        pairs = candidate_pairs(src, other, pad)
        polys, owner, coplanar, isolated = fragments(src, other, pairs, tol,
                                                     cos_tol)
        classes = classify(src, other, polys, owner, coplanar, isolated, tol)
        # A triangle every one of whose fragments is kept is simply kept:
        # its fragments tile it exactly.  This undoes the cuts an infinite
        # plane makes where its triangle never was -- harmless to
        # correctness, but a source of slivers -- and the points neighbours
        # share on its edges come back in the T-junction repair.
        whole = {}
        for i, cls in zip(owner, classes):
            whole[i] = whole.get(i, True) and cls in keep
        kept = {}
        for poly, i, cls in zip(polys, owner, classes):
            if cls in keep:
                kept.setdefault(i, []).append(
                    [np.asarray(p, dtype=np.float64) for p in poly])
        for i in sorted(kept):
            if whole[i]:
                out_tris = [(src.P[i, 0], src.P[i, 1], src.P[i, 2])]
            else:
                out_tris = None
                if len(kept[i]) > 1:
                    out_tris = merge_fragments(kept[i], src.P[i, 0],
                                               src.N[i], tol, flat_tol)
                if out_tris is None:
                    out_tris = []
                    for poly in kept[i]:
                        pts, tris = triangulate_convex(poly, src.N[i],
                                                       flat_tol)
                        out_tris.extend((pts[a], pts[b], pts[c])
                                        for a, b, c in tris)
            for tri in out_tris:
                triangles.append(tri[::-1] if flip else tri)

    if not triangles:
        return np.zeros((0, 3)), np.zeros((0, 3), dtype=np.int64)
    # Weld and repair at the algorithm's own geometric tolerance, so that a
    # vertex counts as lying on an edge only when it genuinely does.  Doing
    # this at yapCAD's coarser quantum let near-miss vertices be "inserted"
    # into edges they were merely close to, bending them out of line.
    V, index = weld([p for tri in triangles for p in tri], quantum=tol)
    F = _drop_collapsed(index.reshape(-1, 3))
    V, F = cluster_weld(V, F, 10.0 * tol)
    V, F = _repair_t_junctions(V, F, 10.0 * tol, flat_tol)
    # Then re-weld once at yapCAD's quantum, so vertex identity in the
    # output is exactly what issolidclosed will compute.
    V, index = weld(V[F].reshape(-1, 3))
    return V, _drop_collapsed(index.reshape(-1, 3))


def _triangulate_labelled(lines):
    """Triangulate a triangle carrying extra points on its edges, exactly.

    ``lines[i]`` is the set of the triangle's original edges (0, 1, 2) that
    ring vertex ``i`` lies on: two for a corner, one for an inserted point.
    Three vertices are collinear exactly when those sets share an edge, so
    straightness is decided combinatorially, with no tolerance at all.

    The geometric test this replaces failed when an inserted point sat a
    couple of tolerances from a corner: measured against the tiny chord to
    that point, the genuine corner's height fell below the threshold, the
    triangle seemed to have only two corners, and it was dropped whole.
    Ears are clipped only where the remainder keeps three true corners, so
    no zero-area triangle is ever emitted and no vertex is lost.  Returns
    index triples into the ring, or ``None`` if no valid ear exists.
    """
    def collinear(i, j, k):
        return bool(lines[i] & lines[j] & lines[k])

    def corners(ring):
        m = len(ring)
        return sum(not collinear(ring[r - 1], ring[r], ring[(r + 1) % m])
                   for r in range(m))

    ring = list(range(len(lines)))
    tris = []
    while len(ring) > 3:
        for k in range(len(ring)):
            p, v, q = ring[k - 1], ring[k], ring[(k + 1) % len(ring)]
            if collinear(p, v, q):
                continue
            rest = ring[:k] + ring[k + 1:]
            if len(rest) == 3 and collinear(*rest):
                continue
            if len(rest) > 3 and corners(rest) < 3:
                continue
            tris.append((p, v, q))
            ring = rest
            break
        else:
            return None
    if collinear(*ring):
        return None
    tris.append(tuple(ring))
    return tris


def cluster_weld(V, F, radius):
    """Merge every pair of vertices closer than ``radius``.

    Grid rounding cannot do this: two points a hair apart can straddle a
    cell boundary and stay distinct.  Instead each vertex is compared with
    those in its own and neighbouring cells of a ``radius``-sized grid, and
    close pairs are joined with union-find; each cluster takes its
    lowest-numbered vertex, so the result is deterministic.

    This has to match the T-junction repair's tolerance.  When the weld was
    finer than the repair, a vertex a few tolerances from a corner survived
    as a separate point, the repair inserted it next to the corner, the two
    near-coincident points stopped counting as a corner, and whole healthy
    triangles were judged degenerate and lost.
    """
    if len(V) == 0:
        return V, F
    parent = list(range(len(V)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    cells = {}
    keys = np.floor(V / radius).astype(np.int64)
    for v, key in enumerate(map(tuple, keys.tolist())):
        cells.setdefault(key, []).append(v)
    offsets = [(dx, dy, dz) for dx in (-1, 0, 1) for dy in (-1, 0, 1)
               for dz in (-1, 0, 1)]
    for key, members in cells.items():
        near = []
        for dx, dy, dz in offsets:
            near.extend(cells.get((key[0] + dx, key[1] + dy, key[2] + dz),
                                  ()))
        near = np.asarray(near)
        for v in members:
            close = near[np.linalg.norm(V[near] - V[v], axis=1) <= radius]
            for u in close.tolist():
                ru, rv = find(u), find(v)
                if ru != rv:
                    parent[max(ru, rv)] = min(ru, rv)
    roots = np.array([find(v) for v in range(len(V))])
    used, remap = np.unique(roots, return_inverse=True)
    return V[used], _drop_collapsed(remap[F])


def _drop_collapsed(F):
    """Remove triangles with a repeated vertex: they cover nothing, and
    removing one from a conforming mesh leaves it conforming."""
    if len(F) == 0:
        return F
    keep = (F[:, 0] != F[:, 1]) & (F[:, 1] != F[:, 2]) & (F[:, 0] != F[:, 2])
    return F[keep]
