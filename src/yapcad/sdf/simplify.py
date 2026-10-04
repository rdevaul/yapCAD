"""Simplify a meshed field, checking every change against the field.

Dual contouring samples uniformly, so a part needs one cell size everywhere
-- the one its finest feature demands.  A threaded nut meshed finely enough
for its 0.16 mm crest flats is meshed just as finely across its flat faces,
where a handful of triangles would do.  Simplification recovers that, and an
SDF makes it unusually well posed: the field says, for any candidate
triangle, how far it strays from the true part.

Edges are collapsed by quadric error (Garland and Heckbert), in rounds.  A
round scores every edge at once, picks greedily, cheapest first, a set of
edges whose neighbourhoods do not overlap, and tests all of them together:

* **topology** -- the link condition, so the mesh stays a closed manifold;
* **orientation** -- no triangle around the merged vertex flips or
  collapses;
* **fidelity** -- the merged vertex is first moved onto the surface (Newton
  steps along the gradient), then the field is evaluated at it and at
  sample points across every triangle that changes.  A collapse is kept
  only if ``|f|`` stays within the tolerance everywhere it is sampled --
  or, where the uniform mesh already strayed further than that, as along
  sharp edges that dual contouring rounds, no further than it already
  did.  The result is never worse than the tolerance or the original.

The rounds repeat until no collapse passes.  Everything is array-at-a-time
and deterministic -- stable sorts, no randomness -- so a simplified preview
regenerates bit-for-bit on the same platform.  Not necessarily across
platforms: the quadric solves and field evaluations can differ in the last
bit between builds of numpy and LAPACK, and a collapse whose error sits
right at the tolerance can then pass on one and fail on the other.  The
results differ in detail, not in quality -- both meet the same bounds.

For an exact field ``|f|`` is the distance to the part.  For an inexact one
it is the field's value, which for hard CSG and the other bounds yapCAD's
fields provide does not exceed the distance; the tolerance is then on that
value.
"""

import math

import numpy as np

from yapcad.sdf.evaluate import evaluate, gradient
from yapcad.sdf.node import SdfError

__all__ = ["simplify_mesh", "simplify_solid"]

#: Most subdivisions per edge of the barycentric sample grid on a triangle.
_MAX_DIVISIONS = 32


def _lattice(n):
    """Barycentric points of a triangle subdivided ``n`` times per edge,
    corners excluded (they are vertices, checked separately or unchanged)."""
    pts = [(i / n, j / n, (n - i - j) / n)
           for i in range(n + 1) for j in range(n + 1 - i)]
    return np.array([q for q in pts if max(q) < 1.0] or [(1 / 3,) * 3])


_LATTICES = {n: _lattice(n) for n in range(1, _MAX_DIVISIONS + 1)}


def _max_error(node, tris, spacing):
    """Worst ``|f|`` sampled across each triangle of ``tris`` ``(T, 3, 3)``
    on a grid no coarser than ``spacing`` -- the original mesh's edge
    length -- so a large merged triangle is checked as densely as the many
    small ones it replaces."""
    longest = np.max(np.linalg.norm(tris - np.roll(tris, 1, axis=1), axis=2),
                     axis=1)
    # Two samples per original edge length, and at least four divisions:
    # a triangle cut across a sharp edge bulges away from it between
    # samples, and two divisions once read 0.095 where the truth was 0.142.
    divisions = np.clip(np.ceil(2.0 * longest / spacing), 4,
                        _MAX_DIVISIONS).astype(int)
    out = np.zeros(len(tris))
    for n in np.unique(divisions).tolist():
        rows = np.nonzero(divisions == n)[0]
        pts = np.einsum("sk,fkd->fsd", _LATTICES[n], tris[rows])
        err = np.abs(evaluate(node, pts.reshape(-1, 3)))
        out[rows] = err.reshape(len(rows), -1).max(axis=1)
    return out


#: A changed triangle's normal must stay within ~78 degrees of its old one.
_MIN_NORMAL_DOT = 0.2


def _face_normals(V, F):
    n = np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]])
    return n, np.linalg.norm(n, axis=1)


def _quadrics(V, F):
    """Area-weighted plane quadrics summed onto vertices: ``(V, 4, 4)``."""
    n, length = _face_normals(V, F)
    ok = length > 0.0
    unit = np.zeros_like(n)
    unit[ok] = n[ok] / length[ok, None]
    d = -np.einsum("ij,ij->i", unit, V[F[:, 0]])
    plane = np.column_stack([unit, d])                       # (F, 4)
    K = plane[:, :, None] * plane[:, None, :] * (0.5 * length)[:, None, None]
    Q = np.zeros((len(V), 4, 4))
    for k in range(3):
        np.add.at(Q, F[:, k], K)
    return Q


def _edges(F):
    e = np.concatenate([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]])
    e.sort(axis=1)
    return np.unique(e, axis=0)


def _csr(keys, values, n):
    """Neighbour lists as CSR: ``values[start[k]:start[k + 1]]``."""
    order = np.argsort(keys, kind="stable")
    start = np.zeros(n + 1, dtype=np.int64)
    np.add.at(start, keys + 1, 1)
    return np.cumsum(start), values[order]


def _candidates(V, Q, E):
    """Quadric-optimal merged positions for each edge, falling back to the
    midpoint when the quadric is singular or its optimum is off the edge."""
    a, b = E[:, 0], E[:, 1]
    q = Q[a] + Q[b]
    mid = 0.5 * (V[a] + V[b])
    A, rhs = q[:, :3, :3], -q[:, :3, 3]
    det = np.linalg.det(A)
    scale = np.einsum("ijj->i", A) ** 3 + 1e-300
    solvable = np.abs(det) > 1e-9 * scale
    pos = mid.copy()
    if solvable.any():
        pos[solvable] = np.linalg.solve(A[solvable],
                                        rhs[solvable, :, None])[:, :, 0]
    span = np.linalg.norm(V[a] - V[b], axis=1)
    far = np.linalg.norm(pos - mid, axis=1) > span
    pos[far] = mid[far]
    return pos


def _project(node, P, steps=2):
    """Move points onto the surface by Newton steps along the gradient."""
    P = P.copy()
    for _ in range(steps):
        f = evaluate(node, P)
        g = gradient(node, P)
        gg = np.einsum("ij,ij->i", g, g)
        ok = gg > 1e-24
        P[ok] -= (f[ok] / gg[ok])[:, None] * g[ok]
    return P


def _quadric_cost(Q, E, P):
    h = np.column_stack([P, np.ones(len(P))])
    q = Q[E[:, 0]] + Q[E[:, 1]]
    return np.einsum("ij,ijk,ik->i", h, q, h)


def _round(node, V, F, tolerance, rejected, spacing, face_err):
    """One round of collapses.

    :param rejected: ``(k, 2)`` edges tested and rejected in earlier rounds
        whose neighbourhoods have not changed since: the reasons still hold,
        so they are not tested again.  This is what keeps the long tail of
        rounds, each collapsing a few dozen edges, cheap.  An edge that was
        merely locked out by a neighbour's collapse is not in it.
    :param face_err: each triangle's sampled ``|f|``, kept up to date across
        rounds: a triangle's error changes only when a collapse rewrites it,
        so it is measured once, not every round it borders a candidate.
    :returns: ``(V, F, collapsed, rejected, face_err)`` for the next round.
    """
    E = _edges(F)
    if len(E) == 0:
        return V, F, 0, rejected, face_err
    Q = _quadrics(V, F)
    P = _candidates(V, Q, E)
    cost = _quadric_cost(Q, E, P)
    order = np.lexsort((np.arange(len(E)), cost))

    nv = len(V)
    vstart, vnbr = _csr(np.concatenate([E[:, 0], E[:, 1]]),
                        np.concatenate([E[:, 1], E[:, 0]]), nv)
    face_of = np.repeat(np.arange(len(F)), 3)
    fstart, ffaces = _csr(F.reshape(-1), face_of, nv)

    # Greedy independent set: an accepted edge locks both endpoints and
    # their neighbours, so no two collapses in a round touch one triangle.
    if len(rejected):
        known = np.isin(E[order, 0] * nv + E[order, 1],
                        rejected[:, 0] * nv + rejected[:, 1])
        order = order[~known]
    locked = np.zeros(nv, dtype=bool)
    chosen = []
    for e in order.tolist():
        a, b = int(E[e, 0]), int(E[e, 1])
        if locked[a] or locked[b]:
            continue
        na = vnbr[vstart[a]:vstart[a + 1]]
        nb = vnbr[vstart[b]:vstart[b + 1]]
        # Link condition: the endpoints may share only the two vertices
        # opposite the edge, or the collapse pinches the surface.
        if len(np.intersect1d(na, nb, assume_unique=True)) != 2:
            continue
        locked[a] = locked[b] = True
        locked[na] = True
        locked[nb] = True
        chosen.append(e)
    if not chosen:
        return V, F, 0, rejected, face_err
    chosen = np.asarray(chosen)
    # Only the chosen candidates are moved onto the surface: projection
    # costs a field and gradient evaluation per point per step.
    P = P.copy()
    P[chosen] = _project(node, P[chosen])

    # Every triangle around each chosen edge, with the merged vertex in
    # place, except the two that contain the edge and vanish.
    owner, faces = [], []
    for k, e in enumerate(chosen.tolist()):
        a, b = int(E[e, 0]), int(E[e, 1])
        around = np.union1d(ffaces[fstart[a]:fstart[a + 1]],
                            ffaces[fstart[b]:fstart[b + 1]])
        faces.append(around)
        owner.append(np.full(len(around), k))
    faces = np.concatenate(faces)
    owner = np.concatenate(owner)
    tri = F[faces]
    a_of, b_of = E[chosen, 0][owner], E[chosen, 1][owner]
    survives = ~((tri == a_of[:, None]).any(1) & (tri == b_of[:, None]).any(1))
    moved = (tri == a_of[:, None]) | (tri == b_of[:, None])
    corners = V[tri]
    new = np.where(moved[:, :, None], P[chosen][owner][:, None, :], corners)

    old_n, _ = _face_normals(V, F[faces])
    new_n = np.cross(new[:, 1] - new[:, 0], new[:, 2] - new[:, 0])
    new_len = np.linalg.norm(new_n, axis=1)
    old_len = np.linalg.norm(old_n, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        dot = np.einsum("ij,ij->i", old_n, new_n) / (old_len * new_len)
    bad_face = survives & ~((new_len > 0.0) & (dot > _MIN_NORMAL_DOT))

    # Fidelity: the new triangles' worst sampled error must be within the
    # tolerance -- or, where dual contouring already left more error than
    # that (along sharp edges, typically), no worse than what was there.
    # Otherwise a sharp edge would freeze every collapse near it.
    # The allowance must not spread: a collapse next to one bad triangle
    # may not let all its new triangles reach that triangle's error, or the
    # bad region grows a ring every round.  So the worst error may not rise,
    # and the excess over the tolerance, integrated over area, may not grow.
    new_face_err = _max_error(node, new[survives], spacing)
    old_face_err = face_err[faces]
    new_err = np.zeros(len(chosen))
    old_err = np.zeros(len(chosen))
    np.maximum.at(new_err, owner[survives], new_face_err)
    np.maximum.at(old_err, owner, old_face_err)
    new_excess = np.zeros(len(chosen))
    old_excess = np.zeros(len(chosen))
    np.add.at(new_excess, owner[survives],
              new_len[survives] * np.maximum(new_face_err - tolerance, 0.0))
    np.add.at(old_excess, owner,
              old_len * np.maximum(old_face_err - tolerance, 0.0))
    vertex_err = np.abs(evaluate(node, P[chosen]))
    allowed = np.maximum(tolerance, old_err)

    ok = ((vertex_err <= allowed) & (new_err <= allowed)
          & (new_excess <= old_excess * (1.0 + 1e-9)))
    ok[np.unique(owner[bad_face])] = False
    rejected = np.concatenate([rejected, E[chosen[~ok]]])
    if not ok.any():
        return V, F, 0, rejected, face_err

    # A rejection stands only while nothing around it changes: forget any
    # that touch a vertex of a triangle around one of these collapses.
    touched = np.zeros(nv, dtype=bool)
    touched[F[faces[ok[owner]]].reshape(-1)] = True
    rejected = rejected[~(touched[rejected[:, 0]] | touched[rejected[:, 1]])]

    applied = ok[owner] & survives
    face_err = face_err.copy()
    face_err[faces[applied]] = new_face_err[applied[survives]]

    keep, drop = E[chosen[ok], 0], E[chosen[ok], 1]
    V = V.copy()
    V[keep] = P[chosen[ok]]
    remap = np.arange(nv)
    remap[drop] = keep
    F = remap[F]
    alive = (F[:, 0] != F[:, 1]) & (F[:, 1] != F[:, 2]) & (F[:, 0] != F[:, 2])
    F, face_err = F[alive], face_err[alive]
    used, F = np.unique(F, return_inverse=True)
    renumber = np.full(nv, -1, dtype=np.int64)
    renumber[used] = np.arange(len(used))
    rejected = renumber[rejected]
    rejected = rejected[(rejected >= 0).all(axis=1)]
    rejected.sort(axis=1)
    return V[used], F.reshape(-1, 3), int(ok.sum()), rejected, face_err


def simplify_mesh(node, vertices, triangles, tolerance, max_rounds=200):
    """Simplify a closed triangle mesh of ``node`` to within ``tolerance``.

    :returns: ``(vertices, triangles)`` arrays.
    """
    try:
        tol = float(tolerance)
    except (TypeError, ValueError):
        raise SdfError("simplify: tolerance must be a number") from None
    if not math.isfinite(tol) or tol <= 0.0:
        raise SdfError("simplify: tolerance must be positive and finite")
    V = np.asarray(vertices, dtype=float)
    F = np.asarray(triangles, dtype=np.int64).reshape(-1, 3)
    rejected = np.zeros((0, 2), dtype=np.int64)
    edges = _edges(F)
    spacing = float(np.median(np.linalg.norm(V[edges[:, 0]] - V[edges[:, 1]],
                                             axis=1))) if len(edges) else 1.0
    if len(F) == 0:
        return V, F
    face_err = _max_error(node, V[F], spacing)
    for _ in range(int(max_rounds)):
        V, F, done, rejected, face_err = _round(node, V, F, tol, rejected,
                                                spacing, face_err)
        if done == 0:
            break
    return V, F


def simplify_solid(solid, tolerance):
    """Simplify every surface of an SDF-authored solid against its field.

    The construction record keeps the field and gains the tolerance, so
    the simplified preview is as reproducible as the original.
    """
    from yapcad.construction import get_construction
    from yapcad.geom3d import solid as make_solid
    from yapcad.geom3d import surface
    from yapcad.sdf.evaluate import normal
    from yapcad.sdf.node import from_construction, \
        meshing_from_construction, to_construction

    from yapcad.sdf.booleans import is_sdf_solid
    if not is_sdf_solid(solid):
        raise SdfError("simplify_solid: the solid is not SDF-authored")
    record = get_construction(solid)
    node = from_construction(record)
    surfaces = []
    for surf in solid[1]:
        V = np.array([p[:3] for p in surf[1]], dtype=float)
        V, F = simplify_mesh(node, V, surf[3], tolerance)
        N = normal(node, V)
        surfaces.append(surface(
            [[float(x), float(y), float(z), 1.0] for x, y, z in V],
            [[float(x), float(y), float(z), 0.0] for x, y, z in N],
            [[int(i), int(j), int(k)] for i, j, k in F]))
    meshing = dict(meshing_from_construction(record) or {})
    if meshing:
        meshing["simplify"] = {"method": "field-qem",
                               "tolerance": float(tolerance)}
    return make_solid(surfaces, [], to_construction(
        node, meshing=meshing or None))
