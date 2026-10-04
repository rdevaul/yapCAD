"""Ground-truth tests for the native mesh boolean engine.

The engine's existing tests check point containment and outward normals,
never whether a result is closed or has the right volume, so open and
double-walled results passed.  These check both, against closed-form
volumes, with ``engine="native"`` forced -- otherwise, with OCC installed,
``prism()`` carries a BREP and the booleans quietly go through OCC.

These began as regression tests for the eight cases the original engine got
right, plus strict expected-failure markers for the thirteen it got wrong --
gaps from independently cut triangles, and coincident faces never classified
as same or opposite.  The rewritten core (:mod:`yapcad.boolean.csg`) passes
all of them, so the markers are gone and they are all regression tests now.
"""

import math

import pytest

from yapcad.brep import _clear_brep_data
from yapcad.geom import point
from yapcad.geom3d import issolidclosed, solid_boolean, translatesolid
from yapcad.geom3d_util import conic, prism, sphere


def mesh_only(solid):
    _clear_brep_data(solid)
    return solid


def moved(solid, delta):
    return mesh_only(translatesolid(solid, point(*delta)))


def signed_volume(solid):
    """Read from the surfaces directly: yapcad.mesh.mesh_view drops faces
    with area <= epsilon, and on small results that measurably changes the
    volume -- enough to break the identities below by ~1e-6."""
    total = 0.0
    for surf in solid[1]:
        v = surf[1]
        for i, j, k in surf[3]:
            a, b, c = v[i], v[j], v[k]
            total += (a[0] * (b[1] * c[2] - b[2] * c[1])
                      - a[1] * (b[0] * c[2] - b[2] * c[0])
                      + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6.0
    return total


def native(a, b, op):
    return solid_boolean(a, b, op, engine="native")


def assert_correct(result, volume):
    assert issolidclosed(result), "result is not a closed solid"
    assert signed_volume(result) == pytest.approx(volume, rel=1e-6,
                                                  abs=1e-9)


# ---------------------------------------------------------------------------
# What the engine gets right
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("op,volume", [
    ("union", 64.0), ("intersection", 8.0), ("difference", 56.0),
])
def test_nested_boxes(op, volume):
    """No triangle crosses the other surface, so every one is classified by
    a single containment test -- the path the octree fallback used to
    replace with a clip against every target plane."""
    assert_correct(native(mesh_only(prism(4, 4, 4)),
                          mesh_only(prism(2, 2, 2)), op), volume)


@pytest.mark.parametrize("op,volume", [
    ("union", 8.0), ("difference", 0.0),
])
def test_a_sphere_inside_a_box(op, volume):
    # geom3d_util.sphere takes a diameter: this one sits wholly inside.
    assert_correct(native(mesh_only(prism(2, 2, 2)),
                          mesh_only(sphere(1.2)), op)
                   if op == "union" else
                   native(mesh_only(sphere(1.2)),
                          mesh_only(prism(2, 2, 2)), op), volume)


def test_disjoint_union_keeps_both():
    assert_correct(native(mesh_only(prism(2, 2, 2)),
                          moved(prism(2, 2, 2), (5.0, 0.0, 0.0)), "union"),
                   16.0)


@pytest.mark.parametrize("op,volume", [
    ("intersection", 0.0), ("difference", 8.0),
])
def test_face_sharing_boxes_except_union(op, volume):
    assert_correct(native(mesh_only(prism(2, 2, 2)),
                          moved(prism(2, 2, 2), (2.0, 0.0, 0.0)), op),
                   volume)


# ---------------------------------------------------------------------------
# Formerly known defects, now fixed
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("op,volume", [
    ("union", 12.175), ("intersection", 3.825), ("difference", 4.175),
])
def test_overlapping_boxes_in_general_position(op, volume):
    """The old engine got the volume but left gaps between triangles."""
    assert_correct(native(mesh_only(prism(2, 2, 2)),
                          moved(prism(2, 2, 2), (0.75, 0.3, 0.2)), op),
                   volume)


def test_face_sharing_union():
    assert_correct(native(mesh_only(prism(2, 2, 2)),
                          moved(prism(2, 2, 2), (2.0, 0.0, 0.0)), "union"),
                   16.0)


@pytest.mark.parametrize("op,volume", [
    ("union", 12.0), ("intersection", 4.0), ("difference", 4.0),
])
def test_boxes_with_coplanar_side_faces(op, volume):
    """Four side faces coincide over the overlap; the old engine's volumes
    were 67% wrong."""
    assert_correct(native(mesh_only(prism(2, 2, 2)),
                          moved(prism(2, 2, 2), (1.0, 0.0, 0.0)), op),
                   volume)


@pytest.mark.parametrize("op", ["union", "intersection", "difference"])
def test_box_and_through_cylinder(op):
    result = native(mesh_only(prism(2, 2, 2)),
                    mesh_only(conic(0.6, 0.6, 3.0,
                                    center=point(0, 0, -1.5))), op)
    assert issolidclosed(result)


@pytest.mark.parametrize("op", ["union", "intersection", "difference"])
def test_overlapping_spheres(op):
    result = native(mesh_only(sphere(1.0)),
                    moved(sphere(1.0), (0.8, 0.1, 0.0)), op)
    assert issolidclosed(result)


# ---------------------------------------------------------------------------
# Volume identities: correctness without any reference engine
# ---------------------------------------------------------------------------
#
# For any correct boolean, vol(A | B) = vol(A) + vol(B) - vol(A & B) and
# vol(A - B) = vol(A) - vol(A & B).  These need nothing but the engine
# itself, so they run in the pure-Python lane too.

IDENTITY_CASES = [
    ("boxes, general position",
     lambda: (mesh_only(prism(2, 2, 2)),
              moved(prism(2, 2, 2), (0.75, 0.3, 0.2)))),
    ("boxes, coincident faces",
     lambda: (mesh_only(prism(2, 2, 2)), moved(prism(2, 2, 2), (1, 0, 0)))),
    ("box and cylinder",
     lambda: (mesh_only(prism(2, 2, 2)),
              mesh_only(conic(0.6, 0.6, 3.0, center=point(0, 0, -1.5))))),
    ("spheres",
     lambda: (mesh_only(sphere(1.0)), moved(sphere(1.0), (0.8, 0.1, 0.0)))),
    ("rotated box and sphere",
     lambda: (mesh_only(prism(2, 2, 2)),
              mesh_only(__import__("yapcad.geom3d", fromlist=["x"])
                        .rotatesolid(moved(sphere(1.6), (0.7, 0.2, 0.4)),
                                     33.0, axis=point(0.6, 0.0, 0.8))))),
]


@pytest.mark.parametrize("name,build", IDENTITY_CASES,
                         ids=[c[0] for c in IDENTITY_CASES])
def test_volume_identities_hold(name, build):
    a, b = build()
    results = {op: native(a, b, op)
               for op in ("union", "intersection", "difference")}
    for op, result in results.items():
        assert issolidclosed(result), f"{op} is not closed"
    va, vb = signed_volume(a), signed_volume(b)
    vu, vi, vd = (signed_volume(results[op])
                  for op in ("union", "intersection", "difference"))
    scale = max(va, vb)
    assert vu == pytest.approx(va + vb - vi, abs=1e-7 * scale)
    assert vd == pytest.approx(va - vi, abs=1e-7 * scale)


def test_results_are_deterministic():
    a, b = mesh_only(sphere(1.0)), moved(sphere(1.0), (0.8, 0.1, 0.0))
    first, second = native(a, b, "union"), native(a, b, "union")
    assert first[1][0][1] == second[1][0][1]
    assert first[1][0][3] == second[1][0][3]


def test_coaxial_cylinders_do_not_shatter():
    """A bore through a bearing seat, as in a printed boss.  Cutting by the
    full plane of every triangle a cap's fan wedge merely straddled diced
    each wedge by dozens of lines: 1,600 input triangles became 32,608 and
    the union took 20 s.  Cuts now need a real triangle-triangle overlap,
    are confined to the intersection segment, and pass-through points on a
    triangle's own edges are dropped."""
    bore = mesh_only(conic(4.2, 4.2, 18.0, center=point(0, 0, -9)))
    seat = mesh_only(conic(11.15, 11.15, 7.2, center=point(0, 0, 0.8)))
    union = native(bore, seat, "union")
    inputs = sum(len(s[3]) for s in bore[1] + seat[1])
    out = sum(len(s[3]) for s in union[1])
    assert issolidclosed(union)
    identity = (signed_volume(bore) + signed_volume(seat)
                - signed_volume(native(bore, seat, "intersection")))
    assert signed_volume(union) == pytest.approx(identity, rel=1e-9)
    assert out < 3 * inputs, f"{inputs} triangles in, {out} out"


def test_identical_solids():
    """Every face coincides with a SAME face: the hardest coplanar case."""
    a, b = mesh_only(prism(2, 2, 2)), mesh_only(prism(2, 2, 2))
    assert_correct(native(a, b, "union"), 8.0)
    assert_correct(native(a, b, "intersection"), 8.0)
    assert native(a, b, "difference")[1] == []


def test_touching_at_an_edge_only():
    """Two cubes sharing only an edge: their union is non-manifold by
    definition -- four faces meet at that edge -- so yapCAD's issolidclosed
    rightly says no.  What must hold is the volume, and that the only edges
    not used exactly twice are the four-way ones on the shared line."""
    from collections import Counter
    a = mesh_only(prism(2, 2, 2))
    b = moved(prism(2, 2, 2), (2.0, 2.0, 0.0))
    union = native(a, b, "union")
    assert signed_volume(union) == pytest.approx(16.0, rel=1e-12)

    def key(p):
        return tuple(round(c / 5e-6) for c in p[:3])

    counts = Counter()
    for surf in union[1]:
        for i, j, k in surf[3]:
            ks = [key(surf[1][n]) for n in (i, j, k)]
            for x, y in ((0, 1), (1, 2), (2, 0)):
                counts[tuple(sorted((ks[x], ks[y])))] += 1
    odd = {e: c for e, c in counts.items() if c != 2}
    assert set(odd.values()) <= {4}
    shared_x, shared_y = round(1.0 / 5e-6), round(1.0 / 5e-6)
    assert all(p[0] == shared_x and p[1] == shared_y
               for e in odd for p in e)
    assert native(a, b, "intersection")[1] == []


def test_the_legacy_keywords_are_accepted():
    a, b = mesh_only(prism(4, 4, 4)), mesh_only(prism(2, 2, 2))
    assert_correct(solid_boolean(a, b, "difference", tol=1e-6, stitch=True,
                                 engine="native"), 56.0)


# ---------------------------------------------------------------------------
# The core's building blocks, and the bugs found while building it
# ---------------------------------------------------------------------------


class TestTriangulateConvex:
    """Each polygon must be covered, keep every vertex, and emit no sliver."""

    @staticmethod
    def check(poly, expected_area):
        import numpy as np
        from yapcad.boolean.csg import triangulate_convex
        pts = [np.asarray(p, dtype=float) for p in poly]
        out, tris = triangulate_convex(pts, np.array([0.0, 0.0, 1.0]), 1e-9)
        used = {i for t in tris for i in t}
        assert used >= set(range(len(poly))), "a vertex was dropped"
        area = 0.0
        for i, j, k in tris:
            twice = np.cross(out[j] - out[i], out[k] - out[i])[2]
            assert twice > 0.0, "a degenerate or flipped triangle"
            area += twice / 2.0
        assert area == pytest.approx(expected_area, rel=1e-12)
        return out, tris

    def test_a_chord_must_not_lie_along_a_straight_chain(self):
        """The first triangulator clipped an ear whose chord ran along a
        straight run of vertices, dropping the collinear remainder."""
        poly = [(0, 0, 0), (1, 0, 0), (2, 0, 0), (2, 1, 0), (0.0, 1.0, 0)]
        self.check(poly, 2.0)

    def test_points_on_every_side(self):
        poly = [(0, 0, 0), (1, 0, 0), (2, 0, 0), (1, 1, 0), (0, 2, 0),
                (0, 1, 0)]
        self.check(poly, 2.0)

    def test_a_sliver_adds_no_point(self):
        """The second triangulator fanned around a centroid, which on a
        sliver sat within tolerance of the long edge and was invisible to
        the T-junction repair."""
        poly = [(0, 0, 0), (1, 0, 0), (2, 0, 0), (2, 1e-6, 0),
                (1, 1e-6, 0), (0, 1e-6, 0)]
        out, _ = self.check(poly, 2e-6)
        assert len(out) == len(poly)

    def test_near_coincident_corners_keep_their_area(self):
        """Two corners a hair apart -- one intersection point computed
        twice -- each read as straight, so the polygon seemed to have two
        corners and the third triangulator dropped it whole (a hole of
        area 2e-6 in an icosahedron minus sphere)."""
        d = 1.2e-9 / math.sqrt(2.0)
        poly = [(0, 0, 0), (1, 0, 0), (1 - d, d, 0), (0, 1, 0)]
        twice = sum(poly[i - 1][0] * poly[i][1] - poly[i][0] * poly[i - 1][1]
                    for i in range(len(poly)))
        self.check(poly, twice / 2.0)

    def test_a_collinear_polygon_has_no_triangles(self):
        import numpy as np
        from yapcad.boolean.csg import triangulate_convex
        pts = [np.array(p, dtype=float) for p in
               ((0, 0, 0), (1, 0, 0), (2, 0, 0))]
        assert triangulate_convex(pts, np.array([0.0, 0.0, 1.0]),
                                  1e-9)[1] == []


def test_winding_numbers_distinguish_inside_from_outside():
    import numpy as np
    from yapcad.boolean.csg import _Mesh, solid_to_mesh, winding_numbers
    box = _Mesh(*solid_to_mesh(mesh_only(prism(2, 2, 2))))
    w = winding_numbers(np.array([[0.0, 0.0, 0.0], [0.99, 0.5, -0.5],
                                  [1.01, 0.0, 0.0], [5.0, 5.0, 5.0]]), box)
    assert w[0] == pytest.approx(1.0) and w[1] == pytest.approx(1.0)
    assert abs(w[2]) < 1e-9 and abs(w[3]) < 1e-9


def test_t_junction_repair_closes_a_cracked_mesh():
    """A cube face split in one triangle but not its neighbour: closed by
    index only once the shared vertex is inserted."""
    import numpy as np
    from collections import Counter
    from yapcad.boolean.csg import _repair_t_junctions, solid_to_mesh
    V, F = solid_to_mesh(mesh_only(prism(2, 2, 2)))
    a, b = F[0][0], F[0][1]
    mid = len(V)
    V = np.vstack([V, (V[a] + V[b]) / 2.0])
    c = F[0][2]
    F = np.vstack([[a, mid, c], [mid, b, c], F[1:]])
    V2, F2 = _repair_t_junctions(V, F, 1e-9, 1e-9)
    counts = Counter()
    for t in F2.tolist():
        for e in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            counts[(min(e), max(e))] += 1
    assert set(counts.values()) == {2}
