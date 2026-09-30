"""Tests for the SDF node DAG and numpy evaluator (design document Phase 1).

The design document asks for two kinds of validation: against closed-form
analytic distances, and against yapCAD's own
:func:`yapcad.geom3d.signedFaceDistance`.  Both are here, alongside tests of
the property that the document says must be right from the first commit —
Lipschitz tracking — and of the DAG serialisation the whole representation
rests on.
"""

import json
import math

import numpy as np
import pytest

from yapcad import sdf
from yapcad.geom3d import signedFaceDistance
from yapcad.geom3d_util import prism


RNG = np.random.default_rng(20260928)


def random_points(n, extent=30.0):
    """Uniform points in a cube of half-width ``extent``."""
    return RNG.uniform(-extent, extent, size=(n, 3))


# ---------------------------------------------------------------------------
# One representative node per registered kind
# ---------------------------------------------------------------------------
#
# Keyed by kind so that the coverage test below fails the moment a new node
# kind is registered without a Lipschitz check, which is the safeguard the
# design document asks for in §4.3.

SAMPLES = {
    "sphere": sdf.sphere(7.0),
    "box": sdf.box((8.0, 12.0, 5.0)),
    "rounded_box": sdf.rounded_box((8.0, 12.0, 5.0), 1.5),
    "cylinder": sdf.cylinder(4.0, 11.0),
    "capsule": sdf.capsule((-3.0, 0.0, -2.0), (4.0, 1.0, 3.0), 2.5),
    "torus": sdf.torus(9.0, 2.5),
    "cone": sdf.cone(6.0, 2.0, 10.0),
    "half_space": sdf.half_space((0.3, 0.5, 0.8), 2.0),
    "gyroid": sdf.gyroid(12.0, 0.4),
    "schwarz_p": sdf.schwarz_p(15.0, 0.5),
    "union": sdf.union(sdf.sphere(6.0), sdf.box(8.0)),
    "intersect": sdf.intersect(sdf.sphere(6.0), sdf.box(8.0)),
    "subtract": sdf.subtract(sdf.box(9.0), sdf.sphere(5.0)),
    "smooth_union": sdf.smooth_union(sdf.sphere(6.0), sdf.box(8.0), 2.0),
    "smooth_intersect": sdf.smooth_intersect(
        sdf.sphere(6.0), sdf.box(8.0), 2.0
    ),
    "smooth_subtract": sdf.smooth_subtract(sdf.box(9.0), sdf.sphere(5.0), 1.5),
    "offset": sdf.offset(sdf.box(8.0), 1.25),
    "shell": sdf.shell(sdf.sphere(6.0), 1.0),
    "transform": sdf.rotate(sdf.box((8.0, 4.0, 2.0)), (0, 0, 1), 30.0),
}


def test_every_registered_kind_has_a_sample():
    """A new kind must come with a sample, so property tests see it."""
    assert set(SAMPLES) == set(sdf.registered_kinds())


# ---------------------------------------------------------------------------
# Analytic exactness
# ---------------------------------------------------------------------------


class TestAnalyticDistances:
    """Primitives claiming ``exact`` must match closed-form distance."""

    def test_sphere(self):
        node = sdf.sphere(7.0)
        pts = random_points(400)
        expected = np.linalg.norm(pts, axis=1) - 7.0
        assert np.allclose(sdf.evaluate(node, pts), expected)

    def test_half_space(self):
        node = sdf.half_space((0.0, 0.0, 2.0), 3.0)
        pts = random_points(400)
        # The constructor normalises the normal, so the plane is z == 3.
        assert np.allclose(sdf.evaluate(node, pts), pts[:, 2] - 3.0)

    def test_torus(self):
        node = sdf.torus(9.0, 2.5)
        pts = random_points(400)
        radial = np.hypot(pts[:, 0], pts[:, 1]) - 9.0
        expected = np.hypot(radial, pts[:, 2]) - 2.5
        assert np.allclose(sdf.evaluate(node, pts), expected)

    def test_capsule_matches_segment_distance(self):
        a = np.array([-3.0, 0.0, -2.0])
        b = np.array([4.0, 1.0, 3.0])
        node = sdf.capsule(a, b, 2.5)
        pts = random_points(400)
        ba = b - a
        t = np.clip(((pts - a) @ ba) / (ba @ ba), 0.0, 1.0)
        expected = np.linalg.norm(pts - (a + t[:, None] * ba), axis=1) - 2.5
        assert np.allclose(sdf.evaluate(node, pts), expected)

    def test_box_exterior_matches_clamped_distance(self):
        """Outside a box, distance is the norm of the axis overshoot."""
        node = sdf.box((8.0, 12.0, 5.0))
        half = np.array([4.0, 6.0, 2.5])
        pts = random_points(600)
        outside = np.any(np.abs(pts) > half, axis=1)
        pts = pts[outside]
        expected = np.linalg.norm(
            np.maximum(np.abs(pts) - half, 0.0), axis=1
        )
        assert np.allclose(sdf.evaluate(node, pts), expected)

    def test_box_interior_is_distance_to_nearest_face(self):
        node = sdf.box((8.0, 12.0, 5.0))
        half = np.array([4.0, 6.0, 2.5])
        pts = RNG.uniform(-1.0, 1.0, size=(400, 3)) * half * 0.95
        expected = -np.min(half - np.abs(pts), axis=1)
        assert np.allclose(sdf.evaluate(node, pts), expected)

    def test_cylinder_surface_is_zero(self):
        node = sdf.cylinder(4.0, 11.0)
        theta = np.linspace(0.0, 2.0 * math.pi, 64, endpoint=False)
        z = RNG.uniform(-5.4, 5.4, size=theta.shape)
        on_wall = np.stack(
            [4.0 * np.cos(theta), 4.0 * np.sin(theta), z], axis=1
        )
        assert np.allclose(sdf.evaluate(node, on_wall), 0.0, atol=1e-12)

    def test_cone_degenerates_to_cylinder_when_radii_match(self):
        cone = sdf.cone(4.0, 4.0, 11.0)
        cyl = sdf.cylinder(4.0, 11.0)
        pts = random_points(400)
        assert np.allclose(
            sdf.evaluate(cone, pts), sdf.evaluate(cyl, pts), atol=1e-9
        )

    def test_exactness_flags_match_the_claim(self):
        for kind in ("sphere", "box", "rounded_box", "cylinder", "capsule",
                     "torus", "cone", "half_space"):
            assert SAMPLES[kind].exact, f"{kind} should be exact"
        for kind in ("gyroid", "schwarz_p", "union", "intersect", "subtract",
                     "smooth_union", "smooth_intersect", "smooth_subtract"):
            assert not SAMPLES[kind].exact, f"{kind} should not be exact"


# ---------------------------------------------------------------------------
# Cross-validation against yapCAD's own triangle distance
# ---------------------------------------------------------------------------


def solid_triangles(solid):
    """Collect every triangle of a yapCAD solid as a list of 3-point lists."""
    tris = []
    for surface in solid[1]:
        verts = surface[1]
        for face in surface[3]:
            tris.append([verts[i] for i in face[:3]])
    return tris


class TestAgainstSignedFaceDistance:
    """The exact box field must agree with yapCAD's triangle distance code."""

    def test_box_matches_mesh_face_distances(self):
        solid = prism(10.0, 16.0, 6.0)
        tris = solid_triangles(solid)
        verts = np.array([v[:3] for tri in tris for v in tri])
        lo = verts.min(axis=0)
        hi = verts.max(axis=0)
        # Derive the box from the mesh rather than assuming prism's axis
        # order, so the test validates the field and not our guess at it.
        assert np.allclose(lo, -hi), "prism is expected to be origin-centred"
        node = sdf.box(tuple(hi - lo))

        pts = random_points(150, extent=14.0)
        field = np.abs(sdf.evaluate(node, pts))
        for p, expected in zip(pts, field):
            probe = [float(p[0]), float(p[1]), float(p[2]), 1.0]
            nearest = min(
                abs(signedFaceDistance(probe, tri)) for tri in tris
            )
            assert nearest == pytest.approx(expected, abs=1e-9)


# ---------------------------------------------------------------------------
# Lipschitz discipline
# ---------------------------------------------------------------------------


class TestLipschitz:
    """``|f(p)| <= L * d(p)`` is the contract; these test what implies it."""

    @pytest.mark.parametrize("kind", sorted(SAMPLES))
    def test_pairwise_lipschitz_inequality_holds(self, kind):
        """The definition, tested directly: ``|f(a)-f(b)| <= L |a-b|``."""
        node = SAMPLES[kind]
        a = random_points(3000, extent=20.0)
        b = a + RNG.normal(scale=2.0, size=a.shape)
        fa = sdf.evaluate(node, a)
        fb = sdf.evaluate(node, b)
        separation = np.linalg.norm(a - b, axis=1)
        allowed = node.lipschitz * separation
        assert np.all(np.abs(fa - fb) <= allowed + 1e-9)

    @pytest.mark.parametrize("kind", sorted(SAMPLES))
    def test_gradient_magnitude_respects_the_constant(self, kind):
        node = SAMPLES[kind]
        pts = random_points(500, extent=20.0)
        grad = sdf.gradient(node, pts, eps=1e-6)
        assert np.all(np.linalg.norm(grad, axis=1) <= node.lipschitz + 1e-4)

    def test_exact_primitives_have_unit_gradient(self):
        """An exact distance field has ``|grad f| == 1`` almost everywhere."""
        for kind in ("sphere", "torus", "capsule", "half_space"):
            node = SAMPLES[kind]
            pts = random_points(300, extent=20.0)
            mag = np.linalg.norm(sdf.gradient(node, pts, eps=1e-6), axis=1)
            assert np.allclose(mag, 1.0, atol=1e-4), kind

    def test_hard_csg_keeps_the_unit_constant(self):
        """min/max of 1-Lipschitz fields stays 1-Lipschitz, if not exact."""
        for kind in ("union", "intersect", "subtract"):
            assert SAMPLES[kind].lipschitz == 1.0

    def test_smooth_blend_preserves_the_constant(self):
        """The polynomial blend's gradient is a convex combination, so L is
        unchanged even though exactness is destroyed."""
        node = SAMPLES["smooth_union"]
        assert node.lipschitz == 1.0
        assert not node.exact

    def test_periodic_surfaces_report_a_constant_above_one(self):
        """The whole point of tracking L: these fields are not metric."""
        w = 2.0 * math.pi / 12.0
        assert SAMPLES["gyroid"].lipschitz == pytest.approx(
            2.0 * math.sqrt(3.0) * w
        )
        assert SAMPLES["gyroid"].lipschitz > 1.0
        w_p = 2.0 * math.pi / 15.0
        assert SAMPLES["schwarz_p"].lipschitz == pytest.approx(
            math.sqrt(3.0) * w_p
        )

    def test_lipschitz_propagates_up_through_combinators(self):
        coarse = sdf.gyroid(12.0, 0.4)
        combined = sdf.union(sdf.sphere(5.0), coarse)
        assert combined.lipschitz == coarse.lipschitz

    def test_safe_step_is_never_longer_than_true_distance(self):
        """A sphere tracer stepping ``|f|/L`` must not tunnel."""
        node = SAMPLES["gyroid"]
        pts = random_points(400, extent=15.0)
        values = sdf.evaluate(node, pts)
        step = sdf.safe_step(node, values)
        assert np.all(step <= np.abs(values) + 1e-12)
        assert np.all(step >= 0.0)


# ---------------------------------------------------------------------------
# Boolean semantics and bounds
# ---------------------------------------------------------------------------


class TestBooleans:

    def test_union_is_inside_when_either_operand_is(self):
        a = sdf.translate(sdf.sphere(5.0), (-3.0, 0.0, 0.0))
        b = sdf.translate(sdf.sphere(5.0), (3.0, 0.0, 0.0))
        node = sdf.union(a, b)
        pts = random_points(800, extent=12.0)
        inside_either = (sdf.evaluate(a, pts) < 0) | (sdf.evaluate(b, pts) < 0)
        assert np.array_equal(sdf.evaluate(node, pts) < 0, inside_either)

    def test_intersect_is_inside_only_when_both_are(self):
        a = sdf.translate(sdf.sphere(5.0), (-3.0, 0.0, 0.0))
        b = sdf.translate(sdf.sphere(5.0), (3.0, 0.0, 0.0))
        node = sdf.intersect(a, b)
        pts = random_points(800, extent=12.0)
        both = (sdf.evaluate(a, pts) < 0) & (sdf.evaluate(b, pts) < 0)
        assert np.array_equal(sdf.evaluate(node, pts) < 0, both)

    def test_subtract_removes_the_tool(self):
        target = sdf.box(10.0)
        tool = sdf.sphere(3.0)
        node = sdf.subtract(target, tool)
        pts = random_points(800, extent=8.0)
        expected = ((sdf.evaluate(target, pts) < 0)
                    & (sdf.evaluate(tool, pts) > 0))
        assert np.array_equal(sdf.evaluate(node, pts) < 0, expected)

    def test_single_operand_union_collapses_to_its_child(self):
        """Canonicalisation, so the DAG does not grow pointless wrappers."""
        child = sdf.sphere(4.0)
        assert sdf.union(child) is child
        assert sdf.intersect(child) is child

    def test_nary_union_matches_repeated_binary_union(self):
        a, b, c = sdf.sphere(4.0), sdf.box(6.0), sdf.cylinder(2.0, 9.0)
        pts = random_points(400)
        assert np.allclose(
            sdf.evaluate(sdf.union(a, b, c), pts),
            sdf.evaluate(sdf.union(sdf.union(a, b), c), pts),
        )

    def test_smooth_union_encloses_the_hard_union(self):
        """A fillet only ever adds material, so its field is never larger."""
        a = sdf.translate(sdf.sphere(5.0), (-4.0, 0.0, 0.0))
        b = sdf.translate(sdf.sphere(5.0), (4.0, 0.0, 0.0))
        pts = random_points(800, extent=12.0)
        hard = sdf.evaluate(sdf.union(a, b), pts)
        soft = sdf.evaluate(sdf.smooth_union(a, b, 2.0), pts)
        assert np.all(soft <= hard + 1e-12)

    def test_smooth_union_converges_to_hard_union_as_radius_shrinks(self):
        a = sdf.translate(sdf.sphere(5.0), (-4.0, 0.0, 0.0))
        b = sdf.translate(sdf.sphere(5.0), (4.0, 0.0, 0.0))
        pts = random_points(400, extent=12.0)
        hard = sdf.evaluate(sdf.union(a, b), pts)
        soft = sdf.evaluate(sdf.smooth_union(a, b, 1e-6), pts)
        assert np.allclose(soft, hard, atol=1e-6)

    def test_shell_is_thin_around_the_original_surface(self):
        node = sdf.shell(sdf.sphere(6.0), 1.0)
        # A point on the original surface sits mid-wall, half a thickness in.
        assert sdf.evaluate(node, [6.0, 0.0, 0.0]) == pytest.approx(-0.5)
        # The wall's own faces are at radius 5.5 and 6.5.
        assert sdf.evaluate(node, [5.5, 0.0, 0.0]) == pytest.approx(0.0)
        assert sdf.evaluate(node, [6.5, 0.0, 0.0]) == pytest.approx(0.0)
        assert sdf.evaluate(node, [0.0, 0.0, 0.0]) == pytest.approx(5.5)

    def test_offset_shifts_the_level_set(self):
        node = sdf.offset(sdf.sphere(5.0), 2.0)
        pts = random_points(300)
        assert np.allclose(
            sdf.evaluate(node, pts), sdf.evaluate(sdf.sphere(7.0), pts)
        )


class TestBounds:

    def test_primitive_bounds_are_tight(self):
        assert sdf.sphere(3.0).bounds == ((-3.0, -3.0, -3.0), (3.0, 3.0, 3.0))
        assert sdf.box((4.0, 6.0, 8.0)).bounds == (
            (-2.0, -3.0, -4.0), (2.0, 3.0, 4.0)
        )
        assert sdf.torus(10.0, 2.0).bounds == (
            (-12.0, -12.0, -2.0), (12.0, 12.0, 2.0)
        )

    def test_union_bounds_cover_both_operands(self):
        node = sdf.union(
            sdf.translate(sdf.sphere(2.0), (-10.0, 0.0, 0.0)),
            sdf.translate(sdf.sphere(2.0), (10.0, 0.0, 0.0)),
        )
        lo, hi = node.bounds
        assert lo[0] == pytest.approx(-12.0)
        assert hi[0] == pytest.approx(12.0)

    def test_subtract_keeps_the_target_bounds(self):
        node = sdf.subtract(sdf.box(10.0), sdf.sphere(50.0))
        assert node.bounds == sdf.box(10.0).bounds

    def test_periodic_fields_are_unbounded(self):
        assert sdf.gyroid(10.0, 0.5).bounds == sdf.UNBOUNDED
        assert sdf.half_space((0, 0, 1)).bounds == sdf.UNBOUNDED

    def test_intersecting_a_lattice_with_a_solid_bounds_it(self):
        """The documented way to bound an infinite lattice."""
        node = sdf.intersect(sdf.box(20.0), sdf.gyroid(5.0, 0.4))
        assert node.bounds == sdf.box(20.0).bounds

    def test_surface_stays_within_reported_bounds(self):
        for kind in ("sphere", "box", "rounded_box", "cylinder", "capsule",
                     "torus", "cone", "union", "intersect", "smooth_union"):
            node = SAMPLES[kind]
            lo, hi = node.bounds
            pts = random_points(2000, extent=40.0)
            outside_box = np.any(
                (pts < np.array(lo)) | (pts > np.array(hi)), axis=1
            )
            # Nothing outside the bounding box may be inside the solid.
            assert np.all(sdf.evaluate(node, pts)[outside_box] > 0), kind


# ---------------------------------------------------------------------------
# Transforms
# ---------------------------------------------------------------------------


class TestTransform:

    def test_translation_moves_the_field(self):
        node = sdf.translate(sdf.sphere(5.0), (10.0, 0.0, 0.0))
        assert sdf.evaluate(node, [10.0, 0.0, 0.0]) == pytest.approx(-5.0)
        assert sdf.evaluate(node, [15.0, 0.0, 0.0]) == pytest.approx(0.0)
        assert node.exact
        assert node.bounds == ((5.0, -5.0, -5.0), (15.0, 5.0, 5.0))

    def test_rotation_preserves_exactness_and_distance(self):
        base = sdf.box((10.0, 2.0, 2.0))
        node = sdf.rotate(base, (0, 0, 1), 90.0)
        assert node.exact
        # The long axis now runs along y.
        assert sdf.evaluate(node, [0.0, 5.0, 0.0]) == pytest.approx(
            0.0, abs=1e-12
        )
        assert sdf.evaluate(node, [5.0, 0.0, 0.0]) == pytest.approx(4.0)

    def test_uniform_scale_stays_exact(self):
        node = sdf.scale(sdf.sphere(1.0), 5.0)
        pts = random_points(300)
        assert node.exact
        assert np.allclose(
            sdf.evaluate(node, pts), sdf.evaluate(sdf.sphere(5.0), pts)
        )

    def test_nonuniform_scale_is_marked_inexact(self):
        node = sdf.scale(sdf.sphere(1.0), (3.0, 1.0, 1.0))
        assert not node.exact
        assert not node.csg_exact

    def test_nonuniform_scale_remains_conservative(self):
        """It may underestimate the distance, but it must never overestimate:
        that is what keeps a ray marcher safe on a stretched primitive."""
        node = sdf.scale(sdf.sphere(1.0), (3.0, 1.0, 1.0))
        # Dense sample of the true zero set: x^2/9 + y^2 + z^2 = 1.
        u = RNG.normal(size=(20000, 3))
        u /= np.linalg.norm(u, axis=1)[:, None]
        surface = u * np.array([3.0, 1.0, 1.0])

        pts = random_points(200, extent=5.0)
        values = sdf.evaluate(node, pts)
        true_distance = np.min(
            np.linalg.norm(pts[:, None, :] - surface[None, :, :], axis=2),
            axis=1,
        )
        assert np.all(np.abs(values) <= true_distance + 1e-9)

    def test_rotation_matches_yapcad_rotation_matrix(self):
        from yapcad.xform import Rotation
        base = sdf.box((10.0, 2.0, 2.0))
        via_helper = sdf.rotate(base, (0, 0, 1), 37.0)
        via_matrix = sdf.transform(base, Rotation([0, 0, 1, 1], 37.0))
        pts = random_points(300)
        assert np.allclose(
            sdf.evaluate(via_helper, pts), sdf.evaluate(via_matrix, pts)
        )

    def test_singular_matrix_is_rejected(self):
        flat = (
            (1.0, 0.0, 0.0, 0.0),
            (0.0, 1.0, 0.0, 0.0),
            (0.0, 0.0, 0.0, 0.0),
            (0.0, 0.0, 0.0, 1.0),
        )
        node = sdf.transform(sdf.sphere(1.0), flat)
        with pytest.raises(sdf.SdfError, match="singular"):
            _ = node.exact

    def test_projective_matrix_is_rejected(self):
        with pytest.raises(sdf.SdfError, match="affine"):
            sdf.transform(sdf.sphere(1.0), (
                (1.0, 0.0, 0.0, 0.0),
                (0.0, 1.0, 0.0, 0.0),
                (0.0, 0.0, 1.0, 0.0),
                (0.1, 0.0, 0.0, 1.0),
            ))


# ---------------------------------------------------------------------------
# CSG-exactness classification (prerequisite for the Phase 3 OCC replay)
# ---------------------------------------------------------------------------


class TestCsgExactness:
    """Structural, and deliberately distinct from field exactness."""

    def test_analytic_primitives_are_csg_exact(self):
        for kind in ("sphere", "box", "rounded_box", "cylinder", "capsule",
                     "torus", "cone", "half_space"):
            assert SAMPLES[kind].csg_exact, kind

    def test_hard_booleans_of_primitives_stay_csg_exact(self):
        for kind in ("union", "intersect", "subtract"):
            assert SAMPLES[kind].csg_exact, kind

    def test_a_hard_boolean_is_csg_exact_but_not_field_exact(self):
        """The two flags mean different things; this is the case that shows
        it, and the case §5.2 relies on for STEP export."""
        node = sdf.union(sdf.sphere(4.0), sdf.box(6.0))
        assert node.csg_exact
        assert not node.exact

    def test_smooth_blends_are_not_csg_exact(self):
        for kind in ("smooth_union", "smooth_intersect", "smooth_subtract"):
            assert not SAMPLES[kind].csg_exact, kind

    def test_periodic_surfaces_are_not_csg_exact(self):
        assert not SAMPLES["gyroid"].csg_exact
        assert not SAMPLES["schwarz_p"].csg_exact

    def test_one_field_only_node_taints_the_whole_tree(self):
        """This is what lets authoring time warn 'will not round-trip to
        STEP' rather than discovering it at export."""
        node = sdf.subtract(
            sdf.box(20.0),
            sdf.intersect(sdf.sphere(8.0), sdf.gyroid(4.0, 0.3)),
        )
        assert not node.csg_exact

    def test_rigid_transforms_preserve_csg_exactness(self):
        moved = sdf.translate(sdf.box(5.0), (1, 2, 3))
        node = sdf.rotate(moved, (0, 1, 0), 45)
        assert node.csg_exact

    def test_offset_and_shell_are_not_csg_exact(self):
        assert not SAMPLES["offset"].csg_exact
        assert not SAMPLES["shell"].csg_exact


# ---------------------------------------------------------------------------
# The DAG and its serialisation
# ---------------------------------------------------------------------------


class TestDag:

    def test_equal_subtrees_are_the_same_node(self):
        assert sdf.sphere(5.0) == sdf.sphere(5.0)
        assert sdf.sphere(5.0).digest == sdf.sphere(5.0).digest
        assert sdf.sphere(5.0).digest != sdf.sphere(5.1).digest

    def test_walk_visits_each_distinct_node_once(self):
        shared = sdf.cylinder(2.0, 10.0)
        node = sdf.union(shared, sdf.translate(shared, (5.0, 0.0, 0.0)))
        kinds = sorted(n.kind for n in sdf.walk(node))
        # cylinder appears once despite two references to it.
        assert kinds == ["cylinder", "transform", "union"]

    def test_shared_subtree_is_stored_once(self):
        shared = sdf.union(sdf.sphere(3.0), sdf.box(4.0))
        node = sdf.union(shared, sdf.translate(shared, (10.0, 0.0, 0.0)))
        doc = sdf.tree_to_json(node)
        # sphere, box, inner union, transform, outer union: five, not seven.
        assert len(doc["nodes"]) == 5

    def test_roundtrip_preserves_the_tree_and_its_field(self):
        node = sdf.subtract(
            sdf.rounded_box((40.0, 40.0, 10.0), 3.0),
            sdf.translate(sdf.cylinder(6.0, 20.0), (10.0, 0.0, 0.0)),
        )
        restored = sdf.tree_from_json(
            json.loads(json.dumps(sdf.tree_to_json(node)))
        )
        assert restored == node
        assert restored.digest == node.digest
        pts = random_points(300)
        assert np.allclose(
            sdf.evaluate(restored, pts), sdf.evaluate(node, pts)
        )

    @pytest.mark.parametrize("kind", sorted(SAMPLES))
    def test_every_kind_roundtrips(self, kind):
        node = SAMPLES[kind]
        restored = sdf.tree_from_json(
            json.loads(json.dumps(sdf.tree_to_json(node)))
        )
        assert restored == node

    def test_document_records_the_format_tag(self):
        doc = sdf.tree_to_json(sdf.sphere(1.0))
        assert doc["format"] == "yapcad-sdf-tree-v1"
        assert doc["root"] in doc["nodes"]

    def test_tree_serialises_compactly(self):
        """The size argument for a tree over a sampled grid (§4.1)."""
        node = sdf.sphere(1.0)
        for i in range(60):
            node = sdf.union(node, sdf.translate(sdf.sphere(1.0), (i, 0, 0)))
        assert len(json.dumps(sdf.tree_to_json(node))) < 32_000

    def test_unknown_format_is_rejected(self):
        doc = sdf.tree_to_json(sdf.sphere(1.0))
        doc["format"] = "yapcad-sdf-tree-v99"
        with pytest.raises(sdf.SdfError, match="unsupported SDF tree format"):
            sdf.tree_from_json(doc)

    def test_dangling_reference_is_rejected(self):
        doc = sdf.tree_to_json(sdf.union(sdf.sphere(1.0), sdf.box(2.0)))
        root = doc["nodes"][doc["root"]]
        root["children"][0] = "ndeadbeefdeadbeef"
        with pytest.raises(sdf.SdfError, match="not present"):
            sdf.tree_from_json(doc)

    def test_cycle_is_rejected(self):
        """A cycle would recurse forever in the evaluator, so it is caught
        at load rather than at evaluation."""
        doc = sdf.tree_to_json(sdf.union(sdf.sphere(1.0), sdf.box(2.0)))
        root_id = doc["root"]
        doc["nodes"][root_id]["children"][0] = root_id
        with pytest.raises(sdf.SdfError, match="cycle"):
            sdf.tree_from_json(doc)

    def test_unknown_node_kind_is_rejected(self):
        doc = sdf.tree_to_json(sdf.sphere(1.0))
        doc["nodes"][doc["root"]]["kind"] = "hyperboloid"
        with pytest.raises(sdf.SdfError, match="unknown SDF node kind"):
            sdf.tree_from_json(doc)

    def test_digest_is_stable_across_processes(self):
        """Package signing depends on this: an identifier that changed run to
        run would churn every package hash (§5.1)."""
        import subprocess
        import sys
        script = (
            "from yapcad import sdf;"
            "print(sdf.union(sdf.sphere(2.5), "
            "sdf.translate(sdf.box(3.0), (1.0, 2.0, 3.0))).digest)"
        )
        runs = {
            subprocess.run(
                [sys.executable, "-c", script],
                capture_output=True, text=True, check=True,
            ).stdout.strip()
            for _ in range(2)
        }
        assert len(runs) == 1
        assert runs.pop() == sdf.union(
            sdf.sphere(2.5), sdf.translate(sdf.box(3.0), (1.0, 2.0, 3.0))
        ).digest


# ---------------------------------------------------------------------------
# The bridge to the Phase 0 construction slot
# ---------------------------------------------------------------------------


class TestConstructionBridge:

    def test_to_and_from_construction_roundtrip(self):
        node = sdf.smooth_union(sdf.sphere(4.0), sdf.box(6.0), 1.0)
        record = sdf.to_construction(node)
        assert record[0] == "sdf"
        assert sdf.from_construction(record) == node

    def test_from_construction_ignores_other_record_kinds(self):
        assert sdf.from_construction(["procedure", "prism(2,2,2)"]) is None
        assert sdf.from_construction([]) is None
        assert sdf.from_construction(None) is None

    def test_record_is_accepted_by_the_construction_module(self):
        from yapcad.construction import SDF, normalize_construction
        record = sdf.to_construction(sdf.sphere(3.0))
        assert record[0] == SDF
        assert normalize_construction(record) == record

    def test_sdf_tree_survives_the_geometry_json_boundary(self):
        """Phase 0 made the slot round-trip; this is the payload it was for."""
        from yapcad.construction import get_construction
        from yapcad.io.geometry_json import (
            geometry_from_json,
            geometry_to_json,
        )

        node = sdf.subtract(sdf.box(10.0), sdf.sphere(4.0))
        solid = sdf.to_solid(node, resolution=12)

        document = json.loads(
            json.dumps(geometry_to_json([solid], units="mm"))
        )
        restored = geometry_from_json(document)

        recovered = sdf.from_construction(get_construction(restored[0]))
        assert recovered == node
        assert recovered.digest == node.digest

    def test_a_solid_cannot_claim_both_sdf_and_brep_authority(self):
        """prism() attaches an analytic BREP, so a prism is BREP-
        authoritative; pinning an SDF record to one is a contradiction, and
        §2 makes authority a single per-solid property."""
        from yapcad.brep import occ_available
        from yapcad.construction import set_construction
        from yapcad.io.geometry_json import geometry_to_json

        if not occ_available():
            pytest.skip("pythonocc-core is not available")

        solid = prism(10.0, 10.0, 10.0)
        set_construction(solid, sdf.to_construction(sdf.box(10.0)))
        with pytest.raises(ValueError, match="both SDF and BREP authority"):
            geometry_to_json([solid], units="mm")


# ---------------------------------------------------------------------------
# The evaluator contract
# ---------------------------------------------------------------------------


class TestEvaluatorContract:
    """Design document §3.1: arrays in, arrays out, from the first commit."""

    def test_n_by_3_in_n_out(self):
        values = sdf.evaluate(sdf.sphere(1.0), np.zeros((17, 3)))
        assert values.shape == (17,)

    def test_yapcad_homogeneous_points_are_accepted(self):
        """yapCAD points carry a w component; it must be ignored, not read
        as a coordinate."""
        plain = sdf.evaluate(sdf.sphere(5.0), [[1.0, 2.0, 3.0]])
        homogeneous = sdf.evaluate(sdf.sphere(5.0), [[1.0, 2.0, 3.0, 1.0]])
        assert plain == pytest.approx(homogeneous)

    def test_single_point_returns_a_scalar(self):
        value = sdf.evaluate(sdf.sphere(5.0), [0.0, 0.0, 0.0])
        assert isinstance(value, float)
        assert value == pytest.approx(-5.0)

    def test_empty_input_is_not_an_error(self):
        values = sdf.evaluate(sdf.sphere(1.0), np.zeros((0, 3)))
        assert values.shape == (0,)

    def test_wrong_component_count_is_rejected(self):
        with pytest.raises(sdf.SdfError, match=r"\(N, 3\) or \(N, 4\)"):
            sdf.evaluate(sdf.sphere(1.0), np.zeros((5, 2)))

    def test_shared_subtree_is_evaluated_once(self):
        """DAG memoisation, not tree expansion: the whole point of sharing."""
        calls = {"n": 0}
        base = sdf.sphere(3.0)
        spec = sdf.node.get_spec("sphere")
        original = spec.backends["numpy"]

        def counting(params, p, children, ev):
            calls["n"] += 1
            return original(params, p, children, ev)

        spec.backends["numpy"] = counting
        try:
            tree = sdf.union(base, sdf.intersect(base, sdf.box(4.0)))
            sdf.evaluate(tree, np.zeros((10, 3)))
        finally:
            spec.backends["numpy"] = original
        assert calls["n"] == 1

    def test_gradient_shape_and_direction(self):
        node = sdf.sphere(5.0)
        grad = sdf.gradient(node, [[3.0, 4.0, 0.0]])
        assert grad.shape == (1, 3)
        # On a sphere the gradient points radially outward.
        assert np.allclose(grad[0], [0.6, 0.8, 0.0], atol=1e-5)

    def test_normal_is_unit_length(self):
        node = SAMPLES["union"]
        normals = sdf.normal(node, random_points(200, extent=15.0))
        assert np.allclose(np.linalg.norm(normals, axis=1), 1.0, atol=1e-4)

    def test_normal_is_zero_where_the_gradient_vanishes(self):
        """The centre of a sphere has no defined normal; report zero rather
        than nan so a mesher can detect and skip it."""
        assert np.allclose(sdf.normal(sdf.sphere(5.0), [0.0, 0.0, 0.0]), 0.0)

    def test_sample_grid_shape_and_coordinates(self):
        node = sdf.sphere(5.0)
        values, (xs, ys, zs) = sdf.sample_grid(node, resolution=(9, 11, 13))
        assert values.shape == (9, 11, 13)
        assert xs[0] == pytest.approx(-5.0) and xs[-1] == pytest.approx(5.0)
        assert len(ys) == 11 and len(zs) == 13
        # Odd counts put a sample exactly on the centre, the deepest point.
        assert values.min() == pytest.approx(-5.0)
        assert values[4, 5, 6] == pytest.approx(-5.0)

    def test_sample_grid_uses_node_bounds_by_default(self):
        values, _ = sdf.sample_grid(sdf.box(10.0), resolution=4)
        # Every corner of the box's own bounding box lies on the surface.
        assert values[0, 0, 0] == pytest.approx(0.0)

    def test_sample_grid_padding_clears_the_surface(self):
        values, _ = sdf.sample_grid(sdf.box(10.0), resolution=4, padding=2.0)
        assert values[0, 0, 0] > 0.0

    def test_sample_grid_refuses_an_unbounded_field(self):
        with pytest.raises(sdf.SdfError, match="unbounded"):
            sdf.sample_grid(sdf.gyroid(5.0, 0.3))

    def test_sample_grid_accepts_explicit_bounds_for_a_lattice(self):
        values, _ = sdf.sample_grid(
            sdf.gyroid(5.0, 0.3),
            bounds=((-10.0, -10.0, -10.0), (10.0, 10.0, 10.0)),
            resolution=12,
        )
        assert values.shape == (12, 12, 12)
        assert values.min() < 0.0 < values.max()


# ---------------------------------------------------------------------------
# Argument validation
# ---------------------------------------------------------------------------


class TestValidation:

    @pytest.mark.parametrize("call", [
        lambda: sdf.sphere(0.0),
        lambda: sdf.sphere(-1.0),
        lambda: sdf.box((1.0, 0.0, 1.0)),
        lambda: sdf.cylinder(1.0, -2.0),
        lambda: sdf.torus(5.0, 0.0),
        lambda: sdf.cone(0.0, 0.0, 5.0),
        lambda: sdf.gyroid(0.0, 0.5),
        lambda: sdf.shell(sdf.sphere(1.0), 0.0),
        lambda: sdf.smooth_union(sdf.sphere(1.0), sdf.box(1.0), 0.0),
        lambda: sdf.half_space((0.0, 0.0, 0.0)),
        lambda: sdf.scale(sdf.sphere(1.0), 0.0),
        lambda: sdf.rotate(sdf.sphere(1.0), (0, 0, 0), 30.0),
    ])
    def test_degenerate_arguments_are_rejected(self, call):
        with pytest.raises(sdf.SdfError):
            call()

    def test_rounded_box_radius_cannot_exceed_the_box(self):
        with pytest.raises(sdf.SdfError, match="exceeds half the small"):
            sdf.rounded_box((10.0, 10.0, 2.0), 3.0)

    def test_non_finite_dimensions_are_rejected(self):
        with pytest.raises(sdf.SdfError):
            sdf.sphere(float("inf"))
        with pytest.raises(sdf.SdfError):
            sdf.sphere(float("nan"))

    def test_subtract_requires_a_tool(self):
        with pytest.raises(sdf.SdfError, match="at least one tool"):
            sdf.subtract(sdf.box(1.0))

    def test_children_must_be_nodes(self):
        with pytest.raises(sdf.SdfError, match="must be SDF nodes"):
            sdf.union(sdf.sphere(1.0), "not a node")

    def test_unsupported_parameter_type_is_rejected_not_coerced(self):
        """Unlike a provenance record, an SDF parameter is load-bearing and
        must never be silently degraded to its repr."""
        with pytest.raises(sdf.SdfError, match="JSON scalar or sequence"):
            sdf.make_node("sphere", {"radius": object()})

    def test_arity_is_enforced(self):
        with pytest.raises(sdf.SdfError, match="at most 1 children"):
            sdf.make_node("shell", {"thickness": 1.0},
                          (sdf.sphere(1.0), sdf.box(1.0)))
