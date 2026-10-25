"""Chamfers of SDF-authored parts, with and without OCC.

Supportless FDM wants bevels, not rounds, so a design that has to print
without supports is written with ``chamfer`` nearly everywhere ``fillet``
would otherwise go -- and until these primitives existed, ``chamfer``
needed a BREP, which meant no SDF-authored part could be chamfered at all.

The shape is pinned two ways: against closed-form volumes derived from the
constraint form, and against OCC's own chamfer-all-edges, which is the
thing the field claims to be a field for.  The field quality is pinned
three ways: exact inside, never longer than the true distance outside (so a
clearance read off it is conservative), and -- where OCC is available to
measure it -- within ten percent of the true distance.
"""

import math

import numpy as np
import pytest

from yapcad import sdf
from yapcad.brep import occ_available
from yapcad.construction import get_construction
from yapcad.dsl.runtime.builtins import get_builtin_registry
from yapcad.dsl.runtime.values import float_val, solid_val
from yapcad.geom3d import issolidclosed, volumeof
from yapcad.geom3d_util import prism

needs_occ = pytest.mark.skipif(not occ_available(),
                               reason="pythonocc-core is not available")


def tree_of(solid):
    return sdf.from_construction(get_construction(solid))


def chamfered_box_volume(a, b, c, d):
    """The box, less a prism along each of the twelve edges, less the
    corners, by inclusion-exclusion over the three bevel directions and the
    corner cut: ``2 d**2`` per unit of edge length, and ``16 d**3 / 3`` put
    back for the eight corners.
    """
    return a * b * c - 2.0 * d * d * (a + b + c) + 16.0 * d ** 3 / 3.0


def chamfered_cylinder_volume(radius, height, d):
    """The cylinder less the two bevel rings, each by Pappus on the
    triangle the bevel cuts out of the profile."""
    ring = math.pi * (radius ** 2 * d
                      - (radius ** 3 - (radius - d) ** 3) / 3.0)
    return math.pi * radius ** 2 * height - 2.0 * ring


# The 26 planes of a chamfered box, written out independently of the
# primitive: six faces, twelve bevels, eight corners.  Their max is the
# exact distance inside any convex polytope, and the sign everywhere.
def plane_max(pts, size, d):
    hx, hy, hz = (v * 0.5 for v in size)
    x, y, z = (np.abs(pts[:, i]) for i in range(3))
    faces = [x - hx, y - hy, z - hz]
    bevels = [(x + y - (hx + hy - d)) / math.sqrt(2.0),
              (y + z - (hy + hz - d)) / math.sqrt(2.0),
              (x + z - (hx + hz - d)) / math.sqrt(2.0)]
    corner = (x + y + z - (hx + hy + hz - 2.0 * d)) / math.sqrt(3.0)
    out = corner
    for plane in faces + bevels:
        out = np.maximum(out, plane)
    return out


BOXES = [((60.0, 40.0, 20.0), 3.0),
         ((20.0, 20.0, 20.0), 5.0),
         ((300.0, 176.0, 10.0), 1.0),
         ((100.0, 10.0, 10.0), 1.0),
         ((8.0, 8.0, 80.0), 3.0),
         ((20.0, 8.0, 8.0), 3.9)]

CYLINDERS = [(10.0, 8.0, 2.0), (18.0, 20.0, 1.5), (4.0, 10.0, 0.5)]


# ---------------------------------------------------------------------------
# the primitives, as fields
# ---------------------------------------------------------------------------


class TestChamferedBoxField:

    def test_the_declared_properties(self):
        node = sdf.chamfered_box((60, 40, 20), 3.0)
        assert node.bounds == ((-30.0, -20.0, -10.0), (30.0, 20.0, 10.0))
        assert node.lipschitz == 1.0
        assert node.csg_exact
        # An intersection of four exact fields: a lower bound outside.
        assert not node.exact

    def test_the_faces_bevels_and_corners_are_where_they_should_be(self):
        node = sdf.chamfered_box((60, 40, 20), 3.0)
        pts = np.array([
            [0, 0, 0],              # centre: half the smallest edge in
            [30, 0, 0],             # on the +x face
            [33, 0, 0],             # three out from it
            [0, 0, 10],             # on the +z face
            [28.5, 18.5, 0],        # on the xy bevel: 28.5 + 18.5 = 47
            [30, 17, 7],            # the corner vertex
            [28.5, 18.5, 8.5],      # the three bevel planes' triple point,
                                    # which the corner cut removes
        ], dtype=float)
        assert sdf.evaluate(node, pts) == pytest.approx(
            [-10.0, 0.0, 3.0, 0.0, 0.0, 0.0, 1.5 / math.sqrt(3.0)],
            abs=1e-12)

    def test_zero_distance_is_a_plain_box(self):
        assert sdf.chamfered_box((10, 10, 10), 0.0) == sdf.box(10.0)

    def test_a_fully_bevelled_edge_is_not_offered_for_replay(self):
        """At distance = half the smallest edge the bevel has no face
        left, which is where rounded_box stops as well."""
        assert sdf.chamfered_box((10.0, 6.0, 4.0), 1.9).csg_exact
        assert not sdf.chamfered_box((10.0, 6.0, 4.0), 2.0).csg_exact

    def test_too_large_a_bevel_is_refused(self):
        with pytest.raises(sdf.SdfError, match="half the smallest edge"):
            sdf.chamfered_box((10.0, 6.0, 4.0), 2.1)

    def test_a_bad_distance_is_refused(self):
        with pytest.raises(sdf.SdfError, match="non-negative"):
            sdf.chamfered_box(10.0, -1.0)
        with pytest.raises(sdf.SdfError, match="must be a number"):
            sdf.chamfered_box(10.0, "wide")

    @pytest.mark.parametrize("size,d", BOXES)
    def test_inside_it_is_the_exact_distance(self, size, d):
        """Inside a convex polytope the distance to the boundary is the
        distance to the nearest face plane, and the ball that reaches it
        stays inside, so the foot is a boundary point."""
        rng = np.random.default_rng(11)
        half = np.array(size, dtype=float) * 0.5
        pts = rng.uniform(-half, half, size=(4000, 3))
        field = sdf.evaluate(sdf.chamfered_box(size, d), pts)
        reference = plane_max(pts, size, d)
        inside = reference < 0.0
        assert inside.sum() > 1000
        assert field[inside] == pytest.approx(reference[inside], abs=1e-9)

    @pytest.mark.parametrize("size,d", BOXES)
    def test_the_sign_is_right_everywhere(self, size, d):
        rng = np.random.default_rng(12)
        half = np.array(size, dtype=float) * 0.5
        pts = rng.uniform(-1.6 * half, 1.6 * half, size=(4000, 3))
        field = sdf.evaluate(sdf.chamfered_box(size, d), pts)
        reference = plane_max(pts, size, d)
        assert np.all(np.sign(field) == np.sign(reference))

    @pytest.mark.parametrize("size,d", BOXES)
    def test_outside_it_never_over_reports(self, size, d):
        """A clearance read off this field has to be safe: no point of the
        part may lie nearer than the field says.  The witnesses are a cloud
        of points known to be inside by the 26-plane reference.
        """
        rng = np.random.default_rng(13)
        half = np.array(size, dtype=float) * 0.5
        cloud = rng.uniform(-half, half, size=(30000, 3))
        cloud = cloud[plane_max(cloud, size, d) < 0.0]
        pts = rng.uniform(-1.7 * half, 1.7 * half, size=(150, 3))
        field = sdf.evaluate(sdf.chamfered_box(size, d), pts)
        outside = field > 0.0
        nearest = np.sqrt(((pts[outside, None, :] - cloud[None, :, :]) ** 2)
                          .sum(axis=2)).min(axis=1)
        assert np.all(field[outside] <= nearest + 1e-9)

    @pytest.mark.parametrize("size,d", BOXES)
    def test_the_field_is_lipschitz_1(self, size, d):
        rng = np.random.default_rng(14)
        half = np.array(size, dtype=float) * 0.5
        a = rng.uniform(-2.0 * half, 2.0 * half, size=(2000, 3))
        b = a + rng.normal(scale=0.35 * half, size=(2000, 3))
        node = sdf.chamfered_box(size, d)
        step = np.abs(sdf.evaluate(node, a) - sdf.evaluate(node, b))
        assert np.all(step <= np.linalg.norm(a - b, axis=1) + 1e-9)


class TestChamferedCylinderField:

    def test_the_declared_properties(self):
        node = sdf.chamfered_cylinder(10.0, 8.0, 2.0)
        assert node.bounds == ((-10.0, -10.0, -4.0), (10.0, 10.0, 4.0))
        assert node.lipschitz == 1.0
        assert node.csg_exact
        # A solid of revolution: the 2D profile distance is the 3D one.
        assert node.exact

    def test_the_side_cap_and_bevel_are_where_they_should_be(self):
        node = sdf.chamfered_cylinder(10.0, 8.0, 2.0)
        pts = np.array([
            [0, 0, 0],          # the axis, half a diameter from the wall
            [10, 0, 0],         # on the side
            [12, 0, 0],         # two out from it
            [0, 0, 4],          # on the cap
            [0, 0, 6],          # two above it
            [9, 0, 3],          # on the bevel: 9 + 3 = 10 + 4 - 2
            [0, 9 + math.sqrt(0.5), 3 + math.sqrt(0.5)],    # one out along
        ], dtype=float)                                     # its normal
        assert sdf.evaluate(node, pts) == pytest.approx(
            [-4.0, 0.0, 2.0, 0.0, 2.0, 0.0, 1.0], abs=1e-12)

    def test_zero_distance_is_a_plain_cylinder(self):
        assert sdf.chamfered_cylinder(5.0, 10.0, 0.0) == sdf.cylinder(5, 10)

    def test_a_fully_bevelled_edge_is_not_offered_for_replay(self):
        assert sdf.chamfered_cylinder(10.0, 8.0, 3.9).csg_exact
        assert not sdf.chamfered_cylinder(10.0, 8.0, 4.0).csg_exact
        assert not sdf.chamfered_cylinder(4.0, 20.0, 4.0).csg_exact

    def test_too_large_a_bevel_is_refused(self):
        with pytest.raises(sdf.SdfError, match="exceeds the radius"):
            sdf.chamfered_cylinder(10.0, 8.0, 4.1)
        with pytest.raises(sdf.SdfError, match="exceeds the radius"):
            sdf.chamfered_cylinder(4.0, 20.0, 4.1)

    @pytest.mark.parametrize("radius,height,d", CYLINDERS)
    def test_the_profile_distance_is_the_exact_distance(self, radius,
                                                        height, d):
        """Every point's distance, checked against the distance to a dense
        sample of the surface of revolution -- which is an upper bound on
        the truth, and an exact field has to meet it."""
        rng = np.random.default_rng(15)
        # The profile boundary: side, both bevels, both caps.
        t = np.linspace(0.0, 1.0, 4001)
        h = 0.5 * height
        rho = np.concatenate([
            np.full_like(t, radius),                 # side
            radius - d * t, radius - d * t,          # bevels
            (radius - d) * t, (radius - d) * t,      # caps
        ])
        z = np.concatenate([
            (2.0 * t - 1.0) * (h - d),
            h - d + d * t, -(h - d + d * t),
            np.full_like(t, h), np.full_like(t, -h),
        ])
        theta = rng.uniform(0.0, 2.0 * math.pi, size=rho.shape)
        surface = np.stack([rho * np.cos(theta), rho * np.sin(theta), z], 1)
        pts = rng.uniform([-1.8 * radius] * 2 + [-1.8 * h],
                          [1.8 * radius] * 2 + [1.8 * h], size=(120, 3))
        field = sdf.evaluate(sdf.chamfered_cylinder(radius, height, d), pts)
        nearest = np.sqrt(
            ((pts[:, None, :] - surface[None, :, :]) ** 2).sum(axis=2)
        ).min(axis=1)
        # Upper bound from a finite sample: never beaten, and met closely.
        assert np.all(np.abs(field) <= nearest + 1e-9)
        assert np.all(np.abs(field) >= nearest - 0.05 * max(radius, h))

    @pytest.mark.parametrize("radius,height,d", CYLINDERS)
    def test_the_field_is_lipschitz_1(self, radius, height, d):
        rng = np.random.default_rng(16)
        box = np.array([2.0 * radius, 2.0 * radius, height])
        a = rng.uniform(-box, box, size=(2000, 3))
        b = a + rng.normal(scale=0.3 * box, size=(2000, 3))
        node = sdf.chamfered_cylinder(radius, height, d)
        step = np.abs(sdf.evaluate(node, a) - sdf.evaluate(node, b))
        assert np.all(step <= np.linalg.norm(a - b, axis=1) + 1e-9)


# ---------------------------------------------------------------------------
# the primitives, as solids
# ---------------------------------------------------------------------------


class TestVolumes:

    @pytest.mark.parametrize("size,d", BOXES[:2])
    def test_the_meshed_box_matches_the_closed_form(self, size, d):
        out = sdf.to_solid(sdf.chamfered_box(size, d), resolution=80)
        assert issolidclosed(out)
        assert volumeof(out) == pytest.approx(
            chamfered_box_volume(*size, d), rel=1e-3)

    @pytest.mark.parametrize("radius,height,d", CYLINDERS[:2])
    def test_the_meshed_cylinder_matches_the_closed_form(self, radius,
                                                         height, d):
        out = sdf.to_solid(sdf.chamfered_cylinder(radius, height, d),
                           resolution=80)
        assert issolidclosed(out)
        assert volumeof(out) == pytest.approx(
            chamfered_cylinder_volume(radius, height, d), rel=1e-3)

    @pytest.mark.parametrize("size,d", [((10, 6, 4), 2.0),
                                        ((10, 10, 10), 5.0),
                                        ((100, 4, 4), 2.0)])
    def test_the_largest_bevel_still_meshes(self, size, d):
        """At half the smallest edge the bevel eats two faces entirely --
        a cube becomes a cuboctahedron.  OCC will not build it, which is
        what ``csg_exact`` says, but the field is still a field."""
        node = sdf.chamfered_box(size, d)
        assert not node.csg_exact
        out = sdf.to_solid(node, resolution=72)
        assert issolidclosed(out)
        assert volumeof(out) == pytest.approx(
            chamfered_box_volume(*size, d), rel=6e-3)

    @needs_occ
    @pytest.mark.parametrize("size,d", BOXES)
    def test_the_replayed_box_is_occs_own_chamfer(self, size, d):
        """The whole claim of the primitive: this field describes the solid
        OCC gets by chamfering every edge of the box, to the last decimal
        of its volume."""
        from yapcad.sdf import occ
        node = sdf.chamfered_box(size, d)
        assert occ.volume(occ.to_occ_shape(node)) == pytest.approx(
            chamfered_box_volume(*size, d), rel=1e-9)

    @needs_occ
    @pytest.mark.parametrize("radius,height,d", CYLINDERS)
    def test_the_replayed_cylinder_is_occs_own_chamfer(self, radius,
                                                       height, d):
        from yapcad.sdf import occ
        node = sdf.chamfered_cylinder(radius, height, d)
        assert occ.volume(occ.to_occ_shape(node)) == pytest.approx(
            chamfered_cylinder_volume(radius, height, d), rel=1e-9)

    @needs_occ
    def test_a_chamfered_primitive_in_a_boolean_still_replays(self):
        from yapcad.sdf import occ
        node = sdf.subtract(sdf.chamfered_box((40, 40, 10), 2.0),
                            sdf.chamfered_cylinder(6.0, 20.0, 1.0))
        assert sdf.csg_blockers(node) == []
        shape = occ.to_occ_shape(node)
        # The plate, less a plain cylindrical hole: the cutter is twice the
        # plate's thickness, so its own bevels fall outside the plate and
        # the hole it leaves has a square entry.
        hole = math.pi * 36.0 * 10.0
        assert occ.volume(shape) == pytest.approx(
            chamfered_box_volume(40, 40, 10, 2.0) - hole, rel=1e-6)


@needs_occ
class TestAgainstOccDistances:
    """OCC measures the true distance from a point to the part, which is
    the only oracle that can tell an exact field from a conservative one.
    """

    @staticmethod
    def skin_of(node):
        from OCC.Core.BRep import BRep_Builder
        from OCC.Core.TopAbs import TopAbs_FACE
        from OCC.Core.TopExp import TopExp_Explorer
        from OCC.Core.TopoDS import TopoDS_Compound, topods
        from yapcad.sdf import occ
        comp = TopoDS_Compound()
        builder = BRep_Builder()
        builder.MakeCompound(comp)
        explorer = TopExp_Explorer(occ.to_occ_shape(node), TopAbs_FACE)
        while explorer.More():
            builder.Add(comp, topods.Face(explorer.Current()))
            explorer.Next()
        return comp

    @staticmethod
    def distances(shape, pts):
        from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
        from OCC.Core.BRepExtrema import BRepExtrema_DistShapeShape
        from OCC.Core.gp import gp_Pnt
        out = []
        for p in pts:
            vertex = BRepBuilderAPI_MakeVertex(
                gp_Pnt(float(p[0]), float(p[1]), float(p[2]))).Shape()
            measure = BRepExtrema_DistShapeShape(shape, vertex)
            measure.Perform()
            out.append(measure.Value())
        return np.array(out)

    @pytest.mark.parametrize("size,d", BOXES)
    def test_the_box_is_conservative_and_within_ten_percent(self, size, d):
        """The under-report is a vertex-region deficit of just under 7% of
        the true distance, whatever the proportions of the box, so a tenth
        is a bound with room in it rather than a tolerance to tune."""
        rng = np.random.default_rng(17)
        node = sdf.chamfered_box(size, d)
        skin = self.skin_of(node)
        half = np.array(size, dtype=float) * 0.5
        pts = rng.uniform(-1.6 * half, 1.6 * half, size=(200, 3))
        field = sdf.evaluate(node, pts)
        truth = self.distances(skin, pts)
        outside = field > 1e-9
        assert np.all(field[outside] <= truth[outside] + 1e-9)
        assert np.all(field[outside] >= 0.9 * truth[outside])
        inside = field < -1e-9
        assert inside.sum() > 10
        assert field[inside] == pytest.approx(-truth[inside], abs=1e-9)

    @pytest.mark.parametrize("radius,height,d", CYLINDERS)
    def test_the_cylinder_is_exact(self, radius, height, d):
        rng = np.random.default_rng(18)
        node = sdf.chamfered_cylinder(radius, height, d)
        skin = self.skin_of(node)
        extent = np.array([1.8 * radius, 1.8 * radius, 1.8 * height])
        pts = rng.uniform(-extent, extent, size=(200, 3))
        field = sdf.evaluate(node, pts)
        truth = self.distances(skin, pts)
        assert np.abs(field) == pytest.approx(truth, abs=1e-9)


class TestSerialisation:

    @pytest.mark.parametrize("node", [
        sdf.chamfered_box((60, 40, 20), 3.0),
        sdf.chamfered_cylinder(18.0, 20.0, 1.5),
    ])
    def test_a_chamfered_primitive_round_trips(self, node):
        assert sdf.tree_from_json(sdf.tree_to_json(node)) == node
        assert sdf.from_construction(sdf.to_construction(node)) == node


# ---------------------------------------------------------------------------
# chamfer(): the substitution
# ---------------------------------------------------------------------------


class TestChamfer:

    def test_a_box_becomes_a_chamfered_box(self):
        assert sdf.chamfer(sdf.box((40, 16, 28)), 3.0) == \
            sdf.chamfered_box((40, 16, 28), 3.0)

    def test_a_cylinder_becomes_a_chamfered_cylinder(self):
        assert sdf.chamfer(sdf.cylinder(18.0, 20.0), 1.5) == \
            sdf.chamfered_cylinder(18.0, 20.0, 1.5)

    def test_a_placed_primitive_keeps_its_place(self):
        """The rover's own idiom: chamfer(translate(box(...)), d)."""
        out = sdf.chamfer(sdf.translate(sdf.box((30, 10, 6)), (5, 0, 2)),
                          1.0)
        assert out == sdf.translate(sdf.chamfered_box((30, 10, 6), 1.0),
                                    (5, 0, 2))

    def test_a_uniform_scale_divides_the_distance(self):
        out = sdf.chamfer(sdf.scale(sdf.box(10.0), 2.0), 3.0)
        assert out == sdf.scale(sdf.chamfered_box(10.0, 1.5), 2.0)

    def test_a_non_uniform_scale_is_refused(self):
        with pytest.raises(sdf.SdfError, match="non-uniform scale"):
            sdf.chamfer(sdf.scale(sdf.box(10.0), (2.0, 1.0, 1.0)), 1.0)

    def test_each_body_of_a_compound_is_chamfered(self):
        pair = sdf.compound(sdf.box(10.0),
                            sdf.translate(sdf.cylinder(3.0, 10.0), (20, 0, 0)))
        assert sdf.chamfer(pair, 1.0) == sdf.compound(
            sdf.chamfered_box(10.0, 1.0),
            sdf.translate(sdf.chamfered_cylinder(3.0, 10.0, 1.0), (20, 0, 0)))

    @pytest.mark.parametrize("node", [
        sdf.sphere(5.0),
        sdf.torus(9.0, 2.5),
        sdf.capsule((-3, 0, 0), (4, 1, 3), 2.5),
        sdf.rounded_box((20, 10, 8), 2.0),
        sdf.rounded_cylinder(10.0, 8.0, 2.0),
    ])
    def test_an_edgeless_primitive_is_left_alone(self, node):
        assert sdf.chamfer(node, 1.0) is node

    def test_a_combined_field_is_refused_with_the_reason(self):
        node = sdf.subtract(sdf.box(20.0), sdf.cylinder(4.0, 30.0))
        with pytest.raises(sdf.SdfError, match="sampled-grid"):
            sdf.chamfer(node, 1.0)

    def test_chamfering_a_chamfer_is_refused(self):
        with pytest.raises(sdf.SdfError, match="already chamfered"):
            sdf.chamfer(sdf.chamfered_box(10.0, 1.0), 0.5)

    def test_a_cone_is_refused(self):
        with pytest.raises(sdf.SdfError, match="cannot bevel"):
            sdf.chamfer(sdf.cone(5.0, 2.0, 10.0), 1.0)

    def test_too_large_a_bevel_names_the_primitive_it_broke(self):
        with pytest.raises(sdf.SdfError, match="smallest edge 2 of a box"):
            sdf.chamfer(sdf.box((10, 4, 20)), 3.0)
        with pytest.raises(sdf.SdfError, match="of a cylinder"):
            sdf.chamfer(sdf.cylinder(10.0, 8.0, ), 5.0)

    @pytest.mark.parametrize("distance", [0.0, -1.0, float("inf")])
    def test_a_bad_distance_is_refused(self, distance):
        with pytest.raises(sdf.SdfError, match="positive and finite"):
            sdf.chamfer(sdf.box(10.0), distance)

    def test_a_distance_that_is_not_a_number_is_refused(self):
        with pytest.raises(sdf.SdfError, match="must be a number"):
            sdf.chamfer(sdf.box(10.0), None)


class TestChamferSolid:

    def test_it_matches_the_closed_form_volume(self):
        solid = sdf.to_solid(sdf.box((40, 16, 28)), resolution=80)
        out = sdf.chamfer_solid(solid, 3.0)
        assert issolidclosed(out)
        assert tree_of(out) == sdf.chamfered_box((40, 16, 28), 3.0)
        assert volumeof(out) == pytest.approx(
            chamfered_box_volume(40, 16, 28, 3.0), rel=2e-3)

    def test_it_refuses_a_mesh_solid(self):
        with pytest.raises(sdf.SdfError, match="not SDF-authored"):
            sdf.chamfer_solid(prism(1, 1, 1), 0.1)


# ---------------------------------------------------------------------------
# the DSL builtin
# ---------------------------------------------------------------------------


def builtin(name):
    return get_builtin_registry().get_function(name).implementation


def test_dsl_chamfer_of_an_sdf_solid_needs_no_occ(monkeypatch):
    import yapcad.brep
    monkeypatch.setattr(yapcad.brep, "occ_available", lambda: False)
    solid = sdf.to_solid(sdf.cylinder(18.0, 20.0), resolution=48)
    out = builtin("chamfer")(solid_val(solid), float_val(1.5)).data
    assert tree_of(out) == sdf.chamfered_cylinder(18.0, 20.0, 1.5)
    assert issolidclosed(out)


def test_dsl_chamfer_of_a_placed_box_needs_no_occ(monkeypatch):
    import yapcad.brep
    monkeypatch.setattr(yapcad.brep, "occ_available", lambda: False)
    node = sdf.translate(sdf.box((300, 176, 10)), (0, 0, 5))
    solid = sdf.to_solid(node, resolution=64)
    out = builtin("chamfer")(solid_val(solid), float_val(1.0)).data
    assert tree_of(out) == sdf.translate(
        sdf.chamfered_box((300, 176, 10), 1.0), (0, 0, 5))


def test_dsl_chamfer_of_a_mesh_solid_explains_itself(monkeypatch):
    import yapcad.brep
    monkeypatch.setattr(yapcad.brep, "occ_available", lambda: False)
    with pytest.raises(RuntimeError, match="mesh-only solid"):
        builtin("chamfer")(solid_val(prism(1, 1, 1)), float_val(0.1))
