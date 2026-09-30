"""Tests for exact CSG replay of SDF trees through OCC (SDF plan, Phase 3).

§5.2's claim is that a tree of analytic primitives, hard booleans and
similarity transforms is isomorphic to a CSG tree, so replaying it through
OCC yields an exact BREP.  "Exact" is tested literally: replayed volumes are
compared with closed-form volumes to near machine precision, not with the
dual-contoured mesh.
"""

import json
import math
import os
import tempfile

import pytest

from yapcad import sdf
from yapcad.brep import occ_available

pytestmark = pytest.mark.skipif(not occ_available(),
                                reason="pythonocc-core is not available")

if occ_available():
    from yapcad.sdf import occ


def rounded_box_volume(a, b, c, r):
    """Minkowski volume: core box, face slabs, edge cylinders, corner balls."""
    A, B, C = a - 2 * r, b - 2 * r, c - 2 * r
    return (A * B * C + 2 * r * (A * B + B * C + A * C)
            + math.pi * r * r * (A + B + C) + 4.0 / 3.0 * math.pi * r ** 3)


EXACT = [
    ("sphere", lambda: sdf.sphere(5.0), 4.0 / 3.0 * math.pi * 125.0),
    ("box", lambda: sdf.box((8.0, 6.0, 4.0)), 192.0),
    ("rounded_box", lambda: sdf.rounded_box((10.0, 6.0, 4.0), 1.0),
     rounded_box_volume(10.0, 6.0, 4.0, 1.0)),
    ("cylinder", lambda: sdf.cylinder(4.0, 10.0), math.pi * 160.0),
    ("cone", lambda: sdf.cone(3.0, 1.0, 10.0), math.pi * 10.0 / 3.0 * 13.0),
    ("pointed cone", lambda: sdf.cone(3.0, 0.0, 10.0), math.pi * 30.0),
    ("cone with equal radii", lambda: sdf.cone(3.0, 3.0, 10.0),
     math.pi * 90.0),
    ("torus", lambda: sdf.torus(9.0, 2.5), 2.0 * math.pi ** 2 * 9.0 * 6.25),
    ("capsule", lambda: sdf.capsule((-3, 0, -2), (4, 1, 3), 2.5),
     math.pi * 6.25 * math.sqrt(75.0) + 4.0 / 3.0 * math.pi * 2.5 ** 3),
    ("degenerate capsule", lambda: sdf.capsule((1, 2, 3), (1, 2, 3), 2.0),
     4.0 / 3.0 * math.pi * 8.0),
    ("cube frame", lambda: sdf.subtract(sdf.box(10.0), sdf.sphere(6.0)),
     1000.0 - 254.0 * math.pi),
    ("union of disjoint", lambda: sdf.union(
        sdf.translate(sdf.box(2.0), (-5, 0, 0)),
        sdf.translate(sdf.box(2.0), (5, 0, 0))), 16.0),
    ("intersection", lambda: sdf.intersect(sdf.box(10.0),
                                           sdf.translate(sdf.box(10.0),
                                                         (5, 5, 5))),
     125.0),
    ("rotated box", lambda: sdf.rotate(
        sdf.translate(sdf.box((8.0, 6.0, 4.0)), (1, 2, 3)), (0, 0, 1), 30.0),
     192.0),
    ("uniformly scaled", lambda: sdf.scale(sdf.sphere(1.0), 3.0),
     4.0 / 3.0 * math.pi * 27.0),
    ("mirrored", lambda: sdf.scale(sdf.cone(3.0, 1.0, 10.0), (-1, -1, -1)),
     math.pi * 10.0 / 3.0 * 13.0),
    ("half-space cut", lambda: sdf.intersect(
        sdf.box(10.0), sdf.half_space((0, 0, 1), 2.0)), 700.0),
    ("tilted half-space", lambda: sdf.intersect(
        sdf.sphere(5.0), sdf.half_space((1, 1, 1), 0.0)),
     4.0 / 3.0 * math.pi * 125.0 / 2.0),
    ("half-space removed", lambda: sdf.subtract(
        sdf.box(10.0), sdf.half_space((0, 0, 1), 2.0)), 300.0),
    ("half-space missing the part", lambda: sdf.intersect(
        sdf.box(10.0), sdf.half_space((0, 0, 1), 50.0)), 1000.0),
    ("half-space under a transform", lambda: sdf.translate(sdf.rotate(
        sdf.intersect(sdf.box(10.0), sdf.half_space((0, 0, 1), 0.0)),
        (1, 0, 0), 90.0), (40, 0, 0)), 500.0),
]


class TestExactReplay:

    @pytest.mark.parametrize("name,build,expected", EXACT,
                             ids=[e[0] for e in EXACT])
    def test_replayed_volume_is_exact(self, name, build, expected):
        shape = occ.to_occ_shape(build())
        assert occ.volume(shape) == pytest.approx(expected, rel=1e-6)

    @pytest.mark.parametrize("name,build,expected", EXACT,
                             ids=[e[0] for e in EXACT])
    def test_replayed_shape_fits_the_field_bounds(self, name, build,
                                                  expected):
        node = build()
        lo, hi = occ.bounding_box(occ.to_occ_shape(node))
        diag = math.dist(*node.bounds)
        for i in range(3):
            assert lo[i] >= node.bounds[0][i] - 1e-6 * diag - 1e-6
            assert hi[i] <= node.bounds[1][i] + 1e-6 * diag + 1e-6

    def test_every_csg_capable_kind_has_an_occ_backend(self):
        """A new analytic kind must come with a replay, or be classified
        non-CSG; either way it cannot silently fail at export time."""
        from yapcad.sdf.node import get_spec, registered_kinds
        from test_sdf import SAMPLES
        for kind in registered_kinds():
            if SAMPLES[kind].csg_exact:
                assert "occ" in get_spec(kind).backends, kind

    def test_the_replay_agrees_with_the_mesh(self):
        """Different routes to the same part should meet: the exact BREP and
        the dual-contoured preview of one tree."""
        node = sdf.subtract(sdf.rounded_box((20.0, 20.0, 8.0), 2.0),
                            sdf.cylinder(4.0, 20.0))
        from yapcad.geom3d import volumeof
        meshed = volumeof(sdf.to_solid(node, resolution=64))
        exact = occ.volume(occ.to_occ_shape(node))
        assert meshed == pytest.approx(exact, rel=5e-3)

    def test_a_shared_subtree_is_replayed_once(self):
        from yapcad.sdf.node import get_spec
        spec = get_spec("sphere")
        original = spec.backends["occ"]
        calls = {"n": 0}

        def counting(*args):
            calls["n"] += 1
            return original(*args)

        spec.backends["occ"] = counting
        try:
            ball = sdf.sphere(2.0)
            occ.to_occ_shape(sdf.union(ball, sdf.intersect(ball,
                                                           sdf.box(3.0))))
        finally:
            spec.backends["occ"] = original
        assert calls["n"] == 1

    def test_the_containment_guard_catches_a_leaking_replay(self):
        """to_occ_shape checks the result against the field's bounds -- the
        one assumption the half-space stand-in rests on.  Make a backend
        misbehave and the guard must fire rather than hand back a bad part."""
        from yapcad.sdf.node import get_spec
        spec = get_spec("sphere")
        original = spec.backends["occ"]
        spec.backends["occ"] = lambda p, c, b, r: original(
            {"radius": p["radius"] * 3.0}, c, b, r)
        try:
            with pytest.raises(sdf.SdfError, match="beyond the SDF tree"):
                occ.to_occ_shape(sdf.sphere(2.0))
        finally:
            spec.backends["occ"] = original


class TestRefusal:

    def test_a_field_only_tree_is_refused_naming_the_cause(self):
        node = sdf.subtract(
            sdf.box(20.0),
            sdf.smooth_union(sdf.sphere(4.0), sdf.cylinder(2.0, 30.0), 1.0))
        with pytest.raises(sdf.SdfError) as excinfo:
            occ.to_occ_shape(node)
        assert "smooth_union" in str(excinfo.value)
        assert "faceted" in str(excinfo.value)

    def test_an_unbounded_tree_is_refused(self):
        with pytest.raises(sdf.SdfError, match="unbounded"):
            occ.to_occ_shape(sdf.half_space((0, 0, 1)))

    def test_a_spindle_torus_is_not_offered_for_replay(self):
        """OCC builds one, calls it valid, and counts the self-overlap twice
        in its volume; the field is the union.  So the classifier says no."""
        assert sdf.torus(9.0, 2.5).csg_exact
        assert not sdf.torus(2.0, 3.0).csg_exact
        assert not sdf.torus(3.0, 3.0).csg_exact

    def test_a_fully_rounded_edge_is_not_offered_for_replay(self):
        """At radius = half the smallest edge the fillet has no face left."""
        assert sdf.rounded_box((10.0, 6.0, 4.0), 1.9).csg_exact
        assert not sdf.rounded_box((10.0, 6.0, 4.0), 2.0).csg_exact


class TestBlockers:

    def test_a_csg_exact_tree_has_none(self):
        assert sdf.csg_blockers(
            sdf.subtract(sdf.box(10.0), sdf.sphere(4.0))) == []

    def test_every_independent_cause_is_named(self):
        """A smooth blend whose operand is a lattice is two reasons, and
        fixing only one would leave the part unreplayable."""
        lattice = sdf.intersect(sdf.box(20.0), sdf.gyroid(4.0, 0.3))
        node = sdf.smooth_union(sdf.sphere(4.0), lattice, 1.0)
        kinds = [n.kind for n in sdf.csg_blockers(node)]
        assert kinds == ["gyroid", "smooth_union"]

    def test_ancestors_that_merely_inherit_are_not_named(self):
        node = sdf.subtract(sdf.box(20.0),
                            sdf.offset(sdf.sphere(4.0), 1.0))
        assert [n.kind for n in sdf.csg_blockers(node)] == ["offset"]

    def test_a_nonuniform_scale_is_a_blocker(self):
        node = sdf.union(sdf.box(4.0),
                         sdf.scale(sdf.sphere(1.0), (3.0, 1.0, 1.0)))
        assert [n.kind for n in sdf.csg_blockers(node)] == ["transform"]

    def test_the_description_is_readable(self):
        text = sdf.describe_blockers(
            sdf.subtract(sdf.box(20.0), sdf.shell(sdf.sphere(5.0), 1.0)))
        assert "shell(thickness=1.0)" in text


class TestToSolidWithBrep:

    def test_brep_true_attaches_a_tagged_derived_brep(self):
        from yapcad.metadata import get_solid_metadata
        node = sdf.subtract(sdf.box(10.0), sdf.sphere(4.0))
        solid = sdf.to_solid(node, resolution=12, brep=True)
        record = get_solid_metadata(solid)["brep"]
        assert record["derivedFrom"] == {"sdf": node.digest}

    def test_the_mesh_is_still_the_sdf_preview(self):
        """brep=True adds a representation; it does not swap the preview."""
        node = sdf.sphere(5.0)
        plain = sdf.to_solid(node, resolution=16)
        with_brep = sdf.to_solid(node, resolution=16, brep=True)
        assert plain[1][0][1] == with_brep[1][0][1]

    def test_brep_true_refuses_a_field_only_tree(self):
        with pytest.raises(sdf.SdfError, match="smooth_union"):
            sdf.to_solid(sdf.smooth_union(sdf.sphere(4.0), sdf.box(5.0), 1.0),
                         resolution=12, brep=True)

    def test_brep_auto_attaches_when_it_can(self):
        from yapcad.brep import has_brep_data
        solid = sdf.to_solid(sdf.box(8.0), resolution=12, brep="auto")
        assert has_brep_data(solid)

    def test_brep_auto_skips_quietly_when_it_cannot(self):
        from yapcad.brep import has_brep_data
        solid = sdf.to_solid(
            sdf.smooth_union(sdf.sphere(4.0), sdf.box(5.0), 1.0),
            resolution=12, brep="auto")
        assert not has_brep_data(solid)

    def test_a_bad_brep_argument_is_refused(self):
        with pytest.raises(sdf.SdfError, match="brep must be"):
            sdf.to_solid(sdf.box(8.0), resolution=12, brep="yes")


def _document(solids):
    from yapcad.io.geometry_json import geometry_to_json
    return json.loads(json.dumps(geometry_to_json(solids, units="mm")))


def _solid_entry(document):
    return next(e for e in document["entities"] if e["type"] == "solid")


class TestDerivedBrepSerialisation:

    def test_the_brep_is_emitted_as_derived(self):
        node = sdf.subtract(sdf.box(10.0), sdf.sphere(4.0))
        reps = _solid_entry(_document(
            [sdf.to_solid(node, resolution=12, brep=True)]))["representations"]
        assert reps["authoritative"] == "sdf"
        assert reps["brep"]["role"] == "derived"
        assert reps["brep"]["derivedFrom"] == "sdf"
        assert reps["brep"]["treeDigest"] == node.digest

    def test_it_round_trips_with_the_brep_intact(self):
        from yapcad.brep import brep_from_solid
        from yapcad.construction import get_construction
        from yapcad.io.geometry_json import geometry_from_json
        node = sdf.subtract(sdf.rounded_box((20.0, 20.0, 8.0), 2.0),
                            sdf.cylinder(4.0, 20.0))
        restored = geometry_from_json(_document(
            [sdf.to_solid(node, resolution=12, brep=True)]))[0]
        assert sdf.from_construction(get_construction(restored)) == node
        assert occ.volume(brep_from_solid(restored).shape) == pytest.approx(
            occ.volume(occ.to_occ_shape(node)), rel=1e-9)

    def test_it_validates_against_the_published_schema(self):
        jsonschema = pytest.importorskip("jsonschema")
        from pathlib import Path
        schema = json.loads(
            (Path(__file__).resolve().parents[1] / "docs" / "schemas"
             / "yapcad-geometry-json-v0.3.schema.json").read_text())
        document = _document([sdf.to_solid(sdf.box(8.0), 12, brep=True)])
        jsonschema.Draft202012Validator(schema).validate(document)

    def test_a_brep_from_a_different_tree_is_rejected_on_read(self):
        from yapcad.io.geometry_json import geometry_from_json
        document = _document([sdf.to_solid(sdf.box(8.0), 12, brep=True)])
        _solid_entry(document)["representations"]["brep"]["treeDigest"] = \
            sdf.box(9.0).digest
        with pytest.raises(ValueError, match="different SDF tree"):
            geometry_from_json(document)

    def test_the_digest_is_checked_against_the_rebuilt_tree(self):
        """Relabelling the tree's root key and the BREP's digest together
        must not fool the check: it compares with the tree as rebuilt."""
        from yapcad.io.geometry_json import geometry_from_json
        document = _document([sdf.to_solid(sdf.box(8.0), 12, brep=True)])
        reps = _solid_entry(document)["representations"]
        tree = reps["sdf"]["tree"]
        tree["nodes"]["nforged"] = tree["nodes"].pop(tree["root"])
        tree["root"] = "nforged"
        reps["brep"]["treeDigest"] = "nforged"
        with pytest.raises(ValueError, match="different SDF tree"):
            geometry_from_json(document)

    def test_a_derived_brep_claiming_authority_is_rejected(self):
        from yapcad.io.geometry_json import geometry_from_json
        document = _document([sdf.to_solid(sdf.box(8.0), 12, brep=True)])
        _solid_entry(document)["representations"]["brep"]["role"] = \
            "authoritative"
        with pytest.raises(ValueError, match="role must be 'derived'"):
            geometry_from_json(document)

    def test_a_derived_brep_needs_an_sdf_to_derive_from(self):
        from yapcad.io.geometry_json import geometry_from_json
        document = _document([sdf.to_solid(sdf.box(8.0), 12, brep=True)])
        reps = _solid_entry(document)["representations"]
        del reps["sdf"]
        reps["authoritative"] = "mesh"
        reps["mesh"]["role"] = "authoritative"
        with pytest.raises(ValueError, match="non-authoritative BREP"):
            geometry_from_json(document)

    def test_an_untagged_brep_beside_a_tree_is_still_refused(self):
        """Phase 2's guard survives: a BREP not replayed from the tree --
        here prism()'s own analytic one -- is a contradiction."""
        from yapcad.construction import set_construction
        from yapcad.geom3d_util import prism
        solid = prism(10.0, 10.0, 10.0)
        set_construction(solid, sdf.to_construction(sdf.box(10.0)))
        with pytest.raises(ValueError, match="not replayed from its SDF"):
            _document([solid])


class TestStepExport:

    def test_a_replayed_part_exports_analytic_step(self):
        from yapcad.io.step import write_step_analytic
        node = sdf.subtract(
            sdf.rounded_box((40.0, 40.0, 10.0), 3.0),
            sdf.translate(sdf.cylinder(6.0, 20.0), (10.0, 0.0, 0.0)))
        solid = sdf.to_solid(node, resolution=16, brep=True)
        handle, path = tempfile.mkstemp(suffix=".step")
        os.close(handle)
        try:
            assert write_step_analytic(solid, path,
                                       fallback_to_faceted=False) is True
            text = open(path).read()
        finally:
            os.unlink(path)
        # The bore and the edge fillets are cylinders, the fillet corners
        # spheres, the faces planes -- analytic surfaces, not facets.
        assert "CYLINDRICAL_SURFACE" in text
        assert "SPHERICAL_SURFACE" in text
        assert "PLANE(" in text
        assert "TRIANGULATED" not in text
