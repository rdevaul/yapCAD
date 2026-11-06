"""Building a design without meshing it, and the refusals that makes safe.

Meshing is about seven eighths of the time to build a field-authored part,
and a design bound for a range-of-motion sweep, a render or a clearance
check never needs the triangles.  Deferring the mesh is therefore mostly a
performance feature -- but it is a *correctness* feature that it refuses
loudly, because an empty surface list is already legal in yapCAD and the
code that consumes solids was written when the only way to get one was to
mean it.  ``issolidclosed`` calls an empty solid closed and ``volumeof``
calls it zero; both are right about nothing at all and wrong about a part
whose mesh has not been built yet.

So the tests here come in two halves: that a deferred part is the same part
once realized, and that everything which would otherwise answer from an
empty surface list says so instead.
"""

import math

import numpy as np
import pytest

from yapcad import sdf
from yapcad.construction import get_construction
from yapcad.geom3d import (
    issolidclosed,
    rotatesolid,
    solid,
    solid_shells,
    solidbbox,
    translatesolid,
    volumeof,
)
from yapcad.geom import point
from yapcad.mesh import mesh_view
from yapcad.sdf.deferred import DeferredMeshError

BLOCK = sdf.box((20.0, 16.0, 10.0))
BORED = sdf.subtract(BLOCK, sdf.cylinder(3.0, 40.0))


def deferred(node=BLOCK, **kwargs):
    with sdf.deferring():
        return sdf.to_solid(node, **kwargs)


class TestDeferring:

    def test_nothing_is_meshed_inside_the_block(self):
        part = deferred()
        assert part[1] == []
        assert sdf.is_deferred(part)

    def test_the_block_only_covers_itself(self):
        inside = deferred()
        outside = sdf.to_solid(BLOCK, resolution=16)
        assert sdf.is_deferred(inside)
        assert not sdf.is_deferred(outside)

    def test_the_flag_overrides_the_block_either_way(self):
        with sdf.deferring():
            assert not sdf.is_deferred(sdf.to_solid(BLOCK, resolution=16,
                                                    defer=False))
        assert sdf.is_deferred(sdf.to_solid(BLOCK, defer=True))

    def test_blocks_nest_and_unwind(self):
        assert not sdf.is_deferring()
        with sdf.deferring():
            assert sdf.is_deferring()
            with sdf.deferring(False):
                assert not sdf.is_deferring()
            assert sdf.is_deferring()
        assert not sdf.is_deferring()

    def test_the_field_comes_through_intact(self):
        part = deferred(BORED)
        assert sdf.is_sdf_solid(part)
        assert sdf.from_construction(get_construction(part)).digest \
            == BORED.digest

    def test_the_plan_records_what_it_would_have_meshed(self):
        plan = sdf.deferred_plan(deferred(resolution=48))
        assert plan["method"] == sdf.DEFERRED
        assert plan["resolution"] == 48

    def test_a_cell_hint_becomes_a_resolution_in_the_plan(self):
        # Recorded in cells, not millimetres: cells survive a transform.
        plan = sdf.deferred_plan(deferred(cell_mm=1.0))
        assert plan["resolution"] == 20

    def test_explicit_bounds_are_recorded(self):
        box = ((-30.0, -30.0, -30.0), (30.0, 30.0, 30.0))
        plan = sdf.deferred_plan(deferred(bounds=box, resolution=8))
        assert plan["bounds"] == [-30.0] * 3 + [30.0] * 3

    def test_asking_for_a_brep_and_no_work_is_refused(self):
        with pytest.raises(sdf.SdfError, match="opposite things"):
            deferred(brep=True)

    def test_an_automatic_brep_is_simply_skipped(self):
        assert sdf.is_deferred(deferred(brep="auto"))

    def test_metadata_survives(self):
        from yapcad.metadata import get_solid_metadata
        part = deferred(metadata={"component": {"id": "PLATE"}})
        assert get_solid_metadata(part)["component"]["id"] == "PLATE"


class TestIsDeferred:

    def test_a_meshed_field_is_not_deferred(self):
        assert not sdf.is_deferred(sdf.to_solid(BLOCK, resolution=16))

    def test_a_genuinely_empty_solid_is_not_deferred(self):
        # An empty CSG result means it: it is not waiting for anything.
        # (``solid([])`` is the spelling; bare ``solid()`` raises.)
        assert not sdf.is_deferred(solid([]))

    def test_something_that_is_not_a_solid_is_not_deferred(self):
        assert not sdf.is_deferred(None)
        assert not sdf.is_deferred(42)

    def test_a_meshed_solid_has_no_plan(self):
        assert sdf.deferred_plan(sdf.to_solid(BLOCK, resolution=16)) is None


class TestRealize:

    def test_it_meshes_in_place_so_other_references_follow(self):
        part = deferred(resolution=24)
        alias = part
        sdf.realize(part)
        assert alias is part
        assert alias[1] and not sdf.is_deferred(alias)

    def test_a_realized_part_is_the_part(self):
        eager = sdf.to_solid(BORED, resolution=32)
        lazy = deferred(BORED, resolution=32)
        sdf.realize(lazy)
        assert volumeof(lazy) == pytest.approx(volumeof(eager), rel=1e-12)

    def test_realizing_twice_is_harmless(self):
        part = deferred(resolution=24)
        sdf.realize(part)
        before = len(part[1])
        sdf.realize(part)
        assert len(part[1]) == before

    def test_realizing_something_already_meshed_leaves_it_alone(self):
        part = sdf.to_solid(BLOCK, resolution=16)
        assert sdf.realize(part) is part

    def test_the_plan_can_be_overridden(self):
        part = deferred(resolution=16)
        sdf.realize(part, resolution=48)
        assert sdf.meshing_from_construction(
            get_construction(part))["resolution"] == 48

    def test_realize_all_counts_what_it_built(self):
        parts = [deferred(resolution=16), deferred(BORED, resolution=16),
                 sdf.to_solid(BLOCK, resolution=16)]
        assert sdf.realize_all(parts) == 2
        assert not any(sdf.is_deferred(p) for p in parts)

    def test_realize_all_meshes_a_shared_part_once(self):
        shared = deferred(resolution=16)
        assert sdf.realize_all([shared, shared, [shared]]) == 1

    def test_realize_all_walks_a_mapping(self):
        assert sdf.realize_all({"a": deferred(resolution=16),
                                "b": deferred(BORED, resolution=16)}) == 2

    def test_realize_all_walks_an_assemblys_geometry(self):
        class FakeAssembly:
            geometry = {"plate": deferred(resolution=16)}

        assembly = FakeAssembly()
        assert sdf.realize_all(assembly) == 1
        assert not sdf.is_deferred(assembly.geometry["plate"])

    def test_a_deferred_marker_with_no_field_is_refused(self, monkeypatch):
        # Defensive: a record that says deferred but yields no tree should
        # say so rather than fail somewhere further down.
        part = deferred(resolution=16)
        monkeypatch.setattr("yapcad.sdf.deferred.from_construction",
                            lambda record: None)
        with pytest.raises(DeferredMeshError, match="no field"):
            sdf.realize(part)


class TestRefusals:
    """The half that makes deferral safe rather than merely fast."""

    def test_volume_refuses_rather_than_reporting_zero(self):
        with pytest.raises(DeferredMeshError, match="realize"):
            volumeof(deferred())

    def test_closure_refuses_rather_than_reporting_closed(self):
        with pytest.raises(DeferredMeshError, match="realize"):
            issolidclosed(deferred())

    def test_shells_refuse_rather_than_reporting_none(self):
        with pytest.raises(DeferredMeshError, match="realize"):
            solid_shells(deferred())

    def test_the_triangle_tap_refuses(self):
        with pytest.raises(DeferredMeshError, match="realize"):
            list(mesh_view(deferred()))

    def test_the_refusal_names_the_operation_and_the_cure(self):
        with pytest.raises(DeferredMeshError) as caught:
            volumeof(deferred())
        message = str(caught.value)
        assert "volumeof" in message
        assert "yapcad.sdf.realize" in message

    def test_a_genuinely_empty_solid_still_answers_as_it_always_did(self):
        empty = solid([])
        assert issolidclosed(empty) is True
        assert volumeof(empty) == 0.0
        assert solid_shells(empty) == ([], [])

    def test_a_meshed_field_answers_normally(self):
        part = sdf.to_solid(BLOCK, resolution=24)
        assert issolidclosed(part)
        assert volumeof(part) > 0.0


class TestBounds:

    def test_the_box_comes_from_the_field_with_no_mesh(self):
        box = solidbbox(deferred())
        assert box[0][:3] == pytest.approx((-10.0, -8.0, -5.0))
        assert box[1][:3] == pytest.approx((10.0, 8.0, 5.0))

    def test_it_agrees_with_the_meshed_part(self):
        eager = solidbbox(sdf.to_solid(BORED, resolution=48))
        lazy = solidbbox(deferred(BORED))
        assert lazy[0][:3] == pytest.approx(eager[0][:3], abs=0.5)
        assert lazy[1][:3] == pytest.approx(eager[1][:3], abs=0.5)

    def test_an_unbounded_deferred_field_says_so(self):
        # Returning an empty box would be read as "no extent" rather than
        # "could not say", and this is the one solid with no mesh to fall
        # back on.
        with pytest.raises(ValueError, match="unbounded"):
            solidbbox(deferred(sdf.half_space((0.0, 0.0, 1.0), 0.0),
                               bounds=((-1.0,) * 3, (1.0,) * 3)))


class TestTransforms:
    """A placed part has to carry its field to where it was placed."""

    def test_translating_moves_the_field(self):
        moved = translatesolid(deferred(), [25.0, 0.0, 0.0, 1.0])
        assert solidbbox(moved)[0][0] == pytest.approx(15.0)
        assert solidbbox(moved)[1][0] == pytest.approx(35.0)

    def test_translating_keeps_the_plan(self):
        # The record is a plan, not a claim about a mesh that exists, so it
        # survives the move -- a placed part is still the same number of
        # cells across.
        moved = translatesolid(deferred(resolution=37), [5.0, 0.0, 0.0, 1.0])
        assert sdf.is_deferred(moved)
        assert sdf.deferred_plan(moved)["resolution"] == 37

    def test_a_meshed_part_still_drops_its_meshing_record(self):
        moved = translatesolid(sdf.to_solid(BLOCK, resolution=16),
                               [5.0, 0.0, 0.0, 1.0])
        assert sdf.meshing_from_construction(get_construction(moved)) is None

    def test_rotating_moves_the_field(self):
        turned = rotatesolid(deferred(), 90.0)
        box = solidbbox(turned)
        assert box[1][0] == pytest.approx(8.0, abs=1e-6)
        assert box[1][1] == pytest.approx(10.0, abs=1e-6)

    def test_a_placed_part_realizes_where_it_was_placed(self):
        moved = translatesolid(deferred(resolution=24), [25.0, 0.0, 0.0, 1.0])
        sdf.realize(moved)
        assert solidbbox(moved)[0][0] == pytest.approx(15.0, abs=0.6)


class TestFieldPredicates:
    """Questions the field answers exactly, and a mesh answers approximately."""

    def test_containment_is_answered_from_the_field(self):
        from yapcad.geom3d import solid_contains_point
        part = deferred()
        assert solid_contains_point(part, point(0.0, 0.0, 0.0))
        assert not solid_contains_point(part, point(50.0, 0.0, 0.0))

    def test_containment_agrees_with_the_meshed_part(self):
        from yapcad.geom3d import solid_contains_point
        lazy, eager = deferred(BORED), sdf.to_solid(BORED, resolution=48)
        for p in (point(0.0, 0.0, 0.0), point(8.0, 6.0, 2.0),
                  point(40.0, 0.0, 0.0)):
            assert solid_contains_point(lazy, p) == \
                solid_contains_point(eager, p), p

    def test_overlapping_parts_are_reported_as_overlapping(self):
        from yapcad.boolean.native import solids_intersect
        here = deferred()
        there = translatesolid(deferred(), [5.0, 0.0, 0.0, 1.0])
        assert solids_intersect(here, there)

    def test_separated_parts_are_not(self):
        from yapcad.boolean.native import solids_intersect
        here = deferred()
        far = translatesolid(deferred(), [500.0, 0.0, 0.0, 1.0])
        assert not solids_intersect(here, far)


class TestBooleans:

    def test_a_boolean_of_deferred_operands_stays_deferred(self):
        with sdf.deferring():
            a = sdf.to_solid(BLOCK, resolution=32)
            b = sdf.to_solid(sdf.cylinder(3.0, 40.0), resolution=32)
            cut = sdf.combine_all([a, b], "difference")
        assert sdf.is_deferred(cut)
        assert cut[1] == []

    def test_the_result_is_the_right_shape_once_realized(self):
        with sdf.deferring():
            a = sdf.to_solid(BLOCK, resolution=32)
            b = sdf.to_solid(sdf.cylinder(3.0, 40.0), resolution=32)
            cut = sdf.combine_all([a, b], "difference")
        sdf.realize(cut)
        truth = 20.0 * 16.0 * 10.0 - math.pi * 9.0 * 10.0
        assert volumeof(cut) == pytest.approx(truth, rel=0.05)

    def test_the_recorded_resolution_still_sizes_the_result(self):
        # Deferred operands carry no mesh to measure, so the plan is the
        # only thing that keeps a boolean from falling back to the default.
        with sdf.deferring():
            a = sdf.to_solid(BLOCK, cell_mm=0.5)
            b = sdf.to_solid(sdf.cylinder(3.0, 40.0), cell_mm=0.5)
            cut = sdf.combine_all([a, b], "difference")
        assert sdf.deferred_plan(cut)["resolution"] >= 40


class TestComposedFields:
    """A tree composed from deferred solids has to carry the plan with it.

    A compound of deferred parts has no surfaces to show for itself.  If the
    new construction record drops the plan, the result is indistinguishable
    from an empty solid -- which is the exact silent answer this whole
    feature exists to prevent, reintroduced one layer up.
    """

    def test_a_compound_of_deferred_parts_stays_marked(self):
        from yapcad.sdf.deferred import plan_from_sources
        a = deferred(resolution=32)
        b = deferred(BORED, resolution=32)
        tree = sdf.compound(
            sdf.from_construction(get_construction(a)),
            sdf.from_construction(get_construction(b)))
        plan = plan_from_sources(tree, [a, b])
        assert plan is not None
        assert plan["method"] == sdf.DEFERRED

    def test_the_plan_is_rescaled_to_the_composed_extent(self):
        # Resolutions are per-solid, so they cannot be carried across to a
        # tree spanning more space; the cell size is what is preserved.
        from yapcad.sdf.deferred import plan_from_sources
        near = deferred(resolution=20)                      # 20 mm / 20 = 1 mm
        far = sdf.translate(BLOCK, (80.0, 0.0, 0.0))
        with sdf.deferring():
            away = sdf.to_solid(far, resolution=20)
        tree = sdf.compound(BLOCK, far)                     # 100 mm across
        plan = plan_from_sources(tree, [near, away])
        assert plan["resolution"] == pytest.approx(100, abs=2)

    def test_nothing_deferred_means_no_plan(self):
        from yapcad.sdf.deferred import plan_from_sources
        meshed = sdf.to_solid(BLOCK, resolution=16)
        assert plan_from_sources(BLOCK, [meshed]) is None

    def test_a_surface_free_field_with_no_record_reads_as_deferred(self):
        # The backstop. A composed record that forgot its plan must not pass
        # for an empty solid.
        from yapcad.sdf.node import to_construction
        naked = solid([], [], to_construction(BLOCK))
        assert sdf.is_deferred(naked)
        with pytest.raises(DeferredMeshError):
            volumeof(naked)

    def test_a_field_meshed_to_nothing_is_not_deferred(self):
        # An empty CSG result says so by carrying the record that produced
        # it, and keeps answering as an empty solid.
        # Cut the whole block away: the bounds survive, the material does
        # not.  (Disjoint operands give *empty bounds*, which to_solid
        # refuses outright -- a different case.)
        nothing = sdf.to_solid(sdf.subtract(BLOCK, sdf.box((60.0, 60.0, 60.0))),
                               resolution=16)
        assert nothing[1] == []
        assert not sdf.is_deferred(nothing)
        assert volumeof(nothing) == 0.0


class TestFingerprint:

    def test_two_different_parts_do_not_share_a_fingerprint(self):
        # Hashing only the triangles would give every mesh-free part the
        # same hash, and a package build dedupes on it.
        from yapcad.package.core import _geometry_fingerprint
        assert _geometry_fingerprint(deferred(BLOCK)) \
            != _geometry_fingerprint(deferred(BORED))

    def test_the_same_part_twice_does(self):
        from yapcad.package.core import _geometry_fingerprint
        assert _geometry_fingerprint(deferred(BORED)) \
            == _geometry_fingerprint(deferred(BORED))


DESIGN = """module deferprobe
command PLATE() -> solid:
    let slab: solid = box(40.0, 30.0, 6.0)
    let hole: solid = translate(cylinder(4.0, 20.0), 10.0, 0.0, -10.0)
    emit difference(slab, hole)
"""


def run(mesh):
    from yapcad.dsl import compile_and_run
    result = compile_and_run(DESIGN, "PLATE", {}, representation="sdf",
                             sdf_cell_mm=1.0, sdf_simplify=False,
                             sdf_mesh=mesh)
    assert result.success, result.error_message
    geometry = result.geometry
    if isinstance(geometry, list) and geometry and geometry[0] != "solid":
        geometry = geometry[0]
    return geometry


class TestFromTheDsl:

    def test_a_design_builds_without_meshing(self):
        part = run(mesh=False)
        assert sdf.is_deferred(part)
        assert part[1] == []

    def test_it_is_the_same_part_as_the_meshed_build(self):
        lazy, eager = run(mesh=False), run(mesh=True)
        sdf.realize(lazy)
        assert volumeof(lazy) == pytest.approx(volumeof(eager), rel=1e-9)

    def test_meshing_is_still_the_default(self):
        assert not sdf.is_deferred(run(mesh=None))

    def test_the_field_is_what_the_renderer_needs(self):
        # The payoff: a part can be looked at without ever being meshed.
        from yapcad import render
        image = render.image(run(mesh=False), view="top", pixels=(64, 48),
                             annotate=False)
        assert image.shape == (48, 64, 3)

    def test_and_what_a_clearance_question_needs(self):
        node = sdf.from_construction(get_construction(run(mesh=False)))
        centre = sdf.evaluate(node, np.array([[0.0, 0.0, 0.0]]))
        assert float(centre[0]) < 0.0          # inside the plate
