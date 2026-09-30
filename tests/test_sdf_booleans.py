"""Booleans between SDF-authored solids stay fields (SDF plan, Phase 4).

Handed two dual-contoured meshes, the native mesh engine returned an open
mesh even for two boxes sharing a face, took tens of seconds at modest
resolution, and at higher resolution timed out or raised on a degenerate
face.  Every yapCAD boolean goes through ``geom3d.solid_boolean``, the DSL's
included, so these tests pin that SDF pairs no longer go near it.
"""

import json
import math

import pytest

from yapcad import sdf
from yapcad.brep import occ_available
from yapcad.construction import get_construction
from yapcad.geom3d import issolid, issolidclosed, solid_boolean, solidbbox, \
    translatesolid, volumeof
from yapcad.geom import point
from yapcad.geom3d_util import prism
from yapcad.sdf import booleans

needs_occ = pytest.mark.skipif(not occ_available(),
                               reason="pythonocc-core is not available")


def tree_of(solid):
    return sdf.from_construction(get_construction(solid))


def meshing_of(solid):
    return sdf.meshing_from_construction(get_construction(solid))


def mesh(tree, resolution=16):
    return sdf.to_solid(tree, resolution=resolution)


A = sdf.box(10.0)
B = sdf.translate(sdf.box(10.0), (10.0, 0.0, 0.0))    # shares A's +x face
C = sdf.sphere(6.0)


@pytest.fixture(autouse=True)
def _no_forced_engine(monkeypatch):
    """A forced engine would, correctly, bypass the field path."""
    monkeypatch.delenv("YAPCAD_BOOLEAN_ENGINE", raising=False)


@pytest.fixture
def mesh_engine_forbidden(monkeypatch):
    """Fail loudly if anything reaches the native mesh engine."""
    import yapcad.boolean.native as native

    def refuse(*args, **kwargs):
        raise AssertionError("an SDF pair reached the native mesh engine")

    monkeypatch.setattr(native, "solid_boolean", refuse)


class TestRouting:

    @pytest.mark.parametrize("operation,build", [
        ("union", sdf.union),
        ("intersection", sdf.intersect),
        ("difference", sdf.subtract),
    ])
    def test_the_result_tree_is_the_field_operation(self, operation, build,
                                                    mesh_engine_forbidden):
        result = solid_boolean(mesh(A), mesh(C), operation)
        assert tree_of(result) == build(A, C)

    def test_difference_is_first_minus_second(self, mesh_engine_forbidden):
        forwards = tree_of(solid_boolean(mesh(A), mesh(C), "difference"))
        backwards = tree_of(solid_boolean(mesh(C), mesh(A), "difference"))
        assert forwards == sdf.subtract(A, C)
        assert backwards == sdf.subtract(C, A)

    def test_the_result_can_be_combined_again(self, mesh_engine_forbidden):
        first = solid_boolean(mesh(A), mesh(C), "difference")
        second = solid_boolean(first, mesh(sdf.sphere(2.0)), "union")
        assert tree_of(second) == sdf.union(sdf.subtract(A, C),
                                            sdf.sphere(2.0))

    def test_an_explicit_engine_still_wins(self, monkeypatch):
        """Backward compatibility: engine= keeps its meaning."""
        import yapcad.boolean.native as native
        monkeypatch.setattr(native, "solid_boolean",
                            lambda *a, **k: "native-was-used")
        assert solid_boolean(mesh(A), mesh(C), "union",
                             engine="native") == "native-was-used"

    def test_the_environment_override_still_wins(self, monkeypatch):
        import yapcad.boolean.native as native
        monkeypatch.setattr(native, "solid_boolean",
                            lambda *a, **k: "native-was-used")
        monkeypatch.setenv("YAPCAD_BOOLEAN_ENGINE", "native")
        assert solid_boolean(mesh(A), mesh(C), "union") == "native-was-used"

    def test_sdf_can_be_named_as_the_engine(self, mesh_engine_forbidden):
        result = solid_boolean(mesh(A), mesh(C), "union", engine="sdf")
        assert tree_of(result) == sdf.union(A, C)

    def test_naming_sdf_for_a_non_sdf_operand_is_refused(self):
        with pytest.raises(sdf.SdfError, match="second operand"):
            solid_boolean(mesh(A), prism(4.0, 4.0, 4.0), "union",
                          engine="sdf")

    def test_mixed_operands_still_go_to_the_existing_engines(self,
                                                             monkeypatch):
        """One field and one mesh is the rest of Phase 4, not this."""
        import yapcad.boolean.native as native
        monkeypatch.setattr(native, "solid_boolean",
                            lambda *a, **k: "native-was-used")
        cube = prism(4.0, 4.0, 4.0)
        from yapcad.brep import _clear_brep_data
        _clear_brep_data(cube)
        assert solid_boolean(mesh(A), cube, "union") == "native-was-used"

    def test_an_unknown_operation_is_refused(self):
        with pytest.raises(ValueError, match="unsupported boolean"):
            booleans.combine_solids(mesh(A), mesh(C), "xor")


class TestTheCasesMeshBooleansFail:
    """Each of these returned an open mesh, or worse, from the mesh engine."""

    def test_boxes_sharing_a_face(self, mesh_engine_forbidden):
        result = solid_boolean(mesh(A, 20), mesh(B, 20), "union")
        assert issolidclosed(result)
        assert volumeof(result) == pytest.approx(2000.0, rel=1e-5)

    def test_a_coplanar_cut(self, mesh_engine_forbidden):
        slab = sdf.translate(sdf.box((10.0, 4.0, 4.0)), (0.0, 0.0, 3.0))
        result = solid_boolean(mesh(A, 20), mesh(slab, 20), "difference")
        assert issolidclosed(result)
        assert volumeof(result) == pytest.approx(840.0, rel=1e-5)

    def test_a_curved_intersection(self, mesh_engine_forbidden):
        result = solid_boolean(mesh(A, 32), mesh(C, 32), "intersection")
        assert issolidclosed(result)
        assert volumeof(result) == pytest.approx(254.0 * math.pi, rel=5e-3)

    def test_it_is_fast(self, mesh_engine_forbidden):
        """The mesh engine took 7-90 s on these; the field path is meshing
        time only."""
        import time
        a, b = mesh(A, 32), mesh(C, 32)
        started = time.time()
        solid_boolean(a, b, "intersection")
        assert time.time() - started < 5.0


class TestResolution:
    """Never coarser than either operand."""

    def _cell(self, solid):
        m = meshing_of(solid)
        lo, hi = tree_of(solid).bounds
        return max(hi[i] - lo[i] for i in range(3)) / m["resolution"]

    def test_the_finer_operand_sets_the_cell(self):
        coarse, fine = mesh(A, 16), mesh(C, 40)
        result = solid_boolean(coarse, fine, "union")
        assert self._cell(result) <= self._cell(fine) + 1e-12

    def test_a_longer_result_gets_more_cells_not_bigger_ones(self):
        result = solid_boolean(mesh(A, 20), mesh(B, 20), "union")
        assert meshing_of(result)["resolution"] == 40

    def test_the_cap_holds(self, monkeypatch):
        monkeypatch.setattr(booleans, "MAX_RESOLUTION", 24)
        tiny = mesh(sdf.translate(sdf.sphere(0.5), (40.0, 0.0, 0.0)), 64)
        result = solid_boolean(mesh(A, 16), tiny, "union")
        assert meshing_of(result)["resolution"] == 24

    def test_a_transformed_operand_is_estimated_from_its_mesh(self):
        """translatesolid drops the meshing record; the cell size then comes
        from the mesh itself, which for dual contouring is about a cell."""
        moved = translatesolid(mesh(A, 20), point(30.0, 0.0, 0.0))
        assert meshing_of(moved) is None
        result = solid_boolean(moved, mesh(C, 20), "union")
        assert issolidclosed(result)
        # 10/20 = 0.5 per cell on the box; the union is ~41 long.
        assert 60 <= meshing_of(result)["resolution"] <= 110

    def test_the_result_is_deterministic(self):
        first = solid_boolean(mesh(A, 20), mesh(C, 20), "difference")
        second = solid_boolean(mesh(A, 20), mesh(C, 20), "difference")
        assert first[1][0][1] == second[1][0][1]
        assert first[1][0][3] == second[1][0][3]

    def test_a_non_manifold_result_is_retried_once_at_double(self,
                                                             monkeypatch):
        calls = []
        real = booleans.to_solid

        def flaky(node, resolution, **kwargs):
            calls.append(resolution)
            if len(calls) == 1:
                raise booleans.NonManifoldMeshError("thin wall")
            return real(node, resolution=resolution, **kwargs)

        monkeypatch.setattr(booleans, "to_solid", flaky)
        solid_boolean(mesh(A, 16), mesh(C, 16), "union")
        # The first attempt is whatever the resolution rule chose (the box's
        # 0.625 cell over the sphere's 12-unit span: 20); the retry doubles it.
        assert len(calls) == 2
        assert calls[1] == 2 * calls[0]

    def test_at_the_cap_the_refusal_stands(self, monkeypatch):
        monkeypatch.setattr(booleans, "MAX_RESOLUTION", 16)

        def always(*args, **kwargs):
            raise booleans.NonManifoldMeshError("thin wall")

        monkeypatch.setattr(booleans, "to_solid", always)
        with pytest.raises(booleans.NonManifoldMeshError):
            solid_boolean(mesh(A, 16), mesh(C, 16), "union")


class TestEdgeCases:

    def test_a_disjoint_intersection_is_an_empty_solid_with_its_tree(self):
        far = mesh(sdf.translate(sdf.box(4.0), (50.0, 0.0, 0.0)))
        result = solid_boolean(mesh(A), far, "intersection")
        assert issolid(result, fast=False)
        assert result[1] == []
        assert tree_of(result).kind == "intersect"

    def test_the_result_round_trips_through_geometry_json(self):
        from yapcad.io.geometry_json import (
            geometry_from_json,
            geometry_to_json,
        )
        result = solid_boolean(mesh(A), mesh(C), "difference")
        document = json.loads(json.dumps(geometry_to_json([result],
                                                          units="mm")))
        assert tree_of(geometry_from_json(document)[0]) == tree_of(result)

    def test_placement_carries_through(self):
        moved = translatesolid(mesh(A), point(100.0, 0.0, 0.0))
        result = solid_boolean(moved, mesh(sdf.translate(C, (100, 0, 0))),
                               "intersection")
        lo, hi = solidbbox(result)[:2]
        assert lo[0] > 90.0 and hi[0] < 110.0


class TestNary:

    def test_combine_all_builds_one_nary_node(self):
        parts = [mesh(sdf.translate(sdf.box(4.0), (6.0 * i, 0.0, 0.0)))
                 for i in range(4)]
        result = booleans.combine_all(parts, "union")
        tree = tree_of(result)
        assert tree.kind == "union" and len(tree.children) == 4
        assert volumeof(result) == pytest.approx(256.0, rel=1e-5)

    def test_difference_subtracts_every_tool_from_the_first(self):
        tools = [mesh(sdf.sphere(3.0)),
                 mesh(sdf.translate(sdf.box(2.0), (4.0, 4.0, 4.0)))]
        result = booleans.combine_all([mesh(A)] + tools, "difference")
        expected = sdf.subtract(A, *[tree_of(t) for t in tools])
        assert tree_of(result) == expected

    def test_combine_all_needs_two_operands(self):
        with pytest.raises(sdf.SdfError, match="at least two"):
            booleans.combine_all([mesh(A)], "union")

    def test_a_later_operand_is_named_when_it_is_not_sdf(self):
        with pytest.raises(sdf.SdfError, match="operand 3"):
            booleans.combine_all([mesh(A), mesh(C), prism(2, 2, 2)], "union")


class TestDslBooleans:
    """The DSL's union/difference/intersection mesh once for SDF operands."""

    def _call(self, name, solids):
        from yapcad.dsl.runtime.builtins import call_builtin
        from yapcad.dsl.runtime.values import solid_val
        return call_builtin(name, [solid_val(s) for s in solids]).data

    def test_union_of_many_is_one_node(self, mesh_engine_forbidden):
        parts = [mesh(sdf.translate(sdf.box(4.0), (6.0 * i, 0.0, 0.0)))
                 for i in range(4)]
        tree = tree_of(self._call("union", parts))
        assert tree.kind == "union" and len(tree.children) == 4

    def test_it_meshes_once(self, monkeypatch):
        calls = []
        real = booleans.to_solid

        def counting(*args, **kwargs):
            calls.append(1)
            return real(*args, **kwargs)

        parts = [mesh(sdf.translate(sdf.box(4.0), (6.0 * i, 0.0, 0.0)))
                 for i in range(5)]
        monkeypatch.setattr(booleans, "to_solid", counting)
        self._call("union", parts)
        assert len(calls) == 1

    def test_difference_with_several_tools(self, mesh_engine_forbidden):
        tools = [mesh(sdf.sphere(3.0)),
                 mesh(sdf.translate(sdf.box(2.0), (4.0, 4.0, 4.0)))]
        result = self._call("difference", [mesh(A)] + tools)
        assert tree_of(result).kind == "subtract"
        assert len(tree_of(result).children) == 3

    def test_intersection(self, mesh_engine_forbidden):
        result = self._call("intersection", [mesh(A, 32), mesh(C, 32)])
        assert volumeof(result) == pytest.approx(254.0 * math.pi, rel=5e-3)

    def test_mixed_operands_keep_the_pairwise_path(self, monkeypatch):
        import yapcad.boolean.native as native
        seen = []
        monkeypatch.setattr(native, "solid_boolean",
                            lambda a, b, op, **k: seen.append(op) or a)
        cube = prism(4.0, 4.0, 4.0)
        from yapcad.brep import _clear_brep_data
        _clear_brep_data(cube)
        self._call("union", [mesh(A), cube])
        assert seen == ["union"]


@needs_occ
class TestDerivedBrep:

    def test_exact_operands_give_an_exact_result(self):
        from yapcad.brep import brep_from_solid
        from yapcad.metadata import get_solid_metadata
        from yapcad.sdf import occ
        a = sdf.to_solid(A, resolution=16, brep=True)
        c = sdf.to_solid(C, resolution=16, brep=True)
        result = solid_boolean(a, c, "intersection")
        marker = get_solid_metadata(result)["brep"]["derivedFrom"]["sdf"]
        assert marker == tree_of(result).digest
        assert occ.volume(brep_from_solid(result).shape) == pytest.approx(
            254.0 * math.pi, rel=1e-6)

    def test_one_inexact_operand_means_no_brep(self):
        from yapcad.brep import has_brep_data
        a = sdf.to_solid(A, resolution=16, brep=True)
        c = sdf.to_solid(C, resolution=16)
        assert not has_brep_data(solid_boolean(a, c, "union"))

    def test_a_non_replayable_result_quietly_has_none(self):
        from yapcad.brep import has_brep_data
        a = sdf.to_solid(A, resolution=16, brep=True)
        blob = sdf.to_solid(sdf.smooth_union(sdf.sphere(3.0), sdf.box(4.0),
                                             1.0), resolution=16)
        assert not has_brep_data(solid_boolean(a, blob, "union"))
