"""Tests for dual contouring and SDF-authoritative solids (Phase 2).

The design document's claims for this phase are specific, so the tests are
too: dual contouring rather than marching cubes *because* it keeps sharp
features; a closed manifold so the rest of yapCAD can consume the result;
and bit-reproducible output, because yapCAD signs packages and a preview
that regenerated differently run to run would churn package hashes.
"""

import json
import math
import subprocess
import sys

import numpy as np
import pytest

from yapcad import sdf
from yapcad.construction import get_construction
from yapcad.geom3d import issolid, issolidclosed, solid_shells, volumeof
from yapcad.io.geometry_json import (
    SCHEMA_ID,
    geometry_from_json,
    geometry_to_json,
)


# ---------------------------------------------------------------------------
# Mesh-level helpers
# ---------------------------------------------------------------------------


def edge_multiplicity(triangles):
    """Unique undirected edges and how many triangles use each."""
    edges = np.concatenate([
        triangles[:, [0, 1]], triangles[:, [1, 2]], triangles[:, [2, 0]]
    ])
    return np.unique(np.sort(edges, axis=1), axis=0,
                     return_counts=True)


def mesh_volume(mesh):
    """Signed volume by the divergence theorem; positive when wound out."""
    a = mesh.vertices[mesh.triangles[:, 0]]
    b = mesh.vertices[mesh.triangles[:, 1]]
    c = mesh.vertices[mesh.triangles[:, 2]]
    return float(np.sum(np.einsum("ij,ij->i", a, np.cross(b, c))) / 6.0)


def euler_characteristic(mesh):
    edges, _ = edge_multiplicity(mesh.triangles)
    return len(mesh.vertices) - len(edges) + len(mesh.triangles)


def assert_closed_manifold(mesh):
    """Every edge shared by exactly two triangles, and none degenerate."""
    _, counts = edge_multiplicity(mesh.triangles)
    assert np.all(counts == 2), (
        f"{int(np.sum(counts != 2))} non-manifold edges"
    )
    tris = mesh.triangles
    assert np.all(tris[:, 0] != tris[:, 1])
    assert np.all(tris[:, 1] != tris[:, 2])
    assert np.all(tris[:, 0] != tris[:, 2])


# ---------------------------------------------------------------------------
# Topology
# ---------------------------------------------------------------------------


class TestTopology:

    @pytest.mark.parametrize("name,node", [
        ("sphere", sdf.sphere(5.0)),
        ("box", sdf.box(8.0)),
        ("rounded_box", sdf.rounded_box((10.0, 6.0, 4.0), 1.0)),
        ("cylinder", sdf.cylinder(4.0, 10.0)),
        ("cone", sdf.cone(5.0, 1.0, 9.0)),
        ("capsule", sdf.capsule((-4, 0, 0), (4, 0, 0), 2.0)),
        ("union", sdf.union(sdf.sphere(4.0), sdf.box(6.0))),
        ("subtract", sdf.subtract(sdf.box(10.0), sdf.cylinder(3.0, 20.0))),
        ("smooth_union", sdf.smooth_union(
            sdf.translate(sdf.sphere(4.0), (-3, 0, 0)),
            sdf.translate(sdf.sphere(4.0), (3, 0, 0)), 2.0)),
        ("shell", sdf.shell(sdf.sphere(6.0), 1.0)),
    ])
    def test_meshes_are_closed_manifolds(self, name, node):
        assert_closed_manifold(sdf.dual_contour(node, resolution=24))

    def test_simply_connected_shapes_have_euler_two(self):
        for node in (sdf.sphere(5.0), sdf.box(8.0), sdf.cylinder(4.0, 10.0)):
            assert euler_characteristic(
                sdf.dual_contour(node, resolution=24)
            ) == 2

    def test_torus_comes_out_genus_one(self):
        mesh = sdf.dual_contour(sdf.torus(9.0, 2.5), resolution=32)
        assert euler_characteristic(mesh) == 0

    def test_a_cube_frame_comes_out_genus_five(self):
        """Box minus a sphere large enough to breach all six faces leaves the
        twelve edges of a frame: V - E + 1 = 5 handles, so chi = -8.  A good
        check that the contouring is not quietly closing off holes."""
        node = sdf.subtract(sdf.box(10.0), sdf.sphere(6.0))
        mesh = sdf.dual_contour(node, resolution=32)
        assert_closed_manifold(mesh)
        assert euler_characteristic(mesh) == -8

    def test_a_hollow_shell_has_two_components(self):
        """A shelled sphere is two nested spheres: chi = 4."""
        node = sdf.shell(sdf.sphere(6.0), 1.0)
        assert euler_characteristic(sdf.dual_contour(node, resolution=32)) == 4

    def test_winding_is_outward(self):
        for node in (sdf.sphere(5.0), sdf.box(8.0), sdf.torus(9.0, 2.5)):
            assert mesh_volume(sdf.dual_contour(node, resolution=24)) > 0.0

    def test_normals_agree_with_the_field_gradient(self):
        mesh = sdf.dual_contour(sdf.sphere(5.0), resolution=24)
        radial = mesh.vertices / np.linalg.norm(mesh.vertices, axis=1)[:, None]
        assert np.allclose(np.sum(mesh.normals * radial, axis=1), 1.0,
                           atol=1e-4)

    def test_an_empty_region_yields_an_empty_mesh(self):
        mesh = sdf.dual_contour(
            sdf.sphere(1.0),
            bounds=((50.0, 50.0, 50.0), (60.0, 60.0, 60.0)),
            resolution=8,
        )
        assert mesh.vertices.shape == (0, 3)
        assert mesh.triangles.shape == (0, 3)


# ---------------------------------------------------------------------------
# Accuracy, and the sharp-feature claim
# ---------------------------------------------------------------------------


class TestAccuracy:

    def test_sharp_edges_are_preserved_exactly(self):
        """The reason §5.1 chose dual contouring over marching cubes.  A box
        aligned to the lattice recovers its volume to floating-point, because
        the QEF puts a vertex exactly on each edge and corner; marching cubes
        would chamfer them and lose volume.

        The residual is around 1e-8 relative, set by the finite-difference
        gradient that supplies the QEF its plane normals -- not by the
        reconstruction, which is exact for a lattice-aligned box."""
        mesh = sdf.dual_contour(sdf.box(8.0), resolution=24)
        assert mesh_volume(mesh) == pytest.approx(512.0, rel=1e-6)

    def test_sharp_edges_survive_an_off_axis_box(self):
        node = sdf.rotate(sdf.box((8.0, 6.0, 4.0)), (0, 0, 1), 30.0)
        mesh = sdf.dual_contour(node, resolution=48)
        assert mesh_volume(mesh) == pytest.approx(192.0, rel=2e-3)

    @pytest.mark.parametrize("node,exact", [
        (sdf.sphere(5.0), 4.0 / 3.0 * math.pi * 125.0),
        (sdf.cylinder(4.0, 10.0), math.pi * 16.0 * 10.0),
        (sdf.torus(9.0, 2.5), 2.0 * math.pi ** 2 * 9.0 * 2.5 ** 2),
    ])
    def test_curved_volumes_converge(self, node, exact):
        assert mesh_volume(
            sdf.dual_contour(node, resolution=48)
        ) == pytest.approx(exact, rel=5e-3)

    def test_accuracy_improves_with_resolution(self):
        exact = 4.0 / 3.0 * math.pi * 125.0
        coarse = abs(
            mesh_volume(sdf.dual_contour(sdf.sphere(5.0), 16)) - exact)
        fine = abs(mesh_volume(sdf.dual_contour(sdf.sphere(5.0), 64)) - exact)
        assert fine < coarse

    def test_vertices_lie_on_the_surface(self):
        node = sdf.subtract(sdf.box(10.0), sdf.sphere(4.0))
        mesh = sdf.dual_contour(node, resolution=32)
        residual = np.abs(sdf.evaluate(node, mesh.vertices))
        spacing = 10.0 / 32.0 * 4.0
        assert np.max(residual) < spacing

    def test_vertices_stay_inside_their_own_cells(self):
        """The QEF is clamped; an unclamped solve can fling a vertex far out
        of its cell on a near-degenerate configuration and tangle the mesh."""
        node = sdf.subtract(sdf.box(9.0), sdf.sphere(5.5))
        mesh = sdf.dual_contour(node, resolution=24)
        lo, hi = node.bounds
        pad = 9.0 / 24.0 * 4.0
        assert np.all(mesh.vertices >= np.array(lo) - pad)
        assert np.all(mesh.vertices <= np.array(hi) + pad)


# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------


class TestDeterminism:
    """Package signing depends on every one of these."""

    def test_repeated_runs_are_bit_identical(self):
        node = sdf.subtract(sdf.rounded_box((20.0, 20.0, 8.0), 2.0),
                            sdf.cylinder(4.0, 20.0))
        first = sdf.dual_contour(node, resolution=24)
        second = sdf.dual_contour(node, resolution=24)
        assert np.array_equal(first.vertices, second.vertices)
        assert np.array_equal(first.normals, second.normals)
        assert np.array_equal(first.triangles, second.triangles)

    def test_output_is_identical_across_processes(self):
        script = (
            "import hashlib, numpy as np;"
            "from yapcad import sdf;"
            "m = sdf.dual_contour(sdf.subtract(sdf.box(10.0), "
            "sdf.sphere(4.0)), resolution=20);"
            "print(hashlib.sha256(m.vertices.tobytes() + "
            "m.triangles.tobytes()).hexdigest())"
        )
        runs = {
            subprocess.run([sys.executable, "-c", script],
                           capture_output=True, text=True,
                           check=True).stdout.strip()
            for _ in range(2)
        }
        assert len(runs) == 1

    def test_lipschitz_pruning_does_not_change_the_result(self):
        """Pruning skips chunks the Lipschitz bound proves the surface cannot
        reach.  If the margin were wrong this would differ, so it is the test
        that makes the optimisation safe to trust."""
        for node in (
            sdf.sphere(5.0),
            sdf.subtract(sdf.box(10.0), sdf.sphere(6.0)),
            sdf.torus(9.0, 2.5),
            sdf.smooth_union(sdf.translate(sdf.sphere(4.0), (-3, 0, 0)),
                             sdf.translate(sdf.sphere(4.0), (3, 0, 0)), 2.0),
        ):
            pruned = sdf.dual_contour(node, resolution=24, prune=True)
            plain = sdf.dual_contour(node, resolution=24, prune=False)
            assert np.array_equal(pruned.vertices, plain.vertices)
            assert np.array_equal(pruned.triangles, plain.triangles)

    def test_chunk_size_does_not_change_the_result(self):
        node = sdf.subtract(sdf.box(10.0), sdf.sphere(6.0))
        base = sdf.dual_contour(node, resolution=24, chunk=64)
        for chunk in (4, 7, 16):
            other = sdf.dual_contour(node, resolution=24, chunk=chunk)
            assert np.array_equal(base.vertices, other.vertices)
            assert np.array_equal(base.triangles, other.triangles)


# ---------------------------------------------------------------------------
# Region selection
# ---------------------------------------------------------------------------


class TestRegion:

    def test_padding_keeps_the_surface_off_the_lattice_boundary(self):
        """Without padding a tangent surface would be clipped by the lattice
        edge and the mesh left open."""
        assert_closed_manifold(sdf.dual_contour(sdf.box(8.0), resolution=16))

    def test_explicit_bounds_are_honoured(self):
        node = sdf.sphere(5.0)
        mesh = sdf.dual_contour(
            node, resolution=16,
            bounds=((-20.0, -20.0, -20.0), (20.0, 20.0, 20.0)),
        )
        assert_closed_manifold(mesh)
        assert mesh_volume(mesh) == pytest.approx(
            4.0 / 3.0 * math.pi * 125.0, rel=0.05
        )

    def test_an_unbounded_field_is_refused_with_a_usable_message(self):
        with pytest.raises(sdf.SdfError, match="unbounded"):
            sdf.dual_contour(sdf.gyroid(5.0, 0.3))

    def test_a_lattice_bounded_by_a_solid_meshes(self):
        """The documented way to get a finite lattice: intersect the infinite
        field with a solid, which gives the tree finite bounds."""
        node = sdf.intersect(sdf.box(20.0), sdf.gyroid(8.0, 0.5))
        mesh = sdf.dual_contour(node, resolution=48)
        assert len(mesh.vertices) > 0
        assert_closed_manifold(mesh)


class TestManifoldGuard:
    """One vertex per cell is dual contouring's structural limit: a cell the
    surface crosses twice cannot represent both sheets.  §5.3 says state the
    limit rather than degrade silently, so it is detected and refused."""

    UNDER_RESOLVED = sdf.intersect(sdf.box(20.0), sdf.gyroid(8.0, 0.5))

    def test_an_under_resolved_wall_is_detected(self):
        mesh = sdf.dual_contour(self.UNDER_RESOLVED, resolution=32)
        assert not sdf.is_manifold(mesh)
        assert sdf.manifold_defects(mesh)["nonmanifold"] > 0

    def test_to_solid_refuses_with_an_actionable_message(self):
        with pytest.raises(sdf.SdfError) as excinfo:
            sdf.to_solid(self.UNDER_RESOLVED, resolution=32)
        message = str(excinfo.value)
        assert "non-manifold" in message
        assert "Raise the resolution" in message

    def test_more_resolution_resolves_it(self):
        result = sdf.to_solid(self.UNDER_RESOLVED, resolution=48)
        assert issolid(result, fast=False)

    def test_the_check_can_be_waived(self):
        result = sdf.to_solid(self.UNDER_RESOLVED, resolution=32, check=False)
        assert issolid(result, fast=False)

    def test_sound_meshes_report_no_defects(self):
        for node in (sdf.sphere(5.0), sdf.box(8.0), sdf.torus(9.0, 2.5)):
            mesh = sdf.dual_contour(node, resolution=24)
            assert sdf.manifold_defects(mesh) == {
                "boundary": 0, "nonmanifold": 0, "coincident": 0
            }

    def test_an_empty_mesh_is_not_a_defect(self):
        mesh = sdf.dual_contour(
            sdf.sphere(1.0), resolution=8,
            bounds=((50.0, 50.0, 50.0), (60.0, 60.0, 60.0)),
        )
        assert sdf.manifold_defects(mesh) == {
            "boundary": 0, "nonmanifold": 0, "coincident": 0
        }

    def test_resolution_must_be_at_least_two(self):
        with pytest.raises(sdf.SdfError, match="resolution"):
            sdf.dual_contour(sdf.sphere(1.0), resolution=1)

    def test_negative_padding_is_refused(self):
        with pytest.raises(sdf.SdfError, match="padding"):
            sdf.dual_contour(sdf.sphere(1.0), padding=-1.0)


# ---------------------------------------------------------------------------
# Conversion to a yapCAD solid
# ---------------------------------------------------------------------------


class TestToSolid:

    def test_result_is_a_well_formed_closed_solid(self):
        result = sdf.to_solid(sdf.sphere(5.0), resolution=24)
        assert issolid(result, fast=False)
        assert issolidclosed(result)

    def test_volume_matches_the_field(self):
        result = sdf.to_solid(sdf.sphere(5.0), resolution=32)
        assert volumeof(result) == pytest.approx(
            4.0 / 3.0 * math.pi * 125.0, rel=1e-2
        )

    def test_volume_is_positive_so_the_winding_reached_yapcad_intact(self):
        for node in (sdf.box(8.0), sdf.torus(9.0, 2.5)):
            assert volumeof(sdf.to_solid(node, resolution=24)) > 0.0

    def test_it_reports_one_outer_shell_and_no_voids(self):
        outer, void = solid_shells(sdf.to_solid(sdf.sphere(5.0), 24))
        assert len(outer) == 1
        assert len(void) == 0

    def test_surface_uses_yapcad_point_and_vector_conventions(self):
        surf = sdf.to_solid(sdf.sphere(5.0), resolution=12)[1][0]
        assert all(len(v) == 4 and v[3] == 1.0 for v in surf[1][:20])
        assert all(len(n) == 4 and n[3] == 0.0 for n in surf[2][:20])
        assert len(surf[1]) == len(surf[2])
        assert surf[4] == [] and surf[5] == []

    def test_the_generating_field_rides_along(self):
        node = sdf.subtract(sdf.box(10.0), sdf.sphere(4.0))
        record = get_construction(sdf.to_solid(node, resolution=12))
        assert sdf.from_construction(record) == node

    def test_meshing_parameters_are_recorded_for_reproducibility(self):
        record = get_construction(sdf.to_solid(sdf.sphere(5.0), resolution=17))
        assert sdf.meshing_from_construction(record) == {
            "method": "dual-contouring",
            "resolution": 17,
        }

    def test_recorded_parameters_regenerate_the_same_mesh(self):
        node = sdf.subtract(sdf.box(10.0), sdf.sphere(4.0))
        first = sdf.to_solid(node, resolution=20)
        params = sdf.meshing_from_construction(get_construction(first))
        second = sdf.to_solid(node, resolution=params["resolution"])
        assert first[1][0][1] == second[1][0][1]
        assert first[1][0][3] == second[1][0][3]

    def test_optional_parameters_are_recorded_only_when_given(self):
        bounds = ((-9.0, -9.0, -9.0), (9.0, 9.0, 9.0))
        record = get_construction(
            sdf.to_solid(sdf.sphere(5.0), 12, bounds=bounds, padding=1.5)
        )
        meshing = sdf.meshing_from_construction(record)
        assert meshing["padding"] == 1.5
        assert meshing["bounds"] == [-9.0, -9.0, -9.0, 9.0, 9.0, 9.0]

    def test_metadata_is_merged(self):
        from yapcad.metadata import get_solid_metadata
        result = sdf.to_solid(sdf.sphere(5.0), 12, metadata={"name": "ball"})
        assert get_solid_metadata(result)["name"] == "ball"


# ---------------------------------------------------------------------------
# Geometry JSON v0.3
# ---------------------------------------------------------------------------


def roundtrip(entities):
    document = json.loads(json.dumps(geometry_to_json(entities, units="mm")))
    return document, geometry_from_json(document)


def solid_entry(document):
    return next(e for e in document["entities"] if e["type"] == "solid")


class TestGeometryJson:

    def test_documents_declare_v03(self):
        document, _ = roundtrip([sdf.to_solid(sdf.sphere(5.0), 12)])
        assert document["schema"] == "yapcad-geometry-json-v0.3"
        assert SCHEMA_ID == "yapcad-geometry-json-v0.3"

    def test_an_sdf_solid_is_sdf_authoritative(self):
        document, _ = roundtrip([sdf.to_solid(sdf.sphere(5.0), 12)])
        reps = solid_entry(document)["representations"]
        assert reps["authoritative"] == "sdf"
        assert reps["mesh"]["role"] == "preview"
        assert reps["mesh"]["generatedFrom"] == "sdf"
        assert reps["mesh"]["method"] == "dual-contouring"
        assert reps["mesh"]["resolution"] == 12

    def test_the_representation_records_the_field_analysis(self):
        node = sdf.smooth_union(sdf.sphere(4.0), sdf.box(6.0), 1.0)
        document, _ = roundtrip([sdf.to_solid(node, 12)])
        record = solid_entry(document)["representations"]["sdf"]
        assert record["role"] == "authoritative"
        assert record["format"] == "yapcad-sdf-tree-v1"
        assert record["lipschitz"] == "bound"
        assert record["lipschitzConstant"] == pytest.approx(
            node.lipschitz)
        assert record["csgExact"] is False
        assert record["bounds"] == [float(v) for v in
                                    (*node.bounds[0], *node.bounds[1])]

    def test_csg_exactness_is_surfaced_for_step_consumers(self):
        """§5.2: a consumer must be able to tell before attempting STEP."""
        csg = sdf.subtract(sdf.box(10.0), sdf.cylinder(3.0, 20.0))
        document, _ = roundtrip([sdf.to_solid(csg, 12)])
        assert solid_entry(document)["representations"]["sdf"]["csgExact"]

        # check=False: this asserts on the recorded classification, and a
        # coarse lattice mesh is beside the point.
        field = sdf.intersect(sdf.box(10.0), sdf.gyroid(4.0, 0.4))
        document, _ = roundtrip([sdf.to_solid(field, 12, check=False)])
        assert not solid_entry(document)["representations"]["sdf"]["csgExact"]

    def test_exact_primitives_report_exact_lipschitz(self):
        document, _ = roundtrip([sdf.to_solid(sdf.sphere(5.0), 12)])
        record = solid_entry(document)["representations"]["sdf"]
        assert record["lipschitz"] == "exact"
        assert record["lipschitzConstant"] == 1.0

    def test_the_tree_is_not_duplicated_in_the_construction_field(self):
        """One encoding per fact: §8.3's lesson from the voids defect."""
        document, _ = roundtrip([sdf.to_solid(sdf.sphere(5.0), 12)])
        assert "construction" not in solid_entry(document)

    def test_the_field_survives_the_round_trip(self):
        node = sdf.subtract(sdf.rounded_box((20.0, 20.0, 8.0), 2.0),
                            sdf.cylinder(4.0, 20.0))
        _, restored = roundtrip([sdf.to_solid(node, 12)])
        recovered = sdf.from_construction(get_construction(restored[0]))
        assert recovered == node
        assert recovered.digest == node.digest

    def test_meshing_parameters_survive_the_round_trip(self):
        original = sdf.to_solid(sdf.sphere(5.0), resolution=19, padding=0.5)
        _, restored = roundtrip([original])
        assert (sdf.meshing_from_construction(get_construction(restored[0]))
                == sdf.meshing_from_construction(get_construction(original)))

    def test_documents_validate_against_the_published_schema(self):
        jsonschema = pytest.importorskip("jsonschema")
        from pathlib import Path
        path = (Path(__file__).resolve().parents[1] / "docs" / "schemas"
                / "yapcad-geometry-json-v0.3.schema.json")
        schema = json.loads(path.read_text(encoding="utf-8"))
        node = sdf.subtract(sdf.box(10.0), sdf.sphere(4.0))
        document, _ = roundtrip([sdf.to_solid(node, 12)])
        jsonschema.Draft202012Validator(schema).validate(document)

    def test_a_tree_serialises_to_a_fraction_of_a_brep_payload(self):
        """§4.1's size argument, measured rather than asserted."""
        node = sdf.subtract(sdf.rounded_box((40.0, 40.0, 10.0), 3.0),
                            sdf.cylinder(6.0, 20.0))
        document, _ = roundtrip([sdf.to_solid(node, 12)])
        record = solid_entry(document)["representations"]["sdf"]
        assert len(json.dumps(record["tree"])) < 2000


class TestGeometryJsonRejections:

    def test_a_doctored_csg_exact_claim_is_rejected(self):
        """The analysis fields cache what the tree already determines, so
        they are verified on read rather than believed -- the same posture as
        the BREP payload hash check."""
        document, _ = roundtrip([
            sdf.to_solid(sdf.intersect(sdf.box(10.0), sdf.gyroid(4.0, 0.4)),
                         12, check=False)
        ])
        solid_entry(document)["representations"]["sdf"]["csgExact"] = True
        with pytest.raises(ValueError, match="csgExact disagrees"):
            geometry_from_json(document)

    def test_a_doctored_lipschitz_constant_is_rejected(self):
        document, _ = roundtrip([sdf.to_solid(sdf.sphere(5.0), 12)])
        reps = solid_entry(document)["representations"]
        reps["sdf"]["lipschitzConstant"] = 4.0
        with pytest.raises(ValueError, match="lipschitzConstant disagrees"):
            geometry_from_json(document)

    def test_a_doctored_lipschitz_kind_is_rejected(self):
        document, _ = roundtrip([sdf.to_solid(sdf.sphere(5.0), 12)])
        solid_entry(document)["representations"]["sdf"]["lipschitz"] = "bound"
        with pytest.raises(ValueError, match="lipschitz kind disagrees"):
            geometry_from_json(document)

    def test_a_corrupt_tree_is_rejected(self):
        document, _ = roundtrip([sdf.to_solid(sdf.sphere(5.0), 12)])
        tree = solid_entry(document)["representations"]["sdf"]["tree"]
        tree["nodes"][tree["root"]]["kind"] = "hyperboloid"
        with pytest.raises(ValueError, match="unknown SDF node kind"):
            geometry_from_json(document)

    def test_a_missing_sdf_record_is_rejected(self):
        document, _ = roundtrip([sdf.to_solid(sdf.sphere(5.0), 12)])
        del solid_entry(document)["representations"]["sdf"]
        with pytest.raises(ValueError, match="missing SDF representation"):
            geometry_from_json(document)

    def test_a_non_authoritative_sdf_record_is_rejected(self):
        document, _ = roundtrip([sdf.to_solid(sdf.sphere(5.0), 12)])
        solid_entry(document)["representations"]["authoritative"] = "mesh"
        solid_entry(document)["representations"]["mesh"]["role"] = (
            "authoritative")
        with pytest.raises(ValueError, match="non-authoritative"):
            geometry_from_json(document)

    def test_a_duplicated_tree_is_rejected(self):
        document, _ = roundtrip([sdf.to_solid(sdf.sphere(5.0), 12)])
        entry = solid_entry(document)
        tree = entry["representations"]["sdf"]["tree"]
        entry["construction"] = ["sdf", tree]
        with pytest.raises(ValueError, match="duplicates its SDF tree"):
            geometry_from_json(document)

    def test_the_mesh_role_must_be_preview(self):
        document, _ = roundtrip([sdf.to_solid(sdf.sphere(5.0), 12)])
        solid_entry(document)["representations"]["mesh"]["role"] = (
            "authoritative")
        with pytest.raises(ValueError, match="mesh role must be 'preview'"):
            geometry_from_json(document)


class TestDownstreamUsability:
    """yapCAD's own predicates are the bar, not our internal edge count.

    ``issolidclosed`` keys edges by vertex *position*, so a mesh that is
    perfectly sound by index can still read as open if two vertices share a
    point -- and then ``volumeof`` refuses it.  These caught exactly that:
    a torus meshed above resolution 24 used to fail here.
    """

    CASES = [
        ("sphere", sdf.sphere(5.0), 32),
        ("box", sdf.box(8.0), 32),
        ("torus", sdf.torus(9.0, 2.5), 32),
        ("torus-fine", sdf.torus(9.0, 2.5), 64),
        ("cone", sdf.cone(5.0, 1.0, 9.0), 32),
        ("capsule", sdf.capsule((-4, 0, 0), (4, 0, 0), 2.0), 32),
        ("cube-frame", sdf.subtract(sdf.box(10.0), sdf.sphere(6.0)), 32),
        ("rotated-box",
         sdf.rotate(sdf.box((8.0, 6.0, 4.0)), (0, 0, 1), 30), 48),
        ("smooth-union", sdf.smooth_union(
            sdf.translate(sdf.sphere(4.0), (-3, 0, 0)),
            sdf.translate(sdf.sphere(4.0), (3, 0, 0)), 2.0), 32),
        ("shell", sdf.shell(sdf.sphere(6.0), 1.0), 40),
        ("lattice", sdf.intersect(sdf.box(20.0), sdf.gyroid(8.0, 0.5)), 48),
    ]

    @pytest.mark.parametrize("name,node,resolution", CASES)
    def test_solids_satisfy_yapcads_own_closure_test(self, name, node,
                                                     resolution):
        result = sdf.to_solid(node, resolution=resolution)
        assert issolidclosed(result), name
        assert volumeof(result) > 0.0

    @pytest.mark.parametrize("name,node,resolution", CASES)
    def test_no_two_vertices_share_a_position(self, name, node, resolution):
        mesh = sdf.dual_contour(node, resolution=resolution)
        assert sdf.manifold_defects(mesh)["coincident"] == 0, name

    def test_meshes_export_as_watertight_stl(self):
        """The Phase 2 promise, end to end: an SDF part goes out through
        yapCAD's existing STL writer with no knowledge of fields."""
        trimesh = pytest.importorskip("trimesh")
        import os
        import tempfile
        from yapcad.io.stl import write_stl

        node = sdf.subtract(sdf.box(10.0), sdf.sphere(6.0))
        result = sdf.to_solid(node, resolution=32)
        handle, path = tempfile.mkstemp(suffix=".stl")
        os.close(handle)
        try:
            write_stl(result, path)
            loaded = trimesh.load_mesh(path)
            assert loaded.is_watertight
            assert loaded.is_winding_consistent
            assert loaded.euler_number == -8
        finally:
            os.unlink(path)

    def test_the_qef_keeps_vertices_in_their_cells(self):
        """Singular-value truncation, rather than uniform damping, is what
        stops a flat cell sliding its vertex along the surface and out of
        the cell -- which was the source of the coincident vertices."""
        node = sdf.torus(9.0, 2.5)
        mesh = sdf.dual_contour(node, resolution=32)
        assert sdf.manifold_defects(mesh)["coincident"] == 0


# ---------------------------------------------------------------------------
# Manifold dual contouring
# ---------------------------------------------------------------------------


def _wedge(degrees):
    """A thin wedge, tilted off the lattice: near its edge the material is
    thinner than any cell, which one vertex per cell cannot represent."""
    import math
    a = math.radians(degrees)
    wedge = sdf.intersect(
        sdf.half_space((math.sin(a / 2), math.cos(a / 2), 0), 0),
        sdf.half_space((math.sin(a / 2), -math.cos(a / 2), 0), 0),
        sdf.translate(sdf.box((10, 10, 4)), (-5.3, 0.17, 0.05)))
    return sdf.rotate(wedge, (0.3, 0.2, 1.0), 17.0)


@pytest.mark.parametrize("resolution", [24, 32])
def test_a_thin_wedge_meshes_manifold(resolution):
    """One vertex per cell gave non-manifold edges here at every
    resolution; one vertex per surface component of a cell does not."""
    from yapcad.sdf.contour import manifold_defects
    from yapcad.sdf.contour import dual_contour
    mesh = dual_contour(_wedge(20.0), resolution=resolution)
    assert not any(manifold_defects(mesh).values())


def test_component_table_is_consistent():
    from yapcad.sdf.contour import _COMPONENT, _COMPONENT_COUNT
    assert _COMPONENT_COUNT[0] == 0 and _COMPONENT_COUNT[255] == 0
    # A single inside corner is one component; two body-diagonal corners
    # (no shared face) are two.
    assert _COMPONENT_COUNT[1] == 1
    assert _COMPONENT_COUNT[(1 << 0) | (1 << 7)] == 2
    # Complementary configurations have the same crossing edges.
    for config in range(256):
        assert np.array_equal(_COMPONENT[config] >= 0,
                              _COMPONENT[255 - config] >= 0)
