"""Tests for the direct manifold3d boolean engine.

Volumes are checked against closed-form values, not against manifold3d run
some other way, so agreement means something.  Included are the cases the
native engine gets wrong (see test_boolean_native.py): general-position
overlaps, shared faces, and coincident side faces.
"""

import warnings

import pytest

from yapcad.boolean import manifold_engine
from yapcad.brep import _clear_brep_data
from yapcad.construction import get_construction
from yapcad.geom import point
from yapcad.geom3d import issolidclosed, solid, solid_boolean, translatesolid
from yapcad.geom3d_util import conic, prism, sphere

pytestmark = pytest.mark.skipif(
    not manifold_engine.is_available(),
    reason="manifold3d with Mesh64 is not installed")


def mesh_only(s):
    _clear_brep_data(s)
    return s


def moved(s, delta):
    return mesh_only(translatesolid(s, point(*delta)))


def signed_volume(s):
    total = 0.0
    for surf in s[1]:
        v = surf[1]
        for i, j, k in surf[3]:
            a, b, c = v[i], v[j], v[k]
            total += (a[0] * (b[1] * c[2] - b[2] * c[1])
                      - a[1] * (b[0] * c[2] - b[2] * c[0])
                      + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6.0
    return total


def run(a, b, op):
    return solid_boolean(a, b, op, engine="manifold")


def assert_correct(result, volume):
    assert issolidclosed(result)
    assert signed_volume(result) == pytest.approx(volume, rel=1e-9, abs=1e-12)


@pytest.fixture(autouse=True)
def _no_forced_engines(monkeypatch):
    monkeypatch.delenv("YAPCAD_BOOLEAN_ENGINE", raising=False)
    monkeypatch.delenv("YAPCAD_MESH_BOOLEAN_ENGINE", raising=False)


class TestCorrectness:

    @pytest.mark.parametrize("op,volume", [
        ("union", 12.175), ("intersection", 3.825), ("difference", 4.175),
    ])
    def test_overlapping_boxes_in_general_position(self, op, volume):
        assert_correct(run(mesh_only(prism(2, 2, 2)),
                           moved(prism(2, 2, 2), (0.75, 0.3, 0.2)), op),
                       volume)

    @pytest.mark.parametrize("op,volume", [
        ("union", 16.0), ("intersection", 0.0), ("difference", 8.0),
    ])
    def test_boxes_sharing_a_face(self, op, volume):
        assert_correct(run(mesh_only(prism(2, 2, 2)),
                           moved(prism(2, 2, 2), (2.0, 0.0, 0.0)), op),
                       volume)

    @pytest.mark.parametrize("op,volume", [
        ("union", 12.0), ("intersection", 4.0), ("difference", 4.0),
    ])
    def test_boxes_with_coincident_side_faces(self, op, volume):
        assert_correct(run(mesh_only(prism(2, 2, 2)),
                           moved(prism(2, 2, 2), (1.0, 0.0, 0.0)), op),
                       volume)

    @pytest.mark.parametrize("op,volume", [
        ("union", 64.0), ("intersection", 8.0), ("difference", 56.0),
    ])
    def test_nested_boxes(self, op, volume):
        assert_correct(run(mesh_only(prism(4, 4, 4)),
                           mesh_only(prism(2, 2, 2)), op), volume)

    @pytest.mark.parametrize("op", ["union", "intersection", "difference"])
    def test_curved_results_are_closed(self, op):
        assert issolidclosed(run(mesh_only(sphere(1.0)),
                                 moved(sphere(1.0), (0.8, 0.1, 0.0)), op))
        assert issolidclosed(run(mesh_only(prism(2, 2, 2)),
                                 mesh_only(conic(0.6, 0.6, 3.0,
                                                 center=point(0, 0, -1.5))),
                                 op))

    def test_a_disjoint_intersection_is_empty(self):
        result = run(mesh_only(prism(2, 2, 2)),
                     moved(prism(2, 2, 2), (10.0, 0.0, 0.0)), "intersection")
        assert result[1] == []

    def test_results_are_deterministic(self):
        a = mesh_only(sphere(1.0))
        b = moved(sphere(1.0), (0.8, 0.1, 0.0))
        first, second = run(a, b, "union"), run(a, b, "union")
        assert first[1][0][1] == second[1][0][1]


class TestConventions:

    def test_construction_names_the_engine(self):
        result = run(mesh_only(prism(2, 2, 2)), mesh_only(prism(1, 1, 1)),
                     "difference")
        assert get_construction(result) == ["boolean", "manifold:difference"]

    def test_faces_are_flat_shaded(self):
        """One normal per triangle, as the other engines produce, so hard
        edges shade as edges rather than being averaged away."""
        surf = run(mesh_only(prism(2, 2, 2)),
                   moved(prism(2, 2, 2), (1.0, 1.0, 1.0)), "union")[1][0]
        for face in surf[3]:
            n0, n1, n2 = (surf[2][i] for i in face)
            assert n0 == n1 == n2

    def test_float64_precision_is_kept(self):
        """Mesh64, not the float32 Mesh: at a 1e6 offset float32 rounds to
        ~0.06, which would move every vertex."""
        far = 1.0e6 + 0.123456789
        a = moved(prism(2, 2, 2), (far, 0.0, 0.0))
        b = moved(prism(2, 2, 2), (far + 1.0, 0.0, 0.0))
        xs = sorted({round(v[0], 12) for v in run(a, b, "union")[1][0][1]})
        assert xs[0] == pytest.approx(far - 1.0, abs=1e-8)
        assert xs[-1] == pytest.approx(far + 2.0, abs=1e-8)

    def test_an_unknown_operation_is_refused(self):
        with pytest.raises(ValueError, match="unsupported"):
            manifold_engine.solid_boolean(mesh_only(prism(1, 1, 1)),
                                          mesh_only(prism(1, 1, 1)), "xor")


def open_box():
    """A cube with its top face removed: not a closed manifold."""
    cube = mesh_only(prism(2, 2, 2))
    return solid(cube[1][:-1], [], [])


class TestRefusal:

    def test_a_non_manifold_operand_is_named(self):
        with pytest.raises(manifold_engine.NotManifoldError,
                           match="second operand"):
            run(mesh_only(prism(2, 2, 2)), open_box(), "union")


class TestAutomaticSelection:
    """With no engine named, mesh-only booleans now prefer manifold3d."""

    def test_it_is_the_default_mesh_engine_when_installed(self):
        result = solid_boolean(mesh_only(prism(2, 2, 2)),
                               moved(prism(2, 2, 2), (1.0, 0.0, 0.0)),
                               "union")
        assert get_construction(result)[1] == "manifold:union"

    def test_the_environment_can_still_force_native(self, monkeypatch):
        monkeypatch.setenv("YAPCAD_MESH_BOOLEAN_ENGINE", "native")
        result = solid_boolean(mesh_only(prism(4, 4, 4)),
                               mesh_only(prism(2, 2, 2)), "difference")
        assert get_construction(result) == ["boolean", "difference"]

    def test_a_refused_input_falls_back_to_native_with_a_warning(self):
        with pytest.warns(RuntimeWarning, match="falling back to the native"):
            result = solid_boolean(mesh_only(prism(4, 4, 4)), open_box(),
                                   "union")
        assert get_construction(result) == ["boolean", "union"]

    def test_without_manifold3d_native_is_used(self, monkeypatch):
        monkeypatch.setattr(manifold_engine, "_manifold", None)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            result = solid_boolean(mesh_only(prism(4, 4, 4)),
                                   mesh_only(prism(2, 2, 2)), "difference")
        assert get_construction(result) == ["boolean", "difference"]

    def test_naming_it_without_manifold3d_explains(self, monkeypatch):
        monkeypatch.setattr(manifold_engine, "_manifold", None)
        with pytest.raises(RuntimeError, match="needs manifold3d"):
            run(mesh_only(prism(2, 2, 2)), mesh_only(prism(1, 1, 1)), "union")
