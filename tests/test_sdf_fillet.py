"""Fillets and compounds of SDF-authored solids, without OCC.

yapRover could not be built without OpenCASCADE because ``fillet`` needed a
BREP and ``compound`` threw away the fields of its operands, sending every
later boolean to the mesh engine.  These pin the field versions of both.
"""

import math

import numpy as np
import pytest

from yapcad import sdf
from yapcad.brep import occ_available
from yapcad.construction import get_construction
from yapcad.dsl.runtime.builtins import get_builtin_registry
from yapcad.dsl.runtime.values import float_val, solid_val
from yapcad.geom import point
from yapcad.geom3d import issolidclosed, solid_boolean, translatesolid, \
    volumeof
from yapcad.geom3d_util import prism

needs_occ = pytest.mark.skipif(not occ_available(),
                               reason="pythonocc-core is not available")


def tree_of(solid):
    return sdf.from_construction(get_construction(solid))


def rounded_box_volume(a, b, c, r):
    x, y, z = a - 2 * r, b - 2 * r, c - 2 * r
    return (x * y * z + 2 * r * (x * y + y * z + x * z)
            + math.pi * r * r * (x + y + z) + 4.0 / 3.0 * math.pi * r ** 3)


def rounded_cylinder_volume(radius, height, e):
    # The core slab, plus two caps: a disc of radius (radius - e) and a
    # quarter-torus round (Pappus: quarter-disc area times its path).
    core = math.pi * radius ** 2 * (height - 2 * e)
    quarter = math.pi * e * e / 4.0
    centroid = radius - e + 4.0 * e / (3.0 * math.pi)
    cap = math.pi * (radius - e) ** 2 * e + 2 * math.pi * centroid * quarter
    return core + 2 * cap


# ---------------------------------------------------------------------------
# rounded_cylinder
# ---------------------------------------------------------------------------


def test_rounded_cylinder_is_an_exact_distance():
    node = sdf.rounded_cylinder(10.0, 8.0, 2.0)
    pts = np.array([[0, 0, 0], [12, 0, 0], [0, 0, 7], [0, 0, 4],
                    [10, 0, 0], [8 + 2 * math.cos(1.0), 0,
                                 2 + 2 * math.sin(1.0)]], dtype=float)
    d = sdf.evaluate(node, pts)
    assert d == pytest.approx([-4.0, 2.0, 3.0, 0.0, 0.0, 0.0], abs=1e-12)
    assert node.exact and node.lipschitz == 1.0 and node.csg_exact


def test_rounded_cylinder_meshes_to_its_volume():
    solid = sdf.to_solid(sdf.rounded_cylinder(65.0, 36.0, 2.5),
                         resolution=128)
    assert issolidclosed(solid)
    assert volumeof(solid) == pytest.approx(
        rounded_cylinder_volume(65.0, 36.0, 2.5), rel=2e-4)


def test_rounded_cylinder_validates_its_radius():
    assert sdf.rounded_cylinder(5.0, 4.0, 0.0).kind == "cylinder"
    with pytest.raises(sdf.SdfError):
        sdf.rounded_cylinder(5.0, 4.0, 2.5)
    with pytest.raises(sdf.SdfError):
        sdf.rounded_cylinder(2.0, 10.0, 2.5)
    # Consuming a whole cap or side is a valid field but not replayable.
    assert not sdf.rounded_cylinder(5.0, 4.0, 2.0).csg_exact


# ---------------------------------------------------------------------------
# fillet on trees
# ---------------------------------------------------------------------------


def test_fillet_of_a_box_is_a_rounded_box():
    out = sdf.fillet(sdf.box((40, 16, 28)), 7.0)
    assert out == sdf.rounded_box((40, 16, 28), 7.0)


def test_fillet_of_a_cylinder_is_a_rounded_cylinder():
    out = sdf.fillet(sdf.cylinder(65.0, 36.0), 2.5)
    assert out == sdf.rounded_cylinder(65.0, 36.0, 2.5)


def test_fillet_follows_rigid_transforms():
    part = sdf.translate(sdf.rotate(sdf.box(10.0), (0, 0, 1), 30.0),
                         (5, 0, 0))
    out = sdf.fillet(part, 1.0)
    assert out.kind == "transform"
    assert out.children[0].children[0] == sdf.rounded_box(10.0, 1.0)


def test_a_uniform_scale_divides_the_radius():
    out = sdf.fillet(sdf.scale(sdf.box(10.0), 2.0), 4.0)
    assert out.children[0] == sdf.rounded_box(10.0, 2.0)


def test_a_non_uniform_scale_is_refused():
    with pytest.raises(sdf.SdfError, match="non-uniform"):
        sdf.fillet(sdf.scale(sdf.box(10.0), (1.0, 2.0, 1.0)), 1.0)


def test_combined_fields_are_refused_not_approximated():
    part = sdf.union(sdf.box(10.0), sdf.translate(sdf.box(10.0), (10, 0, 0)))
    with pytest.raises(sdf.SdfError, match="combined field"):
        sdf.fillet(part, 1.0)


def test_smooth_primitives_are_unchanged():
    for node in (sdf.sphere(3.0), sdf.torus(5.0, 1.0),
                 sdf.rounded_box(10.0, 1.0)):
        assert sdf.fillet(node, 0.5) == node


def test_a_radius_too_large_for_the_part_is_refused():
    with pytest.raises(sdf.SdfError, match="half the smallest edge"):
        sdf.fillet(sdf.box((10, 4, 10)), 2.5)
    with pytest.raises(sdf.SdfError):
        sdf.fillet(sdf.box(10.0), 0.0)


def test_fillet_rounds_each_body_of_a_compound():
    part = sdf.compound(sdf.box(4.0), sdf.translate(sdf.cylinder(1, 2),
                                                    (5, 0, 0)))
    out = sdf.fillet(part, 0.5)
    assert out.kind == "compound"
    assert out.children[0] == sdf.rounded_box(4.0, 0.5)
    assert out.children[1].children[0] == sdf.rounded_cylinder(1, 2, 0.5)


# ---------------------------------------------------------------------------
# fillet on solids
# ---------------------------------------------------------------------------


def test_fillet_solid_matches_the_analytic_volume():
    solid = sdf.to_solid(sdf.box((40, 16, 28)), resolution=80)
    out = sdf.fillet_solid(solid, 7.0)
    assert issolidclosed(out)
    assert tree_of(out) == sdf.rounded_box((40, 16, 28), 7.0)
    assert volumeof(out) == pytest.approx(
        rounded_box_volume(40, 16, 28, 7.0), rel=2e-4)


def test_fillet_solid_refuses_a_mesh_solid():
    with pytest.raises(sdf.SdfError, match="not SDF-authored"):
        sdf.fillet_solid(prism(1, 1, 1), 0.1)


@needs_occ
def test_a_filleted_primitive_still_replays_to_an_exact_brep():
    from yapcad.sdf import occ
    for node, exact in (
            (sdf.rounded_box((40, 16, 28), 7.0),
             rounded_box_volume(40, 16, 28, 7.0)),
            (sdf.rounded_cylinder(65.0, 36.0, 2.5),
             rounded_cylinder_volume(65.0, 36.0, 2.5))):
        assert occ.volume(occ.to_occ_shape(node)) == pytest.approx(
            exact, rel=1e-6)


# ---------------------------------------------------------------------------
# compound
# ---------------------------------------------------------------------------


PAIR = sdf.compound(sdf.box(10.0), sdf.translate(sdf.box(10.0), (12, 0, 0)))


def test_compound_meshes_each_body_separately():
    solid = sdf.to_solid(PAIR, resolution=44)
    assert len(solid[1]) == 2
    assert issolidclosed(solid)
    assert volumeof(solid) == pytest.approx(2000.0, rel=1e-6)


def test_touching_bodies_are_not_fused():
    touching = sdf.compound(sdf.box(10.0),
                            sdf.translate(sdf.box(10.0), (10, 0, 0)))
    assert len(sdf.to_solid(touching, resolution=40)[1]) == 2


def test_a_transformed_compound_remeshes_as_separate_bodies():
    moved = translatesolid(sdf.to_solid(PAIR, resolution=44),
                           point(5, 0, 0))
    assert tree_of(moved).kind == "transform"
    again = sdf.to_solid(tree_of(moved), resolution=44)
    assert len(again[1]) == 2


def test_compound_combines_as_a_union_in_booleans():
    pair = sdf.to_solid(PAIR, resolution=44)
    # A 2 x 2 bar right through both bodies removes 2 x 2 x 10 from each.
    bar = sdf.to_solid(sdf.translate(sdf.box((40, 2, 2)), (6, 0, 0)),
                       resolution=40)
    cut = solid_boolean(pair, bar, "difference")
    assert sdf.is_sdf_solid(cut)
    assert volumeof(cut) == pytest.approx(2000.0 - 2 * 40, rel=1e-3)


def test_compound_survives_serialisation():
    assert sdf.tree_from_json(sdf.tree_to_json(PAIR)) == PAIR


@needs_occ
def test_compound_replays_as_an_occ_compound():
    from yapcad.sdf import occ
    assert occ.volume(occ.to_occ_shape(PAIR)) == pytest.approx(2000.0)


# ---------------------------------------------------------------------------
# DSL builtins
# ---------------------------------------------------------------------------


def builtin(name):
    return get_builtin_registry().get_function(name).implementation


def test_dsl_fillet_of_an_sdf_solid_needs_no_occ(monkeypatch):
    import yapcad.brep
    monkeypatch.setattr(yapcad.brep, "occ_available", lambda: False)
    solid = sdf.to_solid(sdf.cylinder(18.0, 20.0), resolution=48)
    out = builtin("fillet")(solid_val(solid), float_val(1.5)).data
    assert tree_of(out) == sdf.rounded_cylinder(18.0, 20.0, 1.5)
    assert issolidclosed(out)


def test_dsl_fillet_of_a_mesh_solid_explains_itself(monkeypatch):
    import yapcad.brep
    monkeypatch.setattr(yapcad.brep, "occ_available", lambda: False)
    with pytest.raises(RuntimeError, match="mesh-only solid"):
        builtin("fillet")(solid_val(prism(1, 1, 1)), float_val(0.1))


def test_dsl_compound_keeps_the_fields():
    a = sdf.to_solid(sdf.box(10.0), resolution=20)
    b = sdf.to_solid(sdf.translate(sdf.box(10.0), (12, 0, 0)), resolution=20)
    out = builtin("compound")(solid_val(a), solid_val(b)).data
    assert sdf.is_sdf_solid(out)
    assert tree_of(out) == PAIR
    assert len(out[1]) == 2                 # the operands' own meshes


def test_dsl_compound_of_mixed_operands_is_unchanged():
    a = sdf.to_solid(sdf.box(10.0), resolution=20)
    out = builtin("compound")(solid_val(a), solid_val(prism(1, 1, 1))).data
    assert not sdf.is_sdf_solid(out)
