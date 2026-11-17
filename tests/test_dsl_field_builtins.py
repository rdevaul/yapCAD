"""The DSL builtins that exist only as fields.

``capsule`` and ``torus`` are primitives with no mesh construction; the
smooth blends, ``offset_solid``, ``shell_solid``, ``twist_solid`` and the
lattices are arithmetic on signed distance fields.  Each takes and returns
SDF-authored solids, so a design that uses them builds with
``representation="sdf"``.
"""

import math

import pytest

from yapcad import sdf
from yapcad.construction import get_construction
from yapcad.dsl import compile_and_run
from yapcad.geom3d import issolidclosed, solidbbox, volumeof


def run(expr, rep="sdf", **kw):
    source = f"module fields\ncommand PART() -> solid:\n    emit {expr}\n"
    return compile_and_run(source, "PART", {}, representation=rep, **kw)


def part(expr, rep="sdf", **kw):
    result = run(expr, rep, **kw)
    assert result.success, result.error_message
    return result.geometry


def tree(solid):
    return sdf.from_construction(get_construction(solid))


def extent(solid):
    lo, hi = solidbbox(solid)
    return [hi[i] - lo[i] for i in range(3)]


# -- primitives --------------------------------------------------------------

def test_capsule_is_a_swept_sphere():
    solid = part("capsule(point(0.0, 0.0, 0.0), point(10.0, 0.0, 0.0), 2.0)",
                 sdf_cell_mm=0.1)
    assert sdf.is_sdf_solid(solid) and issolidclosed(solid)
    exact = math.pi * 4.0 * 10.0 + 4.0 / 3.0 * math.pi * 8.0
    assert volumeof(solid) == pytest.approx(exact, rel=5e-3)
    assert extent(solid) == pytest.approx([14.0, 4.0, 4.0], abs=0.05)
    assert tree(solid).kind == "capsule"


def test_torus_lies_in_the_xy_plane():
    solid = part("torus(10.0, 2.0)", sdf_cell_mm=0.15)
    assert sdf.is_sdf_solid(solid) and issolidclosed(solid)
    assert volumeof(solid) == pytest.approx(2.0 * math.pi ** 2 * 10.0 * 4.0,
                                            rel=5e-3)
    assert extent(solid) == pytest.approx([24.0, 24.0, 4.0], abs=0.05)


@pytest.mark.parametrize("expr", [
    "capsule(point(0.0, 0.0, 0.0), point(6.0, 0.0, 0.0), 1.5)",
    "torus(6.0, 1.5)",
])
def test_field_primitives_are_fields_in_mesh_mode_too(expr):
    """There is no mesh construction to fall back to, so the default
    representation gets the field as well."""
    solid = part(expr, rep="mesh")
    assert sdf.is_sdf_solid(solid) and issolidclosed(solid)


# -- smooth blends -----------------------------------------------------------

TWO_BALLS = "sphere(5.0), translate(sphere(5.0), 7.0, 0.0, 0.0)"


def test_smooth_union_fills_the_seam():
    hard = part(f"union({TWO_BALLS})", sdf_cell_mm=0.2)
    soft = part(f"smooth_union({TWO_BALLS}, 2.0)", sdf_cell_mm=0.2)
    assert sdf.is_sdf_solid(soft) and issolidclosed(soft)
    assert tree(soft).kind == "smooth_union"
    assert tree(soft).p["radius"] == 2.0
    # The blend adds a fillet in the crease and nothing anywhere else.
    assert volumeof(hard) < volumeof(soft) < volumeof(hard) * 1.05
    assert extent(soft) == pytest.approx(extent(hard), abs=0.05)


def test_smooth_difference_and_intersection_round_the_new_edge():
    block = "box(10.0, 10.0, 10.0), sphere(6.0)"
    cut_hard = part(f"difference({block})", sdf_cell_mm=0.2)
    cut_soft = part(f"smooth_difference({block}, 1.0)", sdf_cell_mm=0.2)
    both_hard = part(f"intersection({block})", sdf_cell_mm=0.2)
    both_soft = part(f"smooth_intersection({block}, 1.0)", sdf_cell_mm=0.2)
    for solid in (cut_soft, both_soft):
        assert sdf.is_sdf_solid(solid) and issolidclosed(solid)
    # Rounding a convex edge removes material from it.
    assert volumeof(cut_soft) < volumeof(cut_hard)
    assert volumeof(both_soft) < volumeof(both_hard)
    assert tree(cut_soft).kind == "smooth_subtract"
    assert tree(both_soft).kind == "smooth_intersect"


def test_a_blend_keeps_the_first_operands_metadata():
    source = '''
module fields
@meta(material="PLA")
command BODY() -> solid:
    emit sphere(5.0)

command PART() -> solid:
    emit smooth_union(BODY(), translate(sphere(4.0), 6.0, 0.0, 0.0), 1.5)
'''
    result = compile_and_run(source, "PART", {}, representation="sdf")
    assert result.success, result.error_message
    from yapcad.metadata import get_solid_metadata
    assert get_solid_metadata(result.geometry).get("material") == "PLA"


# -- offset, shell, twist ----------------------------------------------------

def test_offset_solid_grows_and_shrinks():
    grown = part("offset_solid(box(10.0, 10.0, 10.0), 1.0)", sdf_cell_mm=0.2)
    shrunk = part("offset_solid(box(10.0, 10.0, 10.0), -1.0)", sdf_cell_mm=0.2)
    assert extent(grown) == pytest.approx([12.0] * 3, abs=0.05)
    assert extent(shrunk) == pytest.approx([8.0] * 3, abs=0.05)
    # Growing a box rounds its edges and corners; shrinking keeps them.
    rounded = (1000.0 + 6 * 100.0 + 3 * math.pi * 10.0
               + 4.0 / 3.0 * math.pi)
    assert volumeof(grown) == pytest.approx(rounded, rel=5e-3)
    assert volumeof(shrunk) == pytest.approx(512.0, rel=5e-3)


def test_offset_solid_is_a_uniform_running_clearance():
    """What a print-in-place hinge needs: a bore that clears its pin by the
    same distance everywhere, bevels included."""
    pin = "chamfer(cylinder(4.0, 6.0), 1.0)"
    bore = tree(part(f"offset_solid({pin}, 0.4)"))
    shaft = tree(part(pin))
    on_bore = [[4.4, 0.0, 3.0], [0.0, 0.0, 6.4], [0.0, -4.4, 3.0]]
    assert sdf.evaluate(bore, on_bore) == pytest.approx([0.0] * 3, abs=1e-9)
    assert sdf.evaluate(shaft, on_bore) == pytest.approx([0.4] * 3, abs=1e-9)


def test_shell_solid_straddles_the_surface():
    solid = part("shell_solid(sphere(10.0), 1.0)", sdf_cell_mm=0.25)
    assert sdf.is_sdf_solid(solid) and issolidclosed(solid)
    exact = 4.0 / 3.0 * math.pi * (10.5 ** 3 - 9.5 ** 3)
    assert volumeof(solid) == pytest.approx(exact, rel=1e-2)


def test_twist_solid_turns_sections_by_degrees_per_unit():
    solid = part("twist_solid(box(10.0, 4.0, 20.0), 9.0)", sdf_cell_mm=0.25)
    assert sdf.is_sdf_solid(solid) and issolidclosed(solid)
    # A twist moves material round the axis without adding or removing any.
    assert volumeof(solid) == pytest.approx(800.0, rel=1e-2)
    twisted = tree(solid)
    assert twisted.p["rate"] == pytest.approx(math.radians(9.0))
    # 9 degrees per unit puts the section at z = 10 a quarter turn round:
    # the bar's long axis lies along y there, along x at z = 0.
    probes = [[4.5, 0.0, 0.0], [0.0, 4.5, 0.0],
              [0.0, 4.5, 9.9], [4.5, 0.0, 9.9]]
    inside = [v < 0.0 for v in sdf.evaluate(twisted, probes)]
    assert inside == [True, False, True, False]


# -- lattices ----------------------------------------------------------------

@pytest.mark.parametrize("name", ["gyroid_lattice", "schwarz_p_lattice"])
def test_a_lattice_fills_the_part_and_stays_inside_it(name):
    solid = part(f"{name}(box(16.0, 16.0, 16.0), 8.0, 0.5)")
    assert sdf.is_sdf_solid(solid) and issolidclosed(solid)
    assert 0.1 * 4096.0 < volumeof(solid) < 0.9 * 4096.0
    assert extent(solid) == pytest.approx([16.0] * 3, abs=0.1)
    lattice = tree(solid)
    assert lattice.kind == "intersect"
    assert {child.kind for child in lattice.children} >= {
        name.replace("_lattice", "")}


# -- refusing what is not a field --------------------------------------------

@pytest.mark.parametrize("expr,name", [
    ("smooth_union(box(4.0, 4.0, 4.0), sphere(3.0), 1.0)", "smooth_union"),
    ("offset_solid(box(4.0, 4.0, 4.0), 0.5)", "offset_solid"),
    ("gyroid_lattice(box(8.0, 8.0, 8.0), 4.0, 0.5)", "gyroid_lattice"),
])
def test_a_mesh_operand_is_refused_with_what_to_do(expr, name):
    result = run(expr, rep="mesh")
    assert not result.success
    assert name in result.error_message
    assert "not SDF-authored" in result.error_message
    assert "representation" in result.error_message


def test_the_second_operand_is_named_when_it_is_the_mesh():
    expr = ("smooth_union(sphere(4.0), "
            "extrude(rectangle(4.0, 4.0), 4.0), 1.0)")
    result = run(expr, rep="sdf")
    assert not result.success
    assert "the second solid" in result.error_message


# -- deferred meshing ---------------------------------------------------------

@pytest.mark.parametrize("expr", [
    "capsule(point(0.0, 0.0, 0.0), point(6.0, 0.0, 0.0), 1.5)",
    "smooth_union(sphere(4.0), translate(sphere(4.0), 5.0, 0.0, 0.0), 1.5)",
    "difference(cylinder(7.0, 6.0), "
    "offset_solid(chamfer(cylinder(4.0, 6.0), 1.0), 0.35))",
    "twist_solid(shell_solid(torus(8.0, 2.0), 0.8), 3.0)",
])
def test_field_builtins_defer_their_meshes(expr):
    """With ``sdf_mesh=False`` nothing is meshed on the way through, the
    result carries the same tree and a plan, and realizing it gives the
    part an eager build gives."""
    eager = part(expr, sdf_simplify=False)
    deferred = part(expr, sdf_mesh=False)
    assert sdf.is_deferred(deferred)
    assert not deferred[1]
    assert tree(deferred) == tree(eager)
    sdf.realize_all(deferred)
    assert not sdf.is_deferred(deferred) and issolidclosed(deferred)
    assert volumeof(deferred) == pytest.approx(volumeof(eager), rel=1e-2)
