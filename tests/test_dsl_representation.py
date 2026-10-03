"""Building a DSL design as fields: ``representation="sdf"``.

The same source builds either way; with ``"sdf"`` the primitives,
fasteners and gears are SDF-authored, placed exactly as the mesh primitives
are, so booleans combine fields and ``fillet`` rounds them without OCC.
"""

import math

import pytest

from yapcad import sdf
from yapcad.dsl import compile_and_run
from yapcad.dsl.runtime import representation
from yapcad.geom3d import issolidclosed, solidbbox, volumeof


def run(source, rep="sdf", **kw):
    result = compile_and_run(source, "PART", {}, representation=rep, **kw)
    assert result.success, result.error_message
    geometry = result.geometry
    if isinstance(geometry, list) and geometry and geometry[0] != "solid":
        geometry = geometry[0]
    return geometry


def part(expr):
    return f"module rep\ncommand PART() -> solid:\n    emit {expr}\n"


def bbox(solid):
    lo, hi = solidbbox(solid)
    return [round(v, 6) for v in lo[:3]], [round(v, 6) for v in hi[:3]]


PRIMITIVES = [
    ("box(4.0, 6.0, 8.0)", 4 * 6 * 8),
    ("cylinder(3.0, 5.0)", math.pi * 9 * 5),
    ("cone(3.0, 1.0, 6.0)", math.pi * 6 / 3 * (9 + 3 + 1)),
    ("tube(6.0, 1.0, 4.0)", math.pi * (9 - 4) * 4),
    # sphere() passes its argument to geom3d_util.sphere as a diameter.
    ("sphere(6.0)", 4 / 3 * math.pi * 27),
    ("spherical_shell(8.0, 1.0)", 4 / 3 * math.pi * (64 - 27)),
]


@pytest.mark.parametrize("expr,volume", PRIMITIVES)
def test_sdf_primitives_are_placed_as_the_mesh_ones(expr, volume):
    field = run(part(expr), "sdf", sdf_cell_mm=0.1)
    mesh = run(part(expr), "mesh")
    assert sdf.is_sdf_solid(field)
    assert issolidclosed(field)
    assert volumeof(field) == pytest.approx(volume, rel=5e-3)
    (flo, fhi), (mlo, mhi) = bbox(field), bbox(mesh)
    assert flo == pytest.approx(mlo, abs=0.06)
    assert fhi == pytest.approx(mhi, abs=0.06)


def test_a_filleted_design_builds_as_fields_without_occ(monkeypatch):
    import yapcad.brep
    monkeypatch.setattr(yapcad.brep, "occ_available", lambda: False)
    source = """
module rep
command PART() -> solid:
    let plate: solid = fillet(box(40.0, 30.0, 6.0), 1.0)
    let hole: solid = translate(cylinder(4.0, 10.0), 8.0, 0.0, -5.0)
    emit difference(plate, hole)
"""
    solid = run(source)
    r, a, b, c = 1.0, 40.0, 30.0, 6.0
    x, y, z = a - 2 * r, b - 2 * r, c - 2 * r
    plate = (x * y * z + 2 * r * (x * y + y * z + x * z)
             + math.pi * r * r * (x + y + z) + 4 / 3 * math.pi * r ** 3)
    assert issolidclosed(solid)
    assert volumeof(solid) == pytest.approx(plate - math.pi * 16 * c,
                                            rel=2e-3)


def test_fasteners_and_gears_choose_their_fields():
    nut = run(part('metric_hex_nut("M8")'))
    gear = run(part("miter_gear(teeth=24, outer_module_mm=1.5, "
                    "face_width_mm=8.0, bore_diameter_mm=8.4)"))
    for solid in (nut, gear):
        assert sdf.is_sdf_solid(solid) and issolidclosed(solid)


def test_mesh_remains_the_default():
    assert not sdf.is_sdf_solid(run(part("box(1.0, 1.0, 1.0)"), None))


def test_the_environment_sets_the_default(monkeypatch):
    monkeypatch.setenv("YAPCAD_DSL_REPRESENTATION", "sdf")
    assert sdf.is_sdf_solid(run(part("box(1.0, 1.0, 1.0)"), None))


def test_an_unknown_representation_fails_cleanly():
    result = compile_and_run(part("box(1.0, 1.0, 1.0)"), "PART", {},
                             representation="voxels")
    assert not result.success and "representation" in result.error_message


def test_the_setting_does_not_leak_out_of_a_call():
    run(part("box(1.0, 1.0, 1.0)"), "sdf")
    assert representation.current() == "mesh"
    with representation.using("sdf", 0.5):
        assert representation.is_sdf()
    assert not representation.is_sdf()
