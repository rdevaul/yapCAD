"""Planar profiles and gears as fields (SDF gears and fasteners, step 1).

yapRover's differential uses ``miter_gear``, which built its teeth as an
OCC ruled loft and so could not exist without OpenCASCADE.  The field
version reproduces that construction -- the outer tooth section carried to
the pitch apex -- from a few-hundred-byte spec.
"""

import json
import math

import numpy as np
import pytest

from yapcad import sdf
from yapcad.brep import occ_available
from yapcad.construction import get_construction
from yapcad.dsl.runtime.builtins import get_builtin_registry
from yapcad.dsl.runtime.values import float_val, int_val, string_val
from yapcad.gears.bevel import (
    StraightBevelGearSpec,
    derive_straight_bevel_geometry,
    make_straight_bevel_gear_sdf,
    tooth_section_points,
)
from yapcad.geom3d import issolidclosed, volumeof
from yapcad.sdf import gears as sdf_gears

needs_occ = pytest.mark.skipif(not occ_available(),
                               reason="pythonocc-core is not available")

MITER = StraightBevelGearSpec(teeth=24, mate_teeth=24, outer_module_mm=1.5,
                              face_width_mm=8.0, backlash_mm=0.25,
                              bore_diameter_mm=8.4)


def brute_polygon(points, q):
    """Signed distance to a polygon by checking every segment, with the
    sign from a horizontal even-odd ray: independent of the node's code."""
    pts = np.asarray(points, float)
    a, b = pts, np.roll(pts, -1, axis=0)
    out = []
    for x, y in q[:, :2]:
        e = b - a
        t = np.clip(((x - a[:, 0]) * e[:, 0] + (y - a[:, 1]) * e[:, 1])
                    / (e * e).sum(axis=1), 0, 1)
        d = np.hypot(x - a[:, 0] - t * e[:, 0], y - a[:, 1] - t * e[:, 1])
        inside = False
        for (x0, y0), (x1, y1) in zip(a, b):
            if (y0 > y) != (y1 > y):
                if x < x0 + (y - y0) * (x1 - x0) / (y1 - y0):
                    inside = not inside
        out.append(-d.min() if inside else d.min())
    return np.array(out)


STAR = [(math.cos(k * math.pi / 5) * (5.0 if k % 2 == 0 else 2.0),
         math.sin(k * math.pi / 5) * (5.0 if k % 2 == 0 else 2.0))
        for k in range(10)]


# ---------------------------------------------------------------------------
# polygon
# ---------------------------------------------------------------------------


def test_polygon_is_the_exact_signed_distance():
    rng = np.random.default_rng(3)
    q = rng.uniform(-7, 7, (2000, 3))
    got = sdf.evaluate(sdf.polygon(STAR), q)
    assert got == pytest.approx(brute_polygon(STAR, q), abs=1e-12)


def test_the_symmetry_fold_changes_nothing():
    rng = np.random.default_rng(4)
    q = rng.uniform(-7, 7, (5000, 3))
    q[:3, :2] = 0.0                         # the origin, where angles fail
    full = sdf.evaluate(sdf.polygon(STAR), q)
    folded = sdf.evaluate(sdf.polygon(STAR, symmetry=5), q)
    assert folded == pytest.approx(full, abs=1e-12)
    assert folded[0] < 0.0


def test_polygon_winding_and_closure_do_not_matter():
    a = sdf.polygon(STAR)
    b = sdf.polygon(list(reversed(STAR)) + [STAR[-1]])
    q = np.random.default_rng(5).uniform(-7, 7, (500, 3))
    assert sdf.evaluate(a, q) == pytest.approx(sdf.evaluate(b, q))


def test_polygon_validation():
    with pytest.raises(sdf.SdfError, match="three"):
        sdf.polygon([(0, 0), (1, 0)])
    with pytest.raises(sdf.SdfError, match="no area"):
        sdf.polygon([(0, 0), (1, 0), (2, 0)])
    with pytest.raises(sdf.SdfError, match="symmetric"):
        sdf.polygon(STAR, symmetry=4)
    skewed = list(STAR)
    skewed[0] = (6.0, 0.0)
    with pytest.raises(sdf.SdfError, match="symmetric"):
        sdf.polygon(skewed, symmetry=5)


# ---------------------------------------------------------------------------
# extrude and apex_extrude
# ---------------------------------------------------------------------------


def test_an_extruded_square_is_a_box():
    square = sdf.polygon([(-2, -3), (2, -3), (2, 3), (-2, 3)])
    q = np.random.default_rng(6).uniform(-6, 6, (3000, 3))
    assert sdf.evaluate(sdf.extrude(square, 8.0), q) == pytest.approx(
        sdf.evaluate(sdf.box((4.0, 6.0, 8.0)), q), abs=1e-12)


def test_apex_extrude_scales_the_profile_toward_the_apex():
    square = sdf.polygon([(-4, -4), (4, -4), (4, 4), (-4, 4)])
    cone = sdf.apex_extrude(square, 10.0, 5.0, 10.0)
    # At z = 5.2 the section's half-width is 4 x 5.2 / 10 = 2.08.
    d = sdf.evaluate(cone, np.array([[2.0, 0, 5.5], [1.9, 0, 5.2],
                                     [2.2, 0, 5.2], [0, 0, 7.5]]))
    assert d[0] < 0 and d[1] < 0 < d[2] and d[3] < 0
    solid = sdf.to_solid(cone, resolution=64)
    assert issolidclosed(solid)
    # A frustum of a square pyramid: (h/3)(A1 + A2 + sqrt(A1 A2)).
    assert volumeof(solid) == pytest.approx(5 / 3 * (64 + 16 + 32), rel=2e-3)


def test_apex_extrude_validation():
    square = sdf.polygon([(-1, -1), (1, -1), (1, 1), (-1, 1)])
    with pytest.raises(sdf.SdfError):
        sdf.apex_extrude(square, 10.0, 6.0, 5.0)
    with pytest.raises(sdf.SdfError):
        sdf.apex_extrude(square, 10.0, 0.0, 5.0)


# ---------------------------------------------------------------------------
# straight bevel gear
# ---------------------------------------------------------------------------


def analytic_volume(spec, samples):
    """Every section is the outer one scaled by z / z_outer, so the volume
    is A_outer (z_o^3 - z_i^3) / (3 z_o^2), less the bore."""
    geometry = derive_straight_bevel_geometry(spec)
    pts = np.array(tooth_section_points(spec, geometry, cone_distance="outer",
                                        flank_samples=samples)[:-1])
    x, y = pts[:, 0], pts[:, 1]
    area = 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
    zo, zi = geometry.outer_plane_z_mm, geometry.inner_plane_z_mm
    teeth = area * (zo ** 3 - zi ** 3) / (3 * zo ** 2)
    return teeth - math.pi * (spec.bore_diameter_mm / 2) ** 2 * (zo - zi)


def test_the_gear_tree_is_a_spec_not_an_outline():
    node = sdf.straight_bevel_gear(MITER)
    doc = sdf.tree_to_json(node)
    assert len(json.dumps(doc)) < 1000
    assert sdf.tree_from_json(doc) == node
    assert node.csg_exact


def test_the_yaprover_miter_gear_meshes_closed_to_its_volume():
    solid = make_straight_bevel_gear_sdf(MITER)
    assert issolidclosed(solid)
    assert sdf.is_sdf_solid(solid)
    assert volumeof(solid) == pytest.approx(
        analytic_volume(MITER, sdf_gears.DEFAULT_FLANK_SAMPLES), rel=1e-3)


def test_an_invalid_spec_is_an_sdf_error():
    bad = StraightBevelGearSpec(teeth=4, mate_teeth=24, outer_module_mm=1.5,
                                face_width_mm=8.0)
    with pytest.raises(sdf.SdfError, match="teeth"):
        sdf.straight_bevel_gear(bad)


def test_dsl_miter_gear_needs_no_occ(monkeypatch):
    import yapcad.brep
    monkeypatch.setattr(yapcad.brep, "occ_available", lambda: False)
    miter = get_builtin_registry().get_function("miter_gear").implementation
    solid = miter(int_val(24), float_val(1.5), float_val(8.0),
                  float_val(20.0), float_val(0.25), float_val(8.4),
                  string_val("spherical_involute")).data
    assert issolidclosed(solid)
    tree = sdf.from_construction(get_construction(solid))
    assert tree == sdf.straight_bevel_gear(MITER)


@needs_occ
def test_the_field_matches_the_occ_gear():
    """Field and BREP are the same construction; they differ only between
    flank samples (straight segments against an interpolating spline)."""
    from yapcad.brep import brep_from_solid
    from yapcad.gears.bevel import make_straight_bevel_gear
    from yapcad.sdf import occ
    node = sdf.straight_bevel_gear(MITER)
    brep = brep_from_solid(make_straight_bevel_gear(MITER, flank_samples=17))
    tess = brep.tessellate(deflection=0.02)
    surface = np.array([p[:3] for p in tess[1]])
    assert np.abs(sdf.evaluate(node, surface)).max() < 5e-3
    assert occ.volume(occ.to_occ_shape(node)) == pytest.approx(
        analytic_volume(MITER, 17), rel=1e-3)
