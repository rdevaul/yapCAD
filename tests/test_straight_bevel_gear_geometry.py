"""Pure geometry contract for spherical-involute straight bevel gears."""

import math

import pytest

from yapcad.gears.bevel import (
    StraightBevelGearSpec,
    derive_straight_bevel_geometry,
    make_straight_bevel_gear,
    make_straight_bevel_pair,
    spherical_involute_point,
    tooth_section_points,
)
from yapcad.brep import brep_from_solid, occ_available
from yapcad.geom3d import solidbbox, volumeof


def _miter_spec(**changes):
    values = dict(
        teeth=24,
        mate_teeth=24,
        outer_module_mm=1.5,
        face_width_mm=8.0,
        shaft_angle_deg=90.0,
        pressure_angle_deg=20.0,
        backlash_mm=0.25,
        bore_diameter_mm=8.3,
    )
    values.update(changes)
    return StraightBevelGearSpec(**values)


def test_equal_tooth_miter_geometry_is_45_degrees():
    geometry = derive_straight_bevel_geometry(_miter_spec())
    assert math.degrees(geometry.pitch_cone_angle_rad) == pytest.approx(45.0)
    assert geometry.outer_pitch_radius_mm == pytest.approx(18.0)
    assert geometry.outer_cone_distance_mm == pytest.approx(18.0 * math.sqrt(2.0))
    assert geometry.inner_cone_distance_mm == pytest.approx(
        18.0 * math.sqrt(2.0) - 8.0
    )


def test_unequal_pair_pitch_angles_sum_to_shaft_angle():
    pinion = derive_straight_bevel_geometry(_miter_spec(teeth=16, mate_teeth=32))
    gear = derive_straight_bevel_geometry(_miter_spec(teeth=32, mate_teeth=16))
    assert pinion.pitch_cone_angle_rad + gear.pitch_cone_angle_rad == pytest.approx(
        math.pi / 2.0
    )


def test_spherical_involute_starts_on_base_cone_and_moves_outward():
    geometry = derive_straight_bevel_geometry(_miter_spec())
    start = spherical_involute_point(geometry.base_cone_angle_rad, 0.0)
    later = spherical_involute_point(geometry.base_cone_angle_rad, 0.35)
    assert math.acos(start[2]) == pytest.approx(geometry.base_cone_angle_rad)
    assert math.acos(later[2]) > geometry.base_cone_angle_rad
    assert sum(value * value for value in later) == pytest.approx(1.0)


def test_section_is_closed_planar_and_has_two_flanks_per_tooth():
    spec = _miter_spec()
    geometry = derive_straight_bevel_geometry(spec)
    points = tooth_section_points(spec, geometry, cone_distance="outer", flank_samples=7)
    assert points[0] == pytest.approx(points[-1])
    assert len(points) >= spec.teeth * 14 + 1
    assert max(point[2] for point in points) == pytest.approx(
        min(point[2] for point in points)
    )


@pytest.mark.parametrize("mode", ["octoid", "gleason"])
def test_reserved_core_modes_raise_not_implemented(mode):
    with pytest.raises(NotImplementedError, match=mode):
        derive_straight_bevel_geometry(_miter_spec(generation_type=mode))


@pytest.mark.parametrize(
    "changes",
    [
        {"teeth": 5},
        {"mate_teeth": 0},
        {"outer_module_mm": 0.0},
        {"shaft_angle_deg": 180.0},
        {"pressure_angle_deg": 0.0},
        {"face_width_mm": 30.0},
        {"backlash_mm": 3.0},
        {"bore_diameter_mm": -1.0},
        {"generation_type": "cycloidal"},
    ],
)
def test_invalid_specs_are_rejected(changes):
    with pytest.raises((TypeError, ValueError)):
        derive_straight_bevel_geometry(_miter_spec(**changes))


@pytest.mark.skipif(not occ_available(), reason="requires pythonocc-core")
def test_miter_gear_is_a_valid_brep_with_a_through_bore():
    from OCC.Core.BRepCheck import BRepCheck_Analyzer

    gear = make_straight_bevel_gear(_miter_spec(), flank_samples=7)
    brep = brep_from_solid(gear)
    assert brep is not None
    assert BRepCheck_Analyzer(brep.shape).IsValid()
    assert volumeof(gear) > 0.0
    bbox = solidbbox(gear)
    assert bbox[1][2] - bbox[0][2] == pytest.approx(
        8.0 * math.cos(math.pi / 4.0), abs=0.1
    )


@pytest.mark.skipif(not occ_available(), reason="requires pythonocc-core")
def test_miter_pair_remains_non_intersecting_through_a_tooth_pitch():
    from yapcad.geom3d import solid_boolean

    spec = _miter_spec(
        teeth=8,
        mate_teeth=8,
        outer_module_mm=2.0,
        face_width_mm=4.0,
        backlash_mm=0.2,
        bore_diameter_mm=0.0,
    )
    for angle in (0.0, 5.0, 11.0, 17.0, 30.0, 45.0):
        driver, mate = make_straight_bevel_pair(
            spec, driver_angle_deg=angle, flank_samples=5
        )
        try:
            overlap = volumeof(
                solid_boolean(driver, mate, "intersection", engine="occ")
            )
        except RuntimeError as error:
            if "produced no solids" not in str(error):
                raise
            overlap = 0.0
        assert overlap <= 1.0e-6
