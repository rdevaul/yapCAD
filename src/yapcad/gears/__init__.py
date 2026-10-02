"""Parametric gear geometry."""

from .bevel import (
    StraightBevelGearGeometry,
    StraightBevelGearSpec,
    derive_straight_bevel_geometry,
    make_straight_bevel_gear,
    make_straight_bevel_gear_sdf,
    make_straight_bevel_pair,
    spherical_involute_point,
    tooth_section_points,
)

__all__ = [
    "StraightBevelGearGeometry",
    "StraightBevelGearSpec",
    "derive_straight_bevel_geometry",
    "make_straight_bevel_gear",
    "make_straight_bevel_gear_sdf",
    "make_straight_bevel_pair",
    "spherical_involute_point",
    "tooth_section_points",
]
