"""Placeholder helpers for the DSL fallback blocks."""

from __future__ import annotations

from ..involute_gear import generate_involute_profile


def involute_profile(teeth: int, module_mm: float, pressure_angle_deg: float = 20.0):
    """Return an involute profile as a list of (x, y) tuples."""
    return generate_involute_profile(teeth, module_mm, pressure_angle_deg)
