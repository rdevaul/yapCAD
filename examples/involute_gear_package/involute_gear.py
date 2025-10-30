"""Utility functions for generating involute spur gears."""

from __future__ import annotations

import math
from typing import List, Tuple

from yapcad.geom import point, poly
from yapcad.geom3d import translatesolid, rotatesolid
from yapcad.geom3d_util import extrude, poly2surfaceXY
from yapcad.geom3d import point as point3d

from yapcad.contrib.figgear import make_gear_figure as _figgear_make_gear


def _arc_points(
    radius: float,
    start_angle: float,
    end_angle: float,
    steps: int,
    *,
    clockwise: bool = False,
) -> List[Tuple[float, float]]:
    if steps < 2:
        steps = 2
    if clockwise:
        if end_angle > start_angle:
            end_angle -= 2 * math.pi
    else:
        if end_angle < start_angle:
            end_angle += 2 * math.pi
    samples = []
    for i in range(steps):
        t = start_angle + (end_angle - start_angle) * i / (steps - 1)
        samples.append((radius * math.cos(t), radius * math.sin(t)))
    return samples


def generate_involute_profile(
    teeth: int,
    module_mm: float,
    pressure_angle_deg: float = 20.0,
    *,
    tooth_resolution: int = 12,
    tip_resolution: int = 6,
    root_resolution: int = 5,
) -> List[Tuple[float, float]]:
    if teeth < 6:
        raise ValueError("teeth must be >= 6 for a usable spur gear")
    if module_mm <= 0:
        raise ValueError("module must be positive")

    figure_points, _ = _figgear_make_gear(
        m=module_mm,
        z=teeth,
        alpha_deg=pressure_angle_deg,
        bottom_type="line",
        #bottom_type="spline",
        involute_step=max(module_mm / max(tooth_resolution, 8), 0.05),
        spline_division_num=max(root_resolution * 4, 20),
    )

    profile = [(float(x), float(y)) for x, y in figure_points]
    if profile[0] != profile[-1]:
        profile.append(profile[0])

    area = 0.0
    for (x0, y0), (x1, y1) in zip(profile, profile[1:]):
        area += x0 * y1 - x1 * y0
    if area < 0:
        profile.reverse()
    return profile


def generate_involute_spur(teeth: int, module_mm: float, face_width_mm: float,
                            pressure_angle_deg: float = 20.0) -> Tuple[object, object]:
    outline_pts = generate_involute_profile(teeth, module_mm, pressure_angle_deg)
    outline = [point(x, y) for x, y in outline_pts]
    outline_geom = poly(outline)
    surface, _ = poly2surfaceXY(outline_geom)
    gear_solid = extrude(surface, face_width_mm)
    if len(gear_solid) > 2:
        gear_solid[2] = []
    if len(gear_solid) > 3:
        gear_solid[3] = []
    return surface, gear_solid


def position_gear(solid_geom: list, centre: Tuple[float, float, float], spin_deg: float = 0.0) -> list:
    placed = translatesolid(solid_geom, point3d(*centre))
    if spin_deg:
        placed = rotatesolid(placed, spin_deg, axis=point3d(0, 0, 1), cent=point3d(*centre))
    return placed
