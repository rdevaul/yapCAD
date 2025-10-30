"""
Vendored helpers derived from the MIT-licensed ``figgear`` project.

Original project: https://github.com/chromia/figgear

MIT License
-----------

Copyright (c) chromia

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

This module provides ``make_gear_figure`` with a compatible signature but
implements the limited subset of functionality that yapCAD requires.  The
implementation relies only on the Python standard library to avoid adding
heavy runtime dependencies such as SciPy.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable, List, Sequence, Tuple, Dict

Point = Tuple[float, float]
PointList = List[Point]


@dataclass(frozen=True)
class _GearParameters:
    module: float
    teeth: int
    pressure_angle_deg: float
    involute_step: float
    spline_division_num: int
    bottom_type: str


def _inv(alpha: float) -> float:
    """Return the involute function for ``alpha``."""
    return math.tan(alpha) - alpha


def _catmull_rom(
    control: Sequence[Point],
    *,
    division_num: int,
) -> Iterable[Point]:
    """Yield interpolated points between control[1] and control[2]."""

    if division_num <= 1:
        return []

    p0, p1, p2, p3 = control
    for i in range(1, division_num):
        t = i / division_num
        t2 = t * t
        t3 = t2 * t
        cx = (
            0.5
            * (
                (2 * p1[0])
                + (-p0[0] + p2[0]) * t
                + (2 * p0[0] - 5 * p1[0] + 4 * p2[0] - p3[0]) * t2
                + (-p0[0] + 3 * p1[0] - 3 * p2[0] + p3[0]) * t3
            )
        )
        cy = (
            0.5
            * (
                (2 * p1[1])
                + (-p0[1] + p2[1]) * t
                + (2 * p0[1] - 5 * p1[1] + 4 * p2[1] - p3[1]) * t2
                + (-p0[1] + 3 * p1[1] - 3 * p2[1] + p3[1]) * t3
            )
        )
        yield (cx, cy)


def _add_bottom_points_line(points: PointList, new_points: Sequence[Point]) -> None:
    """Append the interior points of the tooth root as straight segments."""
    points.extend(new_points[1:-1])


def _add_bottom_points_spline(
    points: PointList,
    new_points: Sequence[Point],
    division_num: int,
) -> None:
    """Append interpolated points along the tooth root using a cubic spline."""
    if division_num <= 1:
        _add_bottom_points_line(points, new_points)
        return
    points.extend(_catmull_rom(new_points, division_num=division_num))


def _ensure_closed(points: PointList) -> None:
    if points and points[0] != points[-1]:
        points.append(points[0])


def make_gear_figure(
    m: float,
    z: int,
    alpha_deg: float,
    bottom_type: str,
    **kwargs,
) -> Tuple[PointList, Dict[str, float]]:
    """Generate a 2D involute spur gear profile."""

    if z <= 0:
        raise ValueError("gear must have a positive tooth count")
    if m <= 0:
        raise ValueError("module must be positive")
    bottom_type = bottom_type.lower()
    if bottom_type not in {"spline", "line"}:
        raise ValueError("bottom_type must be 'spline' or 'line'")

    params = _GearParameters(
        module=m,
        teeth=z,
        pressure_angle_deg=alpha_deg,
        involute_step=float(kwargs.get("involute_step", 0.5)),
        spline_division_num=int(kwargs.get("spline_division_num", 50)),
        bottom_type=bottom_type,
    )
    if params.involute_step <= 0:
        raise ValueError("involute_step must be positive")
    if params.spline_division_num <= 0:
        raise ValueError("spline_division_num must be positive")

    alpha = math.radians(params.pressure_angle_deg)
    pitch = params.module * math.pi
    tooth_thickness = pitch / 2.0

    diameter_pitch = params.teeth * params.module
    diameter_addendum = diameter_pitch + 2 * params.module
    diameter_dedendum = diameter_pitch - 2.5 * params.module
    diameter_base = diameter_pitch * math.cos(alpha)

    radius_pitch = diameter_pitch / 2.0
    radius_addendum = diameter_addendum / 2.0
    radius_dedendum = diameter_dedendum / 2.0
    radius_base = diameter_base / 2.0

    angle_per_tooth = 2 * math.pi / params.teeth
    angle_thickness = tooth_thickness / radius_pitch
    inv_at_pitch = _inv(math.acos(radius_base / radius_pitch))
    angle_base = angle_thickness + inv_at_pitch * 2.0
    angle_bottom = angle_per_tooth - angle_base

    cos_bottom = math.cos(-angle_bottom)
    sin_bottom = math.sin(-angle_bottom)

    inv_segments = max(2, int(math.ceil((radius_addendum - radius_base) / params.involute_step)))

    profile: PointList = []

    for tooth_idx in range(params.teeth):
        t = angle_per_tooth * tooth_idx
        cos_t = math.cos(t)
        sin_t = math.sin(t)

        xa = radius_base * cos_t
        ya = radius_base * sin_t
        xb = radius_dedendum * cos_t
        yb = radius_dedendum * sin_t
        xc = xb * cos_bottom - yb * sin_bottom
        yc = xb * sin_bottom + yb * cos_bottom
        xd = xa * cos_bottom - ya * sin_bottom
        yd = xa * sin_bottom + ya * cos_bottom
        base_points = [(xd, yd), (xc, yc), (xb, yb), (xa, ya)]

        if params.bottom_type == "line":
            _add_bottom_points_line(profile, base_points)
        else:
            _add_bottom_points_spline(
                profile, base_points, division_num=params.spline_division_num
            )

        points_inv1: PointList = []
        points_inv2: PointList = []
        cos_inv2 = math.cos(t + angle_base)
        sin_inv2 = math.sin(t + angle_base)

        for segment in range(inv_segments + 1):
            r = radius_base + (radius_addendum - radius_base) * (segment / inv_segments)
            r = max(r, radius_base)
            inv_alpha = _inv(math.acos(radius_base / r))
            x = r * math.cos(inv_alpha)
            y = r * math.sin(inv_alpha)

            x1 = x * cos_t - y * sin_t
            y1 = x * sin_t + y * cos_t
            points_inv1.append((x1, y1))

            x2 = x * cos_inv2 - (-y) * sin_inv2
            y2 = x * sin_inv2 + (-y) * cos_inv2
            points_inv2.append((x2, y2))

        profile.extend(points_inv1)
        profile.extend(reversed(points_inv2))

    _ensure_closed(profile)

    blueprints = {
        "diameter_addendum": diameter_addendum,
        "diameter_pitch": diameter_pitch,
        "diameter_base": diameter_base,
        "diameter_dedendum": diameter_dedendum,
        "radius_addendum": radius_addendum,
        "radius_pitch": radius_pitch,
        "radius_base": radius_base,
        "radius_dedendum": radius_dedendum,
    }

    return profile, blueprints


__all__ = ["make_gear_figure"]
