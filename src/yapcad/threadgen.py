"""Thread profile generator for yapCAD."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable, List

from yapcad.geom import point, epsilon

__all__ = [
    "ThreadProfile",
    "metric_profile",
    "unified_profile",
    "sample_thread_profile",
]


@dataclass(frozen=True)
class ThreadProfile:
    D_nominal: float
    P_pitch: float
    thread_angle: float = 60.0
    crest_flat_ratio: float = 1.0 / 8.0
    root_flat_ratio: float = 1.0 / 4.0
    thread_depth_ratio: float = 0.54127
    handedness: str = "right"
    starts: int = 1
    internal: bool = False
    taper_ratio: float = 0.0

    def lead(self) -> float:
        return self.P_pitch * max(1, self.starts)


def metric_profile(d_nominal_mm: float, pitch_mm: float, *, internal: bool = False) -> ThreadProfile:
    crest = 1.0 / 8.0
    root = 1.0 / 4.0 if not internal else 1.0 / 8.0
    return ThreadProfile(
        D_nominal=d_nominal_mm,
        P_pitch=pitch_mm,
        crest_flat_ratio=crest,
        root_flat_ratio=root,
        thread_depth_ratio=0.54127,
        internal=internal,
    )


def unified_profile(d_nominal_in: float, tpi: float, *, internal: bool = False) -> ThreadProfile:
    pitch = 25.4 / tpi
    crest = 1.0 / 8.0
    root = 1.0 / 4.0 if not internal else 1.0 / 8.0
    depth_ratio = 0.64952 / tpi / pitch
    return ThreadProfile(
        D_nominal=d_nominal_in * 25.4,
        P_pitch=pitch,
        crest_flat_ratio=crest,
        root_flat_ratio=root,
        thread_depth_ratio=depth_ratio,
        internal=internal,
    )


def sample_thread_profile(
    profile: ThreadProfile,
    x_start: float,
    x_end: float,
    theta_deg: float,
    *,
    samples_per_pitch: int = 1,
) -> tuple[List[List[float]], int]:
    if x_end <= x_start:
        raise ValueError("x_end must be greater than x_start")
    if profile.P_pitch <= 0:
        raise ValueError("pitch must be positive")

    points, wrap = _collect_breakpoints(profile, x_start, x_end, samples_per_pitch)
    lead = profile.lead()
    offset = (theta_deg / 360.0) * lead
    if profile.handedness.lower() == "left":
        offset = -offset

    sampled = []
    z_min = x_start
    z_max = x_end
    for base_z in points:
        shifted = base_z + offset
        actual_z = shifted
        if actual_z < z_min:
            actual_z = z_min
        elif actual_z > z_max:
            actual_z = z_max
        radius = _radius_at(profile, actual_z, base_z)
        sampled.append(point(actual_z, radius))
    return sampled, wrap


def _collect_breakpoints(profile: ThreadProfile, x_start: float, x_end: float, samples_per_pitch: int) -> tuple[List[float], int]:
    pitch = profile.P_pitch
    breakpoints = _pitch_breakpoints(profile)
    extra = max(1, samples_per_pitch)
    base_start = math.floor(x_start / pitch) - 1
    base_end = math.ceil(x_end / pitch) + 1

    xs = {x_start, x_end}
    pitch_samples = set()
    for idx in range(base_start, base_end + 1):
        base = idx * pitch
        for bp in breakpoints:
            x = base + bp
            if x_start - 1e-9 <= x <= x_end + 1e-9:
                xs.add(min(max(x, x_start), x_end))
        for seg in range(1, extra):
            frac = seg / extra
            x = base + frac * pitch
            if x_start - 1e-9 <= x <= x_end + 1e-9:
                xs.add(min(max(x, x_start), x_end))
    # samples within single pitch starting at zero
    for base in (0,):
        for bp in breakpoints:
            pitch_samples.add(base + bp)
        for seg in range(1, extra):
            frac = seg / extra
            pitch_samples.add(base + frac * pitch)
        pitch_samples.add(pitch)

    wrap = max(1, len(sorted(pitch_samples)) - 1)

    return sorted(xs), wrap


def _pitch_breakpoints(profile: ThreadProfile) -> Iterable[float]:
    pitch = profile.P_pitch
    crest_half = max(0.0, profile.crest_flat_ratio * pitch / 2.0)
    root_half = max(0.0, profile.root_flat_ratio * pitch / 2.0)
    flank = max(0.0, (pitch - 2 * crest_half - 2 * root_half) / 2.0)
    return [
        0.0,
        crest_half,
        crest_half + flank,
        crest_half + flank + 2 * root_half,
        pitch - crest_half,
        pitch,
    ]


def _radius_at(profile: ThreadProfile, x_actual: float, x_effective: float) -> float:
    pitch = profile.P_pitch
    local = x_effective - math.floor(x_effective / pitch) * pitch
    base_major = (profile.D_nominal + profile.taper_ratio * x_actual) / 2.0
    depth = profile.thread_depth_ratio * pitch

    if profile.internal:
        root = base_major
        crest = base_major + depth
    else:
        crest = base_major
        root = max(crest - depth, 0.0)

    crest_half = max(0.0, profile.crest_flat_ratio * pitch / 2.0)
    root_half = max(0.0, profile.root_flat_ratio * pitch / 2.0)
    flank = max(0.0, (pitch - 2 * crest_half - 2 * root_half) / 2.0)

    x1 = crest_half
    x2 = crest_half + flank
    x3 = x2 + 2 * root_half
    x4 = pitch - crest_half

    if local <= x1:
        return crest
    if local <= x2:
        t = (local - x1) / max(flank, epsilon)
        return crest + t * (root - crest)
    if local <= x3:
        return root
    if local <= x4:
        t = (local - x3) / max(flank, epsilon)
        return root + t * (crest - root)
    return crest
