"""Parametric fastener helpers (threads + hex-cap screws)."""

from __future__ import annotations

import math
from copy import deepcopy
from dataclasses import asdict, dataclass, replace
from typing import Dict, Iterable

from yapcad.geom import arc, epsilon, point
from yapcad.geom_util import geomlist2poly
from yapcad.geom3d import poly2surfaceXY, solid, solid_boolean, translate
from yapcad.geom3d_util import (
    circleSurface,
    conic,
    extrude,
    makeRevolutionThetaSamplingSurface,
)
from yapcad.metadata import add_tags, get_solid_metadata, set_layer
from yapcad.threadgen import ThreadProfile, metric_profile, sample_thread_profile, unified_profile

__all__ = [
    "HexCapScrewSpec",
    "build_hex_cap_screw",
    "metric_hex_cap_screw",
    "unified_hex_cap_screw",
    "metric_hex_cap_catalog",
    "unified_hex_cap_catalog",
]


@dataclass(frozen=True)
class HexCapScrewSpec:
    """Dimensions (all millimeters) for :func:`build_hex_cap_screw`."""

    diameter: float
    thread_length: float
    shank_length: float
    head_height: float
    head_flat_diameter: float
    washer_thickness: float = 0.5
    washer_diameter: float | None = None
    shank_diameter: float | None = None
    thread_arc_samples: int = 180
    thread_samples_per_pitch: int = 6


def build_hex_cap_screw(profile: ThreadProfile, spec: HexCapScrewSpec):
    """Create a watertight solid representing a hex cap screw.

    Args:
        profile: External thread profile describing nominal diameter/pitch.
        spec: Geometry parameters (all millimeters).
    """

    _validate_spec(spec)
    thread = _build_thread(profile, spec)
    components = [thread]
    unthreaded = max(spec.shank_length - spec.thread_length, 0.0)
    if unthreaded > epsilon:
        components.append(_build_shank(spec, unthreaded))

    if spec.washer_thickness > epsilon:
        components.append(_build_washer(spec))

    components.append(_build_hex_head(spec))
    try:
        body = _union_solids(components)
    except Exception:
        body = _stack_solids(components)
    _scrub_surface_octrees(body)
    meta = get_solid_metadata(body, create=True)
    add_tags(meta, ["fastener", "hex_cap_screw"])
    set_layer(meta, "hardware")
    meta["hex_cap_screw"] = {
        "diameter": spec.diameter,
        "thread_length": spec.thread_length,
        "shank_length": spec.shank_length,
        "head_height": spec.head_height,
        "head_flat": spec.head_flat_diameter,
        "washer_thickness": spec.washer_thickness,
        "washer_diameter": spec.washer_diameter or spec.head_flat_diameter,
    }
    return body


def metric_hex_cap_screw(
    size: str,
    length: float,
    *,
    thread_length: float | None = None,
    starts: int = 1,
    thread_arc_samples: int = 180,
    thread_samples_per_pitch: int = 6,
):
    """Return a hex screw for a metric size (e.g. ``'M8'``)."""

    dims = _METRIC_TABLE[size.upper()]
    tl = thread_length or _default_thread_length(length, dims.diameter)
    tl = min(max(tl, epsilon), length)
    profile = metric_profile(dims.diameter, dims.pitch)
    if starts != profile.starts:
        profile = replace(profile, starts=starts)
    washer_diameter, washer_thickness = _normalized_washer_dims(
        head_flat=dims.head_flat,
        base_thickness=dims.washer_thickness,
    )
    spec = HexCapScrewSpec(
        diameter=dims.diameter,
        thread_length=tl,
        shank_length=length,
        head_height=dims.head_height,
        head_flat_diameter=dims.head_flat,
        washer_thickness=washer_thickness,
        washer_diameter=washer_diameter,
        thread_arc_samples=thread_arc_samples,
        thread_samples_per_pitch=thread_samples_per_pitch,
    )
    return build_hex_cap_screw(profile, spec)


def unified_hex_cap_screw(
    size: str,
    length_in: float,
    *,
    thread_length_in: float | None = None,
    starts: int = 1,
    thread_arc_samples: int = 180,
    thread_samples_per_pitch: int = 6,
):
    """Return a hex screw for a UNC/UNF imperial size (e.g. ``'1/4-20'``)."""

    dims = _UNIFIED_TABLE[size.lower()]
    length = length_in * 25.4
    tl = (thread_length_in * 25.4) if thread_length_in is not None else _default_thread_length(length, dims.diameter)
    tl = min(max(tl, epsilon), length)
    profile = unified_profile(dims.diameter / 25.4, dims.tpi)
    if starts != profile.starts:
        profile = replace(profile, starts=starts)
    washer_diameter, washer_thickness = _normalized_washer_dims(
        head_flat=dims.head_flat,
        base_thickness=dims.washer_thickness,
    )
    spec = HexCapScrewSpec(
        diameter=dims.diameter,
        thread_length=tl,
        shank_length=length,
        head_height=dims.head_height,
        head_flat_diameter=dims.head_flat,
        washer_thickness=washer_thickness,
        washer_diameter=washer_diameter,
        thread_arc_samples=thread_arc_samples,
        thread_samples_per_pitch=thread_samples_per_pitch,
    )
    return build_hex_cap_screw(profile, spec)


def _build_thread(profile: ThreadProfile, spec: HexCapScrewSpec):
    samples0, _ = sample_thread_profile(
        profile,
        0.0,
        spec.thread_length,
        0.0,
        samples_per_pitch=spec.thread_samples_per_pitch,
    )
    bottom_r = samples0[0][1]
    top_r = samples0[-1][1]

    def contour(_z0: float, _z1: float, theta: float):
        return sample_thread_profile(
            profile,
            0.0,
            spec.thread_length,
            theta,
            samples_per_pitch=spec.thread_samples_per_pitch,
        )

    surface = makeRevolutionThetaSamplingSurface(
        contour,
        0.0,
        spec.thread_length,
        arcSamples=max(24, spec.thread_arc_samples),
        endcaps=False,
    )
    bottom_cap = circleSurface(point(0, 0, 0), bottom_r, zup=False)
    top_cap = circleSurface(point(0, 0, spec.thread_length), top_r, zup=True)
    return solid([surface, bottom_cap, top_cap])


def _build_shank(spec: HexCapScrewSpec, length: float):
    radius = (spec.shank_diameter or spec.diameter) / 2.0
    return conic(
        radius,
        radius,
        length,
        center=point(0, 0, spec.thread_length),
    )


def _build_washer(spec: HexCapScrewSpec):
    radius = (spec.washer_diameter or spec.head_flat_diameter) / 2.0
    return conic(
        radius,
        radius,
        spec.washer_thickness,
        center=point(0, 0, spec.shank_length),
    )


def _build_hex_head(spec: HexCapScrewSpec):
    base_z = spec.shank_length + spec.washer_thickness
    incircle = spec.head_flat_diameter / 2.0
    circum = incircle / math.cos(math.pi / 6.0)
    pts = []
    for idx in range(6):
        angle = math.pi / 6.0 + idx * (math.pi / 3.0)
        pts.append(point(circum * math.cos(angle), circum * math.sin(angle), 0.0))
    pts.append(pts[0])
    surface, _ = poly2surfaceXY(pts)
    head = extrude(surface, spec.head_height)
    return translate(head, point(0, 0, base_z))


def _default_thread_length(shank_length: float, diameter: float) -> float:
    return max(shank_length - diameter * 0.75, shank_length * 0.65)


def _validate_spec(spec: HexCapScrewSpec):
    if spec.diameter <= epsilon:
        raise ValueError("diameter must be positive")
    if spec.thread_length <= epsilon:
        raise ValueError("thread_length must be positive")
    if spec.thread_length - spec.shank_length > epsilon:
        raise ValueError("thread_length must be <= shank_length")
    if spec.head_height <= epsilon:
        raise ValueError("head_height must be positive")
    if spec.head_flat_diameter <= spec.diameter:
        raise ValueError("head_flat_diameter must exceed the shank diameter")
    if spec.thread_arc_samples < 12:
        raise ValueError("thread_arc_samples must be >= 12")
    if spec.thread_samples_per_pitch < 1:
        raise ValueError("thread_samples_per_pitch must be >= 1")


def _union_solids(solids: Iterable[list]):
    iterator = iter(solids)
    try:
        result = next(iterator)
    except StopIteration:
        raise ValueError("no solids supplied")
    for part in iterator:
        result = solid_boolean(result, part, "union")
    return result


def _stack_solids(solids: Iterable[list]):
    collected = []
    for solid_part in solids:
        if not solid_part[1]:
            continue
        collected.extend([deepcopy(surf) for surf in solid_part[1]])
    if not collected:
        raise ValueError("no surfaces found while stacking solids")
    return solid(collected)


def _scrub_surface_octrees(sld):
    for surf in sld[1]:
        if len(surf) >= 7 and isinstance(surf[6], dict):
            surf[6].pop('_octree', None)
            surf[6].pop('_octree_dirty', None)


@dataclass(frozen=True)
class _MetricDims:
    diameter: float
    pitch: float
    head_height: float
    head_flat: float
    washer_thickness: float
    washer_diameter: float | None = None


@dataclass(frozen=True)
class _UnifiedDims:
    diameter: float  # millimeters
    tpi: float
    head_height: float
    head_flat: float
    washer_thickness: float
    washer_diameter: float | None = None


_METRIC_TABLE: Dict[str, _MetricDims] = {
    "M6": _MetricDims(6.0, 1.0, 4.0, 10.0, 0.4),
    "M8": _MetricDims(8.0, 1.25, 5.3, 13.0, 0.5),
    "M10": _MetricDims(10.0, 1.5, 6.4, 16.0, 0.6),
}

_UNIFIED_TABLE: Dict[str, _UnifiedDims] = {
    "1/4-20": _UnifiedDims(6.35, 20, 4.14, 11.11, 0.5),
    "5/16-18": _UnifiedDims(7.9375, 18, 5.23, 12.7, 0.55),
    "3/8-16": _UnifiedDims(9.525, 16, 6.22, 14.2875, 0.65),
}


def metric_hex_cap_catalog():
    return {name: asdict(dim) for name, dim in _METRIC_TABLE.items()}


def unified_hex_cap_catalog():
    return {name: asdict(dim) for name, dim in _UNIFIED_TABLE.items()}


def _normalized_washer_dims(head_flat: float, base_thickness: float | None):
    diameter = 0.95 * head_flat
    thickness = base_thickness if base_thickness is not None else head_flat * 0.05
    thickness = max(thickness / 2.0, epsilon * 10.0)
    return diameter, thickness
