"""Fastener construction using catalog data.

This module builds fasteners using dimensions from YAML catalogs,
delegating to the core thread generation in yapcad.threadgen.
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import Optional

from yapcad.geom import epsilon
from yapcad.threadgen import ThreadProfile, metric_profile, unified_profile

from .catalog import get_bolt_data, get_nut_data

__all__ = [
    "build_hex_bolt_from_catalog",
    "build_hex_nut_from_catalog",
]


def _make_profile_from_catalog(
    thread_data: dict,
    internal: bool = False,
    handedness: str = "right",
    starts: int = 1,
) -> ThreadProfile:
    """Create a ThreadProfile from catalog thread data.

    Uses the tolerances from the catalog to adjust thread geometry.
    """
    nominal_d = thread_data["nominal_diameter"]
    pitch = thread_data["pitch"]
    tolerances = thread_data.get("tolerances", {})

    # Determine if metric or unified based on presence of 'tpi'
    # (Unified catalog stores tpi, metric stores pitch directly)
    # Both catalogs store pitch in mm for consistency

    # Get tolerance adjustments
    if internal:
        # Internal thread: apply EI (lower deviation, usually 0 for H position)
        ei = tolerances.get("EI", 0.0)
        # For internal threads, the minor diameter limit is what matters
        # TD1 is the tolerance on minor diameter
        td1 = tolerances.get("TD1", 0.0)

        # Adjust nominal diameter for internal thread
        # Basic internal thread has D (major) = nominal
        # The actual cutting diameter accounts for allowance
        adjusted_d = nominal_d + ei

        profile = ThreadProfile(
            D_nominal=adjusted_d,
            P_pitch=pitch,
            thread_angle=60.0,
            crest_flat_ratio=1.0 / 8.0,
            root_flat_ratio=1.0 / 8.0,  # Internal threads have smaller root flat
            thread_depth_ratio=0.57,
            handedness=handedness,
            starts=starts,
            internal=True,
        )
    else:
        # External thread: apply es (upper deviation, negative for g position)
        es = tolerances.get("es", 0.0)
        # Td2 is pitch diameter tolerance, Td is major diameter tolerance
        td2 = tolerances.get("Td2", 0.0)
        td = tolerances.get("Td", 0.0)

        # Adjust nominal diameter for external thread
        # The es deviation is typically negative (undersized)
        adjusted_d = nominal_d + es

        profile = ThreadProfile(
            D_nominal=adjusted_d,
            P_pitch=pitch,
            thread_angle=60.0,
            crest_flat_ratio=1.0 / 8.0,
            root_flat_ratio=1.0 / 4.0,  # External threads have larger root flat
            thread_depth_ratio=0.57,
            handedness=handedness,
            starts=starts,
            internal=False,
        )

    return profile


def _default_thread_length(shank_length: float, diameter: float) -> float:
    """Calculate default thread length per ISO 4014 / ASME B18.2.1.

    For bolts up to 125mm: thread_length = 2*d + 6mm
    For bolts 125-200mm: thread_length = 2*d + 12mm
    For bolts over 200mm: thread_length = 2*d + 25mm

    But never more than shank_length.
    """
    if shank_length <= 125:
        tl = 2 * diameter + 6
    elif shank_length <= 200:
        tl = 2 * diameter + 12
    else:
        tl = 2 * diameter + 25

    # Thread length can't exceed shank length
    return min(tl, shank_length)


def build_hex_bolt_from_catalog(
    thread_series: str,
    size: str,
    length: float,
    tolerance_class: str = "6g",
    thread_length: Optional[float] = None,
    starts: int = 1,
    thread_arc_samples: int = 180,
    thread_samples_per_pitch: int = 6,
    catalog_path: Optional[Path] = None,
):
    """Build a hex bolt using catalog dimensions.

    Args:
        thread_series: "metric_coarse", "unified_coarse", etc.
        size: Size designation from catalog
        length: Shank length in mm
        tolerance_class: Thread tolerance class
        thread_length: Thread length in mm (or computed if None)
        starts: Number of thread starts
        thread_arc_samples: Angular samples for thread surface
        thread_samples_per_pitch: Profile samples per pitch
        catalog_path: Optional custom catalog path

    Returns:
        yapCAD solid representing the bolt
    """
    # Import here to avoid circular dependency
    from yapcad.fasteners_legacy import (
        HexCapScrewSpec,
        build_hex_cap_screw,
    )

    # Get catalog data
    bolt_data = get_bolt_data(thread_series, size, tolerance_class, catalog_path)
    thread_data = bolt_data["thread"]
    head_data = bolt_data["head"]

    # Create thread profile from catalog data
    profile = _make_profile_from_catalog(
        thread_data,
        internal=False,
        handedness="right",
        starts=starts,
    )

    # Calculate thread length
    nominal_d = thread_data["nominal_diameter"]
    if thread_length is None:
        thread_length = _default_thread_length(length, nominal_d)
    thread_length = min(max(thread_length, epsilon), length)

    # Build spec for legacy builder
    spec = HexCapScrewSpec(
        diameter=nominal_d,
        thread_length=thread_length,
        shank_length=length,
        head_height=head_data["head_height"],
        head_flat_diameter=head_data["width_across_flats"],
        washer_thickness=head_data.get("washer_face_thickness", 0.5),
        washer_diameter=head_data.get("washer_face_diameter"),
        shank_diameter=nominal_d,  # Could be slightly larger in reality
        starts=starts,
        thread_arc_samples=thread_arc_samples,
        thread_samples_per_pitch=thread_samples_per_pitch,
    )

    return build_hex_cap_screw(profile, spec)


def build_hex_nut_from_catalog(
    thread_series: str,
    size: str,
    tolerance_class: str = "6H",
    handedness: str = "right",
    starts: int = 1,
    thread_arc_samples: int = 180,
    thread_samples_per_pitch: int = 6,
    catalog_path: Optional[Path] = None,
):
    """Build a hex nut using catalog dimensions.

    Args:
        thread_series: "metric_coarse", "unified_coarse", etc.
        size: Size designation from catalog
        tolerance_class: Thread tolerance class
        handedness: "right" or "left"
        starts: Number of thread starts
        thread_arc_samples: Angular samples for thread surface
        thread_samples_per_pitch: Profile samples per pitch
        catalog_path: Optional custom catalog path

    Returns:
        yapCAD solid representing the nut
    """
    # Import here to avoid circular dependency
    from yapcad.fasteners_legacy import (
        HexNutSpec,
        build_hex_nut,
    )

    # Get catalog data
    nut_data = get_nut_data(thread_series, size, tolerance_class, catalog_path)
    thread_data = nut_data["thread"]
    body_data = nut_data["body"]

    # Create thread profile from catalog data
    profile = _make_profile_from_catalog(
        thread_data,
        internal=True,
        handedness=handedness,
        starts=starts,
    )

    # Build spec for legacy builder
    spec = HexNutSpec(
        diameter=thread_data["nominal_diameter"],
        pitch=thread_data["pitch"],
        width_flat=body_data["width_across_flats"],
        thickness=body_data["thickness"],
        handedness=handedness,
        starts=starts,
        thread_arc_samples=thread_arc_samples,
        thread_samples_per_pitch=thread_samples_per_pitch,
    )

    return build_hex_nut(profile, spec)
