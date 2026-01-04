"""Fastener generation with catalog-based dimensions.

This package provides parametric fastener generation (hex bolts, hex nuts)
with dimensions from YAML catalog files. Bundled catalogs cover common
metric (ISO) and unified (ASME) sizes; users can extend with custom catalogs.

Quick Start:
    >>> from yapcad.fasteners import metric_hex_bolt, metric_hex_nut
    >>> bolt = metric_hex_bolt("M8", 25.0)
    >>> nut = metric_hex_nut("M8")

Available Functions:
    Metric (ISO):
        metric_hex_bolt(size, length, tolerance_class="6g")
        metric_hex_nut(size, tolerance_class="6H")

    Unified (ASME):
        unified_hex_bolt(size, length, tolerance_class="2A")
        unified_hex_nut(size, tolerance_class="2B")

Catalog Customization:
    Set YAPCAD_FASTENER_DATA environment variable to add custom catalog
    directories. These are searched before bundled data.

    Example:
        export YAPCAD_FASTENER_DATA="/path/to/my/fasteners"

    See docs/fastener_catalog_schema.md for YAML format specification.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from .catalog import (
    YAPCAD_FASTENER_DATA,
    THREAD_SERIES,
    load_catalog,
    list_available_sizes,
    list_tolerance_classes,
    get_thread_data,
    get_bolt_data,
    get_nut_data,
    clear_cache,
)

# Re-export legacy API for backward compatibility
from yapcad.fasteners_legacy import (
    HexCapScrewSpec,
    HexNutSpec,
    build_hex_cap_screw,
    build_hex_nut,
    metric_hex_cap_catalog,
    metric_hex_cap_screw,
    metric_hex_nut_catalog,
    unified_hex_cap_catalog,
    unified_hex_cap_screw,
    unified_hex_nut_catalog,
    # Note: metric_hex_nut and unified_hex_nut are redefined below
    # with improved catalog-based API
)

__all__ = [
    # Environment variable
    "YAPCAD_FASTENER_DATA",
    "THREAD_SERIES",
    # Catalog functions
    "load_catalog",
    "list_available_sizes",
    "list_tolerance_classes",
    "get_thread_data",
    "get_bolt_data",
    "get_nut_data",
    "clear_cache",
    # New catalog-based constructors
    "metric_hex_bolt",
    "metric_hex_nut",
    "unified_hex_bolt",
    "unified_hex_nut",
    # Legacy API (backward compatibility)
    "HexCapScrewSpec",
    "HexNutSpec",
    "build_hex_cap_screw",
    "build_hex_nut",
    "metric_hex_cap_catalog",
    "metric_hex_cap_screw",
    "metric_hex_nut_catalog",
    "unified_hex_cap_catalog",
    "unified_hex_cap_screw",
    "unified_hex_nut_catalog",
]


def metric_hex_bolt(
    size: str,
    length: float,
    *,
    tolerance_class: str = "6g",
    thread_length: Optional[float] = None,
    starts: int = 1,
    thread_arc_samples: int = 180,
    thread_samples_per_pitch: int = 6,
    catalog_path: Optional[Path] = None,
):
    """Create a metric hex bolt per ISO 4014/4017.

    Args:
        size: Thread size designation (e.g., "M8", "M10", "M12")
        length: Total shank length in mm
        tolerance_class: Thread tolerance class (default "6g")
            Common classes: "6g" (general), "4g6g" (close), "8g" (coarse)
        thread_length: Override automatic thread length calculation (mm)
        starts: Number of thread starts (default 1)
        thread_arc_samples: Angular resolution for thread generation
        thread_samples_per_pitch: Samples per pitch for thread profile
        catalog_path: Optional path to custom catalog YAML file

    Returns:
        yapCAD solid representing the hex bolt

    Raises:
        KeyError: If size or tolerance_class not found in catalog

    Examples:
        >>> bolt = metric_hex_bolt("M8", 25.0)
        >>> bolt = metric_hex_bolt("M10", 40.0, tolerance_class="4g6g")
        >>> bolt = metric_hex_bolt("M6", 20.0, thread_length=15.0)
    """
    from yapcad.fasteners.builders import build_hex_bolt_from_catalog

    return build_hex_bolt_from_catalog(
        thread_series="metric_coarse",
        size=size,
        length=length,
        tolerance_class=tolerance_class,
        thread_length=thread_length,
        starts=starts,
        thread_arc_samples=thread_arc_samples,
        thread_samples_per_pitch=thread_samples_per_pitch,
        catalog_path=catalog_path,
    )


def metric_hex_nut(
    size: str,
    *,
    tolerance_class: str = "6H",
    handedness: str = "right",
    starts: int = 1,
    thread_arc_samples: int = 180,
    thread_samples_per_pitch: int = 6,
    catalog_path: Optional[Path] = None,
):
    """Create a metric hex nut per ISO 4032.

    Args:
        size: Thread size designation (e.g., "M8", "M10")
        tolerance_class: Thread tolerance class (default "6H" for internal)
            Common classes: "6H" (general), "5H" (close), "7H" (coarse)
        handedness: "right" or "left" hand thread
        starts: Number of thread starts (default 1)
        thread_arc_samples: Angular resolution for thread generation
        thread_samples_per_pitch: Samples per pitch for thread profile
        catalog_path: Optional path to custom catalog YAML file

    Returns:
        yapCAD solid representing the hex nut

    Raises:
        KeyError: If size or tolerance_class not found in catalog

    Examples:
        >>> nut = metric_hex_nut("M8")
        >>> nut = metric_hex_nut("M10", tolerance_class="5H")
    """
    from yapcad.fasteners.builders import build_hex_nut_from_catalog

    return build_hex_nut_from_catalog(
        thread_series="metric_coarse",
        size=size,
        tolerance_class=tolerance_class,
        handedness=handedness,
        starts=starts,
        thread_arc_samples=thread_arc_samples,
        thread_samples_per_pitch=thread_samples_per_pitch,
        catalog_path=catalog_path,
    )


def unified_hex_bolt(
    size: str,
    length: float,
    *,
    tolerance_class: str = "2A",
    thread_length: Optional[float] = None,
    starts: int = 1,
    thread_arc_samples: int = 180,
    thread_samples_per_pitch: int = 6,
    catalog_path: Optional[Path] = None,
):
    """Create a unified (UNC/UNF) hex bolt per ASME B18.2.1.

    Args:
        size: Thread size designation (e.g., "1/4-20", "#10-24", "1/2-13")
        length: Total shank length in inches (will be converted to mm internally)
        tolerance_class: Thread class (default "2A")
            Classes: "1A" (loose), "2A" (general), "3A" (close)
        thread_length: Override automatic thread length in inches
        starts: Number of thread starts (default 1)
        thread_arc_samples: Angular resolution for thread generation
        thread_samples_per_pitch: Samples per pitch for thread profile
        catalog_path: Optional path to custom catalog YAML file

    Returns:
        yapCAD solid representing the hex bolt

    Raises:
        KeyError: If size or tolerance_class not found in catalog

    Examples:
        >>> bolt = unified_hex_bolt("1/4-20", 1.0)  # 1 inch long
        >>> bolt = unified_hex_bolt("1/2-13", 2.0, tolerance_class="3A")
    """
    from yapcad.fasteners.builders import build_hex_bolt_from_catalog

    # Convert inches to mm
    length_mm = length * 25.4
    thread_length_mm = thread_length * 25.4 if thread_length is not None else None

    return build_hex_bolt_from_catalog(
        thread_series="unified_coarse",
        size=size,
        length=length_mm,
        tolerance_class=tolerance_class,
        thread_length=thread_length_mm,
        starts=starts,
        thread_arc_samples=thread_arc_samples,
        thread_samples_per_pitch=thread_samples_per_pitch,
        catalog_path=catalog_path,
    )


def unified_hex_nut(
    size: str,
    *,
    tolerance_class: str = "2B",
    handedness: str = "right",
    starts: int = 1,
    thread_arc_samples: int = 180,
    thread_samples_per_pitch: int = 6,
    catalog_path: Optional[Path] = None,
):
    """Create a unified (UNC/UNF) hex nut per ASME B18.2.2.

    Args:
        size: Thread size designation (e.g., "1/4-20", "#10-24", "1/2-13")
        tolerance_class: Thread class (default "2B" for internal)
            Classes: "1B" (loose), "2B" (general), "3B" (close)
        handedness: "right" or "left" hand thread
        starts: Number of thread starts (default 1)
        thread_arc_samples: Angular resolution for thread generation
        thread_samples_per_pitch: Samples per pitch for thread profile
        catalog_path: Optional path to custom catalog YAML file

    Returns:
        yapCAD solid representing the hex nut

    Raises:
        KeyError: If size or tolerance_class not found in catalog

    Examples:
        >>> nut = unified_hex_nut("1/4-20")
        >>> nut = unified_hex_nut("1/2-13", tolerance_class="3B")
    """
    from yapcad.fasteners.builders import build_hex_nut_from_catalog

    return build_hex_nut_from_catalog(
        thread_series="unified_coarse",
        size=size,
        tolerance_class=tolerance_class,
        handedness=handedness,
        starts=starts,
        thread_arc_samples=thread_arc_samples,
        thread_samples_per_pitch=thread_samples_per_pitch,
        catalog_path=catalog_path,
    )
