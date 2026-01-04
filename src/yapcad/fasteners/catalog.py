"""Fastener catalog loading with bundled data and external override support.

This module provides a flexible catalog system for fastener dimensions:
- Bundled YAML data files for common metric and unified sizes
- Environment variable override for custom data directories
- User config directory support (~/.config/yapcad/fasteners/)
- Explicit path override in API calls

Environment Variables:
    YAPCAD_FASTENER_DATA: Colon-separated (or semicolon on Windows) paths
                          to directories containing custom YAML catalog files.
                          These are searched before bundled data.

Example:
    export YAPCAD_FASTENER_DATA="/path/to/my/fasteners:/another/path"
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional
from functools import lru_cache

import yaml

__all__ = [
    "YAPCAD_FASTENER_DATA",
    "load_catalog",
    "list_available_sizes",
    "list_tolerance_classes",
    "get_thread_data",
    "get_bolt_data",
    "get_nut_data",
    "clear_cache",
]

# Environment variable name for custom data paths
YAPCAD_FASTENER_DATA = "YAPCAD_FASTENER_DATA"

# Bundled data location (relative to this file)
_BUNDLED_DATA_DIR = Path(__file__).parent / "data"

# Supported thread series
THREAD_SERIES = {
    "metric_coarse": "ISO metric coarse thread (M series)",
    "metric_fine": "ISO metric fine thread",
    "unified_coarse": "Unified National Coarse (UNC)",
    "unified_fine": "Unified National Fine (UNF)",
}


def clear_cache() -> None:
    """Clear all cached catalog data.

    Call this if you modify external catalog files and want to reload.
    """
    _get_data_dirs.cache_clear()
    _load_catalog_cached.cache_clear()


@lru_cache(maxsize=None)
def _get_data_dirs() -> tuple[Path, ...]:
    """Return tuple of data directories to search, in priority order.

    Search order:
        1. Directories from YAPCAD_FASTENER_DATA environment variable
        2. User config directory (~/.config/yapcad/fasteners/)
        3. Bundled data directory
    """
    dirs: List[Path] = []

    # 1. Environment variable (highest priority)
    env_path = os.environ.get(YAPCAD_FASTENER_DATA)
    if env_path:
        # Use appropriate path separator for platform
        sep = ";" if sys.platform == "win32" else ":"
        for p in env_path.split(sep):
            p = p.strip()
            if p:
                path = Path(p).expanduser().resolve()
                if path.is_dir():
                    dirs.append(path)

    # 2. User config directory
    if sys.platform == "win32":
        config_base = Path(os.environ.get("APPDATA", "~")).expanduser()
    else:
        config_base = Path.home() / ".config"

    user_config = config_base / "yapcad" / "fasteners"
    if user_config.is_dir():
        dirs.append(user_config)

    # 3. Bundled data (always available as fallback)
    if _BUNDLED_DATA_DIR.is_dir():
        dirs.append(_BUNDLED_DATA_DIR)

    return tuple(dirs)


@lru_cache(maxsize=32)
def _load_catalog_cached(thread_series: str, custom_path_str: Optional[str]) -> Dict[str, Any]:
    """Cached catalog loading (string path for hashability)."""
    custom_path = Path(custom_path_str) if custom_path_str else None
    return _load_catalog_impl(thread_series, custom_path)


def _load_catalog_impl(thread_series: str, custom_path: Optional[Path]) -> Dict[str, Any]:
    """Implementation of catalog loading."""
    filename = f"{thread_series}.yaml"

    if custom_path:
        if not custom_path.exists():
            raise FileNotFoundError(f"Custom catalog not found: {custom_path}")
        return _load_yaml(custom_path)

    for data_dir in _get_data_dirs():
        path = data_dir / filename
        if path.exists():
            return _load_yaml(path)

    searched = [str(d) for d in _get_data_dirs()]
    raise FileNotFoundError(
        f"No catalog found for '{thread_series}'.\n"
        f"Searched directories: {searched}\n"
        f"Available series: {list(THREAD_SERIES.keys())}"
    )


def _load_yaml(path: Path) -> Dict[str, Any]:
    """Load and validate a YAML catalog file."""
    with open(path, 'r', encoding='utf-8') as f:
        data = yaml.safe_load(f)

    if not isinstance(data, dict):
        raise ValueError(f"Invalid catalog format in {path}: expected dict at root")

    # Validate schema version
    schema_version = data.get("schema_version", "1.0")
    if not isinstance(schema_version, str) or not schema_version.startswith("1."):
        raise ValueError(
            f"Unsupported schema version '{schema_version}' in {path}. "
            f"Expected version 1.x"
        )

    # Basic structure validation
    if "sizes" not in data:
        raise ValueError(f"Catalog {path} missing required 'sizes' section")

    # Record source path for debugging
    data["_source_path"] = str(path)

    return data


def load_catalog(
    thread_series: str,
    custom_path: Optional[Path] = None
) -> Dict[str, Any]:
    """Load fastener catalog for a thread series.

    Args:
        thread_series: One of "metric_coarse", "metric_fine",
                       "unified_coarse", "unified_fine"
        custom_path: Optional explicit path to YAML file (overrides search)

    Returns:
        Parsed catalog dictionary containing:
            - schema_version: str
            - standard: str (e.g., "ISO 262")
            - thread_series: str
            - tolerance_classes: dict of class definitions
            - sizes: dict of size data
            - _source_path: str (path catalog was loaded from)

    Raises:
        FileNotFoundError: If no catalog found for the series
        ValueError: If catalog has invalid format

    Search order (unless custom_path specified):
        1. $YAPCAD_FASTENER_DATA directories
        2. ~/.config/yapcad/fasteners/
        3. Bundled data
    """
    if thread_series not in THREAD_SERIES:
        raise ValueError(
            f"Unknown thread series '{thread_series}'. "
            f"Available: {list(THREAD_SERIES.keys())}"
        )

    custom_str = str(custom_path) if custom_path else None
    return _load_catalog_cached(thread_series, custom_str)


def list_available_sizes(
    thread_series: str,
    custom_path: Optional[Path] = None
) -> List[str]:
    """List all available sizes for a thread series.

    Args:
        thread_series: Thread series name
        custom_path: Optional custom catalog path

    Returns:
        List of size strings (e.g., ["M3", "M4", "M5", ...])
    """
    catalog = load_catalog(thread_series, custom_path)
    return sorted(catalog.get("sizes", {}).keys())


def list_tolerance_classes(
    thread_series: str,
    custom_path: Optional[Path] = None
) -> List[str]:
    """List available tolerance classes for a thread series.

    Args:
        thread_series: Thread series name
        custom_path: Optional custom catalog path

    Returns:
        List of tolerance class strings (e.g., ["6g", "6H", "4g6g"])
    """
    catalog = load_catalog(thread_series, custom_path)
    return sorted(catalog.get("tolerance_classes", {}).keys())


def get_thread_data(
    thread_series: str,
    size: str,
    tolerance_class: Optional[str] = None,
    custom_path: Optional[Path] = None
) -> Dict[str, Any]:
    """Get thread dimension data for a specific size.

    Args:
        thread_series: Thread series name
        size: Size designation (e.g., "M8", "1/4-20")
        tolerance_class: Optional tolerance class to validate exists
        custom_path: Optional custom catalog path

    Returns:
        Dict with thread dimensions:
            - nominal_diameter: float (mm)
            - pitch: float (mm)
            - basic_minor_diameter: float (mm)
            - basic_pitch_diameter: float (mm)
            - tolerances: dict of tolerance class data
            - hex_bolt: dict of bolt head dimensions (if present)
            - hex_nut: dict of nut dimensions (if present)

    Raises:
        KeyError: If size or tolerance_class not found
    """
    catalog = load_catalog(thread_series, custom_path)
    sizes = catalog.get("sizes", {})

    # Normalize size string (handle case variations)
    size_upper = size.upper()
    size_lower = size.lower()

    # Try exact match first, then case variations
    if size in sizes:
        size_data = sizes[size]
    elif size_upper in sizes:
        size_data = sizes[size_upper]
    elif size_lower in sizes:
        size_data = sizes[size_lower]
    else:
        available = sorted(sizes.keys())
        raise KeyError(
            f"Size '{size}' not found in {thread_series}.\n"
            f"Available sizes: {available}"
        )

    # Validate tolerance class if specified
    if tolerance_class:
        tolerances = size_data.get("tolerances", {})
        if tolerance_class not in tolerances:
            available = sorted(tolerances.keys())
            raise KeyError(
                f"Tolerance class '{tolerance_class}' not available for {size}.\n"
                f"Available classes: {available}"
            )

    return size_data


def get_bolt_data(
    thread_series: str,
    size: str,
    tolerance_class: str = "6g",
    custom_path: Optional[Path] = None
) -> Dict[str, Any]:
    """Get complete data for creating a hex bolt.

    Args:
        thread_series: Thread series name
        size: Size designation
        tolerance_class: Thread tolerance class (default "6g" for metric)
        custom_path: Optional custom catalog path

    Returns:
        Dict with all data needed to create a bolt:
            - thread: dict with diameter, pitch, tolerances
            - head: dict with width_across_flats, height, etc.

    Raises:
        KeyError: If size not found or no bolt data available
    """
    size_data = get_thread_data(thread_series, size, tolerance_class, custom_path)

    if "hex_bolt" not in size_data:
        raise KeyError(f"No hex bolt data available for {size}")

    tolerances = size_data.get("tolerances", {}).get(tolerance_class, {})

    return {
        "thread": {
            "nominal_diameter": size_data["nominal_diameter"],
            "pitch": size_data["pitch"],
            "basic_minor_diameter": size_data.get("basic_minor_diameter"),
            "basic_pitch_diameter": size_data.get("basic_pitch_diameter"),
            "tolerance_class": tolerance_class,
            "tolerances": tolerances,
        },
        "head": size_data["hex_bolt"],
    }


def get_nut_data(
    thread_series: str,
    size: str,
    tolerance_class: str = "6H",
    custom_path: Optional[Path] = None
) -> Dict[str, Any]:
    """Get complete data for creating a hex nut.

    Args:
        thread_series: Thread series name
        size: Size designation
        tolerance_class: Thread tolerance class (default "6H" for metric internal)
        custom_path: Optional custom catalog path

    Returns:
        Dict with all data needed to create a nut:
            - thread: dict with diameter, pitch, tolerances
            - body: dict with width_across_flats, thickness, etc.

    Raises:
        KeyError: If size not found or no nut data available
    """
    size_data = get_thread_data(thread_series, size, tolerance_class, custom_path)

    if "hex_nut" not in size_data:
        raise KeyError(f"No hex nut data available for {size}")

    tolerances = size_data.get("tolerances", {}).get(tolerance_class, {})

    return {
        "thread": {
            "nominal_diameter": size_data["nominal_diameter"],
            "pitch": size_data["pitch"],
            "basic_minor_diameter": size_data.get("basic_minor_diameter"),
            "basic_pitch_diameter": size_data.get("basic_pitch_diameter"),
            "tolerance_class": tolerance_class,
            "tolerances": tolerances,
            "internal": True,
        },
        "body": size_data["hex_nut"],
    }
