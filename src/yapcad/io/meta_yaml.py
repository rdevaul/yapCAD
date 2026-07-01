"""Sidecar ``.meta.yaml`` emitter for yapCAD STEP/STL exports.

Implements Step 5 of the metadata-namespace-v1.1 migration (RFC
``rfc_assembly_operation_namespaces.md``).

Usage
-----
Call :func:`dump_metadata_yaml` immediately after writing a STEP or STL file::

    write_step(part, "output/my_part.step")
    dump_metadata_yaml(part, "output/my_part.step")   # → my_part.meta.yaml

The sidecar file is placed next to the geometry file by default.  Pass an
explicit *yaml_path* to override.

The function is a no-op (and returns ``None``) when the part carries no
yapCAD metadata — so callers can always call it unconditionally.
"""

from __future__ import annotations

import datetime
import hashlib
import os
from pathlib import Path
from typing import Any, Dict, Optional

try:
    import yaml as _yaml
    _YAML_AVAILABLE = True
except ImportError:
    _YAML_AVAILABLE = False

try:
    import yapcad
    _YAPCAD_VERSION = getattr(yapcad, "__version__", "unknown")
except Exception:
    _YAPCAD_VERSION = "unknown"


def dump_metadata_yaml(
    part: Any,
    geometry_path: str | os.PathLike,
    *,
    yaml_path: Optional[str | os.PathLike] = None,
    source_dsl: Optional[str] = None,
    source_command: Optional[str] = None,
) -> Optional[str]:
    """Write a sidecar ``.meta.yaml`` file next to a STEP or STL export.

    Parameters
    ----------
    part:
        A yapCAD solid, surface, or other geometry object that may carry a
        metadata dict (retrieved via :func:`yapcad.metadata.get_metadata`).
    geometry_path:
        Path to the already-written STEP/STL file.  Used to derive the
        default sidecar path and to compute ``_provenance.source_sha256``
        when *source_dsl* is not provided.
    yaml_path:
        Explicit output path for the YAML sidecar.  When omitted, the file
        is placed alongside *geometry_path* with the suffix replaced by
        ``.meta.yaml``.
    source_dsl:
        Path to the originating DSL file.  Stored in ``_provenance`` and
        used for the SHA-256 hash when the file exists.
    source_command:
        The DSL command name that produced this part (e.g. ``"FORWARD_BULKHEAD"``).

    Returns
    -------
    str or None
        Absolute path of the written YAML file, or ``None`` when no metadata
        was found on *part* (the sidecar is not written in that case).

    Raises
    ------
    ImportError
        If PyYAML is not installed and metadata *is* present (so the caller
        is aware that the sidecar was not written).
    """
    meta = _extract_metadata(part)
    if meta is None:
        return None

    if not _YAML_AVAILABLE:
        raise ImportError(
            "PyYAML is required for .meta.yaml export. "
            "Install it with: pip install pyyaml"
        )

    geo_path = Path(geometry_path).resolve()

    if yaml_path is None:
        # Strip all suffixes, not just the last one, so "part.step" and
        # "part_CURRENT.step" both produce "part.meta.yaml" / "part_CURRENT.meta.yaml".
        stem = geo_path.name
        for suffix in reversed(geo_path.suffixes):
            stem = stem[: -len(suffix)]
        out_path = geo_path.parent / (stem + ".meta.yaml")
    else:
        out_path = Path(yaml_path).resolve()

    provenance = _build_provenance(geo_path, source_dsl, source_command)

    # Merge provenance into a copy of the metadata dict.
    output: Dict[str, Any] = dict(meta)
    output["_provenance"] = provenance

    # Ensure schema string is at the top for readability.
    ordered: Dict[str, Any] = {}
    for key in ("schema", "entityId", "timestamp", "tags", "layer"):
        if key in output:
            ordered[key] = output.pop(key)
    ordered["_provenance"] = output.pop("_provenance")
    ordered.update(output)

    with open(out_path, "w", encoding="utf-8") as f:
        _yaml.dump(ordered, f, default_flow_style=False, allow_unicode=True, sort_keys=False)

    return str(out_path)


def load_metadata_yaml(yaml_path: str | os.PathLike) -> Dict[str, Any]:
    """Load a ``.meta.yaml`` sidecar and return the raw dict.

    Parameters
    ----------
    yaml_path:
        Path to the ``.meta.yaml`` file.

    Returns
    -------
    dict
        The deserialized metadata dictionary.

    Raises
    ------
    ImportError
        If PyYAML is not installed.
    FileNotFoundError
        If *yaml_path* does not exist.
    """
    if not _YAML_AVAILABLE:
        raise ImportError(
            "PyYAML is required to load .meta.yaml files. "
            "Install it with: pip install pyyaml"
        )
    path = Path(yaml_path)
    if not path.exists():
        raise FileNotFoundError(f"Sidecar not found: {path}")
    with open(path, "r", encoding="utf-8") as f:
        return _yaml.safe_load(f) or {}


def sidecar_path_for(geometry_path: str | os.PathLike) -> Path:
    """Return the conventional sidecar path for a geometry file.

    This is a pure path computation — no file I/O is performed.

    Parameters
    ----------
    geometry_path:
        Path to a STEP or STL file.

    Returns
    -------
    pathlib.Path
        The conventional ``.meta.yaml`` path.
    """
    geo_path = Path(geometry_path)
    stem = geo_path.name
    for suffix in reversed(geo_path.suffixes):
        stem = stem[: -len(suffix)]
    return geo_path.parent / (stem + ".meta.yaml")


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _extract_metadata(part: Any) -> Optional[Dict[str, Any]]:
    """Return the metadata dict from *part*, or ``None`` if absent."""
    try:
        from yapcad.metadata import get_metadata  # type: ignore[import]
        meta = get_metadata(part)
        return meta if meta else None
    except (ImportError, Exception):
        pass

    # Fallback: check common attribute names used in yapCAD internals.
    for attr in ("_metadata", "metadata", "meta"):
        val = getattr(part, attr, None)
        if isinstance(val, dict) and val:
            return val

    return None


def _build_provenance(
    geo_path: Path,
    source_dsl: Optional[str],
    source_command: Optional[str],
) -> Dict[str, Any]:
    """Build the ``_provenance`` block for the sidecar."""
    now = datetime.datetime.now(datetime.timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")

    prov: Dict[str, Any] = {
        "yapcad_version": _YAPCAD_VERSION,
        "exported_at": now,
    }

    if source_dsl:
        prov["source_dsl"] = source_dsl
        dsl_path = Path(source_dsl)
        if dsl_path.exists():
            prov["source_dsl_sha256"] = _sha256_file(dsl_path)
    else:
        # Fall back to hashing the geometry file itself for a weak provenance anchor.
        if geo_path.exists():
            prov["geometry_sha256"] = _sha256_file(geo_path)

    if source_command:
        prov["source_command"] = source_command

    return prov


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


__all__ = ["dump_metadata_yaml", "load_metadata_yaml", "sidecar_path_for"]
