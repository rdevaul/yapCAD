"""Validation utilities for `.ycpkg` packages."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Tuple

from yapcad.io.geometry_json import geometry_from_json
from .core import PackageManifest, _compute_hash

# Valid source types per material schema spec
VALID_SOURCE_TYPES = {"standard", "vendor", "custom", "tested"}


def _check_file(path: Path, expected_hash: str | None) -> Tuple[bool, List[str]]:
    messages: List[str] = []
    if not path.exists():
        messages.append(f"ERROR: missing file {path}")
        return False, messages
    actual_hash = _compute_hash(path)
    if expected_hash and expected_hash.lower() != actual_hash.lower():
        messages.append(
            f"ERROR: hash mismatch for {path} (expected {expected_hash}, got {actual_hash})"
        )
        return False, messages
    return True, messages


def _validate_material(mat_id: str, mat_def: Dict[str, Any], strict: bool) -> Tuple[bool, List[str]]:
    """Validate a single material definition per the schema spec."""
    messages: List[str] = []
    ok = True

    # source.type is required
    source = mat_def.get("source")
    if not source:
        messages.append(f"ERROR: material '{mat_id}' missing required 'source' field")
        return False, messages

    source_type = source.get("type")
    if not source_type:
        messages.append(f"ERROR: material '{mat_id}' missing required 'source.type' field")
        ok = False
    elif source_type not in VALID_SOURCE_TYPES:
        messages.append(f"ERROR: material '{mat_id}' has invalid source.type '{source_type}'")
        ok = False
    else:
        # Validate source type-specific fields
        if source_type == "standard":
            std = source.get("standard", {})
            if not std.get("body"):
                messages.append(f"WARNING: material '{mat_id}' (type=standard) missing source.standard.body")
            if not std.get("designation"):
                messages.append(f"WARNING: material '{mat_id}' (type=standard) missing source.standard.designation")

    # Validate visual properties if present
    visual = mat_def.get("visual", {})
    color = visual.get("color")
    if color:
        if not isinstance(color, list) or len(color) != 3:
            messages.append(f"ERROR: material '{mat_id}' visual.color must be [R, G, B] array")
            ok = False
        elif not all(isinstance(c, (int, float)) and 0 <= c <= 1 for c in color):
            messages.append(f"WARNING: material '{mat_id}' visual.color values should be in [0, 1] range")

    metallic = visual.get("metallic")
    if metallic is not None:
        if not isinstance(metallic, (int, float)) or not (0 <= metallic <= 1):
            messages.append(f"WARNING: material '{mat_id}' visual.metallic should be in [0, 1] range")

    roughness = visual.get("roughness")
    if roughness is not None:
        if not isinstance(roughness, (int, float)) or not (0 <= roughness <= 1):
            messages.append(f"WARNING: material '{mat_id}' visual.roughness should be in [0, 1] range")

    # In strict mode, warn about missing engineering properties
    if strict and not mat_def.get("properties"):
        messages.append(f"WARNING: material '{mat_id}' has no physical properties defined")

    return ok, messages


def _validate_materials(materials: Dict[str, Any], strict: bool) -> Tuple[bool, List[str]]:
    """Validate all materials in the manifest."""
    messages: List[str] = []
    overall_ok = True

    if not isinstance(materials, dict):
        return False, ["ERROR: materials section must be a dictionary"]

    for mat_id, mat_def in materials.items():
        if not isinstance(mat_def, dict):
            messages.append(f"ERROR: material '{mat_id}' definition must be a dictionary")
            overall_ok = False
            continue
        ok, mat_messages = _validate_material(mat_id, mat_def, strict)
        messages.extend(mat_messages)
        if not ok:
            overall_ok = False

    return overall_ok, messages


def validate_package(path: Path | str, *, strict: bool = False) -> Tuple[bool, List[str]]:
    """Validate manifest and referenced artefacts.

    Returns (is_valid, messages). Messages are strings with severity prefixes.
    """
    pkg_path = Path(path)
    messages: List[str] = []
    try:
        manifest = PackageManifest.load(pkg_path)
    except Exception as exc:
        return False, [f"ERROR: failed to load manifest: {exc}"]

    data = manifest.data
    if data.get("schema") != "ycpkg-spec-v0.1":
        messages.append(f"ERROR: unsupported package schema {data.get('schema')}")

    # Geometry primary
    try:
        primary_path = manifest.geometry_primary_path()
    except Exception as exc:
        messages.append(f"ERROR: geometry.primary missing: {exc}")
        return False, messages

    ok, file_messages = _check_file(primary_path, data.get("geometry", {}).get("primary", {}).get("hash"))
    messages.extend(file_messages)
    if ok:
        try:
            with primary_path.open("r", encoding="utf-8") as fp:
                doc = json.load(fp)
            geometry_from_json(doc)
        except Exception as exc:
            messages.append(f"ERROR: invalid geometry JSON: {exc}")
            ok = False

    overall_ok = ok and not any(msg.startswith("ERROR") for msg in messages)

    # Derived geometry
    derived_entries = data.get("geometry", {}).get("derived", []) or []
    for entry in derived_entries:
        derived_path = manifest.root / entry["path"]
        ok_entry, msgs = _check_file(derived_path, entry.get("hash"))
        messages.extend(msgs)
        overall_ok = overall_ok and ok_entry

    # Exports and attachments
    for section in ("exports", "attachments"):
        for entry in data.get(section, []) or []:
            file_path = manifest.root / entry["path"]
            ok_entry, msgs = _check_file(file_path, entry.get("hash"))
            messages.extend(msgs)
            overall_ok = overall_ok and ok_entry
            if strict and entry.get("hash") is None:
                messages.append(f"WARNING: {section} entry {entry.get('id')} missing hash")

    # Materials validation
    materials = data.get("materials")
    if materials:
        mat_ok, mat_msgs = _validate_materials(materials, strict)
        messages.extend(mat_msgs)
        overall_ok = overall_ok and mat_ok

    if overall_ok:
        messages.insert(0, f"OK: {pkg_path} passed validation")
    else:
        messages.insert(0, f"FAILED: {pkg_path} has validation errors")
    return overall_ok, messages


__all__ = ["validate_package"]
