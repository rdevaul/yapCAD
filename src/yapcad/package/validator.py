"""Validation utilities for `.ycpkg` packages."""

from __future__ import annotations

import json
from pathlib import Path
from typing import List, Tuple

from yapcad.io.geometry_json import geometry_from_json
from .core import PackageManifest, _compute_hash


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

    if overall_ok:
        messages.insert(0, f"OK: {pkg_path} passed validation")
    else:
        messages.insert(0, f"FAILED: {pkg_path} has validation errors")
    return overall_ok, messages


__all__ = ["validate_package"]
