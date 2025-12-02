"""Core `.ycpkg` packaging helpers."""

from __future__ import annotations

import datetime as _dt
import hashlib
import shutil
import json
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence

from yapcad import __version__ as _yapcad_version
from yapcad.geom3d import issolid, issurface
from yapcad.io.geometry_json import SCHEMA_ID as GEOMETRY_SCHEMA, geometry_from_json, geometry_to_json
from yapcad.metadata import (
    get_solid_metadata,
    get_surface_metadata,
)

PACKAGE_SCHEMA = "ycpkg-spec-v0.1"
MANIFEST_FILENAME = "manifest.yaml"


def _compute_hash(path: Path, algorithm: str = "sha256") -> str:
    h = hashlib.new(algorithm)
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return f"{algorithm}:{h.hexdigest()}"


def _now_iso() -> str:
    return _dt.datetime.now(_dt.timezone.utc).isoformat()


def _collect_tags(entities: Iterable[list]) -> List[str]:
    tags: List[str] = []
    seen = set()
    for entity in entities:
        if issolid(entity):
            meta = get_solid_metadata(entity, create=False) or {}
        elif issurface(entity):
            meta = get_surface_metadata(entity, create=False) or {}
        else:
            continue
        for tag in meta.get("tags", []):
            if tag not in seen:
                tags.append(tag)
                seen.add(tag)
    return tags


def _collect_material_refs(entities: Iterable[list]) -> List[str]:
    """Collect unique material references from entities."""
    refs: List[str] = []
    seen = set()
    for entity in entities:
        if issolid(entity):
            meta = get_solid_metadata(entity, create=False) or {}
        elif issurface(entity):
            meta = get_surface_metadata(entity, create=False) or {}
        else:
            continue
        mat_ref = meta.get("material")
        if mat_ref and mat_ref not in seen:
            refs.append(mat_ref)
            seen.add(mat_ref)
    return refs


def _ensure_subdirs(root: Path) -> None:
    for sub in ("geometry", "metadata", "validation/plans", "validation/results", "exports", "attachments"):
        (root / sub).mkdir(parents=True, exist_ok=True)


def _serialize_geometry(entities: Sequence[list], target: Path, root: Path) -> Dict[str, Any]:
    doc = geometry_to_json(entities)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8") as fp:
        json.dump(doc, fp, indent=2, sort_keys=False)
        fp.write("\n")
    entity_ids = [entry["id"] for entry in doc.get("entities", []) if entry.get("id")]
    return {
        "path": str(target.relative_to(root)),
        "schema": GEOMETRY_SCHEMA,
        "entities": entity_ids,
    }


@dataclass
class PackageManifest:
    """Wrapper around the manifest document."""

    root: Path
    data: Dict[str, Any] = field(default_factory=dict)
    manifest_name: str = MANIFEST_FILENAME

    @property
    def manifest_path(self) -> Path:
        return self.root / self.manifest_name

    @classmethod
    def load(cls, package_path: Path | str) -> "PackageManifest":
        root = Path(package_path)
        manifest_path = root / MANIFEST_FILENAME
        if not manifest_path.exists():
            raise FileNotFoundError(f"manifest not found: {manifest_path}")
        with manifest_path.open("r", encoding="utf-8") as fp:
            data = json.load(fp) if manifest_path.suffix == ".json" else None
        if data is None:
            import yaml  # local import to avoid hard dependency if unused
            with manifest_path.open("r", encoding="utf-8") as fp:
                data = yaml.safe_load(fp) or {}
        return cls(root=root, data=data)

    def save(self) -> None:
        self.data.setdefault("schema", PACKAGE_SCHEMA)
        self.root.mkdir(parents=True, exist_ok=True)
        import yaml

        with self.manifest_path.open("w", encoding="utf-8") as fp:
            yaml.safe_dump(self.data, fp, sort_keys=False)

    def recompute_hashes(self, *, algorithm: str = "sha256") -> None:
        geom = self.data.get("geometry", {})
        for section in ("primary",):
            info = geom.get(section)
            if info:
                path = self.root / info["path"]
                if path.exists():
                    info["hash"] = _compute_hash(path, algorithm)
        for key in ("derived",):
            items = geom.get(key, []) or []
            for info in items:
                path = self.root / info["path"]
                if path.exists():
                    info["hash"] = _compute_hash(path, algorithm)

        for entry_key in ("exports", "attachments"):
            for info in self.data.get(entry_key, []) or []:
                path = self.root / info["path"]
                if path.exists():
                    info["hash"] = _compute_hash(path, algorithm)

    def geometry_primary_path(self) -> Path:
        geom = self.data.get("geometry", {}).get("primary")
        if not geom:
            raise ValueError("manifest missing geometry.primary section")
        return self.root / geom["path"]

    def get_materials(self) -> Dict[str, Any]:
        """Return the materials dictionary from the manifest."""
        return self.data.get("materials", {})

    def get_material(self, material_id: str) -> Optional[Dict[str, Any]]:
        """Get a specific material definition by ID."""
        return self.data.get("materials", {}).get(material_id)


def create_package_from_entities(
    entities: Sequence[list],
    target_dir: Path | str,
    *,
    name: str,
    version: str,
    description: Optional[str] = None,
    author: Optional[str] = None,
    units: Optional[str] = None,
    materials: Optional[Dict[str, Dict[str, Any]]] = None,
    generator: Optional[Dict[str, Any]] = None,
    overwrite: bool = False,
    hash_algorithm: str = "sha256",
) -> PackageManifest:
    if not entities:
        raise ValueError("no entities supplied for packaging")

    root = Path(target_dir)
    if root.exists():
        if not overwrite and any(root.iterdir()):
            raise FileExistsError(f"target directory {root} already exists and is not empty")
    else:
        root.mkdir(parents=True)

    _ensure_subdirs(root)
    primary_path = root / "geometry" / "primary.json"
    geometry_info = _serialize_geometry(entities, primary_path, root)
    geometry_info["hash"] = _compute_hash(primary_path, hash_algorithm)

    tags = _collect_tags(entities)
    material_refs = _collect_material_refs(entities)
    manifest_data: Dict[str, Any] = {
        "schema": PACKAGE_SCHEMA,
        "id": str(uuid.uuid4()),
        "name": name,
        "version": version,
        "description": description or "",
        "created": {
            "timestamp": _now_iso(),
        },
        "generator": generator
        or {
            "tool": "yapCAD",
            "version": _yapcad_version,
        },
        "units": units or "mm",
        "tags": tags,
        "geometry": {
            "primary": geometry_info,
        },
    }
    if author:
        manifest_data["created"]["author"] = author
    # Add materials section if provided or if entities reference materials
    if materials:
        manifest_data["materials"] = materials
    elif material_refs:
        # Entity references materials but none provided - create placeholder entries
        manifest_data["materials"] = {
            ref: {
                "source": {"type": "custom", "custom": {"notes": "Placeholder - define material properties"}},
                "visual": {"color": [0.6, 0.85, 1.0], "metallic": 0.0, "roughness": 0.5},
            }
            for ref in material_refs
        }
    manifest = PackageManifest(root=root, data=manifest_data)
    manifest.save()
    return manifest


def load_geometry(manifest: PackageManifest) -> List[list]:
    primary_path = manifest.geometry_primary_path()
    with primary_path.open("r", encoding="utf-8") as fp:
        doc = json.load(fp)
    return geometry_from_json(doc)


def add_geometry_file(
    manifest: PackageManifest,
    source: Path | str,
    *,
    dest_relative: str | None = None,
    purpose: Optional[str] = None,
    category: str = "derived",
    overwrite: bool = False,
    metadata: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Copy an external geometry file (e.g., STEP/STL) into the package and record it.

    Args:
        manifest: Loaded manifest wrapper.
        source: Path to the external file that should be bundled.
        dest_relative: Optional relative destination path inside the package root.
            Defaults to ``geometry/derived/<source.name>``.
        purpose: Optional description stored alongside the entry.
        category: Manifest section to update. Supported: ``"derived"`` (default), ``"attachments"``.
        overwrite: Allow replacing an existing file at the target location.
        metadata: Additional key/value pairs merged into the manifest entry.

    Returns:
        The manifest entry dictionary that was inserted.
    """

    src_path = Path(source)
    if not src_path.exists():
        raise FileNotFoundError(f"geometry source not found: {src_path}")

    if dest_relative is None:
        if category == "derived":
            dest_relative_path = Path("geometry") / "derived" / src_path.name
        elif category == "attachments":
            dest_relative_path = Path("attachments") / src_path.name
        else:
            dest_relative_path = Path(src_path.name)
    else:
        dest_relative_path = Path(dest_relative)
        if dest_relative_path.is_absolute():
            raise ValueError("dest_relative must be a relative path")

    dest_path = manifest.root / dest_relative_path
    dest_path.parent.mkdir(parents=True, exist_ok=True)

    if dest_path.exists() and not overwrite:
        raise FileExistsError(f"target file already exists: {dest_path}")

    shutil.copy2(src_path, dest_path)

    entry: Dict[str, Any] = {
        "path": str(dest_relative_path.as_posix()),
        "hash": _compute_hash(dest_path),
        "format": src_path.suffix.lstrip(".").lower(),
        "source": {
            "kind": "import",
            "original": str(src_path),
        },
    }
    if purpose:
        entry["purpose"] = purpose
    if metadata:
        entry.update(metadata)

    if category == "derived":
        geometry = manifest.data.setdefault("geometry", {})
        derived = geometry.setdefault("derived", [])
        derived = [item for item in derived if item.get("path") != entry["path"]]
        derived.append(entry)
        geometry["derived"] = derived
    elif category == "attachments":
        attachments = manifest.data.setdefault("attachments", [])
        attachments = [item for item in attachments if item.get("path") != entry["path"]]
        entry.setdefault("id", dest_relative_path.stem)
        attachments.append(entry)
        manifest.data["attachments"] = attachments
    else:
        raise ValueError(f"unsupported category for geometry file: {category}")

    return entry


__all__ = [
    "PACKAGE_SCHEMA",
    "MANIFEST_FILENAME",
    "PackageManifest",
    "create_package_from_entities",
    "load_geometry",
    "_compute_hash",
    "add_geometry_file",
]
