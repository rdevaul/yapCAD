"""Core `.ycpkg` packaging helpers."""

from __future__ import annotations

import datetime as _dt
import hashlib
import shutil
import json
import math
import re
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
ASSEMBLY_PACKAGE_SCHEMA = "ycpkg-spec-v0.2"
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
    subdirectories = (
        "geometry", "geometry/entities", "instances", "metadata",
        "validation/plans", "validation/results", "exports", "attachments",
    )
    for sub in subdirectories:
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

        for component in self.data.get("components", []) or []:
            info = component.get("geometry")
            if info:
                path = self.root / info["path"]
                if path.exists():
                    info["hash"] = _compute_hash(path, algorithm)
        for entry_key in ("assembly", "bom"):
            info = self.data.get(entry_key)
            if info and info.get("path"):
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

    def get_component(self, component_id: str) -> Optional[Dict[str, Any]]:
        """Return a v0.2 component definition by stable ID."""
        return next(
            (item for item in self.data.get("components", [])
             if item.get("id") == component_id),
            None,
        )

    def load_component_geometry(self, component_id: str) -> List[list]:
        """Load canonical local geometry for a v0.2 component."""
        component = self.get_component(component_id)
        if component is None:
            raise KeyError(f"component not found: {component_id}")
        geometry = component.get("geometry")
        if not geometry:
            return []
        with (self.root / geometry["path"]).open("r", encoding="utf-8") as fp:
            return geometry_from_json(json.load(fp))

    def load_assembly_record(self) -> Dict[str, Any]:
        """Load the native v0.2 assembly graph document."""
        entry = self.data.get("assembly") or {}
        if not entry.get("path"):
            raise ValueError("manifest missing assembly.path")
        with (self.root / entry["path"]).open("r", encoding="utf-8") as fp:
            return json.load(fp)

    def load_bom(self) -> Dict[str, Any]:
        """Load the derived v0.2 engineering BOM document."""
        entry = self.data.get("bom") or {}
        if not entry.get("path"):
            raise ValueError("manifest missing bom.path")
        with (self.root / entry["path"]).open("r", encoding="utf-8") as fp:
            return json.load(fp)


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


_COMPONENT_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
_DISPOSITIONS = {"make", "buy", "raw_stock", "consumable"}


def _json_write(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fp:
        json.dump(value, fp, indent=2, sort_keys=False)
        fp.write("\n")


def _datum_record(datum: Any) -> Dict[str, Any]:
    kind = getattr(datum.datum_type, "value", str(datum.datum_type))
    result: Dict[str, Any] = {
        "id": datum.name,
        "kind": kind,
        "origin": [float(value) for value in datum.origin[:3]],
    }
    for field_name in ("direction", "normal", "x_axis", "y_axis"):
        value = getattr(datum, field_name, None)
        if value is not None:
            result[field_name] = [float(item) for item in value[:3]]
    if getattr(datum, "radius", None) is not None:
        result["radius"] = float(datum.radius)
    return result


def _limits_record(limits: Any) -> Optional[Dict[str, float]]:
    if limits is None:
        return None
    result = {}
    for name in ("min_value", "max_value", "max_velocity", "max_effort"):
        value = getattr(limits, name, None)
        if value is not None:
            result[name] = float(value)
    return result or None


def _native_assembly_record(assembly: Any, root_part: Optional[str]) -> Dict[str, Any]:
    parts: Dict[str, Any] = {}
    for instance_name, definition in assembly.parts.items():
        transform = assembly.transforms.get(instance_name)
        parts[instance_name] = {
            "component": definition.component_id or definition.name,
            "transform": (
                [[float(value) for value in row] for row in transform.tolist()]
                if transform is not None else None
            ),
            "datums": [_datum_record(datum) for datum in definition.datums.values()],
        }
    mates = []
    for mate in assembly.mates:
        entry = {
            "id": mate.name,
            "kind": mate.mate_type.value,
            "partA": mate.part_a,
            "datumA": mate.datum_a,
            "partB": mate.part_b,
            "datumB": mate.datum_b,
            "offset": float(mate.offset),
            "angle": float(mate.angle),
        }
        limits = _limits_record(mate.limits)
        if limits:
            entry["limits"] = limits
        mates.append(entry)
    return {
        "schema": "yapcad-assembly-v0.1",
        "name": assembly.name,
        "rootPart": root_part,
        "solved": bool(getattr(assembly, "_solved", False)),
        "parts": parts,
        "mates": mates,
        "jointValues": {
            name: float(value)
            for name, value in getattr(assembly, "_joint_values", {}).items()
        },
        "jointCouplings": [
            coupling.to_dict() for coupling in assembly.joint_couplings
        ],
    }


def _component_record(definition: Any) -> Dict[str, Any]:
    component_id = definition.component_id or definition.name
    if not _COMPONENT_ID_RE.fullmatch(component_id):
        raise ValueError(
            f"invalid component id {component_id!r}; use letters, digits, '.', '_' or '-'"
        )
    disposition = getattr(definition, "disposition", "make")
    if disposition not in _DISPOSITIONS:
        raise ValueError(
            f"invalid component disposition {disposition!r} for {component_id!r}"
        )
    quantity = float(getattr(definition, "quantity_per_instance", 1.0))
    if not math.isfinite(quantity) or quantity <= 0:
        raise ValueError(f"component {component_id!r} quantity must be positive")
    record: Dict[str, Any] = {
        "id": component_id,
        "name": definition.component_name or definition.name,
        "description": definition.description,
        "disposition": disposition,
        "quantityPerInstance": quantity,
        "unit": getattr(definition, "unit", "each"),
    }
    optional = {
        "partNumber": getattr(definition, "part_number", None),
        "revision": getattr(definition, "revision", None),
        "material": getattr(definition, "material", None),
        "manufacturing": getattr(definition, "manufacturing", None),
        "procurement": getattr(definition, "procurement", None),
        "pmi": getattr(definition, "pmi", None),
    }
    record.update({key: value for key, value in optional.items() if value})
    return record


def _geometry_fingerprint(solid: Any) -> str:
    """Hash tessellated local geometry while excluding volatile metadata IDs."""
    def surface_payload(surface: Any) -> Dict[str, Any]:
        return {
            "vertices": surface[1],
            "normals": surface[2],
            "faces": surface[3],
        }

    surfaces = []

    def collect(value: Any) -> None:
        if issurface(value):
            surfaces.append(surface_payload(value))
        elif isinstance(value, (list, tuple)):
            for item in value:
                collect(item)

    collect(solid[1:])
    payload = {"surfaces": surfaces}
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def _derive_bom(
    components: List[Dict[str, Any]], instances: List[Dict[str, Any]],
) -> Dict[str, Any]:
    by_id = {component["id"]: component for component in components}
    quantities = {component_id: 0.0 for component_id in by_id}
    for instance in instances:
        component = by_id[instance["component"]]
        quantities[component["id"]] += float(component["quantityPerInstance"])
    items = []
    for item_number, component in enumerate(components, start=1):
        quantity = quantities[component["id"]]
        if quantity.is_integer():
            quantity = int(quantity)
        item = {
            "item": item_number,
            "component": component["id"],
            "partNumber": component.get("partNumber"),
            "revision": component.get("revision"),
            "description": component["name"],
            "disposition": component["disposition"],
            "quantity": quantity,
            "unit": component["unit"],
        }
        items.append(item)
    return {"schema": "yapcad-bom-v0.1", "generatedFrom": "instances", "items": items}


def create_package_from_assembly(
    assembly: Any,
    target_dir: Path | str,
    *,
    name: str,
    version: str,
    root_part: Optional[str] = None,
    description: Optional[str] = None,
    author: Optional[str] = None,
    units: str = "mm",
    materials: Optional[Dict[str, Dict[str, Any]]] = None,
    generator: Optional[Dict[str, Any]] = None,
    overwrite: bool = False,
    hash_algorithm: str = "sha256",
) -> PackageManifest:
    """Create a v0.2 product-definition package from an Assembly."""
    if not assembly.parts:
        raise ValueError("assembly has no parts")
    root = Path(target_dir)
    if root.exists():
        if not overwrite and any(root.iterdir()):
            raise FileExistsError(f"target directory {root} already exists and is not empty")
    else:
        root.mkdir(parents=True)
    _ensure_subdirs(root)

    primary_path = root / "geometry" / "primary.json"
    primary_info = _serialize_geometry(
        [assembly.compound_geometry(strict=False)], primary_path, root
    )
    primary_info["hash"] = _compute_hash(primary_path, hash_algorithm)

    components: List[Dict[str, Any]] = []
    component_by_id: Dict[str, Dict[str, Any]] = {}
    component_fingerprints: Dict[str, str] = {}
    instances: List[Dict[str, Any]] = []
    for instance_name, definition in assembly.parts.items():
        record = _component_record(definition)
        component_id = record["id"]
        has_geometry = instance_name in assembly.geometry
        if component_id not in component_by_id:
            if not has_geometry and record["disposition"] == "make":
                raise ValueError(f"component {component_id!r} has no geometry")
            if has_geometry:
                component_path = root / "geometry" / "entities" / f"{component_id}.json"
                geometry_info = _serialize_geometry(
                    [assembly.get_part_geometry(instance_name, positioned=False)],
                    component_path, root,
                )
                geometry_info["hash"] = _compute_hash(component_path, hash_algorithm)
                record["geometry"] = geometry_info
                component_fingerprints[component_id] = _geometry_fingerprint(
                    assembly.geometry[instance_name]
                )
            component_by_id[component_id] = record
            components.append(record)
        else:
            existing = component_by_id[component_id]
            for field_name in (
                "name", "description", "disposition", "quantityPerInstance",
                "unit", "partNumber", "revision", "material", "manufacturing",
                "procurement", "pmi",
            ):
                if existing.get(field_name) != record.get(field_name):
                    raise ValueError(
                        f"component {component_id!r} has conflicting {field_name} metadata"
                    )
            existing_has_geometry = component_id in component_fingerprints
            if has_geometry != existing_has_geometry:
                raise ValueError(
                    f"component {component_id!r} is used with conflicting local geometry"
                )
            if has_geometry:
                fingerprint = _geometry_fingerprint(assembly.geometry[instance_name])
                if fingerprint != component_fingerprints[component_id]:
                    raise ValueError(
                        f"component {component_id!r} is used with conflicting local geometry"
                    )
        transform = assembly.transforms[instance_name]
        instances.append({
            "id": instance_name,
            "component": component_id,
            "transform": [[float(value) for value in row] for row in transform.tolist()],
        })

    assembly_path = root / "metadata" / "assembly.json"
    _json_write(assembly_path, _native_assembly_record(assembly, root_part))
    bom_path = root / "metadata" / "bom.json"
    _json_write(bom_path, _derive_bom(components, instances))

    manifest_data: Dict[str, Any] = {
        "schema": ASSEMBLY_PACKAGE_SCHEMA,
        "id": str(uuid.uuid4()),
        "name": name,
        "version": version,
        "description": description or "",
        "created": {"timestamp": _now_iso()},
        "generator": generator or {"tool": "yapCAD", "version": _yapcad_version},
        "units": units,
        "tags": [],
        "product": {"rootAssembly": assembly.name, "lifecycle": "prototype"},
        "geometry": {"primary": primary_info},
        "components": components,
        "instances": instances,
        "assembly": {
            "path": str(assembly_path.relative_to(root)),
            "hash": _compute_hash(assembly_path, hash_algorithm),
            "rootPart": root_part,
        },
        "bom": {
            "path": str(bom_path.relative_to(root)),
            "hash": _compute_hash(bom_path, hash_algorithm),
            "generatedFrom": "instances",
        },
    }
    if author:
        manifest_data["created"]["author"] = author
    if materials:
        manifest_data["materials"] = materials
    else:
        material_refs = sorted({
            component["material"] for component in components
            if component.get("material")
        })
        if material_refs:
            manifest_data["materials"] = {
                material_ref: {
                    "source": {
                        "type": "custom",
                        "custom": {
                            "notes": "Placeholder - define engineering material properties"
                        },
                    },
                    "visual": {
                        "color": [0.6, 0.85, 1.0],
                        "metallic": 0.0,
                        "roughness": 0.5,
                    },
                }
                for material_ref in material_refs
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
    "create_package_from_assembly",
    "ASSEMBLY_PACKAGE_SCHEMA",
    "load_geometry",
    "_compute_hash",
    "add_geometry_file",
]
