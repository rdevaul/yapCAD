"""Metadata helpers for yapCAD geometry objects."""

from __future__ import annotations

import uuid
from typing import Dict, Tuple, Any, Optional, Iterable

from yapcad.geom3d import issolid, issurface

_SURFACE_META_INDEX = 6
_SOLID_META_INDEX = 4
_DEFAULT_SCHEMA = "metadata-namespace-v1.1"

_ROOT_FIELDS = ("schema", "entityId", "timestamp", "tags", "layer")
_SECTION_KEYS = {
    "material",
    "manufacturing",
    "designHistory",
    "constraints",
    "analysis",
    "assembly",
    "operation",
}

# v1.1 enum vocabularies (RFC § 7-8)
_JOINT_KINDS = {"axial", "radial", "mixed", "none"}
_BOLT_RINGS = {"axial", "radial"}
_DATUM_KINDS = {"bolt_circle", "axis", "plane", "point", "edge"}
_SURFACE_KINDS = {"mating", "bearing", "sealing", "datum", "wear"}
_KEEPOUT_KINDS = {"volume", "axis", "ring"}
_OPERATION_KINDS = {"subtract", "intersect", "union"}
_OPERATION_POLICIES = {"strict", "warn", "ignore"}
_FEATURE_KINDS = {
    "access_panel", "vent", "parachute_door", "fastener_through",
    "wire_pass", "channel", "pocket", "other",
}
_FASTENER_THREADS = {"tap", "clearance", "heat-set", "through"}
_FASTENER_HEADS = {"SHCS", "FHCS", "BHCS", "hex", "stud", "none"}


def _validate_enum(value: Any, allowed: set, *, field: str) -> None:
    if value is None:
        return
    if value not in allowed:
        raise ValueError(
            f"invalid {field}={value!r}; expected one of {sorted(allowed)}"
        )


def _initial_root(entity_id: Optional[str] = None) -> Dict[str, Any]:
    data: Dict[str, Any] = {
        "schema": _DEFAULT_SCHEMA,
        "tags": [],
        "layer": "default",
    }
    if entity_id:
        data["entityId"] = entity_id
    return data


def _ensure_root(meta: Dict[str, Any], entity_id: Optional[str] = None) -> Dict[str, Any]:
    # NOTE: when meta is an empty dict, populate it in place rather than
    # returning a detached fresh dict; otherwise downstream namespace
    # writes vanish on garbage-collection.
    if meta is None:
        return _initial_root(entity_id)
    if "schema" not in meta:
        meta["schema"] = _DEFAULT_SCHEMA
    if entity_id and "entityId" not in meta:
        meta["entityId"] = entity_id
    if "tags" not in meta:
        meta["tags"] = []
    if "layer" not in meta or not meta["layer"]:
        meta["layer"] = "default"
    return meta


def _ensure_namespace(meta: Dict[str, Any], namespace: str) -> Dict[str, Any]:
    if namespace not in meta or not isinstance(meta[namespace], dict):
        meta[namespace] = {}
    return meta[namespace]


def _append_tag(meta: Dict[str, Any], tag: str) -> None:
    if not tag:
        return
    tags = meta.setdefault("tags", [])
    if isinstance(tags, list) and tag not in tags:
        tags.append(tag)


def _ensure_metadata(container: list, index: int, create: bool) -> Dict:
    if len(container) <= index:
        if not create:
            return {}
        while len(container) <= index:
            container.append({})
    meta = container[index]
    if not isinstance(meta, dict):
        if not create:
            return {}
        meta = {}
        container[index] = meta
    return meta


def get_surface_metadata(surface: list, create: bool = False) -> Dict:
    if not issurface(surface):
        raise ValueError('invalid surface passed to get_surface_metadata')
    meta = _ensure_metadata(surface, _SURFACE_META_INDEX, create)
    if meta is None:
        return {}
    entity_id = ensure_surface_id(surface) if create else meta.get("entityId")
    root = _ensure_root(meta, entity_id)
    if len(surface) <= _SURFACE_META_INDEX:
        surface.append(root)
    else:
        surface[_SURFACE_META_INDEX] = root
    return root


def set_surface_metadata(surface: list, metadata: Dict) -> list:
    if not isinstance(metadata, dict):
        raise ValueError('surface metadata must be a dictionary')
    if not issurface(surface):
        raise ValueError('invalid surface passed to set_surface_metadata')
    if len(surface) <= _SURFACE_META_INDEX:
        surface.append(metadata)
    else:
        surface[_SURFACE_META_INDEX] = metadata
    return surface


def ensure_surface_id(surface: list) -> str:
    meta = _ensure_metadata(surface, _SURFACE_META_INDEX, True)
    if meta is None:
        meta = {}
        surface.append(meta)
    meta = _ensure_root(meta)
    if 'entityId' not in meta:
        meta['entityId'] = str(uuid.uuid4())
    if 'id' not in meta:
        meta['id'] = meta['entityId']
    surface[_SURFACE_META_INDEX] = meta
    return meta['id']


def set_surface_units(surface: list, units: str) -> list:
    meta = get_surface_metadata(surface, create=True)
    meta['units'] = units
    return surface


def set_surface_origin(surface: list, origin: str) -> list:
    meta = get_surface_metadata(surface, create=True)
    meta['origin'] = origin
    return surface


def get_solid_metadata(sld: list, create: bool = False) -> Dict:
    if not issolid(sld):
        raise ValueError('invalid solid passed to get_solid_metadata')
    meta = _ensure_metadata(sld, _SOLID_META_INDEX, create)
    if meta is None:
        return {}
    entity_id = ensure_solid_id(sld) if create else meta.get("entityId")
    root = _ensure_root(meta, entity_id)
    if len(sld) <= _SOLID_META_INDEX:
        sld.append(root)
    else:
        sld[_SOLID_META_INDEX] = root
    return root


def set_solid_metadata(sld: list, metadata: Dict) -> list:
    if not isinstance(metadata, dict):
        raise ValueError('solid metadata must be a dictionary')
    if not issolid(sld):
        raise ValueError('invalid solid passed to set_solid_metadata')
    if len(sld) <= _SOLID_META_INDEX:
        sld.append(metadata)
    else:
        sld[_SOLID_META_INDEX] = metadata
    return sld


def ensure_solid_id(sld: list) -> str:
    meta = _ensure_metadata(sld, _SOLID_META_INDEX, True)
    if meta is None:
        meta = {}
        sld.append(meta)
    meta = _ensure_root(meta)
    if 'entityId' not in meta:
        meta['entityId'] = str(uuid.uuid4())
    if 'id' not in meta:
        meta['id'] = meta['entityId']
    sld[_SOLID_META_INDEX] = meta
    return meta['id']


def set_solid_context(sld: list, context: Dict) -> list:
    if not isinstance(context, dict):
        raise ValueError('context must be a dictionary')
    meta = get_solid_metadata(sld, create=True)
    meta['context'] = context
    return sld


def set_layer(meta: Dict[str, Any], layer: str) -> Dict[str, Any]:
    if not layer:
        layer = "default"
    root = _ensure_root(meta)
    root["layer"] = layer
    return root


# ---------------------------------------------------------------------------
# Metadata namespace helpers

def add_tags(meta: Dict[str, Any], tags: Iterable[str]) -> Dict[str, Any]:
    root = _ensure_root(meta)
    for tag in tags:
        _append_tag(root, tag)
    return root


def set_material(meta: Dict[str, Any], *, name: Optional[str] = None,
                 standard: Optional[str] = None, grade: Optional[str] = None,
                 density_kg_m3: Optional[float] = None,
                 source: Optional[str] = None) -> Dict[str, Any]:
    root = _ensure_root(meta)
    section = _ensure_namespace(root, "material")
    if name is not None:
        section["name"] = name
    if standard is not None:
        section["standard"] = standard
    if grade is not None:
        section["grade"] = grade
    if density_kg_m3 is not None:
        section["density_kg_m3"] = float(density_kg_m3)
    if source is not None:
        section["source"] = source
    return root


def set_manufacturing(meta: Dict[str, Any], *,
                      process: Optional[str] = None,
                      instructions: Optional[str] = None,
                      fixtures: Optional[Iterable[str]] = None,
                      layers: Optional[Dict[str, Any]] = None,
                      postprocessing: Optional[Iterable[str]] = None) -> Dict[str, Any]:
    root = _ensure_root(meta)
    section = _ensure_namespace(root, "manufacturing")
    if process is not None:
        section["process"] = process
    if instructions is not None:
        section["instructions"] = instructions
    if fixtures is not None:
        section["fixtures"] = list(fixtures)
    if layers is not None:
        section["layers"] = dict(layers)
    if postprocessing is not None:
        section["postprocessing"] = list(postprocessing)
    return root


def add_design_history_entry(meta: Dict[str, Any], *,
                             author: Optional[str] = None,
                             source: Optional[str] = None,
                             context: Optional[str] = None,
                             tools: Optional[Iterable[str]] = None,
                             revision: Optional[str] = None,
                             timestamp: Optional[str] = None,
                             notes: Optional[str] = None) -> Dict[str, Any]:
    root = _ensure_root(meta)
    section = _ensure_namespace(root, "designHistory")
    if author is not None:
        section["author"] = author
    if source is not None:
        section["source"] = source
    if context is not None:
        section["context"] = context
    if tools is not None:
        section["tools"] = list(tools)
    entries = section.setdefault("iterations", [])
    entry: Dict[str, Any] = {}
    if revision is not None:
        entry["revision"] = revision
    if timestamp is not None:
        entry["timestamp"] = timestamp
    if notes is not None:
        entry["notes"] = notes
    if entry:
        entries.append(entry)
    return root


def set_mass_constraint(meta: Dict[str, Any], *,
                        max_kg: Optional[float] = None,
                        target_kg: Optional[float] = None) -> Dict[str, Any]:
    root = _ensure_root(meta)
    section = _ensure_namespace(root, "constraints")
    mass = section.setdefault("mass", {})
    if max_kg is not None:
        mass["max_kg"] = float(max_kg)
    if target_kg is not None:
        mass["target_kg"] = float(target_kg)
    return root


def set_envelope_constraint(meta: Dict[str, Any], *,
                            x_mm: Optional[float] = None,
                            y_mm: Optional[float] = None,
                            z_mm: Optional[float] = None) -> Dict[str, Any]:
    root = _ensure_root(meta)
    section = _ensure_namespace(root, "constraints")
    envelope = section.setdefault("envelope", {})
    if x_mm is not None:
        envelope["x_mm"] = float(x_mm)
    if y_mm is not None:
        envelope["y_mm"] = float(y_mm)
    if z_mm is not None:
        envelope["z_mm"] = float(z_mm)
    return root


def add_compliance(meta: Dict[str, Any], standards: Iterable[str]) -> Dict[str, Any]:
    root = _ensure_root(meta)
    section = _ensure_namespace(root, "constraints")
    listing = section.setdefault("compliance", [])
    for item in standards:
        if item and item not in listing:
            listing.append(item)
    return root


def add_analysis_record(meta: Dict[str, Any], record: Dict[str, Any]) -> Dict[str, Any]:
    root = _ensure_root(meta)
    section = _ensure_namespace(root, "analysis")
    simulations = section.setdefault("simulations", [])
    simulations.append(dict(record))
    return root


# ---------------------------------------------------------------------------
# Assembly namespace (RFC § 7) — v1.1

def get_assembly_metadata(meta: Dict[str, Any], create: bool = False) -> Dict[str, Any]:
    """Return the ``assembly`` namespace dict, creating it when ``create`` is True.

    Returns an empty dict (not attached) if the namespace is absent and
    ``create`` is False, mirroring ``get_*_metadata`` semantics for surfaces
    and solids.
    """
    if not isinstance(meta, dict):
        return {}
    if "assembly" in meta and isinstance(meta["assembly"], dict):
        return meta["assembly"]
    if not create:
        return {}
    root = _ensure_root(meta)
    return _ensure_namespace(root, "assembly")


def set_assembly(meta: Dict[str, Any], *,
                 joint_kind: Optional[str] = None,
                 no_cut: Optional[bool] = None) -> Dict[str, Any]:
    """Set scalar fields on the ``assembly`` namespace. Returns the namespace dict."""
    _validate_enum(joint_kind, _JOINT_KINDS, field="assembly.joint_kind")
    root = _ensure_root(meta)
    section = _ensure_namespace(root, "assembly")
    if joint_kind is not None:
        section["joint_kind"] = joint_kind
    if no_cut is not None:
        section["no_cut"] = bool(no_cut)
    return section


def _build_fastener(spec: Optional[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    if spec is None:
        return None
    if not isinstance(spec, dict):
        raise ValueError("fastener spec must be a dict")
    out: Dict[str, Any] = {}
    if "designation" in spec and spec["designation"] is not None:
        out["designation"] = str(spec["designation"])
    thread = spec.get("thread")
    _validate_enum(thread, _FASTENER_THREADS, field="fastener.thread")
    if thread is not None:
        out["thread"] = thread
    head = spec.get("head")
    _validate_enum(head, _FASTENER_HEADS, field="fastener.head")
    if head is not None:
        out["head"] = head
    return out


def add_bolt_pattern(meta: Dict[str, Any], *,
                     id: str,
                     ring: str,
                     R_mm: float,
                     z_mm: float,
                     count: int,
                     theta0_deg: Optional[float] = None,
                     fastener: Optional[Dict[str, Any]] = None,
                     tolerance_mm: Optional[float] = None) -> Dict[str, Any]:
    """Append a bolt pattern entry to ``assembly.bolt_patterns``. Returns the entry."""
    if not id:
        raise ValueError("bolt_pattern.id is required")
    _validate_enum(ring, _BOLT_RINGS, field="bolt_pattern.ring")
    section = get_assembly_metadata(meta, create=True)
    entry: Dict[str, Any] = {
        "id": str(id),
        "ring": ring,
        "R_mm": float(R_mm),
        "z_mm": float(z_mm),
        "count": int(count),
    }
    if theta0_deg is not None:
        entry["theta0_deg"] = float(theta0_deg)
    fast = _build_fastener(fastener)
    if fast is not None:
        entry["fastener"] = fast
    if tolerance_mm is not None:
        entry["tolerance_mm"] = float(tolerance_mm)
    section.setdefault("bolt_patterns", []).append(entry)
    return entry


def add_datum(meta: Dict[str, Any], *,
              id: str,
              kind: str,
              ring: Optional[str] = None,
              R_mm: Optional[float] = None,
              z_mm: Optional[float] = None,
              direction: Optional[Iterable[float]] = None) -> Dict[str, Any]:
    """Append an explicit datum entry to ``assembly.datums``. Returns the entry."""
    if not id:
        raise ValueError("datum.id is required")
    _validate_enum(kind, _DATUM_KINDS, field="datum.kind")
    if ring is not None:
        _validate_enum(ring, _BOLT_RINGS, field="datum.ring")
    section = get_assembly_metadata(meta, create=True)
    entry: Dict[str, Any] = {"id": str(id), "kind": kind}
    if ring is not None:
        entry["ring"] = ring
    if R_mm is not None:
        entry["R_mm"] = float(R_mm)
    if z_mm is not None:
        entry["z_mm"] = float(z_mm)
    if direction is not None:
        vec = list(direction)
        if len(vec) != 3:
            raise ValueError("datum.direction must have length 3")
        entry["direction"] = [float(c) for c in vec]
    section.setdefault("datums", []).append(entry)
    return entry


def add_surface(meta: Dict[str, Any], *,
                id: str,
                kind: str,
                mate_to: Optional[str] = None,
                finish: Optional[str] = None,
                tolerance_mm: Optional[float] = None) -> Dict[str, Any]:
    """Append a mating/bearing/sealing/datum/wear surface to ``assembly.surfaces``."""
    if not id:
        raise ValueError("surface.id is required")
    _validate_enum(kind, _SURFACE_KINDS, field="surface.kind")
    section = get_assembly_metadata(meta, create=True)
    entry: Dict[str, Any] = {"id": str(id), "kind": kind}
    if mate_to is not None:
        entry["mate_to"] = str(mate_to)
    if finish is not None:
        entry["finish"] = str(finish)
    if tolerance_mm is not None:
        entry["tolerance_mm"] = float(tolerance_mm)
    section.setdefault("surfaces", []).append(entry)
    return entry


def add_keepout(meta: Dict[str, Any], *,
                id: str,
                kind: str,
                shape: Optional[Dict[str, Any]] = None,
                reason: Optional[str] = None) -> Dict[str, Any]:
    """Append a keep-out volume/axis/ring to ``assembly.keepouts``."""
    if not id:
        raise ValueError("keepout.id is required")
    _validate_enum(kind, _KEEPOUT_KINDS, field="keepout.kind")
    section = get_assembly_metadata(meta, create=True)
    entry: Dict[str, Any] = {"id": str(id), "kind": kind}
    if shape is not None:
        if not isinstance(shape, dict):
            raise ValueError("keepout.shape must be a dict")
        entry["shape"] = dict(shape)
    if reason is not None:
        entry["reason"] = str(reason)
    section.setdefault("keepouts", []).append(entry)
    return entry


# ---------------------------------------------------------------------------
# Operation namespace (RFC § 8) — v1.1

def get_operation_metadata(meta: Dict[str, Any], create: bool = False) -> Dict[str, Any]:
    """Return the ``operation`` namespace dict, creating it when ``create`` is True."""
    if not isinstance(meta, dict):
        return {}
    if "operation" in meta and isinstance(meta["operation"], dict):
        return meta["operation"]
    if not create:
        return {}
    root = _ensure_root(meta)
    return _ensure_namespace(root, "operation")


def set_operation(meta: Dict[str, Any], *,
                  kind: str,
                  target_filter: Optional[Iterable[str]] = None,
                  priority: Optional[float] = None,
                  through: Optional[bool] = None,
                  consume: Optional[bool] = None,
                  policy: Optional[str] = None,
                  feature_id: Optional[str] = None,
                  feature_kind: Optional[str] = None,
                  stage: Optional[str] = None) -> Dict[str, Any]:
    """Mark the part as an operation (cutter/intersector/unioner). Returns the namespace."""
    _validate_enum(kind, _OPERATION_KINDS, field="operation.kind")
    _validate_enum(policy, _OPERATION_POLICIES, field="operation.policy")
    _validate_enum(feature_kind, _FEATURE_KINDS, field="operation.feature_kind")
    section = get_operation_metadata(meta, create=True)
    section["kind"] = kind
    if target_filter is not None:
        section["target_filter"] = [str(t) for t in target_filter]
    if priority is not None:
        section["priority"] = float(priority)
    if through is not None:
        section["through"] = bool(through)
    if consume is not None:
        section["consume"] = bool(consume)
    if policy is not None:
        section["policy"] = policy
    if feature_id is not None:
        section["feature_id"] = str(feature_id)
    if feature_kind is not None:
        section["feature_kind"] = feature_kind
    if stage is not None:
        section["stage"] = str(stage)
    return section


__all__ = [
    # entity-level helpers
    'get_surface_metadata',
    'set_surface_metadata',
    'ensure_surface_id',
    'set_surface_units',
    'set_surface_origin',
    'get_solid_metadata',
    'set_solid_metadata',
    'ensure_solid_id',
    'set_solid_context',
    'set_layer',
    # v1.0 namespace helpers
    'add_tags',
    'set_material',
    'set_manufacturing',
    'add_design_history_entry',
    'set_mass_constraint',
    'set_envelope_constraint',
    'add_compliance',
    'add_analysis_record',
    # v1.1 assembly namespace
    'get_assembly_metadata',
    'set_assembly',
    'add_bolt_pattern',
    'add_datum',
    'add_surface',
    'add_keepout',
    # v1.1 operation namespace
    'get_operation_metadata',
    'set_operation',
]
