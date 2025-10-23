"""Metadata helpers for yapCAD geometry objects."""

from __future__ import annotations

import uuid
from typing import Dict, Tuple, Any, Optional, Iterable

from yapcad.geom3d import issolid, issurface

_SURFACE_META_INDEX = 6
_SOLID_META_INDEX = 4
_DEFAULT_SCHEMA = "metadata-namespace-v0.1"

_ROOT_FIELDS = ("schema", "entityId", "timestamp", "tags", "layer")
_SECTION_KEYS = {
    "material",
    "manufacturing",
    "designHistory",
    "constraints",
    "analysis",
}


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
    if not meta:
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


__all__ = [
    'get_surface_metadata',
    'set_surface_metadata',
    'ensure_surface_id',
    'set_surface_units',
    'set_surface_origin',
    'get_solid_metadata',
    'set_solid_metadata',
    'ensure_solid_id',
    'set_solid_context',
    'add_tags',
    'set_material',
    'set_manufacturing',
    'add_design_history_entry',
    'set_mass_constraint',
    'set_envelope_constraint',
    'add_compliance',
    'add_analysis_record',
    'set_layer',
]


__all__ = [
    'get_surface_metadata',
    'set_surface_metadata',
    'ensure_surface_id',
    'set_surface_units',
    'set_surface_origin',
    'get_solid_metadata',
    'set_solid_metadata',
    'ensure_solid_id',
    'set_solid_context',
]
