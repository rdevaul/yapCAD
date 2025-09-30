"""Metadata helpers for yapCAD geometry objects."""

from __future__ import annotations

import uuid
from typing import Dict, Tuple

from yapcad.geom3d import issolid, issurface

_SURFACE_META_INDEX = 6
_SOLID_META_INDEX = 4


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
    return _ensure_metadata(surface, _SURFACE_META_INDEX, create)


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
    meta = get_surface_metadata(surface, create=True)
    if 'id' not in meta:
        meta['id'] = str(uuid.uuid4())
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
    return _ensure_metadata(sld, _SOLID_META_INDEX, create)


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
    meta = get_solid_metadata(sld, create=True)
    if 'id' not in meta:
        meta['id'] = str(uuid.uuid4())
    return meta['id']


def set_solid_context(sld: list, context: Dict) -> list:
    if not isinstance(context, dict):
        raise ValueError('context must be a dictionary')
    meta = get_solid_metadata(sld, create=True)
    meta['context'] = context
    return sld


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
