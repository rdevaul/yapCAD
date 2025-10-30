"""Geometry JSON serialization/deserialization helpers.

Implements the draft schema described in ``docs/geometry_json_schema.md``.
"""

from __future__ import annotations

import uuid
from typing import Any, Dict, Iterable, List, Optional, Tuple, Sequence

from yapcad.geom import (
    arc,
    catmullrom,
    isnurbs,
    iscatmullrom,
    nurbs,
    point,
    vect,
    bbox as bbox2d,
    isarc,
    iscircle,
    isgeomlist,
    isline,
    length,
    line,
    sample,
)
from yapcad.geom3d import issolid, issurface, solidbbox, surfacebbox
from yapcad.metadata import (
    get_solid_metadata,
    get_surface_metadata,
    ensure_solid_id,
    ensure_surface_id,
    set_solid_metadata,
    set_surface_metadata,
    set_layer,
)
SCHEMA_ID = "yapcad-geometry-json-v0.1"


def _float_vec(vec: Iterable[float]) -> List[float]:
    return [float(c) for c in vec]


def _int_vec(vec: Iterable[int]) -> List[int]:
    return [int(c) for c in vec]


def _bbox_or_none(box: Optional[List[List[float]]]) -> Optional[List[float]]:
    if not box:
        return None
    (xmin, ymin, zmin, _), (xmax, ymax, zmax, _) = box
    return [float(xmin), float(ymin), float(zmin), float(xmax), float(ymax), float(zmax)]


def _point_components(pt: Sequence[float]) -> List[float]:
    """Return point components including homogeneous coordinate."""
    x = float(pt[0])
    y = float(pt[1])
    z = float(pt[2]) if len(pt) > 2 else 0.0
    w = float(pt[3]) if len(pt) > 3 else 1.0
    return [x, y, z, w]


def _point_from_components(components: Sequence[float]) -> List[float]:
    if len(components) >= 4:
        return point(float(components[0]), float(components[1]), float(components[2]), float(components[3]))
    if len(components) == 3:
        return point(float(components[0]), float(components[1]), float(components[2]))
    return point(float(components[0]), float(components[1]))


def _serialize_surface(surface: list, metadata_override: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    try:
        surface_id = ensure_surface_id(surface)
    except Exception as exc:
        raise RuntimeError(f"failed to ensure surface metadata for {surface!r}") from exc
    metadata = get_surface_metadata(surface, create=True)
    if metadata_override:
        metadata.update(metadata_override)
    if "layer" not in metadata or not metadata.get("layer"):
        metadata["layer"] = "default"
    verts = surface[1]
    norms = surface[2]
    faces = surface[3]
    try:
        bbox = _bbox_or_none(surfacebbox(surface))
    except Exception:
        bbox = None
    return {
        "id": metadata.get("entityId", surface_id),
        "type": "surface",
        "name": metadata.get("name"),
        "metadata": metadata,
        "boundingBox": bbox,
        "properties": {},
        "vertices": [_float_vec(v) for v in verts],
        "normals": [_float_vec(n) for n in norms],
        "faces": [_int_vec(face) for face in faces],
        "triangulation": {
            "winding": "ccw",
            "topology": "triangle",
        },
    }


def _serialize_solid(solid: list, surface_cache: Dict[str, Dict[str, Any]], metadata_override: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    solid_id = ensure_solid_id(solid)
    metadata = get_solid_metadata(solid, create=True)
    if metadata_override:
        metadata.update(metadata_override)
    if "layer" not in metadata or not metadata.get("layer"):
        metadata["layer"] = "default"
    try:
        bbox = _bbox_or_none(solidbbox(solid))
    except Exception:
        bbox = None

    shell_ids: List[str] = []
    layers_seen = []
    parent_layer = metadata.get("layer")
    for surf in solid[1]:
        if not issurface(surf):
            continue
        surface_meta = get_surface_metadata(surf, create=True)
        surface_override = None
        if metadata_override:
            surface_override = dict(metadata_override)
        if parent_layer and (not surface_meta.get("layer") or surface_meta.get("layer") == "default"):
            surface_override = dict(surface_override or {})
            surface_override["layer"] = parent_layer
        serialized = _serialize_surface(surf, surface_override)
        surface_id = serialized["id"]
        while surface_id in surface_cache:
            new_id = uuid.uuid4().hex
            meta = serialized["metadata"]
            meta["entityId"] = new_id
            meta["id"] = new_id
            serialized["id"] = new_id
            set_surface_metadata(surf, meta)
            surface_id = new_id
        surface_cache[surface_id] = serialized
        shell_ids.append(surface_id)
        layers_seen.append(serialized["metadata"].get("layer", "default"))

    voids: List[List[str]] = []
    if len(solid) > 2:
        for void in solid[2] or []:
            void_ids: List[str] = []
            for surf in void:
                if not issurface(surf):
                    continue
                surface_meta = get_surface_metadata(surf, create=True)
                surface_override = None
                if metadata_override:
                    surface_override = dict(metadata_override)
                if parent_layer and (not surface_meta.get("layer") or surface_meta.get("layer") == "default"):
                    surface_override = dict(surface_override or {})
                    surface_override["layer"] = parent_layer
                serialized = _serialize_surface(surf, surface_override)
                surface_id = serialized["id"]
                while surface_id in surface_cache:
                    new_id = uuid.uuid4().hex
                    meta = serialized["metadata"]
                    meta["entityId"] = new_id
                    meta["id"] = new_id
                    serialized["id"] = new_id
                    set_surface_metadata(surf, meta)
                    surface_id = new_id
                surface_cache[surface_id] = serialized
                void_ids.append(surface_id)
                layers_seen.append(serialized["metadata"].get("layer", "default"))
            voids.append(void_ids)

    if "layer" not in metadata or not metadata.get("layer"):
        unique_layers = [layer for layer in layers_seen if layer]
        if unique_layers:
            first_layer = unique_layers[0]
            if all(layer == first_layer for layer in unique_layers):
                metadata["layer"] = first_layer
            else:
                metadata["layer"] = "default"
        else:
            metadata["layer"] = "default"

    return {
        "id": metadata.get("entityId", solid_id),
        "type": "solid",
        "name": metadata.get("name"),
        "metadata": metadata,
        "boundingBox": bbox,
        "properties": {},
        "shell": shell_ids,
        "voids": voids,
    }


def _polyline_points(sequence: List[float]) -> List[float]:
    return [float(sequence[0]), float(sequence[1])]


def _sample_geometry_element(element: list, min_segments: int = 8) -> List[List[float]]:
    if isline(element):
        return [
            _polyline_points(element[0]),
            _polyline_points(element[1]),
        ]
    if iscatmullrom(element) or isnurbs(element):
        segs = max(min_segments, 32)
    else:
        segs = max(min_segments, int(max(length(element), 1.0)))
    points = []
    for i in range(segs + 1):
        t = i / segs
        p = sample(element, t)
        points.append(_polyline_points(p))
    return points


def _serialize_sketch(geomlist: list, metadata_override: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    poly_vectors: List[List[List[float]]] = []
    primitives: List[Dict[str, Any]] = []

    for element in geomlist:
        if isline(element):
            start = _polyline_points(element[0])
            end = _polyline_points(element[1])
            poly_vectors.append([start, end])
            primitives.append(
                {
                    "kind": "line",
                    "start": start,
                    "end": end,
                }
            )
        elif iscircle(element):
            center = _point_components(element[0])
            radius = float(element[1][0])
            primitive: Dict[str, Any] = {
                "kind": "circle",
                "center": center,
                "radius": radius,
                "orientation": int(element[1][3]),
            }
            if len(element) >= 3:
                primitive["normal"] = _point_components(element[2])
            primitives.append(primitive)
            poly_vectors.append(_sample_geometry_element(element))
        elif isarc(element):
            center = _point_components(element[0])
            radius = float(element[1][0])
            start_angle = float(element[1][1])
            end_angle = float(element[1][2])
            primitive = {
                "kind": "arc",
                "center": center,
                "radius": radius,
                "start": start_angle,
                "end": end_angle,
                "orientation": int(element[1][3]),
            }
            if len(element) >= 3:
                primitive["normal"] = _point_components(element[2])
            primitives.append(primitive)
            poly_vectors.append(_sample_geometry_element(element))
        elif iscatmullrom(element):
            control_points = [_point_components(pt) for pt in element[1]]
            params = dict(element[2])
            primitives.append(
                {
                    "kind": "catmullrom",
                    "points": control_points,
                    "params": params,
                }
            )
            poly_vectors.append(_sample_geometry_element(element))
        elif isnurbs(element):
            control_points = [_point_components(pt) for pt in element[1]]
            params = dict(element[2])
            primitives.append(
                {
                    "kind": "nurbs",
                    "points": control_points,
                    "params": params,
                }
            )
            poly_vectors.append(_sample_geometry_element(element))
        elif isinstance(element, list) and len(element) >= 2 and all(isinstance(pt, list) for pt in element):
            coords = [_polyline_points(pt) for pt in element]
            if coords:
                poly_vectors.append(coords)
                primitives.append({"kind": "polyline", "points": coords})

    box = bbox2d(geomlist)
    if box:
        min_pt, max_pt = box
        bbox = [float(min_pt[0]), float(min_pt[1]), float(max_pt[0]), float(max_pt[1])]
    else:
        bbox = None
    metadata = {
        "schema": "metadata-namespace-v0.1",
        "entityId": str(uuid.uuid4()),
        "tags": [],
        "layer": "default",
    }
    if metadata_override:
        metadata.update(metadata_override)
    return {
        "id": metadata["entityId"],
        "type": "sketch",
        "name": None,
        "metadata": metadata,
        "boundingBox": bbox,
        "polylines": poly_vectors,
        "primitives": primitives,
    }


def geometry_to_json(
    entities: Iterable[list],
    *,
    units: Optional[str] = None,
    generator: Optional[Dict[str, Any]] = None,
    relationships: Optional[List[Dict[str, Any]]] = None,
    attachments: Optional[List[Dict[str, Any]]] = None,
) -> Dict[str, Any]:
    """Serialize solids/surfaces into the geometry JSON document."""
    serialized_entities: List[Dict[str, Any]] = []
    surface_cache: Dict[str, Dict[str, Any]] = {}

    for item in entities:
        metadata_override = None
        entity = item
        if isinstance(item, dict) and 'geometry' in item:
            metadata_override = dict(item.get('metadata', {}))
            entity = item['geometry']
        elif isinstance(item, tuple) and len(item) == 2 and isinstance(item[1], dict):
            entity, metadata_override = item[0], dict(item[1])

        if issolid(entity):
            if metadata_override:
                meta = get_solid_metadata(entity, create=True)
                layer = metadata_override.get("layer")
                if layer:
                    set_layer(meta, layer)
                meta.update({k: v for k, v in metadata_override.items() if k != "layer"})
            solid_entry = _serialize_solid(entity, surface_cache, metadata_override)
            serialized_entities.append(solid_entry)
        elif issurface(entity):
            if metadata_override:
                meta = get_surface_metadata(entity, create=True)
                layer = metadata_override.get("layer")
                if layer:
                    set_layer(meta, layer)
                meta.update({k: v for k, v in metadata_override.items() if k != "layer"})
            surface_entry = _serialize_surface(entity, metadata_override)
            surface_cache.setdefault(surface_entry["id"], surface_entry)
        elif isgeomlist(entity):
            sketch_entry = _serialize_sketch(entity, metadata_override)
            serialized_entities.append(sketch_entry)
        else:
            raise ValueError("unsupported entity type for serialization")

    serialized_entities.extend(surface_cache.values())

    doc: Dict[str, Any] = {
        "schema": SCHEMA_ID,
        "entities": serialized_entities,
    }
    if units:
        doc["units"] = units
    if generator:
        doc["generator"] = generator
    if relationships:
        doc["relationships"] = list(relationships)
    if attachments:
        doc["attachments"] = list(attachments)
    return doc


def _rehydrate_surface(entry: Dict[str, Any]) -> list:
    verts = [list(map(float, v)) for v in entry.get("vertices", [])]
    norms = [list(map(float, n)) for n in entry.get("normals", [])]
    faces = [list(map(int, f)) for f in entry.get("faces", [])]
    surface = ['surface', verts, norms, faces, [], []]
    metadata = entry.get("metadata")
    if metadata:
        set_surface_metadata(surface, metadata)
    return surface


def geometry_from_json(doc: Dict[str, Any]) -> List[list]:
    """Deserialize geometry JSON into yapCAD list structures."""
    if doc.get("schema") != SCHEMA_ID:
        raise ValueError(f"unsupported geometry schema: {doc.get('schema')}")

    entries_by_id: Dict[str, Dict[str, Any]] = {}
    for entry in doc.get("entities", []):
        entry_id = entry.get("id")
        if not entry_id:
            raise ValueError("entity missing id")
        entries_by_id[entry_id] = entry

    surfaces: Dict[str, list] = {}
    solids: List[list] = []

    # First pass: instantiate surfaces explicitly listed
    for entry in doc.get("entities", []):
        if entry.get("type") == "surface":
            surfaces[entry["id"]] = _rehydrate_surface(entry)

    # Second pass: instantiate solids (and any referenced surfaces/voids)
    for entry in doc.get("entities", []):
        if entry.get("type") != "solid":
            continue
        shell_surfaces: List[list] = []
        for sid in entry.get("shell", []):
            surface = surfaces.get(sid)
            if surface is None:
                surf_entry = entries_by_id.get(sid)
                if surf_entry is None:
                    raise ValueError(f"surface {sid} referenced by solid but not provided")
                surface = _rehydrate_surface(surf_entry)
                surfaces[sid] = surface
            shell_surfaces.append(surface)

        voids: List[List[list]] = []
        for void_ids in entry.get("voids", []):
            void_surfaces: List[list] = []
            for sid in void_ids:
                surface = surfaces.get(sid)
                if surface is None:
                    surf_entry = entries_by_id.get(sid)
                    if surf_entry is None:
                        raise ValueError(f"void surface {sid} missing")
                    surface = _rehydrate_surface(surf_entry)
                    surfaces[sid] = surface
                void_surfaces.append(surface)
            voids.append(void_surfaces)

        solid = ['solid', shell_surfaces, voids, []]
        metadata = entry.get("metadata")
        if metadata:
            set_solid_metadata(solid, metadata)
        solids.append(solid)

    if solids:
        return solids
    sketches = [
        entry for entry in doc.get("entities", []) if entry.get("type") == "sketch"
    ]
    if sketches:
        geomlists: List[list] = []
        for sketch in sketches:
            primitives = sketch.get("primitives") or []
            elements: List[list] = []
            for prim in primitives:
                kind = prim.get("kind")
                if kind == "line":
                    start = prim.get("start", [])
                    end = prim.get("end", [])
                    if len(start) >= 2 and len(end) >= 2:
                        p0 = point(float(start[0]), float(start[1]))
                        p1 = point(float(end[0]), float(end[1]))
                        elements.append(line(p0, p1))
                elif kind == "circle":
                    center = _point_from_components(prim.get("center", [0.0, 0.0, 0.0, 1.0]))
                    radius = float(prim.get("radius", 0.0))
                    orientation = int(prim.get("orientation", -1))
                    vect_data = vect(radius, 0.0, 360.0, orientation)
                    normal = prim.get("normal")
                    if normal is not None:
                        elements.append(arc(center, vect_data, _point_from_components(normal)))
                    else:
                        elements.append(arc(center, vect_data))
                elif kind == "arc":
                    center = _point_from_components(prim.get("center", [0.0, 0.0, 0.0, 1.0]))
                    radius = float(prim.get("radius", 0.0))
                    start_angle = float(prim.get("start", 0.0))
                    end_angle = float(prim.get("end", 0.0))
                    orientation = int(prim.get("orientation", -1))
                    vect_data = vect(radius, start_angle, end_angle, orientation)
                    normal = prim.get("normal")
                    if normal is not None:
                        elements.append(arc(center, vect_data, _point_from_components(normal)))
                    else:
                        elements.append(arc(center, vect_data))
                elif kind == "polyline":
                    coords = prim.get("points", [])
                    pts = [point(float(pt[0]), float(pt[1])) for pt in coords if len(pt) >= 2]
                    for idx in range(len(pts) - 1):
                        elements.append(line(pts[idx], pts[idx + 1]))
                elif kind == "catmullrom":
                    ctrl_points = [_point_from_components(pt) for pt in prim.get("points", [])]
                    params = prim.get("params") or {}
                    elements.append(
                        catmullrom(
                            ctrl_points,
                            closed=bool(params.get("closed", False)),
                            alpha=float(params.get("alpha", 0.5)),
                        )
                    )
                elif kind == "nurbs":
                    ctrl_points = [_point_from_components(pt) for pt in prim.get("points", [])]
                    params = prim.get("params") or {}
                    degree = int(params.get("degree", 3))
                    weights = params.get("weights")
                    knots = params.get("knots")
                    elements.append(
                        nurbs(
                            ctrl_points,
                            degree=degree,
                            weights=weights,
                            knots=knots,
                        )
                    )
            if elements:
                geomlists.append(elements)
                continue

            polylines = []
            for poly in sketch.get("polylines", []):
                if len(poly) < 2:
                    continue
                pts = []
                for pt in poly:
                    x, y = pt
                    pts.append(point(x, y))
                segments = []
                for idx in range(len(pts) - 1):
                    segments.append(line(pts[idx], pts[idx + 1]))
                polylines.extend(segments)
            if polylines:
                geomlists.append(polylines)
        if geomlists:
            return geomlists
    return list(surfaces.values())


__all__ = [
    "SCHEMA_ID",
    "geometry_to_json",
    "geometry_from_json",
]
