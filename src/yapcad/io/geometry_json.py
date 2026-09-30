"""Geometry JSON serialization/deserialization helpers.

Implements the draft schema described in ``docs/geometry_json_schema.md``.
"""

from __future__ import annotations

import base64
from copy import deepcopy
import hashlib
import math
import uuid
from typing import Any, Dict, Iterable, List, Optional, Tuple, Sequence

from yapcad.geom import (
    arc,
    catmullrom,
    ellipse,
    isnurbs,
    iscatmullrom,
    isellipse,
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
from yapcad.geom3d import issolid, issurface, solid_shells, solidbbox, surfacebbox
from yapcad.brep import brep_from_solid, occ_available
from yapcad.construction import construction_from_json, construction_to_json
from yapcad.metadata import (
    get_solid_metadata,
    get_surface_metadata,
    ensure_solid_id,
    ensure_surface_id,
    set_solid_metadata,
    set_surface_metadata,
    set_layer,
)
LEGACY_SCHEMA_ID = "yapcad-geometry-json-v0.1"
V0_2_SCHEMA_ID = "yapcad-geometry-json-v0.2"
#: Current schema. Bumped from v0.2 for SDF: extending the ``authoritative``
#: enum is a breaking change for a strict reader, so the id has to move with
#: it (design document §8.2).
SCHEMA_ID = "yapcad-geometry-json-v0.3"
SUPPORTED_SCHEMA_IDS = frozenset({LEGACY_SCHEMA_ID, V0_2_SCHEMA_ID, SCHEMA_ID})

#: Schemas that carry an explicit ``representations`` block on every solid.
REPRESENTATION_SCHEMA_IDS = frozenset({V0_2_SCHEMA_ID, SCHEMA_ID})

#: Serialised SDF tree format, mirrored from :mod:`yapcad.sdf.node` so that
#: reading a document does not require importing the SDF subsystem.
SDF_TREE_FORMAT = "yapcad-sdf-tree-v1"


def _sdf_node_from_tree(tree: Any):
    """Parse and validate a serialised SDF tree.

    Imported lazily: the SDF subsystem pulls in numpy and geom3d, and a
    document with no SDF content should not pay for that.
    """
    from yapcad.sdf.node import tree_from_json
    return tree_from_json(tree)


def _sdf_representation(
    construction: Optional[List[Any]],
) -> Optional[Dict[str, Any]]:
    """Build the ``sdf`` representation record from a construction record.

    Returns ``None`` when the solid was not authored as a field.
    """
    if not isinstance(construction, list) or len(construction) < 2:
        return None
    if construction[0] != "sdf":
        return None
    tree = construction[1]
    node = _sdf_node_from_tree(tree)
    record: Dict[str, Any] = {
        "role": "authoritative",
        "format": SDF_TREE_FORMAT,
        # Whether the field is the true distance, or merely bounds it.
        "lipschitz": "exact" if node.exact else "bound",
        "lipschitzConstant": float(node.lipschitz),
        # Structural: can this tree be replayed through OCC for an exact
        # BREP? Surfaced so a consumer can tell before attempting STEP.
        "csgExact": bool(node.csg_exact),
        "tree": tree,
    }
    lo, hi = node.bounds
    extent = [float(v) for v in (*lo, *hi)]
    if all(math.isfinite(v) for v in extent):
        record["bounds"] = extent
    return record


def _sdf_mesh_parameters(construction: List[Any]) -> Dict[str, Any]:
    """Reproducibility parameters for a mesh generated from a field."""
    params: Dict[str, Any] = {"generatedFrom": "sdf"}
    if len(construction) >= 3 and isinstance(construction[2], dict):
        for key in ("method", "resolution", "padding", "bounds"):
            if key in construction[2]:
                params[key] = deepcopy(construction[2][key])
    return params


def _construction_from_sdf(record: Any, mesh_record: Dict[str, Any],
                           solid_id: str) -> List[Any]:
    """Rebuild the ``['sdf', tree, meshing]`` construction record on read."""
    if not isinstance(record, dict):
        raise ValueError(
            f"solid {solid_id} SDF representation must be an object"
        )
    if record.get("role") != "authoritative":
        raise ValueError(
            f"solid {solid_id} SDF representation role must be 'authoritative'"
        )
    if record.get("format") != SDF_TREE_FORMAT:
        raise ValueError(
            f"solid {solid_id} has unsupported SDF format "
            f"{record.get('format')!r}"
        )
    tree = record.get("tree")
    node = _sdf_node_from_tree(tree)

    # The analysis fields are a cache of what the tree already determines.
    # Verify rather than trust, the same way the BREP payload hash is
    # checked: a document claiming csgExact for a gyroid must be rejected,
    # not believed.
    if bool(record.get("csgExact")) != node.csg_exact:
        raise ValueError(
            f"solid {solid_id} SDF csgExact disagrees with its tree"
        )
    expected_kind = "exact" if node.exact else "bound"
    if record.get("lipschitz") != expected_kind:
        raise ValueError(
            f"solid {solid_id} SDF lipschitz kind disagrees with its tree"
        )
    claimed = record.get("lipschitzConstant")
    if claimed is None or not math.isclose(
        float(claimed), node.lipschitz, rel_tol=1e-9, abs_tol=1e-12
    ):
        raise ValueError(
            f"solid {solid_id} SDF lipschitzConstant disagrees with its tree"
        )

    construction: List[Any] = ["sdf", tree]
    meshing = {
        key: deepcopy(mesh_record[key])
        for key in ("method", "resolution", "padding", "bounds")
        if key in mesh_record
    }
    if meshing:
        construction.append(meshing)
    return construction


def _kernel_record(legacy_brep: Dict[str, Any]) -> Dict[str, str]:
    kernel = legacy_brep.get("kernel")
    if isinstance(kernel, dict) and kernel.get("name"):
        if kernel["name"] != "OpenCASCADE":
            raise ValueError(
                f"unsupported BREP kernel: {kernel['name']!r}"
            )
        return {
            "name": "OpenCASCADE",
            "version": str(kernel.get("version", "unknown")),
        }
    version = "unknown"
    try:  # pragma: no cover - version availability varies by OCC packaging
        import OCC
        version = str(getattr(OCC, "VERSION", version))
    except ImportError:
        pass
    return {"name": "OpenCASCADE", "version": version}


def _explicit_brep_record(legacy_brep: Dict[str, Any],
                          role: str = "authoritative") -> Dict[str, Any]:
    if legacy_brep.get("encoding") != "brep-ascii-base64":
        raise ValueError(
            f"unsupported legacy BREP encoding: {legacy_brep.get('encoding')!r}"
        )
    encoded = legacy_brep.get("data")
    if not isinstance(encoded, str) or not encoded:
        raise ValueError("BREP payload must be a non-empty base64 string")
    try:
        payload = base64.b64decode(encoded, validate=True)
    except Exception as exc:
        raise ValueError("BREP payload is not valid base64") from exc
    record: Dict[str, Any] = {
        "role": role,
        "format": "opencascade-brep",
        "encoding": "base64",
        "payload": encoded,
        "hash": "sha256:" + hashlib.sha256(payload).hexdigest(),
        "kernel": _kernel_record(legacy_brep),
    }
    if legacy_brep.get("modelTolerance") is not None:
        tolerance = float(legacy_brep["modelTolerance"])
        if not math.isfinite(tolerance) or tolerance <= 0:
            raise ValueError("BREP modelTolerance must be finite and positive")
        record["modelTolerance"] = tolerance
    if role == "derived":
        # Which tree this was replayed from, so a reader can tell a derived
        # BREP from a stale or unrelated one without replaying it.
        record["derivedFrom"] = "sdf"
        record["treeDigest"] = legacy_brep["derivedFrom"]["sdf"]
    return record


def _legacy_brep_record(explicit: Dict[str, Any],
                        role: str = "authoritative") -> Dict[str, Any]:
    if not isinstance(explicit, dict):
        raise ValueError("BREP representation must be an object")
    if explicit.get("format") != "opencascade-brep":
        raise ValueError(f"unsupported BREP format: {explicit.get('format')!r}")
    if explicit.get("encoding") != "base64":
        raise ValueError(f"unsupported BREP encoding: {explicit.get('encoding')!r}")
    if explicit.get("role") != role:
        raise ValueError(f"BREP representation role must be {role!r}")
    kernel = explicit.get("kernel")
    if not isinstance(kernel, dict) or kernel.get("name") != "OpenCASCADE":
        raise ValueError("BREP representation requires an OpenCASCADE kernel record")
    encoded = explicit.get("payload")
    if not isinstance(encoded, str) or not encoded:
        raise ValueError("BREP payload must be a non-empty base64 string")
    try:
        payload = base64.b64decode(encoded, validate=True)
    except Exception as exc:
        raise ValueError("BREP payload is not valid base64") from exc
    expected = explicit.get("hash")
    actual = "sha256:" + hashlib.sha256(payload).hexdigest()
    if expected != actual:
        raise ValueError(f"BREP payload hash mismatch: expected {expected!r}, got {actual}")
    legacy: Dict[str, Any] = {
        "encoding": "brep-ascii-base64",
        "data": encoded,
        "hash": actual,
        "kernel": deepcopy(kernel),
    }
    if explicit.get("modelTolerance") is not None:
        tolerance = float(explicit["modelTolerance"])
        if not math.isfinite(tolerance) or tolerance <= 0:
            raise ValueError("BREP modelTolerance must be finite and positive")
        legacy["modelTolerance"] = tolerance
    if role == "derived":
        if explicit.get("derivedFrom") != "sdf":
            raise ValueError("a derived BREP must declare derivedFrom 'sdf'")
        tree_digest = explicit.get("treeDigest")
        if not isinstance(tree_digest, str) or not tree_digest:
            raise ValueError("a derived BREP must name the treeDigest it was "
                             "replayed from")
        legacy["derivedFrom"] = {"sdf": tree_digest}
    return legacy


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

    layers_seen = []
    parent_layer = metadata.get("layer")

    def _emit_surface(surf) -> Optional[str]:
        """Serialize one surface into the cache and return its id."""
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
        layers_seen.append(serialized["metadata"].get("layer", "default"))
        return surface_id

    # yapCAD mesh solids carry cavities as inward-wound shells inside the
    # surface list rather than in a separate slot, so the outer/void split is
    # derived here (see geom3d.solid_shells).  The partition is only emitted
    # when it is unambiguous -- a cavity was actually found -- so documents
    # for ordinary solids are unchanged.
    outer_surfaces = list(solid[1])
    void_groups: List[List[list]] = []
    try:
        derived_outer, derived_voids = solid_shells(solid)
    except Exception:
        derived_outer, derived_voids = None, None
    if derived_voids:
        outer_surfaces = [surf for shell in derived_outer for surf in shell]
        void_groups = derived_voids

    shell_ids: List[str] = []
    for surf in outer_surfaces:
        if not issurface(surf):
            continue
        shell_ids.append(_emit_surface(surf))

    voids: List[List[str]] = []
    for group in void_groups:
        void_ids = [_emit_surface(surf) for surf in group if issurface(surf)]
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

    metadata_record = deepcopy(metadata)
    legacy_brep = metadata_record.pop("brep", None)

    # Provenance from the solid's construction slot. An SDF record is not
    # merely provenance: it is the authoritative definition of the solid,
    # and is promoted into the representations block below.
    construction = construction_to_json(solid)
    sdf_record = _sdf_representation(construction)

    brep_role = "authoritative"
    if sdf_record is not None and legacy_brep:
        # An SDF solid may carry a BREP only as a *derived* representation:
        # one replayed from this very tree.  Anything else -- an analytic
        # BREP from prism(), or one replayed from an earlier version of the
        # tree -- would be a second, disagreeing definition of the solid.
        marker = (legacy_brep.get("derivedFrom") or {}).get("sdf")
        if marker != sdf_record["tree"]["root"]:
            raise ValueError(
                f"solid {solid_id} claims both SDF and BREP authority; "
                "authority is a single per-solid property. Its BREP was not "
                "replayed from its SDF tree -- regenerate it with "
                "yapcad.sdf.to_solid(..., brep=True), or drop one of the two."
            )
        brep_role = "derived"

    if sdf_record is not None:
        authoritative = "sdf"
    elif legacy_brep:
        authoritative = "brep"
    else:
        authoritative = "mesh"

    mesh_record: Dict[str, Any] = {
        "format": "indexed-triangle-set",
        "role": "authoritative" if authoritative == "mesh" else "preview",
    }
    if sdf_record is not None:
        mesh_record.update(_sdf_mesh_parameters(construction))

    representations: Dict[str, Any] = {
        "authoritative": authoritative,
        "mesh": mesh_record,
    }
    if legacy_brep:
        representations["brep"] = _explicit_brep_record(legacy_brep,
                                                        role=brep_role)
    if sdf_record is not None:
        representations["sdf"] = sdf_record

    entry = {
        "id": metadata.get("entityId", solid_id),
        "type": "solid",
        "name": metadata.get("name"),
        "metadata": metadata_record,
        "boundingBox": bbox,
        "properties": {},
        "shell": shell_ids,
        "voids": voids,
        "representations": representations,
    }

    # Omitted entirely when the solid records no construction, so documents
    # for such solids are unchanged. Also omitted for an SDF solid: the tree
    # is already in representations.sdf.tree, and a second copy here would be
    # two encodings of one fact that can drift apart.
    if construction and sdf_record is None:
        entry["construction"] = construction

    return entry


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
        elif isellipse(element):
            center = _point_components(element[1])
            meta = element[2]
            primitive = {
                "kind": "ellipse",
                "center": center,
                "semi_major": float(meta['semi_major']),
                "semi_minor": float(meta['semi_minor']),
                "rotation": float(meta['rotation']),
                "start": float(meta['start']) if meta['start'] != 0 else 0,
                "end": float(meta['end']) if meta['end'] != 360 else 360,
                "normal": _float_vec(meta['normal']),
            }
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
        "schema": "metadata-namespace-v1.1",
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
    schema = doc.get("schema")
    if schema not in SUPPORTED_SCHEMA_IDS:
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

        # A yapCAD mesh solid keeps cavities in its surface list as inward-wound
        # shells; it has no voids slot, and slot 2 is material.  Void surfaces
        # from the document therefore join the shell surfaces, and the split is
        # recovered on demand via geom3d.solid_shells.
        for void_ids in entry.get("voids", []):
            for sid in void_ids:
                surface = surfaces.get(sid)
                if surface is None:
                    surf_entry = entries_by_id.get(sid)
                    if surf_entry is None:
                        raise ValueError(f"void surface {sid} missing")
                    surface = _rehydrate_surface(surf_entry)
                    surfaces[sid] = surface
                shell_surfaces.append(surface)

        try:
            construction = construction_from_json(entry.get("construction"))
        except ValueError as exc:
            raise ValueError(
                f"solid {entry['id']} has a malformed construction record: {exc}"
            ) from exc

        solid = ['solid', shell_surfaces, [], construction]
        metadata = deepcopy(entry.get("metadata") or {})
        if schema in REPRESENTATION_SCHEMA_IDS:
            representations = entry.get("representations")
            if not isinstance(representations, dict):
                raise ValueError(f"solid {entry['id']} missing representations")
            authoritative = representations.get("authoritative")
            # "sdf" joined the enum in v0.3; a v0.2 document claiming it is
            # malformed, which is exactly why the schema id had to bump.
            allowed = {"brep", "mesh"}
            if schema == SCHEMA_ID:
                allowed.add("sdf")
            if authoritative not in allowed:
                raise ValueError(
                    f"solid {entry['id']} has invalid authoritative representation"
                )
            mesh_record = representations.get("mesh")
            if not isinstance(mesh_record, dict):
                raise ValueError(f"solid {entry['id']} missing mesh representation")
            if mesh_record.get("format") != "indexed-triangle-set":
                raise ValueError(f"solid {entry['id']} has invalid mesh format")
            expected_mesh_role = (
                "authoritative" if authoritative == "mesh" else "preview"
            )
            if mesh_record.get("role") != expected_mesh_role:
                raise ValueError(
                    f"solid {entry['id']} mesh role must be {expected_mesh_role!r}"
                )
            if authoritative == "brep":
                if "brep" not in representations:
                    raise ValueError(f"solid {entry['id']} missing BREP representation")
                metadata["brep"] = _legacy_brep_record(representations["brep"])
            elif "brep" in representations and authoritative != "sdf":
                raise ValueError(
                    f"solid {entry['id']} contains non-authoritative BREP representation"
                )
            if authoritative == "sdf":
                if "sdf" not in representations:
                    raise ValueError(
                        f"solid {entry['id']} missing SDF representation"
                    )
                if entry.get("construction") is not None:
                    raise ValueError(
                        f"solid {entry['id']} duplicates its SDF tree in the "
                        "construction field"
                    )
                solid[3] = _construction_from_sdf(
                    representations["sdf"], mesh_record, entry["id"]
                )
                if "brep" in representations:
                    derived = _legacy_brep_record(representations["brep"],
                                                  role="derived")
                    # Compare against the digest of the tree as rebuilt, not
                    # the id the document gives it: a table key is only a
                    # label, and the check is meant to catch drift.
                    rebuilt = _sdf_node_from_tree(solid[3][1]).digest
                    if derived["derivedFrom"]["sdf"] != rebuilt:
                        raise ValueError(
                            f"solid {entry['id']} carries a derived BREP that "
                            "was replayed from a different SDF tree"
                        )
                    metadata["brep"] = derived
            elif "sdf" in representations:
                raise ValueError(
                    f"solid {entry['id']} contains a non-authoritative "
                    "SDF representation"
                )
        if metadata:
            set_solid_metadata(solid, metadata)
        modern = schema in REPRESENTATION_SCHEMA_IDS
        restored_brep = brep_from_solid(solid, refresh=modern)
        if modern and metadata.get("brep") and occ_available():
            if restored_brep is None:
                raise ValueError(f"solid {entry['id']} BREP payload could not be loaded")
            from OCC.Core.BRepCheck import BRepCheck_Analyzer
            if not BRepCheck_Analyzer(restored_brep.shape).IsValid():
                raise ValueError(f"solid {entry['id']} BREP topology is invalid")
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
                elif kind == "ellipse":
                    center = _point_from_components(prim.get("center", [0.0, 0.0, 0.0, 1.0]))
                    semi_major = float(prim.get("semi_major", 1.0))
                    semi_minor = float(prim.get("semi_minor", 1.0))
                    rotation = float(prim.get("rotation", 0.0))
                    start = prim.get("start", 0)
                    end = prim.get("end", 360)
                    normal = prim.get("normal")
                    elements.append(
                        ellipse(
                            center,
                            semi_major,
                            semi_minor,
                            rotation=rotation,
                            start=start,
                            end=end,
                            normal=normal,
                        )
                    )
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
    "LEGACY_SCHEMA_ID",
    "V0_2_SCHEMA_ID",
    "SUPPORTED_SCHEMA_IDS",
    "REPRESENTATION_SCHEMA_IDS",
    "geometry_to_json",
    "geometry_from_json",
]
