# yapCAD Geometry JSON Schema

**Schema ID:** `yapcad-geometry-json-v0.1`  
**Status:** Accepted working version  
**Purpose:** Serialise yapCAD solids, surfaces, assemblies, and associated metadata into a portable JSON document for storage, interchange, or inclusion in `.ycpkg` packages.

---

## 1. Document Structure

```json5
{
  "schema": "yapcad-geometry-json-v0.1",
  "generator": {
    "name": "yapCAD",
    "version": "0.5.1",
    "build": "sha256:…"
  },
  "units": "mm",
  "entities": [
    { ... },   // solids, surfaces, meshes, groups
  ],
  "relationships": [
    { ... }    // optional assembly graph
  ],
  "attachments": [
    { ... }    // external artefacts (STEP/STL manifests)
  ]
}
```

- `schema` (string, required): Version identifier.
- `generator` (object, optional): Tool metadata (name/version/build hash).
- `units` (string, optional): Default length unit (e.g. `"mm"`, `"inch"`).
- `entities` (array, required): Geometry objects described below.
- `relationships` (array, optional): Parent/child links, assembly constraints.
- `attachments` (array, optional): Non-native files referenced by hash/path.

---

## 2. Entity Definitions

Each entity is a JSON object with the following common fields:

| field          | type     | description |
|----------------|----------|-------------|
| `id`           | string   | Stable UUID (should match `entityId` in metadata). |
| `type`         | string   | `"solid"`, `"surface"`, `"mesh"`, `"group"`, `"datum"`, `"sketch"`. |
| `name`         | string   | Optional human-readable label. |
| `metadata`     | object   | Metadata dictionary conforming to `metadata_namespace.md`. |
| `boundingBox`  | number[] | `[xmin, ymin, zmin, xmax, ymax, zmax]` in document units. |
| `properties`   | object   | Derived properties (volume, area). |

Each entity's `metadata` block MUST include the root fields from `metadata_namespace.md`, including `layer` (defaulting to `"default"`). Serialisers SHOULD propagate layer assignments for solids down to child surfaces and sketches so that viewers can offer layer-based visibility controls.

### 2.1 Solids (`type: "solid"`)

```json5
{
  "type": "solid",
  "faces": [
    {
      "surface": "<surface-id>",
      "orientation": 1,         // +1 outer, -1 reversed
      "edges": ["<edge-id>", ...]
    }
  ],
  "shell": ["<surface-id>", ...],  // ordered to define closed surface set
  "voids": [
    ["<surface-id>", ...]
  ]
}
```

### 2.2 Surfaces (`type: "surface"`)

```json5
{
  "type": "surface",
  "vertices": [[x, y, z, 1], ...],
  "normals":  [[nx, ny, nz, 0], ...],
  "faces":    [[i0, i1, i2], ...],   // index into `vertices`
  "triangulation": {
    "winding": "ccw",
    "topology": "triangle"
  }
}
```

### 2.3 Meshes (`type: "mesh"`)

- For imported STL/PLY assets. Same structure as surface but without homogeneous coordinates requirement (`[x, y, z]`).

### 2.4 Groups/Assemblies (`type: "group"`)

```json5
{
  "type": "group",
  "children": ["<entity-id>", ...],
  "transform": [
    [1,0,0,dx],
    [0,1,0,dy],
    [0,0,1,dz],
    [0,0,0,1]
  ]
}
```

---

## 3. Relationships

Relationships capture connections not expressible as simple hierarchy:

```json5
{
  "type": "mated",
  "entities": ["<solid-A>", "<solid-B>"],
  "constraint": {
    "kind": "coincident",
    "faceA": "<face-id>",
    "faceB": "<face-id>",
    "offset": 0.0
  }
}
```

Supported relationship types (initial set):

- `mated`
- `derived-from`
- `superseded-by`
- `references-layer`

---

## 4. Attachments

Attachment entries register external artefacts alongside hashes for integrity.

```json5
{
  "id": "step-export-001",
  "kind": "step",
  "path": "exports/rocket.step",
  "hash": "sha256:…",
  "createdBy": "<entity-id>",
  "metadata": {
    "brep": true,
    "tessellated": true
  }
}
```

---

## 5. Serialization Rules

1. Numeric values MUST be finite doubles. Serialisers must replace `±inf`/`NaN` with errors.
2. Homogeneous coordinates still use yapCAD convention (`w=1` for points, `w=0` for vectors).
3. Vertex indices are zero-based.
4. The root document MAY contain multiple solids or groups; exporters should topologically sort for dependency-free reconstruction.
5. Unknown fields MUST be preserved on round-trip (forward compatibility).

---

## 6. Integration Guidance

- **Geometry export**: implement `to_geometry_json(entity)` that walks solids, surfaces, metadata and writes this schema.
- **Import**: validate `schema` version, then rebuild yapCAD list structures; reattach metadata via helper functions.
- **Manifest use**: `.ycpkg` manifests can reference geometry JSON by path and include matching hashes.
- **Streaming**: allow chunked outputs by splitting `entities` across files and referencing them via `attachments` or manifest entries.

---

## 7. Future Enhancements

- Formal JSON Schema / Pydantic definitions.
- Compression guidelines (e.g. `.json.zst`).
- Support for parametric feature history (link to DSL once available).
- Extend relationships to cover constraint solving results and tolerance stacks.
