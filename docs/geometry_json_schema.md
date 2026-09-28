# yapCAD Geometry JSON Schema

**Schema ID:** `yapcad-geometry-json-v0.3`
**Status:** Implemented draft
**Purpose:** Serialise yapCAD solids, surfaces, assemblies, and associated metadata into a portable JSON document for storage, interchange, or inclusion in `.ycpkg` packages.

---

## 1. Document Structure

```json5
{
  "schema": "yapcad-geometry-json-v0.3",
  "generator": {
    "name": "yapCAD",
    "version": "0.6.1",
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
| `metadata`     | object   | Metadata dictionary conforming to `metadata_namespace.rst`. |
| `boundingBox`  | number[] | `[xmin, ymin, zmin, xmax, ymax, zmax]` in document units. |
| `properties`   | object   | Derived properties (volume, area). |
| `representations` | object | Solid authority and available BREP/mesh forms. |

Each entity's `metadata` block MUST include the root fields from `metadata_namespace.rst`, including `layer` (defaulting to `"default"`). Serialisers SHOULD propagate layer assignments for solids down to child surfaces and sketches so that viewers can offer layer-based visibility controls.

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

`voids` lists interior cavity shells, one array of surface ids per cavity. It
is a **derived, interchange-level** field: an in-memory yapCAD mesh solid has
no voids slot, because a cavity is simply a shell whose faces wind inward, and
`volumeof` already subtracts it correctly. Serialisers compute the split with
`yapcad.geom3d.solid_shells`, which groups surfaces into connected shells by
shared edges and classifies each by the sign of its volume; the field is `[]`
when no cavity is found. Importers that build a mesh solid MUST merge void
surfaces back into the shell surface list rather than into a separate slot, and
may recover the split on demand with the same function.

Version 0.2 separates geometry representation from descriptive metadata. When
a solid originates from an analytic OCC BREP, that BREP is authoritative and
the indexed triangle shell is only a portable preview:

```json5
"representations": {
  "authoritative": "brep",
  "mesh": {
    "format": "indexed-triangle-set",
    "role": "preview"
  },
  "brep": {
    "role": "authoritative",
    "format": "opencascade-brep",
    "encoding": "base64",
    "payload": "R0hJTy4uLg==",
    "hash": "sha256:...",
    "kernel": {
      "name": "OpenCASCADE",
      "version": "7.8.1"
    },
    "modelTolerance": 1e-7
  }
}
```

The hash covers the decoded BREP bytes. Importers MUST verify it before kernel
loading. OCC-enabled yapCAD additionally requires successful rehydration and a
valid topology. Consumers without OCC can still verify the payload and use the
preview mesh. A mesh-only solid instead declares ``authoritative: "mesh"`` and
uses ``role: "authoritative"`` on its mesh record.

The BREP block is not duplicated under ``metadata``. The optional
``modelTolerance`` is expressed in document units and, when present, must be
finite and positive.

#### 2.1.1 SDF representation (v0.3)

Version 0.3 adds ``sdf`` to the ``authoritative`` enum. Extending that enum
is a breaking change for a strict reader, which is why the schema ID moved.
Authority remains a single per-solid property: a solid authored as a field
is SDF-authoritative and its mesh is a regenerated preview, while a solid
imported from STEP stays BREP-authoritative. A document declaring both is
rejected.

```json5
"representations": {
  "authoritative": "sdf",
  "sdf": {
    "role": "authoritative",
    "format": "yapcad-sdf-tree-v1",
    "lipschitz": "bound",
    "lipschitzConstant": 1.0,
    "csgExact": true,
    "bounds": [-10.0, -10.0, -4.0, 10.0, 10.0, 4.0],
    "tree": {
      "format": "yapcad-sdf-tree-v1",
      "root": "n5bb0477e7c48f027",
      "nodes": { "n…": { "kind": "sphere", "params": {"radius": 5.0},
                         "children": [] } }
    }
  },
  "mesh": {
    "format": "indexed-triangle-set",
    "role": "preview",
    "generatedFrom": "sdf",
    "method": "dual-contouring",
    "resolution": 64
  }
}
```

The node table is flat and keyed by content digest rather than nested, so a
shared subtree is stored once and the DAG structure survives the round trip.
Digests are truncated SHA-256 taken over each node's kind, parameters and
its children's digests, and are stable across processes so that a stored
tree does not churn package hashes.

``lipschitz`` distinguishes a field that *is* the signed distance
(``"exact"``) from one that only bounds it (``"bound"`` — hard CSG and
smooth blends). ``lipschitzConstant`` is the ``L`` for which
``|f(p)| <= L * d(p)``, so a ray marcher may safely step ``|f| / L``.
``csgExact`` reports whether the tree holds only analytic primitives, hard
booleans and similarity transforms and can therefore be replayed through
OCC for an exact BREP; ``false`` means the part will not round-trip to STEP
as analytic geometry. ``bounds`` is omitted for an unbounded field.

These three fields are a cache of what the tree already determines.
Importers MUST recompute and verify them rather than trust them, in the same
way the BREP payload hash is verified before kernel loading.

An SDF-authoritative solid MUST NOT also carry a top-level ``construction``
field: the tree is already in ``representations.sdf.tree``, and a second
copy is a second encoding of one fact that can drift. The mesh record's
meshing parameters are recorded so that the preview is exactly reproducible,
which matters because yapCAD signs packages.

#### 2.1.2 Construction provenance (optional)

A solid MAY carry a `construction` record describing how it was produced. It
mirrors the `construction` slot of the in-memory solid
(`['solid', surfaces, material, construction]`) and is a list whose first
element is a non-empty string *kind* tag, followed by kind-specific payload:

```json5
"construction": ["procedure", "yapcad.geom3d_util.sphere(4.0,center=[0,0,0,1],depth=2)"]
"construction": ["boolean", "union"]
"construction": ["boolean", "trimesh:difference"]
```

Kind tags currently emitted are `procedure` (a generating call) and `boolean`
(a CSG operation, optionally prefixed with the engine name). `sdf` is reserved
for the signed-distance-function tree described in `SDF-DESIGN.md`.

The field is **optional and additive**, and is omitted entirely for solids
that record no provenance — so the schema ID does not bump. Producers MUST
emit only JSON-representable payloads; `yapcad.construction` coerces tuples to
lists and falls back to `repr()` for anything else, so serialisation never
fails on provenance alone. Consumers that encounter a present-but-malformed
record (one that is not a list, or whose first element is not a non-empty
string) MUST reject the document rather than silently discard the record.

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
5. Document-level migrators SHOULD preserve unknown fields for forward
   compatibility. Geometry loaders MAY ignore fields they do not understand.
6. BREP payload hashes MUST be checked even when OpenCASCADE is unavailable.
7. Exporters MUST derive STEP and manufacturing STL from authoritative BREP
   when one is present; the preview mesh is not a manufacturing source.

---

## 6. Integration Guidance

- **Geometry export**: implement `to_geometry_json(entity)` that walks solids, surfaces, metadata and writes this schema.
- **Import**: validate `schema` version, then rebuild yapCAD list structures; reattach metadata via helper functions.
- **Manifest use**: `.ycpkg` manifests reference geometry JSON by path and hash.
  The JSON then authenticates its decoded BREP payload independently.
- **Streaming**: allow chunked outputs by splitting `entities` across files and referencing them via `attachments` or manifest entries.

---

## 7. Compatibility

Readers accept ``yapcad-geometry-json-v0.1`` and ``v0.2`` documents. The v0.1
historical ``metadata.brep`` payload is rehydrated when available. Writers
always emit v0.3 and never place new BREP data in generic metadata.

v0.3 differs from v0.2 only by admitting ``"sdf"`` as an ``authoritative``
value, with the accompanying ``sdf`` representation record and the optional
meshing parameters on ``mesh``. A v0.2 document therefore reads unchanged,
and a v0.2 document claiming ``authoritative: "sdf"`` is rejected — which is
the whole reason the ID had to move rather than the enum quietly widening.

The machine-readable schemas are
``docs/schemas/yapcad-geometry-json-v0.3.schema.json`` and, for the previous
revision, ``docs/schemas/yapcad-geometry-json-v0.2.schema.json``.

---

## 8. Future Enhancements

- Compression guidelines (e.g. `.json.zst`).
- Support for parametric feature history (link to DSL once available).
- Extend relationships to cover constraint solving results and tolerance stacks.
