# yapCAD Package (`.ycpkg`) Specification

**Version:** `ycpkg-spec-v0.1`  
**Status:** Draft – initial workflow design  

---

## 1. Goals

The `.ycpkg` package provides a portable container for yapCAD projects, bundling geometry, metadata, validation artefacts, and exports. It is designed to:

- Preserve provenance and traceability (aligns with `metadata-namespace-v0.1`).
- Enable reproducible automation (LLM agents, CI pipelines).
- Support incremental evolution (versioned manifests, optional extensions).

Distribution format is a directory tree or zipped archive. Commands must work on either form.

---

## 2. Directory Layout

```
my_design.ycpkg/
├── manifest.yaml
├── geometry/
│   ├── primary.json          # native yapCAD geometry (schema v0.1)
│   └── derived/…             # optional derived geometry files
├── metadata/
│   ├── requirements.yaml     # requirement statements + trace links
│   └── history.log           # change log (optional)
├── validation/
│   ├── plans/…               # definitions per test/tool
│   └── results/…             # captured outputs (JSON/CSV/etc.)
├── exports/
│   ├── model.step
│   └── model.stl
├── attachments/
│   └── …                     # arbitrary supporting files
└── README.md                 # human guidance (optional)
```

All file references recorded in the manifest must use relative paths within the package.

---

## 3. Manifest (`manifest.yaml`)

High-level structure:

```yaml
schema: ycpkg-spec-v0.1
id: d7d0f8b5-…
name: Rocket Bulkhead
version: 0.1.0
description: >
  Parametric bulkhead design produced via yapCAD.
created:
  timestamp: 2024-10-12T18:04:00Z
  author: rdevaul
generator:
  tool: yapCAD
  version: 0.5.1
  build: sha256:…
units: mm
tags: [prototype, waterjet]

geometry:
  primary:
    path: geometry/primary.json
    hash: sha256:…
    schema: yapcad-geometry-json-v0.1
    entities: ["solid-main", "surface-top"]
    # Each entity's metadata.layer encodes the logical drawing layer (default "default").
  derived:
    - path: geometry/derived/offset.json
      purpose: "shell offset for analysis"
      hash: sha256:…

metadata:
  requirements: metadata/requirements.yaml
  history: metadata/history.log

validation:
  plans:
    - id: mass-check
      path: validation/plans/mass_check.yaml
  results:
    - plan: mass-check
      path: validation/results/mass_check.json
      status: passed
      timestamp: 2024-10-12T19:00:00Z

exports:
  - id: step-main
    kind: step
    path: exports/model.step
    hash: sha256:…
    sourceEntities: ["solid-main"]
  - id: stl-main
    kind: stl
    path: exports/model.stl
    hash: sha256:…

attachments:
  - id: prompt
    path: attachments/design_prompt.txt
    description: "LLM prompt that generated initial geometry"
    hash: sha256:…

provenance:
  parent: null
  revisions:
    - version: 0.0.1
      timestamp: 2024-09-29T12:00:00Z
      notes: Initial concept import
      author: assistant

extensions:
  # Optional vendor-specific data
  com.example.workflow:
    state: approved
```

Key notes:

- `hash` values use lowercase algorithm names (`sha256`, `blake3`, …).
- `geometry.primary.entities` lists IDs found in the JSON geometry file to speed lookup.
- `extensions` is the sanctioned namespace for custom data.

---

## 4. CLI Workflow

### 4.1 Package Creation

`yapcad package create <source.py> -o my_design.ycpkg/ [options]`

Steps:
1. Run user-supplied script to generate geometry.
2. Collect metadata from geometry objects.
3. Emit `geometry/primary.json` via `geometry_to_json`.
4. Optionally run exports (`--export step,stl`).
5. Assemble manifest and write to disk.

### 4.2 Validation Update

`yapcad package validate my_design.ycpkg/ --plan validation/plans/mass_check.yaml`

Steps:
1. Load manifest.
2. Execute specified validation plan (commands defined in YAML).
3. Store outputs under `validation/results/` with status.
4. Update manifest `validation.results`.

### 4.3 Metadata Edit / Sync

`yapcad package annotate my_design.ycpkg/ --material "6061-T6"`

Updates geometry metadata via helper APIs, rewrites `geometry/primary.json`, refreshes manifest hashes.

### 4.4 Export Refresh

`yapcad package export my_design.ycpkg/ --kind step`

Regenerates exports, updates manifest.

---

## 5. API Surface

Implement Python helpers under `yapcad.package`:

- `PackageManifest`: dataclass with `load(path)`, `save()`, `recompute_hashes()`.
- `create_package(source: Callable, target_dir: Path, options)` – drives creation workflow.
- `load_geometry(manifest: PackageManifest)` – returns list of yapCAD entities.
- `update_metadata(manifest, modifier_fn)` – convenience to mutate metadata and persist.
- `run_validation(manifest, plan_id, runner)` – hooks to validation subsystem.

These APIs should be usable by both CLI and programmatic automation.

---

## 6. Hashing Rules

- Use SHA-256 as default; allow override via CLI (`--hash blake3`).
- Hashes computed over raw bytes; store as `sha256:<hex>`.
- When updating a file, manifest must be rewritten with new hash.

---

## 7. Zipped Packages

- `.ycpkg` may be zipped (`zip -r my_design.ycpkg.zip my_design.ycpkg/`).
- CLI must accept `.zip` by transparently mounting to temp dir.
- Manifest should include `packageFormat: directory|zip`.

---

## 8. Future Enhancements

- Digital signatures (`signatures` section) referencing PKI chain.
- Dependency graph for multi-part assemblies.
- Delta packages storing overrides against a base `.ycpkg`.
