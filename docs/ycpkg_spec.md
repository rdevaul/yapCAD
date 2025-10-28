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
│   ├── entities/             # canonical solids/sketches produced by DSL commands
│   └── derived/…             # optional derived geometry files
├── instances/                # placements/replications referencing canonical entities
│   └── *.json
├── metadata/
│   ├── requirements.yaml     # requirement statements + trace links
│   └── history.log           # change log (optional)
├── src/                      # packaged DSL modules
│   └── *.dsl
├── scripts/                  # optional Python helpers invoked by DSL fallback blocks
│   └── *.py
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
    # Sketch entries now capture both sampled polylines and ``primitives``
    # (line/arc/circle/catmullrom/nurbs/polyline records) so round-tripping
    # preserves parametric fidelity for downstream exporters.
  derived:
    - path: geometry/derived/offset.json
      purpose: "shell offset for analysis"
      hash: sha256:…
    - path: geometry/derived/external_gear.step
      purpose: "imported downstream component"
      format: step
      source:
        kind: import
        original: /workspace/library/gear.step
      hash: sha256:…

instances:
  - id: gear_a
    entity: geometry/entities/gear_primary.json
    count: 4
    transforms:
      - instances/gear_a_1.json
      - instances/gear_a_2.json
  - id: bolt_m10
    entity: geometry/entities/bolt_m10.json
    count: 100

source:
  modules:
    - id: involute_gear
      path: src/involute_gear.dsl
      language: yapdsl-v0.1
      hash: sha256:…
      exports: [INVOLUTE_SPUR, INVOLUTE_SPUR2D]
    - id: metric_fasteners
      path: src/metric_fasteners.dsl
      language: yapdsl-v0.1
  runtime:
    python:
      helpers:
        - path: scripts/custom_helpers.py
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

Sketch entities now persist a `primitives` array alongside their sampled `polylines`. Each primitive captures the original parametric definition—`line`, `circle`, `arc`, `catmullrom`, `nurbs`, or explicit `polyline`—so geometry exchanged between yapCAD agents retains full fidelity. Export tools (e.g., `tools/ycpkg_export.py`) can therefore emit native DXF entities instead of approximating curves with short segments.

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
- `add_geometry_file(manifest, source_path, ...)` – copy an external STEP/STL/etc. into `geometry/derived/` (or attachments) and register it in the manifest with hash + provenance.

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
- Layer-aware viewer interactions (already implemented) may evolve into annotated layer libraries.
- DSL compilation pipeline will eventually support signed modules and hashed invocation metadata per entity.
- Canonical-entity instancing will feed BOM generation utilities.
- Digital signatures (`signatures` section) referencing PKI chain.
- Dependency graph for multi-part assemblies.
- Delta packages storing overrides against a base `.ycpkg`.
## Viewer & Validation Quick Reference

- `tools/ycpkg_validate.py <package>` – verifies hashes and geometry JSON.
- `tools/ycpkg_viewer.py <package>` – launches the interactive viewer.
  * 3D mode shows perspective/front/top/side quadrants, supports layer toggles (`1-9` to toggle, `0` reset) and a help overlay (`H`/`F1`).
  * 2D sketches receive the same layer toggles, pan/zoom controls, and help overlay.
│   ├── register/…            # signatures, approvals, etc. (optional)
### 4.5 DSL Compilation & Canonical Entities

- `yapcad dsl compile src/involute_gear.dsl` parses the DSL module, emits canonical solids/sketches under `geometry/entities/`, and records exported commands in the manifest `source.modules` block.
- `yapcad package assemble --dsl module:command=params` resolves DSL invocations, deduplicates canonical geometry, and populates the `instances/` directory plus `manifest.instances` entries.
- Python fallback blocks inside DSL modules must live under `scripts/` and are listed under `source.runtime` for reproducibility.
- `register_instance(entity_id, params, transform)` helpers (TBD) will let higher-level assemblies reuse canonical geometry and populate `manifest.instances`.
