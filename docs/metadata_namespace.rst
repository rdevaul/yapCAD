yapCAD Metadata Namespace Specification
=======================================


**Version:** ``metadata-namespace-v1.0``
**Status:** Stable – yapCAD 1.0  
**Scope:** Metadata dictionaries attached to solids, surfaces, geometry collections, assemblies, and exported artefacts.



1. Overview
-----------

yapCAD embeds metadata alongside geometry objects using dictionary slots (see ``src/yapcad/metadata.py``). This namespace specification formalises those dictionaries so that materials, manufacturing intent, design provenance, and constraints can be shared between authoring tools, automation agents, and downstream exporters.

Every metadata dictionary MUST include:

| key          | type     | description |
|--------------|----------|-------------|
| ``schema``     | string   | Version identifier, e.g. ``metadata-namespace-v0.1``. |
| ``entityId``   | string   | Stable UUID for the geometry entity (use ``ensure_*_id``). |
| ``timestamp``  | string   | Optional ISO 8601 timestamp of last update. |
| ``tags``       | string[] | Free-form labels (e.g. ``["prototype", "revision-A"]``). |
| ``layer``      | string   | Logical drawing layer for the entity (``"default"`` if unspecified). |

Additional sections are nested dictionaries keyed by namespace names described below. Implementations MAY omit an entire section; consumers must treat unknown sections as opaque extension points.



2. Material Namespace (`material`)
----------------------------------


| key            | type          | notes |
|----------------|---------------|-------|
| ``name``         | string        | Human-readable descriptor (``"6061 Aluminum"``, ``"PLA"``). |
| ``standard``     | string        | Reference to published spec (``"ASTM B209-14"``). |
| ``grade``        | string        | Alloy/temper or plastic grade (``"6061-T6"``). |
| ``density_kg_m3``| number        | Density if known. |
| ``source``       | string        | Supplier, stock number. |



3. Manufacturing Namespace (`manufacturing`)
--------------------------------------------


| key               | type     | notes |
|-------------------|----------|-------|
| ``process``         | string   | Primary process (``"waterjet"``, ``"SLA"``, ``"3-axis CNC"``). |
| ``instructions``    | string   | Operator notes or setup text. |
| ``fixtures``        | string[] | Fixture/tooling identifiers. |
| ``layers``          | object   | Optional map relating drawing layers or geometry subsets to process directives (e.g. ``{ "layer0": {"operation": "cut", "draft_deg": 0}, "layer1": {"operation": "drill", "tool": "Ø6mm"} }``). |
| ``postprocessing``  | string[] | e.g. ``["deburr", "anodise-clear"]``. |



4. Design History Namespace (`designHistory`)
---------------------------------------------


| key           | type      | notes |
|---------------|-----------|-------|
| ``author``      | string    | Person or agent responsible. |
| ``source``      | string    | ``"prompt"``, ``"script"``, ``"import"`` etc. |
| ``context``     | string    | Prompt text, script path, or design brief summary. |
| ``tools``       | string[]  | Tools/agents used (``["LLM-gpt4.1", "topology-opt-v2"]``). |
| ``iterations``  | object[]  | Chronological entries (``[{ "revision": "A", "timestamp": "...", "notes": "..."}]``). |



5. Constraint Namespace (`constraints`)
---------------------------------------


| key             | type      | notes |
|-----------------|-----------|-------|
| ``mass``          | object    | ``{ "max_kg": 5.0, "target_kg": 4.2 }``. |
| ``envelope``      | object    | Bounding box or volume limits (``{ "x_mm": 300, "y_mm": 150, "z_mm": 80 }``). |
| ``performance``   | object[]  | Domain-specific constraints (e.g. load, stiffness). |
| ``compliance``    | string[]  | Required certifications (``["ASME Y14.5-2018"]``). |



6. Analysis Namespace (``analysis``)
------------------------------------

Captures analysis intent, solver information, and result references for the entity.

| key             | type      | notes |
|-----------------|-----------|-------|
| ``plans``         | string[]  | IDs referencing ``manifest.validation.plans`` entries executed against this entity. |
| ``loadcases``     | object[]  | Each item describes a load case with fields such as ``id``, ``description``, ``loads`` (forces/pressures/thermal), and ``boundaryConditions``. |
| ``materials``     | object    | Overrides mapping entity regions to material records (e.g. ``{"default": "metadata/materials/al6061.yaml"}``). |
| ``results``       | object[]  | Normalised summaries referencing ``validation/results`` artefacts. Recommended keys: ``plan`` (ID), ``status``, ``metrics`` (dict of key → value), ``summary`` (text), ``path`` (relative to package root). |
| ``safetyFactors`` | object    | Calculated safety factors per load case, e.g. ``{"axial-load": 1.8}``. |
| ``backendHints``  | object    | Optional hints for solver adapters (meshing strategy, element size, preferred backend). |
| ``notes``         | string    | Free-form context for analysts. |



7. Custom Extensions
--------------------


- Extension namespaces MUST be top-level keys using reverse-domain or organisation prefixes, e.g. ``"com.example.additive": {...}``.
- Extension content MUST NOT overwrite standard namespaces.



8. Usage Guidelines
-------------------


1. Call ``ensure_*_id`` before populating metadata to guarantee ``entityId``.
2. When updating metadata, bump ``timestamp`` and keep prior entries in ``designHistory.iterations``.
3. Serialisers SHOULD preserve unknown keys to allow forward compatibility.
4. Exporters SHOULD include the metadata dictionary (or a digest reference) in their manifest entries.



9. Future Work
--------------


- Define JSON Schema files for programmatic validation.
- Add localisation support for human-readable strings (multiple languages).
- Align metadata namespace with planned ``.ycpkg`` manifest structure.
