yapCAD 1.0 Requirements & Roadmap
=================================


Vision
------

Deliver a requirements-driven, provenance-aware design platform where yapCAD projects encapsulate specifications, parametric sources, validation workflows, and exportable geometry within a unified, verifiable package ready for collaboration with CAD/FEA/simulation ecosystems.

Guiding Principles
------------------

- Traceability: every geometry artefact links back to requirements, parametric sources, tests, and toolchain versions.
- Modularity: separate geometry core, viewers, exporters, and validation tooling.
- Openness: use documented schemas and standard formats (STEP, STL, JSON/YAML manifests) for interoperability.
- Automation friendly: support LLM-driven design loops and continuous validation pipelines.

Progress Snapshot (November 2025)
---------------------------------

- ``.ycpkg`` packaging, manifest schema, and CLI tooling implemented (validation & export helpers, metadata tracking, analytic sketch primitives).
- DXF/STEP/STL import and export available; viewers operate on packaged geometry.
- DSL, validation framework, and security/signature features remain in design phase (``docs/dsl_spec.rst``, ``docs/ycpkg_spec.rst``).
- **BREP Integration COMPLETE**: Full native BREP representation with OCC integration (see ``docs/BREP_integration_strategy.md``).
- **Analytic STEP Export COMPLETE**: ``write_step_analytic()`` preserves exact geometric definitions (PLANE, CYLINDER, SPHERE, CONE, etc.).
- **STEP Import COMPLETE**: ``import_step()`` loads STEP files with full BREP topology preservation.
- **STL Import COMPLETE**: ``read_stl()`` / ``import_stl()`` load binary and ASCII STL files with optional vertex deduplication.
- **Advanced curve types**: Parabola and hyperbola primitives implemented.
- **291 regression tests** covering all import/export and BREP workflows.
- Upcoming work: DSL compiler, validation execution layer, automation APIs.

Functional Requirements
-----------------------

1. **Project Packaging**
   - Define a ``.ycpkg`` package (directory or archive) containing manifest, requirements, design sources, validation definitions, results, exports, metadata.
   - Manifest must include identifiers, versioning, dependency hashes, optional signatures.
2. **Specification Management**
   - Structured requirement documents with schemas describing functional, geometric, material, and testing constraints.
   - Revision history and linkage to design iterations.
3. **Parametric Design Representation**
   - Introduce a declarative DSL that compiles into trusted yapCAD macros plus metadata. The DSL covers common parametric features (profiles, sweeps, threads, hole patterns) and can be statically validated.
   - Provide an escape hatch to full Python for advanced cases; mark such sources as “untrusted” unless signed/approved, and capture warnings in project metadata.
   - Track dependencies between parametric parameters and resulting geometry artefacts.
4. **Validation Framework**
   - Test definition language capturing procedures (simulation scripts, unit tests, analysis configs).
   - Run manager recording tool versions, inputs, outputs, and pass/fail criteria; results stored in package.
5. **Export Pipelines**
   - STEP (tessellated and, where possible, analytic BREP) and STL exporters integrated into project workflow.
   - CLI/API to generate exports using project metadata, record references in manifest.
6. **Import Pipelines**
   - Ability to ingest STEP/STL (and future formats) mapping geometry and metadata into project structures.
   - Validate imported data against schemas and mark areas requiring manual confirmation.
7. **Metadata & Provenance**
   - Unique IDs for all geometry entities, change tracking, optional cryptographic signatures.
   - Capture toolchain details (yapCAD version, key libraries, simulation tool versions).
   - Support for author attribution, approvals, and audit trails.
8. **Viewer Decoupling**
   - Repackage the Pyglet viewer (and other visualization utilities) to operate on packaged projects without residing in core geometry modules.
9. **Automation Interfaces**
   - APIs for LLM agents to read/write specs, propose design variants, trigger validation, and compare outcomes.
   - Hooks for CI/CD style workflows (batch validation, export refresh).

Non-Functional Requirements
---------------------------

- **Backward compatibility:** provide migration tools from 0.x projects; maintain 0.x runtime for legacy designs.
- **Performance:** metadata tracking should not significantly degrade geometry operations; exporters must scale to moderately complex assemblies.
- **Security:** hashing/signature features optional but easy to enable; documentation on best practices.
- **Extensibility:** manifest schema versioned; allow third-party extensions (custom validators, exporters).

Roadmap & Milestones
--------------------


**Phase 1 - Shared Foundations (v0.x) [Completed]**
- Mesh/metadata utilities (``yapCADfoundations.md``), viewer refactor, and geometry JSON enhancements delivered.
- STL/STEP exporters available; STL import deferred to Phase 4 analytic work.

**Phase 2 - Project Packaging Prototype [Completed]**
- Manifest schema, ``.ycpkg`` layout, CLIs (``ycpkg_validate``, ``ycpkg_export``), external asset support (``add_geometry_file``) in place.
- Regression tests cover packaging round-trips; migration tooling still TBD for legacy 0.x projects.

**Phase 3 - Parametric DSL & Validation Layer [In Progress]**
- Specifications drafted (``docs/dsl_spec.rst``), but compiler/runtime and validation execution manager not yet implemented.
- Next steps: prototype DSL compiler, define validation schema, integrate with packaging/metadata.
- Validation plan schema (``docs/ycpkg_spec.rst``) and analysis metadata updates published; placeholder analyzer CLI records plan execution summaries pending full solver adapters.

**Phase 4 - Export/Import Expansion [COMPLETE]**
- STEP (faceted and analytic), STL, DXF exports implemented; viewer consumes packaged geometry.
- **STEP import implemented** via ``import_step()`` with full BREP topology preservation.
- **Analytic STEP export implemented** via ``write_step_analytic()`` - preserves exact surface definitions.
- Full BREP kernel with OCC integration complete (see ``docs/BREP_integration_strategy.md``).
- **STL import implemented** via ``read_stl()`` / ``import_stl()`` supporting binary and ASCII formats.

**Phase 5 - Provenance & Security Enhancements [Not Started]**
- Hashing exists for geometry/assets; signatures/approvals still on backlog.

**Phase 6 - Release yapCAD 1.0 [Not Started]**
- Depends on completing Phases 3–5 plus documentation and migration tooling.

Dependencies & Tooling Considerations
-------------------------------------

- **pythonocc-core integrated** via conda environment (``environment.yml``) - provides STEP import/export and OCC boolean engine.
- Hashing/signature libraries (cryptography, nacl) for provenance features.
- Testing infrastructure to run simulations (container support, job orchestration where required).
- **BREP implementation complete**: Native BREP topology with bidirectional OCC conversion (see ``docs/BREP_integration_strategy.md``).

Risks & Mitigations
-------------------

- **Complex migration:** provide automated conversion scripts and dual-loading support during transition.
- **Metadata proliferation:** keep schemas lean, allow opt-in extensions rather than mandatory fields.
- **External tool dependencies:** isolate behind adapters, provide graceful degradation when unavailable.

Open Questions
--------------

- ~~Level of BREP fidelity required for initial STEP release.~~ **RESOLVED**: Full analytic BREP with native topology graph and OCC integration implemented.
- Signature trust model (self-signed vs. PKI integration).
- Integration story for non-visual viewers and headless pipelines.
