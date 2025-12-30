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

Implementation Status (December 2025)
-------------------------------------

This section reflects the **actual current state** of the codebase as of version 0.6.2.

**Implemented and Working:**

- ``.ycpkg`` packaging with manifest schema and CLI tooling (``ycpkg_validate``, ``ycpkg_export``)
- DXF export for 2D geometry (lines, arcs, ellipses, splines, polygons, regions)
- STEP export (tessellated and analytic modes via ``YAPCAD_STEP_FORMAT=analytic``)
- STL export and import (binary and ASCII formats)
- STEP import with BREP topology preservation
- BREP kernel with OCC integration for solid modeling (sweeps, lofts, booleans)
- DSL with lexer, parser, type checker, and runtime interpreter
- DSL CLI (``python -m yapcad.dsl``) with ``check``, ``run``, ``list`` commands
- DSL builtins: primitives, booleans, transforms, sweeps (including adaptive), 2D regions
- DSL conditional expressions (``value if condition else other``)
- DSL list comprehensions (``[expr for var in iterable]`` with optional filter)
- DSL functional combinators: ``union_all``, ``difference_all``, ``intersection_all``, ``sum``, ``product``, ``min_of``, ``max_of``, ``any_true``, ``all_true``, and 2D equivalents
- 2D boolean operations (union, difference, intersection) with proper hole accumulation
- Curve types: ellipse, catmull-rom splines, NURBS, parabola, hyperbola
- Adaptive sweep operations with tangent-tracking profile orientation
- 533+ regression tests across geometry, DSL, import/export, and packaging

**Partially Implemented:**

- DSL packaging integration (``--package`` option works but provenance is basic)
- FEA integration (demo in ``examples/thrust_structure/`` - not production-ready)
- Package validation (basic schema validation, no solver orchestration)

**Not Yet Implemented:**

- Validation test definition language and schema
- Solver adapter framework (generalized FEA/simulation integration)
- Cryptographic signatures and approval workflows (Phase 5)
- Migration tooling for 0.x → 1.0 conversion
- Full provenance chain with audit trails

**Dependencies:**

- **pythonocc-core (OCC)**: Required for BREP operations, STEP import/export, and solid booleans. Install via conda (see ``environment.yml``). Without OCC, yapCAD operates in 2D/tessellation-only mode.
- **ezdxf**: Required for DXF export
- **numpy**: Core dependency for geometry operations

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
   - Provide an escape hatch to full Python for advanced cases; mark such sources as "untrusted" unless signed/approved, and capture warnings in project metadata.
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


**Phase 1 - Shared Foundations (v0.x) [COMPLETE]**

- Mesh/metadata utilities, viewer refactor, and geometry JSON enhancements delivered.
- STL/STEP exporters available.

**Phase 2 - Project Packaging Prototype [COMPLETE]**

- Manifest schema, ``.ycpkg`` layout, CLIs (``ycpkg_validate``, ``ycpkg_export``), external asset support.
- Regression tests cover packaging round-trips.

**Phase 3 - Parametric DSL & Validation Layer [LARGELY COMPLETE]**

*Implemented:*

- Full DSL implementation: lexer, parser, AST, type checker, symbol table, runtime interpreter
- DSL CLI with ``check``, ``run``, ``list`` commands and ``--output``/``--package`` options
- Extensive builtins: primitives, booleans, transforms, sweeps (including adaptive), regions, paths
- 2D geometry: ellipse, catmull-rom, NURBS, polygon, disk, 2D booleans with hole accumulation
- DXF export for 2D geometry visualization
- Basic provenance tracking in DSL execution
- **Conditional expressions**: ``value if condition else other`` syntax (December 2025)
- **List comprehensions**: ``[expr for var in iterable if condition]`` with map and filter
- **Functional combinators** (December 2025):
  - 3D boolean aggregation: ``union_all``, ``difference_all``, ``intersection_all``
  - 2D boolean aggregation: ``union2d_all``, ``difference2d_all``, ``intersection2d_all``
  - Numeric aggregation: ``sum``, ``product``, ``min_of``, ``max_of``
  - Boolean aggregation: ``any_true``, ``all_true``

*Remaining:*

- Validation test definition language (schema specification needed)
- Solver adapter framework (generalize FEA demo)

**Phase 4 - Export/Import Expansion [COMPLETE]**

- STEP export: tessellated (default) and analytic (``YAPCAD_STEP_FORMAT=analytic``)
- STEP import with BREP topology preservation
- STL import/export (binary and ASCII)
- DXF export for 2D geometry
- Full BREP kernel with OCC integration

**Phase 5 - Provenance & Security Enhancements [NOT STARTED]**

- Basic hashing exists for geometry/assets
- Cryptographic signatures not implemented
- Approval workflows not implemented
- **Decision needed**: Scope for 1.0 vs defer to 1.1

**Phase 6 - Release yapCAD 1.0 [NOT STARTED]**

Blockers:

- Phase 3 validation language formalization
- Documentation alignment with implementation
- Decision on Phase 5 scope

Dependencies & Tooling
----------------------

**Required for full functionality:**

- **pythonocc-core** (via conda): BREP operations, STEP import/export, solid booleans. Without this, yapCAD operates in reduced-functionality mode (2D geometry, tessellated solids only).

**Included dependencies:**

- ezdxf: DXF export
- numpy: Core geometry operations
- pyyaml: Manifest parsing

**Optional:**

- pyglet: Interactive viewer
- FEniCSx/Gmsh: FEA validation (demo only)

Risks & Mitigations
-------------------

- **Complex migration:** provide automated conversion scripts and dual-loading support during transition.
- **Metadata proliferation:** keep schemas lean, allow opt-in extensions rather than mandatory fields.
- **External tool dependencies:** isolate behind adapters, provide graceful degradation when unavailable.
- **OCC dependency:** Document clearly that OCC is required for solid modeling; provide meaningful error messages when unavailable.

Open Questions
--------------

- Signature trust model (self-signed vs. PKI integration) - **defer to post-1.0**
- Validation test definition language specification - **needed for 1.0**
- Integration story for non-visual viewers and headless pipelines

Priorities for 1.0 Release
--------------------------

**Required for 1.0:**

1. **Documentation Accuracy** - Ensure all docs reflect actual implementation state (this document, ``dsl_spec.rst``, ``yapBREP.rst``).

2. **DSL Reference Documentation** - Complete builtins reference with examples for all functions.

3. **OCC Dependency Clarity** - Document that OCC is effectively required for solid modeling; update error messages.

**Recommended for 1.0:**

4. **Validation Test Schema** - Define minimal schema for validation tests storable in packages.

5. **DSL Enhancements** - ✅ LARGELY COMPLETE (December 2025):
   - ✅ Conditional expressions implemented
   - ✅ List comprehensions with filter support
   - ✅ Functional combinators (``union_all``, ``sum``, etc.)
   - Remaining: better error messages with source context

6. **Test Coverage** - Ensure new features (2D booleans, DXF export, adaptive sweep, functional combinators) have adequate test coverage.

**Defer to 1.1:**

7. **Phase 5 Security** - Cryptographic signatures and approval workflows.

8. **Migration Tooling** - 0.x → 1.0 package conversion (API not yet stable enough).

9. **Advanced Automation APIs** - Full LLM agent interfaces beyond current CLI.

10. **Additional Import Formats** - IGES, OBJ, 3MF support.
