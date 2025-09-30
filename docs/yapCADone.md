# yapCAD 1.0 Requirements & Roadmap

## Vision
Deliver a requirements-driven, provenance-aware design platform where yapCAD projects encapsulate specifications, parametric sources, validation workflows, and exportable geometry within a unified, verifiable package ready for collaboration with CAD/FEA/simulation ecosystems.

## Guiding Principles
- Traceability: every geometry artefact links back to requirements, parametric sources, tests, and toolchain versions.
- Modularity: separate geometry core, viewers, exporters, and validation tooling.
- Openness: use documented schemas and standard formats (STEP, STL, JSON/YAML manifests) for interoperability.
- Automation friendly: support LLM-driven design loops and continuous validation pipelines.

## Functional Requirements
1. **Project Packaging**
   - Define a `.ycpkg` package (directory or archive) containing manifest, requirements, design sources, validation definitions, results, exports, metadata.
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

## Non-Functional Requirements
- **Backward compatibility:** provide migration tools from 0.x projects; maintain 0.x runtime for legacy designs.
- **Performance:** metadata tracking should not significantly degrade geometry operations; exporters must scale to moderately complex assemblies.
- **Security:** hashing/signature features optional but easy to enable; documentation on best practices.
- **Extensibility:** manifest schema versioned; allow third-party extensions (custom validators, exporters).

## Roadmap & Milestones

### Phase 1 – Shared Foundations (v0.x)
- Implement mesh view utilities, metadata extensions, validation helpers (see `yapCADfoundations.md`).
- Ship STL exporter/importer using new helpers.

### Phase 2 – Project Packaging Prototype
- Introduce manifest schema and `.ycpkg` structure (manifest, design, tests, results, exports).
- Implement loader/saver preserving existing geometry plus new metadata; include migration scripts.
- Add basic CLI for packaging and validation status reporting.

### Phase 3 – Parametric DSL & Validation Layer
- Implement the declarative DSL compiler to Python/macros, including static validation and metadata emission.
- Build validation framework (test definitions, execution records, status propagation).
- Integrate with continuous testing workflows; implement CLI/API for LLM agents.
- Document spec/validation authoring process.

### Phase 4 – Export/Import Expansion
- STEP tessellated exporter leveraging mesh helper; extend to BREP where analytic surfaces available.
- STEP/STL import pipelines mapping to project structures, performing integrity checks.
- Viewer refactor to consume packaged projects.

### Phase 5 – Provenance & Security Enhancements
- Incorporate optional hashing/signature support for specs, geometry, and manifests.
- Provide verification tools and documentation for secure pipelines.

### Phase 6 – Release yapCAD 1.0
- Finalize documentation, migration guides, and API references.
- Update tutorials/examples to use new project model and workflows.
- Tag 1.0 release, deprecate legacy APIs with clear transition plan.

## Dependencies & Tooling Considerations
- Potential third-party STEP libraries (e.g. pythonocc-core) for parsing/advanced export—evaluate licensing and integration cost.
- Hashing/signature libraries (cryptography, nacl) for provenance features.
- Testing infrastructure to run simulations (container support, job orchestration where required).

## Risks & Mitigations
- **Complex migration:** provide automated conversion scripts and dual-loading support during transition.
- **Metadata proliferation:** keep schemas lean, allow opt-in extensions rather than mandatory fields.
- **External tool dependencies:** isolate behind adapters, provide graceful degradation when unavailable.

## Open Questions
- Level of BREP fidelity required for initial STEP release.
- Signature trust model (self-signed vs. PKI integration).
- Integration story for non-visual viewers and headless pipelines.
