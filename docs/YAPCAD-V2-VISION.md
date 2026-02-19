# yapCAD 2.0 — Vision Document

**Status:** Draft v0.2
**Authors:** Rich DeVaul, Jarvis
**Date:** 2026-02-19

---

## One-Liner

yapCAD 2.0 is an **agentic-first engineering platform** — a parametric CAD kernel, DSL, and web workbench designed for AI-mediated design exploration, fabrication, and iteration.

## The Big Shift

yapCAD 1.x was a Python CAD library with a DSL bolted on. yapCAD 2.0 inverts this: the **agent is the primary user**, the DSL is the agent's language for expressing geometry, and the web workbench is where humans supervise, guide, and collaborate.

The default workflow is no longer "human writes DSL → runs script → inspects output." It's:

1. Human describes intent and constraints in natural language
2. Agent explores the design space in parallel (Generative Engineering Loop)
3. Agent proposes candidates that pass requirements-based validation
4. Human reviews, selects, guides — agent iterates
5. Fabrication is downstream of selection, not a separate workflow

This is the GEL philosophy made concrete: **design as parallel search, fabrication as selection, iteration as the norm.**

## Architecture

```
┌─────────────────────────────────────────────────┐
│                 Web Workbench                     │
│  ┌──────────┐  ┌──────────┐  ┌───────────────┐  │
│  │  Monaco   │  │   3D     │  │  Agent Chat   │  │
│  │  DSL      │  │  Viewer  │  │  (WebSocket)  │  │
│  │  Editor   │  │ (Three.js)│  │               │  │
│  └──────────┘  └──────────┘  └───────────────┘  │
└────────────────────┬────────────────────────────┘
                     │ REST + WebSocket
┌────────────────────┴────────────────────────────┐
│              yapCAD Service (FastAPI)             │
│  ┌──────────────────────────────────────────┐   │
│  │  DSL Engine  │  Geometry Kernel  │  Skills │   │
│  │  (parse/eval)│  (OCC BREP)       │  (signed)│   │
│  └──────────────────────────────────────────┘   │
│  ┌──────────────────────────────────────────┐   │
│  │  Assembly  │  Kinematics  │  Collision    │   │
│  └──────────────────────────────────────────┘   │
│  ┌──────────────────────────────────────────┐   │
│  │  Package (.ycpkg)  │  FEA (FEniCS/Calc.) │   │
│  └──────────────────────────────────────────┘   │
└────────────────────┬────────────────────────────┘
                     │
┌────────────────────┴────────────────────────────┐
│           Agent Runtime (OpenClaw)                │
│  Skills: install, geometry2d, cad3d, mechatronics│
│          fea, dfm, printing, ...                  │
│  GEL: parallel exploration, req-based testing     │
└──────────────────────────────────────────────────┘
```

## Core Principles

### 1. Agentic-First

Every API, every data format, every workflow is designed for agent consumption first, human readability second. The DSL is a machine-legible language for parametric geometry. The REST API is the agent's interface. The viewer is the human's window into what the agent is doing.

### 2. Generative Engineering Loop (GEL)

Design is not a linear process. yapCAD 2.0 supports:
- **Parallel exploration:** Multiple design variants evaluated concurrently
- **Requirements as tests:** Constraints expressed as executable validation functions
- **Design families:** Output is a population of valid designs, not a single artifact
- **Iterative refinement:** Agent adjusts parameters, topology, or approach based on validation results
- **Fabrication as selection:** Printing/manufacturing a design doesn't end the exploration

### 3. Watertight by Default

All solid geometry constructors produce manifold, watertight meshes. The OCC BREP kernel is the default backend. Tessellated fallbacks exist but emit warnings. Every solid operation validates output geometry.

### 4. Skills as the Extension Model

yapCAD 2.0 ships with a curated set of **signed agent skills** that provide domain-specific intelligence. Skills are cryptographically signed, sandboxed, and versioned. The skill ecosystem is **open** — anyone can publish signed skills to the marketplace, not just DML. This marketplace will eventually become its own project, open to any signed skill (not limited to yapCAD).

### 5. Packages as Living Design Documents

The `.ycpkg` package format is the unit of design. A package isn't just geometry — it bundles DSL source, requirements, FEA/testing specs, assembly constraints, and (in 2.0) mechatronics linkages, manufacturing hints, and **mini-skills** that guide users through domain-specific design processes.

**Mini-skills** are package-embedded agent instructions that turn a parametric design into an interactive design tool. Examples:
- A **bicycle frame** package with a mini-skill that interviews the user about rider dimensions, terrain, and use case to derive tube angles and lengths
- A **sounding rocket** package with a mini-skill that guides the user through propulsion, structural, and aero requirements to produce a complete multi-stage design
- A **pressure vessel** package with a mini-skill that walks through material selection, working pressure, and safety factor to size wall thickness and dome geometry

This makes yapCAD packages not just reusable parts, but reusable *design processes*.

## Bundled Skills

### `yapcad-install`
**Purpose:** Bootstrap and maintain a yapCAD environment.
- Conda environment setup (pythonocc, trimesh, dependencies)
- Docker container setup (bundled with OpenClaw)
- Version checking and updates
- Environment validation (OCC available? Blender available? Printers configured?)
- First-run tutorial and capability inventory

### `yapcad-geometry2d`
**Purpose:** 2D computational geometry — yapCAD's origin and still essential.
- Path construction, boolean operations on 2D regions
- Bezier/NURBS curve manipulation
- 2D constraint solving
- SVG/DXF import and export
- Profile generation for extrusion, revolution, lofting

### `yapcad-cad3d`
**Purpose:** 3D solid modeling — the core CAD workflow.
- Parametric solid primitives (box, cylinder, sphere, cone, tube, etc.)
- Boolean operations (union, difference, intersection) via OCC
- Revolution solids with guaranteed watertight geometry
- Loft, sweep, helical extrude
- STEP/STL/3MF import and export
- Assembly definition with datums and mate constraints

### `yapcad-mechatronics`
**Purpose:** Mechanism and robotics design.
- Kinematic chain definition and forward/inverse kinematics
- Joint types (revolute, prismatic, spherical, planar)
- Range-of-motion sweeps and collision detection
- COTS component integration (servo datasheets, bearing specs)
- Cable/linkage routing
- Actuator sizing and selection

### `yapcad-fea`
**Purpose:** Finite element analysis and simulation.
- Mesh generation from solid geometry (gmsh integration, tetrahedral/hex)
- Linear static analysis (stress, strain, displacement)
- Thermal analysis
- Modal analysis (natural frequencies)
- Dual solver support: **FEniCS** (existing integration) and **CalculiX** (existing integration)
- Face naming for boundary condition assignment
- Results visualization and design-space sensitivity analysis

*Note: yapCAD already has preliminary FEniCS and CalculiX integration in `src/yapcad/package/analysis/`. The 2.0 skill wraps this existing infrastructure with agent-friendly interfaces.*

### `yapcad-dfm`
**Purpose:** Design for manufacture — bridging design and fabrication.
- Wall thickness analysis
- Overhang detection for additive manufacturing
- Undercut detection for molding/casting
- Tolerance analysis and GD&T
- Material selection guidance
- Cost estimation
- Part segmentation for multi-piece fabrication (from Jeremy's manufacturing module)

### `yapcad-printing`
**Purpose:** 3D printing pipeline — from solid to physical part.
- Slicer integration (OrcaSlicer, PrusaSlicer)
- Printer fleet management (Bambu, Prusa, etc.)
- Print orientation optimization
- Support structure analysis
- Filament/material database
- Print preview rendering
- Job submission and monitoring
- Post-print dimensional validation

## Dependency Architecture

### Required
- **Python 3.11+**
- **conda** (miniconda/miniforge) — manages the native dependency chain
- **pythonocc-core** — OCC BREP kernel (via conda-forge)
- **trimesh** — mesh analysis, import/export, validation
- **FastAPI + uvicorn** — service layer

### Recommended
- **Blender** (headless) — mesh boolean fallback, rendering
- **OpenClaw** — agent runtime (the primary way to use yapCAD 2.0)

### Deprecated (2.0)
- Native Python triangle-mesh boolean engine — replaced by OCC/trimesh:blender
- Direct Python API as primary interface — still supported, but agent-mediated workflow is the default
- Tessellated-only geometry path — OCC BREP is default, tessellation is for visualization only

### Packaging & Distribution

Two supported installation paths:

1. **conda-only:** `conda install -c dml yapcad` — installs the full environment including pythonocc, trimesh, and all dependencies. Recommended for development and direct use.

2. **Docker container (bundled with OpenClaw):** A complete containerized environment with yapCAD service, web workbench, and OpenClaw agent runtime pre-configured. Recommended for deployment and team use. Pull and run — no environment setup required.

Both paths deliver the same capabilities. Docker is the "it just works" option; conda is the "I want to hack on it" option.

## The `.ycpkg` Package Format

The yapCAD package (`.ycpkg`) is the fundamental unit of shareable design. It already supports geometry, DSL source, and FEA/testing specifications. For 2.0, the format extends to include:

### Current Capabilities
- Multi-entity geometry (solids, surfaces, assemblies)
- DSL source with parametric commands
- FEA specifications (boundary conditions, loads, materials)
- Validation tests (requirements as executable checks)
- Package signing (existing `yapcad.package.signing` module)

### 2.0 Extensions
- **Mechatronics metadata:** Kinematic chain definitions, joint limits, actuator specs
- **Manufacturing hints:** Preferred orientation, support strategy, material constraints, fit tolerances
- **Mini-skills:** Embedded agent instructions (natural language + structured prompts) that guide domain-specific design processes
- **Design lineage:** Links to parent designs, variant relationships, exploration history
- **Version control hooks:** Integration with git for checkpoint/merge workflows

### Mini-Skills in Packages

A mini-skill is a structured prompt embedded in a `.ycpkg` that tells an agent how to use the package interactively. It includes:

```yaml
mini_skill:
  name: "Pressure Vessel Designer"
  description: "Interactive design tool for cylindrical pressure vessels"
  interview:
    - question: "What is the working pressure (MPa)?"
      parameter: working_pressure
      type: float
      range: [0.1, 200]
    - question: "What material family?"
      parameter: material
      type: choice
      options: [aluminum, titanium, steel, composite]
    - question: "What is the target internal volume (liters)?"
      parameter: volume
      type: float
      range: [0.1, 10000]
  validation:
    - check: wall_thickness > 2.0mm
    - check: safety_factor >= 2.5
    - check: manifold(outer_shell)
  workflow:
    - step: derive_dimensions
    - step: generate_geometry
    - step: run_fea
    - step: validate_requirements
    - step: propose_to_user
```

This turns any parametric package into an **interactive design wizard** — the agent asks the right questions, fills in the parameters, validates the result, and presents it to the human.

## The Web Workbench

The workbench is where human and agent collaborate. It's a browser application, not a desktop app. It runs locally (no cloud dependency) but is accessible from any device on the network.

### Layout
- **Left panel:** Monaco DSL editor + file management
- **Center:** Three.js 3D viewer (single view + 4-up ortho)
- **Right:** Agent chat panel (WebSocket, geometry-aware)
- **Bottom:** Validation results, print queue, activity log

### Agent Integration
The chat panel isn't a generic chatbot — it's geometry-aware:
- Agent sees the current DSL source, viewport state, and loaded assemblies
- Agent can modify DSL directly (push edits to editor)
- Agent can trigger evaluation, rendering, validation
- Agent proposes design alternatives with visual diffs
- Human can accept, reject, or redirect with natural language
- **Mini-skills** activate when a package with embedded skills is loaded

### Session Authentication
- JWT tokens, generated conversationally (no login page)
- 4-hour TTL, HMAC-SHA256
- Multi-user: each team member gets their own session
- DSL sandbox: `python {}` and `native {}` blocks rejected for untrusted users

### Collaborative Editing

Multiple agents and humans can edit the same design simultaneously. This requires:
- **Version control integration:** Design changes are periodically checkpointed to git
- **Conflict resolution:** Conflicting edits can be merged (if compatible) or kept as parallel explorations in the design space — consistent with the GEL philosophy that divergent designs are features, not bugs
- **Operational transform or CRDT:** For real-time collaborative editing of DSL source (stretch goal — initial implementation may use lock-based turn-taking)
- **Design branching:** Explicit support for forking a design into variants, with the option to merge back or maintain as a design family

The GEL framing makes collaborative editing more natural than in traditional CAD: when two people take a design in different directions, the result is a richer design family, not a merge conflict.

## DSL Evolution for 2.0

### New Primitives
- `lathe(profile_points, steps)` — revolve a 2D polyline, auto-capping
- `dome(radius, height)` — hemispherical/ellipsoidal cap
- `rounded_cylinder(radius, height, fillet_radius)` — cylinder with edge fillets
- `chamfer(solid, edge_selector, distance)` — edge chamfers
- `fillet(solid, edge_selector, radius)` — edge fillets
- `shell(solid, thickness, [faces_to_remove])` — shell a solid (hollow it out)

### Assembly DSL
- `datum(name, type, position, direction)` — define assembly datums
- `mate(part_a.datum, part_b.datum, type)` — constrain parts
- `joint(name, type, axis, limits)` — define kinematic joints
- `pose(assembly, joint_values)` — set joint positions

### Validation DSL
- `require(condition, message)` — assert a design requirement
- `check_manifold(solid)` — verify watertight geometry
- `check_clearance(solid_a, solid_b, min_distance)` — interference check
- `check_printable(solid, printer_profile)` — DFM validation
- `check_stress(solid, load, material, max_stress)` — FEA validation

## Release Plan

### Phase 1: Foundation (current)
- [x] 20 DSL builtins (Phases 1-4)
- [x] Revolution solid disc caps (watertight guarantee)
- [x] Assembly system (datums, mates, kinematics, collision)
- [x] Web viewer with 3-column layout
- [x] REST API service
- [x] FEA integration (FEniCS + CalculiX in `yapcad.package.analysis`)
- [x] Package signing (`yapcad.package.signing`)
- [ ] OCC as default engine with deprecation warnings for fallbacks
- [ ] `lathe()`, `dome()`, `rounded_cylinder()` primitives

### Phase 2: Workbench
- [ ] Monaco DSL editor integration (Phase 1 ✅, refinement needed)
- [ ] JWT session authentication
- [ ] WebSocket chat with geometry context
- [ ] Signed skill loading and sandbox enforcement
- [ ] Multi-file project management

### Phase 3: Skills & Packaging
- [ ] Skill packaging format and signing infrastructure
- [ ] Open skill marketplace (eventually its own project)
- [ ] `yapcad-install` skill
- [ ] `yapcad-cad3d` skill (wraps DSL + geometry operations)
- [ ] `yapcad-printing` skill (wraps existing print pipeline)
- [ ] `yapcad-geometry2d` skill
- [ ] Docker container distribution (yapCAD + OpenClaw bundled)
- [ ] conda package (`conda install -c dml yapcad`)

### Phase 4: GEL Integration & Collaboration
- [ ] Parallel design exploration framework
- [ ] Requirements-as-tests infrastructure
- [ ] Design family management (versioning, comparison, selection)
- [ ] Agent-mediated design review workflow
- [ ] Git-integrated version control for designs
- [ ] Collaborative editing (branching, merging, parallel variants)

### Phase 5: Advanced Skills & Package Extensions
- [ ] `yapcad-mechatronics` skill
- [ ] `yapcad-fea` skill (agent-friendly wrapper around existing solvers)
- [ ] `yapcad-dfm` skill
- [ ] Mini-skill embedding in `.ycpkg` packages
- [ ] `.ycpkg` extensions: mechatronics linkages, manufacturing hints, design lineage
- [ ] External tool integration (FreeCAD, KiCad)
- [ ] Package-embedded design wizards (bicycle frame, rocket, pressure vessel examples)

## What yapCAD 2.0 Is NOT

- **Not a SaaS product.** It runs on your machine, your network. No cloud dependency.
- **Not a Fusion 360 / SolidWorks replacement.** It's a different paradigm — agent-mediated, code-first, exploration-oriented. Use commercial CAD when you need their ecosystem.
- **Not model-specific.** Works with any LLM backend via OpenClaw (Claude, GPT, Qwen, local models). The GEL thesis applies to any capable model.
- **Not just for experts.** The agent lowers the floor — natural language intent → validated geometry. The DSL raises the ceiling — full parametric control for those who want it.

## Decisions Made

1. **Packaging:** Dual-path — conda-only for developers, Docker+OpenClaw for deployment. Both first-class.
2. **FEA solvers:** Support both FEniCS and CalculiX. Both already have preliminary integration in `yapcad.package.analysis`.
3. **Skill marketplace:** Open ecosystem with cryptographic signing. Will become its own project — the "Dark Matter Skill Universe" (see Garrett's `DARK-MATTER-SKILL-UNIVERSE.md`). Marketplace serves both agent skills AND signed `.ycpkg` packages. Reference implementation: `github.com/rdevaul/skill-signer`. **Critical requirement:** yapCAD package signing (`yapcad.package.signing`) must align with the skill-signer infrastructure — one identity, one key, signs both.

   **Multi-repo architecture:** The client maintains a configurable list of signed skill/package repositories. Each repo is an independent source — Gitea (DML internal), GitHub, GitLab, or hypothetically IPFS. The client doesn't care about the backend; it just needs a common protocol for discovery, download, and signature verification. This means:
   - DML team uses internal Gitea repos during development
   - Public releases go to GitHub/GitLab repos
   - Multiple repos can be active simultaneously (like apt sources or conda channels)
   - Git Actions on each repo handle automated security checks on upload
   - Test the full workflow on Gitea first, then expand to public repos
4. **Package format:** Extend existing `.ycpkg` with mechatronics, manufacturing hints, and mini-skills. The package is a living design document, not just a geometry container.
5. **Collaborative editing:** Yes, with git-integrated version control. Conflicts become parallel explorations (consistent with GEL). Real-time CRDT is a stretch goal.

---

*This is a living document. Update it as decisions are made and scope is refined.*
*v0.1 → v0.2: Integrated Rich's answers on packaging (conda + Docker), FEA (both solvers), marketplace (open + signed), .ycpkg extensions (mini-skills), and collaborative editing (git-integrated, GEL-aligned).*
