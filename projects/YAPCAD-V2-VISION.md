# yapCAD 2.0 — Vision Document

**Status:** Draft v0.1
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

yapCAD 2.0 ships with a curated set of **signed agent skills** that provide domain-specific intelligence. Skills are the mechanism for extending the platform — not plugins, not extensions, not arbitrary code. Skills are cryptographically signed, sandboxed, and versioned.

## Bundled Skills

### `yapcad-install`
**Purpose:** Bootstrap and maintain a yapCAD environment.
- Conda environment setup (pythonocc, trimesh, dependencies)
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
- Mesh generation from solid geometry (tetrahedral, hex)
- Linear static analysis (stress, strain, displacement)
- Thermal analysis
- Modal analysis (natural frequencies)
- Integration with external solvers (CalculiX, FEniCS, Elmer)
- Results visualization and design-space sensitivity analysis

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

### Session Authentication
- JWT tokens, generated conversationally (no login page)
- 4-hour TTL, HMAC-SHA256
- Multi-user: each team member gets their own session
- DSL sandbox: `python {}` and `native {}` blocks rejected for untrusted users

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
- [ ] OCC as default engine with deprecation warnings for fallbacks
- [ ] `lathe()`, `dome()`, `rounded_cylinder()` primitives

### Phase 2: Workbench
- [ ] Monaco DSL editor integration (Phase 1 ✅, refinement needed)
- [ ] JWT session authentication
- [ ] WebSocket chat with geometry context
- [ ] Signed skill loading and sandbox enforcement
- [ ] Multi-file project management

### Phase 3: Skills
- [ ] Skill packaging format and signing infrastructure
- [ ] `yapcad-install` skill
- [ ] `yapcad-cad3d` skill (wraps DSL + geometry operations)
- [ ] `yapcad-printing` skill (wraps existing print pipeline)
- [ ] `yapcad-geometry2d` skill

### Phase 4: GEL Integration
- [ ] Parallel design exploration framework
- [ ] Requirements-as-tests infrastructure
- [ ] Design family management (versioning, comparison, selection)
- [ ] Agent-mediated design review workflow

### Phase 5: Advanced Skills
- [ ] `yapcad-mechatronics` skill
- [ ] `yapcad-fea` skill (solver integration)
- [ ] `yapcad-dfm` skill
- [ ] External tool integration (FreeCAD, KiCad, CalculiX)

## What yapCAD 2.0 Is NOT

- **Not a SaaS product.** It runs on your machine, your network. No cloud dependency.
- **Not a Fusion 360 / SolidWorks replacement.** It's a different paradigm — agent-mediated, code-first, exploration-oriented. Use commercial CAD when you need their ecosystem.
- **Not model-specific.** Works with any LLM backend via OpenClaw (Claude, GPT, Qwen, local models). The GEL thesis applies to any capable model.
- **Not just for experts.** The agent lowers the floor — natural language intent → validated geometry. The DSL raises the ceiling — full parametric control for those who want it.

## Open Questions

1. **Packaging:** pip install with conda bootstrap? conda-only? Docker image?
2. **FEA solver:** Bundle CalculiX? Use FEniCS? Both?
3. **Skill marketplace:** DML-only skills, or open ecosystem with signing?
4. **File format:** Should yapCAD define a project format (.yapcad) that bundles DSL + requirements + assembly + print config?
5. **Collaborative editing:** Multiple agents/humans editing the same design simultaneously?

---

*This document will evolve. Update it as decisions are made and scope is refined.*
