# yapCAD Server API Specification

**Status:** Draft
**Purpose:** Agent-first API for yapCAD geometry generation, visualization, and analysis
**Primary users:** LLMs, design agents, GEL automation

---

## Overview

The yapCAD Server enables agentic design workflows where LLMs invoke DSL commands to generate geometry, assemble designs, and receive visual feedback. Every operation produces standardized views so agents can "see" results without parsing raw geometry.

```
┌─────────────┐     ┌─────────────────┐     ┌──────────────┐
│  LLM Agent  │────▶│  yapCAD Server  │────▶│  Web Viewer  │
│  (Claude)   │◀────│  (Python/FastAPI)│◀────│  (React/3JS) │
└─────────────┘     └─────────────────┘     └──────────────┘
        │                   │
        │                   ▼
        │           ┌──────────────┐
        └──────────▶│  FEA / Analysis │
                    └──────────────┘
```

---

## Core Concepts

### Workspace
A workspace is an active design session containing:
- **Entities**: Geometry primitives and solids
- **Assemblies**: Named collections of entities with transforms
- **History**: Undo/redo stack
- **Metadata**: Author, timestamps, constraints

### Package
A yapCAD package (`.ycpkg`) is a serialized workspace that can be saved, loaded, and versioned.

### View Set
Every geometry operation returns a **View Set** — standardized images for agent consumption:
- **Plan** (top-down, XY plane)
- **Front elevation** (XZ plane)
- **Side elevation** (YZ plane)  
- **Perspective** (isometric or user-defined camera)
- **Section** (optional, with specified cutting plane)

---

## API Endpoints

### Workspace Management

#### `POST /workspace/create`
Create a new workspace.

```json
Request:
{
  "name": "thruster_mount_v1",
  "units": "mm",
  "author": "jarvis"
}

Response:
{
  "workspace_id": "ws_abc123",
  "created_at": "2026-02-09T17:34:00Z"
}
```

#### `GET /workspace/{id}`
Get workspace state and entity list.

#### `POST /workspace/{id}/save`
Save workspace to package file.

#### `POST /workspace/load`
Load workspace from package file or URL.

---

### DSL Execution

#### `POST /dsl/execute`
Execute yapCAD DSL commands and return geometry + views.

```json
Request:
{
  "workspace_id": "ws_abc123",
  "commands": [
    "cylinder(r=50, h=100)",
    "translate([0, 0, 50])",
    "difference(box([120, 120, 20]))"
  ],
  "view_options": {
    "plan": true,
    "elevation": true,
    "perspective": true,
    "perspective_angles": [45, 30],
    "resolution": [800, 600]
  }
}

Response:
{
  "success": true,
  "entities_created": ["ent_001", "ent_002"],
  "views": {
    "plan": "data:image/png;base64,...",
    "front_elevation": "data:image/png;base64,...",
    "side_elevation": "data:image/png;base64,...",
    "perspective": "data:image/png;base64,..."
  },
  "bounds": {
    "min": [-60, -60, 0],
    "max": [60, 60, 100]
  },
  "entity_count": 2,
  "history_index": 3
}
```

#### `POST /dsl/validate`
Validate DSL commands without executing (syntax check, constraint check).

```json
Request:
{
  "commands": ["cylinder(r=50, h=100)", "bogus_command()"]
}

Response:
{
  "valid": false,
  "errors": [
    {"line": 2, "message": "Unknown command: bogus_command"}
  ]
}
```

---

### Entity Operations

#### `GET /entities/{workspace_id}`
List all entities in workspace.

#### `GET /entity/{id}`
Get entity details (geometry type, bounds, properties).

#### `POST /entity/{id}/transform`
Apply transform to entity.

```json
Request:
{
  "translate": [10, 0, 0],
  "rotate": [0, 0, 45],
  "scale": 1.0
}
```

#### `DELETE /entity/{id}`
Remove entity from workspace.

#### `POST /entity/{id}/copy`
Duplicate entity with optional transform.

---

### Assembly Operations

#### `POST /assembly/create`
Create named assembly from entities.

```json
Request:
{
  "workspace_id": "ws_abc123",
  "name": "thrust_structure",
  "entities": ["ent_001", "ent_002", "ent_003"],
  "origin": [0, 0, 0]
}
```

#### `POST /assembly/{id}/add`
Add entity to existing assembly.

#### `POST /assembly/{id}/instance`
Create instance of assembly with transform.

---

### Visualization

#### `POST /view/render`
Render current workspace state to images.

```json
Request:
{
  "workspace_id": "ws_abc123",
  "views": ["plan", "perspective"],
  "highlight_entities": ["ent_001"],
  "section_plane": {"origin": [0, 0, 50], "normal": [0, 0, 1]},
  "resolution": [1024, 768],
  "format": "png"
}

Response:
{
  "views": {
    "plan": "data:image/png;base64,...",
    "perspective": "data:image/png;base64,...",
    "section": "data:image/png;base64,..."
  }
}
```

#### `POST /view/animate`
Generate animation frames (for agent to understand motion/mechanism).

```json
Request:
{
  "workspace_id": "ws_abc123",
  "animation": {
    "type": "rotate",
    "entity": "ent_001",
    "axis": [0, 0, 1],
    "range": [0, 360],
    "frames": 8
  }
}
```

---

### Analysis Integration

#### `POST /analysis/bounds`
Get bounding box, volume, surface area, center of mass.

```json
Response:
{
  "bounds": {"min": [...], "max": [...]},
  "volume_mm3": 125000,
  "surface_area_mm2": 15000,
  "center_of_mass": [0, 0, 50]
}
```

#### `POST /analysis/interference`
Check for interference between entities/assemblies.

```json
Request:
{
  "workspace_id": "ws_abc123",
  "entity_a": "ent_001",
  "entity_b": "ent_002"
}

Response:
{
  "interference": true,
  "overlap_volume_mm3": 1250,
  "overlap_bounds": {...}
}
```

#### `POST /analysis/fea`
Submit geometry for FEA analysis (async).

```json
Request:
{
  "workspace_id": "ws_abc123",
  "entity": "ent_001",
  "analysis_type": "static_stress",
  "loads": [
    {"type": "force", "location": [0, 0, 100], "vector": [0, 0, -500], "units": "N"}
  ],
  "constraints": [
    {"type": "fixed", "face": "bottom"}
  ],
  "material": {
    "name": "PC-PBT-CF",
    "E_MPa": 8000,
    "poisson": 0.35,
    "yield_MPa": 80
  }
}

Response:
{
  "job_id": "fea_xyz789",
  "status": "queued",
  "estimated_time_s": 120
}
```

#### `GET /analysis/fea/{job_id}`
Poll FEA job status and retrieve results.

```json
Response:
{
  "status": "complete",
  "results": {
    "max_stress_MPa": 45.2,
    "max_displacement_mm": 0.12,
    "safety_factor": 1.77,
    "stress_plot": "data:image/png;base64,..."
  }
}
```

#### `POST /analysis/print`
Analyze printability (overhangs, support requirements, orientation).

```json
Response:
{
  "printable": true,
  "recommended_orientation": [0, 0, 0],
  "overhang_volume_pct": 3.2,
  "support_volume_mm3": 1200,
  "warnings": ["Thin wall at z=45mm, consider 2mm minimum"]
}
```

---

### Package Management

#### `POST /package/create`
Create package from workspace.

```json
Request:
{
  "workspace_id": "ws_abc123",
  "manifest": {
    "name": "Thrust Structure v1",
    "uuid": "550e8400-e29b-41d4-a716-446655440000",
    "version": "1.0.0",
    "author": "jarvis",
    "description": "Agentic-1 thrust structure"
  },
  "include_views": true
}
```

#### `GET /package/{uuid}/versions`
List all versions of a package.

#### `POST /package/import`
Import package from file, URL, or GEL repository.

---

### GEL Integration

#### `POST /gel/explore`
Start generative exploration from seed design.

```json
Request:
{
  "seed_workspace": "ws_abc123",
  "parameters": [
    {"name": "wall_thickness", "min": 2, "max": 5, "step": 0.5},
    {"name": "rib_count", "min": 4, "max": 12, "step": 2}
  ],
  "constraints": [
    {"type": "max_mass_g", "value": 500},
    {"type": "min_safety_factor", "value": 2.0}
  ],
  "objective": "minimize_mass",
  "max_designs": 50
}

Response:
{
  "exploration_id": "gel_exp_001",
  "status": "running",
  "designs_generated": 0
}
```

#### `GET /gel/{exploration_id}/designs`
Get generated design family.

```json
Response:
{
  "designs": [
    {
      "design_id": "des_001",
      "parameters": {"wall_thickness": 3.0, "rib_count": 8},
      "metrics": {"mass_g": 423, "safety_factor": 2.3},
      "views": {...},
      "pareto_optimal": true
    },
    ...
  ]
}
```

#### `POST /gel/{exploration_id}/select`
Select design(s) from family for further development.

---

## Error Handling

All errors return structured responses:

```json
{
  "success": false,
  "error": {
    "code": "DSL_SYNTAX_ERROR",
    "message": "Unexpected token at line 3",
    "details": {...}
  }
}
```

Error codes:
- `DSL_SYNTAX_ERROR` — Invalid DSL syntax
- `DSL_RUNTIME_ERROR` — DSL execution failed
- `CONSTRAINT_VIOLATION` — Design constraint not met
- `GEOMETRY_INVALID` — Invalid geometry (self-intersecting, etc.)
- `WORKSPACE_NOT_FOUND` — Unknown workspace ID
- `ANALYSIS_FAILED` — FEA or other analysis failed

---

## Authentication

For multi-user/remote deployments:
- Bearer token in `Authorization` header
- Workspace ownership and permissions
- Rate limiting per agent/user

For local single-user (initial deployment):
- No auth required
- Localhost binding only

---

## WebSocket Stream (Optional)

For real-time GEL visualization:

```
WS /stream/{workspace_id}

Messages:
← {"type": "entity_added", "entity": {...}, "views": {...}}
← {"type": "analysis_complete", "job_id": "...", "results": {...}}
→ {"type": "subscribe", "events": ["entity_added", "gel_design"]}
```

---

## Implementation Notes

### Phase 1 — Core API
- Workspace CRUD
- DSL execution with view rendering
- Basic entity operations
- Package save/load

### Phase 2 — Analysis
- Bounds, volume, interference
- Printability analysis
- FEA integration (external solver)

### Phase 3 — GEL
- Parameter exploration
- Constraint evaluation
- Design family management
- Pareto frontier tracking

### Tech Stack
- **Server:** Python + FastAPI
- **Geometry:** yapCAD core
- **Rendering:** Headless Three.js or Cairo for 2D views
- **FEA:** External solver (CalculiX, or cloud API)

---

## Example Agent Session

```
Agent: Create a cylinder with radius 50mm, height 100mm
→ POST /dsl/execute {"commands": ["cylinder(r=50, h=100)"]}
← {views: {plan: ..., perspective: ...}, entities: ["ent_001"]}

Agent: [sees images] Add a 10mm fillet to the top edge
→ POST /dsl/execute {"commands": ["fillet(ent_001, edge='top', r=10)"]}
← {views: {...}}

Agent: Check if this fits in a 320mm print bed
→ POST /analysis/bounds
← {bounds: {max: [50, 50, 100]}} — fits!

Agent: Run stress analysis with 500N downward force
→ POST /analysis/fea {loads: [...], material: "PC-PBT-CF"}
← {job_id: "fea_001", status: "queued"}

Agent: [polls until complete]
← {max_stress: 45MPa, safety_factor: 1.77, stress_plot: ...}

Agent: Safety factor too low. Explore wall thickness 3-6mm
→ POST /gel/explore {parameters: [{name: "wall", min: 3, max: 6}], ...}
← {exploration_id: "gel_001"}

Agent: [polls, reviews design family, selects optimal]
```

---

*Draft: 2026-02-09 — Jarvis*
