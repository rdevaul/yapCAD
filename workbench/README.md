# yapCAD Workbench V2

Interactive web-based CAD workbench for the yapCAD DSL — live 2D sketch editing,
3D solid preview, and parametric design with command composition.

## Architecture

| Layer | Tech | Notes |
|-------|------|-------|
| Frontend | TypeScript + React + Vite | Port 5173 (dev) / 5174 (Caddy proxy) |
| Backend | FastAPI + yapCAD (OCC) | Port 8000 / 8001 (Caddy proxy) |
| Geometry | pythonocc-core, trimesh | BREP booleans + tessellation |
| Comms | WebSocket `/ws/dsl/eval` | REST `/dsl/eval` fallback |

## Quick Start

### Backend (yapCAD service)

The service runs as a LaunchAgent on macOS:

```bash
launchctl load ~/Library/LaunchAgents/ai.dml.yapcad-service.plist
```

Manual start (dev):
```bash
cd ~/Projects/yapCAD
conda run -n yapcad-brep env PYTHONPATH=src \
  YAPCAD_BOOLEAN_ENGINE=occ \
  YAPCAD_TRIMESH_ENGINE=blender \
  python -m uvicorn service.main:app --host 0.0.0.0 --port 8000
```

Logs: `/tmp/yapcad-service.log`

### Frontend

```bash
cd ~/Projects/yapcad-workbench
npm install
npm run dev   # http://localhost:5173
```

### Caddy reverse proxy (Tailscale access)

```bash
caddy run --config Caddyfile
# Workbench: http://100.83.253.120:5174/
# API docs:  http://100.83.253.120:8001/docs
```

## Key Features

### Multi-tab DSL editor
- Load DSL files from the file bar or paste/type in the editor
- Multiple files open simultaneously; each tab has independent state

### Parameter Panel
- Auto-populated from DSL `command` signatures via `/dsl/commands`
- Type-driven widgets: `float` sliders, `int` inputs, `bool` toggles, `point2d` coordinate pairs
- `@ui(widget="point2d", label="...")` annotations enable interactive drag handles
- **▶ Evaluate** — sends current params to backend, renders result
- **✦ Set defaults** — bakes current parameter values into the DSL source as new default expressions

### 2D Sketch viewer
- Live Catmull-Rom spline rendering (centripetal, Barry-Goldman)
- Draggable `point2d` control point handles
- Auto-fit bbox — canvas always centres and scales to content
- Switching to a command with 3D output auto-switches to the 3D viewer

### 3D viewer
- Three.js with BREP tessellation (OCC mesh, deflection 0.3)
- Orbit / pan / zoom controls

### Command composition
DSL commands can call other commands by name with no arguments, inheriting
their current defaults. This is the canonical pattern for dependent geometry:

```dsl
command lathe_solid() -> solid:
    let profile: catmullrom = spline_profile()   # calls sibling command
    let region: region2d = region_from_spline(profile)
    let result: solid = revolve(region, axis, 360.0)
    emit result
```

**Workflow:**
1. Evaluate `spline_profile` → drag control points to shape the profile
2. Click **✦ Set defaults** → `point2d(x, y)` defaults rewritten in source
3. Switch to `lathe_solid` → Evaluate → solid uses updated profile automatically

### WebSocket eval with latency tiers
- First eval latency auto-classifies command as fast / medium / slow
- Fast commands re-eval on every param change (debounced 300 ms)
- Medium/slow commands require explicit Evaluate press

## DSL Conventions

### `@ui` decorator syntax
```dsl
name: type @ui(widget="point2d", label="My Point") = point2d(0.0, 0.0)
```
Goes **after** the type annotation, **before** the `=` default.

### `point2d` coordinate convention for lathe/revolve
- **X** = radius from the rotation axis (must be ≥ 0)
- **Y** = height along the rotation axis

### `revolve()` builtin
```dsl
let region: region2d = region_from_spline(curve)
let axis:   vector   = vector(0.0, 0.0, 1.0)
let solid:  solid    = revolve(region, axis, 360.0)
```
Region is in the XY plane; `revolve` maps X→radius, Y→Z for OCC.

## API Endpoints

| Method | Path | Description |
|--------|------|-------------|
| POST | `/dsl/eval` | Evaluate a command, return geometry JSON |
| POST | `/dsl/commands` | Parse DSL, return command/param signatures |
| POST | `/dsl/callgraph` | Return call graph for a command |
| POST | `/dsl/ui_eval` | Eval with `@ui` param hints |
| WS | `/ws/dsl/eval` | Streaming WebSocket eval |

## Demo Files

`public/spline_demo.dsl` — interactive lathe demo:
- `spline_profile` — 5 draggable `point2d` control points, renders 2D Catmull-Rom curve
- `lathe_solid` — calls `spline_profile()` with no args, revolves profile → solid of revolution

## Key Commits

| Hash | Description |
|------|-------------|
| `46c8b54` (yapCAD) | Phase 1: WebSocket eval, callgraph, `@ui` decorators, scalar emit |
| `b7b4db3` | Phase 2: multi-tab UI, `useWsEval` hook, ParameterPanel |
| `7cdb51b` | `point2d` param type + interactive drag |
| `cee4490` | Client-side Catmull-Rom sampler |
| `6ed5c9d` | Catmull-Rom Barry-Goldman fix |
| `dd04848` (yapCAD) | `revolve()` rewrite — OCC BRepPrimAPI_MakeRevol |
| `26b7294` (yapCAD) | BREP-only solid tessellation on serialize |
| `9cc4059` | Drag hit-test bboxRef fix |
| `4b3b09c` | Lazy eval model + `batchUpdateParams` |
| `cc2fab9` | ✦ Set defaults button + command composition DSL |
| `ad58dac` | Set defaults → editor live update via `setExternalSource` |
