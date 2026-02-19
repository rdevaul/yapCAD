# yapCAD FastAPI Service

FastAPI backend exposing yapCAD DSL execution and geometry operations for the yapCAD Web Viewer and agent workflows.

## Run (dev)

From the yapCAD repo root:

```bash
cd ~/Projects/yapCAD
PYTHONPATH=src conda run -n yapcad-brep \
  env YAPCAD_BOOLEAN_ENGINE=trimesh YAPCAD_TRIMESH_BACKEND=blender \
  uvicorn service.main:app --reload --port 8000
```

CORS is enabled for `http://localhost:5173`.

## Endpoints

- `GET /health`
- `GET /info`
- `GET /builtins`

- `POST /dsl/parse` `{ "source": "..." }`
- `POST /dsl/eval` `{ "source": "...", "command": "MyCmd", "parameters": {}, "format": "json" }`

- `POST /tessellate` `{ "geometry": <geometry-json-doc> }`
- `POST /boolean` `{ "operation": "union"|"difference"|"intersection", "a": <doc>, "b": <doc>, "engine": "trimesh" }`

## Notes

- DSL evaluation is executed in a worker thread and timed out after 30 seconds.
- Geometry JSON schema is `yapcad-geometry-json-v0.1` (see `yapcad.io.geometry_json`).
