#!/bin/bash
# Run the yapCAD service backend (from monorepo root)
# This script must be run from the yapCAD repo root, or the workbench/ subdirectory.
# The backend lives at service/main.py in the yapCAD repo root.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

cd "$REPO_ROOT"

conda run -n yapcad-brep env \
  PYTHONPATH=src \
  YAPCAD_BOOLEAN_ENGINE=occ \
  YAPCAD_MESH_BOOLEAN_ENGINE=trimesh:blender \
  uvicorn service.main:app --host 0.0.0.0 --port 8000
