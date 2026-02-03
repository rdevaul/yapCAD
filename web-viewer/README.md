# yapCAD Web Viewer

Web-based viewer for yapCAD geometry packages (`.ycpkg` files).

## Architecture

- **Frontend**: TypeScript + Three.js loader (`src/yapcad-loader.ts`)
- **Backend**: FastAPI server serving `.ycpkg` packages (`src/backend/`)

## Quick Start

### Backend

```bash
cd web-viewer
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
./run-backend.sh
```

Backend runs on `http://localhost:8000`

### Frontend

```bash
npm install
npm run dev
```

Opens viewer at `http://localhost:5173` (or similar Vite port)

## Package Format

A `.ycpkg` file is a ZIP archive containing:
- `manifest.json` — Package metadata
- `geometry/primary.json` — Main geometry (`yapcad-geometry-json-v0.1`)
- `geometry/*.json` — Optional additional geometry files
- `materials/*.json` — Optional material definitions

## API Endpoints

- `GET /packages` — List available packages
- `GET /packages/{id}` — Get package manifest
- `GET /packages/{id}/geometry/{file}` — Get geometry JSON
- `GET /packages/{id}/thumbnail` — Get package thumbnail

## Status

Early development. Basic Three.js rendering working with test geometry.
