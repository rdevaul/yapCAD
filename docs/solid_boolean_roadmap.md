# Boolean Engine Roadmap

Status snapshot as of this session:

- Native mesh boolean workflow now lives in `yapcad/boolean/native.py`, with `geom3d` re-exporting it for backward compatibility.
- Diagnostics still rely on `tools/validate_mesh.py`; no engine-specific metrics yet.
- External engine integration is underway; a `trimesh` backend is plumbed (requires an installed boolean operator such as Blender/OpenSCAD/Cork).

Near-term tasks:

- [x] **Code extraction** – current boolean implementation now resides in `yapcad/boolean/native.py`, referenced by `geom3d` wrappers.
- [x] **Engine selector UX** – `solid_boolean(..., engine=...)` now routes through the native engine by default, accepts `trimesh` (and optional `trimesh:backend`), and the demo exposes `--engine`; document env flags (`YAPCAD_BOOLEAN_ENGINE`, `YAPCAD_TRIMESH_BACKEND`) for benchmarking.
- [ ] **External prototype** – wrap at least one library boolean (e.g., `trimesh` or PyMeshFix) for STL solids, including geometry conversion helpers.
- [ ] **Benchmark harness** – update the demo CLI and validation scripts to iterate across registered engines, saving STL/lint outputs for comparison.

Open questions:
- Which external kernel offers the best balance of licensing, install footprint, and mesh quality?
- Do we keep the native engine on equal footing (for offline/embedded use), or treat it as a fallback once a third-party backend is available?
- What minimal metrics should every engine report (shell count, single-facet edges, `validate_mesh` scores, render snapshots)?

Once the refactor lands, future work will include:
- Automatic tolerance scaling per engine.
- Capturing validation output in machine-readable form (JSON) to feed regression dashboards.
- Evaluating hybrid workflows (native preprocessing, external boolean, native stitches/cleanup).
