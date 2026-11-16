# tools/

One-off and diagnostic scripts around yapCAD. Library functionality lives
in `src/yapcad`; if a tool grows into something reusable, move it into the
package and leave a thin wrapper (or nothing) here.

## The recommended SDF workflow

For designs built as fields (`--representation sdf`), the pipeline is
**check → preview → mesh once**:

```bash
# 1. static check -- types, imports, native blocks
python -m yapcad.dsl check design.dsl

# 2. previews straight from the field (no meshing at all)
python -m yapcad.render design.dsl PART --view iso -o preview.png
python -m yapcad.render design.dsl PART --view 4up -o sheet.png

# 3. one mesh, at the end, reported for printability
python -m yapcad.dsl run design.dsl PART --representation sdf \
    --sdf-cell 0.2 -o part.stl --require-watertight
```

The run reports triangles, watertightness, connected components and
bounding box on stderr. The same check is available to scripts as
`yapcad.mesh_check.check_mesh` (it accepts a yapCAD solid, a trimesh
mesh, or a mesh file path).

## Tools

| Tool | Status | What it does |
|------|--------|--------------|
| `render_assembly.py` | current | Multi-view PNG of a per-layer DSL assembly. Thin wrapper over `yapcad.render`: builds each layer as a deferred SDF field, one coloured part per layer, 4-up sheet (`--views 4`, default) or one image per named view (`--views 6`). |
| `validate_mesh.py` | current | CAM-oriented lint/heal/slice pipeline over an STL. Starts with the built-in `yapcad.mesh_check` report, then chains `admesh`, `meshfix` and a slicer front-end if installed. |
| `mesh_diagnostics.py` | current | Shared trimesh diagnostics (`MeshDiagnostics`) used by the diagnose tools; watertightness and component counts come from `yapcad.mesh_check`. |
| `diagnose_boolean.py` | current | Boolean-engine diagnostics (cube/sphere/box-hole/tube fixtures) through trimesh. Mesh representation only. |
| `diagnose_primitive.py` | current | Primitive-engine diagnostics through trimesh. Mesh representation only. |
| `analyze_vertices.py` | current | Orphan/duplicate vertex analysis on yapCAD's internal mesh representation. Deferred solids must be `yapcad.sdf.realize`d first. |
| `benchmark_dsl_booleans.py` | current | Times DSL boolean builds; used for the SDF meshing performance notes. |
| `find_func.py` | current | Greps DSL sources for a command definition. |
| `ycpkg_analyze.py`, `ycpkg_export.py`, `ycpkg_validate.py`, `ycpkg_viewer.py` | current | Thin wrappers over `yapcad.package` for inspecting `.ycpkg` packages from the shell. |
| `reference/` | reference | Archived visualization helper scripts (the ancestors of `render_assembly.py`); kept for comparison, not run by anything. |

## SDF notes

- Field-first solids carry their tree and a *plan* to mesh, not
  triangles. Anything here that counts vertices or facets is
  mesh-representation only; realize first.
- `python -m yapcad.dsl run --defer-mesh` is the default for
  `--representation sdf` + STL output; `--no-defer-mesh` restores eager
  per-primitive meshing.
