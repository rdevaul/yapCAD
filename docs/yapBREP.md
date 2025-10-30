# yapCAD BREP Roadmap

This note captures a concrete plan for extending yapCAD’s native data structures so
we can ingest/export analytic STEP models without collapsing them to triangle meshes.
The intent is to evolve yapCAD into a lightweight BREP kernel while keeping today’s
list-based API backwards compatible.

---

## 1. Current State

* **Curves (2D)** – first-class line, arc, circle, Catmull–Rom, NURBS primitives; polygons
  and polylines are stored as lists of points. There is no explicit ellipse/conic type yet.
* **Surfaces (3D)** – planes, spheres, cones, cylinders, prisms etc. exist as “generators”
  that produce tessellated surfaces (`['surface', verts, normals, faces, …]`). Once
  tessellated, the analytic definition is lost; the solid stores metadata about the
  procedure that produced it.
* **Solids** – list-based shells referencing surfaces, plus metadata. No explicit vertex/edge
  topology is stored; the mesh is the only representation.
* **Serialization** – geometry JSON stores solids/surfaces as triangle soup, and sketches as
  polylines + primitives. STEP export is faceted only; STEP import is not implemented.

---

## 2. Goals

1. **Retain analytic fidelity** – curves and surfaces imported from STEP should remain
   parametric (no forced tessellation), with tessellation performed lazily when required.
2. **Explicit topology** – introduce vertices, edges, trims, loops, faces, shells so we can
   mirror STEP-style BREPs.
3. **Backwards compatibility** – existing list-based APIs (`geom*`, `geometry_json`, DSP
   packages) should continue to work; new structures should integrate without breaking them.
4. **Extensible serialization** – `.ycpkg` geometry JSON must round-trip analytic curves and
   surfaces; tessellations become derived assets.
5. **Interoperable export** – STEP export should regenerate analytic entities instead of
   only meshed facets; DXF export for sketches already does this for curves.

---

## 3. Core Data Structure Extensions

### 3.1 Curves (2D/3D)

* **New curve primitives**
  * `ellipse(center, major, minor, axis, span)` / `isellipse`.
  * `conic` for generic quadrics (parabola, hyperbola) – store focal parameters + orientation.
  * 3D curves mirror 2D definitions but with explicit plane/axis metadata, matching STEP
    `CIRCLE`, `ELLIPSE`, `B_SPLINE_CURVE`, `CONIC` etc.
* **Parameter ranges** – store `[u0, u1]` on curve instances when they are used as edges.
* **Tolerance metadata** – attach `curve_meta` with modeling precision, so downstream ops can
  compare within a consistent epsilon.

### 3.2 Surfaces

* Create first-class `planesurface`, `cylindersurface`, `conesurface`, `spheresurface`,
  `tori surface`, and generic `nurbs_surface`. Each stores its parameter domain (`u,v`
  bounds), local axis, and control net (for splines).
* Maintain a tessellation cache per surface (e.g., `surface[TRI_CACHE]`) plus refinement
  options.
* Provide evaluation utilities (`evaluate_surface`, `surface_normal`, `surface_uv_project`)
  to support trimming, intersections, etc.

### 3.3 Topology Graph

Introduce new list-based constructs (or light classes) layered above the existing surface/
solid representation:

| Entity  | Responsibilities                                                            |
|---------|------------------------------------------------------------------------------|
| `brep_vertex` | world-space point + tolerance + incident edges references                |
| `brep_edge`   | references a curve primitive + parameter range + start/end vertices;    |
|               | stores sense (forward/reverse).                                         |
| `brep_trim`   | oriented edge with UV parameterisation for a specific face              |
| `brep_loop`   | ordered list of trims forming outer/inner boundary of a face            |
| `brep_face`   | references a surface primitive + loops + natural orientation            |
| `brep_shell`  | set of faces forming a closed shell (outer or inner)                    |
| `brep_solid`  | collection of shells (outer + optional inner voids)                     |

Each object keeps metadata (`layer`, `tags`, creation history) just like current solids.
Tessellations are derived on demand from these definitions.

---

## 4. Lazy Tessellation Strategy

* Provide a `tessellate(surface, quality)` API that converts analytic surfaces to the
  current triangle format.
* Modify `geometry_to_json` so solids with analytic faces are serialized as analytic
  definitions + optional tessellation cache.
* The viewer/exporters can request tessellation when needed (e.g., to render meshes or
  export STL). For STEP export we use the analytic definitions directly.

---

## 5. Serialization (geometry JSON & STEP)

1. **JSON schema updates**
   * Extend sketch `primitives` to include ellipse/conic entries (already storing lines,
     arcs, splines).
   * Introduce new entity types:
     * `"curve"` (analytic definition + parameters).
     * `"surface"` (existing tessellated surfaces) + `"analyticSurface"` (parametric).
     * `"brep"` object that carries vertices/edges/trims/faces/shells with references.
   * Keep backwards-compatible `polylines` and tessellated surfaces for consumers that
     still expect meshes.
2. **STEP export**
   * Refactor exporter to walk the new BREP structure and emit STEP analytic entities
     (plane, cylinder, B-spline surface, etc.) and topological relationships.
   * For solids that only have tessellations, fall back to faceted BREP (current behaviour).
3. **STEP import**
   * Parse analytic geometry into the new primitives; tessellate lazily for visualization.
   * Support at least the common STEP AP203/AP214 entities used in CAD: planes, cylinders,
     cones, spheres, torus, b-spline surfaces, trimmed surfaces, edges with analytic curves.

---

## 6. API & Package Touch Points

* `geom.py` – add ellipse/conic primitives; extend sampling, length, bbox, intersection
  utilities to understand them.
* `geom3d.py` – new surface constructors; evaluation/tessellation helpers; solids now hold
  `brep_solid` references in metadata.
* `geom3d_util.py` – extrusion/loft/tube functions should emit both analytic surfaces and
  BREP topology (edges, trims) in addition to tessellated surfaces.
* `geometry_json.py` – encode/decode new primitive types; maintain compatibility with
  existing consumers.
* `package/viewer.py` – when given analytic surfaces, tessellate on the fly for display.
* `tools/ycpkg_export.py` – allow exporting either analytic STEP or derived STL/STP meshes.

---

## 7. Implementation Phases

1. **Foundations**
   * Add ellipse/conic primitives (+ tests).
   * Add analytic surface structures & evaluation utilities.
   * Introduce BREP topology containers linked to existing solids.
2. **Serialization**
   * Update geometry JSON reader/writer to store primitives & BREP graph.
   * Ensure `.ycpkg` round-tripping covers new entity types (unit tests).
3. **Tessellation refactor**
   * Centralize tessellation (quality settings, caching).
   * Viewer/export updates to request tessellations lazily.
4. **STEP import/export**
   * Build analytic STEP export pipeline (faces -> trimmed surfaces -> topology).
   * Implement STEP reader that populates the new structures.
5. **Tooling & validation**
   * Extend `ycpkg_export.py` to offer `--formats step-analytic`, `--formats stl`.
   * Add diagnostic tools for watertightness, trimming, topology consistency.

---

## 8. Risks & Mitigations

* **Tolerance management** – adopt a consistent global modeling tolerance, propagate to
  curves/surfaces/vertices. Provide helpers for fuzzy comparisons.
* **Performance** – analytic evaluation and tessellation must be efficient; consider
  caching compiled spline representations (e.g., via NumPy) or delegating to OCC in the
  future.
* **Complex STEP entities** – initial scope should focus on the 80/20 set (planes,
  cylinders, cones, spheres, torus, NURBS). More exotic surfaces (offsets, swept,
  intersection curves) can be mapped via approximation or deferred.
* **Backwards compatibility** – ensure old packages with only tessellated data still load,
  and new packages degrade gracefully when primitives are missing.

---

## 9. Next Steps

1. Prototype ellipse/conic primitives in `geom.py` + tests.
2. Draft analytic surface objects (plane/cylinder/cone) and wire them into extrusion/
   revolution helpers.
3. Define minimal BREP data structures (vertex, edge, face) in a new module (e.g.,
   `yapcad.brep`) and start storing them alongside solids.
4. Update geometry JSON to carry `primitives` and BREP topology; add regression tests.
5. Build a simple STEP analytic exporter to validate the data model.

This staged approach lets us evolve yapCAD incrementally: first by capturing analytic
information during geometry creation, then by enhancing serialization and finally by
supporting full STEP import/export. Once in place, STL import can coexist (as tessellated
geometry), while STEP models benefit from a richer, faithful representation.
