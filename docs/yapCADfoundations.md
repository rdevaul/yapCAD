# yapCAD Shared Foundations

## Purpose
Establish common geometry, metadata, and validation capabilities required to support near-term exporters (STEP, STL) and longer-term design provenance goals while staying compatible with yapCAD 0.x.

## Scope
- Mesh representation utilities
- Surface/solid metadata extensions
- Validation helpers
- Export support primitives
- Testing expectations

## Geometry Canonicalisation
- Provide a `mesh_view(surface)` helper returning a sequence of `(normal, v0, v1, v2)` tuples with consistent right-handed winding and recalculated normals when missing.
- Expose vertex/triangle iterators that accept both single surfaces and solids (flattening shell surfaces where safe) with filtering for duplicates or zero-area tris.
- Standardise conversions from homogeneous 4D points to Cartesian XYZ arrays; define utility to strip/reapply the w-component.

## Metadata Extensions
- Attach optional UUIDs to faces, surfaces, and solids to support traceability through preprocessing (triangulation, boolean ops, export).
- Allow storage of origin identifiers (parametric expression/source file, timestamp) at each surface/solid level.
- Record unit system and default coordinate frame for each solid; exporters fallback to project defaults when absent.
- Introduce lightweight `DesignContext` object capturing tool versions (yapCAD version, dependencies) and optional author info; embed references in surfaces/solids without forcing persistence changes yet.

## Validation Helpers
- Add reusable checks:
  - `is_closed_polygon(points, tol)` to assert loop continuity.
  - `faces_oriented(surface)` verifying triangle orientation matches stored normals.
  - `surface_watertight(surface)` heuristic (matching boundary loops) usable before export.
- Provide warning utilities that exporters/test harnesses can consume without raising exceptions unless requested.

## Export Support Primitives
- Use mesh view helpers inside new STL and STEP writers.
- Provide shared normal computation and vertex deduplication code so exporters stay DRY.
- Define serialization-ready data classes for triangles and loops; ensures binary/ASCII writers use same structures.

## Testing Expectations
- Unit tests covering mesh view normalisation (degenerate triangles, normals missing, reversed loops).
- Golden tests for mesh flattening to guarantee determinism (ordering, orientation).
- Exporters will re-use these helpers, so failures surface early.

## Backwards Compatibility
- Keep current surface/solid structures intact; metadata added via optional fields/defaults.
- Provide migration utilities when saving/loading legacy surfaces lacking UUID/context.

