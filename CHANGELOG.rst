=========
Changelog
=========

Unreleased
==========

Signed distance functions
-------------------------

- New ``yapcad.sdf`` package: signed distance functions as a peer solid
  representation (see ``docs/SDF-DESIGN.md``). Fields are an immutable,
  serialisable node DAG of analytic primitives, hard and smooth booleans,
  offsets, shells, transforms and gyroid/Schwarz-P lattices, evaluated by a
  vectorised numpy backend that tracks each field's exactness and Lipschitz
  bound. ``sdf.to_solid`` dual-contours a field into an ordinary yapCAD solid
  that the rest of yapCAD can consume, deterministically and with sharp edges
  preserved; ``to_solid(..., brep=True)`` replays CSG-expressible trees
  through OpenCASCADE for an exact BREP and analytic STEP export. The package
  needs only numpy; OCC replay is optional. ``examples/sdf_demo.py`` meshes a
  gallery of models to STL, PNG and (with ``--step``) STEP.
- Booleans between two SDF-authored solids now combine their fields instead
  of their meshes: ``geom3d.solid_boolean`` -- and so the DSL's ``union``,
  ``difference`` and ``intersection`` -- returns an exact, closed,
  SDF-authoritative result in well under a second, where the mesh engine
  returned an open mesh, took tens of seconds, or failed. The result is
  meshed no coarser than either operand. An explicit ``engine=`` or
  ``YAPCAD_BOOLEAN_ENGINE`` still takes precedence, and ``engine="sdf"``
  selects the field path. Mixed SDF/mesh or SDF/BREP operands are unchanged.
- ``fillet`` works on SDF-authored solids without OCC. ``sdf.fillet``
  rounds every edge exactly for the analytic primitives -- a box becomes a
  ``rounded_box``, a cylinder the new exact ``rounded_cylinder`` -- through
  similarity transforms and across the bodies of a compound, and the result
  still replays to analytic STEP. The DSL's ``fillet`` uses it for any
  SDF-authored solid. Fillets of combined fields, through non-uniform
  scales, and of mesh-only solids are refused with an explanation; the plan
  for them is in ``docs/SDF-DESIGN.md`` §10.5.
- ``compound`` keeps the fields of SDF-authored operands, as a new
  ``sdf.compound`` node: a union for booleans, separate bodies for meshing
  and OCC replay. Booleans against a compound no longer fall back to the
  mesh engine.
- Straight bevel and miter gears build without OCC. The new
  ``straight_bevel_gear`` field stores the gear spec and evaluates as the
  BREP generator's own construction, the outer tooth section carried to
  the pitch apex; it still replays to that generator's exact BREP. The DSL's
  ``miter_gear`` and ``straight_bevel_gear`` use it when OCC is absent, and
  ``gears.make_straight_bevel_gear_sdf`` exposes it directly. New planar
  field kinds back it: ``polygon`` (exact, with a symmetry fold),
  ``extrude`` and ``apex_extrude``.
- Spur, helical and herringbone gears as fields: ``sdf.spur_gear`` (the
  ``figgear`` profile, extruded, and twisted by the new ``sdf.twist``
  operator) and ``gears.make_involute_gear_sdf``. The DSL's
  ``herringbone_gear`` uses it when OCC is absent.
- ``metric_hex_nut`` and ``unified_hex_nut`` take ``representation="sdf"``
  for a closed, SDF-authoritative nut with a helical ``thread`` field whose
  profile is ``threadgen``'s own. The swept mesh, which is not closed,
  remains the default.
- ``metric_hex_bolt`` and ``unified_hex_bolt`` take ``representation="sdf"``
  for a closed bolt with an external ``thread`` field (``sdf.hex_bolt``).
- DSL designs can be built as SDF fields: ``compile_and_run`` and
  ``package_from_dsl`` take ``representation="sdf"`` (also
  ``yapcad.dsl run --representation sdf`` and
  ``YAPCAD_DSL_REPRESENTATION=sdf``). The solid primitives, fasteners and
  gears then produce SDF-authored solids placed as their mesh versions are,
  so an existing design builds as fields, without OCC, unchanged.
- ``sdf.simplify_solid`` simplifies a meshed field by quadric edge collapse,
  checking every collapse against the field: a closed manifold result, no
  flipped triangles, and the field's value on every changed triangle within
  the tolerance, or no worse than the uniform mesh already was. At 0.01 mm
  a threaded M8 nut keeps 4% of its triangles and a miter gear 6%, both
  with lower 99th-percentile error than the uniform mesh. The tolerance is
  recorded with the meshing parameters.
- Dual contouring places one vertex per surface component of a cell
  (Manifold Dual Contouring), so acute edges -- helical gear tips, lattice
  walls -- no longer mesh non-manifold. Ordinary meshes are unchanged.
- Symmetric ``polygon`` fields evaluate three times faster.
- ``geom3d`` solid transforms now carry an SDF-authored solid's tree along
  with its mesh; previously ``translatesolid`` and friends moved the mesh and
  left the defining tree where it was.
- ``sdf.to_solid`` raises ``sdf.NonManifoldMeshError`` (an ``SdfError``
  subclass) when it refuses a mesh, so callers can catch that case alone.

Mesh booleans
-------------

- The native mesh boolean engine is rewritten (``yapcad.boolean.csg``) and
  now returns closed, volume-correct solids: 24 of 24 benchmark cases (was
  8) and 1,800 randomised boxes/spheres/cylinders/icosahedra cases against a
  manifold3d reference, where the old engine left holes, dropped slivers or
  timed out. Triangles are split by the planes of the triangles they cross,
  fragments are classified by generalized winding number (coplanar faces by
  orientation), each kept triangle is re-triangulated without losing a
  vertex, and the result is welded and T-junction-repaired into a
  conforming mesh. It is pure numpy and deterministic; a 10k-triangle sphere
  union takes about two seconds. ``native.solid_boolean`` keeps its
  signature; ``tol`` and ``stitch`` are accepted and ignored.
- The native engine no longer shatters triangles. Cuts now require a real
  triangle-triangle overlap (the Moller interval test), are confined to the
  pieces the intersection segment touches, and pass-through vertices on a
  triangle's own edges are dropped when its fragments are merged. Two
  coaxial ~800-triangle cylinders unioned into 32,608 triangles in 20 s;
  they now give 3,904 in about a second, and a full yapRover build with
  the native engine completes.
- New ``manifold`` boolean engine (``yapcad.boolean.manifold_engine``) that
  calls manifold3d directly, without trimesh, in float64. Install it with
  ``pip install 'yapCAD[manifold]'``. When manifold3d is installed it is now
  the default mesh boolean engine, ahead of ``native``: on a 24-case
  benchmark it returned 24 closed, correct results where the old ``native``
  engine managed 8, and it is faster than the rewritten one.
  ``YAPCAD_MESH_BOOLEAN_ENGINE=native`` restores the old default. An
  operand that is not a closed 2-manifold is refused with
  ``NotManifoldError`` when the engine is named explicitly, and falls back
  to ``native`` with a warning otherwise.

Geometry, provenance and serialisation
--------------------------------------

- Geometry JSON is now ``yapcad-geometry-json-v0.3``. It adds ``"sdf"`` as an
  ``authoritative`` representation, with the tree stored in
  ``representations.sdf`` and an optional BREP marked ``role: "derived"``.
  v0.1 and v0.2 documents still read; writers emit v0.3, so a reader that
  checks the schema ID strictly must be updated to accept it.
- Solid construction provenance (``['procedure', call]``,
  ``['boolean', operation]``, ``['sdf', tree]``) now survives geometry
  JSON as an optional ``construction`` field. The new ``yapcad.construction``
  module defines the record format and provides accessors, normalisation and
  the JSON round-trip.
- New ``yapcad.geom3d.solid_shells``
  groups a solid's surfaces into connected shells by shared edges and
  classifies each as outer or void by the sign of its volume. yapCAD mesh
  solids carry cavities as inward-wound shells in the surface list -- a hollow
  cube already reports the correct ``volumeof`` -- so the partition is derived
  on demand rather than stored in a slot that could drift out of sync with the
  geometry.

Fixes
-----

- Fix the DSL's ``sphere(radius)`` making a sphere of half the documented
  radius: it passed its argument to ``geom3d_util.sphere``, which takes a
  diameter. **Designs that call ``sphere`` in the DSL now get spheres twice
  the size they did**; halve the argument to keep the old geometry. No DSL
  source in this repository or in yapRover calls it. The Python API's
  ``geom3d_util.sphere`` is unchanged and still takes a diameter.
- Fix ``geom3d.solid()`` filing the construction record in the material
  slot. It assigned its list arguments by an emptiness test rather than by
  position, so the near-universal ``solid(surfaces, [], construction)``
  idiom stored the record as material and left construction empty; every
  ``geom3d_util`` procedure record and boolean engine record was affected,
  and serialised solids gained a spurious ``"voids": [[], []]``.
- Fix ``geometry_from_json`` writing a document's void surfaces into the
  solid's material slot (``yapcad.geometry._retessellate_brep_solid`` and
  ``service/core/tessellator`` shared the misreading). Void surfaces now
  rejoin the shell list, where mesh solids keep cavities, and the
  interchange-level ``voids`` field is computed from ``solid_shells`` on
  write.
- Fix scaling about a centre point. ``geom.scale`` and ``geom3d.scalesurface``
  composed the transform in the wrong order and scaled about ``-cent`` rather
  than ``cent``; both BREP representations already scaled about ``cent``, so a
  centred scale of a BREP solid left its mesh and BREP disagreeing. Code that
  passed a non-origin ``cent`` gets different (correct) results.
- ``geom3d.scalesolid`` now drops a solid's BREP after a non-uniform scale,
  which the BREP cannot follow. It previously stayed attached and unscaled,
  and was then written out and exported as the solid's authoritative
  geometry.
- ``Geometry.mirror`` and ``Geometry.scale`` on solids now delegate to
  ``geom3d.mirrorsolid`` / ``scalesolid``, as ``translate`` and ``rotate``
  already did. ``mirror`` previously did nothing to a mesh-only solid,
  ``scale`` produced an invalid one, and on BREP solids both discarded the
  solid's metadata.

Packaging, testing and documentation
------------------------------------

- Test the pure-Python package on Python 3.11 and 3.14, the supported-version
  endpoints, and require Python 3.11 or newer. Python 3.10 support is retired
  ahead of its October 2026 end-of-life.
- ``jsonschema`` joins the ``tests`` extra, so schema validation runs by
  default, and the conda environment files request the ``tests`` extra by
  its correct name.
- Add ``docs/SDF-DESIGN.md``, planning native signed-distance-function support:
  node-DAG representation, Lipschitz-bound tracking, dual-contouring meshing,
  CSG-tree replay into OCC for exact BREP, and mixed-authority booleans.

Version 1.1.0 (2026-07-05)
==========================

This release introduces the **metadata v1.1** system — declarative ``@meta`` /
``@ui`` decorators, an assembly/operation metadata namespace, and a mechatron
assembly-integration pipeline — alongside a formalized **pure-python core /
BREP optional-dependency** packaging split.

what's new:
-----------

  - **@meta / @ui decorators**: DSL commands can now carry declarative metadata
    via ``@meta(...)`` (design/assembly/operation namespaces) and presentation
    hints via ``@ui(...)``. Metadata is validated at check time and normalized
    into the emitted solid's metadata dict.

  - **@meta checker warnings (W310 / W311 / W312)**: The DSL checker now emits
    structured warnings for unknown or malformed ``@meta`` keys, so metadata
    typos are caught before execution.

  - **AST-level MetadataTransform**: ``@meta`` decorators are validated and
    normalised at the AST/transform layer, giving consistent metadata semantics
    across the interpreter and exporters.

  - **.meta.yaml sidecar exporter**: A new ``yapcad.io.meta_yaml`` module emits a
    ``.meta.yaml`` sidecar next to exported geometry, capturing the full
    metadata namespace (including applied assembly operations) for downstream
    tools. Exporter hooks write it automatically.

  - **Assembly resolver registry**: New ``yapcad.assembly.resolver`` provides
    ``BasicResolver`` and ``StagedResolver`` plus a resolver registry, which
    apply ``operation.kind="subtract"`` (and related) cutters against target
    parts in priority order during assembly resolution — letting a DSL declare
    process-aware cutters that are applied at assembly load rather than baked
    into each part.

  - **Mechatron assembly integration (Phases 2-4)**: New process-aware mechatron
    export path and assembly builtins bridge yapCAD assemblies to a mechatron
    ``graph.json`` canonical model. This includes off-axis (non-±Z) feature
    datums so radial joints carry correct per-hole bolt-axis orientation.

  - **Assembly load cases + bolt patterns**: ``yapcad.assembly.load_case`` adds
    ``LoadCase`` / ``LoadAttach`` / ``BoltPattern`` dataclasses, and
    ``Assembly`` gains ``add_load_case`` / ``get_load_case`` and
    ``add_bolt_pattern`` / ``get_bolt_pattern`` registries.

  - **FEA setup bridge**: New ``yapcad.package.analysis.mechatron_fea_setup``
    reads bolt patterns and load cases from a mechatron ``graph.json`` and
    produces structured FEA inputs — per-bolt world coordinates, bolt-spring
    stiffness lookup, and resolved load-attach descriptors — with a
    per-hole-datum kinematic path plus a legacy PCD/access-direction fallback.

  - **Pure-python core / BREP optional-dependency split**: yapCAD's core is now
    formally pure Python and installs from PyPI with **no compiled
    dependencies**. Solid-modeling (BREP) features remain optional and layer on
    top of OpenCASCADE via ``pythonocc-core`` (distributed through conda-forge,
    not pip). New ``[project.optional-dependencies]`` extras ``brep`` and
    ``full`` document the OCC requirement; ``yapcad.has_brep()`` reports whether
    OCC is available so downstream code can branch cleanly; and every
    ``yapcad.brep`` OCC import degrades gracefully when OCC is absent. A
    ``requires_occ`` pytest marker lets the pure-python test lane run via
    ``pytest -m "not requires_occ"``.

  - **First PyPI release of the pure-Python core**: yapCAD is now published to
    PyPI (``pip install yapcad``) as a pure-Python wheel with no compiled
    dependencies. BREP/boolean/STEP features still require ``pythonocc-core``
    from conda-forge. To make the limitation obvious, importing yapCAD without
    OCC present now emits a one-time, suppressible
    ``YapcadBrepUnavailableWarning`` describing the reduced-functionality state
    and the conda upgrade path.

Version 1.0.1 (never released; these changes shipped in 1.1.0)
===============================================================

what's new:
-----------

  - **Assembly, Kinematics, Collision, and Viewer Packages**: Added
    ``yapcad.assembly``, ``yapcad.kinematics``, ``yapcad.collision``, and
    ``yapcad.viewer`` modules for datum-driven multi-body assemblies:

    - Datum and mate abstractions for feature-based part placement
    - Constraint validation and mate solving for common mechanical joints
    - Kinematic chain support for articulated assemblies
    - Collision checks with interface volumes for designed overlaps
    - Optional VTK viewer and API server for assembly inspection

  - **Manufacturing Post-Processing**: New ``yapcad.manufacturing`` package for
    splitting swept structures into manufacturable segments:

    - Beam segmentation at specified or computed cut points
    - Interior connector specifications and mating relationships
    - STEP export support for individual segments
    - Assembly manifest generation for segmented designs

  - **Watertight Revolution Solids**: ``makeRevolutionSolid()`` now adds disc caps
    for closed profile revolutions, improving watertightness for downstream mesh
    validation and export.

  - **3D Text Support**: New ``yapcad.text3d`` module provides 3D text generation for
    labeling parts and creating extruded or engraved text:

    - ``text_solid()`` - Creates extruded 3D text as a solid
    - ``engrave_text()`` - Cuts text into a solid surface using boolean difference
    - ``text_to_polygons()`` - Converts text to 2D polygons for custom workflows
    - Supports TrueType/OpenType fonts (via freetype-py) or built-in block font
    - Auto-detects Arial font or falls back to simple printable block characters

  - **Fillet and Chamfer Operations**: New BREP edge operations in ``yapcad.brep``
    for rounding and beveling edges (requires pythonocc-core):

    - ``fillet_all_edges()`` - Apply rounded fillets to all edges of a solid
    - ``chamfer_all_edges()`` - Apply beveled chamfers to all edges of a solid
    - ``fillet_edges()`` - Apply fillets to specific selected edges
    - ``chamfer_edges()`` - Apply chamfers to specific selected edges
    - Integrates with DSL via ``fillet()``, ``chamfer()``, ``fillet_edges()``,
      and ``chamfer_edges()`` builtins

  - **Edge Selection Predicates**: New ``yapcad.brep_edge_select`` module provides
    sophisticated edge selection for selective fillet/chamfer operations:

    - ``select_vertical_edges()`` - Select edges parallel to the Z axis
    - ``select_horizontal_edges()`` - Select edges perpendicular to the Z axis
    - ``select_edges_by_direction()`` - Select edges parallel to an arbitrary direction
    - ``select_edges_by_length()`` - Select edges within a length range
    - ``select_edges_at_z()`` - Select edges at a specific Z height
    - ``select_edges_in_z_range()`` - Select edges within a Z coordinate range
    - ``select_edges_crossing_z()`` - Select vertical edges spanning a Z value
    - ``select_top_edges()`` - Select edges at the maximum Z height
    - ``select_bottom_edges()`` - Select edges at the minimum Z height
    - ``select_edges_near_point()`` - Select edges near a target point
    - ``select_edges_in_cylinder()`` - Select edges within a cylindrical region
    - ``filter_linear_edges()`` - Filter to keep only straight edges
    - ``filter_curved_edges()`` - Filter to keep only curved edges
    - Set operations: ``union_edges()``, ``intersect_edges()``, ``subtract_edges()`` for combining selections
    - ``edge_info()`` - Query detailed geometric properties of edges
    - DSL bindings currently expose direction, length, Z-position, set-operation,
      and selective fillet/chamfer helpers.

  - **Helical Extrusion**: New ``helical_extrude()`` function in ``yapcad.geom3d_util``
    creates smooth helical/twisted extrusions using high-resolution lofting. Ideal for
    helical gears, twisted columns, and spiral features. Requires pythonocc-core.

  - **Pattern Functions**: New pattern generation functions for creating arrays of geometry:

    - ``radial_pattern()`` in ``yapcad.geom_util`` - Creates circular patterns of 2D geometry
    - ``linear_pattern()`` in ``yapcad.geom_util`` - Creates linear arrays of 2D geometry
    - ``radial_pattern_solid()`` in ``yapcad.geom3d_util`` - Creates circular patterns of 3D solids
    - ``linear_pattern_solid()`` in ``yapcad.geom3d_util`` - Creates linear arrays of 3D solids
    - ``radial_pattern_surface()`` in ``yapcad.geom3d_util`` - Creates circular patterns of surfaces
    - ``linear_pattern_surface()`` in ``yapcad.geom3d_util`` - Creates linear arrays of surfaces

  - **OCC Helix Helper**: New ``make_occ_helix()`` function creates mathematically exact
    helix curves using OpenCascade's 2D parametric curve on cylindrical surface technique.
    Used internally by ``helical_extrude()`` but also available for advanced users.

  - **Bezier Curve Support**: New Bezier curve primitives in ``yapcad.geom`` and
    ``yapcad.spline`` modules for creating smooth parametric curves:

    - ``bezier()`` - Constructor for Bezier curves of any degree (linear, quadratic, cubic, etc.)
    - ``bezier_point()`` - Evaluate Bezier curve at parameter t using De Casteljau's algorithm
    - ``bezier_tangent()`` - Compute tangent vector at any point on the curve
    - ``bezier_curve()`` - Sample Bezier curve as a polyline for rendering
    - ``isbezier()`` - Type predicate for Bezier curve objects
    - Full integration with ``sample()``, ``length()``, ``center()``, and ``bbox()`` functions

  - **B-spline Curve Support**: New B-spline curve primitives for complex smooth curves
    with local control in ``yapcad.geom`` and ``yapcad.spline`` modules:

    - ``bspline()`` - Constructor for open or closed B-spline curves with configurable degree
    - ``bspline_point()`` - Evaluate B-spline curve at parameter t
    - ``bspline_tangent()`` - Compute tangent vector using numerical differentiation
    - ``bspline_curve()`` - Sample B-spline curve as a polyline
    - ``isbspline()`` - Type predicate for B-spline curve objects
    - Support for both open and closed curves via ``closed=True`` parameter
    - Full integration with yapCAD geometry functions (``sample()``, ``length()``, ``center()``, ``bbox()``)

  - **Curve-based 3D Operations**: New functions in ``yapcad.geom3d_util`` for creating
    solids from curved paths:

    - ``loft()`` - Create smooth solids by transitioning between 2D profiles at different Z heights
    - ``sweep_along_path()`` - Extrude 2D profiles along Bezier or B-spline path curves
    - Support for organic shapes, tapered forms, and complex swept geometries

  - **Expanded DSL Builtins**: Added DSL access for hollow primitives
    (``tube()``, ``conic_tube()``, ``spherical_shell()``), ``dodecahedron()``,
    path utilities, text placement helpers, advanced gear helpers, and pattern
    generation.

  - **Documentation Improvements**: All new functions include comprehensive docstrings
    with parameter descriptions, return values, usage examples, and notes following
    Sphinx documentation standards.

chore:
------

  - Removed internal project files from ``projects/`` so private robot, hanger,
    and jig designs can live in separate internal repositories.

Version 1.0.0rc2 (2026-01-04)
=============================

- DSL static verifiability: ``while`` loops are removed and execution has
  resource limits, so a program cannot run without bound.
- Emacs major mode for the DSL, a fastener catalog system, and a DSL
  fastener example.
- Documentation and roadmap updates for the 1.0.0 release candidates.

Version 1.0.0rc1 (2025-12-30)
=============================

**First release candidate for yapCAD 1.0**

This release represents the culmination of the yapCAD 1.0 roadmap, delivering a
requirements-driven, provenance-aware design platform with full solid modeling,
packaging, and validation capabilities.

what's new:
-----------

  - **Validation Schema Implementation**: Complete validation test schema with
    ``yapcad.package.analysis.schema`` module providing ``ValidationKind``,
    ``ResultStatus``, ``ComparisonOp`` enums and ``validate_plan()``,
    ``validate_result()`` APIs. Supports geometric, measurement, structural,
    thermal, CFD, multiphysics, and assembly test kinds.

  - **Package Signing**: Provisional package signing with GPG and SSH keys.
    ``sign_package()``, ``verify_package()``, ``list_signatures()`` APIs for
    cryptographic verification of package integrity and authorship.

  - **DSL Enhancements**: Conditional expressions (``value if cond else other``),
    nested list comprehensions (``[f(x,y) for x in xs for y in ys if cond]``),
    and functional combinators (``union_all``, ``difference_all``, ``sum``,
    ``product``, ``min_of``, ``max_of``, ``any_true``, ``all_true``).

  - **Documentation Reorganization**: Historical planning documents moved to
    ``docs/historical/``. All specification documents updated to v1.0 status.

  - **Test Suite**: 584 passing tests covering geometry, DSL, import/export,
    packaging, validation schema, and fasteners.

specifications finalized:
-------------------------

  - ``ycpkg-spec-v1.0`` - Package format specification
  - ``metadata-namespace-v1.0`` - Metadata dictionary specification
  - ``validation-schema-v1.0`` - Validation test schema
  - ``ycpkg-signing-v1.0`` - Package signing specification
  - ``yapcad-geometry-json-v1.0`` - Geometry JSON schema
  - Material schema specification v1.0

deferred to 1.1:
----------------

  - Multi-signature approval workflows
  - Delegation and authority chains
  - Revocation lists
  - Full audit trails
  - Solver adapter framework (generalized FEA/simulation integration)

Version 0.6.3 (2025-12-17)
==========================

what's new:
-----------

  - **DSL 2D Geometry Features**: New curve types (ellipse, catmull-rom splines,
    NURBS), 2D regions (polygon, disk), and 2D boolean operations (union2d,
    difference2d, intersection2d) with proper hole accumulation for chained
    operations.
  - **DXF Export**: 2D geometry can now be exported to DXF format via
    ``--output foo.dxf`` for visualization in CAD programs. Supports lines,
    arcs, ellipses, splines, polylines, and regions.
  - **Curve Operations**: New ``sample_curve()`` and ``curve_length()`` builtins
    for sampling points along curves and measuring arc length.
  - **OCC BREP Kernel Complete**: Full bidirectional conversion between native
    yapCAD BREP representation and OpenCascade shapes. Box, cylinder, sphere,
    and cone primitives achieve 100% volume fidelity in round-trip testing.
    All seven development phases documented in ``docs/BREP_integration_strategy.md``
    are now complete.
  - **Adaptive Sweep Operations**: ``sweep_adaptive()`` and ``sweep_adaptive_hollow()``
    with tangent-tracking profile orientation and ruled lofting.
  - **Materials & Fasteners**: Added material property schema (``docs/material_schema_spec.md``)
    supporting density, color, finish, and physical properties. Complete metric and
    unified fastener catalogs with proper thread geometry.
  - **Documentation Accuracy**: Updated ``docs/yapCADone.rst`` roadmap to clearly
    separate implemented features from planned work.
  - **Viewer Clipping Planes**: Added X/Y/Z clipping plane toggles to the package
    viewer for inspecting interior geometry. Press X, Y, or Z to cycle through
    off/+/- states; press C to clear all planes. Essential for examining screw/nut
    fit and internal features.

bug fixes:
----------

  - Fixed ``difference2d`` to properly accumulate holes in chained operations.
    Previously, subtracting multiple holes would lose earlier holes due to
    structure flattening.
  - Fixed ``isinsideXY`` corner case workaround in 2D booleans where test rays
    passing through polygon vertices caused incorrect inside/outside detection.
  - Removed unsupported SOLID entities from DXF output by using ``setup=False``
    in ezdxf initialization, improving FreeCAD compatibility.

Known problems
--------------

- Analytic STEP export preserves BREP data but complex surfaces may tessellate.

Version 0.6.1 (2025-10-30)
==========================

what's new:
-----------

  - **Documentation Polish**: Converted packaging/DSL/BREP specifications to
    reStructuredText, added the roadmap snapshot to the Sphinx toctree, and
    rebased all release references on the refreshed docs set. Sphinx builds
    now complete without missing toctree warnings.
  - **Analysis Planning**: Documented validation plan schema, expanded the
    metadata namespace for analysis annotations, and introduced the
    ``ycpkg_analyze`` CLI plus plan-loading APIs. The bundled CalculiX backend
    now generates an axisymmetric plate approximation, invokes ``ccx`` when
    available, and records maximum axial deflection for acceptance checks.
  - **Viewer Performance**: Layered triangle meshes are cached as pyglet vertex
    lists, eliminating per-frame immediate-mode uploads and noticeably improving
    pan/rotate responsiveness in the four-view window.

Known problems
--------------

- Analytic STEP export and full BREP kernel support remain in planning (`docs/yapBREP.rst`).
- DSL compiler/validator tooling tracked in `docs/dsl_spec.rst` remains under development.

Version 0.6.0 (2025-10-28)
==========================

what's new:
-----------

  - **Package Format**: `.ycpkg` manifests now capture canonical entities, reusable
    instances, metadata layers, and JSON geometry version `0.6.0` for reliable
    round-tripping between authoring, export, and viewer tooling.
  - **Involute Gear Toolkit**: Vendored the MIT-licensed `figgear` generator
    (`yapcad.contrib.figgear`) so gear examples and tests no longer require external
    clones. Canonical gear packages embed both 2D profiles and extruded solids with
    provenance metadata.
  - **Viewer Enhancements**: Added four-view layout updates (grids, lighting,
    layer toggles, help overlay) and fixed revolution solid cap normals along with
    trackpad gesture support.
  - **DXF / Spline Upgrades**: Native exporter writes analytic `CIRCLE`/`ARC`
    elements, spline loops remain watertight, and extrusion tests cover spline-based
    perimeters and holes.
  - **Documentation Refresh**: Roadmap snapshot (`docs/yapCADone.rst`), BREP plan
    (`docs/yapBREP.rst`), and README/Sphinx front matter updated for the 0.6 series.

Known problems
--------------

- Analytic STEP export and full BREP kernel support remain in planning (`docs/yapBREP.rst`).
- DSL compiler/validator tooling tracked in `docs/dsl_spec.rst` remains under development.

Version 0.5.1 (2025-10-14)
==========================

what's new:
-----------

  - **3D Boolean Operations Fixes**: Complete overhaul of solid boolean operations
    with robust normal orientation and interior triangle filtering.

    - Fixed sphere union normal orientation issues by filtering interior overlap triangles
    - Added quality-based filtering for degenerate sliver triangles (aspect ratio checks)
    - Implemented containment-based filtering to remove artifacts in overlap regions
    - All primitive tests now pass with correct watertight geometry

  - **2D Boolean Operations Fixes**: Resolved crash when performing boolean operations
    on ``Circle`` and other single-geometry primitives.

    - Fixed geometry wrapping in ``Boolean._prepare_geom()`` to handle unwrapped arc format
    - Added comprehensive regression tests for 2D boolean operations

  - **Primitive Improvements**: Enhanced reliability of 3D geometric primitives.

    - Fixed ``conic()`` primitive to generate proper watertight solids
    - Fixed ``tube()`` primitive normal orientation and end cap connectivity
    - All 9 core primitives (box, sphere, cylinder, cone, tube, etc.) validated as watertight

  - **Modular Boolean Engine Architecture**: Separated boolean operations into
    ``yapcad.boolean.native`` module for better maintainability.

    - Support for multiple boolean engine backends (native, trimesh:manifold, trimesh:blender)
    - Engine selection via ``solid_boolean(..., engine='native')`` parameter
    - Environment variable support (``YAPCAD_BOOLEAN_ENGINE``, ``YAPCAD_TRIMESH_BACKEND``)

  - **Test Suite Improvements**: Enhanced test coverage and reliability.

    - 106 tests passing (up from 99 in v0.5.0)
    - Added boolean regression test suite
    - Improved solid topology tests with better error reporting

Known problems
--------------

- Incomplete documentation for some advanced 3D features.
- STEP export currently supports tessellated geometry; analytical BREP support planned for 1.0.

Version 0.5.0 (2024-09-30)
==========================

what's new:
-----------

  - Adds shared geometry utilities and metadata helpers.
  - Introduces STL export (`yapcad.io.stl`) plus tests.
  - Provides `examples/rocket_demo.py` showing a full 3D workflow.
  - Updates documentation with 3D-focused imagery and instructions.

Known problems
--------------

- Incomplete documentation, though this is improving.

Version 0.4.0 (Development)
============================

what's new:
-----------

- **Testing Infrastructure Overhaul**: Completely redesigned test execution system
  to properly support both automated and interactive visual tests.

  - Added comprehensive pytest markers: ``@pytest.mark.visual`` for interactive tests
  - Created ``run_visual_tests.py`` and ``run_visual_tests_venv.sh`` for isolated
    visual test execution using subprocess isolation
  - Enhanced test discovery using AST parsing to automatically find decorated visual tests
  - Fixed visual test termination issues that were causing pytest to exit prematurely
  - Updated all test documentation with clear separation between non-visual and visual testing

- **3D Geometry Enhancements**: Merged advanced 3D surface representation and
  geometry system improvements from development branch, including enhanced
  ``Geometry`` class architecture and improved computational geometry operations.

Known problems
--------------

- Incomplete documentation, especially outside the ``yapcad.geom`` module.
- Occasional problems with complex boolean operations.
- Incomplete functionality around 3D modeling.

Version 0.3.1
=============

what's new:
-----------

- Added Read the Docs configuration and ``docs/requirements.txt`` so hosted
  builds use a consistent environment.
- Updated README instructions for building documentation and running tests.
- Follow-up to 0.3.0 (no functional code changes).

Known problems
--------------

- Incomplete documentation, especially outside the ``yapcad.geom`` module.
- Occasional problems with complex boolean operations.
- Incomplete functionality around 3D modeling.

Version 0.3.0
=============

what's new:
-----------

- Require Python 3.10+ and align dependency metadata with current
  interpreter and library versions.
- Pin pyglet to 1.x rendering backend and add fallback
  guards to every OpenGL-enabled example so they degrade gracefully on
  systems without a working pyglet/Cocoa stack.
- Sphinx documentation now builds even when optional themes are
  missing, and `sphinx-apidoc` no longer depends on ``pkg_resources``.

Known problems
--------------

- Incomplete documentation, especially outside the ``yapcad.geom`` module.
- Occasional problems with complex boolean operations.
- Incomplete functionality around 3D modeling.

Version 0.2.0
=============

what's new:
-----------

- First announced version of **yapCAD**. Yay!

- Added new ``boxcut`` example, showing a fully worked (if simple)
  parametric design system.

- Additional documentation updates and minor bugfixes.

Known problems
--------------

- Our `yapCAD readthedocs`_ documentation is missing the expanded
  documentation from submodules, which is a problem since much of
  **yapCAD**'s documentation is in the form of docstrings in the
  source.  I'm working on getting this sorted out.  In the mean time,
  you may want to build a local copy of the documentation as described
  in the main ``README`` file.   Or, checkout and read the source.

- Incomplete documentation, especially outside the ``yapcad.geom`` module.

- Occasional problems with complex boolean operations.  A bug in the
  ``intersectXY`` method of the ``Boolean`` class.

- Incomplete functionality around 3D modeling

- Inconsistent inclusion of licensing boilerplate, other minor
  formatting issues.

Version 0.1.5
=============

what's new:
-----------

- Pre-release, heading towards V0.2.x

- Restructuring for package release

- Lots more documentation (still incomplete)

- Fixes to package configuration

Known problems
--------------

- Incomplete documentation, especially outside the ``yapcad.geom`` module.

- Occasional problems with complex boolean operations

- Incomplete functionality around 3D modeling

- Inconsistent inclusion of licensing boilerplate
  

.. _yapCAD readthedocs: https://yapcad.readthedocs.io/en/latest/index.html
