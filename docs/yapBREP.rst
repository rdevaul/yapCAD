yapCAD BREP Implementation
==========================


This document describes yapCAD's BREP (Boundary Representation) implementation,
which provides native solid modeling through integration with OpenCascade (OCC).

Overview
--------

yapCAD supports two modes of solid geometry:

1. **Tessellated Mode** (default, no dependencies) - Solids are represented as
   triangle meshes. Suitable for visualization and STL export.

2. **BREP Mode** (requires pythonocc-core) - Solids carry native OCC BREP data
   with exact geometric definitions. Required for STEP import/export and
   high-fidelity boolean operations.

Implementation Status
---------------------

**Complete:**

* Full OCC BREP wrapper classes (``BrepSolid``, ``BrepFace``, ``BrepEdge``, ``BrepVertex``)
* BREP-aware primitives: ``box()``, ``sphere()``, ``cylinder()``, ``cone()``, ``prism()``
* Solid operations: ``extrude()``, ``tube()``, ``makeRevolutionSolid()``, ``makeLoftSolid()``
* Adaptive sweep: ``sweep_adaptive()``, ``sweep_adaptive_hollow()``
* Boolean operations via OCC when BREP data present
* STEP import with topology preservation
* STEP export (tessellated and analytic modes)
* STL import/export
* Transform propagation to BREP data

**Curves (2D/3D):**

* Line, arc, circle, ellipse (full and partial arcs)
* Catmull-Rom splines (open and closed)
* NURBS curves
* Parabola, hyperbola primitives

**Surfaces:**

* Planar, cylindrical, conical, spherical surfaces via OCC
* Tessellation on demand with quality controls
* Surface evaluation and normal computation

Installation
------------

BREP functionality requires pythonocc-core, installed via conda::

    conda env create -f environment.yml
    conda activate yapcad-brep

Without OCC, yapCAD operates in tessellation-only mode with reduced functionality:

* No STEP import/export
* Boolean operations use mesh-based algorithms (lower fidelity)
* Sweep operations produce tessellated results only

Architecture
------------

BREP Wrapper Classes
~~~~~~~~~~~~~~~~~~~~

Located in ``yapcad/brep.py``:

* ``BrepSolid`` - Wraps OCC TopoDS_Solid with lazy tessellation
* ``BrepFace`` - Surface face with trim curves
* ``BrepEdge`` - Edge with curve geometry
* ``BrepVertex`` - Vertex with tolerance

These classes provide bidirectional conversion::

    # yapCAD solid -> OCC shape
    from yapcad.brep import brep_from_solid
    occ_shape = brep_from_solid(solid).shape

    # OCC shape -> yapCAD solid
    from yapcad.brep import solid_from_brep
    solid = solid_from_brep(occ_shape)

BREP Attachment
~~~~~~~~~~~~~~~

Solids store BREP data in metadata under ``'brep'`` key::

    solid = box(10, 20, 30)
    brep = solid[5].get('brep')  # BrepSolid instance if OCC available

Boolean Operations
~~~~~~~~~~~~~~~~~~

Environment variables control engine selection:

* ``YAPCAD_BOOLEAN_ENGINE`` - Force engine: ``native``, ``trimesh``, ``occ``
* ``YAPCAD_MESH_BOOLEAN_ENGINE`` - Mesh fallback when OCC unavailable

Auto-selection logic:

1. If both operands have BREP data → use OCC
2. Otherwise → use mesh-based engine

STEP Import/Export
------------------

Import
~~~~~~

::

    from yapcad.brep import import_step

    # Returns yapCAD Geometry with attached BREP
    geometry = import_step("model.step")

    # Access the underlying solid
    solid = geometry.geom

Export
~~~~~~

Two modes available:

**Tessellated (default)**::

    from yapcad.geom3d import write_solid_step
    write_solid_step(solid, "output.step")

**Analytic** (preserves exact geometry)::

    # Via environment variable
    import os
    os.environ['YAPCAD_STEP_FORMAT'] = 'analytic'
    write_solid_step(solid, "output.step")

    # Or via write_step_analytic directly
    from yapcad.brep import write_step_analytic
    write_step_analytic(solid, "output.step")

Adaptive Sweep Operations
-------------------------

The ``sweep_adaptive()`` function creates solids by sweeping a profile along
a path with tangent-tracking orientation::

    from yapcad.geom3d_util import sweep_adaptive

    solid = sweep_adaptive(
        profile,           # region2d for cross-section
        spine,             # path3d with line/arc segments
        angle_threshold_deg=5.0,  # threshold for new section
        frame_mode='minimal_twist'  # or 'frenet', 'custom'
    )

For hollow profiles (pipes)::

    solid = sweep_adaptive(
        outer_profile,
        spine,
        inner_profiles=[inner_profile],  # one or more voids
        angle_threshold_deg=5.0
    )

Edge Treatment: Fillets and Chamfers
-------------------------------------

Apply rounded (fillet) or beveled (chamfer) edges to BREP solids::

    from yapcad.brep import (
        fillet_all_edges, chamfer_all_edges,
        fillet_edges, chamfer_edges,
        brep_from_solid
    )
    from yapcad.geom3d_util import box

    # Create a box solid
    my_box = box(20, 20, 10)
    brep_solid = brep_from_solid(my_box)

    # Apply fillet (rounded edges) to all edges
    filleted = fillet_all_edges(brep_solid, radius=2.0)

    # Apply chamfer (beveled edges) to all edges
    # Creates symmetric 45° chamfers
    chamfered = chamfer_all_edges(brep_solid, distance=1.5)

    # Apply fillet to specific edges only
    selected_edges = [edge1, edge2]  # BrepEdge objects
    partial_fillet = fillet_edges(brep_solid, selected_edges, radius=1.0)

    # Apply chamfer to specific edges
    partial_chamfer = chamfer_edges(brep_solid, selected_edges, distance=1.0)

**Notes:**

* Requires pythonocc-core (BREP environment)
* Returns a new ``BrepSolid`` instance with modified edges
* Edges that cannot be filleted/chamfered (too small, geometric constraints) are skipped
* For selective edge operations, use edge selection helpers (see below)
* See ``examples/fillet_chamfer_demo.dsl`` for DSL usage examples

Edge Selection Functions
-------------------------

The ``yapcad.brep_edge_select`` module provides functions to select edges from BREP
solids based on geometric criteria. This is particularly useful for applying selective
fillets or chamfers to specific edges rather than all edges of a solid.

Basic Selection Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~

**Select by Direction:**

::

    from yapcad.brep_edge_select import (
        select_vertical_edges,
        select_horizontal_edges,
        select_edges_by_direction
    )

    # Select edges parallel to Z axis (vertical edges)
    vertical = select_vertical_edges(brep_solid, tolerance_deg=1.0)

    # Select edges perpendicular to Z axis (horizontal edges)
    horizontal = select_horizontal_edges(brep_solid, tolerance_deg=1.0)

    # Select edges parallel to a custom direction
    edges_45deg = select_edges_by_direction(brep_solid,
                                            direction=(1, 1, 0),
                                            tolerance_deg=1.0)

**Select by Z Position:**

::

    from yapcad.brep_edge_select import (
        select_top_edges,
        select_bottom_edges,
        select_edges_at_z,
        select_edges_in_z_range,
        select_edges_crossing_z
    )

    # Select edges at the top of the solid (maximum Z)
    top_edges = select_top_edges(brep_solid, tolerance=0.001)

    # Select edges at the bottom of the solid (minimum Z)
    bottom_edges = select_bottom_edges(brep_solid, tolerance=0.001)

    # Select edges at a specific Z height
    mid_edges = select_edges_at_z(brep_solid, z_value=5.0, tolerance=0.001)

    # Select edges within a Z range
    range_edges = select_edges_in_z_range(brep_solid,
                                          z_min=2.0, z_max=8.0,
                                          tolerance=0.001)

    # Select edges that cross a specific Z height (vertical edges spanning Z)
    crossing = select_edges_crossing_z(brep_solid, z_value=5.0, tolerance=0.001)

**Select by Length:**

::

    from yapcad.brep_edge_select import select_edges_by_length

    # Select edges within a length range
    long_edges = select_edges_by_length(brep_solid, min_length=10.0)
    short_edges = select_edges_by_length(brep_solid, max_length=5.0)
    mid_edges = select_edges_by_length(brep_solid,
                                       min_length=5.0,
                                       max_length=15.0)

**Select by Position:**

::

    from yapcad.brep_edge_select import (
        select_edges_near_point,
        select_edges_in_cylinder
    )

    # Select edges near a point (based on edge midpoint)
    near_origin = select_edges_near_point(brep_solid,
                                          target_point=(0, 0, 5),
                                          max_distance=2.0)

    # Select edges within a cylindrical region (useful for holes/bosses)
    around_hole = select_edges_in_cylinder(brep_solid,
                                           center=(10, 10, 0),
                                           radius=5.0,
                                           axis=(0, 0, 1))

Filtering Functions
~~~~~~~~~~~~~~~~~~~

::

    from yapcad.brep_edge_select import (
        filter_curved_edges,
        filter_linear_edges,
        get_all_edges
    )

    # Get all edges from a solid
    all_edges = get_all_edges(brep_solid)

    # Filter to only curved (non-linear) edges
    curved = filter_curved_edges(all_edges)

    # Filter to only linear (straight) edges
    linear = filter_linear_edges(all_edges)

Set Operations on Edge Lists
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Combine, intersect, or subtract edge selections::

    from yapcad.brep_edge_select import (
        union_edges,
        intersect_edges,
        subtract_edges
    )

    # Combine multiple selections (removes duplicates)
    combined = union_edges(vertical_edges, top_edges, bottom_edges)

    # Find edges common to all selections
    common = intersect_edges(horizontal_edges, long_edges, top_edges)

    # Remove edges from a selection
    all_but_bottom = subtract_edges(all_edges, bottom_edges)

Edge Information
~~~~~~~~~~~~~~~~

Get detailed information about an edge::

    from yapcad.brep_edge_select import edge_info

    info = edge_info(some_edge)
    # Returns dict with:
    #   - length: Edge length
    #   - endpoints: ((x1, y1, z1), (x2, y2, z2))
    #   - midpoint: (x, y, z)
    #   - direction: (dx, dy, dz) for linear edges, None for curved
    #   - is_linear: True if edge is straight
    #   - is_vertical: True if parallel to Z axis
    #   - is_horizontal: True if perpendicular to Z axis

Complete Example: Selective Filleting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here's a complete example showing how to apply fillets only to specific edges::

    from yapcad.brep import brep_from_solid, fillet_edges
    from yapcad.brep_edge_select import (
        select_vertical_edges,
        select_top_edges,
        intersect_edges
    )
    from yapcad.geom3d_util import box

    # Create a box
    my_box = box(20, 20, 10)
    brep = brep_from_solid(my_box)

    # Select only the vertical edges at the top of the box
    vertical = select_vertical_edges(brep)
    top = select_top_edges(brep)
    top_vertical = intersect_edges(vertical, top)

    # Apply fillet only to these edges
    filleted_box = fillet_edges(brep, top_vertical, radius=1.0)

**Notes:**

* All selection functions require pythonocc-core
* Selection functions return lists of ``BrepEdge`` objects
* Edge selections can be combined using set operations
* Tolerance parameters control how precisely edges must match criteria
* Edge selections avoid duplicates (each unique edge appears once)

Helical Extrusion
-----------------

Create smooth helical extrusions (twisted features) for gears, columns, and spiral geometry::

    from yapcad.geom3d_util import helical_extrude
    from yapcad.geom import line, point

    # Define a square profile centered at origin
    square = [
        line(point(-5, -5), point(5, -5)),
        line(point(5, -5), point(5, 5)),
        line(point(5, 5), point(-5, 5)),
        line(point(-5, 5), point(-5, -5))
    ]

    # Create twisted column: 20mm tall, 90° total twist
    twisted_column = helical_extrude(
        square,
        height=20,
        twist_angle_deg=90,
        segments=64  # Higher = smoother (recommended: 64+)
    )

**Parameters:**

* ``profile`` - region2d (closed loop of 2D curves) centered near origin
* ``height`` - Extrusion height along Z-axis
* ``twist_angle_deg`` - Total rotation angle (positive = counterclockwise from +Z)
* ``segments`` - Number of lofting sections (default 64, use 64+ for smooth surfaces)

**Notes:**

* Requires pythonocc-core (uses lofting for smooth surfaces)
* Zero twist produces simple extrusion
* Ideal for helical gears (use with involute gear profiles)

API Reference
-------------

Key Functions
~~~~~~~~~~~~~

``yapcad.brep``:

* ``import_step(path)`` - Import STEP file
* ``write_step_analytic(solid, path)`` - Export analytic STEP
* ``brep_from_solid(solid)`` - Extract BREP wrapper
* ``solid_from_brep(shape)`` - Create solid from OCC shape
* ``fillet_all_edges(brep_solid, radius)`` - Round all edges
* ``chamfer_all_edges(brep_solid, distance)`` - Bevel all edges
* ``fillet_edges(brep_solid, edges, radius)`` - Round specific edges
* ``chamfer_edges(brep_solid, edges, distance)`` - Bevel specific edges

``yapcad.brep_edge_select``:

* ``get_all_edges(brep_solid)`` - Get all unique edges as BrepEdge objects
* ``select_vertical_edges(brep_solid, tolerance_deg)`` - Select edges parallel to Z
* ``select_horizontal_edges(brep_solid, tolerance_deg)`` - Select edges perpendicular to Z
* ``select_edges_by_direction(brep_solid, direction, tolerance_deg)`` - Select by direction
* ``select_edges_by_length(brep_solid, min_length, max_length)`` - Select by length range
* ``select_edges_at_z(brep_solid, z_value, tolerance)`` - Select edges at Z height
* ``select_edges_in_z_range(brep_solid, z_min, z_max, tolerance)`` - Select edges in Z range
* ``select_edges_crossing_z(brep_solid, z_value, tolerance)`` - Select edges spanning Z
* ``select_top_edges(brep_solid, tolerance)`` - Select edges at maximum Z
* ``select_bottom_edges(brep_solid, tolerance)`` - Select edges at minimum Z
* ``select_edges_near_point(brep_solid, point, max_distance)`` - Select by proximity
* ``select_edges_in_cylinder(brep_solid, center, radius, axis)`` - Select in cylindrical region
* ``filter_curved_edges(edges)`` - Keep only curved edges
* ``filter_linear_edges(edges)`` - Keep only linear edges
* ``edge_info(edge)`` - Get detailed edge properties
* ``union_edges(*edge_lists)`` - Combine edge lists (remove duplicates)
* ``intersect_edges(*edge_lists)`` - Find common edges
* ``subtract_edges(base_edges, edges_to_remove)`` - Remove edges from list

``yapcad.geom3d``:

* ``box(l, w, h)`` - Rectangular prism with BREP
* ``sphere(center, radius)`` - Sphere with BREP
* ``cylinder(center, axis, radius, height)`` - Cylinder with BREP
* ``cone(center, axis, radius1, radius2, height)`` - Cone/frustum with BREP

``yapcad.geom3d_util``:

* ``extrude(profile, direction, height)`` - Linear extrusion
* ``makeRevolutionSolid(profile, axis, angle)`` - Revolution
* ``makeLoftSolid(profiles)`` - Loft through profiles
* ``sweep_adaptive(profile, spine, ...)`` - Adaptive sweep
* ``helical_extrude(profile, height, twist_angle_deg, ...)`` - Helical/twisted extrusion
* ``radial_pattern_solid(solid, count, ...)`` - Circular array of solids
* ``linear_pattern_solid(solid, count, spacing)`` - Linear array of solids
* ``radial_pattern_surface(surf, count, ...)`` - Circular array of surfaces
* ``linear_pattern_surface(surf, count, spacing)`` - Linear array of surfaces

``yapcad.geom_util``:

* ``radial_pattern(geometry, count, ...)`` - Circular array of 2D geometry
* ``linear_pattern(geometry, count, spacing)`` - Linear array of 2D geometry

``yapcad.text3d``:

* ``text_solid(text, height, depth, ...)`` - Create 3D extruded text
* ``engrave_text(target, text, position, normal, ...)`` - Engrave text into solid
* ``text_width(text, height, spacing, font)`` - Calculate text width for positioning

Testing
-------

BREP tests are in ``tests/test_brep.py`` and ``tests/test_step_import.py``.
Run with::

    PYTHONPATH=./src pytest tests/test_brep.py -v

Limitations
-----------

* Complex surfaces (NURBS, offset) may tessellate during operations
* Boolean operations on complex geometry may fail; check results
* Performance depends on OCC installation quality

Future Work
-----------

* Improved NURBS surface support
* Direct BREP editing operations
* Better error reporting for failed operations
* Integration with DSL for BREP-specific operations
