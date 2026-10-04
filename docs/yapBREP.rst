yapCAD BREP Implementation
==========================


This document describes yapCAD's BREP (Boundary Representation) implementation,
which provides native solid modeling through integration with OpenCascade (OCC).

Overview
--------

yapCAD solids are triangle meshes that can carry a more precise
definition: an OpenCASCADE BREP, described here, or a signed distance field
(see :doc:`sdf_guide`). :doc:`representations` compares the three.

With ``pythonocc-core`` installed, primitives carry an exact BREP, booleans
between BREP solids use OCC, and STEP import and analytic STEP export are
available. Without it, solids are meshes or fields.

Implementation Status
---------------------

**Complete:**

* Full OCC BREP wrapper classes (``BrepSolid``, ``BrepFace``, ``BrepEdge``, ``BrepVertex``)
* BREP-aware primitives in ``yapcad.geom3d_util``: ``prism()``, ``conic()``
  (cylinders and cones), ``sphere()``, ``tube()``
* Solid operations: ``extrude()``, ``makeRevolutionSolid()``, ``makeLoftSolid()``
* Adaptive sweep: ``sweep_adaptive()``, hollow with ``inner_profiles=``
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

Without OCC there is no STEP import or analytic STEP export, booleans use
the mesh or field engines (see :doc:`representations`), and
``sweep_adaptive()`` raises. See :doc:`installation`.

Architecture
------------

BREP Wrapper Classes
~~~~~~~~~~~~~~~~~~~~

Located in ``yapcad/brep.py``:

* ``BrepSolid`` - Wraps OCC TopoDS_Solid with lazy tessellation
* ``BrepFace`` - Surface face with trim curves
* ``BrepEdge`` - Edge with curve geometry
* ``BrepVertex`` - Vertex with tolerance

Converting between yapCAD solids and OCC shapes::

    from yapcad.brep import BrepSolid, attach_brep_to_solid, brep_from_solid
    from yapcad.geom3d import solid
    from yapcad.geom3d_util import prism

    # yapCAD solid -> OCC shape
    occ_shape = brep_from_solid(prism(10, 20, 30)).shape

    # OCC shape -> yapCAD solid: tessellate it, and keep the BREP attached
    brep = BrepSolid(occ_shape)
    part = solid([brep.tessellate()])
    attach_brep_to_solid(part, brep)

BREP Attachment
~~~~~~~~~~~~~~~

A solid keeps its BREP, serialised, in its metadata under ``'brep'``.
Read it with ``brep_from_solid``, which returns ``None`` for a solid
without one::

    from yapcad.brep import brep_from_solid, has_brep_data
    from yapcad.geom3d_util import prism

    part = prism(10, 20, 30)
    has_brep_data(part)            # True when OCC is installed
    brep = brep_from_solid(part)   # a BrepSolid, or None

Boolean Operations
~~~~~~~~~~~~~~~~~~

``yapcad.geom3d.solid_boolean`` chooses an engine per call. An explicit
``engine=`` argument, or the ``YAPCAD_BOOLEAN_ENGINE`` environment variable,
forces one of ``native``, ``trimesh``, ``occ``, ``manifold`` or ``sdf``.
Otherwise:

1. If both operands were authored as signed distance fields → combine the
   fields (``yapcad.sdf.booleans``); the result stays SDF-authoritative.
2. If both operands have BREP data and OCC is available → use OCC.
3. Otherwise → a mesh engine: ``YAPCAD_MESH_BOOLEAN_ENGINE`` if set, else
   ``manifold`` when manifold3d is installed, else ``native``.

The ``manifold`` engine calls manifold3d directly (``pip install
'yapCAD[manifold]'``), in float64. It is robust on the cases mesh booleans
find hard -- shared and coincident faces, general-position overlaps -- and
refuses an operand that is not a closed 2-manifold rather than guess; in
automatic selection such an operand falls back to ``native`` with a
warning. The ``native`` engine (``yapcad.boolean.csg``) needs only numpy:
it splits triangles by the planes of the triangles they cross, classifies
fragments by generalized winding number, and repairs the result into a
closed, conforming mesh, including for shared and coincident faces. It is
slower than ``manifold`` -- seconds rather than milliseconds at ten
thousand triangles -- and is the fallback when manifold3d is absent.

STEP Import/Export
------------------

Import
~~~~~~

``import_step`` returns one ``Geometry`` per solid in the file, each
wrapping a ``BrepSolid``::

    from yapcad.brep import attach_brep_to_solid
    from yapcad.geom3d import solid
    from yapcad.io.step_importer import import_step

    parts = []
    for geometry in import_step("model.step"):
        part = solid([geometry.surface()])        # its tessellation
        attach_brep_to_solid(part, geometry.geom)  # and its exact BREP
        parts.append(part)

Export
~~~~~~

Two modes available:

**Faceted**, from any solid's mesh::

    from yapcad.io.step import write_step
    write_step(part, "output.step")

**Analytic**, from the solid's BREP, falling back to faceted when it has
none (the return value says which you got)::

    from yapcad.io.step import write_step_analytic
    analytic = write_step_analytic(part, "output.step")

From the DSL command line, ``YAPCAD_STEP_FORMAT=analytic`` selects the
analytic writer for ``yapcad.dsl run ... -o part.step``.

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
    from yapcad.geom3d_util import prism

    # Create a box solid
    my_box = prism(20, 20, 10)
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
    from yapcad.geom3d_util import prism

    # Create a box
    my_box = prism(20, 20, 10)
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

* ``brep_from_solid(solid)`` - The solid's ``BrepSolid``, or ``None``
* ``attach_brep_to_solid(solid, brep)`` - Attach a ``BrepSolid`` to a solid
* ``has_brep_data(solid)`` / ``occ_available()`` - What is present
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

``yapcad.io.step`` and ``yapcad.io.step_importer``:

* ``import_step(path)`` - One ``Geometry`` per solid in a STEP file
* ``write_step(solid, path)`` - Faceted STEP
* ``write_step_analytic(solid, path)`` - Analytic STEP from the BREP

``yapcad.geom3d_util`` (with a BREP when OCC is installed):

* ``prism(length, width, height, center=...)`` - Box, centred on ``center``
* ``conic(base_r, top_r, height, center=...)`` - Cylinder or cone, standing on ``center``
* ``sphere(diameter, center=...)`` - Sphere (note: a diameter)
* ``tube(outer_diameter, wall_thickness, length)`` - Hollow cylinder
* ``extrude(surface, distance, direction=...)`` - Linear extrusion of a surface
* ``makeRevolutionSolid(contour, zStart, zEnd, steps, ...)`` - Solid of revolution
* ``makeLoftSolid(lower_loop, upper_loop)`` - Loft between two loops
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

BREP tests are in ``tests/test_brep.py``, ``tests/test_brep_edge_select.py``,
``tests/test_step_importer.py`` and ``tests/test_io_step.py``. Run with::

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
* Promoting BREP solids to fields, for booleans between BREP and SDF parts
