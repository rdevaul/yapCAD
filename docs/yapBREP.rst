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

API Reference
-------------

Key Functions
~~~~~~~~~~~~~

``yapcad.brep``:

* ``import_step(path)`` - Import STEP file
* ``write_step_analytic(solid, path)`` - Export analytic STEP
* ``brep_from_solid(solid)`` - Extract BREP wrapper
* ``solid_from_brep(shape)`` - Create solid from OCC shape

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
