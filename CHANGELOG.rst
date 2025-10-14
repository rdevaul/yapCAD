=========
Changelog
=========

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
