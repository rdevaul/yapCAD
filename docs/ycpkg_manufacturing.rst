Package Manufacturing Exports
=============================

Assembly-aware packages can export canonical component-local fabrication
models alongside their portable geometry and product structure. Run BREP and
STEP operations from an OpenCASCADE-enabled environment.


Component-Local Models
----------------------

An assembly-aware package can export each unique fabricated component once in
its canonical local coordinate frame. Repeated instances do not create
duplicate files, and purchased or raw-stock components are excluded by
default. Manufacturing STL is tessellated from the live OCC BREP, and strict
STEP exports that BREP analytically. Both therefore run during package
creation rather than from the portable display mesh stored in geometry JSON.

.. code-block:: bash

   conda run -n yapcad-brep python -m yapcad.dsl run \
     examples/product.dsl BUILD_ASSEMBLY \
     --package build/product.ycpkg --name product --version 0.2.0 \
     --component-export stl --component-export step --force

Files are written below ``exports/components/<component-id>/`` and registered
in both the manifest export list and their component definition. STEP export
through this package path is strict: missing analytic BREP is an error rather
than a faceted STEP or display-mesh STL fallback.


Strict Direct Exports
---------------------

For a direct command output, ``--strict-step`` requires analytic BREP and
disables faceted STEP fallback. ``--strict-stl`` tessellates the authoritative
BREP, removes zero-area facets, and requires a watertight mesh. These modes are
appropriate for manufacturing gates; ordinary display and preview workflows
may continue to use the portable mesh fallback.


Manufacturing Boundary
----------------------

Component export follows the engineering component definitions in the
package. Assembly-view compounds that overlap at intended interfaces are not
single printable components and can correctly fail strict STL validation.
Split such geometry into its actual fabricated component definitions before
export. Manufacturing exports are prototype fabrication inputs, not released
production drawings or semantic GD&T.
