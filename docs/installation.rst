============
Installation
============

yapCAD needs Python 3.11 or newer. The core -- 2D geometry, mesh solids,
signed distance functions, the DSL, packages and the assembly system --
installs from PyPI with pip and needs only numpy and a few small
libraries. OpenCASCADE, for exact BREP solids and analytic STEP, is
optional and comes from conda-forge.


Choosing an install
===================

=========================================  =======================================================
You want                                   Install
=========================================  =======================================================
Meshes, SDF parts, the DSL, STL export     ``pip install yapCAD``
Fast, robust mesh booleans                 ``pip install 'yapCAD[manifold]'``
Strict watertightness checks and repair    ``pip install 'yapCAD[meshcheck]'``
BREP solids, fillets of BREPs, STEP        the conda environment below
Development and tests                      ``pip install -e '.[tests]'`` from a checkout
=========================================  =======================================================

Extras combine: ``pip install 'yapCAD[manifold,meshcheck]'``.

Without OpenCASCADE, yapCAD is still a complete solid modeller for
SDF-authored parts: fields combine exactly, ``fillet`` rounds field
primitives, fasteners and gears have field versions, and a DSL design
builds as fields with ``representation="sdf"`` (see :doc:`sdf_guide`).
Mesh booleans run on the native engine, or on manifold3d when the
``manifold`` extra is installed.

.. note::

   The ``brep`` and ``full`` extras name ``pythonocc-core``, but it is
   not published on PyPI, so pip cannot install them. Use conda for OCC.


With OpenCASCADE (conda)
========================

From a checkout:

.. code-block:: bash

    conda env create -f environment.yml
    conda activate yapcad-brep

``environment.yml`` installs ``pythonocc-core`` 7.7.2 from conda-forge,
yapCAD in editable mode with its test extra, trimesh, gmsh and
manifold3d. ``environment-full.yml`` adds FEniCSx for finite-element
analysis.

Check what is available:

.. code-block:: python

    from yapcad.brep import occ_available
    print(occ_available())      # True in the conda environment

Without OCC, importing yapCAD prints a ``YapcadBrepUnavailableWarning``
that lists what needs it. Everything else works.


Development
===========

.. code-block:: bash

    python -m venv .venv && source .venv/bin/activate
    pip install -e '.[tests]'
    pytest

The suite runs with or without OCC; tests that need it are skipped and
carry the ``requires_occ`` marker, so ``pytest -m "not requires_occ"``
runs the pure-Python set explicitly. CI tests the pure-Python install on
Python 3.11 and 3.14 and the OCC environment on 3.12.


Environment variables
=====================

=================================  ======================================================================
``YAPCAD_BOOLEAN_ENGINE``          force a boolean engine: ``sdf``, ``occ``, ``manifold``, ``native``,
                                   ``trimesh`` (see :doc:`representations`)
``YAPCAD_MESH_BOOLEAN_ENGINE``     the engine for mesh operands when none is forced
``YAPCAD_TRIMESH_BACKEND``         the backend trimesh uses when the ``trimesh`` engine is chosen
``YAPCAD_DSL_REPRESENTATION``      ``mesh`` (default) or ``sdf``: what DSL primitives produce
``YAPCAD_DSL_SDF_CELL_MM``         the cell size DSL primitives are meshed at in ``sdf`` mode
``YAPCAD_SDF_SIMPLIFY``            ``0`` keeps finished SDF parts' uniform meshes
``YAPCAD_DSL_RECURSION_LIMIT``     depth limit for DSL command-to-command calls (default 100)
``YAPCAD_STEP_FORMAT``             ``faceted`` (default) or ``analytic``: STEP from ``yapcad.dsl run``
``YAPCAD_FASTENER_DATA``           extra fastener catalog directories, searched before the bundled data
``YAPCAD_ASSEMBLY_RESOLVER``       assembly mate resolver: ``basic``, ``staged`` or an entry point
=================================  ======================================================================
