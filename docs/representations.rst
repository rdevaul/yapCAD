===================================
Representations and boolean engines
===================================

Every yapCAD solid is a triangle mesh you can draw, measure and export.
Behind the mesh, a solid can carry one of two more precise definitions:
an OpenCASCADE BREP, or a signed distance field. Which definition a solid
has decides how booleans, fillets and exports treat it. This page explains
the three, which one is in charge of a given solid, and how yapCAD picks a
boolean engine.


Three representations
=====================

**Mesh.** Triangles, nothing more. Every solid has one, so every exporter,
renderer and measurement works on every solid. A mesh is only as accurate
as its triangles, and booleans between meshes are the hardest operation in
the library to make robust.

**BREP** (boundary representation). Exact surfaces -- planes, cylinders,
splines -- bounded by exact edges, held by OpenCASCADE. Primitives built
with OCC available carry one; so do STEP imports and OCC booleans and
fillets. A BREP gives analytic STEP export and exact fillets on any edge.
It needs ``pythonocc-core`` (see :doc:`installation`).

**Signed distance field** (SDF). A tree of primitives and operations that
gives the distance to the surface at any point (see :doc:`sdf_guide`).
Booleans, offsets, smooth blends, lattices and fillets of primitives are
exact and cannot fail. Fields need only numpy, and a field made of
primitives and hard booleans replays through OCC into an exact BREP.


Authority
=========

A solid's *authoritative* representation is the one that defines it; the
others are derived from it and can be regenerated.

* A solid built from a field is **SDF-authoritative**. Its mesh is a
  preview; its tree rides in the solid's ``construction`` record with the
  meshing parameters, so the preview can be regenerated exactly. A BREP
  replayed from the tree is marked as derived.
* A solid with an OCC BREP and no field is **BREP-authoritative**; its mesh
  is a tessellation.
* Anything else is a plain mesh.

Transforms keep a solid's authority: moving a field-authored solid moves
its tree too. Geometry JSON (``yapcad-geometry-json-v0.3``) records which
representation is authoritative, so a package reads back with the same
authority it was written with.

.. code-block:: python

    from yapcad import sdf
    from yapcad.brep import has_brep_data
    from yapcad.geom3d_util import prism

    field_part = sdf.to_solid(sdf.box(10.0), resolution=16)
    mesh_part = prism(10, 10, 10)

    sdf.is_sdf_solid(field_part)    # True
    sdf.is_sdf_solid(mesh_part)     # False
    has_brep_data(mesh_part)        # True when OCC is installed


Boolean engines
===============

:func:`yapcad.geom3d.solid_boolean` -- and so the DSL's ``union``,
``difference`` and ``intersection`` -- chooses an engine for each call:

1. **An engine you name wins.** ``engine=`` on the call, or the
   ``YAPCAD_BOOLEAN_ENGINE`` environment variable: one of ``sdf``, ``occ``,
   ``manifold``, ``native`` or ``trimesh``.
2. **Two fields combine as fields** (``sdf``). The result is exact, keeps
   its tree, and is meshed no coarser than either operand.
3. **Two BREPs use OpenCASCADE** (``occ``) when it is installed. If OCC
   fails, the call falls through to a mesh engine.
4. **Otherwise a mesh engine**: ``YAPCAD_MESH_BOOLEAN_ENGINE`` if set, else
   ``manifold`` when manifold3d is installed, else ``native``.

The mesh engines:

``manifold``
    manifold3d, called directly in double precision
    (``pip install 'yapCAD[manifold]'``). Fast and robust, including on
    shared and coincident faces. It refuses an operand that is not a
    closed 2-manifold; in automatic selection such a call falls back to
    ``native`` with a warning.

``native``
    yapCAD's own engine (``yapcad.boolean.csg``), in pure numpy. It splits
    triangles by the triangles they cross, classifies the pieces by
    generalised winding number, and repairs the result into a closed,
    conforming mesh. It is correct on the same cases as manifold3d --
    checked against it on 1,800 randomised cases -- but slower, and it is
    always available.

``trimesh``
    trimesh's boolean backends (Blender, OpenSCAD and others), selected
    with ``YAPCAD_TRIMESH_BACKEND``. Kept for compatibility.

A boolean between a field and a mesh, or a field and a BREP, goes to a
mesh engine and returns a mesh: the field is lost. Keep a design in one
representation where you can; the DSL's ``representation="sdf"`` switch
builds a whole design as fields.


Which to use
============

* **Mechanical parts with exact geometry and STEP deliverables**: BREP,
  with OCC installed -- or fields made of primitives and hard booleans,
  which replay to the same exact BREP.
* **Parts that need blends, lattices, offsets or shells, or no OCC**:
  fields.
* **Imported STL and other scanned or tessellated input**: meshes, with
  the ``manifold`` extra for booleans.
