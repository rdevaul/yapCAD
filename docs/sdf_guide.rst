====================================
Signed distance functions: the guide
====================================

A signed distance function (SDF) describes a solid by a field: at every
point in space it gives the distance to the part's surface, negative inside
and positive outside. yapCAD treats fields as a third way to author solids,
alongside meshes and OpenCASCADE BREPs, and builds them with nothing but
numpy.

Fields are worth having because the operations that are fragile on meshes
and BREPs are arithmetic on fields. A union is a ``min``, a difference a
``max``, a fillet a rounded primitive or a smooth blend. They cannot fail,
have no coplanar-face or sliver cases, and stay exact. A part authored as a
field is meshed when it has to be, at whatever resolution the job needs,
and the field travels with the mesh so it can be meshed again.

This guide covers authoring fields in Python and in the DSL, turning them
into solids, and getting them out as STL, STEP and packages. The design
reasoning, including the representation and the Lipschitz discipline,
is in :doc:`SDF-DESIGN`.

.. contents:: On this page
   :local:
   :depth: 1


A first part
============

.. code-block:: python

    from yapcad import sdf
    from yapcad.geom3d import issolidclosed, volumeof
    from yapcad.io.stl import write_stl

    # A 40 x 40 x 10 mm plate with rounded edges, one bore and two holes.
    plate = sdf.subtract(
        sdf.rounded_box((40.0, 40.0, 10.0), 3.0),
        sdf.translate(sdf.cylinder(6.0, 20.0), (10.0, 0.0, 0.0)),
        sdf.translate(sdf.cylinder(3.0, 20.0), (-12.0, -12.0, 0.0)),
        sdf.translate(sdf.cylinder(3.0, 20.0), (-12.0, 12.0, 0.0)),
    )

    solid = sdf.to_solid(plate, resolution=96)
    assert issolidclosed(solid)
    print(round(volumeof(solid), 1))        # about 13,640 mm^3
    write_stl(solid, "plate.stl")

``plate`` is a *field*: an immutable tree of nodes. ``sdf.to_solid`` meshes
it into an ordinary yapCAD solid, and everything in yapCAD that takes a
solid takes this one: booleans, transforms, exporters, packages. The solid
also carries the tree, so it is *SDF-authoritative*: yapCAD knows the mesh
is a preview of the field and can regenerate it.

Every node reports what it is:

.. code-block:: python

    plate.bounds        # ((-20, -20, -5), (20, 20, 5))
    plate.exact         # False: hard CSG bounds the distance, not equals it
    plate.lipschitz     # 1.0: |f| never changes faster than distance
    plate.csg_exact     # True: replayable through OCC as exact CSG
    sdf.evaluate(plate, [[0, 0, 0], [30, 0, 0]])   # [-4.0, 10.0]


Building fields
===============

**Primitives** are centred on the origin.

=====================================  ==================================================
``sphere(radius)``                     sphere
``box(size)``                          box; ``size`` is a scalar or ``(x, y, z)``
``rounded_box(size, radius)``          box with every edge rounded, same outer extent
``cylinder(radius, height)``           capped cylinder on the z axis
``rounded_cylinder(r, height, edge)``  cylinder with both circular edges rounded
``cone(r_bottom, r_top, height)``      frustum or cone on the z axis
``capsule(start, end, radius)``        a segment swept by a sphere
``torus(major, minor)``                torus in the xy plane
``half_space(normal, offset)``         everything on one side of a plane (unbounded)
``gyroid(period, thickness)``          gyroid lattice sheet (unbounded)
``schwarz_p(period, thickness)``       Schwarz P lattice sheet (unbounded)
=====================================  ==================================================

**Combinators.**

====================================  =================================================
``union(*parts)``                     the region any part occupies
``intersect(*parts)``                 the region every part occupies
``subtract(target, *tools)``          ``target`` less every tool
``compound(*parts)``                  a union for booleans; separate bodies when meshed
``smooth_union(a, b, radius)``        a union blended over ``radius``
``smooth_intersect(a, b, radius)``    an intersection blended over ``radius``
``smooth_subtract(a, b, radius)``     a difference blended over ``radius``
``offset(part, distance)``            grow (or, negative, shrink) a part
``shell(part, thickness)``            a wall of ``thickness`` straddling the surface
====================================  =================================================

**Transforms**: ``translate(part, delta)``, ``rotate(part, axis, degrees)``,
``scale(part, factor)`` and ``transform(part, matrix)``. A field moved by
``geom3d.translatesolid`` and the other solid transforms carries its tree
along, so the definition always matches the mesh.

**Profiles** for parts built from a 2D outline:

* ``polygon(points, symmetry=1)`` -- the exact distance to a closed 2D
  polygon, constant along z. Pass the rotational symmetry (a gear's tooth
  count) and evaluation tests only the nearby segments.
* ``extrude(profile, height)`` -- extrude a profile about ``z = 0``.
* ``apex_extrude(profile, z_ref, z_lo, z_hi)`` -- the cone over a profile
  toward the origin, as a bevel gear's teeth are built.
* ``twist(part, rate, fold_z=None)`` -- rotate each section by ``rate * z``
  radians, or rise to a fold and back for a herringbone.

Lattices are infinite, so intersect them with a solid before meshing:

.. code-block:: python

    cube = sdf.intersect(sdf.box(24.0), sdf.gyroid(9.0, 0.55))


Exactness and the Lipschitz bound
---------------------------------

A field is *exact* when its value is the true distance to the surface. The
primitives are; most combinations are not, because ``min`` and ``max`` of
distances only bound the distance away from the surface. What every yapCAD
field guarantees is a Lipschitz bound ``L``: ``|f|`` never changes faster
than ``L`` times the distance moved. That bound is what lets the mesher
skip empty space safely and a ray tracer step without tunnelling, so every
node computes it. Domain operators such as ``twist`` raise it. Nothing in
the API needs it set by hand.


Turning fields into solids
==========================

``sdf.to_solid(node, resolution=64)`` dual-contours the field over its
bounds. ``resolution`` is the number of cells across the longest axis, so
the cell size is the bounds' longest extent divided by it. Two rules
choose it:

* **A cell must be smaller than the thinnest feature**, or the surface
  crosses one cell twice. A wall a cell thick usually survives; a thread
  crest 0.16 mm wide needs cells of about 0.14 mm.
* **Sharp edges are kept.** Each vertex is placed by solving for the point
  that best fits the surface normals in its cell, so a box's edges stay
  sharp at any resolution.

Where the surface does cross a cell twice -- the thin wedge at an acute
edge, a lattice wall -- the cell gets one vertex per piece of surface
(Manifold Dual Contouring), so the mesh stays manifold. ``to_solid`` checks
the result and raises ``sdf.NonManifoldMeshError`` rather than return a
broken solid. Raise the resolution if it does.

Meshing uniformly makes more triangles than a part needs: an M8 nut fine
enough for its threads is just as fine across its flat faces.
``sdf.simplify_solid`` merges triangles wherever the field confirms the
merged surface stays within a tolerance of the true part -- by default a
twentieth of the cell the part was meshed at:

.. code-block:: python

    from yapcad.fasteners import metric_hex_nut

    uniform = metric_hex_nut("M8", representation="sdf", simplify=False)
    light = sdf.simplify_solid(uniform)    # about 97,000 -> 3,500 triangles
    assert issolidclosed(light)

The simplified mesh is closed, no less accurate than the uniform one, and
records the tolerance alongside the other meshing parameters.

**Finished parts are simplified by default**: the fastener and gear
builders, the solids a field-mode DSL call returns, and the parts of a
packaged field-mode assembly. Intermediate results -- every boolean -- are
not, because simplifying costs about twenty times the meshing and those
meshes are usually discarded. Pass ``simplify=False`` (fasteners, gears),
``sdf_simplify=False`` (``compile_and_run``, ``package_from_dsl``),
``--no-simplify`` on the DSL command line, or set ``YAPCAD_SDF_SIMPLIFY=0``
to keep the uniform meshes. ``sdf.simplify_finished`` applies the same rule
to your own parts.


Booleans
========

When both operands of :func:`yapcad.geom3d.solid_boolean` are
SDF-authoritative, the boolean combines their trees and meshes the result,
no coarser than either operand was meshed. The result stays a field, so it
can be combined, filleted or remeshed again. The DSL's ``union``,
``difference`` and ``intersection`` all go through ``solid_boolean``.

.. code-block:: python

    from yapcad.geom3d import solid_boolean

    block = sdf.to_solid(sdf.box(20.0), resolution=40)
    ball = sdf.to_solid(sdf.sphere(12.0), resolution=48)
    hollowed = solid_boolean(block, ball, "difference")
    assert sdf.is_sdf_solid(hollowed)

An operand that is a mesh or a BREP sends the boolean to the mesh engines
instead (see :doc:`representations`); the result is then a mesh.


Fillets and compounds
=====================

``sdf.fillet(node, radius)`` rounds every edge, as OCC's fillet-all-edges
does, and is exact for the analytic primitives: a box becomes a
``rounded_box`` and a cylinder a ``rounded_cylinder``. Both still replay
through OCC, so a filleted part keeps its analytic STEP export. The DSL's
``fillet`` uses it for any SDF-authoritative solid; no OCC is needed.

.. code-block:: python

    rounded = sdf.fillet(sdf.box((40.0, 16.0, 28.0)), 7.0)
    assert rounded == sdf.rounded_box((40.0, 16.0, 28.0), 7.0)

Fillet primitives before combining them. A fillet of a combined field is
refused, because rounding each operand would also round edges the boolean
removes; a correct general fillet waits on sampled fields.

``sdf.compound`` keeps several bodies separate, as the DSL's ``compound``
does: they combine as a union in booleans but mesh as distinct bodies.


Fasteners and gears
===================

The fastener and gear generators build fields on request:

.. code-block:: python

    from yapcad.fasteners import metric_hex_bolt, metric_hex_nut
    from yapcad.gears import make_involute_gear_sdf, make_straight_bevel_gear_sdf
    from yapcad.gears.bevel import StraightBevelGearSpec

    nut = metric_hex_nut("M8", representation="sdf")
    bolt = metric_hex_bolt("M8", 30.0, representation="sdf")
    miter = make_straight_bevel_gear_sdf(StraightBevelGearSpec(
        teeth=24, mate_teeth=24, outer_module_mm=1.5, face_width_mm=8.0,
        bore_diameter_mm=8.4))
    herringbone = make_involute_gear_sdf(20, 1.5, 8.0, helix_angle_deg=25.0,
                                         herringbone=True)

Threads are true helices whose profile comes from the same definition as
the mesh fasteners. A bevel gear's field is the BREP generator's own
construction, and replays to that generator's exact BREP when OCC is
installed. Each part's tree is a few hundred bytes, because the nodes store
the specification -- teeth, module, pitch -- rather than the geometry.


Fields in the DSL
=================

An existing DSL design builds as fields with one setting:

.. code-block:: python

    from yapcad.dsl import compile_and_run

    source = '''
    module plate
    command PART() -> solid:
        let blank: solid = fillet(box(40.0, 30.0, 6.0), 1.0)
        emit difference(blank, translate(cylinder(4.0, 10.0), 8.0, 0.0, -5.0))
    '''
    result = compile_and_run(source, "PART", {}, representation="sdf")
    assert result.success and sdf.is_sdf_solid(result.geometry)

In ``"sdf"`` mode ``box``, ``cylinder``, ``sphere``, ``cone``, ``tube`` and
``spherical_shell`` produce fields placed exactly as their mesh versions
are; the fasteners and gears produce their field versions; and booleans,
``fillet`` and ``compound`` keep fields. Other primitives -- region
extrusions, sweeps, lofts -- still produce meshes. Set it with:

* ``compile_and_run(..., representation="sdf", sdf_cell_mm=0.5)``
* ``package_from_dsl(..., representation="sdf")``
* ``python -m yapcad.dsl run design.dsl PART --representation sdf``
* ``YAPCAD_DSL_REPRESENTATION=sdf`` (and ``YAPCAD_DSL_SDF_CELL_MM``) for
  tools that do not pass it on.

``sdf_cell_mm`` sets the cell size primitives are meshed at; a primitive
is never meshed coarser than a quarter of its thinnest extent. The yapRover
release design builds this way with no OCC: every structural part comes
out closed and within 0.1% of the exact volume OCC computes.


Getting parts out
=================

**STL**: ``yapcad.io.stl.write_stl(solid, path)``, as for any solid.

**STEP**: a tree of primitives, hard booleans, rigid transforms and the
rounded primitives is *CSG-exact*: it can be replayed through
OpenCASCADE into an exact BREP, and exported as analytic STEP. Ask for it
when meshing, then export:

.. code-block:: python

    from yapcad.io.step import write_step_analytic

    exact = sdf.to_solid(plate, resolution=96, brep="auto")
    # True: analytic STEP.  Without OCC there is no BREP to replay, so it
    # writes a faceted STEP and returns False.
    write_step_analytic(exact, "plate.step")

``brep=True`` requires a replay and raises if it is impossible;
``"auto"`` replays when it can. Smooth blends, lattices, offsets and the
other non-CSG nodes cannot replay. ``sdf.csg_blockers(node)`` lists the
nodes responsible before anyone tries an export.

**Packages and geometry JSON**: an SDF-authoritative solid serialises with
its tree as the authoritative representation
(``yapcad-geometry-json-v0.3``), its mesh as a preview, and any replayed
BREP marked as derived. Reading it back restores the field.


Limits
======

* **Mixed booleans** of a field with a mesh or BREP give a mesh, not a
  field. Promoting meshes and BREPs to fields is planned.
* **Fillets** of combined fields and of mesh-only solids are refused for
  now; both wait on sampled fields.
* **Meshing is uniform.** A part with fine detail is meshed finely
  everywhere, which simplification then undoes; adaptive meshing would
  avoid building the large mesh at all.
* **Meshing intermediate results**: every boolean meshes its result, even
  when the next step discards it.

:doc:`SDF-DESIGN` tracks each of these.
