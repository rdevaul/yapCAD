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
``max``, a fillet or chamfer a finished primitive or a smooth blend. They
cannot fail, have no coplanar-face or sliver cases, and stay exact. A part
authored as a field is meshed when it has to be, at whatever resolution the
job needs, and the field travels with the mesh so it can be meshed again.

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
``chamfered_box(size, distance)``      box with every edge bevelled 45 degrees
``cylinder(radius, height)``           capped cylinder on the z axis
``rounded_cylinder(r, height, edge)``  cylinder with both circular edges rounded
``chamfered_cylinder(r, height, d)``   cylinder with both circular edges bevelled
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
(Manifold Dual Contouring), so the mesh stays manifold.

When the cell turns out to be too coarse for the part -- the mesh comes back
empty, or dual contouring cannot close it -- ``to_solid`` halves the cell and
meshes again, up to 256 cells across, and says so:

.. code-block:: text

    RuntimeWarning: to_solid: dual contouring found no surface at cell 2 mm
    (resolution 48), so it refined to cell 1 mm (resolution 96), where the
    mesh came out closed with 992 vertices. Pass cell_mm=1 to mesh it in one
    pass.

That matters most for booleans, which can leave a feature thinner than
anything either operand had -- a cut that leaves a 2 mm remnant of a 12 mm
gusset -- where no rule reading the operands could have known. The cell the
warning reports is the coarsest one that *found* the surface, not one that
measures it well: a remnant a cell and a bit thick meshes, and its volume is
still rough. Pass it back as ``cell_mm`` to skip the extra pass, or pass
something finer if the dimension matters:

.. code-block:: python

    part = sdf.to_solid(tree, cell_mm=1.0)       # no refinement pass
    exact = sdf.to_solid(tree, cell_mm=0.25)     # and resolved, not just found

``cell_mm`` is the cell size in model units, rounded down to a whole number
of cells across the longest axis, so the cell is never coarser than asked.
``sdf.combine_all(solids, op, cell_mm=...)`` passes one to a boolean, and
the cell a part was meshed at rides in its meshing record, so a later
boolean on a refined part starts from the cell that worked.

A field that holds nothing is not refined for: when the Lipschitz bound
proves the region is empty -- ``sdf.proves_empty`` is the test --
the result is the empty solid yapCAD allows for an empty CSG result, and it
still carries its tree. Only a region that cannot be proved empty is
refined, which is why an empty boolean costs one pass rather than an
escalation. ``refine=False`` restores the old behaviour: ``to_solid`` then
raises ``sdf.NonManifoldMeshError``, or an ``SdfError`` naming the cell it
tried, rather than mesh again.

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


Fillets, chamfers and compounds
===============================

``sdf.fillet(node, radius)`` rounds every edge, as OCC's fillet-all-edges
does, and is exact for the analytic primitives: a box becomes a
``rounded_box`` and a cylinder a ``rounded_cylinder``. Both still replay
through OCC, so a filleted part keeps its analytic STEP export. The DSL's
``fillet`` uses it for any SDF-authoritative solid; no OCC is needed.

.. code-block:: python

    rounded = sdf.fillet(sdf.box((40.0, 16.0, 28.0)), 7.0)
    assert rounded == sdf.rounded_box((40.0, 16.0, 28.0), 7.0)

``sdf.chamfer(node, distance)`` is the same thing with a 45 degree bevel in
place of the round, as OCC's chamfer-all-edges does -- the finish a
supportless FDM part wants, since a bevel prints over thin air and a
fillet's tangent does not. A box becomes a ``chamfered_box``, a cylinder a
``chamfered_cylinder``, and both replay through OCC as well.

.. code-block:: python

    bevelled = sdf.chamfer(sdf.cylinder(18.0, 20.0), 1.5)
    assert bevelled == sdf.chamfered_cylinder(18.0, 20.0, 1.5)

A ``chamfered_cylinder`` is an exact field; a ``chamfered_box`` is the
intersection of four exact fields, so like any hard CSG it is exact inside
and on the surface and a Lipschitz-1 *lower bound* outside -- within 7% of
the true distance at worst, which makes a clearance measured from it
conservative rather than optimistic.

Fillet or chamfer primitives before combining them. Either one applied to a
combined field is refused, because finishing each operand would also finish
edges the boolean removes; the general case waits on sampled fields.

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
``fillet``, ``chamfer`` and ``compound`` keep fields. Other primitives -- region
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

**A picture**: :mod:`yapcad.render` sphere-traces the field straight to a
PNG, so looking at a part needs neither a mesh nor an external viewer:

.. code-block:: python

    from yapcad import render

    render.four_up(plate, path="plate.png")

A planar section is one ``evaluate`` call over a plane -- no marching, no
mesh -- and shows clearance outside the material and wall thickness inside
it. See :doc:`render_guide`.



Building without meshing
========================

Meshing is the expensive half of field authoring -- about seven eighths of
the time to build a part -- and a great deal of what a design does with a
part afterwards never needs the triangles. Placing it by a datum mate
transforms its tree; a range-of-motion sweep measures clearance against the
tree; :mod:`yapcad.render` traces the tree. All of that works on a part
whose mesh has never been built.

``sdf_mesh=False`` builds the fields and defers every mesh:

.. code-block:: python

    result = compile_and_run(source, "BUILD", {}, representation="sdf",
                             sdf_mesh=False)

or, around any field code:

.. code-block:: python

    with sdf.deferring():
        plate = sdf.to_solid(field)          # no triangles are made

What comes back carries the field and a *plan* -- the cell it would have
meshed at, recorded in cells across the longest axis so that it survives a
transform. :func:`yapcad.sdf.realize` builds the mesh in place when
something wants it, and :func:`yapcad.sdf.realize_all` does a whole design:

.. code-block:: python

    sdf.realize(plate)                # this part
    sdf.realize_all(result.geometry)  # all of them, simplified as usual

Deferring is not only postponement. A design meshes every intermediate
boolean on the way through and throws those meshes away; deferring meshes
only what survives. On a 105-node bracket, building eagerly took 10.6 s,
while deferring and then realizing the one part that mattered took 5.2 s --
for a bit-identical result.

What refuses, and why
---------------------

An empty surface list is already legal in yapCAD: it is what an empty CSG
result looks like. So the code that consumes solids was written when the
only way to get one was to mean it, and it answers accordingly --
``issolidclosed`` calls an empty solid closed, ``volumeof`` calls it zero,
``solid_contains_point`` puts every point outside. Those are reasonable
answers about nothing at all, and wrong about a part whose mesh has not been
built yet.

So a deferred solid never answers from an absent mesh. Either the field
answers the question, or the call refuses and names the cure:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - asked for
     - what happens
   * - ``solidbbox``, ``center``, ``bbox``
     - answered from the tree's bounds, no meshing
   * - ``solid_contains_point``, ``solids_intersect``
     - answered from the field, which is *exact* where ray casting
       against triangles has to pick a tolerance
   * - a transform, a datum mate
     - moves the field and keeps the plan
   * - ``volumeof``, ``issolidclosed``, ``solid_shells``
     - ``DeferredMeshError`` naming :func:`yapcad.sdf.realize`
   * - STL, STEP, a mesh boolean, a viewer
     - the same refusal, raised once in ``mesh_view``

A genuinely empty solid is unaffected: it is not waiting for anything, and
still answers as it always did.

Caveats
-------

* Deferral is opt-in, and meshing remains the default.
* ``brep=True`` and deferral ask for opposite things, and say so.
* A deferred part is not simplified until it is realized;
  :func:`~yapcad.sdf.realize_all` simplifies as the eager path does, while
  :func:`~yapcad.sdf.realize` leaves the mesh alone unless asked.


Limits
======

* **Mixed booleans** of a field with a mesh or BREP give a mesh, not a
  field. Promoting meshes and BREPs to fields is planned.
* **Fillets** of combined fields and of mesh-only solids are refused for
  now; both wait on sampled fields.
* **Meshing is uniform.** A part with fine detail is meshed finely
  everywhere, which simplification then undoes; adaptive meshing would
  avoid building the large mesh at all. Refinement is uniform for the same
  reason: a thin feature anywhere takes the whole part, or the whole
  compound, to a finer cell.
* **Refinement finds a feature; it does not resolve it.** It stops at the
  first cell that meshes, which can still measure a thin remnant roughly.
  A dropped *part* of a mesh -- some other feature missed while the rest
  came out fine -- is not detected at all; only a wholly empty or unclosed
  mesh is.
* **Meshing intermediate results**: every boolean meshes its result, even
  when the next step discards it.

:doc:`SDF-DESIGN` tracks each of these.
