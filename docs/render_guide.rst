Rendering fields: headless engineering images
=============================================

Agentic CAD work is gated by whether the tool can see what it built. Until
now, getting a picture out of yapCAD meant meshing the part or generating
BREP surfaces and handing the result to an external viewer. :mod:`yapcad.render`
removes that step for SDF-authored geometry: it sphere-traces the field in
numpy and writes a PNG, with no window, no GPU, no mesh and no new
dependency.

.. code-block:: python

   from yapcad import render, sdf

   part = sdf.subtract(sdf.chamfered_box((40.0, 30.0, 6.0), 1.0),
                       sdf.translate(sdf.cylinder(4.0, 20.0),
                                     (10.0, 0.0, 0.0)))
   render.image(part, view="iso", path="part.png")
   render.four_up(part, path="part_4up.png")
   render.section_view(part, axis="z", where=0.0, path="cut.png")

or from the command line, against a DSL design:

.. code-block:: console

   $ python -m yapcad.render design.dsl PLATE --view 4up -o plate.png
   PLATE: field built in 0.9s, 4up drawn in 1.1s at 640x480 -> plate.png

Where this sits
---------------

``docs/SDF-DESIGN.md`` section 6 plans two rendering tiers: mesh the field and
hand the triangles to an existing viewer (Tier 0), or emit GLSL and
sphere-trace it in a browser (Tier 1). Both put a human in front of a screen.
This is the tier in between, for a reader that is a *program*, and it is the
cheapest of the three: it carries none of the shader emitter's risk of
drifting away from the numpy backend, because it **is** the numpy backend.
:func:`yapcad.sdf.evaluate` and :func:`yapcad.sdf.normal` are the whole of it,
so ``smooth_union`` cannot come to mean two different things.

Sphere tracing works here for the same reason the mesher does. Hard CSG
produces a conservative *lower bound* on distance rather than an exact one,
so a march that advances by ``|f| / L`` can never step through a surface:
``exact=False`` costs iterations, never correctness.

Views
-----

The default projection is orthographic, because that is the projection
engineering drawings are read in -- a millimetre is the same number of pixels
everywhere, which is what makes the scale bar honest. Perspective is
available for a presentation view and is labelled *not to scale*.

.. list-table::
   :header-rows: 1
   :widths: 14 20 20 46

   * - ``view``
     - screen right
     - screen up
     - looks along
   * - ``iso``
     - --
     - ``+Z``-ish
     - ``(-1, -1, -1)``, the eye above ``+X +Y +Z``
   * - ``front``
     - ``+X``
     - ``+Z``
     - ``+Y``
   * - ``back``
     - ``-X``
     - ``+Z``
     - ``-Y``
   * - ``left``
     - ``-Y``
     - ``+Z``
     - ``+X``
   * - ``right``
     - ``+Y``
     - ``+Z``
     - ``-X``
   * - ``top``
     - ``+X``
     - ``+Y``
     - ``-Z``
   * - ``bottom``
     - ``+X``
     - ``-Y``
     - ``+Z``

:meth:`~yapcad.render.Camera.fit` frames the model's box rather than a sphere
around it, so a long part fills the frame instead of floating in it. Pass a
direction vector in place of a name for anything else.

Sections are not traced
-----------------------

Cutting a field is a single evaluation of the tree over a plane of points.
There is no marching, no meshing and no boolean -- the slice *is* the field
restricted to a plane -- so a section is the cheapest view there is, and it
carries information a mesh renderer cannot give you:

* outside the material the value is the **clearance** to the nearest
  surface, so the exterior of the slice is a clearance map;
* inside, its magnitude is the local **half-thickness**;
* the zero contour is the cut outline.

.. code-block:: python

   cut = render.section(part, axis="z", where=4.0)
   cut.area_mm2          # the material the plane passes through
   cut.distance[y, x]    # signed millimetres at that pixel

Scenes, parts and colour
------------------------

A scene is a list of parts, and each is traced against **its own** field and
composited by nearest depth. That is not only how parts get individual
colours and a legend: it is the thing that makes an assembly affordable.
Evaluation throughput falls as a tree grows -- a five-node part evaluates
about seven times faster per point than a hundred-node one -- so the cost of
a scene should grow with the number of parts in it, not with the size of one
combined tree.

.. code-block:: python

   render.image({"CRADLE": cradle, "PIN": pin}, view="iso", path="asm.png")

Colour is assigned to *groups* -- a material, a subassembly, whatever you key
on -- never to individual parts, and there are eight slots. A ninth group
folds into a neutral tone rather than inventing a hue that no reader could
separate from the eight above it. By default a two-to-eight-part scene
colours by part, because that is exactly the range where per-part colour
helps; past that, identity is carried by the part-boundary lines the
compositor draws and by the legend. Pass ``colour_by="group"``, ``"part"`` or
``"none"`` to decide for yourself.

What makes it an engineering image
----------------------------------

A reader that is a program cannot measure a pixel, hover a tooltip or orbit
the model, so everything it needs is drawn into the frame: the title, the
model's bounding box, a scale bar whose length is a round number of
millimetres, a labelled axis triad (an axis pointing at the eye is parked
beside the origin rather than drawn on top of it), a legend when there is
more than one group, and, on a section, the plane's coordinate and a
clearance key.

The line work matters more than the lighting. A shaded surface says roughly
which way a face points; a drawn crease says exactly where one face stops and
the next begins, which is what makes a chamfer, a gusset or a bolt hole
legible at a few hundred pixels. Four kinds of line are found together: the
silhouette against the background, the boundary between two parts, a depth
step within one part, and a crease where the normal turns sharply.

What it costs
-------------

Measured at 512x512 on one laptop core, against parts from the yapRover
design:

.. list-table::
   :header-rows: 1
   :widths: 34 14 14 38

   * - part
     - nodes
     - render
     - note
   * - a plate, a tub
     - 5
     - 0.1 s
     - the whole view
   * - a bracket with gussets and bolt holes
     - 105
     - 2.8 s
     - 0.76 M field evaluations
   * - the same, as a section
     - 105
     - 0.4 s
     - one ``evaluate`` call
   * - the same, four-up
     - 105
     - 12.5 s
     - three views and a section

Two things carry that, and both are on by default: rays are clipped to the
field's bounding box before marching, and each iteration evaluates only the
rays still in flight. Rays finish at wildly different times -- a silhouette
ray takes far longer than one landing on a flat face -- so without
compaction every iteration pays for the slowest ray in the frame. Turning
each off, on that 105-node part:

.. list-table::
   :header-rows: 1
   :widths: 18 18 22 42

   * - clipping
     - compaction
     - march
     - field evaluations
   * - on
     - on
     - 1.8 s
     - 0.76 M
   * - off
     - on
     - 5.4 s
     - 2.37 M
   * - on
     - off
     - 68 s
     - 36.96 M
   * - off
     - off
     - 79 s
     - 41.42 M

The bigger cost is usually next door. Building that 105-node tree through
the DSL took 11.7 s, against 2.8 s to render it, because the DSL meshes every
primitive and every boolean on the way through -- meshes this renderer never
looks at. Raising ``--cell`` makes the build cheaper and does not affect the
image at all, since the render uses the field rather than the mesh.

Extending it: the backend seam
------------------------------

Each part is traced into a :class:`~yapcad.render.PartBuffer` -- a depth and
a normal, and nothing about the rest of the scene -- and
:func:`~yapcad.render.composite` merges those by nearest depth into a
:class:`~yapcad.render.GBuffer` that also records which part won each pixel.
Shading, line work and annotation read only the merged buffer.

Keep that seam. A second backend that rasterises triangles from a tessellated
BREP or an imported mesh can fill the same :class:`~yapcad.render.PartBuffer`
and composite against field parts, and the silhouette where the two meet
falls out for free, because the edge pass never asks where a pixel came from.
Unions and placements -- which is most of an assembly -- need nothing more
than that. A boolean *across* representations is the case the compositor
cannot help with, and it is already resolved before the renderer sees it:
the modeller has promoted it to a mesh by then.

Nearest-wins is decided with a tolerance rather than exactly. Datum-mated
parts share a face, and an exact comparison between two depths that differ
only by the float noise a transform chain leaves behind decides each pixel
independently, which reads as speckle down every flush joint. Note that
quantising depth into bins does not fix this -- two depths either side of a
bin edge still split -- which is why the comparison is a tolerance.

Limits
------

* **SDF-authored geometry only.** A STEP import or a mesh-only part has no
  field to trace until Phase 4 promotion lands. :func:`~yapcad.render.as_parts`
  says so rather than rendering nothing.
* **Grazing rays stop short.** The march stops where the field is within the
  hit tolerance of zero, which bounds the *field*, not the depth. Near a
  silhouette the field changes slowly along the ray, so the reported depth
  can be short by much more than the tolerance. It always errs toward the
  camera, never through the surface.
* **A feature thinner than the hit tolerance can be missed** by a
  near-tangential ray -- the same minimum-feature-size story as meshing.
* **These are images for seeing, not for measuring.** Numbers belong to
  :mod:`yapcad.rom` and the clearance APIs.
* **Not interactive.** At a tenth of a second to a few seconds a frame this
  is fine for a program and wrong for a human orbiting a model; that is what
  the browser tier in section 6 is for.
* **Lettering is capitals only.** The built-in 5x7 face has no lowercase and
  :func:`~yapcad.render.raster.text` upper-cases its input, which is how
  drawings are lettered anyway.

Reference
=========

.. autoclass:: yapcad.render.Camera
   :members:
.. autoclass:: yapcad.render.Part
.. autoclass:: yapcad.render.PartBuffer
   :members:
.. autoclass:: yapcad.render.GBuffer
   :members:
.. autoclass:: yapcad.render.Section
   :members:
.. autoclass:: yapcad.render.Theme
   :members:

.. autofunction:: yapcad.render.image
.. autofunction:: yapcad.render.four_up
.. autofunction:: yapcad.render.section_view
.. autofunction:: yapcad.render.render
.. autofunction:: yapcad.render.render_part
.. autofunction:: yapcad.render.composite
.. autofunction:: yapcad.render.shade
.. autofunction:: yapcad.render.edges
.. autofunction:: yapcad.render.section
.. autofunction:: yapcad.render.as_parts
.. autofunction:: yapcad.render.scene_bounds
.. autofunction:: yapcad.render.parts_from_dsl
.. autofunction:: yapcad.render.annotate
.. autofunction:: yapcad.render.grid
.. autofunction:: yapcad.render.write_png
