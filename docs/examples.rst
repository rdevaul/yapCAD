========
Examples
========

Two example scripts show signed distance fields at work. Both run on a
plain ``pip install yapCAD``, with no OpenCASCADE; both write STL meshes,
PNG previews and a JSON summary. Other examples are listed in the README.


Parts built as fields
=====================

``examples/sdf_parts_demo.py`` builds the parts a mechanism needs, the ways
a yapCAD user would build them, and reports each part's triangle count
before and after simplification and its mass in PETG.

.. code-block:: bash

    python examples/sdf_parts_demo.py                  # all parts, to build/sdf-parts
    python examples/sdf_parts_demo.py --only m8-nut    # one part
    python examples/sdf_parts_demo.py --tolerance 0    # skip simplification
    python examples/sdf_parts_demo.py --list

.. list-table::
   :widths: 40 60
   :class: borderless

   * - .. image:: images/sdf/yaprover-wheel.jpg
          :width: 260px
          :alt: the yapRover wheel built as a signed distance field
     - **yaprover-wheel.** The yapRover wheel, written in the DSL exactly as
       the rover's design writes it and built with
       ``representation="sdf"``. Its rounded rim is an exact field fillet,
       so no OCC is needed. 466,032 triangles uniform, 37,416 simplified;
       231 g in PETG.
   * - .. image:: images/sdf/gyroid-wheel.jpg
          :width: 260px
          :alt: the wheel with a gyroid lattice web
     - **gyroid-wheel.** The same wheel with a gyroid lattice in place of the
       spokes, blended into the rim and hub with smooth unions: a part that
       is a few lines as a field and very hard as a BREP. The lattice is
       curved everywhere, so simplification keeps 40% of its triangles;
       329 g in PETG.
   * - .. image:: images/sdf/m8-nut.jpg
          :width: 260px
          :alt: an M8 hex nut with a helical internal thread
     - **m8-nut.** An ISO 4032 M8 nut from
       ``metric_hex_nut("M8", representation="sdf")``: a hexagonal prism less
       a helical thread field whose profile is the mesh fastener's own.
       96,988 triangles uniform, 3,450 simplified.
   * - .. image:: images/sdf/m8x30-bolt.jpg
          :width: 260px
          :alt: an M8x30 hex bolt
     - **m8x30-bolt.** An M8x30 hex bolt: external thread, shank, washer face
       and head. 207,584 triangles uniform, 8,882 simplified.
   * - .. image:: images/sdf/miter-gear.jpg
          :width: 260px
          :alt: a 24-tooth miter gear
     - **miter-gear.** The yapRover differential's 24-tooth, module 1.5
       straight bevel gear: the outer tooth section carried to the pitch
       apex, the same construction as the OCC generator, whose exact BREP it
       replays to. 118,192 triangles uniform, 5,262 simplified.
   * - .. image:: images/sdf/herringbone-gear.jpg
          :width: 260px
          :alt: a 20-tooth herringbone gear
     - **herringbone-gear.** A 20-tooth involute gear twisted 25° each way
       from mid face. Its tooth tips' acute edges are why the mesher places a
       vertex per surface component of a cell. 218,128 triangles uniform,
       5,566 simplified.


Field primitives and operators
==============================

``examples/sdf_demo.py`` is a gallery of what fields are made of. Each
model reports its volume, whether its mesh is a closed manifold, whether
its field is an exact distance, its Lipschitz bound, and whether it
replays as exact CSG. ``--step`` also writes analytic STEP for every model
that does, when OCC is installed.

.. code-block:: bash

    python examples/sdf_demo.py                       # to build/sdf-demo
    python examples/sdf_demo.py --step                # and analytic STEP
    python examples/sdf_demo.py --only cube-frame --no-render

.. list-table::
   :widths: 33 33 34
   :class: borderless

   * - .. figure:: images/sdf/cube-frame.jpg
          :width: 200px

          **cube-frame**: a box less a sphere that breaches all six faces
     - .. figure:: images/sdf/box.jpg
          :width: 200px

          **box**: sharp edges kept exactly by the vertex solve
     - .. figure:: images/sdf/sphere.jpg
          :width: 200px

          **sphere**: the smoothness baseline
   * - .. figure:: images/sdf/torus.jpg
          :width: 200px

          **torus**: genus one
     - .. figure:: images/sdf/bracket-plate.jpg
          :width: 200px

          **bracket-plate**: hard CSG, replayable to analytic STEP
     - .. figure:: images/sdf/smooth-union.jpg
          :width: 200px

          **smooth-union**: a field-native blend
   * - .. figure:: images/sdf/shell-cutaway.jpg
          :width: 200px

          **shell-cutaway**: a 2 mm wall cut open by a half space
     - .. figure:: images/sdf/gyroid-lattice.jpg
          :width: 200px

          **gyroid-lattice**: an infinite lattice bounded by a box
     - .. figure:: images/sdf/lattice-in-a-part.jpg
          :width: 200px

          **lattice-in-a-part**: a solid skin over a lattice core

See :doc:`sdf_guide` for how these are built.
