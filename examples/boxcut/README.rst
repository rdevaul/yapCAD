boxcut example
==============

An example procedural design system in **yapCAD** for creating
squeeze-fit box designs cut from flat sheets of material

how to use
----------

install **yapCAD**, and invoke the design system as follows:

::

   python3 boxcut.py <length> <width> <height> [options]

where options includes:

::

   -h, --help                         # show help message
   -u, --units {mm,inch}              # units for parameters, defult is mm
   -o, --outname <dxf file name base> # default is boxout
   -t, --thick <sheet thickness>      # default is 3.175 mm, or 0.125"
   -k, --kerf <kerf>                  # default is 0.397, or 1/64"
   -g, --gl                           # do a 3D visualization in openGL

For example, let’s say we want to create a box that is 6" in length, 4"
in width, and 3" in height. We are cutting quarter-inch stock acrylic
(the default for ``boxcut``) on a laser cutter without built-in kerf
compensation, so we need to specify the kerf, which happens to be about
0.5 mm.

First, let’s decide what name to give the output design file — let’s
call this design boxout6x4x3. We don’t need to specify material
thickness, since we are using the default. Finally, we want to see the
result visualized in 3D, so we will need to select that option.

Ok, we are ready run boxcut as follows:

::

   python3 boxcut 6.0 4.0 3.0 -u inch --outname boxout6x4x3 --kerf 0.5 --gl

This will cause an OpenGL window to pop up with a visualization of the
box, including instructions for manipulating the view. Hit ``esc`` to
exit.

The design file ``boxout6x4x3.dxf`` will also be produced. When you open
this file in your favorite DXF file viewer, you will see the profile for
each of the six pieces tiled and labeled.

The top surface of each piece as cut represent the outside surface of
each piece as assembled, except for the bottom piece which is right-side
up.
