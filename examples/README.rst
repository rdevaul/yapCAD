yapCAD Examples
===============

A collection of example designs for `yapCAD <../README.rst>`__

Example List
------------

-  `boxcut <./boxcut>`__ — parametric design system for creating
   squeeze-fit rectangular boxes from laser-cut acrylic or similar flat
   sheet material. Specify desired length, width, and height on the
   command line and an appropriate DXF file will be created. Compensates
   for laser-cutter kerf, optionally visualizes the resulting box in 3D.

-  `example1.py <./example1.py>`__ — simple drawing example, uses pyglet
   OpenGL interactive rendering. Click and drag to rotate drawing, zoom
   in and out with the up and down arrow keys, quit with ESC.

-  `example2.py <./example2.py>`__ — simple computational geometry
   example demonstrarting computing the intersection of lines and lines,
   and lines and arcs. No rendered output.

-  `example3.py <./example3.py>`__ — combines drawing and computational
   geometry. Creates circles, computes tangent lines to circles, renders
   the output as a DXF drawing.

-  `example4.py <./example4.py>`__ — shows the use of the ``Polygon``
   class to create, draw, and sample iso-distance lines offset from a
   central convex polygon. Uses sampled isodistance line to position
   regularly spaced (in terms of distance along that line) circles.
   Produces DXF rendered output.

-  `example5.py <./example5.py>`__ — another combined computational
   geometry and DXF drawing example, this time using “simple”
   vertex-list polygons rather than the full ``Polygon`` class. Shows
   the creation, sampling, and drawing of such polygons.

-  `example6.py <./example6.py>`__ — a combined computational geometry
   and DXF drawing example that reproduces the computation and output of
   `example5.py <./example5.py>`__ but with rounded polygons using the
   ``Polygon`` class.

-  `example7.py <./example7.py>`__ — draw and compute intersections
   between polylines, arcs, and lines.

-  `example8.py <./example8.py>`__ — create “flowers” with random
   rounded petals using the Polygon geometry generation class and mirror
   them around two axes.

-  `example8-gl.py <./example8-gl.py>`__ — same as previous, but render
   interactively using pyglet and OpenGL

-  `example9.py <./example9.py>`__ — three-dimensional geometry creation
   demonstration, rendered interactively using pyglet and OpenGL. Click
   and drag to rotate drawing, toggle lighting by hitting the ‘l’ key,
   zoom in and out with the up and down arrow keys, quit with ESC.

-  `example10.py <./example10.py>`__ — demo of inside testing for a
   variety of geometry. Generates 60 random geometric figures and 2,000
   random points, performs inside testing on all points with respect to
   all geometry. Draws inside points in aqua, outside in red. Optionally
   render interactively with pyglet OpenGL or DXF by way of command-line
   arguments. Run ``python3 example10.py help`` to see options.

-  `example11.py <./example11.py>`__ — demo of using boolean operations
   (intersection, union, difference) with yapCAD Polygon() instances to
   create simple combined shapes. Optionally render interactively with
   pyglet OpenGL or DXF by way of command-line arguments. Run
   ``python3 example11.py help`` to see options.

-  `example12.py <./example11.py>`__ — more complex demo of using
   boolean operations (intersection, union, difference) with yapCAD
   Polygon() instances to create combined shapes. Optionally render
   interactively with pyglet OpenGL or DXF by way of command-line
   arguments. Run ``python3 example12.py help`` to see options.
