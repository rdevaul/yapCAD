# yapCAD Examples
A collection of example designs for [yapCAD](../README.md)

## Example List

* [example1.py](./example1.py) &mdash; simple DXF drawing example
* [example2.py](./example2.py) &mdash; simple computational geometry example demonstrarting computing the intersection of lines and lines, and lines and arcs.  No rendered output.
* [example3.py](./example3.py) &mdash; combines DXD drawing and computational geometry. Creates circles, computes tangent lines to circles, renders the output as a DXF drawing.
* [example4.py](./example4.py) &mdash; shows the use of the `Polygon` class to create, draw, and sample iso-distance lines offset from a central convex polygon. Uses sampled isodistance line to position regularly spaced (in terms of distance along that line) circles. Produces DXF rendered output.
* [example5.py](./example5.py) &mdash; another combined computational geometry and DXF drawing example, this time using "simple" vertex-list polygons rather than the full `Polygon` class.  Shows the creation, sampling, and drawing of such polygons. 
* [example6.py](./example6.py) &mdash; a combined computational geometry and DXF drawing example that reproduces the computation and output of [example5.py](./example5.py) but with rounded polygons using the `Polygon` class. 
