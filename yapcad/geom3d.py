## geom3d, enhanced two- and three-dimensional geometry support for yapCAD
## started on Thu Oct  1 13:52:45 PDT 2020
## Richard W. DeVaul

from yapcad.geom import *

## the geometric representations of point, line, arc, poly, and
## geomlist provided by yapcad.geom are suitable for representing
## zero- or one-dimensional figures embedded in a three-dimensional
## space.  And while there is no requirement that the representations
## provided by yapcad.geom (such as an arc, line, or polyline) lie in
## the XY plane, many of the functions in yapcad.geom are explicitly
## intended to perform computational geometry operations on XY-planar
## entities. 

## Further, while a closed figure described by yapcad.geom may
## implicitly bound a two-dimensional face, there is little support
## for working with two-dimensional surfaces provided by that module.
## There is no direct support of any kind for working with
## three-dimemnsional volumes in yapcad.geom.

## in this module, we define the concept of a parametric
## two-dimensional surface and three-dimensional volume, and provide
## implicit and explicit geometry operations for working with them.

## The goal of the geom3d.yapcad module is to allow for the
## construction of two-dimensinal surfaces and three-dimensional
## geometry for the purposes of modeling, computational geometry, and
## rendering.  Specifically, we wish to support the following:

## (1) Support the implicit representation of three-dimensional
## geometry, and the performance of constructive solid geometry
## operations on this implicit geometry (union, intersection,
## difference) to produce more complex implicit three dimensional
## forms. 

## (2) Support the implicit representation of two dimensional
## surfaces, such as a planar surface specified by three points, or
## the surface of a three-dimensional object like a sphere, and allow
## for computational geometry operations on these surfaces, such as
## intersection operations, to produce explicit one-dimensional
## objects, such as lines, arcs, etc.

## (3) support the conversion of an implicit two-dimensional surface
## to an explicit, teselated triangluar geometry that may be easily
## rendered using a conventional 3D graphics rendering pipline, such
## as OpenGL

## (4) Support for the conversion of implicit three-dimenaional
## constructive solid geometry into an explicit, contiguous closed
## surface representation using the marching cubes algortihm, or any
## other user-specified conversion algoritm, for the purposes of
## interactive 3D rendering and conversion to 3D CAM formats, such as
## STL.

## Structure for 2D/3D geometry:

## geom3d in tuple( <type>, <transform>, <primary representation>,
##                  <sampled representation>, <rendering hints> )
## where:
## type in tuple ('surface' | 'solid' | 'transform', <subtype>)

