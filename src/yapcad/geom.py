## foundational computational geometry library for yapCAD
## Born on 29 July, 2020
## Copyright (c) 2020 Richard DeVaul

# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation files
# (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
# BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""foundational computational geometry library for **yapCAD**

====================
OVERVIEW
====================

The yapcad.geom module provides foundational computational geometry
functions for **yapCAD**.  This includes operations for computing
intersections, bounding boxes, inside-outside testing, *etc.*, on both
simple and compound geometry.

=====================================
yapcad.geom geometric representations
=====================================

The yapcad.geom module includes constants, scalar operations, vector
operations, operations on non-compound geometry elements such as
points, lines and arcs, as well as operations on compound geometry
elements such as multi-vertex polylines and polygons, and arbitrary
geometry lists.

constants
=========

yapcad.geom provides the "constants" ``epsilon`` and ``pi2`` (2*pi).
Redefine these at your peril.

scalars
=======

Scalar numbers in **yapCAD** are ordinary Python3 ``int`` or ``float``
numbers. This means that in general you will have ordinary
double-precision floating point dynamic range and precision for your
computational geometry operations.  These precision limitations are
reflected in the empirically-chosen value of ``epsilon`` of 5E-6.

vectors
=======

vectors are defined as a list of four numbers, i.e. ``[x,y,z,w]``

When vectors in **yapCAD** are interpreted as three-dimensional
coordinates (as opposed to lists of parameters) the "extra" w
coordinate is a normalization factor. This approach is sometimes
referred to as generalized homogeneous coordinates or projective
coordinates, and was popular with computer graphicists in the
1990s. See https://en.wikipedia.org/wiki/Homogeneous_coordinates

In this interpretation, the w coordinate is a normalization
coordinate, and is either 1.0 or all coordinates are assumed to be
scaled by 1/w.  The inclusion of the w coordinate allows the general
use of transformation matrices for affine transforms.  

There is a yapcad.geom convenience function, ``vect()``, that will make
a vector out of just about any plausible set of arguments.  Unspecified
w values are set to 1, and unspecified z values are set to 0.

points
======

Ordinary **yapCAD** geometry is assumed to lie in the w=1 hyperplane, and
when a vector lies in that plane it is usually interpreted as
transformable coordinate, or point. Because we
expect w=1 unless certain types of affine transforms have been
applied, most **yapCAD** geometry operations do not look at the w
coordinate, and operate as though w=1.  

Any vector that lies in the w>0 half-space is a valid coordinate, or
point. And in **yapCAD**, ``[x,y,z,w]`` corresponds to the same 3D
coordinate as ``[x/w,y/w,z/w,1]``

yapcad.geom provides the ``point()`` convenience function that operates
much like ``vect()``, except that ``point()`` enforces that w is always
greater than zero.

All of the following are examples of a ``yapcad.geom`` point: ::

   pnt1 = point(0,0)
   pnt2 = point(2.0,-2.0,5.0)
   pnt3 = vect(0,0,0,1)
   pnt4 = [1.0, 2.0, 3.0, 1.0]
 
lines
=====

Lines are Python3 lists of two points, e.g. the following are all lines: ::

   line1 = line(point(0,0),point(1,1))
   line2 = line(point(-5,-5,-5),point(5,5,5))
   line3 = [point(x1,y1,z1),point(x2,y2,z2)]

Lines are paramaterized over the interval `0 <= u <= 1` where `u=0`
corresponds to the first point, and `u=1` corresponds to the
second. Points on the interval between 0 and 1 are considered to be
"inside" the line segment, and points outside the interval are still
on the line, but not inside the line segment.  This distinction is
important for intersection calculations, among other operations.

Lines support all standard ``yapcad.geom`` computational geometry
operations. (see below)

*bounding boxes*
----------------

One special type of ``yapcad.geom`` line is a 3D bounding box.  A
bounding box is a line that spans the "lower bottom left" to "upper
top right" of a figure, *e.g.* ``bbx =
[[xmin,ymin,zmin,1],[xmax,ymax,zmax,1]]``. Bounding boxes are used
internally to speed up inside testing and for a variety of other
purposes.

arcs
====

In yapcad.geom, an arc or circle is defined by a center, a radius, a
start angle, an end angle, and a normal that specifies the plane of
the arc/circle.  Here are some examples of ``yapcad.geom`` arcs: ::

  circle1 = arc(point(5.0,5.0),5.0) # full circle implied
  circle2 = arc(point(5.0,5.0),5.0,0,360) # special values for start and end
  notCirc = arc(point(5.0,5.0),5.0,90.0,450.0) # zero-length arc due to mod 360.0
  semicirc1 = arc(point(0,0),10.0,90.0,270.0)
  arc2 = arc(point(-1,1,10),1.0,-45.0,45.0,point(0,0,1),samplereverse=True)

Angles are specified in degrees and are right-handed, which is to say
a positive angle specifies a counter-clockwise sweep.

A full circle is specified through special values of `start` and
`end`; If and only if ``start == 0 and end == 360`` (both integer values)
then the arc is a full circle.  **NOTE:** Other values for start and
end that span a 360 degree difference will actually produce a
zero-length arc, which is probably not what you want.

Arcs are paramaterized over an interval `0 <= u <=1`, where `u=0`
corresponds to the point on the arc at the start angle, and `u=1`
corresponds to the point on the arc at the end angle (unless
``samplereverse==True``, see below).  Parameters less than zero or
greater than one correspond to point on the circle outside the angular
interval between start and end.  **NOTE:** In the special case of a
circle, the point at `u=0` and `u=1` is the same.  Parameters which
correspond to angles greater than 360 degrees or less than zero
degrees "wrap around" to the corresponding angle modulus 360.

An arc is is represented as a list of one or two vectors and a
quasivector, *e.g.* ``[ center, [radius, start, end, -1 ], <normal>]``.
Most of the time, we assume that arcs lie in the x-y plane, which is
to say **yapCAD** assumes ``normal = [0,0,1]`` if it's not
specified. **NOTE:** At present, **yapCAD** doesn't support non-unit-z
normal vectors for arcs, though it will not stop you from specifying
one.

To mark that the second list element is a quasivector, we set the `w`
component to a negative value.  Since negative `w` values should never
exist for points in our projective geometry system this should be a
robust convention.  For ordinary arcs, `w=-1`.  For sample-reversed arcs
(left-handed arcs, where `u=0` begins at the end angle, and `u=1` is
the start angle) `w=-2`

**NOTE:** Sample-reverse arcs are left-handed with respect to the
sampling parameter, not with respect to the specification of angles.
This means that two arcs that are geometrically the same except for
sampling order will have the beneficial property of having the same
`start` and `end` values.

Arcs support all standard ``yapcad.geom`` computational geometry
operations. (see below)

polylines/polygons
==================

A polyline is specified as a list of three or more points which define
a series of continuous line segments.  A polygon is a list of four or
more coplanar points, where the first and last points are coincident,
`i.e.` ``dist(points[0],points[-1]) < epsilon``.  Here are some
examples: ::

  plyline1 =[point(-5,-5),point(5,-5),point(0,5),point(0,10)]
  plygon1 =[point(-5,-5),point(5,-5),point(0,5),point(-5,-5)]

Polylines are parameterized in much the same way as lines.  Values `0
<= u <=1` correspond to points along the polyline or polygon between
the start and end points.  In the case of the non-closed polyline,
values of `u<0` correspond to points on the line in the negative `u`
parameter space of the first line segment, and values of `u>1`
correspond to points on the line in the `u>1` space of the last line
segment.

Polygons are parameterized in the same way, with the exception that,
like circles, the value at `u=0` and `u=1` are the same, and `u`
values outside the `[1,0]` interval "wrap around" (are subject to ``u
% 1.0``) and thus always return a point on the line segments that make
up the polygon.

Polygons with positive area are defined in right-hand
(counterclockwise) order.  A polygon with left-hand point order has
negative area, and implies a hole in a larger surface.

Polylines and polygons support all standard ``yapcad.geom``
computational geometry operations. (see below)

geometry lists
===============

Geometry lists are pretty much exactly what they sound like, which is
to say a Python3 list of other ``yapcad.geom`` elements.  There is
nothing enforcing continuity of elements in a geometry list, so if
continuity, coplanar elements, a closed representation, *etc.* are
required, these constraints must be checked or enforced by the
functions that operate on them. 

Even though there is no inherent continuity constraint, geometry lists
do support element-ordered sampling, unsampling, and intersection, as
well as bounding box calculation and even inside-outside testing,
though the results of such testing will only make sense for geometry
lists composed of closed elements.

Geometry lists are paramaterized a bit like polylines, except that
there is no guarantee of continuity.

======================
COMPUTATIONAL GEOMETRY
======================

The primary purpose of the ``yapcad.geom`` module is to support
computational geometry operations on two-dimensional figures in
three-dimensional space.  

All of the ``yapcad.geom`` geometry representations, or *figures*,
described here support the following operations:

- ``length(x)`` -- return the scalar length of figure ``x``

- ``center(x)`` -- return the center point of figure ``x`` 

- ``bbox(x)`` -- return the three-dimensional bounding box of ``x``

- ``sample(x,u)`` -- parameterizing figure ``x`` over the closed
  interval `[0,1]`, return the point on ``x`` corresponding to the
  parameter ``u``.

- ``unsample(x,p)`` -- given a point ``p`` and figure ``x``,
  parameterizing ``x`` over the closed interval `[0,1]`, return the
  parameter `u` resulting in the point closest to ``p``, or ``False``
  if no point on ``x`` lies within `epsilon` of ``p``.  **FIXME:**
  Unsampling of geometry lists is currently unimplemented.

- ``segment(x,u1,u2)`` -- given a figure ``x`` and two scalar
  parameters ``u1 != u2 and u1 >= 0 and u2 >= 0 and u1 <= 1 and u2 <=
  1``, return a new figure sliced from the original figure spanning
  the specified sampling interval.

- ``isinsideXY(x,p)`` -- for figure ``x`` and point ``p`` that lie in
  the same XY plane, determine whether ``p`` lies within epsilon of
  the interior of figure ``x``.  In the case where ``x`` is a point,
  line, or polyline where there is no two-dimensional interior
  defined, determine if ``p`` lies on or within epsilon of the
  `0 <= u <=1` parameter domain of figure ``x``.

- ``intersectXY(g1,g2,inside=True,params=False)`` -- Compute the
  intersections of the two figures ``g1`` and ``g2`` that lie in the
  same XY plane.  If parameter ``inside == True``, then the
  intersection is only considered valid if it lies within the parameter
  domain `0 <= u <= 1` for each figure.  Return a list of intersection
  points, or the boolean value ``False`` if there are none.  If
  ``params==True``, then return a list of intersection parameter
  values for each figure, or ``False`` if there are no intersections.

In addition to the above, ``yapcad.geom`` supports following affine
transformation operations for all figures:

- ``scale(x,sx=1.0,sy=False,sz=False,cent=point(0,0))`` -- scale the
  figure ``x``, either uniformly or by specified factors in `x`, `y`, and
  `z`.  Returns a new scaled figure.

- ``translate(x,delta)`` -- translate the figure ``x`` by the specified
  ``delta`` vector.  Returns a new translated figure.

- ``rotate(x,ang,cent=point(0,0),axis=point(0,0,1.0))`` -- rotate the
  figure ``x`` by angle ``ang`` degrees about point ``cent`` and axis
  ``axis``.  Returns a new rotated figure.

- ``mirror(x,plane)`` -- mirror the figure ``x`` by the axis-aligned,
  origin-intersecting plane specified by the string ``plane in
  ['xy','yz','xz']``.  **FIXME:** Allow for the specification of an
  arbitrary plane described by an intersection point and a normal, or
  by three points in space. Return a new mirrored figure.

- ``transform(x,m)`` -- apply the arbitrary transformation specified
  by matrix ``m`` (see the documentation for yapcad.xform) to figure
  ``x``.  Return a new transformed figure.

specialized computational geometry operations
===============================================

*operations on points*
----------------------

Points are defined as vectors that lie in a positive, non-zero
hyperplane, i.e. `[x, y, z, w]` such that w > 0.  points are
distinguished from quasivectors, such as the parameters to an
arc, which don't transform as vectors.

Quasivectors, by contrast, lie in a w < 0 hyperplane.  For example,
by convention right-handed arc parameters lie in the quasivector w=-1
hyperplane.

The ``yapcad.geom`` module provides a relatively full set of vector
operations for adding, subtracting, computing inner and outer
products, scaling, *etc.*, including versions that treat
``yapcad.geom`` points/vectors as 3 vectors and versions that treat
them as 4 vectors.

Other notable computational geometry operations on points:

- ``barycentricXY(a,p1,p2,p3):`` -- Given a point ``a`` and three
  vertices ``p1``, ``p2``, and ``p3``, all of which fall into the same
  XY plane, compute the barycentric coordinates lam1, lam2, and lam3
  and return them as a list.
    

*operations on lines*
---------------------

There are several line-specific computational geometry operations, notably:

- ``linePointXY(l,p,inside=True,distance=False,params=False)`` -- For
  a point ``p`` and a line ``l`` that lie in the same XY plane,
  compute the point on ``l`` that is closest to ``p``, and return
  that point. If ``inside`` is true, then return the closest distance
  point between the point and the line segment. If ``distance`` is
  true, return the closest distance, not the point. If ``params`` is
  true, return the sampling parameter value of the closest point.

- ``linePointXYDist(l,p,inside=True)`` -- Convenience function
  wrapping fast point-line distance calculation using
  ``linePointXY()``

*operations on arcs*
--------------------

In addition to the standard operations, there are several arc-specific
computational geometry operations, notably:

- ``lineArcIntersectXY(l,c,inside=True,params=False)`` -- Value safe
  function to compute the intersection of a line ``l`` and an arc
  ``c`` that lie in the same XY plane.  If ``inside == True`` only
  return intersections that lie within the line and arc segment,
  otherwise treat the line as infinite and the arc as a circle.
  Return a list of intersection points, or if ``params == True``
  return a list of intersection parameters instead.

- ``circleCircleTangentsXY(c0,c1)`` -- Value-safe function to compute
  tangent lines to two coplanar circles ``c0`` and ``c1`` lying in the
  same XY plane.  Function will either return two lines or ``False``,
  if the center of the circles are too close together.

*operations on polylines and polygons*
--------------------------------------

In addition to the standard operations, there are several polyline-
and polygon-specific computational geometry operations, notably:

- ``samplepoly(a,u)`` -- Sample the poly ``a`` at parameter ``u``,
  where u is mapped proportionally across the length of all line
  segments, and return the resulting point.  If ``a`` is a closed
  polygon, then sample at ``u % 1.0``, such that samples outside the
  `[0,1]` interval "wrap around".  If ``a`` is not closed and ``u<0``,
  then samples are drawn from the first line in the ``u<0`` paramter
  space of that line.  If ``a`` is not closed and ``u>1``, then
  samples are drawn from the last line in the ``u>1`` paramter space
  of that line.

"""

from math import *
import mpmath as mpm
import copy
import yapcad.xform as xform

## constants
#epsilon=0.0000001
epsilon=0.000005
pi2 = 2.0*pi

## operations on scalars
## -----------------------

    
## utility function to determine if argument is a "real" python
## number, since booleans are considered ints (True=1 and False=0 for
## integer arithmetic) but 1 and 0 are not considered boolean

def isgoodnum(n):
    """ determine if an argument is actually a scalar number, and not boolean
    """
    return (not isinstance(n,bool)) and isinstance(n,(int,float))

## utilty function to determine if scalars a and b are the same to
## within epsilon
def close(a,b):
    """ are two scalars the same within epsilon
    """
    return abs(a-b) < epsilon


## operations on vectors
## ------------------------


def vect(a=False,b=False,c=False,d=False):
    """Convenience function for making a homogeneous coordinates 4 vector
from practically anything
    """
    r = [0,0,0,1]
    if isgoodnum(a):
        r[0]=a
        if isgoodnum(b):
            r[1]=b
            if isgoodnum(c):
                r[2]=c
                if isgoodnum(d):
                    r[3]=d
    elif isinstance(a,(tuple,list)):
        for i in range(min(4,len(a))):
            x=a[i]
            if isgoodnum(x):
                r[i]=x
    return r

## determine if two vectors are the same, to within epsilon
def vclose(a,b):
    return close(mag(sub(a,b)),0)
    
## check to see if argument is a proper vector for our purposes
def isvect(x):
    """
    check to see if argument is a proper vector for our purposes
    """
    return isinstance(x,list) and len(x) == 4 and isgoodnum(x[0]) and isgoodnum(x[1]) and isgoodnum(x[2]) and isgoodnum(x[3])
    
## R^3 -> R^3 functions: ignore w component
## ------------------------------------------------
def add(a,b):
    """ 3 vector, `a + b`"""
    return [a[0]+b[0],a[1]+b[1],a[2]+b[2],1.0]

def sub(a,b):
    """ 3 vector, `a - b`"""
    return [a[0]-b[0],a[1]-b[1],a[2]-b[2],1.0]

def scale3(a,c):
    """ 3 vector, vector ''a'' times scalar ``c``, `a * c`"""
    return [a[0]*c,a[1]*c,a[2]*c,1.0]

## component-wise 3vect multiplication
def mul(a,b):
    """ component-wise 3 vector multiplication"""
    return [a[0]*b[0],a[1]*b[1],a[2]*b[2],1.0]

## NOTE: this function assumes that a lies in the x,y plane.  If this
## is not the case, the results are bogus.
def orthoXY(a):
    """compute an orthogonal vector to vector ``a`` which lies in an XY
plane by crossing `[a1, a2, a3]` with `[0, 0, 1]`"""

    return [ a[1], -a[0], 0, 1.0 ]

## Compute the cross generalized product of a x b, assuming that both
## fall into the w=1 hyperplane
def cross(a,b):
    """Compute the cross generalized product of a x b, assuming that both
    fall into the w=1 hyperplane

    """
    
    return [ a[1]*b[2] - a[2]*b[1],
             a[2]*b[0] - a[0]*b[2],
             a[0]*b[1] - a[1]*b[0],
             1.0 ]

## R^4 -> R^4 functions: operate on w component
def add4(a,b):
    """ 4 vector `a + b`"""
    return [a[0]+b[0],a[1]+b[1],a[2]+b[2],a[3]+b[3]]

def sub4(a,b):
    """ 4 vector `a - b`"""
    return [a[0]-b[0],a[1]-b[1],a[2]-b[2],a[3]-b[3]]

def scale4(a,c):
    """ 4 vector ``a`` times scalar ``c``"""
    return [a[0]*c,a[1]*c,a[2]*c,a[3]*c]

## Homogenize, or project back to the w=1 plane by scaling all values
## by w
def homo(a):
    """Homogenize, or project back to the w=1 plane by scaling all values by w"""
    return [ a[0]/a[3],
             a[1]/a[3],
             a[2]/a[3],
             1 ]

## R^3 -> R functions -- ignore w component
## ----------------------------------------
def dot(a,b):
    """ 3 vector ``a`` dot ``b`` """
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]

def mag(a):
    """ compute the magnitude of 3 vector ``a``"""
    return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2])

def dist(a,b):  # compute distance between two points a & b
    """ compute the euclidean disgtance between two 3 vector points ``a`` and ``b``"""
    return mag(sub(a,b))

## R^3 -> bool functions
## ---------------------

# does point p lie inside 3D bounding box bbox
def isinsidebbox(bbox,p):
    """ does point ``p`` lie inside 3D bounding box ``bbox``?"""
    return p[0] >= bbox[0][0] and p[0] <= bbox[1][0] and\
        p[1] >= bbox[0][1] and p[1] <= bbox[1][1] and\
        p[2] >= bbox[0][2] and p[2] <= bbox[1][2]

# utility function to determine if a list of points lies in the specified
# cardinal plane, one of XY, YZ, XZ
def isCardinalPlanar(plane="xy",points=[]):
    """ utility function to determine if a list of points lies in the specified
# cardinal plane, one of XY, YZ, XZ"""
    if plane=="xy":
        idx = 2
    elif plane == "yz":
        idx = 0
    elif plane == "xz":
        idx = 1
    else:
        raise ValueError('bad plane specification: {}'.format(plane))
        return False

    if len(points) < 2:
        return True

    ref = points[0][idx]
    for i in range(1,len(points)):
        if abs(points[i][idx]-ref) > epsilon:
            return False
    return True

## conveniience functions
def isXYPlanar(points=[]):
    """do all points in list lie in XY plane?"""
    return isCardinalPlanar("xy",points)

def isYZPlanar(points=[]):
    """do all points in list lie in YZ plane?"""
    return isCardinalPlanar("yz",points)

def isXZPlanar(points=[]):
    """do all points in list lie in XZ plane?"""
    return isCardinalPlanar("xz",points)


## R^4 -> R functions 
## ----------------------------------------
def dot4(a,b):
    """ 4 vect dot product"""
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]+a[3]*b[3]

def mag4(a):
    """ 4 vect magnitude calculation"""
    return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]+a[3]*a[3])

def dist4(a,b):  # compute distance between two points a & b
    """ 4 vector distance calculation"""
    return mag4(sub4(a,b))


## misc operations
## -----------------------------------------

## function to deep-copy geometry, which is to say lists containing
## non-zero-length lists or 4-vectors.  Returns False if a isn't valid
## geometry list

# def deepcopy(a):
#     if isvect(a):
#         return vect(a)
#     elif isinstance(a,list):
#         if len(a) > 0:
#             c = list(map(deepcopy,a))
#             for i in c:
#                 if not i:
#                     return False
#             return c
#         else:
#             return False
#     else:
#         return False

deepcopy = copy.deepcopy

# pretty printing string formatter for vectors, lines, and polygons.
# You can use this anywhere you use str(), since it will fall back to
# str() if the argument isn't a yapCAD vector, line, or polygon.
def vstr(a):
    """ utility function for recursively checking and formatting lists
    """
    def _isallvect(foo):
        if not isinstance(foo,list):
            return False
        if len(foo)==1:
            return isinstance(foo[0],list) and \
                (isvect(foo[0]) or _isallvect(foo[0]))
        else:
            return (isvect(foo[0]) or _isallvect(foo[0])) and \
                _isallvect(foo[1:])
    def _makestr(foo):
        if len(foo) ==1:
            return vstr(foo[0])
        else:
            return vstr(foo[0]) + ", " + _makestr(foo[1:])
    # if it's not a list, it's not a vector, line, or poly    
    if not isinstance(a,list):
        return str(a) # fall back to default string formatting
    # if it is a vector, format it to leave out extraneous
    # coordinates.  NOTE: this means that 3 vectors that happen to
    # fall into the z=0 plane will be formatted as though they were 2
    # vectors.
    if isvect(a):
        if abs(a[3]-1.0) > epsilon: # not in w=1
            return "[{}, {}, {}, {}]".format(a[0],a[1],a[2],a[3])
        elif abs(a[2]) > epsilon: # not in z=0
            return "[{}, {}, {}]".format(a[0],a[1],a[2])
        else: # in x-y plane
            return "[{}, {}]".format(a[0],a[1])
    # if it is a list of vectors (line or poly) format it appropriately
    elif len(a)>0 and _isallvect(a):
        return "["+_makestr(a)+"]"
    else: #it's none of those things, fall back to default string
          #formatting
        return str(a)


## COMPUTATIONAL GEOMETRY
## ======================
## operations on points
## --------------------

## points are defined as vectors that lie in a positive, non-zero
## hyperplane, i.e. [x, y, z, w] such that w > 0.  points are
## distinguished from pesuedovectors, such as the parameters to an
## arc, which don't transform as vectors.

## quasivectors, by contrast, should lie in a w < 0 hyperplane.  For
## example, by convention arc parameters lie in the quasivector w=-1
## hyperplane

## since points are vectors, we don't have a lot of special operations
## on points, except to explicitly create them and to test for them.

def point(x=False,y=False,z=False,w=False):
    """Point creation from point or scalars"""
    if ispoint(x):
        return deepcopy(x)
    r = [0,0,0,1]
    if isgoodnum(x):
        r[0]=x
        if isgoodnum(y):
            r[1]=y
            if isgoodnum(z):
                r[2]=z
                if isgoodnum(w):
                    r[3]=w
    if r[3] > 0:
        return r
    else:
        raise ValueError('bad w argument to point()')


def ispoint(x):
    """ is it a point?"""
    if isvect(x) and x[3] > 0.0:
        return True
    return False

## the length of a point is by definition zero
def pointlength(x):
    """ an odd, completionist function that always returns zero"""
    return 0.0

## the center of a point is.... the same point
def pointcenter(x):
    """ the center of a point is.... the same point"""
    return point(x)

## this only makes sense in the context of generaling computational
## geometry operations
def samplepoint(x,u):
    """whatever parameter you pass to sampling a point, you get the same point"""
    return point(x)

## compute 3D bounding box of point, which is 2 epsilon on a side.  Note,
## this guarantees that two points that are wtihin epsilon of
## eachother will be within eiach other's bbox
def pointbbox(x):
    """compute 3D bounding box of point, which is 2 epsilon on a side.  Note,
    this guarantees that two points that are wtihin epsilon of
    eachother will be within each other's bbox
    """
    ee = point(epsilon,epsilon,epsilon)
    return [sub(x,ee),add(x,ee)]

## inside testing for points.  Only true if points are the same
## within epsilon
def isinsidepointXY(x,p):
    """ inside testing for points.  Only true if points are the same
    within epsilon"""
    
    return dist(x,p) < epsilon

## operations on lines
## --------------------

## lines are defined as lists of two points, i.e.  [point(x1,
## y1),point(x2, y2)].  Lines must lie in a positive hyperplane,
## viz. w>0

## make a line, copying points, value-safe
def line(p1,p2=False):
    """Value-safe line creation"""
    if isline(p1):
        return deepcopy(p1)
    elif ispoint(p1) and ispoint(p2):
        return [ point(p1), point(p2) ]
    else:
        raise ValueError('bad values passed to line()')

## is it a line?
def isline(l):
    """ is it a line? """
    return isinstance(l,list) and len(l) == 2 \
        and ispoint(l[0]) and ispoint(l[1])

## return the length of a line
def linelength(l):
    """ return the length of a line"""
    return dist(l[0],l[1])

## return the center of a line
def linecenter(l):
    """ return the center of a line"""
    return scale3(add(l[0],l[1]),0.5)

## Sample a parameterized line.  Values 0 <= u <= 1.0 will fall within
## the line segment, values u < 0 and u > 1 will fall outside the line
## segment.

def sampleline(l,u):
    """Sample a parameterized line ``l``.  Values `0 <= u <= 1.0` will
    fall within the line segment, values `u < 0` and `u > 1` will fall
    outside the line segment.

    """

    p1=l[0]
    p2=l[1]
    p = 1.0-u
    return add(scale3(p1,p),scale3(p2,u))

def segmentline(l,u1,u2):
    """slice out a parameterized segment from a line and return this as a new line segment"""
    p1=sampleline(l,u1)
    p2=sampleline(l,u2)
    return [p1,p2]

## function to "unsample" a line -- given a point on a line, provide
## the corresponding parametric value.  Return False if the distance
## of the point from the line is greater than epsilon

def unsampleline(l,p):
    """
    function to "unsample" a line -- given a point on a line, provide
    the corresponding parametric value.  Return False if the distance
    of the point from the line is greater than epsilon
    """

    v1 = sub(l[1],l[0])
    v2 = sub(p,l[0])
    z = cross(v1,v2)            # magnitude of cross product is zero
                                # if vectors parallel
    if mag(z) > epsilon:        #  not parallel case
        return False
    len1 = mag(v1)
    len2 = mag(v2)
    if dot(v1,v2) > 0:
        return len2/len1
    else:
        return -len2/len1
        
## compute 3D line bounding box
def linebbox(l):
    """ compute 3D line bounding box"""
    p1=l[0]
    p2=l[1]
    return [ point(min(p1[0],p2[0]),min(p1[1],p2[1]),min(p1[2],p2[2])),
             point(max(p1[0],p2[0]),max(p1[1],p2[1]),max(p1[2],p2[2])) ]

## inside testing for line -- only true of point lies on line, to
## within epsilon

def isinsidelineXY(l,p):
    """
    inside testing for line -- only true of point lies on line, to within epsilon.
    """

    return linePointXY(l,p,distance=True) < epsilon

## Compute the intersection of two lines that lie in the same x,y plane
def lineLineIntersectXY(l1,l2,inside=True,params=False):
    """Compute the intersection of two lines that lie in the same XY
    plane. **NOTE:** It's usually preferable to use the generalized
    ``intersectXY()`` function than this non-type-safe, line-specific one.
    """

    x1=l1[0][0]
    y1=l1[0][1]
    z1=l1[0][2]
    
    x2=l1[1][0]
    y2=l1[1][1]
    z2=l1[1][2]

    x3=l2[0][0]
    y3=l2[0][1]
    z3=l2[0][2]
    
    x4=l2[1][0]
    y4=l2[1][1]
    z4=l2[1][2]

    ## check for x,y planar consistency
    if abs(z2-z1) > epsilon or abs(z3-z1) > epsilon or abs(z4-z1) > epsilon:
        raise ValueError('lines not in same x-y plane')

    ## do lines intersect anywhere?
    denom=(x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)
    if denom*denom < epsilon:
        return False

    ## the lines do intersect, so let's see if they intersect
    ## inside both line segments
    t = ((x1-x3)*(y3-y4) - (y1-y3)*(x3-x4))/denom
    u = -1 * ((x1-x2)*(y1-y3) - (y1-y2)*(x1-x3))/denom

    ## return the paramater space intersection
    if params:
        return [t,u]
    
    ## do we care about falling inside the line segments? if so,
    ## check that the intersection falls within
    if inside and ( t < 0.0 or t > 1.0 or u < 0.0 or u > 1.0):
        return False

    return [x1 + t*(x2-x1), y1+t*(y2-y1), z1, 1.0]

## for a point and a line that lie in the x,y plane, compute the
## closest distance point on the line to the point, and return that
## point. If inside is true, then return the closest distance point
## between the point and the line segment.  If distance is true,
## return the distance, not the point.
##
## NOTE: in the case where inside is False, it's faster to compute the
## distance with linepoint(... inside=False, distance=True) than it is
## to compute the intersection point, so that is prefereable to the
## two-step operation of calling linePoint() and then dist() on the
## result.

def linePointXY(l,p,inside=True,distance=False,params=False):

    """
    For a point ``p`` and a line ``l`` that lie in the same XY plane,
    compute the point on ``l`` that is closest to ``p``, and return
    that point. If ``inside`` is true, then return the closest distance
    point between the point and the line segment. If ``distance`` is
    true, return the closest distance, not the point. If ``params`` is
    true, return the sampling parameter value of the closest point.

    """
    a=l[0]
    b=l[1]
    # check for degenerate case of zero-length line
    abdist = dist(a,b)
    if abdist < epsilon:
        #raise ValueError('zero-length line passed to linePointXY')
        print('zero-length line passed to linePointXY')
        return False

    if distance and params:
        raise ValueError('incompatible distance and params parameters passed to linePointXY')

    x0=p[0]
    y0=p[1]
    z0=p[2]
    x1=a[0]
    y1=a[1]
    z1=a[2]
    x2=b[0]
    y2=b[1]
    z2=b[2]

    ## check to see if all three points lie in the same x,y plane
    if not isXYPlanar([p,a,b]):
        raise ValueError('non-XY points in linePointXY call')
        return false
    # if abs(z1-z0) > epsilon or abs(z2-z0) > epsilon:
    #    return False

    linedist = abs( ((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1)/abdist)

    ## this is the fast case:
    if not inside and distance:
        return linedist
    
    ## find out where the intersection between the original line and a
    ## line defined by the point and an orthogonal direction vector
    ## is.  We do this by constructing two direction vectors
    ## orthogonal to the orgiginal line scaled by the line distance,
    ## and adding them to the point in question.  Assuming that the
    ## line distance is not zero, only one of these constructed points
    ## will fall on the line

    ## compute unit direction vector for original line
    dir = sub(b,a)
    dir = scale3(dir,1.0/mag(dir))

    ## compute two orthogonal direction vectors of length linedist
    ordir1 = scale3(orthoXY(dir),linedist)
    ordir2 = scale3(ordir1, -1.0)
    
    ## there are two possible intersection points
    pi1 = add(p,ordir1)
    pi2 = add(p,ordir2)

    ## compute distances
    d1pa = dist(a,pi1)
    d1pb = dist(pi1,b)
    d1 =  d1pa+d1pb # "triangle" with pi1

    d2pa = dist(a,pi2)
    d2pb = dist(pi2,b)
    d2 =  d2pa+d2pb # "triangle" with pi2

    ## the shortest "triangle" distance will signal the point that
    ## is actually on the line, even if that point falls outside
    ## the a,b line interval
    
    if params or not inside: # if we don't care about being inside the
                             # line segment
        if d1 <= d2:
            if distance:
                return d1
            elif params:
                return d1pb/abdist
            else:
                return pi1
        else:
            if distance:
                return d2
            elif params:
                return d2pb/abdist
            else:
                return pi2
        
    
    ## if the closest point on the line to point p lies between
    ## the endpoints of the line, then either d1 or d2 will equal
    ## abdist.  IF neither do, then we know that the closest point lies
    ## outside the endpoints

    if abs(d1-abdist) < epsilon:
        if distance:
            return linedist
        else:
            return pi1

    if  abs(d2-abdist) < epsilon:
        if distance:
            return linedist
        else:
            return pi2

    ## closest point is outside the interval.  That means that the
    ## distance from point p to whichever endpoint is smaller is the
    ## closest distance

    d3 = dist(a,p)
    d4 = dist(b,p)

    if d3 < d4:
        if distance:
            return d3
        else:
            return a
    else:
        if distance:
            return d4
        else:
            return b

## convenience function for fast distance calc
def linePointXYDist(l,p,inside=True):
    """
    Convenience function wrapping point-line distance calculation using ``linePointXY()``
    """
    return linePointXY(l,p,inside,distance=True)

## functions that operate on 2D arcs/circles
##------------------------------------------------------------

## an arc or circle is defined by a center, a radius, a start angle,
## an end angle, and a normal that specifies the plane of the
## arc/cicle.

## angles are in degress, and specify a counter-clockwise sweep.  If
## start = 0 and end = 360 (both integer) the arc is a circle.  Other
## values that create a 360 degree difference may produce an
## infentesimal gap.

## an arc is defined as a list of one or two vectors and a
## quasivector: [ center, [r, s, e], <normal>].  Most of the time, we
## assume that arcs lie in the x-y plane, which is to say that
## normal = [0,0,1] if it's not specified.

## to mark that the second list element is a quasivector, we set the
## w component to a negative value.  Since negative w values should
## never exist for vectors in our projective geometry system this
## should be a robust convention.  For ordinary arcs, w=-1.  For
## sample-reversed arcs where u=0.0 begins at the end angle, and u=1.0
## is the start agle, w=-2

## make an arc, copying points, value-safe
## NOTE: if start and end are not specified, a full circle is created
def arc(c,rp=False,sn=False,e=False,n=False,samplereverse=False):
    """
    Construct an arc by copying an eisting arc or specifying a center ``c`` and various optional parameters.

    """
    if isarc(c):
        return deepcopy(c)
    elif ispoint(c):
        cen = point(c)
        w=-1
        if samplereverse:
            w=-2
        if isvect(rp):
            psu = deepcopy(rp)
            r=psu[0]
            if r < 0:
                raise ValueError('negative radius not allowed for arc')
            psu[3]=w
            if not sn:
                return [ cen, psu ]
            else:
                if not ispoint(sn) or abs(mag(sn)-1.0) > epsilon:
                    raise ValueError('bad (non-unitary) plane vector for arc')
                return [ cen, psu, point(sn) ]
        elif isgoodnum(rp):
            r = rp
            if r < 0:
                raise ValueError('negative radius not allowed for arc')
            start = 0
            end = 360
            if isgoodnum(sn) and isgoodnum(e):
                start = sn
                end = e
            psu = vect(r,start,end,w)
            if not n:
                return [ cen, psu ]
            else:
                if not ispoint(n) or abs(mag(n)-1.0) > epsilon:
                    raise ValueError('bad (non-unitary) plane vector for arc')
                return [ cen, psu, point(n) ]

    raise ValueError('bad arguments passed to arc()')
        
## function to determine if argument can be interpreted as a valid
## arc.

def isarc(a):
    """ is it an arc? """
    if not isinstance(a,list):
        return False
    n = len(a)
    if n < 2 or n > 3:
        return False
    if not (ispoint(a[0]) and isvect(a[1])):
        return False
    if a[1][3] not in (-1,-2):           # is psuedovector marked?
        return False
    r =a[1][0]                
    if r < 0:
        # is radius less than zero? if so, not a valid arc
        return False
    if n == 3 and ( not ispoint(a[2]) or abs(mag(a[2])-1.0) > epsilon):
        # if plane-definition vector is included but is non-unitary,
        # it's not a valid arc
        return False
    
    return True

## test to see if an arc is a circle.  NOTE that due to modulus
## arithmetic, it's not possible to specify an arc with a full 360
## degees range, as this would "wrap around" to an arc of 0 degrees
## range.  To solve this, we use a special integer convention for
## start and end to signal a true full circle

def iscircle(a):
    """ is it a circle? """
    if isarc(a):
        start=a[1][1] 
        end=a[1][2]
        ## these are special, integer values that flag a true full
        ## circle.
        if start==0 and end==360:
            return True
    else:
        return False

## function to return the midpoint of an arc
def arccenter(c):
    """ 
    return the midpoint of an arc.  In the special case of a complete circle, 
    return the center of the circle.
    """
    start=c[1][1] 
    end=c[1][2]
    if start == 0 and end == 360:
        return c[0]
    else:
        return samplearc(c,0.5)
    
## function to return the length of an arc
def arclength(c):
    """return scalar length of an arc"""
    r=c[1][0]
    start=c[1][1]
    end=c[1][2]
    if start == 0 and end == 360:
        return r * pi2
    else:
        start = start % 360.0
        end = end % 360.0

    if start > end:
        end = end+360.0
    d = pi2*r
    l = d*(end-start)/360.0
    return l

## Sample the arc over the start-end angle interval
def samplearc(c,u,polar=False):
    """sample the arc ``c`` at parameter ``u`` and return the resulting point.  If ``polar`` is true, return polar coordinates of angle and radius with cartesian center."""
    p=c[0]
    r=c[1][0]
    start=c[1][1]
    end=c[1][2]
    w=c[1][3]
    if w == -2: # if arc is flagged as sample-reversed, flip u
        u=1.0-u
    if start != 0 and end != 360:
        start = start % 360.0
        end = end % 360.0
        if end < start:
            end += 360.0
    if len(c) == 3:
        norm = c[2]
        if dist(norm,vect(0,0,1)) > epsilon:
            raise NotImplementedError('non x-y plane arc sampling not yet supported')
    angle = ((end-start)*u+start)%360.0
    radians = angle*pi2/360.0
    if polar: # return polar coordinates with cartesian center
        return [ p,r,angle]
    
    q = scale3(vect(cos(radians),sin(radians)),r)
    return add(p,q)

# Given an arc paramaterized on a 0,1 interval, return a new arc with
# the same center and radius spanning the interval u1,u2

def segmentarc(c,u1,u2):
    """
    Given an arc paramaterized on a 0,1 interval, return a new arc with
    the same center and radius spanning the interval u1,u2
    """

    pol1=samplearc(c,u1,polar=True)
    pol2=samplearc(c,u2,polar=True)
    sr= (c[1][3] == -2)
    if sr:
        return arc(pol1[0],pol1[1],pol2[2],pol1[2],samplereverse=True)
    else:
        return arc(pol1[0],pol1[1],pol1[2],pol2[2])

## unsample the arc, which is to day given a point that is on the
## circle, return it's corresponding sample parameter, or False if the
## point doesn't lie on the circle
def unsamplearc(c,p):
    """
    unsample the arc ``c`` at point ``p``, which is to day given a point that is on the
    circle, return it's corresponding sample parameter, or False if the
    point is more than `epsilon` from the circle.
    """
    x = sub(p,c[0])
    r=c[1][0]
    start=c[1][1]
    end=c[1][2]
    # if close(end-start,360.0): #rotated circles may look like this
    #     start = 0
    #     end = 360
        
    if start != 0 and end != 360:
        start = start % 360.0
        end = end % 360.0
        if end < start:
            end += 360.0    
    if close(start,end):
        # degenerate, zero-length arc
        return False
    if len(c) == 3:
        norm = c[2]
        if dist(norm,vect(0,0,1)) > epsilon:
            raise NotImplementedError('non x-y plane arc unsampling not yet supported')
    if abs(mag(x)-r) > epsilon:
        return False            # point is not on arc circle
    ang = (atan2(x[1],x[0]) % pi2)*360/pi2
    if end > 360.0 and ang <= end-360.0:
        ang = ang + 360.0
    u = (ang-start)/(end-start)
    if c[1][3] == -2: #samplereverse
        u = 1.0-u
    return u
    

def arcbbox(c):
    """return bounding box for arc"""
    if iscircle(c):
        rr=point(c[1][0],c[1][0])
        return [sub(c[0],rr),add(c[0],rr)]
    else:
        pp = []
        for i in range(5):
            u = i/4
            pp.append(samplearc(c,u))
        return polybbox(pp)

def isinsidearcXY(c,p):
    """
    Determine if point ``p`` lies on (within `epsilon of) the arc ``c``.  In the special
    case where ``c`` is a circle, determine if point ``p`` lies
    anywhere on or inside the circle.

    """

    x = c[0]
    r = c[1][0]
    if dist(x,p) > r:
        return False
    if iscircle(c):
        return True
    start = c[1][1]%360.0
    end = c[1][2]%360.0
    if end < start:
        end+= 360.0
    p2 = sub(p,x)
    ang = (atan2(p2[1],p2[0]) % pi2)*360/pi2

    if end <= 360.0:
        return (ang >= start and ang <= end)
    else:
        return ang >= start or ang <= (end-360.0)
        

## Intersection functions for arcs and circles
## we can compute the intersection of an arc and a line, or an arc and
## an arc, when these all lie in the same plane.

## arc-arc intersection calculation, non-value-safe version
def _arcArcIntersectXY(c1,c2,inside=True,params=False):
    """non-value-safe function to compute the intersetion of two arcs"""
    x1=c1[0]
    x2=c2[0]
    r1=c1[1][0]
    r2=c2[1][0]

    # check for sample reverse condition
    sr1 = c1[1][3]==-2
    sr2 = c2[1][3]==-2

    ## first check for non-intersection due to distance between the
    ## centers of the arcs, treating both arcs as circles for the moment

    d=dist(x1,x2) #calculate the distance d between circle centers

    if d > r1+r2:
        return False # too far, no possiblity of intersection

    if ( r1> r2 and d < r1-r2) or (r2 >= r1 and d < r2-r1):
        return False # too close, little arc is fully inside bigger arc

    if d < epsilon:
        return False # circle centers too close for stable calculation

    ## OK, we are in the goldilocks zone of intersection.  this means
    ## that if boh arcs are cicles or if inside=False we are
    ## guaranteed one or two intersections.  Calculate those
    ## intersections and then test to see if they fall between start
    ## and end of the respective arcs

    ## we start by calculating the distance id of the intersection plane
    ## from the center of arc 1, knowing that by definition id <= r1

    ## Math: consider the triangle with side lengths r1, r2, and d,
    ## where d is the previously calculated distance between arc
    ## centers.  Consider the two right triangles with side lengths
    ## r1, id, h, and r2, h, (d-id).  We know that:
    ## id^2 + h^2 = r1^2, (d-id)^2 + h^2 = r2^2
    ## solving both for h2 and then substituting, this means:
    ## r1^2 - id^2 = r2^2 - (d-id)^2
    ## collecting terms and solving for id produces:
    ## id = (r1^2-r2^2 + d^2)/2d

    id = (r1*r1 - r2*r2 + d*d)/(2 * d)

    ## compute the point on the line connecting the two arc centers
    ## that is id away from the first arc

    v1 = scale3(sub(x2,x1),1.0/d) # unitary direction vector pointing
                                 # from x1 to x2
    v2 = scale3(v1,id) # point on line between two circles in
                      # coordinate space centered at x1

    ## compute direction vector o orthgonal to v1 -- the line that
    ## intersects point v2 and v2+o will pass through our intersection
    ## points

    o = orthoXY(v1)
    
    ## now, transform back into world coordinates and calculate the
    ## intersection of this line with either of our arcs, treating
    ## them as circles for now

    l = [add(v2,x1),add(add(v2,o),x1)]

    s = _lineArcIntersectXY(l,c1,False)

    ## as a sanity check, do the same with the other arc.  Results
    ## should be within epsilon
    #ss = _lineArcIntersectXY(l,c2,False)
    #foo = list(map(lambda x, y: dist(x,y) < epsilon,s,ss))
    #print("sanity check: " , foo)

    if not s or len(s) == 0:
        raise ValueError('no computed intersections, something is wrong')

    if not inside and not params:
        return s
    
    ## jump back to arc1 and arc2 space and check angles

    s1 = list(map(lambda x: sub(x,x1),s))
    s2 = list(map(lambda x: sub(x,x2),s))

    ## compute start and end angles for arcs
    start1=c1[1][1]
    end1=c1[1][2]
    if not (start1 == 0 and end1 == 360):
        start1 = start1 % 360.0
        end1 = end1 % 360.0
        if end1 < start1:
            end1 = end1 + 360.0
        
    start2=c2[1][1]
    end2=c2[1][2]
    
    if not (start2 == 0 and end2 == 360):
        start2 = start2 % 360.0
        end2 = end2 % 360.0
        if end2 < start2:
            end2 = end2 + 360.0
        

    ## check each intersection against angles for each arc.  
    ss = []
    uparam1 = []
    uparam2 = []
    for i in range(len(s)):
        p1 =s1[i]
        p2 =s2[i]
        ang1 = (atan2(p1[1],p1[0]) % pi2)*360.0/pi2
        ang2 = (atan2(p2[1],p2[0]) % pi2)*360.0/pi2

        if params:
            u1 = 0
            u2 = 0
            if end1 <= 360.0 or ang1 >= start1 or \
               ( end1 > 360.0 and ang1 > end1-360.0):
                u1 = (ang1-start1)/(end1-start1)
                if sr1:
                    u1 = 1.0-u1
            elif end1 > 360.0:
                u1 = (ang1+360.0-start1)/(end1-start1)
                if sr1:
                    u1 = 1.0-u1
            uparam1 = uparam1 + [ u1 ]
                
            if end2 <= 360.0 or ang2 >= start2 or \
               ( end2 > 360.0 and ang2 > end1-360.0):
                u2 = (ang2-start2)/(end2-start2)
                if sr2:
                    u2 = 1.0-u2
            elif end2 > 360.0:
                u2 = (ang2+360.0-start2)/(end2-start2)
                if sr2:
                    u2 = 1.0-u2
            uparam2 = uparam2 + [ u2]
                
        else:
            good = False
            ## check angle against first arc
            if end1 <= 360.0 and ang1 >= start1 and ang1 <= end1:
                good = True
            elif end1 > 360.0 and (ang1 >= start1 or ang1<= end1-360.0):
                good = True

            ## check angle against second arc
            if end2 <= 360.0 and  ang2 >= start2 and ang2 <= end2:
                good = good and True
            elif end2 > 360.0 and (ang2 >= start2 or ang2<= end2-360.0):
                good = good and True
            else:
                good = False

            ## only add instersection to the list if both checks were passed
            if good:
                ss = ss + [ s[i] ]
                
    if not params and len(ss) == 0:
        return False
    else:
        if params:
            return [uparam1,uparam2]
        else:
            return ss

## value-safe wrapper for arc-arc intersection function
def arcArcIntersectXY(c1,c2,inside=True,params=False):
    """
    Value-safe function to compute the intersection of two arcs that
    lie in the same XY plane.  If ``inside == True``, only count
    intersections within the arc segments.  If ``inside == False``,
    then treat the arcs as circles for intersection
    calculation. Return a list of the intersection points, or if
    ``params==True`` return a list of the intersection parameters for
    each arc.

    """
    
    for c in [c1,c2]:
        if len(c) == 3:
            norm = c[2]
            if dist(norm,vect(0,0,1)) > epsilon:
                raise ValueError('arc passed to lineArcIntersectXY does not lie in x-y plane')
    if not isXYPlanar([c1[0],c2[0]]):
        raise ValueError('arcs passed to arcArcIntersectXY do not lie in same x-y plane')
    return _arcArcIntersectXY(c1,c2,inside,params)
    

## non-value-safe line-arc intersection function
def _lineArcIntersectXY(l,c,inside=True,params=False):
    """
    non-value-safe function for computing the intersection of a line and an arc
    """
    x=c[0]
    r=c[1][0]
    mpr=mpm.mpf(r)
    
    # is the arc a full circle?
    circle = False
    if c[1][1] == 0 and c[1][2] == 360:
        circle = True
        
    start=c[1][1] % 360.0
    end=c[1][2] %360.0

    ## what is the shortest distance between the line and the center
    ## of the arc?  If that is greater than r, then there is no
    ## intersection
    dst = linePointXYDist(l,x,inside and not params)
    if dst > r+epsilon:
        return False

    ## start by treating the arc as a circle.  At this point we know
    ## we have one or two intersections within the line segment,
    ## though perhaps none within the arc segment, which we will test
    ## for later
    
    ## transform points so arc is located at the origin
    p0=sub(l[0],x)
    p1=sub(l[1],x)
    
    ## solve for b in:  | b*p0 + (1-b)*p1 | = r
    ## let V= p0-p1, P=p1
    ##     | b*V + P |^2 = r^2
    ##       b^2(Vx^2 + Vy^2) + 2b(VxPx+VyPy) + Px^2 + Py^2 - r^2 = 0
    ## let a = Vx^2 + Vy^2,
    ##     b = 2*(VxPx + VyPy)
    ##     cc = Px^2 + Py^2 - r^2
    ## b0 = ( -b + sqrt(b^2 - 4ac) )/ 2a
    ## b1 = ( -b - sqrt(b^2 - 4ac) )/ 2a
    
    V = sub(p0,p1)
    P = p1
    #a = V[0]*V[0]+V[1]*V[1]
    mpV0 = mpm.mpf(V[0])
    mpV1 = mpm.mpf(V[1])
    mpP0 = mpm.mpf(P[0])
    mpP1 = mpm.mpf(P[1])
    a = mpV0*mpV0+mpV1*mpV1
    mpepsilon = mpm.mpf(epsilon)
    if mpm.fabs(a) < mpepsilon*mpepsilon:
        print('degenerate line in lineArcIntersectXY')
        raise ValueError('bad!')
        return False
    # b = 2*(V[0]*P[0]+V[1]*P[1])
    b = 2*(mpV0*mpP0+mpV1*mpP1)
    #cc = P[0]*P[0]+P[1]*P[1]-r*r
    cc = mpP0*mpP0+mpP1*mpP1-mpr*mpr
    d = b*b-4*a*cc
    ## Check to see if we are within epsilon, scaled by the length of the line
    if mpm.fabs(d) < mpm.sqrt(a)*2*mpepsilon: # one point of intersection
        b0 = -b/(2*a)
        b1 = False
    elif d < 0:
        print("value of d: ",d,"  value of sqrt(a)*epsilon",sqrt(a)*epsilon)
        raise ValueError("imaginary solution to circle line intersection -- shouldn't happen here")
    else: # two points of intersection
        b0 = (-b + mpm.sqrt(d))/(2*a)
        b1 = (-b - mpm.sqrt(d))/(2*a)

    # use computed parameters to calculate solutions, still in
    # circle-at-origin coordinates
    s = [ add(scale3(V,float(b0)),p1) ]
    if b1:
        s = s + [ add(scale3(V,float(b1)),p1) ]

    if not inside or circle or params:              # transform back into world
                                          # coordinates
        pp = list(map(lambda q: add(q,x),s))
        if params:
            uu1 = []
            uu2 = []
            for i in range(len(pp)):
                uu1 = uu1 + [ unsampleline(l,pp[i]) ]
                uu2 = uu2 + [ unsamplearc(c,pp[i]) ]
            return [uu1, uu2]
        else:
            return pp

    ## see if any of the intersections we've found lie between
    ## start and end of the arc
    
    ss = []
    for i in s:
        ang = (atan2(i[1],i[0]) % pi2)*360.0/pi2

        if end > start and ang >= start and ang <= end:
            ss = ss + [ add(x,i) ]
        elif end < start and (ang >= start or ang<= end):
            ss = ss + [ add(x,i) ]

    if len(ss) == 0:
        return False
    return ss

## value-safe wrapper for line-arc intersection function
def lineArcIntersectXY(l,c,inside=True,params=False):
    """
    Value safe function to compute the intersection of a line ``l``
    and an arc ``c`` that lie in the same XY plane.  If ``inside ==
    True`` only return intersections that lie within the line and arc
    segment, otherwise treat the line as infinite and the arc as a
    circle.  Return a list of intersection points, or if ``params ==
    True`` return a list of intersection parameters instead.

    """
    
    if len(c) == 3:
        norm = c[2]
        if dist(norm,vect(0,0,1)) > epsilon:
            raise ValueError('arc passed to lineArcIntersectXY does not lie in x-y plane')
    points = l + [ c[0] ]
    if not isXYPlanar(points):
        raise ValueError('line and circle passed to lineArcIntersectXY do not all lie in same x-y plane')
    return _lineArcIntersectXY(l,c,inside,params)


## function to compute tangent lines to two coplanar circles lying in
## an x-y plane.  Function will either return two lines or False, if
## the center of the circles are too close together.

def _circleCircleTangentsXY(c1,c2):
    """ non-value-safe function to compute tangent lines to two coplanar circles lying in
    an x-y plane.  Function will either return two lines or False, if
    the center of the circles are too close together.
    """

    a = c1[1][0]
    b = c2[1][0]
    if a>b:
        bigIsOne=True
        bigC = c1
        smallC = c2
    else:
        bigIsOne=False
        bigC = c2
        smallC = c1
    ## Consdier the triangle created by the center of the small
    ## circle, the center of the large circle, and the point at the 90
    ## degree intersection of the line from the center of the small
    ## circle to the radian of the tangent point on the large circle.
    ## This is a right triangle with one leg of length d (distance of
    ## centers), one leg of length bigR-smallR, and one leg of unknown
    ## length, beta. theta is the angle formed by d and beta, which is
    ## also the angle of one of the the tangent lines, the other being
    ## -theta.
    ## 
    ## we will calulate theta as follows:
    ## beta^2 - (r2-r1)^2 = d^2
    ## beta = sqrt( d^2 - (r2-r1)^2 )
    ## theta = atan ((r2-r1)/beta)
    
    r1 = smallC[1][0]
    r2 = bigC[1][0]

    d = dist(c1[0],c2[0])
    mpd = mpm.mpf(d)
    dr = r2-r1
    mpdr = mpm.mpf(dr)

    if d <= dr: #centers too close
        raise ValueError('circleCircleTangentsXY: centers of circles too close')
    
    beta = mpm.sqrt( mpd*mpd - mpdr*mpdr)
    theta = float(mpm.atan2(dr,beta))

    ## now, figure out the angle created by the center of the large
    ## circle with respect to the small circle
    dd = sub(bigC[0],smallC[0])
    phi = atan2(dd[1],dd[0])

    ## the two lines have angle phi+theta, and phi-theta.  The
    ## intersection point of these lines is at the point on the circle
    ## phi+theta+90', and phi-theta-90'
    gamma1 = phi+theta+pi/2
    gamma2 = phi-theta-pi/2
    n1 = point(cos(gamma1),sin(gamma1))
    n2 = point(cos(gamma2),sin(gamma2))
    p1 = add(scale3(n1,r1),smallC[0])
    p2 = add(scale3(n1,r2),bigC[0])
    p3 = add(scale3(n2,r1),smallC[0])
    p4 = add(scale3(n2,r2),bigC[0])

    l1 = l2 = []
    if bigIsOne:
        l1=line(p2,p1)
        l2=line(p4,p3)
    else:
        l1 = line(p1,p2)
        l2 = line(p3,p4)

    return [l1,l2]

## value safe wrapper
def circleCircleTangentsXY(c0,c1):
    """
    Value-safe function to compute tangent lines to two coplanar circles ``c0`` and ``c1`` lying in
    the same XY plane.  Function will either return two lines or ``False``, if
    the center of the circles are too close together.
    """
    
    if not iscircle(c0) or not iscircle(c1):
        raise ValueError('non circles passed to circleCirlceTangentsXY')
    if len(c0) == 3:
        norm = c0[2]
        if dist(norm,vect(0,0,1)) > epsilon:
            raise ValueError('first circle passed to circleCircleTangentsXY does not lie in an x-y plane')
    if len(c1) == 3:
        norm = c1[2]
        if dist(norm,vect(0,0,1)) > epsilon:
            raise ValueError('second circle passed to circleCircleTangentsXY does not lie in an x-y plane')
    points = [ c0[0], c1[0]]
    if not isXYPlanar(points):
        raise ValueError('circles passed to circleCircleTangentsXY do not all lie in same x-y plane')
    return _circleCircleTangentsXY(c0,c1)

def pointCircleTangentsXY(p,c1):
    c0 = arc(p,epsilon)
    return circleCircleTangentsXY(c0,c1)

## functions that compute barycentric coordinates for 2D triangles
## -----------------------------------------------------------

## given a point a and three vertices, all of which fall into the x-y
## plane, compute the barycentric coordinates lam1, lam2, and lam3
def _barycentricXY(a,p1,p2,p3):
    """
    given a point a and three vertices, all of which fall into the XY
    plane, compute the barycentric coordinates lam1, lam2, and lam3
    """
    x=a[0]
    y=a[1]

    x1=p1[0]
    y1=p1[1]

    x2=p2[0]
    y2=p2[1]

    x3=p3[0]
    y3=p3[1]

    denom=(y2-y3)*(x1-x3) + (x3-x2)*(y1-y3)

    if abs(denom) < epsilon:
        raise ValueError('degenerate triangle in barycentricXY call')
    
    lam1 = ((y2-y3)*(x-x3) + (x3-x2)*(y-y3))/denom
    lam2 = ((y3-y1)*(x-x3) + (x1-x3)*(y-y3))/denom
    lam3 = 1.0-(lam1+lam2)

    return [lam1,lam2,lam3]

def barycentricXY(a,p1,p2,p3):
    """Given a point ``a`` and three vertices ``p1``, ``p2``, and ``p3``,
    all of which fall into the same XY plane, compute the barycentric
    coordinates lam1, lam2, and lam3

    """
    if not isXYPlanar([a,p1,p2,p3]):
        raise ValueError('non-XY points in barycentricXY call')
    return _barycentricXY(a,p1,p2,p3)

## Functions that operate on simple polylines/polygons
## ---------------------------------------------------

## A simple polyline/polygon is a list of three or more points, which
## define the verticies of the line.  If the last point and the first
## point are coincident, i.e. dist(points[0],points[-1]) < epsilon,
## then the line will be interpreted as a closed polygon.

## polygons with positive area are defined in right-hand
## (counterclockwise) order

## polylines can be sampled, can be intersected with arcs, lines, and
## other polylines, and when closed allow for inside-outside
## testing. polylines also have bounding boxes.

def poly(*args):
    """ Make a polyline or polygon by copying an existing poly or from a list of points or individual points
    """
    if len(args) == 0 or len(args) == 2:
        raise ValueError('bad number of arguments {} passed to poly()'.format(len(args)))
    if len(args) == 1:
        if ispoly(args[0]):
            return deepcopy(args[0])
        else:
            raise VauleError('non-poly list passed to poly()')
    # args is of length 3 or greater.  Check to see if args are points
    a = list(args)
    b = list(filter(lambda x: not ispoint(x),a))
    if len(b) > 0:
        raise ValueError('non-point arguments to poly(): {} '.format(b))
    return deepcopy(a)

def ispoly(a):
    """is ``a`` a poly?"""
    return isinstance(a,list) and len(a) > 2 and \
        len(list(filter(lambda x: not ispoint(x),a))) == 0

## Note: this does not test for all points lying in the same plane
def ispolygon(a):
    """is ``a`` a polygon?"""
    return ispoly(a) and dist(a[0],a[-1]) < epsilon

## if it's a polyline, sample the middle.  If it's a polygon, compute
## the barycentric cemter with equal vertex weighting.
def polycenter(a):
    """
    Return the midpoint of polyline ``a``.  In the case where ``a`` is
    a closed polygon, return the barycentric center with equal vertex
    weighting.

    """
    if ispolygon(a):
        p=point(0,0)
        l= len(a)-1
        for i in range(l):
            p = add(p,a[i])
        return scale3(p,1.0/l)
    else:
        return samplepoly(a,0.5)

## we heart poly bboxes
def polybbox(a):
    """Compute the bounding box of polyline/polygon ``a``"""
    if len(a) == 0:
        return False
    elif len(a) == 1:
        return pointbbox(a[0])
    else:
        minx = maxx = a[0][0]
        miny = maxy = a[0][1]
        for i in range(1,len(a)):
            x=a[i][0]
            y=a[i][1]
            if x < minx:
                minx =x
            elif x > maxx:
                maxx = x
            if y < miny:
                miny = y
            elif y > maxy:
                maxy = y
        return [ point(minx,miny),point(maxx,maxy)]

## only valid for closed polylines.  Count the intersections for a
## line drawn from point to test to a point outside the bounding
## box. Inside only if the number of intersections is odd.
def isinsidepolyXY(a,p):
    """
    Determine if point ``p`` lies on (within `epsilon of) the polyline ``a``.  In the 
    case where ``a`` is a closed polygon, determine if point ``p`` lies anywhere within
    the polygon by even-odd line intersection method.

    """
    closed=False

    if len(a) > 2 and dist(a[0],a[-1]) < epsilon:
        closed = True

    ## if not closed, use "unsample" test to determine if ``p`` lies
    ## on the polyline
    if not closed:
        if unsamplepoly(a,p) == False:
            return False
        else:
            return True
    ## poly is closed polygon
    bb = polybbox(a)
    ## do quick bounding box check
    if not isinsidebbox(bb,p):
        return False
    ## inside the bounding box, do intersection testing
    p2 = add([1,1,0,1],bb[1])
    if vclose(p2,p): ## did we randomly pick an outside point near the
                     ## test point?
        p2 = sub(bb[0],[1,1,0,1])
    l = line(p,p2)

    pp = intersectSimplePolyXY(l,a)
    if pp == False:
        return False
    return len(pp) % 2 == 1

## is it a polygon with all points in the same x-y plane
def ispolygonXY(a):
    """is ``a`` an XY-coplanar polygon?"""
    return ispolygon(a) and isXYPlanar(a)

## utility function used by both samplepoly() and unsamplepoly()
def __lineslength(a):
    lines=[]
    lengths=[]
    length=0
    for i in range(1,len(a)):
        p0=  a[i-1]
        p1=  a[i]
        lines.append(line(p0,p1))
        l = dist(p0,p1)
        lengths.append(l)
        length += l
    return lines,lengths,length

## map the interval 0 to 1 to the total length of all poly segments
## and return the sample point corresponding to the parameter u

def samplepoly(a,u):
    """
    Sample the poly ``a`` at parameter ``u`` and return the resulting
    point.  If ``a`` is a closed polygon, then sample at ``u % 1.0``.
    If ``a`` is not closed and ``u<0``, then samples are drawn from
    the first line in the ``u<0`` paramter space of that line.  If
    ``a`` is not closed and ``u>1``, then samples are drawn from the
    last line in the ``u>1`` paramter space of that line.

    """
    
    if not ispoly(a):
        raise ValueError('non-poly passed to samplepoly')

    closed=False

    lines,lengths,length = __lineslength(a)
    if len(a) > 2 and dist(a[0],a[-1]) < epsilon:
        closed = True

    if closed:
        u=u%1.0
        
    dst = u * length;
    d = 0
    if u < 1.0:
        for i in range(len(lines)):
            l=lengths[i]
            if dst <= d+l:
                uu = 1.0 - (d+l-dst)/l
                return sampleline(lines[i],uu)
            else:
                d+=l
    else:
        uu = (dst-length+lengths[-1])/lengths[-1]
        return sampleline(lines[-1],uu)

def segmentpoly(a,u1,u2):
    """slice out a parameterized segment from polyline or polygon ``a``
    and return this as a new polyline"""

    if not ispoly(a):
        raise ValueError('non-poly passed to segmentpoly')
    closed=False
    if len(a) > 2 and dist(a[0],a[-1]) < epsilon:
        closed = True

    if closed:
        u1%= 1.0
        u2%= 1.0
        if u2 < u1:
            return segmentpoly(a,u1,1.0-epsilon) + segmentpoly(a,0.0,u2)
        
    else:
        if u1 < 0 or u2 < 0 or u1 > 1 or u2 > 1:
            raise ValueError('parameters fall outside 0,1 interval: {},{}'.format(u1,u2))
    if u2 < u1:
        a = deepcopy(a).reverse()
        u1 = 1.0-u1
        u2 = 1.0-u2
        
    lines,lengths,length = __lineslength(a)
    sgl = segmentgeomlist(lines,u1,u2)

    ply = []
    for l in sgl:
        ply.append(l[0])
    ply.append(sgl[-1][1])

    return ply

    
        
## anlogous to the unsampleline() and unsamplearc() functions, given a
## point on a poly, return the corresponding sample parameter, or
## False if the point is more than epsilon away from any poly line
## segment

def unsamplepoly(a,p):
    """
    anlogous to the unsampleline() and unsamplearc() functions, given a
    point on a poly, return the corresponding sample parameter, or
    False if the point is more than epsilon away from any poly line
    segment
    """
    if not ispoly(a):
        raise ValueError('non-poly passed to unsamplepoly')
    closed = False
    if len(a) > 2 and dist(a[0],a[-1]) < epsilon:
        closed = True
    lines,lengths,length = __lineslength(a)
    if len(lines) == 1:
        return unsampleline(lines[0],p)
    else:
        uu1 = unsampleline(lines[0],p)
        if not isinstance(uu1,bool) and uu1 < 1.0 and\
           (not closed or uu1 >= 0.0):
            return uu1*lengths[0]/length
        dst = lengths[0]
        if len(lines) > 2:
            for i in range(1,len(lines)-1):
                uu = unsampleline(lines[i],p)
                if not isinstance(uu,bool) and uu >= 0.0 and uu < 1.0:
                    return (uu*lengths[i]+dst)/length
                else:
                    dst = dst+lengths[i]
        uu2 = unsampleline(lines[-1],p)
        
        if not isinstance(uu2,bool) and uu2 >= 0.0 and\
           (not closed or uu2 <= 1.0):
            return (uu2*lengths[-1]+dst)/length
        else:
            return False


## ccompute the intersection between non-compound geometric element g,
## and poly a.
def intersectSimplePolyXY(g,a,inside=True,params=False):
    if not (ispoly(a) and isXYPlanar(a)):
        raise ValueError('non-XY-planar or bad poly argument to intersectSimplePOlyXY: {}'.format(vstr(p)))
    closed = False
    ARC=False
    LINE=False
    if len(a) > 2 and dist(a[0],a[-1]) < epsilon:
        closed = True
    if isline(g):
        LINE=True
    elif isarc(g):
        ARC=True
    else:
        raise ValueError('bad non-line or non-arc argument passed to intersectSimplePolyXY: {}'.format(vstr(g)))
    pnts = []
    uu1s = []
    uu2s = []
    lines, lengths, leng = __lineslength(a)
    if len(lines) == 1:
        if LINE:
            return lineLineIntersectXY(lines[0],g,inside,params)
        else:
            return lineArcIntersectXY(lines[0],g,inside,params)
    dst = 0.0
    if len(lines) > 2:
        for i in range(len(lines)):
            if LINE:
                uu = lineLineIntersectXY(lines[i],g,params=True)
                if not isinstance(uu,bool) and \
                   (((closed or (i > 0 and i < len(lines)-1)) and \
                     uu[0] >= 0.0 and uu[0] <= 1.0) or\
                    (not closed and i == 0 and uu[0] <= 1.0) or\
                    (not closed and i == len(lines)-1 and uu[0] >= 0.0)):
                    if params:
                        uu1s.append((uu[0]*lengths[i]+dst)/leng)
                        uu2s.append(uu[1])
                    else:
                        if (not inside) or \
                           (uu[0] >= 0.0 and uu[0] <= 1.0) and \
                           (uu[1] >= 0.0 and uu[1] <= 1.0):
                            pnts.append(sampleline(g,uu[1]))
            elif ARC:
                uu = lineArcIntersectXY(lines[i],g,params=True)
                if not isinstance(uu,bool):
                    for j in range(len(uu[0])):
                        if (((closed or (i > 0 and i < len(lines)-1)) and \
                             uu[0][j] >= 0.0 and uu[0][j] <= 1.0) or\
                            (not closed and i == 0 and uu[0][j] <= 1.0) or\
                            (not closed and \
                             i == len(lines)-1 and uu[0][j] >= 0.0)):
                            if params:
                                uu1s.append((uu[0][j]*lengths[i]+dst)/leng)
                                uu2s.append(uu[1][j])
                            else:
                                if (not inside) or \
                                   (uu[0][j] >= 0.0 and uu[0][j] <= 1.0) and \
                                   (uu[1][j] >= 0.0 and uu[1][j] <= 1.0):
                                    pnts.append(samplearc(g,uu[1][j]))

            else:
                raise ValueError('unknown geometry type -- should never happen here')
            dst = dst+lengths[i]

    if params:
        if len(uu1s) > 0:
            return [ uu2s, uu1s]
        else:
            return False
    else:
        if len(pnts) > 0:
            return pnts
        else:
            return False
        
    
def polycenter(a):
    return samplepoly(a,0.5)

def polylength(a):
    if len(a) < 2:
        return 0.0
    l = 0.0
    for i in range(1,len(a)):
        l += dist(a[i-1],a[i])
    return l

## geometry lists ----------------------- functions for lists of
## geometry.  There is no guarantee of continutiy for geometry
## elements, only that they are all valid geometry representationsm,
## or other geometry lists.

def isgeomlist(a):
    if not isinstance(a,list):
        return False
    b = list(filter(lambda x: not (ispoint(x) or isline(x) \
                                   or isarc(x) or ispoly(x) \
                                   or isgeomlist(x)),a))
    return not len(b) > 0


def __geomlistlength(gl):
    leng=0.0
    lengths=[]
    for g in gl:
        l = length(g)
        lengths.append(l)
        leng=leng+l
    return lengths,leng

## map the interval 0 to 1 to the total length of all geometry list
## elements and return the sample point corresponding to the parameter
## u.  Note that points elements are ignored for the purpose of
## sampling.

def samplegeomlist(gl,u):
    if not isgeomlist(gl):
        raise ValueError('non-geomlist passed to samplegeomlist')
    
    lengths,leng = __geomlistlength(gl)
        
    dst = u * leng
    d = 0
    if u < 1.0:
        for i in range(len(gl)):
            l=lengths[i]
            if dst <= d+l:
                uu = 1.0 - (d+l-dst)/l
                return sample(gl[i],uu)
            else:
                d+=l
    else:
        uu = (dst-leng+lengths[-1])/lengths[-1]
        return sample(gl[-1],uu)

## function to reverse a geometry list (presumably contiguous) for
## sampling purposes
def reverseGeomList(gl):
    # print("flippy")
    rgl = []
    
    for g in gl:
        if ispoint(g):
            rgl.insert(0,point(g))
        elif isline(g):
            rgl.insert(0,line(g[1],g[0]))
        elif isarc(g):
            c = deepcopy(g)
            c[1][3] = -2 # set samplereverse flag
            rgl.insert(0,c)
        elif ispoly(g):
            ply = deepcopy(g).reverse()
            rgl.insert(0,ply)
        elif isgeomlist(g):
            rgl.insert(0,reverseGeomList(g))
        else:
            raise ValueError("don't know what to do with {}".format(g))
    return rgl
    
## given a geometry list paramaterized over a 0,1 interval, return a
## geometry list corresponding to the interval 0 <= u1 <= u2 <=1.0
def segmentgeomlist(gl,u1,u2,closed=False,reverse=False):
    if closed:
        u1 %= 1.0
        u2 %= 1.0
        if u2 < u1:
            if reverse:
                return segmentgeomlist(gl,0.0,u2,closed=False,reverse=True) +\
                    segmentgeomlist(gl,u1,1.0-epsilon,closed=False,reverse=True)
            else:
                return segmentgeomlist(gl,u1,1.0-epsilon,closed=False) \
                    + segmentgeomlist(gl,0.0,u2,closed=False)
    if u1 < 0 or u1 > u2 or u2 < 0 or u2 > 1.0:
        raise ValueError('bad parameters {} and {} passed to segmentgeomlist'.format(u1,u2))
    lengths, leng = __geomlistlength(gl)

    STARTED = False
    DONE= False
    dst1 = u1 * leng
    dst2 = u2 * leng
    d = 0
    rgl = []
    def _segment(g,u1,u2):
        if isline(g):
            return segmentline(g,u1,u2)
        elif isarc(g):
            return segmentarc(g,u1,u2)
        elif ispoly(g):
            return segmentpoly(g,u1,u2)
        elif isgeomlist(g):
            return segmentgeomlist(g,u1,u2)
        else:
            raise NotImplementedError("don't know how to segment {}".format(g))
        
    for i in range(len(gl)):
        l=lengths[i]
        if dst1 <= d+l:
            if not STARTED:
                uu = 1.0 - (d+l-dst1)/l
                uu2 = 1.0
                if dst2 < d+l:
                    uu2 = 1.0 - (d+l-dst2)/l
                    DONE=True
                rgl.append(_segment(gl[i],uu,uu2))
                STARTED=True
                                      
            elif dst2 <= d+l:
                uu = 1.0 - (d+l-dst2)/l
                rgl.append(_segment(gl[i],0.0,uu))
                DONE=True
            else:
                rgl.append(deepcopy(gl[i]))
        if DONE:
            if reverse:
                rgl = reverseGeomList(rgl)
            return rgl
        d+=l
    if reverse:
        rgl = reverseGeomList(rgl)
    return rgl


def geomlistbbox(gl):
    if not isgeomlist(gl):
        raise ValueError('non geomlist passed to geomlistbbox: {}'.format(gl))
    ply=[]
    for g in gl:
        if ispoint(g):
            ply.append(g)
        elif isline(g) or ispoly(g):
            ply = ply+g
        elif isarc(g):
            ply = ply + arcbbox(g)
        elif isgeomlist(g):
            ply = ply + geomlistbbox(g)
        else:
            raise ValueError('bad thing in geomlist passed to geomlistbbox -- should never happen')
    return polybbox(ply)

## determine if a point lies inside closed regions of a geometry list.
## only valid for geometry lists with only closed regions.  Count the
## intersections for a line drawn from point to test to a point
## outside the bounding box. Inside only if the number of
## intersections is odd.
def isinsidegeomlistXY(a,p):

    bb = geomlistbbox(a)
    if not isinsidebbox(bb,p):
        return False
    p2 = add([1,1,0,1],bb[1])
    if vclose(p2,p): ## did we randomly pick an outside point near the
                     ## test point?
        p2 = sub(bb[0],[1,1,0,1])
    l = line(p,p2)

    pp = intersectGeomListXY(l,a)
    if pp == False:
        return False
    return len(pp) % 2 == 1


## determint if the contents of a geometry list lie in the same x-y
## plane
def isgeomlistXYPlanar(gl):
    if not isgeomlist(gl):
        return False
    pp = []
    for g in gl:
        if ispoint(g):
            pp.append(g)
        elif isline(g) or ispoly(g):
            pp = pp + g
        elif isarc(g):
            pp.append(g[0])
            if len(g) > 2 and not vclose(g[2],point(0,0,1)):
                return False
        elif isgeomlist(g):
            if not isgeomlistXYPlanar(g):
                return False
    return isXYPlanar(pp)

## compute the intersection between geometric element g, and geometry
## list gl.  NOTE: this function does not impose continuity
## requirements on the geometry list, and point elements are ignored
## for intersection testing.

def intersectGeomListXY(g,gl,inside=True,params=False):
    if not isgeomlistXYPlanar(gl + [ g ]):
        raise ValueError('non-XY-planar or geometry arguments to intersectSimpleGeomListXY: {}'.format(vstr(gl + [ g])))
    gTypes = ('point','line','arc','poly','glist')
    def typeString(geom):
        if ispoint(geom):
            return 'point'
        elif isline(geom):
            return 'simple'
        elif isarc(geom):
            return 'simple'
        elif ispoly(geom):
            return 'poly'
        elif isgeomlist(g):
            return 'glist'
        else:
            return False
    gtype = typeString(g)
    if not gtype or gtype == 'point':
        raise ValueError('bad geometry argument passed to intersectGeomListXY: {}'.format(vstr(g)))
    pnts = []
    uu1s = []
    uu2s = []
    lengths, leng = __geomlistlength(gl)
    if len(gl) == 1:
        if gtype == 'simple':
            return _intersectSimpleXY(g,gl[0],inside,params)
        if gtype == 'poly':
            r = intersectSimplePolyXY(gl[0],g,inside,params)
            if params:
                return [ r[1], r[0] ]
            else:
                return r
        if gtype == 'glist':
            r = intersectGeomListXY(gl[0],g,inside,params)
            if params:
                return [ r[1], r[0] ]
            else:
                return r
        else:
            raise ValueError('bad gtype in intersectGeomListXY, this should never happen')
        
    dst = 0.0
    if len(gl) > 2:
        for i in range(len(gl)):
            g2 = gl[i]
            uu = []
            gtype2 = typeString(g2)
            if not gtype2:
                raise ValueError('This should never happen: bad contents of geometry list after checks in intersectGeomListXY')
            if gtype2 == 'point':
                continue
            elif gtype == 'simple' and  gtype2 == 'simple':
                try:
                    uu = _intersectSimpleXY(g,g2,params=True)
                except ValueError:
                    print('simple intersection problem with: ',vstr(g),' and ',vstr(g2))
                    raise
            elif gtype == 'simple' and  gtype2 == 'poly':
                uu = intersectSimplePolyXY(g,g2,params=True)
            elif gtype == 'poly' and gtype2 == 'simple':
                z = intersectSimplePolyXY(g2,g,params=True)
                if not isinstance(z,bool):
                    uu = [z[1],z[0]]
            elif gtype == 'poly' and gtype2 == 'poly':
                zz1 = []
                zz2 = []
                for i in range(1,len(g)):
                    zz = intersectSimplePolyXY(line(g[i-1],g[i]),
                                               g2, params=True)
                    if not isinstance(zz,bool):
                        zz1.append(zz[0])
                        zz2.append(zz[1])
                if len(zz1) > 0:
                    uu = [ zz1,zz2 ]
            elif gtype == 'simple' and gtype2 == 'glist':
                uu = intersectGeomListXY(g,g2,params=True)
            elif gtype == 'glist' and gtype2 == 'simple':
                z = intersectGeomListXY(g2,g,params=True)
                if not isinstance(z,bool):
                    uu = [z[1],z[0]]
            elif gtype == 'glist' and gtype2 == 'glist':
                zz1 = []
                zz2 = []
                for gg in g:
                    zz = intersectGeomListXY(gg,g2,params=True)
                    if not isinstance(z,bool):
                        zz1.append(zz[0])
                        zz2.append(zz[1])
                if len(zz1) > 0:
                    uu = [zz1,zz2 ]
            else:
                raise ValueError('Aaaaahhhhh.... this should never happen')

            ## if there was at least one intersection
            if not isinstance(uu,bool) and len(uu) > 0:
                for j in range(len(uu[0])):
                    if params:
                        if ((not inside) or \
                            (uu[0][j] >= 0.0 and uu[0][j] <= 1.0)) and \
                            (uu[1][j] >= 0.0 and uu[1][j] <= 1.0):
                            uu1s.append(uu[0][j])
                            uu2s.append((uu[1][j]*lengths[i]+dst)/leng)
                    else:
                        if ((not inside) or \
                           (uu[0][j] >= 0.0 and uu[0][j] <= 1.0)) and \
                           (uu[1][j] >= 0.0 and uu[1][j] <= 1.0):
                            pnts.append(sample(g,uu[0][j]))
            dst = dst+lengths[i]

    if params:
        if len(uu1s) > 0:
            return [ uu1s, uu2s]
        else:
            return False
    else:
        if len(pnts) > 0:
            return pnts
        else:
            return False

    
        
## functions on trangles -- a subset of polys
## ------------------------------------------

## a triangle is a list of three points, or four points if the first
## and the last are the same.
def istriangle(a):
    return ispoly(a) and ( len(a) == 3 or \
                           len(a) == 4 and dist(a[0],a[3]) < epsilon )

## given a triangle defined as a list of three points and an
## additional test point, determine if that test point lies inside or
## outside the triangle, assuming that all points lie in the x-y plane

def _isInsideTriangleXY(a,poly): #no value checking version
    p1=poly[0]
    p2=poly[1]
    p3=poly[2]
    bary = _barycentricXY(a,p1,p2,p3)
    if bary[0] >= 0.0 and bary[0] <= 1.0 and \
       bary[1] >= 0.0 and bary[1] <= 1.0 and \
       bary[2] >= 0.0 and bary[2] <= 1.0:
        return True
    return False
        

# value-safe wrapper for triangle inside testing
def isInsideTriangleXY(a,poly):
    if not istriangle(poly) or not isXYPlanar(poly + [a ]):
        print("isInsideTriangleXY -- a: ",vstr(a)," poly: ",vstr(poly))
        print("istriangle(a): ",istriangle(a)," isXYPlanar(poly + [a ]) :",isXYPlanar(poly + [a ]))
        print("ispoly(a): ",ispoly(a))
        raise ValueError('bad triangle or non-XY points in insidetriangleXY call')
    return _isInsideTriangleXY(a,poly)

## given a convex polygon defined as a list of three or more points
## and a test point, determine if that test point lies inside or
## ouside the polygon, assuming that all points lie in the x-y plane

def _isInsideConvexPolyXY(a,poly): #no value checking version
    if len(poly) == 3:
        return _isInsideTriangleXY(a,poly)
    else:
        return _isInsideTriangleXY(a,poly[0:3]) or \
               _isInsideConvexPolyXY(a,poly[1:len(poly)])

# value-safe wrapper for convex polygon inside testing function
def isInsideConvexPolyXY(a,poly):
    if len(poly) < 3 or not isXYPlanar(poly + [ a]):
        raise ValueError('bad poly length or non-XY points in insidepolyXY call')
    return _isInsideConvexPolyXY(a,poly)

## Generalized computational geometry functions
## ----------------------------------------

## Functions that operate on points, lines, arcs, and polys,
## determining the nature of their arguments and generalizing the
## operations of sampling, finding centers, computing bounding boxes,
## and finding intersections.

def length(x):

    """
    Return the scalar length of figure x.
    """
    if ispoint(x):
        # return pointlength(x):
        return 0.0
    elif isline(x):
        return linelength(x)
    elif isarc(x):
        return arclength(x)
    elif ispoly(x):
        return polylength(x)
    elif isgeomlist(x):
        l = 0.0
        for g in x:
            l += length(g)
        return l
    else:
        raise ValueError("inappropriate type for length(): ".format(x))

def center(x):
    """
    Return the point corresponding to the center of figure x.
    """
    if ispoint(x):
        # return pointcenter(x)
        return point(x)
    elif isline(x):
        return linecenter(x)
    elif isarc(x):
        return arccenter(x)
    elif ispoly(x):
        return polycenter(x)
    elif isgeomlist(x):
        pl = []
        for g in x:
            pl.append(center(g))
        return polycenter(pl)
    else:
        raise ValueError("inappropriate type for center(): ",format(x))
    
def bbox(x):
    """
    Given a figure x, return the three-dimensional bounding box of the figure.
    """
    if ispoint(x):
        return pointbbox(x)
    elif isline(x):
        return linebbox(x)
    elif isarc(x):
        return arcbbox(x)
    elif ispoly(x):
        return polybbox(x)
    elif isgeomlist(x):
        return geomlistbbox(x)
    else:
        raise ValueError("inappropriate type for bbox(): ",format(x))

## 
def sample(x,u):
    """
    Given a figure x and a parameter u, return the point on the figure
    corresponding to the specified sampling parameter.

    """
    if ispoint(x):
        # return pointsample(x)
        return point(x)
    elif isline(x):
        return sampleline(x,u)
    elif isarc(x):
        return samplearc(x,u)
    elif ispoly(x):
        return samplepoly(x,u)
    elif isgeomlist(x):
        return samplegeomlist(x,u)
    else:
        raise ValueError("inappropriate type for sample(): " + str(x))

def unsample(x,p):
    """
    Invert the sampling operation: return the parameter corresponding
    to the closest point on the figure to p as long as the distance is
    less than epsilon.

    """
    if ispoint(x):
        if vclose(x,p):
            return 0.0
        else:
            return False
    elif isline(x):
        return unsampleline(x,p)
    elif isarc(x):
        return unsamplearc(x,p)
    elif ispoly(x):
        return unsamplepoly(x,p)
    elif isgeomlist(x):
        raise NotImplementedError('unsampling geometry lists currently not supported')
    else:
        raise ValueError("inappropriate type for unasample(): "+str(x))

def segment(x,u1,u2):
    """ given a figure x, create a new figure spanning the specified interval in the original figure
    """
    if not (isgoodnum(u1) and isgoodnum(u2)) or close(u1,u2) or u1<0 or u2 < 0 or u1 > 1 or u2 > 1:
        raise ValueError('bad parameter arguments passed to segment: '+str(u1)+', '+str(u2))
    if ispoint(x):
        return deepcopy(x)
    elif isline(x):
        return segmentline(x,u1,u2)
    elif isarc(x):
        return segmentarc(x,u1,u2)
    elif ispoly(x):
        return segmentpoly(x,u1,u2)
    elif isgeomlist(x):
        return segmentgeomlist(x,u1,u2)
    else:
        raise ValueError("inappropriate figure type for segment(): "+str(x))
    

    
    
def isinsideXY(x,p):
    if ispoint(x):
        return isinsidepointXY(x,p)
    elif isline(x):
        return isinsidelineXY(x,p)
    elif isarc(x):
        return isinsidearcXY(x,p)
    elif ispoly(x):
        return isinsidepolyXY(x,p)
    elif isgeomlist(x):
        return isinsidegeomlistXY(x,p)
    else:
        raise ValueError("bad thing passed to inside: {}".format(x))

def translate(x,delta):
    if ispoint(x):
        return add(x,delta)
    elif isline(x):
        return [add(x[0],delta),add(x[1],delta)]
    elif isarc(x):
        a = deepcopy(x)
        a[0] = add(a[0],delta)
        return a
    elif ispoly(x):
        np = []
        for p in x:
            np.append(add(p,delta))
        return np
    elif isgeomlist(x):
        ngl = []
        for g in x:
            ngl.append(translate(g,delta))
        return ngl
    else:
        raise ValueError("don't know how to translate {}".format(x))

def scale(x,sx=1.0,sy=False,sz=False,cent=point(0,0),mat=False):
    if sy == False and sz == False:
        sy = sz = sx
    if vclose(point(sx,sy,sz),point(1.0,1.0,1.0)):
        return deepcopy(x)
    if not mat:
        if vclose(cent,point(0,0,0)):
            mat = xform.Scale(sx,sy,sz)
        else:
            mat = xform.Translation(cent,inverse=True)
            mat = mat.mul(xform.Scale(sx,sy,sz))
            mat = mat.mul(xform.Translation(cent))

    if ispoint(x):
        return mat.mul(x)
    elif isline(x):
        return line(mat.mul(x[0]),
                    mat.mul(x[1]))
    elif isarc(x):
        c = arc(x)
        if not close(sx,sy):
            raise ValueError('asymmetric scaling of XY-plane arcs not allowed')
        c[0] = mat.mul(x[0])
        c[1][0] *= sx
        return c

    elif ispoly(x):
        ply = []
        for p in x:
            ply.append(mat.mul(p))
        return ply
    
    elif isgeomlist(x):
        gl = []
        for g in x:
            gl.append(scale(g,sx,sy,sz,cent,mat))
        return gl
    else:
        raise VauleError("don't know how to scale ",vstr(x))
    

    
def transform(x,m):
    if not isinstance(m,xform.Matrix):
        raise ValueError('bad transformation matrix passed to transform')
    if ispoint(x):
        return m.mul(x)

    elif isline(x):
        return [m.mul(x[0]),m.mul(x[1])]

    elif isarc(x):
        raise NotImplementedError('not yet implemented arc transformation')

    elif ispoly(x):
        rval = []
        for p in x:
            np = m.mul(p)
            rval.append(np)
        return rval
    
    elif isgeomlist(x):
        gl = []
        for g in x:
            gl.append(transform(g,m))
        return gl
    else:
        raise ValueError('bad object passed to transform')
    

# Generalized rotation function.  Note that mat should only be set during
# recursive processing of geometry lists.

def rotate(x,ang,cent=point(0,0),axis=point(0,0,1.0),mat=False):
    if close(ang,0.0):
        return deepcopy(x)
    if not mat: # if matrix isn't pre-specified, calculate it
        if vclose(cent,point(0,0,0)):
            mat = xform.Rotation(axis,ang)
        else:
            mat = xform.Translation(cent,inverse=True)
            mat = mat.mul(xform.Rotation(axis,ang))
            mat = mat.mul(xform.Translation(cent))

    # arcs are wierd, since we will have to deal with a non-trivial
    # change of basis function to handle the interpretation of "start"
    # and "end" if the axis of rotation isn't the z axis.
    if ispoint(x):
        return mat.mul(x)
    elif isline(x):
        return line(mat.mul(x[0]),
                    mat.mul(x[1]))
    elif isarc(x):
        if not vclose(axis,point(0,0,1.0)):
            raise NotImplementedError('rotation of arcs out of XY plane not yet implemented')
        c = arc(x)
        c[0] = mat.mul(x[0])
        if not iscircle(c):
            c[1][1] += ang
            c[1][2] += ang
        return c

    elif isgeomlist(x):
        gl = []
        for g in x:
            gl.append(rotate(g,ang,cent,axis,mat))
        return gl
    
    else:
        return transform(x,mat)
    

## generalized geometry mirror function

def mirror(x,plane):
    flip=point(1,1,1)
    if plane == 'xz':
        flip[1]= -1
    elif plane == 'yz':
        flip[0]= -1
    elif plane == 'xy':
        flip[2]= -1
    else:
        raise ValueError('bad reflection plane passed to mirror')

    if ispoint(x):
        return point(mul(x,flip))
    elif isarc(x):
        a2=arc(x)
        a2[0] = mul(x[0],flip)
        start = x[1][1]
        end  = x[1][2]
        if not (start == 0 and end == 360):
            ## mirror the arc segment
            ps = point(cos(start*pi2/360.0),sin(start*pi2/360.0))
            pe = point(cos(end*pi2/360.0),sin(end*pi2/360.0))
            ps = mul(ps,flip)
            pe = mul(pe,flip)
            end = (atan2(ps[1],ps[0])%pi2)*360/pi2
            start = (atan2(pe[1],pe[0])%pi2)*360/pi2
            a2[1][1]=start
            a2[1][2]=end
        return a2
    elif ispoly(x):
        ply = []
        for p in x:
            ply.append(mul(p,flip))
        return ply
    elif isgeomlist(x):
        r = []
        for xx in x:
            r.append(mirror(xx,plane))
        return r
    else:
        raise ValueError('bad thing in list passed to mirror: {}'.format(x))

    
    
## Non-value-safe intersection calculation for non-compound geometic
## elements.

## If params is true, return a list of two lists of intersection
## parameters. If not, return a list of intersection points, or False
## if no intersections

def _intersectSimpleXY(g1,g2,inside=True,params=False):
    g1line = True
    g2line = True
    if isarc(g1):
        g1line=False
    if isarc(g2):
        g2line=False
    if g1line and g2line:
        r = lineLineIntersectXY(g1,g2,inside,params)
        if isinstance(r,bool):
            return False
        elif params:
            return [ [ r[0] ],[ r[1] ]]
        else:
            return [ r ]
    elif g1line and not g2line:
        return lineArcIntersectXY(g1,g2,inside,params)
    elif g2line and not g1line:
        r = lineArcIntersectXY(g2,g1,inside,params)
        if r and params:
            return [ r[1],r[0]]
        else:
            return r
    elif not (g2line or g1line):
        return arcArcIntersectXY(g1,g2,inside,params)
    raise ValueError('this should never happen, something wrong in _intersectSimpleXY')

    
## Value-safe simple wrapper for calculation of intersection of
## non-compound geometric elements
def intersectSimpleXY(g1,g2,inside=True,params=False):
    if not (isline(g1) or isarc(g1)) \
       or not (isline(g2) or isarc(g2)):
        raise ValueError('bad geometry passed to intersectSimpleXY')
    if not isgeomlistXYPlanar([g1,g2]):
        raise ValueError('geometry not in same XY plane in intersectSimpleXY')

    return _intersectSimpleXY(g1,g2,inside,params)


## is the argument a "simple" (non-compound) geometry object
def issimple(g):
    if ispoint(g) or isline(g) or isarc(g):
        return True
    else:
        return False

## grand unified XY plane intersection calculation, for any two simple
## or compound geometry instances.  The true swiss-army-knife of
## intersection.  Booo-ya. 

def intersectXY(g1,g2,inside=True,params=False):
    if ispoint(g1) or ispoint(g2):
        return False
    elif issimple(g1):
        if issimple(g2):
            return intersectSimpleXY(g1,g2,inside,params)
        elif ispoly(g2):
            return intersectSimplePolyXY(g1,g2,inside,params)
        elif isgeomlist(g2):
            return intersectGeomListXY(g1,g2,inside,params)
        else:
            raise ValueError('bad thing passed to intersectXY, should never get here')
    elif issimple(g2):
        rr = []
        if ispoly(g1):
            rr = intersectSimplePolyXY(g2,g1,inside,params)
        elif isgeomlist(g1):
            rr = intersectGeomListXY(g2,g1,inside,params)
            
        if rr == False:
            return False
        if params:
            return [ rr[1],rr[0] ]
        else:
            return rr
    elif ispoly(g2):
        if ispoly(g1):
            closed = False
            if len(g1) > 2 and dist(g1[0],g1[-1]) < epsilon:
                closed= True
            uu1s = uu2s = []
            pnts = []
            lines, lengths, leng =  __lineslength(g1)
            if lines == []:
                return False
            if len(lines) == 1:
                return intersectSimplePolyXY(lines[0],g2,inside,params)
            dst = 0.0
            for i in range(len(lines)):
                uu = intersectSimplePolyXY(lines[i],g2,params=True)
                if not uu == False:
                    for j in range(len(uu[0])):
                        if (((closed or (i > 0 and i < len(lines)-1)) and \
                             uu[0][j] >= 0.0 and uu[0][j] <= 1.0) or\
                            (not closed and i == 0 and uu[0][j] <= 1.0) or\
                            (not closed and i == len(lines)-1 and\
                             uu[0][j] >= 0.0)):
                            if params:
                                uu1s.append((uu[0][j]*lengths[i]+dst)/leng)
                                uu2s.append(uu[1][j])
                            else:
                                if (not inside) or \
                                   (uu[0][j] >= 0.0 and uu[0][j] <= 1.0) and \
                                   (uu[1][j] >= 0.0 and uu[1][j] <= 1.0):
                                    pnts.append(sample(g2,uu[1][j]))
            if params:
                if len(uu1s) > 0:
                    return [ uu2s, uu1s]
                else:
                    return False
            else:
                if len(pnts) > 0:
                    return pnts
                else:
                    return False
        elif isgeomlist(g1):
            rr = intersectGeomListXY(g2,g1,inside,params)
            if rr == False:
                return False
            if params:
                return [ rr[1],rr[0] ]
            else:
                return rr
        else:
            raise ValueError('bad thing, this should never happen')

    elif isgeomlist(g2):
        return intersectGeomListXY(g1,g2,inside,params)
    else:
        raise ValueError('very bad thing, this should never happen')
    
        
            
            
        
    

## UNIT TESTS
## ========================================
## check to see if we have been invoked on the command line
## if so, run some tests

if __name__ == "__main__":
    print("------------------------------------------")
    print("yapCAD geometry tests")
    print("------------------------------------------")
    print("light-weight unit tests for geom.py module")
    a = point(5,0)
    b = point(0,5)
    c = point(-3,-3)
    d = point(1,1)
    e = point(10,10)
    # print("---> point creation and testing")
    # print("some points: a:" + vstr(a) + ", b:"
    #       + vstr(b) + ", c:" + vstr(c) + ", d:" + vstr(d) + ", e:" + vstr(e))
    # print("ispoint(a): ",ispoint(a))
    # print("ispoint(point(a)): ",ispoint(point(a)))
    # print("ispoint([1,2]: ",ispoint([1,2]))
    # print("ispoint(vect(1,2)): ",ispoint(vect(1,2)))
    # print("ispoint(vect(1,2,3,4)): ",ispoint(vect(1,2,3,4)))
    # print("ispoint(vect(1,2,3,-1)): ",ispoint(vect(1,2,3,-1)))
    
    print("---> basic vector operations tests")
    # print("mag a: " + str(mag(a)))
    # print("add(a,b): " + vstr(add(a,b)))
    # print("sub(a,b): " + vstr(sub(a,b)))
    # print("mag(sub(a,b)): " + str(mag(sub(a,b))))
    # print("mag(sub(a,b)) == sqrt(50): " + str(mag(sub(a,b))==sqrt(50.0)))

    print("---> line creation and testing")
    l1 = [a,b]
    # print("l1 = [a,b] -- l1:",vstr(l1))
    # l11 = line(a,b)
    # print("l11 = line(a,b) -- l1:",vstr(l11))
    l2 = [c,d]
    # l22 = line(l2)
    # print("l2 = [c,d], l22 = line(l2) -- l22: ",vstr(l22))
    l3 = line(c,e)
    # print("l3 = line(c,e), isline(l3) : ",isline(l3))
    # print("a {}, isline(a): {}".format(vstr(a),isline(a)))
    
    
    print("---> vector and geometry copying tests")
    foo = [a,b,l3,l2,d]
    print("foo: ",vstr(foo))
    print("deepcopy(foo)",vstr(deepcopy(foo)))
    bar = [a,b,[1,2],l2,l3]
    print("bar: ",vstr(bar))
    print("expect False: deepcopy(bar)",vstr(deepcopy(bar)))
        
    print("---> line-line intersection tests")
    # print("l1:" + vstr(l1) + ", l2:" + vstr(l2) +", l3:" + vstr(l3))

    # int0 = lineLineIntersectXY(l1,l1)
    # int1 = lineLineIntersectXY(l1,l2,False)
    # int2 = lineLineIntersectXY(l1,l2,True)
    # int3 = lineLineIntersectXY(l1,l3,True)

    # print("expect False: lineLineIntersectXY(l1,l1): " + vstr(int0))
    # print("expect [2.5, 2.5]: lineLineIntersectXY(l1,l2,False): " + vstr(int1))
    # print("expect False: lineLineIntersectXY(l1,l2,True): " + vstr(int2))
    # print("expect [2.5, 2.5]: lineLineIntersectXY(l1,l3,True): " + vstr(int3))
    
    # print("linePointXY(l1,vect(0,0)): "
    #       + vstr(linePointXY(l1,vect(0,0))))
    # print("linePointXYDist(l1,vect(0,0)) == sqrt(12.5): "
    #       + vstr(abs(linePointXYDist(l1,vect(0,0))-sqrt(12.5))<epsilon))
    # print("linePointXY(l1,vect(0,10),False): "
    #       + vstr(linePointXY(l1,vect(vect(0,10)),False)))
    # print("linePointXY(l1,vect(0,10),True): "
    #       + vstr(linePointXY(l1,vect(0,10),True)))
    # print("linePointXY(l1,vect(10,0),False): "
    #       + vstr(linePointXY(l1,vect(10,0),False)))
    # print("linePointXY(l1,vect(10,0),True): "
    #       + vstr(linePointXY(l1,vect(10,0),True)))

    print("---> arc creation and testing")
    # arc1=[vect(2.5,2.5),vect(2.5,90.0,270.0,-1)]
    # print("arc1=[vect(2.5,2.5),vect(2.5,90.0,270.0,-1)], arc1: ",vstr(arc1))
    # arc11=arc(vect(2.5,2.5),2.5,90.0,270.0)
    # print("arc11=arc(vect(2.5,2.5),2.5,90.0,270.0), arc11: ",vstr(arc11))
    # print("isarc(arc1): {}  isarc(arc11): {}".format(isarc(arc1),isarc(arc11)))
    # arc12=arc(arc11)
    # print("arc12=arc(arc11), arc12: {}, isarc(arc12): {}".format(vstr(arc12),isarc(arc12)))
    # try:
    #     print("try creating an arc with a negative radius, should raise ValueError")
    #     print("arc(vect(0,0),-2): ",arc(vect(0,0),-2))
    # except ValueError as err:
    #     print('got expected result:',err)
    # print("--> line-arc disambiguation")
    # print("l1: ",vstr(l1)," arc1: ",vstr(arc1))
    # print("isline(l1): {} isline(arc1): {}".format(isline(l1),isline(arc1)))
    # print("isarc(l1): {} isarc(arc1): {}".format(isarc(l1),isarc(arc1)))
    # print("---> arc-line intersection tests")
    arc1=[vect(2.5,2.5),vect(2.5,90.0,270.0)]
    print("arc1: {}".format(vstr(arc1)))
    print("l1: {}".format(vstr(l1)))
    l2[1]=vect(0,0)
    print("l2: {}".format(vstr(l2)))
    int4 = lineArcIntersectXY(l1,arc1,False)
    int5 = lineArcIntersectXY(l1,arc1,True)
    int6 = lineArcIntersectXY([vect(0,5),vect(5,5)],arc1,True)
    int7 = lineArcIntersectXY(l2,arc1,True)
    int8 = lineArcIntersectXY(l2,arc1,False)
    print("lineArcIntersectXY(l1,arc1,False): {}".format(vstr(int4)))
    print("lineArcIntersectXY(l1,arc1,True): {}".format(vstr(int5)))
    print("lineArcIntersectXY([vect(0,5),vect(5,5)],arc1,True): {}".format(vstr(int6)))
    print("lineArcIntersectXY(l2,arc1,False): {}".format(vstr(int7)))
    print("lineArcIntersectXY(l2,arc1,True): {}".format(vstr(int8)))

    print("---> circle-circle tangent testing")
    circ1 = arc(point(5,5),5)
    circ2 = arc(point(-5,5),7.5)
    circ3 = arc(point(0,0),1)

    tl1 = circleCircleTangentsXY(circ1,circ2)
    tl2 = circleCircleTangentsXY(circ2,circ1)
    tl3 = circleCircleTangentsXY(circ3,circ2)

    print("circ1: ",vstr(circ1))
    print("circ2: ",vstr(circ2))
    print("circ3: ",vstr(circ3))
    
    print("circleCircleTangentsXY(circ1,circ2) :", vstr(tl1))
    print("circleCircleTangentsXY(circ2,circ1) :", vstr(tl2))
    print("circleCircleTangentsXY(circ3,circ2) :", vstr(tl3))

    print("---> arc-arc intersection tests")
    arc2=[vect(4.0,2.5),vect(2.5,90.0,270.0)]
    print("arc1: {}".format(vstr(arc1)))
    print("arc2: {}".format(vstr(arc2)))

    int9 = arcArcIntersectXY(arc1,arc2,False)
    int10 = arcArcIntersectXY(arc1,arc2,True)
    int11 = arcArcIntersectXY(arc2,arc1,True)
    print("arcArcIntersectXY(arc1,arc2,False):",vstr(int9))
    print("arcArcIntersectXY(arc1,arc2,True):",vstr(int10))
    print("arcArcIntersectXY(arc2,arc1,True):",vstr(int11))

    
    print("---> do some planar point testing")
    # points in the x-y plane
    p1 = point(2.5,0)
    p2 = point(5,5)
    p3 = point(0,5)
    p4 = point(2.5,10)

    # points in the x,z plane
    p5 = point(1,0,-5)
    p6 = point(10,0,10)
    
    # points in the y-z plane
    p7 = point(0,2,-1)
    p8 = point(0,10,10)

    print('expect True: isCardinalPlanar("xy",[p1,p2,p3,p4]) : {}'.format(
        isCardinalPlanar("xy",[p1,p2,p3,p4])))
    print('expect False: isCardinalPlanar("xz",[p1,p2,p3,p4]) : {}'.format(
        isCardinalPlanar("xz",[p1,p2,p3,p4])))
    print('expect False: isCardinalPlanar("yz",[p1,p2,p3,p4]) : {}'.format(
        isCardinalPlanar("yz",[p1,p2,p3,p4])))

    print('expect True: isCardinalPlanar("xz",[p1,p5,p6]) : {}'.format(
        isCardinalPlanar("xz",[p1,p5,p6])))

    print('expect True: isCardinalPlanar("yz",[p3,p7,p8]) : {}'.format(
        isCardinalPlanar("yz",[p3,p7,p8])))

    try:
        print('deliberate bad plane specification, should raise ValueError')
        print('isCardinalPlanar("FOO!",[p1,p2,p3,p4]) : {}'.format(
            isCardinalPlanar("FOO!",[p1,p2,p3,p4])))
    except ValueError as err:
        print('got expected result:',err)


    print("---> convex polygon inside testing")
    tri1 = [p1,p2,p3]
    poly1 = [p1,p2,p4,p3]
    p = point(2.5,1.0)
    q = point(2.5,5.0)
    r = point(2.5,7.0)
    s= point(-10,-10)
    t= point(10,10)
    print ("p: {}, q: {}, r: {}, s: {}, t: {}".format(vstr(p), vstr(q),
                                                      vstr(r), vstr(s),
                                                      vstr(t)))
    print ("tri1: {}".format(vstr(tri1)))
    print ("ispoly(tri1): ",ispoly(tri1))
    print ("istriangle(tri1): ",istriangle(tri1))
    print("expect True: isInsideTriangleXY(p,tri1): {}"\
          .format(isInsideTriangleXY(p,tri1)))
    print("expect True: isInsideTriangleXY(q,tri1): {}"\
          .format(isInsideTriangleXY(q,tri1)))
    print("expect False: isInsideTriangleXY(r,tri1): {}"\
          .format(isInsideTriangleXY(r,tri1)))
    print("expect False: isInsideTriangleXY(s,tri1): {}"\
          .format(isInsideTriangleXY(s,tri1)))
    print("expect False: isInsideTriangleXY(t,tri1): {}"\
          .format(isInsideTriangleXY(t,tri1)))
    print("inside poly testing")
    print ("tri1: {}".format(vstr(tri1)))
    print ("poly1: {}".format(vstr(poly1)))
    print("expect True: isInsideConvexPolyXY(p,tri1): {}"\
          .format(isInsideConvexPolyXY(p,tri1)))
    print("expect True: isInsideConvexPolyXY(q,tri1): {}"\
          .format(isInsideConvexPolyXY(q,tri1)))
    print("expect False: isInsideConvexPolyXY(r,tri1): {}"\
          .format(isInsideConvexPolyXY(r,tri1)))
    
    print("expect True: isInsideConvexPolyXY(p,poly1): {}"\
          .format(isInsideConvexPolyXY(p,poly1)))
    print("expect True: isInsideConvexPolyXY(q,poly1): {}"\
          .format(isInsideConvexPolyXY(q,poly1)))
    print("expect True: isInsideConvexPolyXY(r,poly1): {}"\
          .format(isInsideConvexPolyXY(r,poly1)))
    print("expect False: isInsideConvexPolyXY(s,tri1): {}"\
          .format(isInsideConvexPolyXY(s,tri1)))
    print("expect False: isInsideConvexPolyXY(t,tri1): {}"\
          .format(isInsideConvexPolyXY(t,tri1)))
    
    print("done!")
