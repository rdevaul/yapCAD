## yapCAD geometry-generating superclass
## =====================================

## Copyright (c) 2020 Richard W. DeVaul
## Copyright (c) 2020 yapCAD contributors
## All rights reserved

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

"""object-oriented computational geometry figure classes for **yapCAD**

===============
Overview
===============

The ``yapcad.geometry`` module provides the ``Geometry`` class. At its
most basic level, the ``Geometry`` class wrapps ``yapcad.geom``
figures and caches properties, such as the figure's bounding box,
length, center, etc., speeding certain computational geometry
operations.

Why would I use ``yapcad.geometry`` vs. ``yapcad.geom``
=======================================================

You might perfer the convenience and simplicity of the object-oriented
interface as a way to access the underlying power of the
``yapcad.geom`` and associated modules.  In addition, in many cases it
is actually more efficient to use the ``Geometry`` class wrappers for
figures rather than the ``yapcad.geom`` figures themselves because of
the memoizing and caching feautres provided.

For examle, Geometry class provides wrappers around the
``yapcad.geom3d poly2surface()`` function for the triangulation of
figures into triangle mesh surfaces, which is cached. This, in turn,
provides the foundations for extrusions and the construction of 3D
triangluated objects.

Likewise, for compound figures above a certain complexity threshold,
the Geometry class implements an internal quadtree decomposition of
the figure to speed intersection testing.  This is done in a lazy way,
which is to say that the quadtree is constructed the first time
intersection calculation is requested, and persists for as long as the
figure's geometry remains unchanged.  *NOTE: Fixeme, quadtree-based
intersection not yet implemented*

object oriented vs. functional approch
--------------------------------------

The ``yapcad.geom`` module generally takes a functional programming
approach, minimizing side effects.  This comes at the cost of some
redundant computation, as the representation of ``yapcad.geom``
figures are generally quite minimal, and properties such as centers,
bounding boxes, lengths, **etc.**, must be recomputed each time they
are requested, unless the user has made their own explicit cache. 

By contrast, the ``Geometry`` class caches these properties, and can
provide higher efficiency for certain types of complex operations.  In
addition, derived classes provide representations and functionality
absent from the underlying functional representations, such as
"growable" polygons, boolean operations, *etc.*.

yapcad.geometry conveneince functions
=====================================

Convenience functions, such as ``Line()``, ``Arc()``, *etc.*, operate
analogously to their uncapitalized counterparts in ``yapcad.geom``,
creating ``Geometry`` class instances wrapping the corresponding
figure.

"""

from copy import deepcopy
from yapcad.geom import *
from yapcad.geom_util import *
from yapcad.geom3d import *

class Geometry:
    """generalized computational geometry base class, also acts as a
    wrapper around yapcad.geom elements.  

    Using Geometry subclass qinstances can be more effiient than
    workng with the corresponding yapcad.geom elements, as the
    Geometry instance uses lazy evaluation and caching to speed
    repeated evaluations of various properties.

    """

    def __repr__(self):
        return f"Geometry({self.__elem})"

    def __init__(self,a=False):
        self.__update=True # do we need to update geom?
        self.__elem=[] # internal geometry
        self.__length=0.0
        self.__sampleable = False
        self.__intersectable = False
        self.__continuous = False
        self.__closed = False
        self.__center = None
        self.__bbox = None
        self.__surface = None # surface representation
        self.__surface_ang = -1 # surface parameter
        self.__surface_len = -1 # surface parameter
        self.__derived = False  # set to true for derived geometry subclasses
        if a != False:
            if ispoint(a):
                self.__elem=  deepcopy(a) 
                self.__sampleable=True
            elif isline(a) or isarc(a):
                self.__elem= deepcopy(a)
                self.__sampleable=True
                self.__intersectable=True
                self.__continuous=True
                if iscircle(a):
                    self.__closed = True
            elif ispoly(a):
                self.__elem=deepcopy(a)
                self.__sampleable=True
                self.__intersectable=True
                self.__continuous=True
                self.__closed=ispolygon(a)
            elif isgeomlist(a):
                self.__elem=deepcopy(a)
                self.__sampleable=True
                self.__intersectable=True
                self.__continuous=iscontinuousgeomlist(a)
                self.__closed=isclosedgeomlist(a)

            elif isinstance(a,Geometry):
                self.__elem = deepcopy(a.elem)
                self.__sampleable = a.issampleable()
                self.__intersectable= a.isintersectable()
                self.__continuous= a.iscontinuous()
                self.__closed = a.isclosed()
            else:
                raise ValueError(f'bad argument to Geometry class constructor: {a}')


    @property
    def sampleable(self):
        return self.__sampleable

    def _setSampleable(self,bln):
        self.__sampleable=bln

    def issampleable(self):
        """is the figure sampleable?"""
        return self.__sampleable

    @property
    def intersectable(self):
        return self.__intersectable

    def isintersectable(self):
        """is the figure intersectable?"""
        return self.__intersectable

    @property
    def continuous(self):
        return self.__continuous

    def iscontinuous(self):
        """is the figure C0 continuous over the interval [0,1]"""
        return self.__continuous

    @property
    def closed(self):
        return self.__closed

    def _setClosed(self,bln):
        self.__closed=bln

    def isclosed(self):
        """is the figure C0 continuous over the [0,1] interval and
        is ``self.sample(0.0)`` within epsilon of ``self.sample(1.0)``
        """
        return self.__closed

    @property
    def derived(self):
        return self.__derived

    def _setDerived(self,bln):
        self.__derived = bln

    def isderived(self):
        """is this an instance of a derived geometry subclass that
        computes self.geom from self.elem?  Set only in constructors
        for derived geometry subclasses
        """
        return self.__derived

    @property
    def length(self):
        """return length of figure"""
        if self.update:
            self._updateInternals()
        return self.__length

    def _setLength(self,l):
        self.__length=l

    @property
    def center(self):
        """return center of figure"""
        if self.update:
            self._updateInternals()
        return self.__center

    def _setCenter(self,c):
        self.__center=c

    @property
    def bbox(self):
        """return 3D bounding box of figure"""
        if self.update:
            self._updateInternals()
        return self.__bbox

    def _setBbox(self,bbx):
        self.__bbox=bbx

    def isinsideXY(self,p):
        """determine if a point is inside a figure.  In the case of non-closed
        figures, such as lines, determine if the point lies within
        epsilon of one of the lines of the figure.
        """
        if self.update:
            self._updateInternals()
        return isinsideXY(self.geom,p)


    def translate(self,delta):
        """apply a translation to figure"""
        self._setElem(translate(self.__elem,delta))
        self._setUpdate(True)

    def scale(self,sx=1.0,sy=False,sz=False,cent=point(0,0,0)):
        """apply a scaling to a figure"""
        self._setElem(scale(self.__elem,sx,sy,sz,cent))
        self._setUpdate(True)

    def rotate(self,ang,cent=point(0,0,0),axis=point(0,0,1.0)):
        """apply a rotation to a figure"""
        self._setElem(rotate(self.__elem,ang,cent,axis))
        self._setUpdate(True)

    def mirror(self,plane,keepSign=True):
        """apply a mirror operation to a figure.  Currently, the following
        values of "plane" are allowed: 'xz', 'yz', xy'.  Generalized
        arbitrary reflection plane specification will be added in the
        future.

        If ``keepSign == True`` (default) the sign of the area will be
        maintained, meaning that if ``mirror`` is applied to a
        right-handed closed figure (a figure with positive area) the
        resulting mirrored figure will also have positive area.  This
        is probably what you want, unless you are specifically turning
        a face into a hole, or vice-versa.

        """
        nelm = mirror(self.__elem,plane)
        if keepSign:
            nelm = reverseGeomList(nelm)
        self._setElem(nelm)
        self._setUpdate(True)

    def transform(self,m):
        """apply an arbitrary transformation to a figure, as specified by a
        transformation matrix.
        """
        self._setElem(transform(self.elem,m))
        self._setUpdate(True)
        

    # one underscore to make this easily overridable in subclasses
    def _updateInternals(self):
        """update internals: set basic attributes based on geom"""
        if self.update:
            self._setUpdate(False)
            if self.__elem == []:
                self.__length = 0.0
                self.__center = None
                self.__bbox = None
            else:
                self.__length = length(self.__elem)
                self.__center = center(self.__elem)
                self.__bbox = bbox(self.__elem)

        return

    ## more properties
    
    @property
    def update(self):
        return self.__update

    def _setUpdate(self,bln):
        """method to set status of __update, accessible from derived classes"""
        self.__update = bln

    @update.setter
    def update(self,bln):
        """property to flag update.  Ignores right-hand side of expression, always sets update flag to True"""
        self._setUpdate(True)

    @property 
    def elem(self):
        return self.__elem

    def _setElem(self,e):
        self.__elem = e

    @elem.setter
    def elem(self,e):
        self._setElem(e)

    @property
    def geom(self):
        """return yapcad.geom representation of figure"""
        if self.update:
            self._updateInternals()
        return deepcopy(self.__elem)

    def sample(self,u):
        """If the figure is sampleable, given a parameter u, return the point
        on the figure corresponding to the specified sampling
        parameter.

        """
        if not self.__sampleable:
            raise ValueError('figure is not sampleable')
        if self.update:
            self._updateInternals()
            
        gl = self.geom

        if len(gl) == 1:
            return sample(gl[0],u)
        else:
            return sample(gl,u)

    def unsample(self,p):
        """
        Invert the sampling operation: return the parameter corresponding
        to the closest point on the figure to p as long as the distance is
        less than epsilon.
        
        """
        if not self.__sampleable:
             raise ValueError('figure is not sampleable, or unsampleable for that matter')
        if self.update:
            self._updateInternals()

        gl = self.geom

        if len(gl) == 1:
            return unsample(gl[0],p)
        else:
            return unsample(gl,p)

    def segment(self,u1,u2,reverse=False):
        gl = self.geom
        if gl == []:
            raise ValueError('empty Boolean, segment not defined')
        return segmentgeomlist(gl,u1,u2,closed=self.closed,reverse=reverse)

    def intersectXY(self,g,inside=True,params=False):
        """given two XY-coplanar figures, this figure and ``g``,
        calculate the intersection of these two figures, and return a
        list of intersection points, or False if none.  

        ``g`` must be an instance of ``IntersectGeometry``, or a
        ``yapcad.geom`` figure.

        If ``inside == True``, only return intersections that are
        within the ``0 <= u <= 1.0`` interval for both figures.  If
        ``params == True``, instead of returning a list of points,
        return two lists corresponding to the sampling parameter value
        of the intersections corresponding to each figure.

        """
        if not self.isintersectable():
            raise ValueError(f'this instance {self} not intersectable')
        if self.update:
            self._updateInternals()
        if isinstance(g,Geometry):
            if g.isintersectable():
                g = g.geom
            else:
                raise ValueError(f'Geometry object {g} not intersectable')
        elif isgeomlist([g]):
            pass
        else:
            raise ValueError(f'bad thing passed to intersectXY: {g}')
        
        return intersectXY(self.geom,g,inside,params)

    
    def surface(self,minang = 5.0, minlen = 0.5):
        """
        triangulate a closed polyline or geometry list, return a surface.
        the ``minang`` parameter specifies a minimum angular resolution
        for sampling arcs, and ``minleng`` specifies a minimum distance
        between sampled points. 
        ``surface = ['surface',vertices,normals,faces]``, where:
            ``vertices`` is a list of ``yapcad.geom`` points,
            ``normals`` is a list of ``yapcad.geom`` points of the same length as ``vertices``,
            and ``faces`` is the list of faces, which is to say lists of three indices that 
            refer to the vertices of the triangle that represents each face.
        """
        if not self.isclosed():
            raise ValueError("non-closed figure has no surface representation")
        if self.update:
            self._updateInternals()
        if (self.__surface and
            close(self.__surface_ang,minang) and
            close(self.__surface_len,minlen)):
            return self.__surface
        self.__surface_ang = minang
        self.__surface_len = minlen
        geo = self.geom
        if len(geo) == 0:
            return []
        if ispolygon(geo):
            ply = geo
            holes = []
        else:
            ply, holes = geomlist2poly_with_holes(geo, minang, minlen)

        self.__surface,bnd = poly2surface(ply, holepolys=holes, checkclosed=False)

        return self.__surface

## Utility functions

def Point(x=False,y=False,z=False,w=False):
    return Geometry(point(x,y,z,w))

def Line(p1,p2=False):
    return Geometry(line(p1,p2))

def Arc(c,rp=False,sn=False,e=False,n=False,samplereverse=False):
    return Geometry(arc(c,rp,sn,e,n,samplereverse))

# def Poly(*args):
#     ply = poly(*args)
#     # print(f'ply: {ply}')
#     return Geometry(ply)

def Figure(*args):
    return Geometry(list(*args))


                    
