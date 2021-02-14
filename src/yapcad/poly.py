## Derived Geometry clases for yapCAD
from yapcad.geom import *
from yapcad.geometry import *

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

"""
Derived Geometry Classes
========================

These subclasses of Geometry compute a figure based on the contents
of their elements list, rather than simply providing a wrapper around
the contents of their elements.  

The Polygon class is an intepolating closed figure class that will
interpolate lines between each point or open figure element in its
element list.  It will also treat circles as rounded convex conrners
of the specified radius.  

By treating circles as rounded corners, polygon instances can be
"grown" to create new figures that are a fixed distance larger, which
is very useful for offsetting drill holes from the edge of a boundary,
etc.  """

class Polygon(Geometry):
    """Multi-element closed figure derived geometry class"""

    def __init__(self,a=None):
        super().__init__()
        self._setDerived(True)
        self._setClosed(True)
        self.__lengths=0
        self.__outline=[]
        
        if isinstance(a,Polygon):
            self._setElem(deepcopy(a.elem))
            self._setUpdate(True)
        elif isgeomlist(a):
            for i in a:
                if ispoint(i) or isarc(i):
                    self.elem.append(i)
                elif isline(i):
                    pass
                elif ispoly(i):
                    self._setElem(self.elem + i)
                else:
                    raise ValueError("bad argument to Polygon constructor")

        elif a != None:
            raise ValueError("bad argument to Polygon constructor")

        elm = self.elem

        l = len(elm)
        if l > 2:
            # if first and last element are identicial points or circles,
            # remove them
            if ((ispoint(elm[0]) and ispoint(elm[-1]) and 
                 dist(elm[0],elm[-1]) < epsilon) or 
                (iscircle(elm[0]) and iscircle(elm[-1]) and 
                 abs(elm[0][1][0] - elm[-1][1][0]) < epsilon)):
                self.elem.pop()
        
    def __repr__(self):
        return 'Polygon({})'.format(vstr(self.elem))

    ## compute barycentric center of figure by equally weighting the
    ## points.  If last and first points are the same, ignore the
    ## last point.
    def _calcCenter(self):
        l = len(self.elem)
        if l == 0:
            return None
        elif l == 1:
            return center(self.elem[0]) # center is sole point
        else:
            
            if dist(center(self.elem[0]),
                    center(self.elem[-1])) < epsilon:
                l -= 1

            
            p = center(self.elem[0])
            for i in range(1,l):
                p = add(center(self.elem[i]),p)

            return scale3(p,1/l)

    
    def _updateInternals(self):
        if self.update:
            self._setUpdate(False)
            self._makeoutline()
            self._setCenter(self._calcCenter())
            self._setBbox(bbox(self.__outline))

    def addPoint(self,element):
        if ispoint(element):
            self._setUpdate(True)  # flag that we need to recalculate stuff
            self.elem.append(deepcopy(element))
        else:
            raise ValueError('attempt to add a non point to Polyline')

    ## set a point to a specified value
    def setPoint(self,i,p):
        if not ispoint(p):
            raise ValueError('bad point passed to Polyline.setPoint(): '.format(p))
        elem = self.elem
        if i < len(elem):
            elem[i]=point(p)
        elif i == len(elem):
            elem.append(p)
        else:
            raise ValueError('index out of range in PolylinesetPoint(): '.format(i))
        self._setElem(elem)
        self._setUpdate(True)

    ## add another drawing element
    def addLine(self,element):
        if isline(element):
            self._setUpdate(True)  # flag that we need to recalculate stuff
            self.elem.append(deepcopy(element))
        else:
            raise ValueError('attempt to add a non point, line or arc to poly')
    
    ## add another drawing element
    def addArc(self,element):
        if isarc(element):
            self._setUpdate(True)  # flag that we need to recalculate stuff
            self.elem.append(deepcopy(element))
        else:
            raise ValueError('attempt to add a non point, line or arc to poly')
                
    def remove(self,element):
        self._setUpdate(True)
        self.elem.remove(element)
                
    ## function to take the elements in the elem[] list and
    ## construct the full outline.  For example, consider a list of
    ## three points in the elem[] list -- this would result in the
    ## construction of three lines that results in a closed curve.
    ## Similarly, three circles in the elem list would result in a
    ## "rounded triangle" composed of the circle-tangent lines joined
    ## by three arcs.  Elements that are explicitly specified lines or
    ## non-circular arcs are joined to adjacent elements by lines

    def _makeoutline(self):
        # import pdb; pdb.set_trace()
        def _calclength():
            l = 0
            ll = []
            length = 0
            for o in self.__outline:
                if isarc(o):
                    length=arclength(o)
                elif isline(o): # line
                    length=linelength(o)
                else:
                    #print("self.__outline: ",vstr(self.__outline))
                    #print("o: ",vstr(o))
                    raise ValueError("bad element in outline list for _calclength()")
                l += length
                ll.append(length)
            self._setLength(l)
            self.__lengths=ll

        def _handleCircle(e0,e1,e2):
            ## get the two tangent lines from circle e1 to the circle
            ## e2
            ll = circleCircleTangentsXY(e1,e2)
            l=[]
            x0=center(e0)
            x1=center(e1)
            x2=center(e2)
            #print("x0: ",vstr(x0)," x1: ",vstr(x1)," x2: ",vstr(x2))
            #print("ll: ",vstr(ll))
            v1= sub(x1,x0)
            v2= sub(x2,x1)
            r0 = cross(v1,v2)
            x3 = linecenter(ll[0])
            x4 = linecenter(ll[1])
            #print("x3: ",vstr(x3),"x4: ",vstr(x4))
            #print("v1: ",vstr(v1),"v2: ",vstr(v2))
            v3 = sub(x3,x2)
            v4 = sub(x4,x2)
            r1 = cross(v2,v3)
            r2 = cross(v2,v4)
            #print("v3: ",vstr(v3)," v4: ",vstr(v4))
            #print("r0: ", r0[2], " r1: ",r1[2]," r2: ",r2[2])

            if r0[2] >= 0:
                #print("left-hand turn")
                if r1[2] >= 0:
                    l = ll[1]
                else:
                    l = ll[0]
            else:
                #print("right hand turn")
                if r2[2] >= 0:
                    l = ll[0]
                else:
                    l = ll[1]
                        
            self.__outline.append(line(l))
            
        def _fromPointAdd(e0,p1,e2):
            if ispoint(e2): #r1 is a point -- simplest case
                #ll = dist(p1,e2)
                self.__outline.append(line(p1,e2))
            elif isline(e2):
                self.__outline.append(line(p1,e2[0]))
                self.__outline.append(line(e2))
            elif isarc(e2) and not iscircle(e2):
                p = samplearc(e2,1*epsilon)
                self.__outline.append(line(p1,p))
                self.__outline.append(arc(e2))
            elif iscircle(e2): 
                _handleCircle(e0,arc(p1,1*epsilon),e2)
                self.__outline.append(arc(e2))
            else:
                raise ValueError('bad object in element list')
            
        self.__outline=[]
        self._setLength(0.0)
        self.__lengths=[]

        if len(self.elem) < 2:
            ## if one element, the outline is the element
            if len(self.elem) == 1:
                self.__outline=deepcopy(self.elem)
                _calclength()
            ## if there are fewer than one elem, there is nothing to do
            return
        ## we have work to do.  Construct all lines and circle-tangent
        ## lines first, inserting explicit line and arc elements.
        ## Then go back and convert circles to arcs.  Finally compute
        ## total length and length list for use in sample() function.

        
        for i in range(len(self.elem)):
            if i == 0:
                e0 = self.elem[-1]
            else:
                e0 = self.elem[i-1]
            e1 = self.elem[i]
            if i == len(self.elem)-1: #last item
                e2 = self.elem[0]
            else:
                e2 = self.elem[i+1]
            ## work through element types
            if ispoint(e1): #e1 is a point
                _fromPointAdd(e0,e1,e2)
            elif isline(e1):
                _fromPointAdd(e0,e1[0],e2)
            elif isarc(e1) and not iscircle(e1):
                p1 = samplearc(e1,1.0)
                _fromPointAdd(e0,p1,e2)
            elif iscircle(e1):
                c2 = []
                if ispoint(e2):
                    c2 = arc(e2,1*epsilon)
                if iscircle(e2):
                    c2 = e2
                else:
                    p = sample(e2,0.0)
                    c2 = arc(p,1*epsilon)
                _handleCircle(e0,e1,c2)
                
                if not ispoint(e2):
                    self.__outline.append(deepcopy(e2))

        ## OK, now go back through and replace full circles with arcs
        ## and catch any intersecting lines due to non-convex
        ## curvature
        for i in range(1,len(self.__outline)):
            e0 = self.__outline[i-1]
            e1 = self.__outline[i]
            if i == len(self.__outline)-1: #last item
                e2 = self.__outline[0]
            else:
                e2 = self.__outline[i+1]
            if isline(e1) and isline(e2):
                pi = lineLineIntersectXY(e1,e2,inside=True)
                if pi == False:
                    # print("not expected: adjacent lines don't intersect. We will fix that")
                    pi = lineLineIntersectXY(e1,e2,inside=False)
                    if pi == False:
                        # raise ValueError('parallel adjacent lines -- no clue')
                        print("odd parallel? adjacent line condition")
                        continue
                    else:
                        self.__outline[i] = line(e1[0],pi)
                        self.__outline[(i+1)%len(self.__outline)] = line(pi,e2[1])
            if iscircle(e1):
                if not isline(e0) or not isline(e2):
                    raise ValueError('circle not bracketed by lines')
                pp0 = lineArcIntersectXY(e0,e1,False)
                pp1 = lineArcIntersectXY(e2,e1,False)
                # these should be tangent lines, so exactly one
                # intersection each.  sometimes we get two that are
                # just a bit more than epsilon apart, which is OK.
                if len(pp0) < 1 or len(pp1) < 1:
                    raise ValueError('bad line-circle intersection in poly outline calculation')
                # we are assuming that elements are ordered in a
                # counter-clockwise fashion
                p0=sub(pp0[0],e1[0])
                p1=sub(pp1[0],e1[0])
                start = (atan2(p0[1],p0[0]) % pi2) * 360.0/pi2
                end = (atan2(p1[1],p1[0]) % pi2) * 360.0/pi2
                newarc = arc(e1[0],e1[1][0],start,end)
                self.__outline[i]=newarc

        self.__outline = cullZeroLength(self.__outline)
                
        _calclength()

    # sample function uses cached length and lengths[], otherwise
    # behavior is identical to superclass
    def sample(self,u):
        if self.update:
            self._updateInternals()
        if len(self.__outline) == 0:
            raise ValueError('no geometry to sample, empty poly')
        dist = (u % 1.0) * self.length
        d = 0
        for i in range(len(self.__outline)):
            l=self.__lengths[i]
            if dist <= d+l:
                u = 1.0 - (d+l-dist)/l
                e = self.__outline[i]
                return sample(e,u)
            else:
                d=d+l
                
    @property
    def geom(self):
        if self.update:
            self._updateInternals()
        return deepcopy(self.__outline)
                

    def grow(self,r):
        if close(r,0.0):
            return
        elif r < 0:
            raise ValueError('negative growth values not valid')
        def _growgeom(elem):
            ne = []
            for e in elem:
                if ispoint(e):
                    ne.append(arc(e,r))
                elif isline(e):
                    ne.append(arc(e[0],r))
                    ne.append(arc(e[1],r))
                elif isarc(e):
                    a = arc(e)
                    a[1][0] += r
                    ne.append(a)
                elif ispoly(e):
                    for p in e:
                        ne.append(arc(p,r))
                elif isgeomlist(e):
                    ne += _growgeom(e)
            return ne
        self._setElem(_growgeom(self.elem))
        self._setUpdate(True)

        
    def shrink(self,r):
        raise NotImplementedError('shrink operation not yet implemented')
        return False

    
## Utility functions


## function to process a geometry list and remove zero-length
## elements, including points, zero-length lines, and zero-length
## arcs.  Will recursively process geometry lists.

def cullZeroLength(geom):
    gl = []
    for g in geom:
        if ispoint(g):
            pass
        elif (isline(g) or isarc(g) or ispoly(g)) \
             and length(g) > epsilon:
            gl.append(deepcopy(g))
        elif isgeomlist(g) and length(g) > epsilon:
            gl.append(cullZeroLength(g))
        else: # zero-length or unknown element
            print("zero length element ",vstr(g)," culled")
    return gl



## make a Polygon instance containing a circle as two half-arcs

# def Circle(center=point(0,0,0),radius=1.0):
#     poly = Polygon()
#     poly.addArc(arc(center,radius,0.0,180.0))
#     poly.addArc(arc(center,radius,180.0,360.0))
#     return poly

def Circle(center=point(0,0,0),radius=1.0):
    return Geometry(arc(center,radius,0,360))

## make a rectangle with the specified width and height

def Rect(width,height,center=point(0,0,0)):
    w=width/2.0
    h=height/2.0
    p0=add(point(-w,h),center)
    p1=add(point(-w,-h),center)
    p2=add(point(w,-h),center)
    p3=add(point(w,h),center)
    poly = Polygon()
    poly.addPoint(p0)
    poly.addPoint(p1)
    poly.addPoint(p2)
    poly.addPoint(p3)
    return poly

## make rounded rectangle with specified width, height, and chamfer,
## centered at the origin by default

def RoundRect(width,height,chamf,center=point(0,0,0)):
    cr=chamf/2.0
    wid=width-chamf
    hei=height-chamf
    w=wid/2.0
    h=hei/2.0
    p0=add(point(-w,h),center)
    p1=add(point(-w,-h),center)
    p2=add(point(w,-h),center)
    p3=add(point(w,h),center)
    c0=arc(p0,cr)
    c1=arc(p1,cr)
    c2=arc(p2,cr)
    c3=arc(p3,cr)
    poly = Polygon()
    poly.addArc(c0)
    poly.addArc(c1)
    poly.addArc(c2)
    poly.addArc(c3)
    return poly
