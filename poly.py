# pure geometry polyline and polygon support for yapCAD
from geom import *

## the goal is to describe open and closed multi-element geometric
## figures that are specified in terms of point, line, and arc
## primitives.  For example, a rounded rectangle or rounded triangle.
## These figures can be sampled as though they are continuous curves
## and can be subject to operations like "grow" or "intersect."  The
## pure geometry will be the foundation of creating drawable
## representations.

class poly:
    """simple base class for multi-element geometric figures"""
    elements=[]
    closed=True

    def __repr__(self):
        return 'an abstract multi-element geometric figure'

    ## add another drawing element
    def addPoint(self,element):
        if isvect(element):
            self.elements.append(element)
        else:
            raise ValueError('attempt to add a non point, line or arc to poly')
    
    ## add another drawing element
    def addLine(self,element):
        if isline(element):
            self.elements.append(element)
        else:
            raise ValueError('attempt to add a non point, line or arc to poly')
    
    ## add another drawing element
    def addArc(self,element):
        if isarc(element):
            ## mark the psuedovector with a w=-1 tag to indicate that
            ## this is an arc, not a line
            self.elements.append(element)
        else:
            raise ValueError('attempt to add a non point, line or arc to poly')

    ## compute barycentric center of figure by equally weighting the
    ## center of all of the elements
    def center(self):
        p = point(0,0)
        for e in self.elements:
            if ispoint(e):
                p = add(p,e)
            elif isline(e):
                p = add(p,linecenter(e))
            elif iscircle(e):
                p = add(p,e[0])
            elif isarc(e):
                p = add(p,arccenter(e))
            else:
                raise ValueError('bad element in poly elements')
        if len(self.elements) > 0:
            p = scale(p,1.0/length(self.elements))
        return p
                
    def remove(self,element):
        self.elements.remove(element)
                
    ## function to take the elements in the elements[] list and
    ## construct the full outline.  For example, consider a list of
    ## three points in the elements[] list -- this would result in the
    ## construction of three lines that results in a closed curve.
    ## Similarly, three circles in the elements list would result in a
    ## "rounded triangle" composed of the circle-tangent lines joined
    ## by three arcs.  Elements that are explicitly specified lines or
    ## non-circular arcs are joined to adjacent elements by lines
    
    def makeoutline(self):
        def _calclength():
            l = 0
            ll = []
            length = 0
            for o in self.outline:
                if isarc(o) and o[1][3] == -1:  # arc with "marked" psuedovector
                    length=arclength(o)
                elif isline(o): # line
                    length=linelength(o)
                else:
                    raise ValueError("bad element in outline list for _calclength()")
                l += length
                ll.append(length)
            self.length = l
            self.lengths=ll
            
        def _fromPointAdd(p0,e1,e2):
            if isvect(e1): #r1 is a point -- simplest case
                ll = dist(p0,e1)
                self.outline.append(line(p0,e1))
            elif isline(e1):
                self.outline.append(line(p0,e1[0]))
                self.outline.append(line(e1))
            elif isarc(e1) and not iscircle(e1):
                p = samplearc(e1,0.0)
                self.outline.append(line(p0,p))
                self.outline.append(arc(e1))
            elif iscircle(e1): ## complicated
                c = self.center() # center of the current figure
                ## get the two tangent lines from point p0 to the
                ## circle e1
                ll = pointCircleTangentsXY(p0,e1)
                d1 = dist(c,linecenter(ll[0]))
                d2 = dist(c,linecenter(ll[1]))
                l=[]
                if d1 < d2:
                    l = ll[1]
                else:
                    l = ll[0]
                self.outline.append(line(l))
                self.outline.append(arc(e1))
            else:
                raise ValueError('bad object in element list')
            
        self.outline=[]
        self.length=0
        self.lengths=[]

        if len(self.elements) < 2:
            ## if one element, the outline is the element
            if len(self.elements) == 1:
                self.outline=deepcopy(self.elements)
                self._calclength()
            ## if there are fewer than one elements, there is nothing to do
            return

        ## we have work to do.  Construct all lines and circle-tangent
        ## lines first, inserting explicit line and arc elements.
        ## Then go back and convert circles to arcs.  Finally compute
        ## total length and length list for use in sample() function.

        for i in range(1,len(self.elements)):
            e0 = self.elements[i-1]
            e1 = self.elements[i]
            if i == len(self.elements)-1: #last item
                e2 = self.elements[0]
            else:
                e2 = self.elements[i+1]
            ## work through element types
            if isvect(e0): #e0 is a point
                _fromPointAdd(e0,e1,e2)
            elif isline(e0):
                _fromPointAdd(e0[1],e1,e2)
            elif isarc(e0) and not iscircle(e0):
                p0 = samplearc(e0,1.0)
                _fromPointAdd(p0,e1,e2)
            elif iscircle(e0):
                c = self.center()
                ll = circleCircleTangentsXY(e0,e1)
                d1 = dist(c,linecenter(ll[0]))
                d2 = dist(c,linecenter(ll[1]))
                l = []
                if d1 < d2:
                    l = ll[1]
                else:
                    l = ll[0]
                self.outline.append(line(l))
                self.outline.append(arc(e1)) # note: this appends the
                                             # full circle.  We will
                                             # fix it in post. :)

        ## OK, now go back through and replace full circles with arcs
        ## and catch any intersecting lines due to non-convex
        ## curvature
        for i in range(1,len(self.outline)):
            e0 = self.outline[i-1]
            e1 = self.outline[i]
            if i == len(self.outline)-1: #last item
                e2 = self.outline[0]
            else:
                e2 = self.outline[i+1]
            if iscircle(e1):
                if not isline(e0) or not isline(e2):
                    raise ValueError('circle not bracketed by lines')
                pp0 = lineArcIntersectXY(e0,e1,False)
                pp1 = lineArcIntersectXY(e2,e1,False)
                # these should be tangent lines, so exactly one
                # intersection each.  If not, we have a problem
                if len(pp0) != 1 or len(pp1) != 1:
                    raise ValueError('bad line-circle intersection in poly outline calculation')
                # we are assuming that elements are ordered in a
                # counter-clockwise fashion
                p0=sub(pp0[0],e1[0])
                p1=sub(pp1[0],e1[0])
                start = (atan2(p0[1],p0[0]) % pi2) * 360.0/pi2
                end = (atan2(p1[1],p1[0]) % pi2) * 360.0/pi2
                if end < start:
                    raise NotImplementedError("not handling non-convex case for circle-tangent lines in poly outline yet")
                newarc = arc(e1[0],e1[1][0],start,end)
                self.outline[i]=newarc
                
        _calclength()
        
    def sample(self,u):
        if len(self.elements) == 0:
            return vect(0,0)
        dist = (u % 1.0) * self.length
        d = 0
        for e in self.outline:
            l=0
            if isline(e):
                l=linelength(e)
            elif isarc(e):
                l=arclength(e)
            else:
                raise ValueError('bad element in poly outline list {}'.format(e))
            if dist <= d+l:
                uu = (d+l-dist)/l
                if isline(e):
                    return sampleline(e,uu)
                else:
                    return samplearc(e,uu)
            else:
                d=d+l
                
    def bounding(self):
        return vect(0,0)

    def is_inside(self,p):
        return False

    def grow(self,r):
        return False

    def shrink(self,r):
        return False

    def lineIntersect(self,l,inside=True):
        return False

    def arcIntersect(self,c,inside=True):
        return False
    
    def polyIntersectPoints(self,other):
        return False

    def polyIntersect(self,other):
        return False

    
