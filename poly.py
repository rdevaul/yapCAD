# pure geometry polyline and polygon support for yapCAD
from geom import *
from geometry import *

## The goal is to describe open and closed multi-element geometric
## figures that are specified in terms of point, line, and arc
## primitives.  For example, a "chamfered" line with arcs instead of
## sharp line intersection corners, or rounded polygons like rounded
## rectangles or rounded triangles.  These figures can be sampled as
## though they are continuous curves and can be subject to operations
## like "grow" or "intersect."

## These classes can be seen as extensions of the pure geometry
## list-based poly() and polygon() respresentations that represent
## figures with line segments



class Polyline(SampleGeometry):
    """Generalized multi-element open figure class"""

    def __init__(self,a=False):
        self._elem=[]
        self._length=0.0
        self._lengths=[]
        self._lines=[]
        self._update=True
        self._closed=False
        self._center=point(0,0,0)
        self._bbox=line(point(-epsilon,-epsilon),
                        point(epsilon,epsilon))
        if isinstance(a,Polyline):
            self._elem = deepcopy(a._elem)
            self._updateInternals()
        elif isinstance(a,list):
            for i in a:
                if ispoint(i): 
                    self._elem.append(i)
                else:
                    raise ValueError("bad argument to Polygon constructor")

    def __repr__(self):
        return 'Polyline({})'.format(vstr(self._elem))

    ## add another drawing element
    def addPoint(self,element):
        if ispoint(element):
            self._update=True  # flag that we need to recalculate stuff
            self._elem.append(element)
        else:
            raise ValueError('attempt to add a non point to Polyline')

    ## set a point to a specified value
    def setPoint(self,i,p):
        if not ispoint(p):
            raise ValueError('bad point passed to Polyline.setPoint(): '.format(p))
        if i < len(self._elem):
            self._elem[i]=point(p)
        elif i == len(self._elem):
            self._elem.append(p)
        else:
            raise ValueError('index out of range in PolylinesetPoint(): '.format(i))
        self._update=True
        
    ## return a copy of the elem list
    def getElem(self):
        return deepcopy(self._elem)
    
    def _updateInternals(self):
        if self._update:
            self._updateCenter()
            self._updateLines()
            self._update=False


    ## compute barycentric center of figure by equally weighting the
    ## points.  If last and first points are the same, ignore the
    ## last point.
    def _updateCenter(self):
        l = len(self._elem)
        if l == 0:
            return
        elif l == 1:
            self._center = center(self._elem[0]) # center is sole point
            e = point(epsilon,epsilon)            # make bounding box
            self._bbox = line(sub(self._center,e),
                             add(self._center,e))
        else:
            
            if dist(center(self._elem[0]),
                    center(self._elem[-1])) < epsilon:
                l -= 1

            p = point(0,0,0)
            for i in range(l):
                p = add(center(self._elem[i]),p)

            self._center = scale(p,1/l)
        return

    ## return the center of the figure.  If necessary, recompute that
    def getCenter(self):
        if self._update:
            self._updateInternals()
        return self._center


    def getLength(self):
        if self._update:
            self._updateInternals()
        return self._length
    
    ## function to take the points in the elem[] list and build a
    ## list of the individual lines, line lengths and total length to
    ## facilitate sampling
    
    def _updateLines(self):
        for i in range(1,len(self._elem)):
            p0=  center(self._elem[i-1])
            p1=  center(self._elem[i])
            self._lines.append(line(p0,p1))
            l = dist(p0,p1)
            self._lengths.append(l)
            self._length += l
        if len(self._elem) > 2 and dist(center(self._elem[0]),
                                            center(self._elem[-1])) < epsilon:
            self._closed = True

    def sample(self,u):
        dist = u * self._length;
        d = 0
        if self._closed:
            u = u%1.0
        if u < 1.0:
            for i in range(len(self._lines)):
                l=self._lengths[i]
                if dist <= d+l:
                    uu = 1.0 - (d+l-dist)/l
                    return sampleline(self._lines[i],uu)
                else:
                    d+=l
        else:
            uu = (dist-self._length+self._lengths[-1])/self._lengths[-1]
            return sampleline(self._lines[-1],uu)
        

class Polygon(Polyline):
    """Generalized multi-element closed figure class"""

    def __init__(self,a=False):
        super().__init__(a)
        self.outline=[]

        elm = self._elem

        l = len(elm)
        if l > 2:
            # if first and last element are identicial points or circles,
            # remove them
            if ispoint(elm[0]) and ispoint(elm[-1]) and \
               dist(elm[0],elm[-1]) < epsilon or \
               iscircle(elm[0]) and iscircle(elm[-1]) and \
               abs(elm[0][1][0] - elm[-1][1][0]) < epsilon:
                self._elem.pop()
        
    def __repr__(self):
        return 'Polygon({})'.format(vstr(self._elem))

    ## add another drawing element
    def addLine(self,element):
        if isline(element):
            self._elem.append(element)
        else:
            raise ValueError('attempt to add a non point, line or arc to poly')
    
    ## add another drawing element
    def addArc(self,element):
        if isarc(element):
            ## mark the psuedovector with a w=-1 tag to indicate that
            ## this is an arc, not a line
            self._elem.append(element)
        else:
            raise ValueError('attempt to add a non point, line or arc to poly')

    ## compute barycentric center of figure by equally weighting the
    ## center of all of the elem
    def center(self):
        p = point(0,0)
        for e in self._elem:
            if ispoint(e):
                p = add(p,e)
            elif isline(e):
                p = add(p,linecenter(e))
            elif iscircle(e):
                p = add(p,e[0])
            elif isarc(e):
                p = add(p,arccenter(e))
            else:
                raise ValueError('bad element in poly elem')
        if len(self._elem) > 0:
            p = scale(p,1.0/len(self._elem))
        return p
                
    def remove(self,element):
        self._elem.remove(element)
                
    ## function to take the elements in the elem[] list and
    ## construct the full outline.  For example, consider a list of
    ## three points in the elem[] list -- this would result in the
    ## construction of three lines that results in a closed curve.
    ## Similarly, three circles in the elem list would result in a
    ## "rounded triangle" composed of the circle-tangent lines joined
    ## by three arcs.  Elements that are explicitly specified lines or
    ## non-circular arcs are joined to adjacent elements by lines
    
    def makeoutline(self):
        def _calclength():
            l = 0
            ll = []
            length = 0
            for o in self.outline:
                if isarc(o):
                    length=arclength(o)
                elif isline(o): # line
                    length=linelength(o)
                else:
                    raise ValueError("bad element in outline list for _calclength()")
                l += length
                ll.append(length)
            self._length = l
            self._lengths=ll
            
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
                c = self.getCenter() # center of the current figure
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
        self._length=0
        self._lengths=[]

        if len(self._elem) < 2:
            ## if one element, the outline is the element
            if len(self._elem) == 1:
                self.outline=deepcopy(self._elem)
                self._calclength()
            ## if there are fewer than one elem, there is nothing to do
            return
        ## we have work to do.  Construct all lines and circle-tangent
        ## lines first, inserting explicit line and arc elements.
        ## Then go back and convert circles to arcs.  Finally compute
        ## total length and length list for use in sample() function.

        
        for i in range(len(self._elem)):
            if i == 0:
                e0 = self._elem[-1]
            else:
                e0 = self._elem[i-1]
            e1 = self._elem[i]
            if i == len(self._elem)-1: #last item
                e2 = self._elem[0]
            else:
                e2 = self._elem[i+1]
            ## work through element types
            if isvect(e0): #e0 is a point
                _fromPointAdd(e0,e1,e2)
            elif isline(e0):
                _fromPointAdd(e0[1],e1,e2)
            elif isarc(e0) and not iscircle(e0):
                p0 = samplearc(e0,1.0)
                _fromPointAdd(p0,e1,e2)
            elif iscircle(e0):
                c = self.getCenter()
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
                if False and end < start:
                    raise NotImplementedError("not handling non-convex case for circle-tangent lines in poly outline yet")
                newarc = arc(e1[0],e1[1][0],start,end)
                self.outline[i]=newarc
                
        _calclength()
        
    def sample(self,u):
        if len(self.outline) == 0:
            return vect(0,0)
        dist = (u % 1.0) * self._length
        d = 0
        for i in range(len(self.outline)):
            l=self._lengths[i]
            if dist <= d+l:
                uu = (d+l-dist)/l
                e = self.outline[i]
                if isline(e):
                    return sampleline(e,1.0-uu)
                else:
                    return samplearc(e,1.0-uu)
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

    
