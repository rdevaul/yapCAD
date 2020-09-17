## yapCAD boolen operation support

from yapcad.geom import *
from yapcad.poly import *

class Boolean(IntersectGeometry):
    """Boolean operations on Polygons"""

    types = ('union','intersection','difference')
    
    def __init__(self,type='union',polys=[]):
        for p in polys:
            if not isinstance(p,Polygon):
                raise ValueError('non-poly passed to Boolean(): {}'.format(p))
            if not type in self.types:
                raise ValueError('invalid type passed to Boolean(): {}'.format(tpe))
            self._elem=list(polys)
            self._type=type


    def union_geom(self,g1,g2):
        bbox1 = g1.bbox()
        bbox2 = g2.bbox()
        inter = intersectXY(g1.geom(),g2.geom(),params=True)

        if inter == False: # disjoint, but one could be inside the other
            if isinsidebbox(bbox1,bbox2[0]) and isinsidebbox(bbox1,bbox2[1]):
                ## g2 is inside g1
                return g1.geom()
            elif isinsidebbox(bbox2,bbox1[0]) and isinsidebbox(bbox2,bbox1[1]):
                ## g1 is indside g2
                return g2.geom()
            else:
                return [g1.geom(),g2.geom()]
        if len(inter[0]) == 1 and len(inter[1] == 1):
        ## single point of intersection: return a geometry list with both
            return [g1, g2]
        ## There are two or more points of intersection.
        if len(inter[0]) == 2:
            ## for g1, we have two possiblities for the "outside"
            ## segment: u1 --> u2 or u2 --> u1.  Likewise, the same is
            ## true for g2.  If we sample a point (u1+u2)/2 on g1 and
            ## it falls outside g2, then that is the segment to
            ## choose.  Likewise for the other.
            g1s = inter[0][0]
            g1e = inter[0][1]
            g2s = inter[1][0]
            g2e = inter[1][1]

            p=g1.sample((g1s+g1e)/2)
            if g2.inside(p):
                g1s = g1e
                g1e = inter[0][0]

            p=g2.sample((g2s+g2e)/2)
            if g1.inside(p):
                g2s = g1e
                g2e = inter[1][0]

            seg1 = g1.segment(g1s,g1e)
            seg2 = g2.segment(g2s,g2e)
            return seg1 + seg2
            
    def geom(self):
        if self._type == 'union' and len(self._elem)==2:
            return self.union_geom(self._elem[0],self._elem[1])
        else:
            raise NotImplementedError("don't know how to do {} yet for {} polygons".format(self._type,len(self._elem)))
        
        
