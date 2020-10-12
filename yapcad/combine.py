## yapCAD boolen operation support

from yapcad.geom import *
from yapcad.poly import *

class Boolean(IntersectGeometry):
    """Boolean operations on Polygons"""

    types = ('union','intersection','difference')
    
    def __init__(self,type='union',polys=[]):
        for p in polys:
            if not ( isinstance(p,Polygon) or isinstance(p,Boolean) ):
                raise ValueError('non-poly or non-boolean passed to Boolean(): {}'.format(p))
            if not type in self.types:
                raise ValueError('invalid type passed to Boolean(): {}'.format(tpe))
            self._elem=list(polys)
            self._type=type
            self._update=True
            self._outline=[]


    def combine_geom(self,g1,g2):
        bbox1 = g1.bbox()
        bbox2 = g2.bbox()
        inter = intersectXY(g1.geom(),g2.geom(),params=True)

        ## utility to perform combination on one "segment"
        def cmbin(g1,g2,itr):
            g1s = itr[0][0]
            g1e = itr[0][1]
            g2s = itr[1][0]
            g2e = itr[1][1]

            seg = []
            if g1e < g1s:
                g1e+=1.0
            if g2e < g2s:
                g2e+=1.0

            p1=g1.sample(((g1s+g1e)/2)%1.0)
            # seg.append(arc(p1,0.2))
            p2=g2.sample(((g2s+g2e)/2)%1.0)
            # seg.append(arc(p2,0.4))

            ZLEN1=close(g1s,g1e)
            ZLEN2=close(g2s,g2e)

            if ZLEN1 and ZLEN2:
                print ('both segments zero length')
                return []
            elif ZLEN2 and not ZLEN1:
                print ('zero length segment 2')
                if self._type=='union' or self._type=='difference':
                    return g1.segment(g1s,g1e)
                else: #intersection
                    return []
            elif ZLEN1 and not ZLEN2:
                print ('zero length segment 1')
                if self._type=='union':
                    return g2.segment(g2s,g2e)
                else: # difference or intersection
                    return []
            
            if self._type == 'union':
                if g2.isinside(p1):
                    seg += g2.segment(g2s,g2e)
                else:
                    seg += g1.segment(g1s,g1e)                    
            elif self._type == 'intersection':
                if g2.isinside(p1):
                    seg += g1.segment(g1s,g1e)
                else:
                    seg += g2.segment(g2s,g2e)
            elif self._type == 'difference':
                if g2.isinside(p1):
                    pass
                else:
                    seg += g1.segment(g1s,g1e)
                    seg += g2.segment(g2s,g2e)
            return seg

        ## utility function to sort intersections into non-decreasing
        ## order
        def rsort(il):
            nl = []
            rl = []
            rr = []
            for i in range(len(il[0])):
                nl.append([il[0][i],il[1][i]])
            nl.sort(key=lambda x: x[0])
            for i in range(len(nl)):
                rl.append(nl[i][0])
                rr.append(nl[i][1])
            return [rl,rr]

        if inter == False: # disjoint, but one could be inside the other
            if not bbox1 and bbox2:
                if self._type=='union':
                    return g2.geom()
                else:
                    return []
            if bbox1 and not bbox2:
                if self._type=='union' or self._type=='difference':
                    return g1.geom()
            if not bbox1 and not bbox2:
                return []
            if isinsidebbox(bbox1,bbox2[0]) and isinsidebbox(bbox1,bbox2[1]):
                ## g2 is inside g1
                if self._type == 'union':
                    return g1.geom()
                elif self._type == 'intersection':
                    return g2.geom()
                else: #difference, g2 is a hole in g1
                    return [g1.geom(),g2.geom()]
            elif isinsidebbox(bbox2,bbox1[0]) and isinsidebbox(bbox2,bbox1[1]):
                ## g1 is indside g2
                if self._type == 'union':
                    return g2.geom()
                elif self._type == 'intersection':
                    return g1.geom()
                else: #difference, g2 has eaten g1
                    return []
            else: # disjoint
                if self._type == 'union':
                    return [g1.geom(),g2.geom()]
                elif self._type == 'difference':
                    return g1.geom()
                else: #intersection
                    return []
        if len(inter[0]) == 1 and len(inter[1] == 1):
        ## single point of intersection:
            if self._type == 'union':
                return [g1.gom(), g2.geom()]
            elif self._type == 'difference':
                return g1.geom()
            else: #intersection
                return []
        ## There are two or more points of intersection.
        inter = rsort(inter)
        # print("rsort: ",vstr(inter))
        
        if len(inter[0]) % 2 == 0 or True:
            r = []
            for i in range(1,len(inter[0])+1):
                r += cmbin(g1,g2,[[inter[0][i-1],
                                   inter[0][i%len(inter[0])]],
                                  [inter[1][i-1],
                                   inter[1][i%len(inter[1])]]])
            return r
#        else:
#            raise ValueError('odd number of intersections -- bailing')

    def bbox(self):
        return bbox(self.geom())

    def getCenter(self):
        gl = self.geom()
        if gl == []:
            raise ValueError('empty Boolean, no center')
        return center(gl)

    def getLength(self):
        gl = self.geom()
        if gl == []:
            raise ValueError('empty Boolean, no length')
        return length(gl)

    def segment(self,u1,u2):
        gl = self.geom()
        if gl == []:
            raise ValueError('empty Boolean, segment not defined')
        return segmentgeomlist(gl,u1,u2,closed=True)

    def sample(self,u):
        gl = self.geom()
        if gl == []:
            raise ValueError('empty Boolean, sample not defined')
        return sample(gl,u)

    def isinside(self,p):
        gm = self.geom()
        if gm == []:
            raise ValueError('empty Boolean, inside not defined')
        bb = bbox(gm)
        p2 = add([1,1,0,1],bb[1])
        l = line(p,p2)

        pp = intersectGeomListXY(l,gm)
        if pp == False:
            return False
        return len(pp) % 2 == 1

    def grow(self,r):
        if close(r,0.0):
            return
        elif r < 0:
            raise ValueError('negative growth not allowed')

        if self._type in ['union','intersection']:
            for p in self._elem:
                p.grow(r)
            self._update = True
        else:
            raise NotImplementedError("Don't have grow support for {} yet".format(self._type))
        

    def geom(self):
        if self._update:
            if len(self._elem)==2:
                #if len(self._elem)==2 and \
                #self._type == 'union' or self._type == 'intersection':
                self._outline = self.combine_geom(self._elem[0],self._elem[1])
                self._update = False
            else:
                raise NotImplementedError("don't know how to do {} yet for {} polygons".format(self._type,len(self._elem)))
        return deepcopy(self._outline)
        
        
