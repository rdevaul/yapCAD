## yapCAD boolen operation support

from yapcad.geom import *
from yapcad.poly import *

combineDebugGL=[]

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
        try :
            inter = intersectXY(g1.geom(),g2.geom(),params=True)
        except ValueError:
            print("had a problem intersecting following geometries:")
            print("g1.geom(): ",g1.geom())
            print("g2.geom(): ",g2.geom())
            raise

        if inter != False and (inter[0] == False or inter[1] == False):
            raise ValueError('bad intersection list: ',inter)
        # Given parameter values that define two potential sub-arcs of
        # a figure, determine which sub-arc is valid by looking for
        # intersections between these parameter values. In the
        # condition where u2 is smaller than u1, it's possible that
        # this represents an arc in the counter-clockwise direction
        # between u2 and u1.  It also could represent a clockwise arc
        # that "wraps around" from u1 to (u2+1).  We will check for
        # values x1, x2 in ilist that are u2 < x1 < u1 and u1 < x2 <
        # (u2+1.0).  The existance of x1 or x2 rules out the
        # corresponding arc.  If neither x1 nor x2 exists, we bias
        # towards the counter-clockwise arc.

        def between(u1,u2,ilist):
            if len(ilist) < 2:
                raise ValueError('bad ilist')
            
            if u1 > u2 :
                if len(ilist) == 2:
                    return u1,u2+1.0,False
                x1s = list(filter(lambda x: x > u2 and x < u1,ilist))
                x2s = list(filter(lambda x: x > u1 or x < u2,ilist))
                l1 = len(x1s)
                l2 = len(x2s)
                #print("u1: ",u1," u2: ",u2," ilist: ",ilist," x1s: ",x1s," x2s: ",x2s)
                if l1 > 0 and l2 > 0:
                    print('WARNING: intersections on both sides')
                elif l1 > l2:
                    print("AA")
                    if True or self._type == 'union':
                        return u1,u2+1.0,False
                    else:
                        return u2,u1,True
                else:
                    print("A-")
                    return u2,u1,True

            else:
                if len(ilist) == 2:
                    return u1,u2,False
                x1s = list(filter(lambda x: x > u1 and x < u2,ilist))
                x2s = list(filter(lambda x: x > u2 or x < u1,ilist))
                l1 = len(x1s)
                l2 = len(x2s)
                #print("u1: ",u1," u2: ",u2," ilist: ",ilist," x1s: ",x1s," x2s: ",x2s)
                if l1 > 0 and l2 > 0:
                    print('WARNING: intersections on both sides')
                elif l1 > l2:
                    print("BB")

                    if True or self._type == 'union':
                        return u2,u1+1,False
                    else:
                        return u1,u2,False
                else:
                    print("B-")
                    return u1,u2,False
            

        ## utility to perform combination on one "segment"
        def cmbin(g1,g2,itr):
            g1s = itr[0][0]
            g1e = itr[0][1]
            g2s = itr[1][0]
            g2e = itr[1][1]

            seg = []
            ZLEN1=close(g1s,g1e)
            ZLEN2=close(g2s,g2e)

            g1reverse=False
            g2reverse=False
            
            if True or self._type == 'difference':
                if g1e < g1s:
                    g1e+=1.0
                #g1s,g1e,g1reverse = between(g1s,g1e,inter[0])
                g2s,g2e,g2reverse = between(g2s,g2e,inter[1])
                g2reverse = False
            else:
                if g1e < g1s:
                    g1e+=1.0
                if g2e < g2s:
                    g2e += 1.0

            p1=g1.sample(((g1s+g1e)/2)%1.0)

            p1inside=0
            for i in range(5):
                u = (i+1)/6.0
                p = g1.sample((u*g1e+(1.0-u)*g1s)%1.0)
                if g2.isinside(p):
                    p1inside=p1inside+1

            p2inside = 0
            for i in range(5):
                u = (i+1)/6.0
                p = g2.sample((u*g2e+(1.0-u)*g2s)%1.0)
                if g1.isinside(p):
                    p2inside=p2inside+1
                
            if p1inside > 0 and p2inside > 0:
                print("warning: inside test succeeded for both p1s and p2s: ",
                      p1inside," ",p2inside)

            if p1inside == 0 and p2inside == 0:
                print("warning: inside test failed for both p1s and p2s")
                
            p2=g2.sample(((g2s+g2e)/2)%1.0)

            if ZLEN1 and ZLEN2:
                print ('both segments zero length')
                return []
            elif ZLEN2 and not ZLEN1:
                print ('zero length segment 2')
                if self._type=='union':
                    return g1.segment(g1s,g1e)
                elif self._type=='difference':
                    return 
                else: #intersection
                    return []
            elif ZLEN1 and not ZLEN2:
                print ('zero length segment 1')
                if self._type=='union':
                    if g2e < g2s:
                        g2e += 1.0
                    return g2.segment(g2s,g2e)
                else: # difference or intersection
                    return []
            
            if self._type == 'union':
                #if g2.isinside(p1):
                if p1inside > p2inside:
                    # if g2e < g2s:
                    #     g2e += 1.0
                    seg += g2.segment(g2s,g2e)
                else:
                    seg += g1.segment(g1s,g1e)                    
            elif self._type == 'intersection':
                #if g2.isinside(p1):
                if p1inside > p2inside:
                    seg += g1.segment(g1s,g1e)
                else:
                    # if g2e < g2s:
                    #     g2e += 1.0
                    #seg += g2.segment(g2s,g2e,reverse=g2reverse)
                    seg += g2.segment(g2s,g2e)
            elif self._type == 'difference':
                s = []
                #if g2.isinside(p1):
                if p1inside > p2inside:
                    pass
                else:
                    # print("rsort: ",vstr(inter))
                    seg += g1.segment(g1s,g1e)
                    # print("g2s: ",g2s," g2e: ",g2e," g2reverse: ",g2reverse)
                    s = g2.segment(g2s,g2e,reverse=g2reverse)
                    s = reverseGeomList(s)
                    seg += s
                if len(inter[0]) > 2:
                    combineDebugGL.append(s)
                    # print("seg: ",vstr(seg))

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

        if inter == False: # disjoint, but bounding boxes might be
                           # null or one could be inside the other
            if not bbox1 and bbox2: # g1 is empty, but g2 contains geometry
                if self._type=='union':
                    return g2.geom()
                else:
                    return []
            if bbox1 and not bbox2: # g2 is empty, but g1 isn't
                if self._type=='union' or self._type=='difference':
                    return g1.geom()
            if not bbox1 and not bbox2: # no geometry at all
                return []
            ## OK, no intersection but it is possible that one profile
            ## could be inside the other.  Do fast bounding box checks
            ## before doing intersection-based checking.
            if isinsidebbox(bbox1,bbox2[0]) and isinsidebbox(bbox1,bbox2[1]) \
               and g1.isinside(g2.sample(0.0)): # g2 is inside g1
                ## g2 is inside g1
                if self._type == 'union':
                    return g1.geom()
                elif self._type == 'intersection':
                    return g2.geom()
                else: #difference, g2 is a hole in g1
                    return [g1.geom(),g2.geom()]
            elif isinsidebbox(bbox2,bbox1[0]) and isinsidebbox(bbox2,bbox1[1]) \
                 and g2.isinside(g1.sample(0.0)): # g1 is inside g2
                ## g1 is indside g2
                if self._type == 'union':
                    return g2.geom()
                elif self._type == 'intersection':
                    return g1.geom()
                else: #difference, g2 has eaten g1
                    return []
            else: # g1 and g2 are disjoint
                if self._type == 'union':
                    return [g1.geom(),g2.geom()]
                elif self._type == 'difference':
                    return g1.geom()
                else: #intersection
                    return []
        if len(inter[0]) == 1 and len(inter[1]) == 1:
        ## single point of intersection:
            if self._type == 'union':
                return [g1.geom(), g2.geom()]
            elif self._type == 'difference':
                return g1.geom()
            else: #intersection
                return []
        ## There are two or more points of intersection.
        inter = rsort(inter)
        #print("rsort: ",vstr(inter))

        if len(inter[0]) %2 != 0:
            print("WARNING: odd number of intersections (",len(inter[0]),", unpredictable behavior may result")
        r = []
        for i in range(1,len(inter[0])+1):
            r += cmbin(g1,g2,[[inter[0][i-1],
                               inter[0][i%len(inter[0])]],
                              [inter[1][i-1],
                               inter[1][i%len(inter[1])]]])
        return r

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

    def segment(self,u1,u2,reverse=False):
        gl = self.geom()
        if gl == []:
            raise ValueError('empty Boolean, segment not defined')
        return segmentgeomlist(gl,u1,u2,closed=True,reverse=reverse)

    def mirror(self,plane,poly=False):
        b = deepcopy(self)
        b._elem = []
        b._update=True
        for p in self._elem:
            p2 = p.mirror(plane,poly=True)
            b._elem.append(p2)
        if poly:
            return b
        return b.geom()

    def rotate(self,angle,cent=point(0,0,0),axis=point(0,0,1),poly=False):
        b = deepcopy(self)
        b._elem = []
        b._update=True
        for p in self._elem:
            p2 = p.rotate(angle,cent,axis,poly=True)
            b._elem.append(p2)
        if poly:
            return b
        return b.geom()

    def scale(self,sx,sy=False,sz=False,cent=point(0,0),poly=False):
        b = deepcopy(self)
        b._elem = []
        b._update=True
        for p in self._elem:
            p2 = p.scale(sx,sy,sz,cent,poly=True)
            b._elem.append(p2)
        if poly:
            return b
        return b.geom()                        

    def translate(self,delta,poly=False):
        b = deepcopy(self)
        b._elem = []
        b._update=True
        for p in self._elem:
            p2 = p.translate(delta,poly=True)
            b._elem.append(p2)
        if poly:
            return b
        return b.geom()                        
            
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
                self._outline = self.combine_geom(self._elem[0],self._elem[1])
                self._outline = cullZeroLength(self._outline)
                self._update = False
            else:
                raise NotImplementedError("don't know how to do {} yet for {} polygons".format(self._type,len(self._elem)))
        return deepcopy(self._outline)
        
        
