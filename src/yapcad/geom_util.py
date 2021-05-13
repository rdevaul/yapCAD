## utilty functions companion to yapcad.geom
## Born on 26 December, 2020
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

from yapcad.geom import *
import math
import random

""" 
Less-essential yapCAD utility functions to support the creation and
operations on yapcad.geom figures.

"""

def randomPoints(bbox,numpoints):
    """Given a 3D bounding box and a number of points to generate, 
    return a list of uniformly generated random points within the 
    bounding box"""
    
    points = []
    minx = bbox[0][0]
    maxx = bbox[1][0]
    miny = bbox[0][1]
    maxy = bbox[1][1]
    minz = bbox[0][2]
    maxz = bbox[1][2]
    rangex = maxx-minx
    rangey = maxy-miny
    rangez = maxz-minz
    for i in range(numpoints):
        points.append(point(random.random()*rangex+minx,
                            random.random()*rangey+miny,
                            random.random()*rangez+minz))
    return points

def randomCenterInBox(bbox,r):
    """given a bounding box and a radius, generate a random center point
    such that a circle with the specified radius and center point falls
    completely inside the box

    """
    
    minx = bbox[0][0]+r
    maxx = bbox[1][0]-r
    miny = bbox[0][1]+r
    maxy = bbox[1][1]-r
    minz = bbox[0][2]+r
    maxz = bbox[1][2]-r
    rangex = maxx-minx
    rangey = maxy-miny
    rangez = maxz-minz
    x = point(random.random()*rangex+minx,
              random.random()*rangey+miny,
              random.random()*rangez+minz)
    return x


def randomArc(bbox,minr=0.0,maxr=10.0,circle=False):
    """given a bounding box and a minimum and maximum radius, generate an
    arc with random center and radius that falls within the bounding
    box.  If ``circle==False``, use randomly generated start and end
    angles, otherwise generate only full circles

    """
    radr = maxr-minr
    r = random.random()*radr+minr
    start = 0
    end = 360
    if not circle:
        start = random.random()*360.0
        end = start + random.random()*360.0

    x = randomCenterInBox(bbox,r)
    return arc(x,r,start,end)

def randomPoly(bbox,numpoints=10,minr = 1.0,maxr = 10.0):
    """given a bounding box, a number of vertices, and a minimum and
    maximum radius, generate a simple polygon that completely lies
    inside the bounding box, whose vertices are evenly spaced by angle
    around a center point, with randomly chosen radius between
    ``minr`` and ``maxr``
    """

    angles = []
    rads = []
    ang = 0.0
    rr = maxr-minr
    for i in range(numpoints):
        a = random.random()
        angles.append(ang)
        ang = ang + a
        # print("a: ",a," ang: ",ang)
        rads.append(random.random()*rr+minr)

    sf = pi2/ang

    points = []
    x = randomCenterInBox(bbox,maxr)
    
    for i in range(numpoints):
        p = [cos(angles[i]*sf)*rads[i],
             sin(angles[i]*sf)*rads[i],0,1]
        points.append(add(p,x))

    return points + [ points[0] ]

def makeLineSpiral(center, turnRad, # radius after one full turn
                   turns, # number of turns
                   dstep = 10.0, # sampling resolution in degrees
                   minlen = 0.25): # minimum distance between points
    """given a center point, the increase in radius per turn, the number
    of turns, and the angular resolution of the approximation,
    generate a yapcqad.geom poly approximation of the spiral

    """

    # make a spiral of points
    spiral = []
    rstep = turnRad*dstep/360.0

    def addpoint(p,lastpoint):
        if not lastpoint or dist(p,lastpoint) > minlen:
            lastpoint = p
            spiral.append(p)
        return lastpoint

    lastpoint = None
    for i in range(round(360*turns/dstep)):
        ang = i * dstep*pi2/360.0
        r = i * rstep
        p = add(center,
                point(math.cos(ang)*r,math.sin(ang)*r))
        lastpoint = addpoint(p,lastpoint)
    return spiral

def makeArcSpiral(center, turnRad, # radius after one full turn
                  turns, # number of turns
                  dstep = 45): # sampling resolution in degrees
    """given a center point, the increase in radius per turn, the number
    of turns, and the angular resolution of the approximation,
    generate an approximation of the spiral using circular arcs as
    segments instead of straightlines. Return a yapcqad.geom geomlist
    approximation of the spiral

    """
    spPoints = makeLineSpiral(center,turnRad,turns,dstep)

    arcs=[]
    for i in range(1,len(spPoints)):
        p0 = spPoints[i-1]
        p1 = spPoints[i]
        direct = sub(p1,p0)
        orth = scale3(orthoXY(direct),-1.1)
        mid = scale3(add(p0,p1),0.5)
        cent = add(mid,orth)
        r= dist(cent,p0)
        r2 = dist(cent,p1)
        assert close(r,r2)
        circ = arc(cent,r)
        u1 = unsamplearc(circ,p0)
        u2 = unsamplearc(circ,p1)
        a = segmentarc(circ,u1,u2)
        arcs.append(a)
    return arcs

def polarSampleArc(c,ang,inside=True):
    """given an arc ``c`` , sample that arc at ``ang`` degrees from the
    start of the arc, proceeding clockwise.  If ``samplereverse`` is
    true, sample the arc at ``-ang`` degrees from the end, proceeding
    counterclockwise.  If ``c`` is not a complete circle and the
    specified ``ang`` falls outside the closed [start,end] interval
    and ``inside == True``, return ``False``, otherwise return the
    sampled value.

    """
    p = c[0]
    r = c[1][0]
    start=c[1][1]
    end=c[1][2]
    samplereverse = (c[1][3] == -2)
    circle = (start == 0 and end == 360)

    if not circle:
        start = start % 360.0
        end = end % 360.0
        if end < start:
            end += 360.0
        if inside and ang > (end - start):
            return False
    
    if samplereverse:
        start=end
        ang *= -1

    srad = (start+ang)*pi2/360.0
    q = scale3(vect(cos(srad),sin(srad)),r)
    return add(p,q)


def triarea(p1,p2,p3):
    """
    utility function to return the area of a triangle
    """
    v1=sub(p2,p1)
    v2=sub(p3,p1)
    cp = cross(v1,v2)
    return mag(cp)/2

def geomlist2poly(gl,minang=5.0,minlen=0.25,checkcont=False):

    """
    convert a continuous geometry list that contains only lines and
    arcs into a poly representation that contains only points. Arcs
    are sampled at a minimum resolution of ``minang`` degrees, and
    resulting points are guaranteed to be no closer than ``minlen``
    apart.  Returns a valid poly.

    """

    if checkcont and not iscontinuousgeomlist(gl):
        raise ValueError('non-continuous geometry list passed to geomlist2poly')
    
    ply = []
    lastpoint = None
    def addpoint(p,lastpoint):
        if not lastpoint or dist(p,lastpoint) > minlen:
            lastpoint = p
            ply.append(p)
        return lastpoint

            
    for e in gl:
        if ispoint(e):
            pass  # ignore points in source geometry list
        elif isline(e):
            p0 = e[0]
            p1 = e[1]
            lastpoint = addpoint(p0,lastpoint)
            lastpoint = addpoint(p1,lastpoint)
        elif isarc(e):
            firstpoint = None
            for ang in range(0,360,round(minang)):
                p = polarSampleArc(e,float(ang),inside=True)
                if not p:
                    break
                else:
                    if not firstpoint:
                        firstpoint = p
                    lastpoint = addpoint(p,lastpoint)
            if iscircle(e):
                lastpoint = addpoint(firstpoint,lastpoint)
        elif ispoly(e):
            for p in e:
                lastpoint = addpoint(p,lastpoint)
        elif isgeomlist(e):
            pts = geomlist2poly(e,minang,minlen,checkcont)
            for p in pts:
                lastpoint = addpoint(p,lastpoint)
        else:
            raise ValueError(f'bad object in list passed to geomlist2poly: {e}')
    return ply



def combineglist(g1,g2,operation):
    """function to perform set operations on geometry lists.

    ``g1`` and ``g2`` are closed figure geometry lists.
    ``operation`` is one of ``["union","intersection","difference"]``
    
    result is a potentially zero-length geometry list representing the
    result, which may or may not be simply connected, but should
    always be a closed figure.

    """

    def poly2lines(pol):
        l = []
        for i in range(1,len(pol)):
            l.append([pol[i-1],pol[i]])
        return l

    if ispoly(g1):
        g1 = poly2lines(g1)

    if ispoly(g2):
        g2 = poly2lines(g2)
        
    # if not (isclosedgeomlist(g1) and isclosedgeomlist(g2)):
    #     raise ValueError("bad arguments passed to combineglist")

    if not operation in ["union","intersection","difference"]:
        raise ValueError("bad operation specified for combineglist")
    
    bbox1 = bbox(g1)
    bbox2 = bbox(g2)
    try :
        inter = intersectXY(g1,g2,params=True)
    except ValueError:
        print("had a problem intersecting following geometries:")
        print("g1: ",g1)
        print("g2: ",g2)
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

            if l1 > l2:
                print("AA")
                if True or operation == 'union':
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
                
            if l1 > l2:
                print("BB")

                if True or operation == 'union':
                    return u2,u1+1,False
                else:
                    return u1,u2,False
            else:
                print("B-")
                return u1,u2,False
            

    ## utility to perform combination on one "segment"
    def cmbin(g1,g2,itr):
        if not isgeomlist(g1):
            g1 = [ g1 ]

        if not isgeomlist(g2):
            g2 = [ g2 ]
        
        g1s = itr[0][0]
        g1e = itr[0][1]
        g2s = itr[1][0]
        g2e = itr[1][1]
        
        seg = []
        ZLEN1=close(g1s,g1e)
        ZLEN2=close(g2s,g2e)

        g1reverse=False
        g2reverse=False
            
        if True or operation == 'difference':
            if g1e < g1s:
                g1e+=1.0
            #g1s,g1e,g1reverse = between(g1s,g1e,inter[0])
            #import pdb ; pdb.set_trace()
            g2s,g2e,g2reverse = between(g2s,g2e,inter[1])
            g2reverse = False
        else:
            if g1e < g1s:
                g1e+=1.0
            if g2e < g2s:
                g2e += 1.0

        #p1=sample(g1,((g1s+g1e)/2)%1.0)

        p1inside=0
        for i in range(5):
            u = (i+1)/6.0
            p = sample(g1,(u*g1e+(1.0-u)*g1s)%1.0)
            if isinsideXY(g2,p):
                p1inside=p1inside+1

        p2inside = 0
        for i in range(5):
            u = (i+1)/6.0
            p = sample(g2,(u*g2e+(1.0-u)*g2s)%1.0)
            if isinsideXY(g1,p):
                p2inside=p2inside+1
                
        if p1inside > 0 and p2inside > 0:
            print("warning: inside test succeeded for both p1s and p2s: ",
                  p1inside," ",p2inside)

        if p1inside == 0 and p2inside == 0:
            print("warning: inside test failed for both p1s and p2s")
                
        #p2=sample(g2,((g2s+g2e)/2)%1.0)

        if ZLEN1 and ZLEN2:
            print ('both segments zero length')
            return []
        elif ZLEN2 and not ZLEN1:
            print ('zero length segment 2')
            if operation=='union':
                return segmentgeomlist(g1,g1s,g1e,closed=True)
            elif operation=='difference':
                return []
            else: #intersection
                return []
        elif ZLEN1 and not ZLEN2:
            print ('zero length segment 1')
            if operation=='union':
                if g2e < g2s:
                    g2e += 1.0
                return segmentgeomlist(g2,g2s,g2e,closed=True)
            else: # difference or intersection
                return []
            
        if operation == 'union':
            #if isinsideXY(g2,p1):
            if p1inside > p2inside:
                # if g2e < g2s:
                #     g2e += 1.0
                seg += segmentgeomlist(g2,g2s,g2e,closed=True)
            else:
                seg += segmentgeomlist(g1,g1s,g1e,closed=True)
        elif operation == 'intersection':
            #if isinsideXY(g2,p1):
            if p1inside > p2inside:
                seg += segmentgeomlist(g1,g1s,g1e,closed=True)
            else:
                # if g2e < g2s:
                #     g2e += 1.0
                #seg += segmentgeomlist(g2,g2s,g2e,reverse=g2reverse)
                seg += segmentgeomlist(g2,g2s,g2e,closed=True)
        elif operation == 'difference':
            s = []
            #if isinsideXY(g2,p1):
            if p1inside > p2inside:
                pass
            else:
                # print("rsort: ",vstr(inter))
                seg += segmentgeomlist(g1,g1s,g1e,closed=True)
                # print("g2s: ",g2s," g2e: ",g2e," g2reverse: ",g2reverse)
                s = segmentgeomlist(g2,g2s,g2e,closed=True,reverse=g2reverse)
                s = reverseGeomList(s)
                seg += s
            if len(inter[0]) > 2:
                pass
                #combineDebugGL.append(s)
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
            if operation=='union':
                return g2
            else:
                return []
        if bbox1 and not bbox2: # g2 is empty, but g1 isn't
            if operation=='union' or operation=='difference':
                return g1
        if not bbox1 and not bbox2: # no geometry at all
            return []
        ## OK, no intersection but it is possible that one profile
        ## could be inside the other.  Do fast bounding box checks
        ## before doing intersection-based checking.
        if isinsidebbox(bbox1,bbox2[0]) and isinsidebbox(bbox1,bbox2[1]) \
           and isinsideXY(g1,sample(g2,0.0)): # g2 is inside g1
            ## g2 is inside g1
            if operation == 'union':
                return g1
            elif operation == 'intersection':
                return g2
            else: #difference, g2 is a hole in g1
                return [g1,g2]
        elif isinsidebbox(bbox2,bbox1[0]) and isinsidebbox(bbox2,bbox1[1]) \
             and isinsideXY(g2,sample(g1,0.0)): # g1 is inside g2
            ## g1 is indside g2
            if operation == 'union':
                return g2
            elif operation == 'intersection':
                return g1
            else: #difference, g2 has eaten g1
                return []
        else: # g1 and g2 are disjoint
            if operation == 'union':
                return [g1,g2]
            elif operation == 'difference':
                return g1
            else: #intersection
                return []
    if len(inter[0]) == 1 and len(inter[1]) == 1:
    ## single point of intersection:
        if operation == 'union':
            return [g1, g2]
        elif operation == 'difference':
            return g1
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
