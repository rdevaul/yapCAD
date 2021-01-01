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

def issurface(s,fast=True):
    """
    Check to see if ``s`` is a valid surface.
       ``surface = ['surface',vertices,normals,indices]``, where:
            ``vertices`` is a list of ``yapcad.geom`` points,
            ``normals`` is a list of ``yapcad.geom`` points of the same length as ``vertices``,
            and ``faces`` is the list of faces, which is to say lists of three indices that 
            refer to the vertices of the triangle that represents each face.
    """
    if not isinstance(s,list) or s[0] != 'surface' or len(s) != 4:
        return False
    if fast:
        return True
    else:
        verts=s[1]
        norms=s[2]
        faces=s[3]
        if (not ispoly(verts) or
            not ispoly(norms) or
            len(verts) != len(norms)):
            return False
        l = len(verts)
        if (len(lit(filter(lambda x: not len(x) == 3, faces))) > 0):
            return False
        if (len(list(filter(lambda x: not (isinstance(x,int) or
                                           x < 0 or x >= l),
                            reduce((lambda x,y: x + y),faces)))) > 0):
            return False
        return True

def normfunc(tri):
    v1 = sub(tri[1],tri[0])
    v2 = sub(tri[2],tri[1])
    d = cross(v1,v2)
    n = scale3(d,1.0/mag(d))
    return n,n,n
    
def addTri2Surface(tri,s,check=False,nfunc=normfunc):
    """
    Add triangle ``tri`` (a list of three points) to a surface ``s``,
    returning the updated surface.  *NOTE:* There is no enforcement of
    contiguousness or coplainarity -- this function will add any triangle.
    """

    def addVert(p,n,vrts,nrms):
        for i in range(len(vrts)):
            if vclose(p,vrts[i]):
                return i,vrts,nrms
        vrts.append(p)
        nrms.append(n)
        return len(vrts)-1,vrts,nrms

    if check and (not issurface(s) or not istriangle(tri)):
        raise ValueError(f'bad arguments to addTri2Surface({tri},{s})')
    
    vrts = s[1]
    nrms = s[2]
    faces = s[3]

    n1,n2,n3 = nfunc(tri)
    i1,vrts,nrms = addVert(tri[0],n1,vrts,nrms)
    i2,vrts,nrms = addVert(tri[1],n2,vrts,nrms)
    i3,vrts,nrms = addVert(tri[2],n3,vrts,nrms)
    faces.append([i1,i2,i3])

    return ['surface',vrts,nrms,faces]

def poly2surface(ply,minlen=0.5,minarea=0.0001,checkclosed=False,box=None):
    """
    given an XY-coplanar polygon, return the triangulated surface
    representation of that polygon. If ``checkclosed`` is true, make
    sure ``ply`` is a vaid, closed, XY-coplanar polygon.  if ``box``
    exists, use it as the bounding box.
        ``surface = ['surface',vertices,normals,indices]``, where:
            ``vertices`` is a list of ``yapcad.geom`` points,
            ``normals`` is a list of ``yapcad.geom`` points of the same length as ``vertices``,
            and ``indices`` is the list of indices (three at a time) that represent the
            indices of each of the tiangles that make up the surface
    """
    
    if checkclosed and not ispolygonXY(ply):
        raise ValueError('ply is not an XY coplanar closed polygon')

    if not box:
        box = bbox(ply)

    def triarea(p1,p2,p3):
        v1=sub(p2,p1)
        v2=sub(p3,p2)
        cp = cross(v1,v2)
        return cp[2]/2
    
    def leftTurn(p1,p2,p3):
        v1=sub(p2,p1)
        v2=sub(p3,p2)
        cp = cross(v1,v2)
        return (cp[2] > epsilon)

    outpoint = add(box[1],[1,1,0,1])
    def insideTest(p,poly=ply):
        l = line(p,outpoint)
        pp = intersectSimplePolyXY(l,poly)
        if pp == False:
            return False
        return len(pp) % 2 == 1
    
    ome = 1.0-epsilon
    def selfIntersect(l,bndry):
        for i in range(1,len(bndry)+1):
            uu = lineLineIntersectXY(l,line(bndry[i-1],
                                            bndry[i%len(bndry)]),params=True)
            if (not isinstance(uu,bool) and
                uu[0] > epsilon and uu[0] < ome and
                uu[1] > epsilon and uu[1] < ome):
                return True
        return False


    def trimsharp(bndry):
        if len(bndry) < 3:
            return bndry, False
        print(f"len bndry: {len(bndry)}")
        nodel = True
        for i in range(1,len(bndry)+1):
            delit  = True
            while len(bndry) > 3 and delit:
                p1 = bndry[(i-2)%len(bndry)]
                p2 = bndry[(i-1)%len(bndry)]
                p3 = bndry[i%len(bndry)]
                v1=sub(p2,p1)
                v2=sub(p3,p2)
                cp = cross(v1,v2)
                dp = dot(v1,v2)
                if abs(cp[2]/2) < minarea and dp < 0:
                    del bndry[(i-1)%len(bndry)]
                    nodel = False
                else:
                    delit = False
        print(f"len bndry: {len(bndry)}")
        return bndry, not nodel
    
    def subdivide(i,bndry):
        if len(bndry) < 2 or i < 1 or i > len(bndry):
            raise ValueError('bad index or list for subdivision')
        p0 = bndry[i-1]
        p1 = bndry[i%len(bndry)]
        if dist(p0,p1) < minlen:
            return bndry,False
        np = scale3(add(p0,p1),0.5)
        bndry.insert(i,np)
        return bndry,True
        
    vrts=[]
    nrms=[]
    faces=[]
    surf=['surface',[],[],[]]

    bndry=deepcopy(ply[0:-1]) # last point is redundant

    cnt = 0
    addcnt = 0
    while (len(bndry) > 2):
        cnt+= 1
        if cnt > 1000:
            assert False
        l = len(bndry)
        added = False
        for i in range(l-2):
            testl = line(bndry[i],
                         bndry[i+2])
            testp = scale3(add(bndry[i],
                               bndry[i+2]),0.5)
            area = triarea(bndry[i],
                           bndry[i+1],
                           bndry[i+2])
            if (area > epsilon and
                not selfIntersect(testl, bndry) and
                insideTest(testp)):
                if area > minarea:
                    surf = addTri2Surface(bndry[i:i+3],surf,
                                          nfunc=lambda x: ([0,0,1],
                                                           [0,0,1],
                                                           [0,0,1]))
                addcnt+=1
                
                del bndry[i+1]
                added = True
                break

        bndry, trimmed = trimsharp(bndry)

        if not added and not trimmed and len(bndry) > 2:
            div = False
            skip = False
            for i in range(1,len(bndry)+1):
                if skip:
                    print(f"i : {i} len(bndry): {len(bndry)}")
                    skip=False
                else:
                    bndry,diddiv = subdivide(i,bndry)
                    div = div or diddiv
                    if diddiv:
                        print(f"did div at {i}")
                        skip = True

                
            if not div:
                print("incomplete triangulation")
                print(f"cnt: {cnt}, addcnt: {addcnt}")
                print(f"bndry len: {len(bndry)}, surface triangles: {len(surf[3])}")
                print(f"bndry: {vstr(bndry)}")
                return surf,bndry
                #raise ValueError('unable to finish triangulating surface')
        

    return surf,bndry

def surf2lines(surf):
    """
    convert a surface representation to a non-redundant set of lines
    for line-based rendering purposes
    """
    
    drawn = []

    verts = surf[1]
    norms = surf[2]
    faces = surf[3]

    lines = []

    def inds2key(i1,i2):
        if i1 <= i2:
            return f"{i1}-{i2}"
        else:
            return f"{i2}-{i1}"

    def addLine(i1,i2,lines):
        key = inds2key(i1,i2)
        if not key in drawn:
            lines.append(line(verts[i1],
                              verts[i2]))
            drawn.append(key)
        return lines
    
    for f in faces:
        lines = addLine(f[0],f[1],lines)
        lines = addLine(f[1],f[2],lines)
        lines = addLine(f[2],f[0],lines)

    return lines
    

