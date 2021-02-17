## geom3d_util, additional 3D geometry support for yapCAD
## started on Mon Feb 15 20:23:57 PST 2021 Richard W. DeVaul

from yapcad.geom import *
from yapcad.geom_util import *
from yapcad.xform import *
from yapcad.geom3d import *

import math

"""
==================================================
Utility functions to support 3D geometry in yapCAD
==================================================

This module is mostly a collection of 
parametric soids and surfaces, and supporting functions.

"""


def sphere2cartesian(lat,lon,rad):
    """
    Utility function to convert spherical polar coordinates to
    cartesian coordinates for a sphere centered at the origin.
        ``lat`` -- latitude
        ``lon`` -- longitude
        ``rad`` -- sphere radius

    returns a ``yapcad.geom`` point
    """
    if lat == 90:
        return [0,0,rad,1]
    elif lat == -90:
        return [0,0,-rad,1]
    else:
        latr = (((lat+90)%180)-90)*pi2/360.0
        lonr = (lon%360)*pi2/360.0
    
        smallrad = math.cos(latr)*rad
        z = math.sin(latr)*rad
        x = math.cos(lonr)*smallrad
        y = math.sin(lonr)*smallrad
        return [x,y,z,1]

## icosohedron-specific function, generate intial geometry
def makeIcoPoints(center,radius):
    """
    Procedure to generate the verticies of an icosohedron with the specified
    ``center`` and ``radius``.
    """
    points = []
    normals = []
    p = sphere2cartesian(90,0,radius)
    n = scale3(p,1.0/mag(p))
    n[3] = 0.0
    points.append(p)
    normals.append(n)
    for i in range(10):
        sgn = 1
        if i%2 == 0:
            sgn =-1
        lat = math.atan(0.5)*360.0*sgn/pi2
        lon = i*36.0
        p = sphere2cartesian(lat,lon,radius)
        n = scale3(p,1.0/mag(p))
        n[3] = 0.0
        points.append(p)
        normals.append(n)

    p = sphere2cartesian(-90,0,radius)
    n = scale3(p,1.0/mag(p))
    n[3] = 0.0
    points.append(p)
    normals.append(n)
    return list(map( lambda x: add(x,center),points)),normals

# face indices for icosahedron
icaIndices = [ [1,11,3],[3,11,5],[5,11,7],[7,11,9],[9,11,1],
               [2,1,3],[2,3,4],[4,3,5],[4,5,6],[6,5,7],[6,7,8],[8,7,9],[8,9,10],[10,9,1],[10,1,2],
               [0,2,4],[0,4,6],[0,6,8],[0,8,10],[0,10,2] ]


def addVertex(nv,nn,verts,normals):
    """
    Utility function that takes a vertex and associated normal and a
    list of corresponding vertices and normals, and returns an index
    that corresponds to the given vertex/normal.  If the vertex
    doesn't exist in the list, the lists are updated to include it and
    the corresponding normal.

    returns the index, and the (potentiall updated) lists
    """
    for i in range(len(verts)):
        if vclose(nv,verts[i]):
            return i,verts,normals
    verts.append(nv)
    normals.append(nn)
    return len(verts)-1,verts,normals

def subdivide(f,verts,normals,rad):
    """
    Given a face (a list of three vertex indices), a list of vertices,
    normals, and a radius, subdivide that face into four new faces and
    update the lists of vertices and normals accordingly.

    return the updated vertex and normal lists, and a list of the four
    new faces.
    """

    ind1 = f[0]
    ind2 = f[1]
    ind3 = f[2]
    v1 = verts[ind1]
    v2 = verts[ind2]
    v3 = verts[ind3]
    n1 = normals[ind1]
    n2 = normals[ind2]
    n3 = normals[ind3]
    va = add(v1,v2)
    vb = add(v2,v3)
    vc = add(v3,v1)
    ma = rad/mag(va)
    mb = rad/mag(vb)
    mc = rad/mag(vc)
    va = scale3(va,ma)
    vb = scale3(vb,mb)
    vc = scale3(vc,mc)
    
    na = add4(n1,n2)
    na = scale4(na,1.0/mag(na))
    nb = add4(n2,n3)
    nb = scale4(nb,1.0/mag(nb))
    nc = add4(n3,n1)
    nc = scale4(nc,1.0/mag(nc))

    inda,verts,normals = addVertex(va,na,verts,normals)
    indb,verts,normals = addVertex(vb,nb,verts,normals)
    indc,verts,normals = addVertex(vc,nc,verts,normals)

    f1 = [ind1,inda,indc]
    f2 = [inda,ind2,indb]
    f3 = [indb,ind3,indc]
    f4 = [inda,indb,indc]

    return verts,normals, [f1,f2,f3,f4]

# make the sphere, return a surface representation
def sphereSurface(diameter,center=point(0,0,0),depth=2):
    rad = diameter/2
    verts,normals = makeIcoPoints(center,rad)
    faces = icaIndices
    
    for i in range(depth):
        ff = []
        for f in faces:
            verts, norms, newfaces = subdivide(f,verts,normals,rad)
            ff+=newfaces
        faces = ff
                      
    
    return ['surface',verts,normals,faces,[],[]]

# make sphere, return solid representation
def sphere(diameter,center=point(0,0,0),depth=2):
    call = f"yapcad.geom3d_util.sphere({diameter},{center},{depth})"
    return solid( [ sphereSurface(diameter,center,depth)],
                  [],['procedure',call] )
                    

def rectangularPlane(length,width,center=point(0,0,0)):
    """ return a rectangular surface with the normals oriented in the
    positive z direction """
    c = center
    l = point(length,0,0)
    w = point(0,width,0)
    p0 = add(c,(scale3(add(l,w),-0.5)))
    p1 = add(p0,l)
    p2 = add(p1,w)
    p3 = add(p0,w)
    n = vect(0,0,1,0)
    
    surf = surface( [p0,p1,p2,p3],[n,n,n,n],
                    [[0,1,2],
                     [2,3,0]])

    return surf

# make a rectangular prism from six surfaces, return a solid
def prism(length,width,height,center=point(0,0,0)):
    """make a rectangular prism solid composed of six independent 
    faces"""
    call = f"yapcad.geom3d_util.prism({length},{width},{height},{center})"

    l2 = length/2
    w2 = width/2
    h2 = height/2
    
    topS = rectangularPlane(length,width,point(0,0,h2))
    bottomS = rotatesurface(topS,180,axis=point(1,0,0))
    frontS = rectangularPlane(length,height,point(0,0,w2))
    frontS = rotatesurface(frontS,90,axis=point(1,0,0))
    backS = rotatesurface(frontS,180)
    rightS = rectangularPlane(height,width,point(0,0,l2))
    rightS = rotatesurface(rightS,90,axis=point(0,1,0))
    leftS = rotatesurface(rightS,180)

    sol = solid([topS,bottomS,frontS,backS,rightS,leftS],
                [],
                ['procedure',call])
    return sol

def circleSurface(center,radius,angr=10,zup=True):
    """make a circular surface centered at ``center`` lying in the XY
    plane with normals pointing in the positive z direction if ``zup
    == True``, negative z otherwise"""

    if angr < 1 or angr > 45:
        raise ValueError('angular resolution must be between 1 and 45 degrees')
    
    samples = round(360.0/angr)
    angr = 360.0/samples
    basep=[center]
    
    for i in range(samples):
        theta = i*angr*pi2/360.0
        pp = [math.cos(theta)*radius,math.sin(theta)*radius,0.0,1.0]
        pp = add(pp,center)
        basep.append(pp)       

    
    basef=[]

    rng = range(1,len(basep))
    if not zup:
        rng = reversed(rng)

    ll = len(basep)-1
    for i in rng:
        face = []
        if zup:
            face = [0,i,1+i%ll]
        else:
            face = [0,1+i%ll,i]
            
        basef.append(face)

    z=-1
    if zup:
        z=1
    n = vect(0,0,z,0)
    basen= [ n ] * len(basep)

    return surface(basep,basen,basef)

def conic(baser,topr,height, center=point(0,0,0),angr=10):

    """Make a conic frustum splid, center is center of first 'base'
    circle, main axis aligns with positive z.  This function can be
    used to make a cone, a conic frustum, or a cylinder, depending on
    the parameters.

         ``baser`` is the radius of the base, must be greater than
         zero (epsilon).

         ``topr`` is the radius of the top, may be zero or
         positive. If ``0 <= topr < epsilon``, then the top
         is treated as a single point.

         ``height`` is distance from base to top, must be greater than
         epsilon.
    
         ``center`` is the location of the center of the base.

         ``angr`` is the requested angular resolution in degrees for
         sampling circles.  Actual angular resolution will be
         ``360/round(360/angr)``

    """
    call = f"yapcad.geom3d_util.conic({baser},{topr},{height},{center},{angr})"
    if baser < epsilon:
        raise ValueError('bad base radius for conic')
    base = arc(center,baser)

    toppoint = False
    if topr < 0:
        raise ValueError('bad top radius for conic')
    if topr < epsilon:
        toppoint = True

    if height < epsilon:
        raise ValueError('bad height in conic')

    baseS = circleSurface(center,baser,zup=False)
    baseV = baseS[1]
    ll = len(baseV)
        
    if not toppoint:
        topS = circleSurface(add(center,point(0,0,height)),
                             topr,zup=True)
        topV = topS[1]
        cylV = baseV[1:] + topV[1:]
        ll = ll-1
        baseN = []
        topN = []
        cylF = []
        for i in range(ll):
            p0 = cylV[(i-1)%ll]
            p1 = cylV[(i+1)%ll]
            p2 = cylV[ll+i]

            cylF.append([i,(i+1)%ll,ll+(i+1)%ll])
            cylF.append([i,ll+(i+1)%ll,ll+i])

            pp,n0 = tri2p0n([p0,p1,p2])
            
            baseN.append(n0)
            topN.append(n0)

        cylN = baseN+topN

        cylS = surface(cylV,cylN,cylF)

        return solid([baseS,cylS,topS],[],
                     ['procedure',call])
    else:
        topP = add(center,point(0,0,height))
        conV = [ topP ] + baseV
        ll = len(conV)
        conN = [[0,0,1,0]]
        conF = []

        for i in range(1,ll):
            p0= conV[0]
            p1= conV[(i-1)%ll]
            p2= conV[(i+1)%ll]

            conF.append([0,i,(i+1)%ll])
            pp,n0 = tri2p0n([p0,p1,p2])

            conN.append(n0)

        conS = surface(conV,conN,conF)

        return solid([baseS,conS],[],
                     ['procedure',call])
