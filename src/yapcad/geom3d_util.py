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

vertexHash = {}
_vertexHash_owner = None
def addVertex(nv,nn,verts,normals):
    """
    Utility function that takes a vertex and associated normal and a
    list of corresponding vertices and normals, and returns an index
    that corresponds to the given vertex/normal.  If the vertex
    doesn't exist in the list, the lists are updated to include it and
    the corresponding normal.

    returns the index, and the (potentiall updated) lists
    """
    global vertexHash, _vertexHash_owner
    owner_id = id(verts)
    if _vertexHash_owner != owner_id or len(verts) == 0:
        vertexHash = {}
        _vertexHash_owner = owner_id

    found = False
    # Normalize tiny values to avoid "-0.00" != "0.00" hash key mismatch
    x = 0.0 if abs(nv[0]) < epsilon else nv[0]
    y = 0.0 if abs(nv[1]) < epsilon else nv[1]
    z = 0.0 if abs(nv[2]) < epsilon else nv[2]
    vkey = f"{x:.2f}{y:.2f}{z:.2f}"
    if vkey in vertexHash:
        found = True
        inds = vertexHash[vkey]
        valid_inds = []
        for i in inds:
            if i >= len(verts):
                continue
            if vclose(nv,verts[i]):
                return i,verts,normals
            valid_inds.append(i)
        vertexHash[vkey] = valid_inds
            
    verts.append(nv)
    normals.append(nn)
    i = len(verts)-1
    if found:
        vertexHash[vkey] += [ i ]
    else:
        vertexHash[vkey] = [ i ]
    return i,verts,normals

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
    ## subdivision works only when center is origin, so we add
    ## any center offset after we subdivide
    verts,normals = makeIcoPoints(point(0,0,0),rad)
    faces = icaIndices

    for i in range(depth):
        ff = []
        for f in faces:
            verts, norms, newfaces = subdivide(f,verts,normals,rad)
            ff+=newfaces
        faces = ff
                      
    if not vclose(center,point(0,0,0)):
        verts = list(map(lambda x: add(x,center),verts))
    return ['surface',verts,normals,faces,[],[]]

# make sphere, return solid representation
def sphere(diameter,center=point(0,0,0),depth=2):
    call = f"yapcad.geom3d_util.sphere({diameter},center={center},depth={depth})"
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

    surfaces = [topS,bottomS,frontS,backS,rightS,leftS]
    if not vclose(center,[0,0,0,1]):
        surfaces = list(map(lambda x: translatesurface(x,center),surfaces))
    
    sol = solid(surfaces,
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
        # Only use perimeter vertices (baseV[1:]), skip the center point
        conV = [ topP ] + baseV[1:]
        ll = len(conV)
        # Initialize all normals to a default value
        conN = [[0,0,1,0] for _ in range(ll)]
        conF = []

        # ll = apex(1) + perimeter vertices(36) = 37
        # Perimeter vertex indices: 1 to ll-1
        num_perimeter = ll - 1
        for i in range(1, ll):
            p0 = conV[0]  # apex
            # Wrap indices within perimeter range [1, ll-1]
            prev_idx = ((i - 2) % num_perimeter) + 1
            next_idx = (i % num_perimeter) + 1
            p1 = conV[prev_idx]
            p2 = conV[next_idx]

            try:
                pp, n0 = tri2p0n([p0, p1, p2])
            except ValueError:
                # Skip degenerate faces near the apex
                continue

            conF.append([0, i, next_idx])
            conN[i] = n0

        conS = surface(conV,conN,conF)

        return solid([baseS,conS],[],
                     ['procedure',call])

def makeRevolutionSurface(contour,zStart,zEnd,steps,arcSamples=36):
    """
    Take a countour (any function z->y mapped over the interval

     ``zStart`` and ``zEnd`` and produce the surface of revolution
     around the z axis.  Sample ``steps`` contours of the function,
    which in turn are turned into circles sampled `arcSamples`` times.
    """

    sV=[]
    sN=[]
    sF=[]
    zRange = zEnd-zStart
    zD = zRange/steps

    degStep = 360.0/arcSamples
    radStep = pi2/arcSamples

    # Pre-compute cos/sin values to avoid floating point errors at the seam
    # Explicitly ensure that index 0 uses exact values
    angle_cos = []
    angle_sin = []
    for i in range(arcSamples):
        if i == 0:
            angle_cos.append(1.0)
            angle_sin.append(0.0)
        else:
            angle = i * radStep
            angle_cos.append(math.cos(angle))
            angle_sin.append(math.sin(angle))

    # Check if we need pole caps
    r_start = contour(zStart)
    r_end = contour(zEnd)
    need_start_cap = r_start < epsilon * 10
    need_end_cap = r_end < epsilon * 10

    # Add pole vertices if needed
    start_pole_idx = None
    end_pole_idx = None
    if need_start_cap:
        pole_point = [0.0, 0.0, zStart, 1.0]
        pole_normal = [0.0, 0.0, -1.0, 0.0]  # Points downward for bottom pole
        start_pole_idx, sV, sN = addVertex(pole_point, pole_normal, sV, sN)

    if need_end_cap:
        pole_point = [0.0, 0.0, zEnd, 1.0]
        pole_normal = [0.0, 0.0, 1.0, 0.0]  # Points upward for top pole
        end_pole_idx, sV, sN = addVertex(pole_point, pole_normal, sV, sN)

    for i in range(steps):
        z = i*zD+zStart
        r0 = contour(z)
        r1 = contour(z+zD)

        # Handle pole caps
        if i == 0 and need_start_cap:
            # Create triangular faces from pole to first ring
            if r1 < epsilon:
                r1 = epsilon
            for j in range(arcSamples):
                a1_idx = j
                a2_idx = (j+1) % arcSamples

                pp1 = [angle_cos[a1_idx]*r1, angle_sin[a1_idx]*r1, z+zD, 1.0]
                pp2 = [angle_cos[a2_idx]*r1, angle_sin[a2_idx]*r1, z+zD, 1.0]

                try:
                    _, n = tri2p0n([sV[start_pole_idx], pp1, pp2])
                except ValueError:
                    continue

                k1, sV, sN = addVertex(pp1, n, sV, sN)
                k2, sV, sN = addVertex(pp2, n, sV, sN)
                sF.append([start_pole_idx, k1, k2])
            continue

        if i == steps - 1 and need_end_cap:
            # Create triangular faces from last ring to pole
            if r0 < epsilon:
                r0 = epsilon
            for j in range(arcSamples):
                a1_idx = j
                a2_idx = (j+1) % arcSamples

                p1 = [angle_cos[a1_idx]*r0, angle_sin[a1_idx]*r0, z, 1.0]
                p2 = [angle_cos[a2_idx]*r0, angle_sin[a2_idx]*r0, z, 1.0]

                try:
                    _, n = tri2p0n([p1, sV[end_pole_idx], p2])
                except ValueError:
                    continue

                k1, sV, sN = addVertex(p1, n, sV, sN)
                k2, sV, sN = addVertex(p2, n, sV, sN)
                sF.append([k1, end_pole_idx, k2])
            continue

        # Regular quad strips for non-pole sections
        if r0 < epsilon:
            r0 = epsilon
        if r1 < epsilon:
            r1 = epsilon

        for j in range(arcSamples):
            # Use pre-computed values with proper wrapping
            a0_idx = (j-1) % arcSamples
            a1_idx = j
            a2_idx = (j+1) % arcSamples

            p0 = [angle_cos[a0_idx]*r0, angle_sin[a0_idx]*r0, z, 1.0]
            p1 = [angle_cos[a1_idx]*r0, angle_sin[a1_idx]*r0, z, 1.0]
            p2 = [angle_cos[a2_idx]*r0, angle_sin[a2_idx]*r0, z, 1.0]

            pp1 = [angle_cos[a1_idx]*r1, angle_sin[a1_idx]*r1, z+zD, 1.0]
            pp2 = [angle_cos[a2_idx]*r1, angle_sin[a2_idx]*r1, z+zD, 1.0]

            try:
                p,n = tri2p0n([p0,p2,pp1])
            except ValueError:
                # Skip degenerate faces
                continue

            k1,sV,sN = addVertex(p1,n,sV,sN)
            k2,sV,sN = addVertex(p2,n,sV,sN)
            k3,sV,sN = addVertex(pp2,n,sV,sN)
            k4,sV,sN = addVertex(pp1,n,sV,sN)
            sF.append([k1,k2,k3])
            sF.append([k1,k3,k4])
        
    return surface(sV,sN,sF)

def contour(poly,distance,direction, samples,scalefunc= lambda x: (1,1,1)):
    """take a closed polygon and apply a scaling function defined on the
    interval 0,1 that returns a tuple of x,y,z scaling values. For
    each of ``samples`` number of samples, translate in ``direction``
    direction and scale the contour according to ``scalefunc()``,
    producing a surface."""
    if not ispolygon(poly):
        raise ValueError('invalid polygon passed to contour')
    if samples < 2:
        raise ValueError('number of samples must be 2 or greater')
    if not close(mag(direction),1.0):
        raise ValueError('bad direction vector')
    if distance <= epsilon:
        raise ValueError('bad distance passed to contour')

    raise NotImplemented("this function is a work in progress")

    p0 = poly
    p1 = []

    u = 0.0
    scle = scalefunc(u)
    sctx = Scale(scle[0],scle[1],scle[2])
    ply = list(map(lambda p: sctx.mul(p),poly))
    surf = poly2surface(ply)
    s1 = reversesurface(surf)

    u = 1.0
    scle = scalefunc(u)
    sctx = Scale(scle[0],scle[1],scle[2])
    ply2 = list(map(lambda p:  sctx.mul(p),poly))
    surf2 = poly2surface(ply2)
    s2 = translatesurface(surf2,scale4(direction,distance))
    
    vrts = surf[1]
    nrms = surf[2]
    facs = surf[3]
    bndr = surf[4]

    vrts1 = list(map(lambda i: vrts[i],bndr))
    for ii in range(1,samples):
        u = ii/samples
        
        stripF = []
        #vrts2 = list(map(lambda i:
        for i in range(len(bndr)):
            j0 = bndry1[(i-1)%len(bndry1)]
            j1 = bndry1[i]
            j2 = bndry1[(i+1)%len(bndry1)]
            j3 = j2+len(s2[1])
            j4 = j1+len(s2[1])
            p0 = stripV[j0]
            p1 = stripV[j2]
            p2 = stripV[j3]
            try:
                pp,n0 = tri2p0n([p0,p1,p2])
            except ValueError:
                # bad face, skip
                continue
            stripN[j1]=n0
            stripN[j4]=n0
            stripF.append([j1,j2,j3])
            stripF.append([j1,j3,j4])
            
            


def extrude(surf,distance,direction=vect(0,0,1,0)):

    """ Take a surface and extrude it in the specified direction to
    create a solid.  Return the solid. """
    call = f"yapcad.geom3d_util.extrude({surf},{distance},{direction})"

    if not issurface(surf):
        raise ValueError('invalid surface passed to extrude')

    if distance <= epsilon:
        raise ValueError('bad distance passed to extrude')

    s1 = translatesurface(surf,scale4(direction,distance))
    s2 = reversesurface(surf)

    loops = []
    if s2[4]:
        loops.append(list(s2[4]))
    loops.extend([list(loop) for loop in s2[5] if loop])

    stripV = s2[1] + s1[1]  # vertices for the edge strips
    stripN = [vect(0, 0, 1, 0)] * len(stripV)  # placeholder normals
    stripF: list[list[int]] = []
    offset = len(s2[1])

    def _project(idx):
        point3 = stripV[idx]
        return point3[0], point3[1]

    if surf[3]:
        try:
            tri = surf[3][0]
            basis = tri2p0n([surf[1][tri[0]], surf[1][tri[1]], surf[1][tri[2]]], basis=True)
            if basis:
                forward = basis[2]

                def _project(idx):  # type: ignore[redefinition]
                    vec = forward.mul(stripV[idx])
                    return vec[0], vec[1]
        except Exception:  # pragma: no cover
            pass

    def loop_area(loop):
        if len(loop) < 3:
            return 0.0
        area = 0.0
        for i in range(len(loop)):
            x0, y0 = _project(loop[i])
            x1, y1 = _project(loop[(i + 1) % len(loop)])
            area += x0 * y1 - x1 * y0
        return area / 2.0

    if loops:
        outer_area = loop_area(loops[0])
        if outer_area < 0:
            loops[0].reverse()
            outer_area = -outer_area
        outer_sign = 1 if outer_area >= 0 else -1
        for loop in loops[1:]:
            if loop_area(loop) * outer_sign > 0:
                loop.reverse()

    for bndry in loops:
        if len(bndry) < 2:
            continue
        for i in range(len(bndry)):
            j0 = bndry[(i - 1) % len(bndry)]
            j1 = bndry[i]
            j2 = bndry[(i + 1) % len(bndry)]
            j3 = j2 + offset
            j4 = j1 + offset
            p0 = stripV[j0]
            p1 = stripV[j2]
            p2 = stripV[j3]
            try:
                pp, n0 = tri2p0n([p0, p1, p2])
            except ValueError:
                continue
            stripN[j1] = n0
            stripN[j4] = n0
            stripF.append([j1, j2, j3])
            stripF.append([j1, j3, j4])

    #import pdb ; pdb.set_trace()
    strip = surface(stripV,stripN,stripF)

    return solid([s2,strip,s1],
                 [],
                 ['procedure',call])



def _loft_surface(lower_loop, upper_loop, invert=False):
    """Create a surface connecting two loops.

    Args:
        lower_loop: List of points forming the lower loop (must be open, not closed)
        upper_loop: List of points forming the upper loop (must be open, not closed)
        invert: If True, reverse the face winding order

    Returns:
        A yapCAD surface connecting the two loops with triangle strips
    """

    if not lower_loop or not upper_loop:
        raise ValueError('invalid loops passed to loft surface')
    lower = [point(p) for p in lower_loop]
    upper = [point(p) for p in upper_loop]
    if len(lower) != len(upper):
        raise ValueError('loop length mismatch in loft surface')

    vertices = lower + upper
    normals = [[0, 0, 1, 0] for _ in vertices]
    faces = []
    count = len(lower)

    for idx in range(count):
        j0 = idx
        j1 = (idx + 1) % count
        j2 = j1 + count
        j3 = idx + count

        tri1 = [j0, j1, j2]
        tri2 = [j0, j2, j3]
        if invert:
            tri1 = list(reversed(tri1))
            tri2 = list(reversed(tri2))

        try:
            _, normal = tri2p0n([vertices[tri1[0]],
                                 vertices[tri1[1]],
                                 vertices[tri1[2]]])
        except ValueError:
            continue

        for vid in {tri1[0], tri1[1], tri1[2], tri2[0], tri2[1], tri2[2]}:
            normals[vid] = normal

        faces.append(tri1)
        faces.append(tri2)

    return surface(vertices, normals, faces)


        

    

def _circle_loop(center_xy, radius, minang):
    arc_geom = [arc(point(center_xy[0], center_xy[1]), radius)]
    loop = geomlist2poly(arc_geom, minang=minang, minlen=0.0)
    if not loop:
        raise ValueError('failed to generate circle loop')
    return loop


def tube(outer_diameter, wall_thickness, length,
         center=None, *, base_point=None, minang=5.0, include_caps=True):
    """Create a cylindrical tube solid.

    ``base_point`` (or legacy ``center`` argument) identifies the base of the
    cylindrical wall, i.e. the plane where ``z == base_point[2]``.
    """

    if base_point is not None and center is not None:
        raise ValueError('specify only base_point (preferred) or center, not both')
    if base_point is None:
        base_point = center if center is not None else point(0, 0, 0)
    if len(base_point) < 3:
        raise ValueError('base_point must contain x, y, z components')
    base_point = point(base_point)

    if wall_thickness <= epsilon:
        raise ValueError('wall thickness must be positive')

    outer_radius = outer_diameter / 2.0
    inner_radius = outer_radius - wall_thickness
    if inner_radius <= epsilon:
        raise ValueError('wall thickness too large for tube')

    base_z = base_point[2]
    center_xy = (base_point[0], base_point[1])

    base_loop_xy = _circle_loop(center_xy, outer_radius, minang)
    inner_loop_xy = list(reversed(_circle_loop(center_xy, inner_radius, minang)))

    base_surface, _ = poly2surfaceXY(base_loop_xy, holepolys=[inner_loop_xy])
    base_surface = reversesurface(base_surface)
    base_surface = translatesurface(base_surface, point(0, 0, base_z))

    top_surface, _ = poly2surfaceXY(base_loop_xy, holepolys=[inner_loop_xy])
    top_surface = translatesurface(top_surface, point(0, 0, base_z + length))

    outer_base = [point(p[0], p[1], base_z, 1.0) for p in base_loop_xy[:-1]]
    inner_base = [point(p[0], p[1], base_z, 1.0) for p in inner_loop_xy[:-1]]
    outer_top = [point(p[0], p[1], base_z + length, 1.0) for p in base_loop_xy[:-1]]
    inner_top = [point(p[0], p[1], base_z + length, 1.0) for p in inner_loop_xy[:-1]]

    outer_side = _loft_surface(outer_base, outer_top)
    inner_side = _loft_surface(inner_top, inner_base, invert=True)

    call = f"yapcad.geom3d_util.tube({outer_diameter}, {wall_thickness}, {length}, base_point={base_point})"
    surfaces = [base_surface, outer_side, top_surface, inner_side]
    if not include_caps:
        surfaces = [outer_side, inner_side]
    return solid(surfaces,
                 [],
                 ['procedure', call])


def conic_tube(bottom_outer_diameter, top_outer_diameter, wall_thickness,
               length, center=None, *, base_point=None, minang=5.0, include_caps=True):
    """Create a conic tube with varying outer diameter.

    ``base_point`` (or ``center`` legacy argument) marks the axial base of the
    frustum (the larger-diameter end when stacked)."""

    if base_point is not None and center is not None:
        raise ValueError('specify only base_point (preferred) or center, not both')
    if base_point is None:
        base_point = center if center is not None else point(0, 0, 0)
    base_point = point(base_point)

    if wall_thickness <= epsilon:
        raise ValueError('wall thickness must be positive')

    r0_outer = bottom_outer_diameter / 2.0
    r1_outer = top_outer_diameter / 2.0
    r0_inner = r0_outer - wall_thickness
    r1_inner = r1_outer - wall_thickness
    if r0_inner <= epsilon or r1_inner <= epsilon:
        raise ValueError('wall thickness too large for conic tube')

    base_z = base_point[2]
    center_xy = (base_point[0], base_point[1])

    base_outer_loop = _circle_loop(center_xy, r0_outer, minang)
    base_inner_loop = list(reversed(_circle_loop(center_xy, r0_inner, minang)))

    top_outer_loop = _circle_loop(center_xy, r1_outer, minang)
    top_inner_loop = list(reversed(_circle_loop(center_xy, r1_inner, minang)))

    base_surface, _ = poly2surfaceXY(base_outer_loop, holepolys=[base_inner_loop])
    base_surface = reversesurface(base_surface)
    base_surface = translatesurface(base_surface, point(0, 0, base_z))

    top_surface, _ = poly2surfaceXY(top_outer_loop, holepolys=[top_inner_loop])
    top_surface = translatesurface(top_surface, point(0, 0, base_z + length))

    outer_base = [point(p[0], p[1], base_z, 1.0) for p in base_outer_loop[:-1]]
    outer_top = [point(p[0], p[1], base_z + length, 1.0) for p in top_outer_loop[:-1]]
    inner_base = [point(p[0], p[1], base_z, 1.0) for p in base_inner_loop[:-1]]
    inner_top = [point(p[0], p[1], base_z + length, 1.0) for p in top_inner_loop[:-1]]

    outer_side = _loft_surface(outer_base, outer_top)
    inner_side = _loft_surface(inner_top, inner_base, invert=True)

    call = ("yapcad.geom3d_util.conic_tube("
            f"{bottom_outer_diameter}, {top_outer_diameter}, {wall_thickness}, {length}, base_point={base_point})")
    surfaces = [base_surface, outer_side, top_surface, inner_side]
    if not include_caps:
        surfaces = [outer_side, inner_side]
    return solid(surfaces,
                 [],
                 ['procedure', call])


def spherical_shell(outer_diameter, wall_thickness,
                    solid_angle=4 * math.pi, center=point(0, 0, 0),
                    *, minang=5.0, steps=24):
    """Create a spherical shell or cap defined by a solid angle."""

    if wall_thickness <= epsilon:
        raise ValueError('wall thickness must be positive')

    outer_radius = outer_diameter / 2.0
    inner_radius = outer_radius - wall_thickness
    if inner_radius <= epsilon:
        raise ValueError('wall thickness too large for spherical shell')

    solid_angle = max(min(solid_angle, 4 * math.pi), epsilon)
    cos_theta = 1.0 - solid_angle / (2.0 * math.pi)
    cos_theta = max(-1.0, min(1.0, cos_theta))
    theta = math.acos(cos_theta)

    arc_samples = max(12, int(round(360.0 / minang)))
    lat_steps = max(4, int(math.ceil(theta / math.radians(minang))))

    cz = center[2]

    def _sphere_contour(radius):
        def contour(z):
            dz = z - cz
            dz = max(min(dz, radius), -radius)
            return math.sqrt(max(radius * radius - dz * dz, 0.0))
        return contour

    z_start_outer = cz + outer_radius * math.cos(theta)
    z_start_inner = cz + inner_radius * math.cos(theta)
    z_end_outer = cz + outer_radius
    z_end_inner = cz + inner_radius

    outer_surface = makeRevolutionSurface(_sphere_contour(outer_radius),
                                          z_start_outer, z_end_outer,
                                          max(steps, lat_steps),
                                          arcSamples=arc_samples)
    outer_surface = translatesurface(outer_surface, point(center[0], center[1], 0))

    inner_surface = makeRevolutionSurface(_sphere_contour(inner_radius),
                                          z_start_inner, z_end_inner,
                                          max(steps, lat_steps),
                                          arcSamples=arc_samples)
    inner_surface = translatesurface(inner_surface, point(center[0], center[1], 0))
    inner_surface = reversesurface(inner_surface)

    surfaces = [outer_surface, inner_surface]

    if theta < math.pi - epsilon:
        r_outer_ring = outer_radius * math.sin(theta)
        r_inner_ring = inner_radius * math.sin(theta)

        base_outer = [point(center[0] + p[0], center[1] + p[1], cz + outer_radius * math.cos(theta), 1.0)
                      for p in _circle_loop((0, 0), r_outer_ring, minang)[:-1]]
        base_inner = [point(center[0] + p[0], center[1] + p[1], cz + inner_radius * math.cos(theta), 1.0)
#                      for p in list(reversed(_circle_loop((0, 0), r_inner_ring, minang)))[:-1]]
                      for p in _circle_loop((0, 0), r_inner_ring, minang)[:-1]]

        conic_surface = _loft_surface(base_outer, base_inner, invert=False)
        surfaces.append(conic_surface)

    call = ("yapcad.geom3d_util.spherical_shell("
            f"{outer_diameter}, {wall_thickness}, {solid_angle}, center={center})")
    return solid(surfaces,
                 [],
                 ['procedure', call])


def stack_solids(solids, *, axis='z', start=0.0, gap=0.0, align='center'):
    """Return translated copies of ``solids`` stacked along an axis."""

    if not solids:
        return []

    axis = axis.lower()
    if axis not in ('x', 'y', 'z'):
        raise ValueError('axis must be one of x, y, or z')

    axis_idx = {'x': 0, 'y': 1, 'z': 2}[axis]
    other_idx = [i for i in range(3) if i != axis_idx]

    placed = []
    cursor = start
    reference = None
    pending_gap = 0.0

    for entry in solids:
        if isinstance(entry, str):
            directive = entry.strip().lower()
            if directive.startswith('space:'):
                try:
                    value = float(directive.split(':', 1)[1])
                except ValueError as exc:
                    raise ValueError(f'bad spacing directive {entry}') from exc
                pending_gap += value
                continue
            raise ValueError(f'unsupported directive {entry!r} in stack_solids')

        solid_obj = entry
        bbox = solidbbox(solid_obj)
        length = bbox[1][axis_idx] - bbox[0][axis_idx]
        if length < epsilon:
            raise ValueError('solid has zero length along stacking axis')

        cursor += pending_gap
        pending_gap = 0.0

        translation = [0.0, 0.0, 0.0]
        translation[axis_idx] = cursor - bbox[0][axis_idx]

        if reference is None:
            if align == 'center':
                reference = [
                    (bbox[0][idx] + bbox[1][idx]) / 2.0 for idx in other_idx
                ]
            elif align == 'min':
                reference = [bbox[0][idx] for idx in other_idx]
            elif align == 'max':
                reference = [bbox[1][idx] for idx in other_idx]
            else:
                raise ValueError('align must be center, min, or max')

        if align == 'center':
            for ref_val, idx in zip(reference, other_idx):
                translation[idx] = ref_val - (bbox[0][idx] + bbox[1][idx]) / 2.0
        elif align == 'min':
            for ref_val, idx in zip(reference, other_idx):
                translation[idx] = ref_val - bbox[0][idx]
        elif align == 'max':
            for ref_val, idx in zip(reference, other_idx):
                translation[idx] = ref_val - bbox[1][idx]

        placed.append(translatesolid(solid_obj, vect(translation[0], translation[1], translation[2], 0)))
        cursor += length + gap

    return placed
