## geom3d, enhanced two- and three-dimensional geometry support for yapCAD
## started on Thu Oct  1 13:52:45 PDT 2020
## Richard W. DeVaul

from yapcad.geom import *
from yapcad.geom_util import *
from functools import reduce

"""
==========================================================
geom3d -- functional 3D geometry representation for yapCAD
==========================================================

The geometric representations of point, line, arc, poly, and
geomlist provided by ``yapcad.geom`` are suitable for representing
zero- or one-dimensional figures embedded in a three-dimensional
space.  And while there is no requirement that the representations
provided by ``yapcad.geom`` (such as an arc, line, or polyline) lie in
the XY plane, many of the functions in yapcad.geom are explicitly
intended to perform computational geometry operations on XY-planar
entities. 

Further, while a closed figure described by ``yapcad.geom`` may
implicitly bound a two-dimensional face, there is no direct support
for working with two-dimensional surfaces provided by that module.
There is no direct support of any kind for working with
three-dimemnsional volumes in ``yapcad.geom``.

In this module we specify representations for two-dimensional surfaces
and bounded three-dimensional volumes, and provide tools for working
with them in implicit, parametric form as well as in explicit,
triangulated form.

The goal of the ``geom3d.yapcad module`` is to allow for the
construction of two-dimensinal surfaces and three-dimensional
geometry for the purposes of modeling, computational geometry, and
rendering.  Specifically, we wish to support the following:

(1) Support the implicit representation of three-dimensional
geometry, and the performance of constructive solid geometry
operations on this implicit geometry (union, intersection,
difference) to produce more complex implicit three dimensional
forms. 

(2) Support the implicit representation of two dimensional
surfaces, such as a planar surface specified by three points, or
the surface of a three-dimensional object like a sphere, and allow
for computational geometry operations on these surfaces, such as
intersection operations, to produce explicit one-dimensional
objects, such as lines, arcs, etc.

(3) support the conversion of an implicit two-dimensional surface
to an explicit, teselated triangluar geometry that may be easily
rendered using a conventional 3D graphics rendering pipline, such
as OpenGL

(4) Support for the conversion of implicit three-dimenaional
constructive solid geometry into an explicit, contiguous closed
surface representation using the marching cubes algortihm, or any
other user-specified conversion algoritm, for the purposes of
interactive 3D rendering and conversion to 3D CAM formats, such as
STL.

Structures for 2D/3D geometry:

``surface = ['surface',vertices,normals,faces,boundary,holes]``, where:

           ``vertices`` is a list of ``yapcad.geom`` points,

           ``normals`` is a list of ``yapcad.geom`` points of the same
           length as ``vertices``,

           ``faces`` is the list of faces, which is to say lists
           of three indices that refer to the vertices of the triangle
           that represents each face,

           ``boundary`` is a list of indices for the vertices that form
           the outer perimeter of the surface.

           ``holes`` is a (potentially zero-length) list of lists of
           holes, each of which is a non-zero list of three or more
           indices of vertices that form the perimeter of any holes in
           the surface.

A surface has an inside and an outside.  To deterine which side of a
surface a point lies on, find the closest face and determine if the
point lies on the positive or negative side of the face.

``solid = ['solid', bbox, surfaces ]``, where:

    ``bbox`` is the bounding box for the solid

    ``surfaces`` is a list of surfaces with contiguous boundaries that
    completely encloses an interior space.

"""

def signedPlaneDistance(p,p0,n):
    """Given a point on the plane ``p0``, and the unit-length plane
    normal ``n``, determine the signed distance to the plane.
    """
    return dot(sub(p,p0),n)

def tri2p0n(face):
    """Given ``face``, a non-degenrate poly of length 3, return the center
    point and normal of the face, otherwise known as the Hessian
    Normal Form: https://mathworld.wolfram.com/HessianNormalForm.html

    """
    p1=face[0]
    p2=face[1]
    p3=face[2]
    p0 = scale3(add(face[0],add(face[1],face[2])),1.0/3.0)
    v1=sub(p2,p1)
    v2=sub(p3,p2)
    c= cross(v1,v2)
    m = mag(c)
    if m < epsilon:
        raise ValueError('degenerate face in tri2p0n')
    n= scale3(c,1.0/m)
    return [p0,n]

def signedFaceDistance(p,face):
    """given a test point ``p`` and a three-point face ``face``, determine
    the minimum signed distance of the test point from the face.  Any
    point lying in the positive normal direction or on the surface of
    the plane will result in a zero or positive distace.  Any point
    lying in the negative normal direction from the surface of the
    plane will result in a negative distance.

    """
    p0,n = tuple(tri2p0n(face))
    d = sub(p,face[0])
    m = -dot(d,n) # negative of distance of p from plane
    a = add(p,scale3(n,m)) # projection of point p into plane

    # create a coordinate system based on the first face edge and the
    # plane normal
    v1 = sub(face[1],face[0])
    v2 = sub(face[2],face[0])
    vx = scale3(v1,1.0/mag(v1))
    vy = scale3(cross(vx,n),-1)
    vz = n
    p1=[0,0,0,1]
    p2=[mag(v1),0,0,1]
    p3=[dot(v2,vx),dot(v2,vy),0,1]
    aa= sub(a,face[0])
    aa=[dot(aa,vx),dot(aa,vy),0,1]

    #barycentric coordinates for derermining if projected point falls
    #inside face
    lam1,lam2,lam3 = tuple(barycentricXY(aa,p1,p2,p3))
    inside = ( (lam1 >= 0.0 and lam1 <= 1.0) and
               (lam2 >= 0.0 and lam2 <= 1.0) and
               (lam3 >= 0.0 and lam3 <= 1.0) )

    ind = 0 # in-plane distance is zero if projected point inside
            # triangle
    if not inside:
        d1 = linePointXY([p1,p2],aa,distance=True)
        d2 = linePointXY([p2,p3],aa,distance=True)
        d3 = linePointXY([p3,p1],aa,distance=True)
        ind = min(d1,d2,d3) #in-plane distance is smallest distance from each edge
            
    if close(m,0): # point lies in plane
        return ind
    else:
        dist = sqrt(m*m+ind*ind) # total distance is hypotenuse of
                                 # right-triangle of in-plane and
                                 # out-plane distance
        return copysign(dist,-1*m)
    

def triTriIntersect(t1,t2,inside=True,inPlane=False,plane1=None):
    """Function to compute the intersection of two triangles.  Returns
    ``False`` if no intersection, a line (a list of two points) if the
    planes do not overlap and there is a linear intersection, and a
    polygon (list of three or more points) if the triangles are
    co-planar and overlap.
    
    If ``inside == True`` (default) return line-segment or poly
    intersection that falls inside both bounded triangles, otherwise
    return linear intersection of two planes, or False if planes are
    degenerate.

    If ``inPlane==True``, return the intersection as a poly in the
    planar coordinate system implied by ``t1``, or in the planar
    coordinate system specified by ``plane1``

    If ``plane1`` is not ``False``, it should be planar coordinate
    system in the form of ``[p0,x,y,z]``, where ``p0`` specifies the
    origin in world coordinates, and ``x``, ``y``, and ``z`` are
    orthonormal basis vectors in which ``z`` is the plane normal.

    NOTE: when ``plane1`` is True and ``inPlane`` is False, it is
    assumed that ``plane1`` is a previously-computed planar basis that
    is coplanar with ``t1`` and will be used to speed computation. 

    """
    return False
    
def issurface(s,fast=True):

    """
    Check to see if ``s`` is a valid surface.
    """
    def filterInds(inds,verts):
        l = len(verts)
        if l < 3:
            return False
        return (len(list(filter(lambda x: not (isinstance(x,int) or
                                               x < 0 or x >= l),
                                inds))) == 0)
    
    if not isinstance(s,list) or s[0] != 'surface' or len(s) != 6:
        return False
    if fast:
        return True
    else:
        verts=s[1]
        norms=s[2]
        faces=s[3]
        boundary= s[4]
        holes= s[5]
        if (not ispoly(verts) or
            not ispoly(norms) or
            len(verts) != len(norms)):
            return False
        l = len(verts)
        if (len(list(filter(lambda x: not len(x) == 3, faces))) > 0):
            return False
        if not filterInds(reduce( (lambda x,y: x + y),faces),verts):
            return False
        if not filterInds(boundary,verts):
            return False
        if len(holes)>0:
            for h in holes:
                if not filterInds(h,verts):
                    return False
        return True

def normfunc(tri):
    """
    utility funtion to compute normals for a flat facet triangle
    """
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
    boundary = s[4]
    holes = s[5]

    n1,n2,n3 = nfunc(tri)
    i1,vrts,nrms = addVert(tri[0],n1,vrts,nrms)
    i2,vrts,nrms = addVert(tri[1],n2,vrts,nrms)
    i3,vrts,nrms = addVert(tri[2],n3,vrts,nrms)
    faces.append([i1,i2,i3])

    return ['surface',vrts,nrms,faces,boundary,holes]

def poly2surface(ply,holepolys=[],minlen=0.5,minarea=0.0001,
                 checkclosed=False,box=None):
    """Given ``ply``, an XY-coplanar polygon, return the triangulated
    surface representation of that polygon. If ``holepolys`` is not
    the empty list, treat each polygon in that list as a hole in
    ``ply``.  If ``checkclosed`` is true, make sure ``ply`` and all
    members of ``holepolys`` are a vaid, closed, XY-coplanar polygons.
    if ``box`` exists, use it as the bounding box.

    """
    
    if checkclosed:
        polys = holepolys + [ply]
        if not isgeomlistXYPlanar(polys):
            raise ValueError('non-XY-coplanar arguments')
        for p in polys:
            if not ispolygonXY(p):
                raise ValueError(f'{p} is not a closed polygon')

    if len(holepolys) != 0:
        raise NotImplementedError('not supporting surfaces with holes, yet')
    
    if not box:
        box = bbox(ply)
    
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
        #print(f"len bndry: {len(bndry)}")
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
        #print(f"len bndry: {len(bndry)}")
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

    def makeboundary(poly,vertices,normals):
        bndry = []
        i = len(vertices)
        for p in poly:
            vertices.append(point(p))
            normals.append([0,0,1,1])
            bndry.append(i)
            i+=1
        return bndry,vertices,normals
    
    vrts=[]
    nrms=[]
    faces=[]
    boundary=[]
    holes=[]

    bndry=deepcopy(ply[0:-1]) # last point is redundant
    # make the perimeter
    
    boundary,vrts,nrms = makeboundary(bndry,vrts,nrms)
    for h in holepolys:
        hole,vrts,nrms = makeboundary(h[0:-1],vrts,nrms)
        holes.append(hole)

    surf=['surface',vrts,nrms,faces,boundary,holes]
    cnt = 0
    addcnt = 0
    while (len(bndry) > 2):
        cnt+= 1
        if cnt > 1000:
            assert False
            return surf, bndry
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
                                          nfunc=lambda x: ([0,0,1,1],
                                                           [0,0,1,1],
                                                           [0,0,1,1]))
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
                    #print(f"i : {i} len(bndry): {len(bndry)}")
                    skip=False
                else:
                    bndry,diddiv = subdivide(i,bndry)
                    div = div or diddiv
                    if diddiv:
                        #print(f"did div at {i}")
                        skip = True

                
            if not div:
                print("incomplete triangulation")
                print(f"cnt: {cnt}, addcnt: {addcnt}")
                print(f"bndry len: {len(bndry)}, surface triangles: {len(surf[3])}")
                print(f"bndry: {vstr(bndry)}")
                return surf,bndry
                #raise ValueError('unable to finish triangulating surface')
        

    return surf,bndry

def surfArea(surf):
    """
    given a surface, return the surface area
    """
    area = 0.0
    vertices = surf[1]
    faces= surf[3]
    for f in faces:
        area += triarea(vertices[f[0]],
                        vertices[f[1]],
                        vertices[f[2]])

    return area

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
    

