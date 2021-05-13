## geom3d, enhanced two- and three-dimensional geometry support for yapCAD
## started on Thu Oct  1 13:52:45 PDT 2020
## Richard W. DeVaul

from yapcad.geom import *
from yapcad.geom_util import *
from yapcad.xform import *
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

Structures for 2D/3D geometry
=============================

surfaces
--------

``surface = ['surface',vertices,normals,faces,boundary,holes]``, where:

           ``vertices`` is a list of ``yapcad.geom`` points,

           ``normals`` is a list of ``yapcad.geom`` direction vectors
           of the same length as ``vertices``,

           ``faces`` is the list of faces, which is to say lists
           of three indices that refer to the vertices of the triangle
           that represents each face,

           ``boundary`` is a list of indices for the vertices that form
           the outer perimeter of the surface, or [] if the surface
           has no boundary, such as that of a torus or sphere

           ``holes`` is a (potentially zero-length) list of lists of
           holes, each of which is a non-zero list of three or more
           indices of vertices that form the perimeter of any holes in
           the surface.

A surface has an inside and an outside.  To deterine which side of a
surface a point lies on, find the closest face and determine if the
point lies on the positive or negative side of the face.

Solids
------

To represent a completely bounded space (a space for which any point
can be unambiguously determined to be inside or outside) there is the
``solid`` representation.  A solid is composed of zero or more
surfaces, and may be completely empty, as empty solids are legal
products of constructive solid geometry operations like intersection
and difference.

The gurantee, which is not enforced by the representation, is that for
any point inside the bounding box of the solid, and for any point
chosen outside the bounding box of the solid, a line drawn between the
two will have either even (or zero), or odd point intersections (as
opposed to tangent intersections), regardless of the choice of the
outside-the-bounding-box point.

Solids have optional associated metadata about the material properties
of the solid and about how it was constructed.  

For example, material properties might include OpenGL-type rendering
data, mechanical properties, or a reference to a material dictionary
that includes both.

Construction meta-data might include the model-file and material-file
name from which the geometry was loaded, the polygon from which the
solid was extruded (and associated extrusion parameters), or the
function call with parameters for algorithmically-generated geometry. 

``solid = ['solid', surfaces, material, construction ]``, where:

           ``surfaces`` is a list of surfaces with contiguous boundaries 
           that completely encloses an interior space,

           ``material`` is a list of domain-specific representation of
           the material properties of the solid, which may be empty.
           This information may be used for rendering or simulating
           the properties of the solid.

           ``construction`` is a list that contains information about
           how the solid was constructed, and may be empty


Assembly
--------

Assemblies are lists of elements, which are solids or assemblies, in
which each list has an associated geometric transformation.

``assembly = ['assembly', transform, elementlist]``, where:

            ``transform = [xformF, xformR]``, a pair of forward and inverse
              transformation matricies such that ``xformF`` * ``xformR`` = 
              the identity matrix.  

            ``elementlist = [element0, elememnt1, ... ]``, in which each
              list element is either a valid ``solid`` or ``assembly``.

"""

def signedPlaneDistance(p,p0,n):
    """Given a point on the plane ``p0``, and the unit-length plane
    normal ``n``, determine the signed distance to the plane.
    """
    return dot(sub(p,p0),n)

def tri2p0n(face,basis=False):
    """Given ``face``, a non-degenrate poly of length 3, return the center
    point and normal of the face, otherwise known as the Hessian
    Normal Form: https://mathworld.wolfram.com/HessianNormalForm.html

    In addition, if ``basis==True``, calculate an orthnormal basis
    vectors implied by the triangle, with the x' vector aligned with
    the p1-p2 edge, the z' vector as the nornal vector, and the y'
    vector as -1 * x' x z', and return the transformation matrix that
    will transform a point in world coordinates to a point in the
    orthonormal coordinate system with the origin at the center
    point. Return that transformation matrix and its inverse.

    """
    # import pdb ; pdb.set_trace()
    p1=face[0]
    p2=face[1]
    p3=face[2]
    p0 = scale3(add(p1,add(p2,p3)),1.0/3.0)
    v1=sub(p2,p1)
    v2=sub(p3,p2)
    c= cross(v1,v2)
    m = mag(c)
    if m < epsilon:
        raise ValueError('degenerate face in tri2p0n')
    n= scale3(c,1.0/m)
    n[3] = 0.0 #direction vectors lie in the w=0 hyperplane
    if not basis:
        return [p0,n]
    else:
        # compute orthonormal basis vectors
        x=scale3(v1,1.0/mag(v1))
        z=n
        y=scale3(cross(x,z),-1)
        # direction vectors lie in the w=0 hyperplane
        x[3] = y[3] = z[3] = 0.0
        # build the rotation matrix
        T = [x,
             y,
             z,
             [0,0,0,1]]
        rm = Matrix(T)
        # build translation matrix
        tm = Translation(scale3(p0,-1))
        # make composed matrix
        forward = rm.mul(tm)
        # make inverse matrix
        rm.trans=True
        inverse =  Translation(p0).mul(rm)
        return [p0,n,forward,inverse]
    

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

def linePlaneIntersect(lne,plane="xy",inside=True):
    """Function to calculate the intersection of a line and a plane.

    ``line`` is specified in the usual way as two points.  ``plane``
    is either specified symbolicly as one of ``["xy","yz","xz"]``, a
    list of three points, or a planar coordinate system in the form of
    ``[p0,n,forward,reverse]``, where ``p0`` specifies the origin in
    world coordinates, ``n`` is the normal (equivalent to the ``z``
    vector), and ``forward`` and ``reverse`` are transformation
    matricies that map from world into local and local into world
    coordinates respectively.

    Returns ``False`` if the line and plane do not intersect, or if
    ``inside==True`` and the point of intersection is outside the line
    interval.  Returns the point of intersection otherwise.

    NOTE: if plane is specified as three points, then setting
    ``inside=True`` will also force a check to see if the intersection
    point falls within the specified triangle.

    """

    def lineCardinalPlaneIntersect(lne,idx,inside=True):
        if close(lne[0][idx]-lne[1][idx],0.0): #degenerate
            return False
        # (1-u)*l[0][idx] + u*l[1][idx] = 0.0
        # u*(l[1][idx]-l[0][idx]) = l[0][idx]
        # u = l[0][idx]/(l[1][idx]-l[0][idx])
        u = lne[0][idx]/(lne[0][idx]-lne[1][idx])
        if inside and (u < 0.0 or u > 1.0):
            return False
        else:
            return sampleline(lne,u)

    # is the plane specified symbolicaly?
    trangle = False
    idx = -1
    if plane=="xy":
        idx=2
    elif plane=="yz":
        idx=0
    elif plane=="xz":
        idx=1

    if idx > -1:
        return lineCardinalPlaneIntersect(lne,idx,inside)
    else:
        if istriangle(plane):
            triangle = True
            tri = plane
            plane = tri2p0n(plane,basis=True)
            tri2 = [ plane[2].mul(tri[0]),
                     plane[2].mul(tri[1]),
                     plane[2].mul(tri[2]) ]
        else:
            raise ValueError('non-plane passed to linePlaneIntersect')
        # otherwise assume that plane is a valid planar basis

        # transform into basis with plane at z=0
        l2 = [plane[2].mul(lne[0]),plane[2].mul(lne[1])]
        p = lineCardinalPlaneIntersect(l2,2,inside)
        if not p:
            return False
        else:
            if triangle and not isInsideTriangleXY(p,tri2):
                return False
            return plane[3].mul(p)
        

def triTriIntersect(t1,t2,inside=True,inPlane=False,basis=None):

    """Function to compute the intersection of two triangles.  Returns
    ``False`` if no intersection, a line (a list of two points) if the
    planes do not overlap and there is a linear intersection, and a
    polygon (list of three or more points) if the triangles are
    co-planar and overlap.
    
    If ``inside == True`` (default) return line-segment or poly
    intersection that falls inside both bounded triangles, otherwise
    return a line segment that lies on the infinite linear
    intersection of two planes, or False if planes are degenerate.

    If ``inPlane==True``, return the intersection as a poly in the
    planar coordinate system implied by ``t1``, or in the planar
    coordinate system specified by ``basis``

    If ``basis`` is not ``False``, it should be planar coordinate
    system in the form of ``[p0,n,forward,reverse]``, where ``p0``
    specifies the origin in world coordinates, ``n`` is the normal
    (equivalent to the ``z`` vector), and ``forward`` and ``reverse``
    are transformation matricies that map from world into local and
    local into world coordinates respectively.

    NOTE: when ``basis`` is True and ``inPlane`` is False, it is
    assumed that ``basis`` is a planar basis computed by tri2p0n
    coplanar with ``t1``.

    """
    if not basis:
        #create basis from t1
        basis = tri2p0n(t1,basis=True)
    p01,n1,tfor,tinv = tuple(basis)

    p02,n2 = tuple(tri2p0n(t2,basis=False))

    #transform both triangles into new coordinate system
    t1p = list(map(lambda x: tfor.mul(x),t1))
    t2p = list(map(lambda x: tfor.mul(x),t2))

    # check for coplanar case
    if (abs(t2p[0][2]) <= epsilon and abs(t2p[1][2]) <= epsilon
        and abs(t2p[2][2]) <= epsilon):
        if not inside:
            if inPlane:
                return t2p
            else:
                return t2
        else: # return poly that is in-plane intersection
            intr = combineglist(t1p,t2p,'intersection')
            if len(intr) < 1: #no intersection
                return False
            else:
                if inPlane:
                    return intr
                else:
                    return transform(intr,tinv)
    # not coplanar, check to see if planes are parallel
    if vclose(n1,n2):
        return False #yep, degenerate

    if inside:
        # check to see if t2p lies entirely above or below the z=0 plane
        if ((t2p[0][2] > epsilon and t2p[1][2] > epsilon and t2p[2][2] > epsilon)
            or
            (t2p[0][2] < -epsilon and t2p[1][2] < -epsilon and
             t2p[2][2] < -epsilon)):
            return False
        # linear intersection.  Figure out which two of three lines
        # cross the z=0 plane

    # this should work whether or not the intersection is
    # inside t2.
    ip1 = linePlaneIntersect([t2p[0],t2p[1]],"xy",False)
    ip2 = linePlaneIntersect([t2p[1],t2p[2]],"xy",False)
    ip3 = linePlaneIntersect([t2p[2],t2p[0]],"xy",False)

    a=ip1
    b=ip2
    if not a:
        a=ip3
    if not b:
        b=ip3
    if inPlane:
        return [a,b]
    else:
        return [tinv.mul(a),tinv.mul(b)]
    
def surface(*args):
    """given a surface or a list of surface parameters as arguments,
    return a conforming surface representation.  Checks arguments
    for data-type correctness.

    """
    if args==[]:
        # empty surface
        return ['surface',[],[],[],[],[] ]
    if len(args) == 1:
        # one argument, produce a deep copy of surface (it that is
        # what it is)
        if issurface(args[0],fast=False):
            return deepcopy(args[0])
    if len(args) > 2 and len(args) < 6:
        vrts = args[0]
        nrms = args[1]
        facs = args[2]
        bndr = []
        hle = []
        if (isinstance(vrts,list) and isinstance(nrms,list)
            and isinstance(facs,list) and
            len(vrts) == len(nrms)):
            if len(args) > 3 and isinstance(args[3],list):
                bndr = args[3]
            if len(args) == 5 and isinstance(args[4],list):
                hle = args[4]
            surf = ['surface',vrts,nrms,facs,bndr,hle]
            if issurface(surf,fast=False):
                return surf
    raise ValueError('bad arguments to surface')

def surfacebbox(s):
    """return bounding box for surface"""
    if not issurface(s):
        raise ValueError('bad surface passed to surfacebbox')
    return polybbox(s[1])
            
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
            not isdirectlist(norms) or
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

## save pointers to the yapcad.geom transformation functions
geom_rotate = rotate
geom_translate = translate
geom_scale = scale
geom_mirror = mirror

def rotatesurface(s,ang,cent=point(0,0,0),axis=point(0,0,1.0),mat=False):
    """ return a rotated copy of the surface"""
    if close(ang,0.0):
        return deepcopy(s)
    if not mat: # if matrix isn't pre-specified, calculate it
        if vclose(cent,point(0,0,0)):
            mat = xform.Rotation(axis,ang)
        else:
            mat = xform.Translation(cent)
            mat = mat.mul(xform.Rotation(axis,ang))
            mat = mat.mul(xform.Translation(cent,inverse=True))
    s2 = deepcopy(s)
    #import pdb ; pdb.set_trace()
    s2[1] = geom_rotate(s2[1],ang,cent,axis,mat)
    s2[2] = geom_rotate(s2[2],ang,cent,axis,mat)
    return s2

def translatesurface(s,delta):
    """ return a translated copy of the surface"""
    if vclose(delta,point(0,0,0)):
        return deepcopy(s)
    s2 = deepcopy(s)
    for i in range(len(s2[1])):
        s2[1][i] = add(s2[1][i],delta)
    return s2

def mirrorsurface(s,plane):
    """return a mirrored version of a surface.  Currently, the following
    values of "plane" are allowed: 'xz', 'yz', xy'.  Generalized
    arbitrary reflection plane specification will be added in the
    future.

    Note that this is a full surface geometry reflection, and not
    simply a normal reverser.
    """
    s2 = deepcopy(s)
    s2[1] = mirror(s[1],plane)
    s2[2] = mirror(s[2],plane)
    s2[3] = list(map(lambda x: [x[0],x[2],x[1]],s2[3]))
    return s2
    
def reversesurface(s):
    """ return a nomal-reversed copy of the surface """
    s2 = deepcopy(s)
    s2[2] = list(map(lambda x: scale4(x,-1.0),s2[2]))
    s2[3] = list(map(lambda x: [x[0],x[2],x[1]],s2[3]))
    return s2

def solid(*args):
    """given a solid or a list of solid parameters as arguments,
    return a conforming solid representation.  Checks arguments
    for data-type correctness.

    """
    if args==[] or (len(args) == 1 and args[0] == []):
        # empty solid, which is legal because we must support
        # empty results of CSG operations, etc.
        return ['solid',[],[],[] ]

    # check for "copy constructor" case
    if len(args) == 1 and issolid(args[0],fast=False):
        # one argument, it's a solid.  produce a deep copy of solid
        return deepcopy(args[0])
    
    # OK, step through arguments
    if len(args) > 0 and len(args) < 4:
        surfaces = []
        material = []
        construction = []
        
        if isinstance(args[0],list): # is 1st argument a surface list?
            for s in args[0]:
                if not issurface(s):
                    raise ValueError('bad arguments to solid')
        else:
            raise ValueError('bad arguments to solid')

        surfaces = args[0]
        if len(args) > 1:
            if isinstance(args[1],list):
                material = args[1]
            else:
                raise ValueError('bad material arguments to solid')
        if len(args) == 3:
            if isinstance(args[2],list):
                construction = args[2]
            else:
                raise ValueError('bad construction arguments to solid')

        return ['solid', surfaces, material, construction]
        
    raise ValueError('bad arguments to solid')
            
    
    
def issolid(s,fast=True):

    """
    Check to see if ``s`` is a solid.  NOTE: this function only determines
    if th data structure is correct, it does not verify that the collection
    of surfaces completely bounds a volume of space without holes
    """

    if not isinstance(s,list) or s[0] != 'solid' or len(s) != 4:
        return False
    if fast:
        return True
    else:
                         
        for surface in s[1]:
            if not issurface(surface,fast=fast):
                return False
        if not (isinstance(s[2],list) and isinstance(s[3],list)):
            return False
        return True

def solidbbox(sld):
    if not issolid(sld):
        raise ValueError('bad argument to solidbbox')

    box = []
    for s in sld[1]:
        box = surfacebbox(s + box)

    return box

def translatesolid(x,delta):
    if not issolid(x):
        raise ValueError('bad solid passed to translatesolid')
    s2 = deepcopy(x)
    surfs = []
    for s in x[1]:
        surfs.append(translatesurface(s,delta))
    s2[1] = surfs    
    return s2

def rotatesolid(x,ang,cent=point(0,0,0),axis=point(0,0,1.0),mat=False):
    if not issolid(x):
        raise ValueError('bad solid passed to rotatesolid')
    s2 = deepcopy(x)
    surfs=[]
    for s in x[1]:
        surfs.append(rotatesurface(s,ang,cent=cent,axis=axis,mat=mat))
    s2[1] = surfs
    return s2

def mirrorsolid(x,plane,preserveNormal=True):
    if not issolid(x):
        raise ValueError('bad solid passed to mirrorsolid')
    s2 = deepcopy(x)
    surfs=[]
    for s in x[1]:
        surf = mirrorsurface(s,plane)
        if preserveNormal and False:
            surf = reversesurface(surf)
        surfs.append(surf)
    s2[1] = surfs
    return s2

def normfunc(tri):
    """
    utility funtion to compute normals for a flat facet triangle
    """
    v1 = sub(tri[1],tri[0])
    v2 = sub(tri[2],tri[1])
    d = cross(v1,v2)
    n = scale3(d,1.0/mag(d))
    n[3] = 0.0 # direction vectors lie in the w=0 hyperplane
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


def surfacearea(surf):
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
    
def poly2surface(ply,holepolys=[],minlen=0.5,minarea=0.0001,
                 checkclosed=False,basis=None):

    """Given ``ply``, a coplanar polygon, return the triangulated surface
    representation of that polygon and its boundary.  If ``holepolys``
    is not the empty list, treat each polygon in that list as a hole
    in ``ply``.  If ``checkclosed`` is true, make sure ``ply`` and all
    members of ``holepolys`` are a vaid, closed, coplanar polygons.
    if ``box`` exists, use it as the bounding box.

    if ``basis`` exists, use it as the planar coordinate basis to
    transform the poly into the z=0 plane.

    Returns surface and boundary

    """

    if len(ply) < 3:
        raise ValueError(f'poly must be at least length 3, got {len(ply)}')

    if not basis:
        v0 = sub(ply[1],ply[0])
        v1 = None
        for i in range(2,len(ply)):
            v1 = sub(ply[i],ply[1])
            if mag(cross(v0,v1)) > epsilon:
                break
        if not v1:
            raise ValueError(f'degenerate poly passed to poly2surface')
        basis = tri2p0n([ply[0],ply[1],ply[i]],basis=True)

    ply2 = list(map(lambda x: basis[2].mul(x),ply))
    holes2 = []
    for holes in holepolys:
        holes2.append(list(map(lambda x: basis[2].mul(x),hole)))

    surf,bnd = poly2surfaceXY(ply2,holes2,minlen,minarea,checkclosed)

    verts2 = list(map(lambda x: basis[3].mul(x),surf[1]))
    norm2 = list(map(lambda x: basis[3].mul([x[0],x[1],x[2],0]),surf[2]))
    bnd2 = list(map(lambda x: basis[3].mul(x),bnd))

    surf[1]=verts2
    surf[2]=norm2

    return surf,bnd2

def poly2surfaceXY(ply,holepolys=[],minlen=0.5,minarea=0.0001,
                   checkclosed=False,box=None):
    """Given ``ply``, an XY-coplanar polygon, return the triangulated
    surface representation of that polygon. If ``holepolys`` is not
    the empty list, treat each polygon in that list as a hole in
    ``ply``.  If ``checkclosed`` is true, make sure ``ply`` and all
    members of ``holepolys`` are a vaid, closed, XY-coplanar polygons.
    if ``box`` exists, use it as the bounding box.

    Returns surface and boundary
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
            normals.append([0,0,1,0])
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
                                          nfunc=lambda x: ([0,0,1,0],
                                                           [0,0,1,0],
                                                           [0,0,1,0]))
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

### updated, surface- and solid-aware generalized geometry functions

# length -- scalar length doesn't make sense for sufface or solid

geom_center = center
def center(x):
    """Return the point corresponding to the center of surface, solid, or
    figure x.

    """
    if issurface(x):
        box = surfacebbox(x)
        return scale3(add(box[0],box[1]),0.5)
    elif issolid(x):
        box = solidbbox(x)
        return scale3(add(box[0],box[1]),0.5)
    else:
        return geom_center(x)

geom_bbox = bbox
def bbox(x):
    """Given a figure, surface, or solid x, return the three-dimensional
    bounding box of that entity."""
    if issolid(x):
        return solidbbox(x)
    elif issurface(x):
        return surfacebbox(x)
    else:
        return geom_bbox(x)

# sample -- doesn't make sense for surface or solid
# unsample -- doesn't make sense for surface or solid
# segment -- doesn't make sense for surface or solid
# isnsideXY -- doesn't make sense for a suface or solid

def translate(x,delta):
    """ return a translated version of the surface, solid, or figure"""
    if issolid(x):
        return translatesolid(x,delta)
    elif issurface(x):
        return translatesurface(x,delta)
    else:
        return geom_translate(x,delta)
  
def rotate(x,ang,cent=point(0,0),axis=point(0,0,1.0),mat=False):
    """ return a rotated version of the surface, solid, or figure"""
    if issolid(x):
        return rotatesolid(x,ang,cent=cent,axis=axis,mat=mat)
    elif issurface(x):
        return rotatesurface(x,ang,cent=cent,axis=axis,mat=mat)
    else:
        return geom_rotate(x,ang,cent=cent,axis=axis,mat=mat)
  

def mirror(x,plane):
    """
    return a mirrored version of a figure.  Currently, the following
    values of "plane" are allowed: 'xz', 'yz', xy'.  Generalized
    arbitrary reflection plane specification will be added in the
    future.

    NOTE: this operation will reverse the sign of the area of ``x`` if
    x is a closed polyline or geometry list
    """
    if issolid(x):
        return mirrorsolid(x,plane)
    elif issurface(x):
        return mirrorsurface(x,plane)
    else:
        return geom_mirror(x,plane)
  
                     
