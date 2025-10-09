## geom3d, enhanced two- and three-dimensional geometry support for yapCAD
## started on Thu Oct  1 13:52:45 PDT 2020
## Richard W. DeVaul

from yapcad.geom import *
from yapcad.geom_util import *
from yapcad.xform import *
from functools import reduce
from yapcad.octtree import NTree

from yapcad.triangulator import triangulate_polygon


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

Topology Analysis
-----------------

Two key functions are provided for analyzing solid topology:

``issolidclosed(solid)`` -- Verifies that a solid is topologically closed
by ensuring every edge is shared by exactly two faces across all surfaces.
This is essential for determining if a solid properly encloses a volume
without gaps or holes. Returns True if closed, False otherwise.

``volumeof(solid)`` -- Calculates the volume enclosed by a closed solid
using the divergence theorem. Requires that the solid be topologically
closed (verified by calling ``issolidclosed()``). Returns the volume
as a non-negative floating point number.

Example usage::

    from yapcad.geom3d_util import prism
    from yapcad.geom3d import issolidclosed, volumeof

    cube = prism(2, 2, 2)
    if issolidclosed(cube):
        vol = volumeof(cube)  # Returns 8.0


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
    if len(args) >= 3 and len(args) <= 6:
        vrts = args[0]
        nrms = args[1]
        facs = args[2]
        bndr = []
        hle = []
        metadata = None

        if not (isinstance(vrts,list) and isinstance(nrms,list)
                and isinstance(facs,list) and len(vrts) == len(nrms)):
            raise ValueError('bad arguments to surface')

        extras = list(args[3:])
        for item in extras:
            if isinstance(item, dict):
                if metadata is not None:
                    raise ValueError('multiple metadata dictionaries passed to surface')
                metadata = item
            elif isinstance(item, list):
                if not bndr:
                    bndr = item
                elif not hle:
                    hle = item
                else:
                    raise ValueError('too many list arguments passed to surface')
            else:
                raise ValueError('bad arguments to surface')

        surf = ['surface',vrts,nrms,facs,bndr,hle]
        if metadata is not None:
            if not isinstance(metadata, dict):
                raise ValueError('surface metadata must be a dict')
            surf.append(metadata)
        if issurface(surf,fast=False):
            return surf
    raise ValueError('bad arguments to surface')

def surfacebbox(s):
    """return bounding box for surface"""
    if not issurface(s):
        raise ValueError('bad surface passed to surfacebbox')
    return polybbox(s[1])
    
def _ensure_surface_metadata_dict(s):
    if not isinstance(s, list) or len(s) < 4 or s[0] != 'surface':
        raise ValueError('bad surface passed to _ensure_surface_metadata_dict')
    if len(s) == 4:
        s.extend([[], [], {}])
        return s[6]
    if len(s) == 5:
        s.extend([[], {}])
        return s[6]
    if len(s) == 6:
        s.append({})
        return s[6]
    metadata = s[6]
    if metadata is None:
        metadata = {}
    elif not isinstance(metadata, dict):
        metadata = {'legacy_metadata': metadata}
    s[6] = metadata
    return metadata


def _triangle_bbox(tri, tol):
    mins = [min(pt[i] for pt in tri) - tol for i in range(3)]
    maxs = [max(pt[i] for pt in tri) + tol for i in range(3)]
    return [point(mins[0], mins[1], mins[2]),
            point(maxs[0], maxs[1], maxs[2])]


def _triangle_plane_intersection_points(tri, plane, tol):
    n, d = plane
    points = []
    for idx in range(3):
        p0 = tri[idx]
        p1 = tri[(idx + 1) % 3]
        dir_vec = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]]
        denom = n[0] * dir_vec[0] + n[1] * dir_vec[1] + n[2] * dir_vec[2]
        if abs(denom) < tol:
            continue
        numer = d - (n[0] * p0[0] + n[1] * p0[1] + n[2] * p0[2])
        t = numer / denom
        if t < -tol or t > 1.0 + tol:
            continue
        t_clamped = max(0.0, min(1.0, t))
        pt = _lerp_point(p0, p1, t_clamped)
        if abs(n[0] * pt[0] + n[1] * pt[1] + n[2] * pt[2] - d) <= tol * 10.0:
            points.append(pt)
    return _unique_points(points, tol) if points else []


def _candidate_planes_for_triangle(tri, target, tri_plane, tol):
    box = solidbbox(target)
    if box:
        center = point((box[0][0] + box[1][0]) / 2.0,
                        (box[0][1] + box[1][1]) / 2.0,
                        (box[0][2] + box[1][2]) / 2.0)
        extent = max(box[1][0] - box[0][0],
                     box[1][1] - box[0][1],
                     box[1][2] - box[0][2])
    else:
        center = point(0, 0, 0)
        extent = 0.0

    epsilon_dot = max(extent * 1e-6, 1e-6)

    tri_bbox = _triangle_bbox(tri, tol * 2.0)
    seen = set()
    planes = []
    snap_candidates = []
    for surf in target[1]:
        tree = surface_octree(surf)
        if tree is None:
            candidates = list(_iter_triangles_from_surface(surf))
        else:
            elems = tree.getElements(tri_bbox)
            candidates = []
            for elem in elems:
                if isinstance(elem, list):
                    candidates.append(elem)
                elif isinstance(elem, tuple):
                    candidates.append(elem[0])
                else:
                    candidates.append(elem)
            if not candidates:
                candidates = list(_iter_triangles_from_surface(surf))
        for cand in candidates:
            plane = _triangle_plane(cand)
            n, d = plane
            centroid = point(
                (cand[0][0] + cand[1][0] + cand[2][0]) / 3.0,
                (cand[0][1] + cand[1][1] + cand[2][1]) / 3.0,
                (cand[0][2] + cand[1][2] + cand[2][2]) / 3.0,
            )
            vec = point(centroid[0] - center[0],
                        centroid[1] - center[1],
                        centroid[2] - center[2])
            dot_sign = n[0] * vec[0] + n[1] * vec[1] + n[2] * vec[2]
            if dot_sign > epsilon_dot:
                sense = -1
            elif dot_sign < -epsilon_dot:
                sense = 1
            else:
                center_eval = _plane_eval(plane, center)
                sense = -1 if center_eval <= 0 else 1
            plane_with_sense = (n, d, sense)
            key = _plane_key(plane_with_sense, tol)
            if key not in seen:
                seen.add(key)
                planes.append(plane_with_sense)
            snap_candidates.extend(_triangle_plane_intersection_points(cand, tri_plane, tol))
    unique_snap = _unique_points(snap_candidates, max(tol * 10.0, 1e-6)) if snap_candidates else []
    return planes, unique_snap


_DEFAULT_RAY_TOL = 1e-7


def invalidate_surface_octree(s):
    meta = _ensure_surface_metadata_dict(s)
    meta.pop('_octree', None)
    meta['_octree_dirty'] = True


def _bbox_overlap(box_a, box_b, tol):
    return not (
        box_a[1][0] < box_b[0][0] - tol or
        box_a[0][0] > box_b[1][0] + tol or
        box_a[1][1] < box_b[0][1] - tol or
        box_a[0][1] > box_b[1][1] + tol or
        box_a[1][2] < box_b[0][2] - tol or
        box_a[0][2] > box_b[1][2] + tol
    )


def _segment_bbox(p0, p1, pad):
    return [
        point(min(p0[0], p1[0]) - pad, min(p0[1], p1[1]) - pad, min(p0[2], p1[2]) - pad),
        point(max(p0[0], p1[0]) + pad, max(p0[1], p1[1]) + pad, max(p0[2], p1[2]) + pad),
    ]


def _plane_eval(plane, p):
    n, d = plane
    return dot(n, p) - d


def _segment_plane_intersection(p1, p2, plane, tol):
    n, d = plane
    direction = sub(p2, p1)
    denom = dot(n, direction)
    if abs(denom) < tol:
        return point((p1[0] + p2[0]) * 0.5, (p1[1] + p2[1]) * 0.5, (p1[2] + p2[2]) * 0.5)
    t = (d - dot(n, p1)) / denom
    t = max(0.0, min(1.0, t))
    return point(p1[0] + direction[0] * t,
                 p1[1] + direction[1] * t,
                 p1[2] + direction[2] * t)


def _dedupe_polygon(poly, tol):
    if not poly:
        return []
    deduped = [poly[0]]
    for pt in poly[1:]:
        if dist(pt, deduped[-1]) > tol:
            deduped.append(pt)
    if len(deduped) > 2 and dist(deduped[0], deduped[-1]) <= tol:
        deduped[-1] = deduped[0]
        deduped.pop()
    return deduped


def _clip_polygon_against_plane(poly, plane, tol, keep_inside=True):
    if not poly:
        return []
    n, d, sense = plane
    clipped = []
    prev = poly[-1]
    prev_eval = _plane_eval((n, d), prev)
    if sense <= 0:
        prev_inside = prev_eval <= tol if keep_inside else prev_eval >= -tol
    else:
        prev_inside = prev_eval >= -tol if keep_inside else prev_eval <= tol
    for curr in poly:
        curr_eval = _plane_eval((n, d), curr)
        if sense <= 0:
            curr_inside = curr_eval <= tol if keep_inside else curr_eval >= -tol
        else:
            curr_inside = curr_eval >= -tol if keep_inside else curr_eval <= tol
        if curr_inside:
            if not prev_inside:
                clipped.append(_segment_plane_intersection(prev, curr, (n, d), tol))
            clipped.append(point(curr))
        elif prev_inside:
            clipped.append(_segment_plane_intersection(prev, curr, (n, d), tol))
        prev = curr
        prev_eval = curr_eval
        prev_inside = curr_inside
    clipped = _dedupe_polygon(clipped, tol)
    if not keep_inside and len(clipped) >= 3:
        clipped = list(reversed(clipped))
    return clipped


def _split_polygon_by_plane(poly, plane, tol):
    if not poly:
        return [], []

    n, d, sense = plane
    evals = [_plane_eval((n, d), p) for p in poly]
    max_eval = max(evals)
    min_eval = min(evals)

    if sense <= 0:
        if max_eval <= tol:
            return [point(p) for p in poly], []
        if min_eval >= -tol:
            return [], [[point(p) for p in poly]]
    else:
        if min_eval >= -tol:
            return [point(p) for p in poly], []
        if max_eval <= tol:
            return [], [[point(p) for p in poly]]

    inside = _clip_polygon_against_plane(poly, plane, tol, keep_inside=True)
    outside = _clip_polygon_against_plane(poly, plane, tol, keep_inside=False)
    outside_polys = [outside] if outside else []
    return inside, outside_polys


def _split_polygon_by_planes(poly, planes, tol):
    inside_polys = [poly]
    outside_polys = []
    for plane in planes:
        next_inside = []
        for current in inside_polys:
            inside, outside = _split_polygon_by_plane(current, plane, tol)
            outside_polys.extend(outside)
            if inside:
                next_inside.append(inside)
        inside_polys = next_inside
        if not inside_polys:
            break
    return inside_polys, outside_polys


def _triangulate_polygon(poly, reference_normal=None):
    if len(poly) < 3:
        return []
    if len(poly) == 3:
        tri = [point(poly[0]), point(poly[1]), point(poly[2])]
        if reference_normal is not None:
            v01 = sub(tri[1], tri[0])
            v02 = sub(tri[2], tri[0])
            if dot(cross(v01, v02), reference_normal) < 0:
                tri[1], tri[2] = tri[2], tri[1]
        return [tri]

    anchor = point(poly[0])
    triangles = []
    for i in range(1, len(poly) - 1):
        tri = [anchor, point(poly[i]), point(poly[i + 1])]
        if reference_normal is not None:
            v01 = sub(tri[1], tri[0])
            v02 = sub(tri[2], tri[0])
            if dot(cross(v01, v02), reference_normal) < 0:
                tri[1], tri[2] = tri[2], tri[1]
        triangles.append(tri)
    return triangles


def _triangle_plane(tri):
    p0, n = tri2p0n(tri)
    d = dot(n, p0)
    return n, d


def _plane_key(plane, tol):
    n, d, sense = plane
    scale = 1.0 / tol
    return (round(n[0] * scale), round(n[1] * scale), round(n[2] * scale), round(d * scale), sense)


def _sub3(a, b):
    return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]


def _dot3(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross3(a, b):
    return [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]


def _mag3(v):
    return math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])


def _normalize3(v, tol):
    m = _mag3(v)
    if m < tol:
        return [0.0, 0.0, 0.0]
    return [v[0] / m, v[1] / m, v[2] / m]


def _lerp_point(p0, p1, t):
    return point(
        p0[0] + (p1[0] - p0[0]) * t,
        p0[1] + (p1[1] - p0[1]) * t,
        p0[2] + (p1[2] - p0[2]) * t,
    )



def _point_in_triangle(pt, tri, tol):
    a = [tri[0][0], tri[0][1], tri[0][2]]
    b = [tri[1][0], tri[1][1], tri[1][2]]
    c = [tri[2][0], tri[2][1], tri[2][2]]
    p = [pt[0], pt[1], pt[2]]

    v0 = _sub3(c, a)
    v1 = _sub3(b, a)
    v2 = _sub3(p, a)

    dot00 = _dot3(v0, v0)
    dot01 = _dot3(v0, v1)
    dot02 = _dot3(v0, v2)
    dot11 = _dot3(v1, v1)
    dot12 = _dot3(v1, v2)

    denom = dot00 * dot11 - dot01 * dot01
    if abs(denom) < tol:
        return False
    inv = 1.0 / denom
    u = (dot11 * dot02 - dot01 * dot12) * inv
    v = (dot00 * dot12 - dot01 * dot02) * inv
    return u >= -tol and v >= -tol and (u + v) <= 1.0 + tol


def _segment_triangle_intersection_param(p0, p1, tri, tol):
    n, d = _triangle_plane(tri)
    dir_vec = _sub3([p1[0], p1[1], p1[2]], [p0[0], p0[1], p0[2]])
    denom = _dot3(n, dir_vec)
    if abs(denom) < tol:
        return None
    numer = d - (n[0] * p0[0] + n[1] * p0[1] + n[2] * p0[2])
    t = numer / denom
    if t < -tol or t > 1.0 + tol:
        return None
    t_clamped = max(0.0, min(1.0, t))
    intersection = _lerp_point(p0, p1, t_clamped)
    if not _point_in_triangle(intersection, tri, tol):
        return None
    return (t_clamped, intersection)


def _collect_segment_intersections(p0, p1, solid, tol):
    results = []
    query_box = _segment_bbox(p0, p1, tol * 5)
    for surf in solid[1]:
        tree = surface_octree(surf)
        if tree is None:
            candidates = list(_iter_triangles_from_surface(surf))
        else:
            elems = tree.getElements(query_box)
            candidates = []
            for elem in elems:
                if isinstance(elem, list):
                    candidates.append(elem)
                elif isinstance(elem, tuple):
                    candidates.append(elem[0])
                else:
                    candidates.append(elem)
        for tri in candidates:
            res = _segment_triangle_intersection_param(p0, p1, tri, tol)
            if res is not None:
                results.append(res)
    return results


def _unique_points(points, tol):
    unique = []
    for pt in points:
        candidate = point(pt)
        if not any(dist(candidate, existing) <= tol for existing in unique):
            unique.append(candidate)
    return unique


def _points_to_polygon(points, reference_normal, tol):
    unique = _unique_points(points, tol)
    if len(unique) < 3:
        return []

    centroid = point(
        sum(p[0] for p in unique) / len(unique),
        sum(p[1] for p in unique) / len(unique),
        sum(p[2] for p in unique) / len(unique),
    )

    basis_u = None
    for pt in unique:
        vec = _sub3([pt[0], pt[1], pt[2]], [centroid[0], centroid[1], centroid[2]])
        if _mag3(vec) > tol:
            basis_u = _normalize3(vec, tol)
            break
    if basis_u is None:
        return []

    ref = _normalize3([reference_normal[0], reference_normal[1], reference_normal[2]], tol)
    basis_v = _cross3(ref, basis_u)
    if _mag3(basis_v) < tol:
        basis_v = _cross3(basis_u, ref)
    basis_v = _normalize3(basis_v, tol)

    def _angle(pt):
        vec = _sub3([pt[0], pt[1], pt[2]], [centroid[0], centroid[1], centroid[2]])
        x = _dot3(vec, basis_u)
        y = _dot3(vec, basis_v)
        return math.atan2(y, x)

    ordered = sorted(unique, key=_angle)
    compact = [ordered[0]]
    for pt in ordered[1:]:
        if dist(compact[-1], pt) > tol:
            compact.append(pt)
    if dist(compact[0], compact[-1]) <= tol:
        compact[-1] = compact[0]
        compact = compact[:-1]
    if len(compact) < 3:
        return []

    v0 = _sub3([compact[1][0], compact[1][1], compact[1][2]], [compact[0][0], compact[0][1], compact[0][2]])
    v1 = _sub3([compact[2][0], compact[2][1], compact[2][2]], [compact[1][0], compact[1][1], compact[1][2]])
    normal = _cross3(v0, v1)
    if _mag3(normal) < tol:
        return []
    if _dot3(normal, [reference_normal[0], reference_normal[1], reference_normal[2]]) < 0:
        compact.reverse()

    return compact


def _orient_triangle(tri, reference_normal, tol):
    pts = [point(tri[0]), point(tri[1]), point(tri[2])]
    v01 = sub(pts[1], pts[0])
    v02 = sub(pts[2], pts[0])
    normal = cross(v01, v02)
    area_mag = mag(normal)
    if area_mag < tol:
        return None
    ref = reference_normal
    if ref is not None:
        ref_vec = [ref[0], ref[1], ref[2]]
        if dot(normal, ref_vec) < 0:
            pts[1], pts[2] = pts[2], pts[1]
    return pts


def _clip_triangle_against_solid(tri, target, tol, reference_normal):
    polygon = [point(tri[0]), point(tri[1]), point(tri[2])]
    centroid = point(
        (tri[0][0] + tri[1][0] + tri[2][0]) / 3.0,
        (tri[0][1] + tri[1][1] + tri[2][1]) / 3.0,
        (tri[0][2] + tri[1][2] + tri[2][2]) / 3.0,
    )

    tri_plane = _triangle_plane(tri)
    planes, snap_candidates = _candidate_planes_for_triangle(tri, target, tri_plane, tol)

    inside_polys = []
    if planes:
        inside_raw, _ = _split_polygon_by_planes(polygon, planes, tol)
        for poly in inside_raw:
            cleaned = _dedupe_polygon([point(p) for p in poly], tol)
            if len(cleaned) >= 3:
                inside_polys.append(cleaned)

    if not inside_polys:
        if solid_contains_point(target, centroid, tol=tol):
            inside_polys = [polygon]
        else:
            oriented = _orient_triangle(tri, reference_normal, tol)
            if oriented:
                return [], [oriented], False
            return [], [], False

    _, _, to_local, to_world = tri2p0n(tri, basis=True)

    if snap_candidates:
        snap_tol = max(5e-2, tol * 1e6)

        def _snap_point(pt):
            best = snap_tol
            best_candidate = None
            for cand in snap_candidates:
                d = dist(pt, cand)
                if d < best:
                    best = d
                    best_candidate = cand
            return point(best_candidate) if best_candidate is not None else point(pt)

        snapped = []
        for poly in inside_polys:
            snapped_poly = [_snap_point(pt) for pt in poly]
            snapped_poly = _dedupe_polygon(snapped_poly, tol)
            if len(snapped_poly) >= 3:
                snapped.append(snapped_poly)
        if snapped:
            inside_polys = snapped

    def _to_local(pt):
        loc = to_local.mul(pt)
        return (loc[0], loc[1])

    def _from_local(xy):
        local_pt = point(xy[0], xy[1], 0.0)
        world_pt = to_world.mul(local_pt)
        return point(world_pt[0], world_pt[1], world_pt[2])

    outer_loop = [_to_local(pt) for pt in polygon]
    inside_loops = [[_to_local(pt) for pt in poly] for poly in inside_polys]

    inside_triangles = []
    for loop in inside_loops:
        tri2d_list = triangulate_polygon(loop)
        for tri2d in tri2d_list:
            tri3d = [_from_local(pt) for pt in tri2d]
            oriented = _orient_triangle(tri3d, reference_normal, tol)
            if oriented:
                inside_triangles.append(oriented)

    holes = inside_loops if inside_loops else None
    outside_triangles = []
    tri2d_list = triangulate_polygon(outer_loop, holes=holes)
    for tri2d in tri2d_list:
        tri3d = [_from_local(pt) for pt in tri2d]
        oriented = _orient_triangle(tri3d, reference_normal, tol)
        if oriented:
            outside_triangles.append(oriented)

    split = bool(inside_triangles) and bool(outside_triangles)
    return inside_triangles, outside_triangles, split




def stitch_open_edges(triangles, tol):
    """Experimental: attempt to close open edge loops by triangulating them.

    Parameters
    ----------
    triangles : Iterable
        A sequence of triangle coordinate lists ``[[x, y, z, w], ...]``
        describing a mesh with potential boundary edges.
    tol : float
        Tolerance used when deduplicating vertices and detecting shared
        edges. Typically reuse ``_DEFAULT_RAY_TOL``.

    Returns
    -------
    list
        A new list of triangles with any stitched faces appended.
    """
    if not triangles:
        return triangles

    oriented_edges = []
    canonical_counts = {}
    base_triangles = []
    for tri in triangles:
        try:
            n = _triangle_normal(tri)
        except ValueError:
            continue
        pts = [point(tri[0]), point(tri[1]), point(tri[2])]
        base_triangles.append(pts)
        edges = ((pts[0], pts[1]), (pts[1], pts[2]), (pts[2], pts[0]))
        for start, end in edges:
            canonical = _canonical_edge_key(start, end)
            canonical_counts[canonical] = canonical_counts.get(canonical, 0) + 1
            oriented_edges.append((start, end, n))

    boundary = []
    for start, end, normal in oriented_edges:
        if canonical_counts[_canonical_edge_key(start, end)] == 1:
            boundary.append((start, end, normal))

    if not boundary:
        return base_triangles

    def edge_id(start, end):
        return (_point_to_key(start), _point_to_key(end))

    adjacency = {}
    for start, end, normal in boundary:
        adjacency.setdefault(_point_to_key(start), []).append((start, end, normal))

    used = set()
    loops = []
    for start, end, normal in boundary:
        eid = edge_id(start, end)
        if eid in used:
            continue
        loop = []
        normals = []
        current_start = start
        current_end = end
        current_normal = normal
        start_key = _point_to_key(current_start)
        while True:
            loop.append(point(current_start))
            normals.append(current_normal)
            used.add(edge_id(current_start, current_end))
            next_key = _point_to_key(current_end)
            if next_key == start_key:
                break
            candidates = adjacency.get(next_key)
            if not candidates:
                loop = []
                break
            next_edge = None
            for candidate in candidates:
                cid = edge_id(candidate[0], candidate[1])
                if cid not in used:
                    next_edge = candidate
                    break
            if next_edge is None:
                loop = []
                break
            current_start, current_end, current_normal = next_edge
        if loop and len(loop) >= 3:
            loops.append((loop, normals))

    if not loops:
        return base_triangles

    stitched = list(base_triangles)
    for loop_points, normals in loops:
        polygon = _unique_points(loop_points, tol)
        if len(polygon) < 3:
            continue
        avg = [0.0, 0.0, 0.0]
        for n in normals:
            avg[0] += n[0]
            avg[1] += n[1]
            avg[2] += n[2]
        avg = _normalize3(avg, tol)
        if _mag3(avg) < tol:
            avg = normals[0]
        poly = _points_to_polygon(polygon, avg, tol)
        if not poly:
            continue
        for tri in _triangulate_polygon(poly, avg):
            try:
                normal = _triangle_normal(tri)
            except ValueError:
                continue
            if _dot3(normal, [avg[0], avg[1], avg[2]]) < 0:
                tri = [tri[0], tri[2], tri[1]]
                try:
                    normal = _triangle_normal(tri)
                except ValueError:
                    continue
            stitched.append(tri)
    return stitched

def stitch_solid(sld, tol=_DEFAULT_RAY_TOL):
    """Experimental helper that runs :func:`stitch_open_edges` on a solid."""

    if not issolid(sld, fast=False):
        raise ValueError('invalid solid passed to stitch_solid')

    triangles = list(_iter_triangles_from_solid(sld))
    stitched = stitch_open_edges(triangles, tol)
    if not stitched:
        return sld

    surface_result = _surface_from_triangles(stitched)
    if surface_result:
        return solid([surface_result], [], ['utility', 'stitch_open_edges'])
    return sld



def _boolean_fragments(source, target, tol):
    outside_tris = []
    inside_tris = []
    inside_overlap = []

    for tri in _iter_triangles_from_solid(source):
        reference_normal = _triangle_normal(tri)
        inside_parts, outside_parts, split = _clip_triangle_against_solid(
            tri, target, tol, reference_normal
        )
        if inside_parts:
            inside_tris.extend(inside_parts)
        if outside_parts:
            outside_tris.extend(outside_parts)
        if split and inside_parts:
            inside_overlap.extend(inside_parts)
    return outside_tris, inside_tris, inside_overlap



def _ray_triangle_intersection(origin, direction, triangle, tol=_DEFAULT_RAY_TOL):
    v0, v1, v2 = triangle
    e1 = sub(v1, v0)
    e2 = sub(v2, v0)
    h = cross(direction, e2)
    a = dot(e1, h)
    if abs(a) < tol:
        return None
    f = 1.0 / a
    s = sub(origin, v0)
    u = f * dot(s, h)
    if u < -tol or u > 1.0 + tol:
        return None
    q = cross(s, e1)
    v = f * dot(direction, q)
    if v < -tol or u + v > 1.0 + tol:
        return None
    t = f * dot(e2, q)
    if t < -tol:
        return None
    w = 1.0 - u - v
    if w < -tol:
        return None
    hit = point(origin[0] + direction[0] * t,
                origin[1] + direction[1] * t,
                origin[2] + direction[2] * t)
    return t, hit, (u, v, w)


def surface_octree(s, rebuild=False):
    meta = _ensure_surface_metadata_dict(s)
    tree = meta.get('_octree')
    if rebuild or meta.get('_octree_dirty') or tree is None:
        verts = s[1]
        faces = s[3]
        if not faces:
            meta['_octree'] = None
            meta['_octree_dirty'] = False
            return None
        extent = 0.0
        for axis in range(3):
            vals = [v[axis] for v in verts]
            extent = max(extent, max(vals) - min(vals))
        mindim = max(extent / 32.0, epsilon)
        tree = NTree(n=8, mindim=mindim)
        for face in faces:
            if len(face) != 3:
                continue
            tri = [verts[face[0]], verts[face[1]], verts[face[2]]]
            tree.addElement(tri)
        tree.updateTree()
        meta['_octree'] = tree
        meta['_octree_dirty'] = False
    return meta['_octree']


def _surface_from_triangles(triangles):
    if not triangles:
        return None
    verts = []
    normals = []
    faces = []
    for tri in triangles:
        v01 = sub(tri[1], tri[0])
        v02 = sub(tri[2], tri[0])
        if mag(cross(v01, v02)) < epsilon:
            continue
        n = _triangle_normal(tri)
        indices = []
        for pt in tri:
            verts.append(point(pt))
            normals.append([n[0], n[1], n[2], 0.0])
            indices.append(len(verts) - 1)
        faces.append(indices)
    surface_obj = ['surface', verts, normals, faces, [], []]
    _ensure_surface_metadata_dict(surface_obj)
    invalidate_surface_octree(surface_obj)
    return surface_obj


def _iter_triangles_from_surface(surf):
    verts = surf[1]
    faces = surf[3]
    for face in faces:
        if len(face) != 3:
            continue
        yield [verts[face[0]], verts[face[1]], verts[face[2]]]


def _iter_triangles_from_solid(sld):
    for surf in sld[1]:
        yield from _iter_triangles_from_surface(surf)


def _group_hits(hits, tol):
    if not hits:
        return []
    hits.sort(key=lambda x: x[0])
    groups = [[hits[0]]]
    for hit in hits[1:]:
        if abs(hit[0] - groups[-1][-1][0]) <= tol:
            groups[-1].append(hit)
        else:
            groups.append([hit])
    return groups


def _triangle_normal(tri):
    _, n = tri2p0n(tri)
    return n
def solid_contains_point(sld, p, tol=_DEFAULT_RAY_TOL):
    if not issolid(sld, fast=False):
        raise ValueError('invalid solid passed to solid_contains_point')
    if not ispoint(p):
        raise ValueError('invalid point passed to solid_contains_point')

    surfaces = sld[1]
    if not surfaces:
        return False

    box = solidbbox(sld)
    if box:
        expanded = [
            point(box[0][0] - tol, box[0][1] - tol, box[0][2] - tol),
            point(box[1][0] + tol, box[1][1] + tol, box[1][2] + tol),
        ]
        if not isinsidebbox(expanded, p):
            return False
        extent = max(box[1][0] - box[0][0],
                     box[1][1] - box[0][1],
                     box[1][2] - box[0][2])
        ray_length = max(extent * 1.5, 1.0)
    else:
        ray_length = 1.0

    directions = [
        vect(1, 0, 0, 0),
        vect(-1, 0, 0, 0),
        vect(0, 1, 0, 0),
        vect(0, -1, 0, 0),
        vect(0, 0, 1, 0),
        vect(0, 0, -1, 0),
    ]

    seen_inside = False
    for direction in directions:
        far_point = point(p[0] + direction[0] * ray_length,
                          p[1] + direction[1] * ray_length,
                          p[2] + direction[2] * ray_length)
        query_box = _segment_bbox(p, far_point, tol * 5)
        hits = []
        for surf in surfaces:
            tree = surface_octree(surf)
            if tree is None:
                candidates = list(_iter_triangles_from_surface(surf))
            else:
                elems = tree.getElements(query_box)
                candidates = []
                for elem in elems:
                    if isinstance(elem, list):
                        candidates.append(elem)
                    elif isinstance(elem, tuple):
                        candidates.append(elem[0])
                    else:
                        candidates.append(elem)
                if not candidates:
                    candidates = list(_iter_triangles_from_surface(surf))
            for tri in candidates:
                result = _ray_triangle_intersection(p, direction, tri, tol=tol)
                if result is None:
                    continue
                t_hit, hit_point, _ = result
                if t_hit < -tol or t_hit > ray_length + tol:
                    continue
                if t_hit <= tol:
                    return True
                normal = _triangle_normal(tri)
                sign = -1 if dot(normal, direction) > 0 else 1
                hits.append((t_hit, hit_point, sign))
        groups = _group_hits(hits, tol)
        parity = 0
        for group in groups:
            sign_sum = sum(hit[2] for hit in group)
            if sign_sum == 0:
                continue
            parity ^= 1
        if parity == 0:
            return False
        seen_inside = True
    return seen_inside


def solids_intersect(a, b, tol=_DEFAULT_RAY_TOL):
    if not issolid(a, fast=False) or not issolid(b, fast=False):
        raise ValueError('invalid solid passed to solids_intersect')
    if not a[1] or not b[1]:
        return False

    box_a = solidbbox(a)
    box_b = solidbbox(b)
    if box_a and box_b and not _bbox_overlap(box_a, box_b, tol):
        return False

    center_a = point((box_a[0][0] + box_a[1][0]) / 2.0,
                     (box_a[0][1] + box_a[1][1]) / 2.0,
                     (box_a[0][2] + box_a[1][2]) / 2.0)
    center_b = point((box_b[0][0] + box_b[1][0]) / 2.0,
                     (box_b[0][1] + box_b[1][1]) / 2.0,
                     (box_b[0][2] + box_b[1][2]) / 2.0)
    if solid_contains_point(a, center_b, tol=tol) or solid_contains_point(b, center_a, tol=tol):
        return True

    for tri_a in _iter_triangles_from_solid(a):
        for tri_b in _iter_triangles_from_solid(b):
            if triTriIntersect(tri_a, tri_b):
                return True
    return False


def solid_boolean(a, b, operation, tol=_DEFAULT_RAY_TOL, *, stitch=False):
    if operation not in {'union', 'intersection', 'difference'}:
        raise ValueError(f'unsupported solid boolean operation {operation!r}')

    outside_a, inside_a, overlap_a = _boolean_fragments(a, b, tol)
    outside_b, inside_b, overlap_b = _boolean_fragments(b, a, tol)

    if operation == 'union':
        result_tris = outside_a + outside_b
        if not result_tris:
            result_tris = inside_a if inside_a else inside_b
    elif operation == 'intersection':
        result_tris = inside_a + inside_b
        if not result_tris:
            result_tris = overlap_a + overlap_b
    else:  # difference
        reversed_inside_b = [[tri[0], tri[2], tri[1]] for tri in inside_b]
        if not reversed_inside_b and overlap_b:
            reversed_inside_b = [[tri[0], tri[2], tri[1]] for tri in overlap_b]
        result_tris = outside_a + reversed_inside_b

    if stitch:
        result_tris = stitch_open_edges(result_tris, tol)

    surface_result = _surface_from_triangles(result_tris)
    if surface_result:
        return solid([surface_result], [], ['boolean', operation])
    return solid([], [], ['boolean', operation])



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
    
    if not isinstance(s,list) or s[0] != 'surface' or len(s) not in (6,7):
        return False
    if fast:
        return True
    else:
        verts=s[1]
        norms=s[2]
        faces=s[3]
        boundary= s[4]
        holes= s[5]
        metadata = s[6] if len(s) == 7 else None
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
        if metadata is not None and not isinstance(metadata, dict):
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
    if len(args) >= 1 and len(args) <= 4:
        if not isinstance(args[0],list):
            raise ValueError('bad arguments to solid')
        for srf in args[0]:
            if not issurface(srf):
                raise ValueError('bad arguments to solid')

        surfaces = args[0]
        material = []
        construction = []
        metadata = None

        for item in args[1:]:
            if isinstance(item, dict):
                if metadata is not None:
                    raise ValueError('multiple metadata dictionaries passed to solid')
                metadata = item
            elif isinstance(item, list):
                if material == []:
                    material = item
                elif construction == []:
                    construction = item
                else:
                    raise ValueError('too many list arguments passed to solid')
            else:
                raise ValueError('bad arguments to solid')

        sld = ['solid', surfaces, material, construction]
        if metadata is not None:
            sld.append(metadata)
        return sld

    raise ValueError('bad arguments to solid')
            
    
    
def issolid(s,fast=True):

    """
    Check to see if ``s`` is a solid.  NOTE: this function only determines
    if th data structure is correct, it does not verify that the collection
    of surfaces completely bounds a volume of space without holes
    """

    if not isinstance(s,list) or s[0] != 'solid' or len(s) not in (4,5):
        return False
    if fast:
        return True
    else:
                         
        for surface in s[1]:
            if not issurface(surface,fast=fast):
                return False
        if not (isinstance(s[2],list) and isinstance(s[3],list)):
            return False
        if len(s) == 5 and not isinstance(s[4], dict):
            return False
        return True

def solidbbox(sld):
    if not issolid(sld):
        raise ValueError('bad argument to solidbbox')

    box = []
    for surf in sld[1]:
        sb = surfacebbox(surf)
        if not box:
            box = sb
        else:
            box = [point(min(box[0][0], sb[0][0]),
                         min(box[0][1], sb[0][1]),
                         min(box[0][2], sb[0][2])),
                   point(max(box[1][0], sb[1][0]),
                         max(box[1][1], sb[1][1]),
                         max(box[1][2], sb[1][2]))]

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

def _point_to_key(p):
    """
    Convert a point to a hashable key for edge/vertex identification.
    Uses rounded coordinates to handle floating point precision issues.
    """
    # Round to a reasonable precision to handle floating point comparison
    return (round(p[0] / epsilon) * epsilon,
            round(p[1] / epsilon) * epsilon,
            round(p[2] / epsilon) * epsilon)

def _canonical_edge_key(p1, p2):
    """
    Create a canonical edge key from two points.
    Returns tuple of keys in sorted order so (p1,p2) and (p2,p1) map to same edge.
    """
    k1 = _point_to_key(p1)
    k2 = _point_to_key(p2)
    return (min(k1, k2), max(k1, k2))

def issolidclosed(x):
    """
    Check if solid x is topologically closed.

    A solid is closed if and only if every edge is shared by exactly two
    faces across all surfaces. This ensures no holes or gaps exist in the
    solid's boundary.

    The function analyzes face adjacency by:
    1. Building a global edge map using vertex positions (not indices)
    2. Counting how many faces share each edge
    3. Verifying that every edge is shared by exactly 2 faces

    Args:
        x: A solid data structure

    Returns:
        True if the solid is topologically closed, False otherwise

    Raises:
        ValueError: if x is not a valid solid

    Example:
        >>> from yapcad.geom3d_util import prism, sphere
        >>> cube = prism(2, 2, 2)
        >>> issolidclosed(cube)
        True
    """
    # First verify this is a valid solid
    if not issolid(x, fast=False):
        raise ValueError('invalid solid passed to issolidclosed')

    surfaces = x[1]

    # Empty solid is trivially closed
    if not surfaces:
        return True

    # Build a global edge map across all surfaces
    # Key: canonical edge tuple (point_key1, point_key2)
    # Value: count of faces that share this edge
    global_edge_count = {}

    for surf_idx, surf in enumerate(surfaces):
        faces = surf[3]
        vertices = surf[1]

        # Process each face in this surface
        for face_idx, face in enumerate(faces):
            if len(face) != 3:
                raise ValueError(f'non-triangular face in surface {surf_idx}, face {face_idx}')

            # Get the three vertex positions
            p0 = vertices[face[0]]
            p1 = vertices[face[1]]
            p2 = vertices[face[2]]

            # Extract the three edges of this triangular face
            edges = [
                _canonical_edge_key(p0, p1),
                _canonical_edge_key(p1, p2),
                _canonical_edge_key(p2, p0)
            ]

            # Count this face's contribution to each edge
            for edge in edges:
                if edge not in global_edge_count:
                    global_edge_count[edge] = 0
                global_edge_count[edge] += 1

    # Check that every edge is shared by exactly 2 faces
    for edge, count in global_edge_count.items():
        if count != 2:
            return False

    return True

def volumeof(x):
    """
    Calculate the volume enclosed by a solid.

    Uses the divergence theorem to compute volume from the surface triangulation.
    For each triangular face with vertices (p0, p1, p2), the signed volume
    contribution is: V_i = (1/6) * dot(p0, cross(p1-p0, p2-p0))

    The total volume is the sum of absolute values of all face contributions.

    Args:
        x: A solid data structure (must be topologically closed)

    Returns:
        float: The volume of the solid (always non-negative)

    Raises:
        ValueError: if x is not a valid solid or is not closed

    Example:
        >>> from yapcad.geom3d_util import prism
        >>> cube = prism(2, 2, 2)
        >>> abs(volumeof(cube) - 8.0) < 0.001
        True
    """
    # Verify this is a valid, closed solid
    if not issolid(x, fast=False):
        raise ValueError('invalid solid passed to volumeof')

    if not issolidclosed(x):
        raise ValueError('solid must be topologically closed to compute volume')

    surfaces = x[1]

    # Handle empty solid
    if not surfaces:
        return 0.0

    total_volume = 0.0

    # Accumulate signed volume contributions from all faces
    for surf in surfaces:
        vertices = surf[1]
        faces = surf[3]

        for face in faces:
            if len(face) != 3:
                raise ValueError('non-triangular face encountered')

            # Get the three vertices of this face and normalise to w=1 if needed
            p0 = point(vertices[face[0]])
            p1 = point(vertices[face[1]])
            p2 = point(vertices[face[2]])

            # Compute vectors from p0 to other vertices
            v1 = sub(p1, p0)
            v2 = sub(p2, p0)

            # Signed volume contribution: (1/6) * dot(p0, cross(v1, v2))
            # This is the volume of the tetrahedron formed by the origin
            # and the three face vertices
            cross_product = cross(v1, v2)
            signed_volume = dot(p0, cross_product) / 6.0

            total_volume += signed_volume

    # Return absolute value (orientation might cause negative result)
    return abs(total_volume)

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
    for hole in holepolys:
        holes2.append(list(map(lambda x: basis[2].mul(x), hole)))

    surf,bnd = poly2surfaceXY(ply2,holes2,minlen,minarea,checkclosed)

    verts2 = list(map(lambda x: basis[3].mul(x),surf[1]))
    norm2 = list(map(lambda x: basis[3].mul([x[0],x[1],x[2],0]),surf[2]))
    bnd2 = list(map(lambda x: basis[3].mul(x),bnd))

    surf[1]=verts2
    surf[2]=norm2

    return surf,bnd2

def poly2surfaceXY(ply,holepolys=[],minlen=0.5,minarea=0.0001,
                   checkclosed=False,box=None):
    """Given ``ply``, return a triangulated XY surface (holes supported)."""

    if checkclosed:
        polys = holepolys + [ply]
        if not isgeomlistXYPlanar(polys):
            raise ValueError('non-XY-coplanar arguments')
        for p in polys:
            if not ispolygonXY(p):
                raise ValueError(f'{p} is not a closed polygon')

    if not box:
        box = bbox(ply)

    def _normalize_loop(poly):
        pts = [point(p) for p in poly]
        if pts and dist(pts[0], pts[-1]) <= epsilon:
            pts = pts[:-1]
        return pts

    outer_loop = _normalize_loop(ply)
    if len(outer_loop) < 3:
        raise ValueError('degenerate polygon passed to poly2surfaceXY')

    hole_loops = [_normalize_loop(loop) for loop in holepolys]

    triangles = triangulate_polygon([(p[0], p[1]) for p in outer_loop],
                                    [[(q[0], q[1]) for q in loop]
                                     for loop in hole_loops])

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

    boundary,vrts,nrms = makeboundary(outer_loop,vrts,nrms)
    for loop in hole_loops:
        hole,vrts,nrms = makeboundary(loop,vrts,nrms)
        holes.append(hole)

    surf=['surface',vrts,nrms,faces,boundary,holes]

    def _signed_triangle(tri):
        (x1,y1),(x2,y2),(x3,y3) = tri
        return ((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1)) / 2.0

    outer_area = sum(outer_loop[i][0]*outer_loop[(i+1)%len(outer_loop)][1]
                     - outer_loop[(i+1)%len(outer_loop)][0]*outer_loop[i][1]
                     for i in range(len(outer_loop)))
    orientation = 1 if outer_area >= 0 else -1

    for tri in triangles:
        tri_points = [point(x, y, 0, 1) for (x, y) in tri]
        area = _signed_triangle(tri)
        if area * orientation < 0:
            tri_points[1], tri_points[2] = tri_points[2], tri_points[1]
            area = -area
        if abs(area) > minarea:
            surf = addTri2Surface(tri_points,surf,
                                  nfunc=lambda x: ([0,0,1,0],
                                                   [0,0,1,0],
                                                   [0,0,1,0]))

    return surf,[]

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
  
                     
