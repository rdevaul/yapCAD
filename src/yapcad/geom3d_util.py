## geom3d_util, additional 3D geometry support for yapCAD
## started on Mon Feb 15 20:23:57 PST 2021 Richard W. DeVaul

from yapcad.geom import *
from yapcad.geom_util import *
from yapcad.xform import *
from yapcad.geom3d import *
from yapcad.brep import attach_brep_to_solid, occ_available, BrepSolid

import math

"""
==================================================
Utility functions to support 3D geometry in yapCAD
==================================================

This module is mostly a collection of
parametric soids and surfaces, and supporting functions.

"""


def adaptive_arc_segments(radius, chord_error=0.1, min_segments=12, max_segments=360):
    """
    Calculate the optimal number of arc segments for a given radius
    to achieve a target maximum chord error (deviation from true circle).

    This ensures that larger objects get more segments for smooth surfaces,
    while small features don't waste polygons.

    :param radius: radius of the circle/arc in mm
    :param chord_error: maximum acceptable chord error in mm (default 0.1mm)
    :param min_segments: minimum number of segments (default 12)
    :param max_segments: maximum number of segments (default 360)
    :returns: optimal number of segments as an integer

    The chord error is the distance between the actual circle and the
    straight line segment connecting two adjacent points. For a circle
    of radius r and angle θ between segments:
        chord_error = r * (1 - cos(θ/2))

    Solving for θ:
        θ = 2 * arccos(1 - chord_error/r)
        segments = 2π / θ

    Example:
        - 5mm radius: ~31 segments (for 0.1mm error)
        - 170mm radius: ~183 segments (for 0.1mm error)
    """
    if radius <= 0 or chord_error <= 0:
        return 36  # fallback default

    # Prevent math domain error when error >= radius
    # (entire radius fits within error tolerance)
    ratio = min(chord_error / radius, 0.999999)

    # Calculate angle subtended by each segment
    theta_rad = 2 * math.acos(1 - ratio)

    # Calculate number of segments for full circle
    segments = int(math.ceil(2 * math.pi / theta_rad))

    # Clamp to reasonable bounds
    return max(min_segments, min(max_segments, segments))


def adaptive_angr_from_radius(radius, chord_error=0.1, min_angr=1.0, max_angr=10.0):
    """
    Calculate angular resolution (degrees per segment) adaptively based on radius.

    This is the inverse of adaptive_arc_segments, providing the angular resolution
    for functions that take 'angr' parameter instead of segment count.

    :param radius: radius of the circle/arc in mm
    :param chord_error: maximum acceptable chord error in mm (default 0.1mm)
    :param min_angr: minimum angular resolution in degrees (default 1.0)
    :param max_angr: maximum angular resolution in degrees (default 10.0)
    :returns: optimal angular resolution in degrees
    """
    segments = adaptive_arc_segments(radius, chord_error)
    angr = 360.0 / segments
    return max(min_angr, min(max_angr, angr))


def sphere2cartesian(lat,lon,rad):
    """
    Convert spherical polar coordinates to Cartesian coordinates for a
    sphere centered at the origin.

    :param lat: latitude in degrees
    :param lon: longitude in degrees
    :param rad: sphere radius
    :returns: ``yapcad.geom`` point in homogeneous coordinates
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

    # Attach analytic surface definition for precise operations
    from yapcad.analytic_surfaces import sphere_surface
    analytic = sphere_surface(center, rad)
    metadata = {'analytic_surface': analytic}

    return ['surface',verts,normals,faces,[],[], metadata]

# make sphere, return solid representation
def sphere(diameter,center=point(0,0,0),depth=2):
    call = f"yapcad.geom3d_util.sphere({diameter},center={center},depth={depth})"
    sld = solid([sphereSurface(diameter, center, depth)],
                [], ['procedure', call])
    if occ_available():
        try:
            from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeSphere
            from OCC.Core.gp import gp_Pnt

            c = center
            center_point = gp_Pnt(float(c[0]), float(c[1]), float(c[2]))
            shape = BRepPrimAPI_MakeSphere(center_point, float(diameter) / 2.0).Shape()
            attach_brep_to_solid(sld, BrepSolid(shape))
        except Exception:  # pragma: no cover - OCC optional
            pass
    return sld


def oblate_spheroid(equatorial_diameter, oblateness, center=point(0,0,0), depth=3):
    """Create an oblate spheroid (ellipsoid with two equal equatorial radii).

    An oblate spheroid is an ellipsoid with a < b = c (polar radius < equatorial).
    The oblateness (flattening) f = (a - c) / a where a is equatorial, c is polar.
    For Earth, f ≈ 0.00335; for Mars, f ≈ 0.00648.

    :param equatorial_diameter: diameter at the equator (X and Y axes)
    :param oblateness: geometric oblateness/flattening (0 = sphere, higher = flatter)
    :param center: center point of the spheroid
    :param depth: subdivision depth for mesh (default 2)
    :returns: yapCAD solid representing the oblate spheroid
    """
    equatorial_radius = equatorial_diameter / 2.0
    # polar_radius = equatorial_radius * (1 - oblateness)
    polar_radius = equatorial_radius * (1.0 - oblateness)

    call = f"yapcad.geom3d_util.oblate_spheroid({equatorial_diameter},{oblateness},center={center},depth={depth})"

    # Create mesh by scaling sphere mesh in Z direction
    # Start with unit sphere mesh, then scale
    unit_surf = sphereSurface(2.0, point(0,0,0), depth)  # unit sphere (radius 1)
    verts = unit_surf[1]
    normals = unit_surf[2]
    faces = unit_surf[3]

    # Scale vertices: X,Y by equatorial_radius, Z by polar_radius
    scaled_verts = []
    scaled_normals = []
    for v in verts:
        sv = [v[0] * equatorial_radius + center[0],
              v[1] * equatorial_radius + center[1],
              v[2] * polar_radius + center[2],
              1]
        scaled_verts.append(sv)

    # Normals need to be transformed by inverse transpose of scale matrix
    # For diagonal scale (sx, sy, sz), inverse transpose is (1/sx, 1/sy, 1/sz)
    for n in normals:
        sn = [n[0] / equatorial_radius,
              n[1] / equatorial_radius,
              n[2] / polar_radius,
              0]
        # Normalize
        mag_n = math.sqrt(sn[0]**2 + sn[1]**2 + sn[2]**2)
        if mag_n > 0:
            sn = [sn[0]/mag_n, sn[1]/mag_n, sn[2]/mag_n, 0]
        scaled_normals.append(sn)

    surf = ['surface', scaled_verts, scaled_normals, faces, [], []]
    sld = solid([surf], [], ['procedure', call])

    if occ_available():
        try:
            from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeSphere
            from OCC.Core.gp import gp_Pnt, gp_Trsf, gp_GTrsf
            from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_GTransform

            # Create unit sphere at origin, then apply non-uniform scale
            shape = BRepPrimAPI_MakeSphere(gp_Pnt(0, 0, 0), 1.0).Shape()

            # Apply non-uniform scaling using general transformation
            gtrsf = gp_GTrsf()
            gtrsf.SetValue(1, 1, equatorial_radius)  # X scale
            gtrsf.SetValue(2, 2, equatorial_radius)  # Y scale
            gtrsf.SetValue(3, 3, polar_radius)       # Z scale
            gtrsf.SetValue(1, 4, center[0])          # X translation
            gtrsf.SetValue(2, 4, center[1])          # Y translation
            gtrsf.SetValue(3, 4, center[2])          # Z translation

            transformer = BRepBuilderAPI_GTransform(shape, gtrsf, True)
            if transformer.IsDone():
                scaled_shape = transformer.Shape()
                attach_brep_to_solid(sld, BrepSolid(scaled_shape))
        except Exception:  # pragma: no cover - OCC optional
            pass

    return sld


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
    surf[4] = [0, 1, 2, 3]
    surf[5] = []
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
    if occ_available():
        try:
            from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
            from OCC.Core.gp import gp_Pnt
            cx, cy, cz = center[0], center[1], center[2]
            corner = gp_Pnt(cx - l2, cy - w2, cz - h2)
            shape = BRepPrimAPI_MakeBox(corner, float(length), float(width), float(height)).Shape()
            attach_brep_to_solid(sol, BrepSolid(shape))
        except Exception:
            pass
    return sol

def circleSurface(center,radius,angr=None,zup=True,chord_error=0.1):
    """make a circular surface centered at ``center`` lying in the XY
    plane with normals pointing in the positive z direction if ``zup
    == True``, negative z otherwise

    :param center: center point of circle
    :param radius: radius of circle
    :param angr: angular resolution in degrees (if None, use adaptive resolution)
    :param zup: if True, normal points +Z, else -Z
    :param chord_error: max chord error for adaptive resolution (default 0.1mm)
    """

    if angr is None:
        # Use adaptive resolution based on radius
        angr = adaptive_angr_from_radius(radius, chord_error)
    elif angr < 1 or angr > 45:
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

    surf = surface(basep,basen,basef)
    surf[4] = list(range(1, len(basep)))
    surf[5] = []
    return surf

def conic(baser,topr,height, center=point(0,0,0),angr=None,chord_error=0.1):

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
         sampling circles. If None (default), adaptive resolution is used
         based on the larger of baser and topr. Actual angular resolution
         will be ``360/round(360/angr)``

         ``chord_error`` is the maximum chord error in mm for adaptive
         resolution (default 0.1mm), ignored if angr is specified.

    """
    # Use adaptive resolution based on larger radius if angr not specified
    if angr is None:
        max_radius = max(baser, topr if topr >= epsilon else baser)
        angr = adaptive_angr_from_radius(max_radius, chord_error)

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

    baseS = circleSurface(center,baser,angr=angr,zup=False)
    baseV = baseS[1]
    ll = len(baseV)

    if not toppoint:
        topS = circleSurface(add(center,point(0,0,height)),
                             topr,angr=angr,zup=True)
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

        result = solid([baseS,cylS,topS],[],
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

        result = solid([baseS,conS],[],
                       ['procedure',call])

    if occ_available():
        try:
            brep_shape = _make_conic_brep(baser, topr, height, center)
            if brep_shape is not None:
                attach_brep_to_solid(result, BrepSolid(brep_shape))
        except Exception:
            pass

    return result

def _make_revolution_brep(contour, zStart, zEnd, steps):
    """Build a BREP solid of revolution for a simple r(z) contour."""
    if not occ_available():
        return None
    try:
        from OCC.Core.gp import gp_Pnt, gp_Ax1, gp_Dir  # type: ignore
        from OCC.Core.BRepBuilderAPI import (  # type: ignore
            BRepBuilderAPI_MakeWire,
            BRepBuilderAPI_MakeEdge,
            BRepBuilderAPI_MakeFace,
            BRepBuilderAPI_MakeSolid,
        )
        from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeRevol  # type: ignore
        from OCC.Core.TopExp import TopExp_Explorer  # type: ignore
        from OCC.Core.TopAbs import TopAbs_SOLID  # type: ignore
        from OCC.Core import TopoDS  # type: ignore
        from OCC.Core.TopoDS import topods  # type: ignore
    except Exception:
        return None

    def _safe_radius(z):
        try:
            r = float(contour(z))
        except Exception:
            return None
        return max(abs(r), epsilon)

    zs = [zStart + (zEnd - zStart) * i / float(max(1, steps)) for i in range(steps + 1)]
    samples = []
    for z in zs:
        r = _safe_radius(z)
        if r is None:
            return None
        samples.append((r, z))

    # Build a closed polyline: profile along r(z), then axis back to start.
    pts = [gp_Pnt(r, 0.0, z) for (r, z) in samples]
    axis_top = gp_Pnt(0.0, 0.0, samples[-1][1])
    axis_bottom = gp_Pnt(0.0, 0.0, samples[0][1])

    wire_maker = BRepBuilderAPI_MakeWire()
    for a, b in zip(pts[:-1], pts[1:]):
        wire_maker.Add(BRepBuilderAPI_MakeEdge(a, b).Edge())
    wire_maker.Add(BRepBuilderAPI_MakeEdge(pts[-1], axis_top).Edge())
    wire_maker.Add(BRepBuilderAPI_MakeEdge(axis_top, axis_bottom).Edge())
    wire_maker.Add(BRepBuilderAPI_MakeEdge(axis_bottom, pts[0]).Edge())
    wire = wire_maker.Wire()

    face_builder = BRepBuilderAPI_MakeFace(wire)
    if not face_builder.IsDone():
        return None
    face = face_builder.Face()

    axis = gp_Ax1(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(0.0, 0.0, 1.0))
    revol = BRepPrimAPI_MakeRevol(face, axis, pi2)
    shape = revol.Shape()

    # Prefer an explicit solid if present
    exp = TopExp_Explorer(shape, TopAbs_SOLID)
    if exp.More():
        return topods.Solid(exp.Current())

    # If we got a shell, try to promote it to a solid
    try:
        shell = TopoDS.topods_Shell(shape)
        solid_builder = BRepBuilderAPI_MakeSolid(shell)
        if solid_builder.IsDone():
            return solid_builder.Solid()
    except Exception:
        pass

    try:
        return topods.Solid(shape)
    except Exception:
        return None


def makeRevolutionSurface(contour,zStart,zEnd,steps,arcSamples=None,*,chord_error=0.1,return_brep=False):
    """
    Generate a surface of revolution by sampling a contour function.

    :param contour: callable mapping ``z`` to a radial distance
    :param zStart: lower bound for the ``z`` interval
    :param zEnd: upper bound for the ``z`` interval
    :param steps: number of contour samples between ``zStart`` and ``zEnd``
    :param arcSamples: number of samples around the revolution arc.
                       If None (default), adaptive resolution is used.
    :param chord_error: maximum chord error for adaptive resolution (default 0.1mm)
    :returns: ``['surface', vertices, normals, faces]`` list representing the surface
    """

    sV=[]
    sN=[]
    sF=[]
    zRange = zEnd-zStart
    zD = zRange/steps

    # Use adaptive resolution if arcSamples not specified
    if arcSamples is None:
        # Sample the contour to find maximum radius
        max_radius = 0
        for i in range(steps + 1):
            z = i * zD + zStart
            r = contour(z)
            max_radius = max(max_radius, r)
        arcSamples = adaptive_arc_segments(max_radius, chord_error)

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
                    _, n = tri2p0n([sV[start_pole_idx], pp2, pp1])
                except ValueError:
                    continue

                k1, sV, sN = addVertex(pp1, n, sV, sN)
                k2, sV, sN = addVertex(pp2, n, sV, sN)
                sF.append([start_pole_idx, k2, k1])
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
                    _, n = tri2p0n([p1, p2, sV[end_pole_idx]])
                except ValueError:
                    continue

                k1, sV, sN = addVertex(p1, n, sV, sN)
                k2, sV, sN = addVertex(p2, n, sV, sN)
                sF.append([k1, k2, end_pole_idx])
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
        
    surf = surface(sV,sN,sF)

    if not return_brep:
        return surf

    brep_shape = None
    try:
        brep_shape = _make_revolution_brep(contour, zStart, zEnd, steps)
    except Exception:
        brep_shape = None
    return surf, brep_shape

def makeRevolutionThetaSamplingSurface(contour, zStart, zEnd, arcSamples=360,
                                       endcaps=False, degrees=True, *,
                                       return_brep=False):
    """Generate a surface of revolution using a theta-dependent contour.

    If ``return_brep`` is True and the contour is axisymmetric (identical for
    all theta and ``wrap_shift == 0``), a native BREP solid is also returned;
    otherwise the second return value is ``None``.
    """

    sV = []
    sN = []
    sF = []
    zRange = zEnd - zStart
    if zRange <= 0:
        raise ValueError('zEnd must be greater than zStart')

    base_theta = 0.0 if degrees else 0.0
    profile_zero_raw = contour(zStart, zEnd, base_theta)
    profile_zero, wrap_shift = _normalize_theta_profile(profile_zero_raw)
    steps = len(profile_zero)
    if steps < 2:
        raise ValueError('contour function must return at least two samples')
    wrap_shift = int(max(0, min(wrap_shift, steps - 1)))

    arcSamples = max(3, int(arcSamples))
    degStep = 360.0 / arcSamples
    radStep = pi2 / arcSamples

    angle_cos = []
    angle_sin = []
    for i in range(arcSamples):
        theta_val = degStep * i if degrees else radStep * i
        angle_cos.append(math.cos(math.radians(theta_val) if degrees else theta_val))
        angle_sin.append(math.sin(math.radians(theta_val) if degrees else theta_val))
    angle_cos[0] = 1.0
    angle_sin[0] = 0.0

    need_start_cap = False
    need_end_cap = False
    start_pole_idx = None
    end_pole_idx = None

    if endcaps:
        r_start = profile_zero[0][1]
        r_end = profile_zero[-1][1]
        need_start_cap = r_start < epsilon * 10
        need_end_cap = r_end < epsilon * 10
        if need_start_cap:
            pole_point = [0.0, 0.0, profile_zero[0][0], 1.0]
            pole_normal = [0.0, 0.0, -1.0, 0.0]
            start_pole_idx, sV, sN = addVertex(pole_point, pole_normal, sV, sN)
        if need_end_cap:
            pole_point = [0.0, 0.0, profile_zero[-1][0], 1.0]
            pole_normal = [0.0, 0.0, 1.0, 0.0]
            end_pole_idx, sV, sN = addVertex(pole_point, pole_normal, sV, sN)

    profiles = [profile_zero]
    axisymmetric = True

    for i in range(1, arcSamples):
        theta_val = degStep * i if degrees else radStep * i
        prof_raw = contour(zStart, zEnd, theta_val)
        prof, wrap_val = _normalize_theta_profile(prof_raw)
        if len(prof) != steps:
            raise ValueError('contour returned inconsistent sample counts for theta sweep')
        if wrap_val and wrap_shift == 0:
            wrap_shift = int(max(0, min(wrap_val, steps - 1)))
        if axisymmetric and wrap_val != wrap_shift:
            axisymmetric = False
        if axisymmetric:
            for a, b in zip(prof, profile_zero):
                if abs(a[0] - b[0]) > 1e-9 or abs(a[1] - b[1]) > 1e-6:
                    axisymmetric = False
                    break
        profiles.append(prof)

    for ang_idx in range(arcSamples):
        next_idx = (ang_idx + 1) % arcSamples
        cos_a = angle_cos[ang_idx]
        sin_a = angle_sin[ang_idx]
        cos_b = angle_cos[next_idx]
        sin_b = angle_sin[next_idx]
        prof_a = profiles[ang_idx]
        prof_b = profiles[next_idx]
        apply_shift = wrap_shift if (wrap_shift > 0 and next_idx == 0) else 0

        for j in range(steps - 1):
            z0, r0 = prof_a[j]
            z1, r1 = prof_a[j + 1]

            idx_b0 = j
            idx_b1 = j + 1
            if apply_shift:
                idx_b0 = j + apply_shift
                idx_b1 = j + 1 + apply_shift
                if idx_b1 >= steps:
                    continue

            z2, r2 = prof_b[idx_b1]
            z3, r3 = prof_b[idx_b0]

            v0 = [cos_a * r0, sin_a * r0, z0, 1.0]
            v1 = [cos_a * r1, sin_a * r1, z1, 1.0]
            v2 = [cos_b * r2, sin_b * r2, z2, 1.0]
            v3 = [cos_b * r3, sin_b * r3, z3, 1.0]

            try:
                _, n0 = tri2p0n([v1, v0, v2])
            except ValueError:
                continue

            k0, sV, sN = addVertex(v0, n0, sV, sN)
            k1, sV, sN = addVertex(v1, n0, sV, sN)
            k2, sV, sN = addVertex(v2, n0, sV, sN)
            k3, sV, sN = addVertex(v3, n0, sV, sN)

            sF.append([k0, k2, k1])
            sF.append([k0, k3, k2])

    if endcaps:
        if need_start_cap and start_pole_idx is not None:
            for ang_idx in range(arcSamples):
                next_idx = (ang_idx + 1) % arcSamples
                prof = profiles[ang_idx]
                prof_next = profiles[next_idx]
                r_curr = prof[0][1]
                r_next = prof_next[0][1]
                v_curr = [angle_cos[ang_idx] * r_curr, angle_sin[ang_idx] * r_curr, prof[0][0], 1.0]
                v_next = [angle_cos[next_idx] * r_next, angle_sin[next_idx] * r_next, prof_next[0][0], 1.0]
                try:
                    _, n_cap = tri2p0n([v_next, v_curr, sV[start_pole_idx]])
                except ValueError:
                    continue
                i_curr, sV, sN = addVertex(v_curr, n_cap, sV, sN)
                i_next, sV, sN = addVertex(v_next, n_cap, sV, sN)
                sF.append([start_pole_idx, i_next, i_curr])
        if need_end_cap and end_pole_idx is not None:
            for ang_idx in range(arcSamples):
                next_idx = (ang_idx + 1) % arcSamples
                prof = profiles[ang_idx]
                prof_next = profiles[next_idx]
                r_curr = prof[-1][1]
                r_next = prof_next[-1][1]
                v_curr = [angle_cos[ang_idx] * r_curr, angle_sin[ang_idx] * r_curr, prof[-1][0], 1.0]
                v_next = [angle_cos[next_idx] * r_next, angle_sin[next_idx] * r_next, prof_next[-1][0], 1.0]
                try:
                    _, n_cap = tri2p0n([v_curr, v_next, sV[end_pole_idx]])
                except ValueError:
                    continue
                i_curr, sV, sN = addVertex(v_curr, n_cap, sV, sN)
                i_next, sV, sN = addVertex(v_next, n_cap, sV, sN)
                sF.append([end_pole_idx, i_curr, i_next])

    surf = surface(sV, sN, sF)

    if not return_brep:
        return surf

    brep_shape = None
    if wrap_shift == 0 and axisymmetric:
        def _contour_r(z_val):
            # linear interpolate over profile_zero
            if z_val <= profile_zero[0][0]:
                return profile_zero[0][1]
            if z_val >= profile_zero[-1][0]:
                return profile_zero[-1][1]
            for idx in range(len(profile_zero) - 1):
                z0, r0 = profile_zero[idx]
                z1, r1 = profile_zero[idx + 1]
                if z0 <= z_val <= z1 or z1 <= z_val <= z0:
                    t = (z_val - z0) / (z1 - z0) if abs(z1 - z0) > 1e-12 else 0.0
                    return r0 + t * (r1 - r0)
            return profile_zero[-1][1]

        try:
            brep_shape = _make_revolution_brep(_contour_r, zStart, zEnd, steps)
        except Exception:
            brep_shape = None

    return surf, brep_shape


def makeRevolutionSolid(contour, zStart, zEnd, steps, arcSamples=None, chord_error=0.1, metadata=None):
    """
    Build a solid of revolution around the Z axis. When pythonocc-core is
    available, a native BREP is attached; otherwise we fall back to the
    tessellated representation.

    :param arcSamples: number of samples around the revolution arc.
                       If None (default), adaptive resolution is used.
    :param chord_error: maximum chord error for adaptive resolution (default 0.1mm)
    """
    surf, brep_shape = makeRevolutionSurface(
        contour, zStart, zEnd, steps, arcSamples=arcSamples, chord_error=chord_error, return_brep=True
    )
    call = f"yapcad.geom3d_util.makeRevolutionSolid(contour,{zStart},{zEnd},{steps},{arcSamples})"
    construction = ['procedure', call]
    if metadata is not None and isinstance(metadata, dict):
        sld = solid([surf], [], construction, metadata)
    else:
        sld = solid([surf], [], construction)
    if brep_shape is not None:
        try:
            attach_brep_to_solid(sld, BrepSolid(brep_shape))
        except Exception:
            pass
    return sld


def _normalize_theta_profile(samples):
    wrap = 0
    data = samples
    if isinstance(samples, tuple):
        if len(samples) != 2:
            raise ValueError('contour tuple must be (points, wrap)')
        data, wrap = samples
    if not isinstance(data, (list, tuple)) or not data:
        raise ValueError('contour must return a sequence of samples')

    normalized = []
    for sample in data:
        if ispoint(sample):
            z_val = float(sample[0])
            r_val = abs(float(sample[1]))
        elif isinstance(sample, (list, tuple)) and len(sample) >= 2:
            z_val = float(sample[0])
            r_val = abs(float(sample[1]))
        else:
            raise ValueError('invalid contour sample; expected point or [x,y] pair')
        normalized.append((z_val, r_val))
    return normalized, int(wrap)


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
    if stripF:
        strip = surface(stripV, stripN, stripF)
    else:
        strip = ['surface', [], [], [], [], []]

    result = solid([s2,strip,s1],
                   [],
                   ['procedure',call])
    if occ_available():
        try:
            brep_shape = _extrude_brep_shape(s2[1], loops, distance, direction)
            if brep_shape is not None:
                attach_brep_to_solid(result, BrepSolid(brep_shape))
        except Exception:
            pass
    return result


def _extrude_brep_shape(vertices, loops, distance, direction):
    if not loops:
        return None
    try:
        from OCC.Core.BRepBuilderAPI import (
            BRepBuilderAPI_MakeWire,
            BRepBuilderAPI_MakeEdge,
            BRepBuilderAPI_MakeFace,
        )
        from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakePrism
        from OCC.Core.gp import gp_Pnt, gp_Vec
    except ImportError:
        return None

    def _loop_points(loop):
        return [
            (float(vertices[idx][0]), float(vertices[idx][1]), float(vertices[idx][2]))
            for idx in loop
        ]

    def _wire_for_points(points):
        if len(points) < 2:
            return None
        writer = BRepBuilderAPI_MakeWire()
        for i in range(len(points)):
            p0 = gp_Pnt(*points[i])
            p1 = gp_Pnt(*points[(i + 1) % len(points)])
            if p0.Distance(p1) <= epsilon:
                continue
            edge = BRepBuilderAPI_MakeEdge(p0, p1).Edge()
            writer.Add(edge)
        return writer.Wire()

    outer = _wire_for_points(_loop_points(loops[0]))
    if outer is None:
        return None
    face_builder = BRepBuilderAPI_MakeFace(outer, False)
    for hole in loops[1:]:
        wire = _wire_for_points(_loop_points(hole))
        if wire is not None:
            face_builder.Add(wire)
    face = face_builder.Face()
    dir_vec = scale4(direction, distance)
    prism_vec = gp_Vec(float(dir_vec[0]), float(dir_vec[1]), float(dir_vec[2]))
    prism = BRepPrimAPI_MakePrism(face, prism_vec, True).Shape()
    try:
        from OCC.Core.TopoDS import topods
        prism = topods.Solid(prism)
    except Exception:
        pass
    return prism


def _make_conic_brep(baser, topr, height, center):
    try:
        from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeCone, BRepPrimAPI_MakeCylinder
        from OCC.Core.gp import gp_Ax2, gp_Pnt, gp_Dir
    except ImportError:
        return None
    axis = gp_Ax2(gp_Pnt(float(center[0]), float(center[1]), float(center[2])),
                  gp_Dir(0.0, 0.0, 1.0))
    if abs(baser - topr) < epsilon:
        if baser < epsilon:
            return None
        return BRepPrimAPI_MakeCylinder(axis, float(baser), float(height)).Shape()
    return BRepPrimAPI_MakeCone(axis, float(baser), float(topr), float(height)).Shape()


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
    def _dedupe(loop):
        cleaned = []
        for pt in loop:
            if not cleaned or mag(sub(pt, cleaned[-1])) > epsilon:
                cleaned.append(point(pt))
        return cleaned

    lower = _dedupe(lower_loop)
    upper = _dedupe(upper_loop)
    if len(lower) != len(upper) or len(lower) < 3:
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

    if not faces:
        raise ValueError('loft surface produced no faces; check input loops')

    return surface(vertices, normals, faces)


def _loft_brep(lower_loop, upper_loop):
    if not occ_available():
        return None
    try:
        from OCC.Core.gp import gp_Pnt  # type: ignore
        from OCC.Core.BRepBuilderAPI import (  # type: ignore
            BRepBuilderAPI_MakeWire,
            BRepBuilderAPI_MakeEdge,
            BRepBuilderAPI_MakeFace,
        )
        from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections  # type: ignore
        from OCC.Core.TopoDS import topods  # type: ignore
    except Exception:
        return None

    def _wire_from_loop(loop):
        builder = BRepBuilderAPI_MakeWire()
        for i in range(len(loop)):
            p0 = loop[i]
            p1 = loop[(i + 1) % len(loop)]
            e = BRepBuilderAPI_MakeEdge(
                gp_Pnt(float(p0[0]), float(p0[1]), float(p0[2])),
                gp_Pnt(float(p1[0]), float(p1[1]), float(p1[2])),
            ).Edge()
            builder.Add(e)
        return builder.Wire()

    w1 = _wire_from_loop(lower_loop)
    w2 = _wire_from_loop(upper_loop)
    loft = BRepOffsetAPI_ThruSections(True, False, 1.0e-6)
    loft.AddWire(w1)
    loft.AddWire(w2)
    loft.Build()
    if not loft.IsDone():
        return None
    shape = loft.Shape()
    try:
        return topods.Solid(shape)
    except Exception:
        return shape


def makeLoftSolid(lower_loop, upper_loop, *, metadata=None):
    """
    Create a solid loft between two planar loops (matching vertex counts).
    Attempts to attach a native BREP when pythonocc-core is available,
    otherwise falls back to the tessellated representation.
    """
    surf = _loft_surface(lower_loop, upper_loop)
    call = f"yapcad.geom3d_util.makeLoftSolid(lower_loop, upper_loop)"
    construction = ['procedure', call]
    if metadata is not None and isinstance(metadata, dict):
        sld = solid([surf], [], construction, metadata)
    else:
        sld = solid([surf], [], construction)

    try:
        brep_shape = _loft_brep(lower_loop, upper_loop)
        if brep_shape is not None:
            attach_brep_to_solid(sld, BrepSolid(brep_shape))
    except Exception:
        pass

    return sld


        

    

def _circle_loop(center_xy, radius, minang=None, chord_error=0.1):
    """Generate a circle loop with adaptive or fixed resolution.

    :param center_xy: (x, y) center coordinates
    :param radius: radius of circle
    :param minang: minimum angular resolution in degrees. If None, use adaptive.
    :param chord_error: maximum chord error for adaptive resolution (default 0.1mm)
    """
    if minang is None:
        minang = adaptive_angr_from_radius(radius, chord_error)

    arc_geom = [arc(point(center_xy[0], center_xy[1]), radius)]
    loop = geomlist2poly(arc_geom, minang=minang, minlen=0.0)
    if not loop:
        raise ValueError('failed to generate circle loop')
    return loop


def tube(outer_diameter, wall_thickness, length,
         center=None, *, base_point=None, minang=None, chord_error=0.1, include_caps=True):
    """Create a cylindrical tube solid.

    ``base_point`` (or legacy ``center`` argument) identifies the base of the
    cylindrical wall, i.e. the plane where ``z == base_point[2]``.

    :param minang: minimum angular resolution in degrees. If None, use adaptive.
    :param chord_error: maximum chord error for adaptive resolution (default 0.1mm)
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

    # Use adaptive resolution based on outer radius if not specified
    if minang is None:
        minang = adaptive_angr_from_radius(outer_radius, chord_error)

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
    result = solid(surfaces,
                   [],
                   ['procedure', call])
    if occ_available():
        try:
            from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeCylinder
            from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut
            from OCC.Core.gp import gp_Pnt, gp_Dir, gp_Ax2

            axis = gp_Ax2(gp_Pnt(float(base_point[0]), float(base_point[1]), float(base_point[2])),
                          gp_Dir(0.0, 0.0, 1.0))
            cyl_outer = BRepPrimAPI_MakeCylinder(axis, float(outer_radius), float(length)).Shape()
            cyl_inner = BRepPrimAPI_MakeCylinder(axis, float(inner_radius), float(length)).Shape()
            brep_shape = BRepAlgoAPI_Cut(cyl_outer, cyl_inner).Shape()
            attach_brep_to_solid(result, BrepSolid(brep_shape))
        except Exception:
            pass
    return result


def conic_tube(bottom_outer_diameter, top_outer_diameter, wall_thickness,
               length, center=None, *, base_point=None, minang=None, chord_error=0.1, include_caps=True):
    """Create a conic tube with varying outer diameter.

    ``base_point`` (or ``center`` legacy argument) marks the axial base of the
    frustum (the larger-diameter end when stacked).

    :param minang: minimum angular resolution in degrees. If None, use adaptive.
    :param chord_error: maximum chord error for adaptive resolution (default 0.1mm)
    """

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

    # Use adaptive resolution based on larger outer radius if not specified
    if minang is None:
        max_radius = max(r0_outer, r1_outer)
        minang = adaptive_angr_from_radius(max_radius, chord_error)

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
    result = solid(surfaces,
                   [],
                   ['procedure', call])
    if occ_available():
        try:
            from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeCone
            from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut
            from OCC.Core.gp import gp_Ax2, gp_Pnt, gp_Dir
            axis = gp_Ax2(gp_Pnt(float(base_point[0]), float(base_point[1]), float(base_point[2])),
                          gp_Dir(0.0, 0.0, 1.0))
            outer = BRepPrimAPI_MakeCone(axis, float(r0_outer), float(r1_outer), float(length)).Shape()
            inner = BRepPrimAPI_MakeCone(axis, float(r0_inner), float(r1_inner), float(length)).Shape()
            brep_shape = BRepAlgoAPI_Cut(outer, inner).Shape()
            attach_brep_to_solid(result, BrepSolid(brep_shape))
        except Exception:
            pass
    return result


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
    result = solid(surfaces,
                   [],
                   ['procedure', call])

    if occ_available():
        try:
            from OCC.Core.gp import gp_Ax2, gp_Pnt, gp_Dir  # type: ignore
            from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeSphere  # type: ignore
            from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut  # type: ignore
            axis = gp_Ax2(gp_Pnt(float(center[0]), float(center[1]), float(center[2])),
                          gp_Dir(0.0, 0.0, 1.0))
            outer_shape = BRepPrimAPI_MakeSphere(axis, float(outer_radius), float(theta)).Shape()
            inner_shape = BRepPrimAPI_MakeSphere(axis, float(inner_radius), float(theta)).Shape()
            brep_shape = BRepAlgoAPI_Cut(outer_shape, inner_shape).Shape()
            attach_brep_to_solid(result, BrepSolid(brep_shape))
        except Exception:
            pass

    return result


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


# ============================================================================
# Path3D Sampling and Adaptive Sweep Functions
# ============================================================================

def _normalize_vector(v):
    """Normalize a 3D vector to unit length."""
    import math
    mag = math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    if mag < 1e-10:
        return [0, 0, 1]  # Default to Z-up for degenerate case
    return [v[0]/mag, v[1]/mag, v[2]/mag]


def _cross_product(a, b):
    """Compute cross product of two 3D vectors."""
    return [
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    ]


def _dot_product(a, b):
    """Compute dot product of two 3D vectors."""
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]


def _vector_subtract(a, b):
    """Subtract two 3D vectors."""
    return [a[0]-b[0], a[1]-b[1], a[2]-b[2]]


def _vector_add(a, b):
    """Add two 3D vectors."""
    return [a[0]+b[0], a[1]+b[1], a[2]+b[2]]


def _vector_scale(v, s):
    """Scale a 3D vector."""
    return [v[0]*s, v[1]*s, v[2]*s]


def _lerp_point(p1, p2, t):
    """Linear interpolation between two points."""
    return [
        p1[0] + t*(p2[0]-p1[0]),
        p1[1] + t*(p2[1]-p1[1]),
        p1[2] + t*(p2[2]-p1[2])
    ]


def _angle_between_vectors(v1, v2):
    """Compute angle in degrees between two unit vectors."""
    import math
    dot = _dot_product(v1, v2)
    dot = max(-1.0, min(1.0, dot))  # Clamp for numerical stability
    return math.degrees(math.acos(dot))


def _compute_line_tangent(start, end):
    """Compute unit tangent for a line segment."""
    direction = _vector_subtract(end, start)
    return _normalize_vector(direction)


def _compute_arc_tangent(center, point, normal):
    """Compute unit tangent for an arc at a given point.

    The tangent is perpendicular to the radius and lies in the arc plane.
    """
    # Radius vector from center to point
    radius_vec = _vector_subtract(point, center)
    # Tangent is perpendicular to radius, in the plane defined by normal
    # tangent = normal × radius (gives direction along arc)
    tangent = _cross_product(normal, radius_vec)
    return _normalize_vector(tangent)


def _sample_line_segment(start, end, num_samples):
    """Sample a line segment uniformly.

    Returns list of (point, tangent, parameter) tuples.
    """
    tangent = _compute_line_tangent(start, end)
    samples = []
    for i in range(num_samples + 1):
        t = i / num_samples
        point = _lerp_point(start, end, t)
        samples.append((point, tangent, t))
    return samples


def _sample_arc_segment(center, start, end, normal, num_samples):
    """Sample an arc segment uniformly.

    Returns list of (point, tangent, parameter) tuples.
    """
    import math

    # Compute arc angle
    r_start = _vector_subtract(start, center)
    r_end = _vector_subtract(end, center)
    radius = math.sqrt(_dot_product(r_start, r_start))

    # Normalize radius vectors
    r_start_n = _normalize_vector(r_start)
    r_end_n = _normalize_vector(r_end)

    # Compute angle between start and end
    dot = _dot_product(r_start_n, r_end_n)
    dot = max(-1.0, min(1.0, dot))
    arc_angle = math.acos(dot)

    # Determine rotation direction using cross product with normal
    cross = _cross_product(r_start_n, r_end_n)
    if _dot_product(cross, normal) < 0:
        arc_angle = 2*math.pi - arc_angle

    samples = []
    for i in range(num_samples + 1):
        t = i / num_samples
        angle = t * arc_angle

        # Rotate r_start around normal by angle using Rodrigues' formula
        cos_a = math.cos(angle)
        sin_a = math.sin(angle)
        r_rot = _vector_add(
            _vector_scale(r_start_n, cos_a),
            _vector_add(
                _vector_scale(_cross_product(normal, r_start_n), sin_a),
                _vector_scale(normal, _dot_product(normal, r_start_n) * (1 - cos_a))
            )
        )

        point = _vector_add(center, _vector_scale(r_rot, radius))
        tangent = _compute_arc_tangent(center, point, normal)
        samples.append((point, tangent, t))

    return samples


def _sample_path3d(path3d, samples_per_segment=20):
    """Sample a path3d at regular intervals.

    Args:
        path3d: Dict with 'segments' list of line/arc segments
        samples_per_segment: Number of samples per segment

    Returns:
        List of (point, tangent, global_param) tuples where global_param
        is in [0, 1] across the entire path.
    """
    segments = path3d.get('segments', [])
    if not segments:
        return []

    all_samples = []
    total_segments = len(segments)

    for seg_idx, seg in enumerate(segments):
        seg_type = seg.get('type', 'line')

        if seg_type == 'line':
            start = seg['start']
            end = seg['end']
            seg_samples = _sample_line_segment(start, end, samples_per_segment)
        elif seg_type == 'arc':
            center = seg['center']
            start = seg['start']
            end = seg['end']
            normal = seg.get('normal', [0, 0, 1])
            seg_samples = _sample_arc_segment(center, start, end, normal, samples_per_segment)
        else:
            continue

        # Adjust parameters to global range
        for point, tangent, local_t in seg_samples:
            # Skip first point of non-first segments to avoid duplicates
            if seg_idx > 0 and local_t == 0:
                continue
            global_t = (seg_idx + local_t) / total_segments
            all_samples.append((point, tangent, global_t))

    return all_samples


def _adaptive_sample_path3d(path3d, angle_threshold_deg=5.0, samples_per_segment=50):
    """Sample path adaptively, emitting samples when tangent changes exceed threshold.

    Args:
        path3d: Dict with 'segments' list
        angle_threshold_deg: Angle change in degrees that triggers new sample
        samples_per_segment: Dense sampling rate for angle detection

    Returns:
        List of (point, tangent, global_param) tuples at profile locations.
    """
    # Get dense samples
    dense_samples = _sample_path3d(path3d, samples_per_segment)
    if len(dense_samples) < 2:
        return dense_samples

    # Always include start
    result = [dense_samples[0]]
    last_emitted_tangent = dense_samples[0][1]

    for i in range(1, len(dense_samples) - 1):
        point, tangent, param = dense_samples[i]
        angle = _angle_between_vectors(last_emitted_tangent, tangent)

        if angle >= angle_threshold_deg:
            result.append((point, tangent, param))
            last_emitted_tangent = tangent

    # Always include end
    result.append(dense_samples[-1])

    return result


def _slerp(v1, v2, t):
    """Spherical linear interpolation between two unit vectors."""
    import math

    dot = _dot_product(v1, v2)
    dot = max(-1.0, min(1.0, dot))

    # If vectors are very close, use linear interpolation
    if dot > 0.9995:
        result = _vector_add(
            _vector_scale(v1, 1 - t),
            _vector_scale(v2, t)
        )
        return _normalize_vector(result)

    theta = math.acos(dot)
    sin_theta = math.sin(theta)

    s1 = math.sin((1 - t) * theta) / sin_theta
    s2 = math.sin(t * theta) / sin_theta

    return _vector_add(
        _vector_scale(v1, s1),
        _vector_scale(v2, s2)
    )


def _interpolate_up_samples(up_samples, param):
    """Interpolate user-provided up vectors at a given parameter.

    Args:
        up_samples: List of [t, [ux, uy, uz]] pairs, sorted by t
        param: Parameter in [0, 1] to interpolate at

    Returns:
        Normalized up vector at param
    """
    if not up_samples:
        return [0, 0, 1]  # Default Z-up

    # Sort by parameter
    sorted_samples = sorted(up_samples, key=lambda x: x[0])

    # Find bracketing samples
    for i, (t, v) in enumerate(sorted_samples):
        if t >= param:
            if i == 0:
                return _normalize_vector(v)
            # Interpolate between i-1 and i
            t0, v0 = sorted_samples[i-1]
            t1, v1 = sorted_samples[i]
            if abs(t1 - t0) < 1e-10:
                return _normalize_vector(v0)
            local_t = (param - t0) / (t1 - t0)
            return _slerp(_normalize_vector(v0), _normalize_vector(v1), local_t)

    # Past end, use last sample
    return _normalize_vector(sorted_samples[-1][1])


def _compute_minimal_twist_frame(tangent, prev_up):
    """Compute a coordinate frame that minimizes twist from previous frame.

    Args:
        tangent: Unit tangent vector (profile normal direction)
        prev_up: Previous frame's up vector

    Returns:
        (right, up) unit vectors forming a frame with tangent
    """
    # Project prev_up onto plane perpendicular to tangent
    dot = _dot_product(prev_up, tangent)
    up = _vector_subtract(prev_up, _vector_scale(tangent, dot))

    # Handle degenerate case (prev_up parallel to tangent)
    mag = (_dot_product(up, up)) ** 0.5
    if mag < 1e-10:
        # Pick arbitrary perpendicular
        if abs(tangent[2]) < 0.9:
            up = _cross_product(tangent, [0, 0, 1])
        else:
            up = _cross_product(tangent, [1, 0, 0])
        up = _normalize_vector(up)
    else:
        up = _vector_scale(up, 1.0/mag)

    right = _cross_product(tangent, up)
    return right, up


def _compute_frenet_frame(tangent, tangent_derivative):
    """Compute Frenet frame from tangent and its derivative.

    Args:
        tangent: Unit tangent vector
        tangent_derivative: Rate of change of tangent (curvature direction)

    Returns:
        (right, up) unit vectors, or None if degenerate
    """
    # Binormal = T × T'
    binormal = _cross_product(tangent, tangent_derivative)
    mag = (_dot_product(binormal, binormal)) ** 0.5

    if mag < 1e-10:
        return None  # Degenerate (straight segment)

    binormal = _vector_scale(binormal, 1.0/mag)
    # Normal = B × T
    normal = _cross_product(binormal, tangent)

    return binormal, normal


def sweep_adaptive(profile, spine, *, inner_profiles=None,
                   angle_threshold_deg=5.0,
                   frame_mode='minimal_twist',
                   up_samples=None,
                   ruled=True,
                   metadata=None):
    """Sweep a profile along a path with adaptive tangent tracking.

    The profile normal tracks the path tangent. New profile sections are
    generated whenever the tangent direction changes by more than the
    threshold angle. Uses BRepOffsetAPI_ThruSections to loft between sections.

    Args:
        profile: yapCAD region2d for outer boundary (in XY plane, centered at origin)
        spine: path3d dict with line/arc segments
        inner_profiles: Optional region2d or list of region2d for inner voids.
                        Supports multiple voids (e.g., split pipe for heat exchanger).
                        Creates hollow solids by boolean subtraction.
        angle_threshold_deg: Angle change (degrees) that triggers new section (default 5.0)
        frame_mode: One of:
            - 'minimal_twist': Profile 'up' stays consistent (default)
            - 'frenet': Natural Frenet frame (may twist at inflections)
            - 'custom': Use up_samples for interpolated orientation
        up_samples: For frame_mode='custom', list of [t, [ux, uy, uz]] pairs
                    where t in [0,1] is path parameter and [ux,uy,uz] is up vector.
                    Vectors are normalized automatically.
        ruled: If True (default), use ruled surfaces that preserve straight edges.
               If False, use smooth interpolation between sections.
        metadata: Optional metadata dict

    Returns:
        yapCAD solid created by lofting between adapted sections
    """
    if not occ_available():
        raise RuntimeError("sweep_adaptive requires pythonocc-core")

    from OCC.Core.gp import gp_Pnt, gp_Dir, gp_Vec, gp_Ax2, gp_Trsf, gp_Circ
    from OCC.Core.BRepBuilderAPI import (
        BRepBuilderAPI_MakeWire, BRepBuilderAPI_MakeEdge,
        BRepBuilderAPI_MakeFace, BRepBuilderAPI_Transform
    )
    from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections
    from OCC.Core.GC import GC_MakeArcOfCircle
    from OCC.Core.GProp import GProp_GProps
    from OCC.Core.BRepGProp import brepgprop
    from OCC.Core.TopoDS import topods
    import math

    def _build_wire_from_region2d(region, reverse=False):
        """Build an OCC wire from a yapCAD region2d in the XZ plane."""
        wire_builder = BRepBuilderAPI_MakeWire()
        segments = list(region)
        if reverse:
            segments = segments[::-1]

        for seg in segments:
            if isline(seg):
                p1 = seg[0]
                p2 = seg[1]
                if reverse:
                    p1, p2 = p2, p1
                # Profile in XZ plane (Y=0)
                start = gp_Pnt(p1[0], 0, p1[1] if len(p1) > 1 else 0)
                end = gp_Pnt(p2[0], 0, p2[1] if len(p2) > 1 else 0)
                edge = BRepBuilderAPI_MakeEdge(start, end).Edge()
                wire_builder.Add(edge)
            elif isarc(seg):
                center = seg[0]
                params = seg[1]
                radius = params[0]
                start_ang = math.radians(params[1])
                end_ang = math.radians(params[2])
                if reverse:
                    start_ang, end_ang = end_ang, start_ang
                cx, cy = center[0], center[1] if len(center) > 1 else 0
                start_pt = gp_Pnt(cx + radius * math.cos(start_ang), 0, cy + radius * math.sin(start_ang))
                end_pt = gp_Pnt(cx + radius * math.cos(end_ang), 0, cy + radius * math.sin(end_ang))
                arc_center = gp_Pnt(cx, 0, cy)
                circ = gp_Circ(gp_Ax2(arc_center, gp_Dir(0, 1, 0)), radius)
                arc_maker = GC_MakeArcOfCircle(circ, start_pt, end_pt, True)
                if arc_maker.IsDone():
                    edge = BRepBuilderAPI_MakeEdge(arc_maker.Value()).Edge()
                    wire_builder.Add(edge)
        return wire_builder.Wire()

    def _transform_wire_to_frame(wire, position, tangent, right, up):
        """Transform wire from XZ plane to specified coordinate frame.

        Wire is built in XZ plane with:
        - X axis = profile X (will become 'right')
        - Z axis = profile Y (will become 'up')
        - Y axis = profile normal (will become 'tangent'/forward)
        """
        from OCC.Core.gp import gp_Ax3

        # Build transformation matrix
        # From: X=(1,0,0), Y=(0,1,0), Z=(0,0,1)
        # To: X=right, Y=tangent, Z=up
        trsf = gp_Trsf()

        # Set rotation part using Ax3 (SetDisplacement requires Ax3)
        ax_from = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 1, 0), gp_Dir(1, 0, 0))  # Y forward, X right
        ax_to = gp_Ax3(
            gp_Pnt(position[0], position[1], position[2]),
            gp_Dir(tangent[0], tangent[1], tangent[2]),
            gp_Dir(right[0], right[1], right[2])
        )
        trsf.SetDisplacement(ax_from, ax_to)

        transformer = BRepBuilderAPI_Transform(wire, trsf, True)
        return topods.Wire(transformer.Shape())

    # Normalize inner_profiles to list
    if inner_profiles is None:
        inner_profiles_list = []
    elif isinstance(inner_profiles, list) and inner_profiles and isinstance(inner_profiles[0], list):
        # Check if it's a list of regions (each region is a list of segments)
        # A segment is also a list, so check deeper
        if inner_profiles and inner_profiles[0] and isline(inner_profiles[0][0]):
            # It's a list of region2d's
            inner_profiles_list = inner_profiles
        else:
            # Single region2d
            inner_profiles_list = [inner_profiles]
    else:
        # Single region2d
        inner_profiles_list = [inner_profiles]

    # Build profile wires in XZ plane
    outer_wire_template = _build_wire_from_region2d(profile)
    inner_wire_templates = [_build_wire_from_region2d(inner, reverse=True)
                           for inner in inner_profiles_list]

    # Get adaptive samples
    samples = _adaptive_sample_path3d(spine, angle_threshold_deg)
    if len(samples) < 2:
        raise ValueError("Path too short for adaptive sweep")

    # Initialize frame tracking
    if frame_mode == 'custom' and up_samples:
        initial_up = _interpolate_up_samples(up_samples, 0)
    else:
        initial_up = [0, 0, 1]  # Default Z-up

    # Check if path is closed (first and last sample at same position)
    first_pt = samples[0][0]
    last_pt = samples[-1][0]
    dist_sq = sum((a - b) ** 2 for a, b in zip(first_pt, last_pt))
    is_closed_path = dist_sq < 1e-6  # Within 1 micron

    # Helper to build a lofted solid from wires at sample locations
    def _build_loft_from_wire_template(wire_template, samples_list, ruled_mode, is_closed):
        """Build lofted solid from wire template at sample locations.

        For closed paths, reuses the first wire to close the loft instead of
        creating a duplicate wire at the same position (which causes mesh errors).
        """
        from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections

        loft = BRepOffsetAPI_ThruSections(True, ruled_mode, 1.0e-6)

        prev_up_local = initial_up
        first_wire = None

        # For closed paths, skip the last sample (we'll reuse first wire instead)
        num_samples = len(samples_list) - 1 if is_closed else len(samples_list)

        for i in range(num_samples):
            point, tangent, param = samples_list[i]

            # Compute frame based on mode
            if frame_mode == 'frenet' and i > 0:
                # Approximate derivative from previous/next tangent
                if i < len(samples_list) - 1:
                    _, next_tangent, _ = samples_list[i + 1]
                    tangent_deriv = _vector_subtract(next_tangent, tangent)
                else:
                    _, prev_tangent, _ = samples_list[i - 1]
                    tangent_deriv = _vector_subtract(tangent, prev_tangent)

                frame = _compute_frenet_frame(tangent, tangent_deriv)
                if frame is None:
                    right, up_local = _compute_minimal_twist_frame(tangent, prev_up_local)
                else:
                    right, up_local = frame
            elif frame_mode == 'custom' and up_samples:
                target_up = _interpolate_up_samples(up_samples, param)
                right, up_local = _compute_minimal_twist_frame(tangent, target_up)
            else:
                right, up_local = _compute_minimal_twist_frame(tangent, prev_up_local)

            prev_up_local = up_local

            # Transform wire to frame
            wire_transformed = _transform_wire_to_frame(wire_template, point, tangent, right, up_local)
            loft.AddWire(wire_transformed)

            # Save first wire for closing
            if i == 0:
                first_wire = wire_transformed

        # For closed paths, add the first wire again to properly close the loft
        if is_closed and first_wire is not None:
            loft.AddWire(first_wire)

        loft.Build()
        if not loft.IsDone():
            raise RuntimeError("Adaptive sweep loft failed")

        return loft.Shape()

    # Build outer loft
    outer_shape = _build_loft_from_wire_template(outer_wire_template, samples, ruled, is_closed_path)

    # If we have inner profiles, build inner lofts and subtract
    if inner_wire_templates:
        from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut

        shape = outer_shape
        for inner_template in inner_wire_templates:
            inner_shape = _build_loft_from_wire_template(inner_template, samples, ruled, is_closed_path)
            cut_op = BRepAlgoAPI_Cut(shape, inner_shape)
            cut_op.Build()
            if not cut_op.IsDone():
                raise RuntimeError("Boolean cut for hollow profile failed")
            shape = cut_op.Shape()
    else:
        shape = outer_shape

    # Normalize shape: if it's a Compound containing exactly one Solid, extract that Solid
    # This ensures proper behavior in subsequent boolean operations
    from OCC.Core.TopExp import TopExp_Explorer
    from OCC.Core.TopAbs import TopAbs_SOLID
    from OCC.Core.TopoDS import topods

    if shape.ShapeType() == 0:  # Compound
        exp = TopExp_Explorer(shape, TopAbs_SOLID)
        solids = []
        while exp.More():
            solids.append(topods.Solid(exp.Current()))
            exp.Next()
        if len(solids) == 1:
            shape = solids[0]

    # Compute volume for verification
    props = GProp_GProps()
    brepgprop.VolumeProperties(shape, props)
    volume = abs(props.Mass())

    # Create yapCAD solid with BREP
    construction = ['procedure', 'sweep_adaptive']
    sld = solid([], [], construction)

    # Attach BREP
    from yapcad.brep import attach_brep_to_solid, BrepSolid
    brep_solid = BrepSolid(shape)
    attach_brep_to_solid(sld, brep_solid)

    # Attach metadata if provided
    if metadata:
        from yapcad.metadata import get_solid_metadata
        meta = get_solid_metadata(sld, create=True)
        meta.update(metadata)

    return sld


def extrude_region2d(profile, height, *, direction=None, metadata=None):
    """
    Extrude a 2D region (list of line/arc segments) to create a solid.

    This uses OCC's BRepPrimAPI_MakePrism to create the extruded solid.
    The profile should be a closed 2D region (list of line/arc segments in XY plane).

    Args:
        profile: A yapCAD region2d (list of 2D curve segments forming a closed loop)
        height: Extrusion distance
        direction: Optional direction vector [dx, dy, dz]. Defaults to [0, 0, 1] (Z-up).
        metadata: Optional dict of metadata to attach

    Returns:
        A yapCAD solid representing the extruded shape
    """
    if not occ_available():
        raise RuntimeError("extrude_region2d requires pythonocc-core")

    from OCC.Core.gp import gp_Pnt, gp_Dir, gp_Ax2, gp_Circ, gp_Vec
    from OCC.Core.BRepBuilderAPI import (
        BRepBuilderAPI_MakeWire, BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace
    )
    from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakePrism
    from OCC.Core.GC import GC_MakeArcOfCircle
    from OCC.Core.GProp import GProp_GProps
    from OCC.Core.BRepGProp import brepgprop
    import math

    if direction is None:
        direction = [0, 0, 1]

    # Build wire from region2d in XY plane
    wire_builder = BRepBuilderAPI_MakeWire()

    for seg in profile:
        if isline(seg):
            p1 = seg[0]
            p2 = seg[1]
            # Keep in XY plane (Z=0)
            start = gp_Pnt(p1[0], p1[1] if len(p1) > 1 else 0, 0)
            end = gp_Pnt(p2[0], p2[1] if len(p2) > 1 else 0, 0)
            edge = BRepBuilderAPI_MakeEdge(start, end).Edge()
            wire_builder.Add(edge)
        elif isarc(seg):
            center = seg[0]
            params = seg[1]
            radius = params[0]
            start_ang = math.radians(params[1])
            end_ang = math.radians(params[2])
            cx, cy = center[0], center[1] if len(center) > 1 else 0
            start_pt = gp_Pnt(cx + radius * math.cos(start_ang), cy + radius * math.sin(start_ang), 0)
            end_pt = gp_Pnt(cx + radius * math.cos(end_ang), cy + radius * math.sin(end_ang), 0)
            arc_center = gp_Pnt(cx, cy, 0)
            circ = gp_Circ(gp_Ax2(arc_center, gp_Dir(0, 0, 1)), radius)
            arc_maker = GC_MakeArcOfCircle(circ, start_pt, end_pt, True)
            if arc_maker.IsDone():
                edge = BRepBuilderAPI_MakeEdge(arc_maker.Value()).Edge()
                wire_builder.Add(edge)

    if not wire_builder.IsDone():
        raise RuntimeError("Failed to build wire from region2d")

    outer_wire = wire_builder.Wire()

    # Create face from wire
    face_maker = BRepBuilderAPI_MakeFace(outer_wire, True)
    if not face_maker.IsDone():
        raise RuntimeError("Failed to build face from wire")

    profile_face = face_maker.Face()

    # Create extrusion vector
    extrude_vec = gp_Vec(direction[0] * height, direction[1] * height, direction[2] * height)

    # Create prism (extrusion)
    prism = BRepPrimAPI_MakePrism(profile_face, extrude_vec, True)
    prism.Build()

    if not prism.IsDone():
        raise RuntimeError("Failed to create prism")

    shape = prism.Shape()

    # Create BREP wrapper and tessellate to get yapCAD surface
    from yapcad.brep import BrepSolid, attach_brep_to_solid
    brep = BrepSolid(shape)
    surface = brep.tessellate()

    # Create solid from the tessellated surface
    sld = solid([surface], [], ['procedure', f'extrude_region2d(height={height})'])

    # Attach BREP data to the solid
    attach_brep_to_solid(sld, brep)

    if metadata:
        from yapcad.metadata import get_solid_metadata
        meta = get_solid_metadata(sld, create=True)
        meta.update(metadata)

    return sld


def sweep_profile_along_path(profile, spine, *, inner_profile=None, metadata=None):
    """
    Sweep a 2D profile along a 3D path (spine) to create a solid.

    This uses OCC's BRepOffsetAPI_MakePipe to create the swept solid.
    The profile should be a closed 2D region (list of line/arc segments in XZ plane).
    The spine is a path3d dict containing line and arc segments.

    Args:
        profile: A yapCAD region2d (list of 2D curve segments forming a closed loop)
                 This is the OUTER boundary of the profile.
        spine: A path3d dict with format:
            {
                'type': 'path3d',
                'segments': [
                    {'type': 'line', 'start': [x,y,z], 'end': [x,y,z]},
                    {'type': 'arc', 'center': [x,y,z], 'start': [x,y,z], 'end': [x,y,z], 'normal': [nx,ny,nz]},
                    ...
                ]
            }
        inner_profile: Optional yapCAD region2d for the inner boundary (hole).
                       If provided, creates a hollow profile (like a box girder).
        metadata: Optional dict of metadata to attach

    Returns:
        A yapCAD solid representing the swept shape
    """
    if not occ_available():
        raise RuntimeError("sweep_profile_along_path requires pythonocc-core")

    from OCC.Core.gp import gp_Pnt, gp_Dir, gp_Ax2, gp_Circ, gp_Vec
    from OCC.Core.BRepBuilderAPI import (
        BRepBuilderAPI_MakeWire, BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace
    )
    from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_MakePipe, BRepOffsetAPI_MakePipeShell
    from OCC.Core.GC import GC_MakeArcOfCircle
    from OCC.Core.GProp import GProp_GProps
    from OCC.Core.BRepGProp import brepgprop
    import math

    def _build_wire_from_region2d(region, reverse=False):
        """Build an OCC wire from a yapCAD region2d in the XZ plane.

        Args:
            region: List of line/arc segments forming a closed loop
            reverse: If True, reverse the winding order (for inner holes)
        """
        wire_builder = BRepBuilderAPI_MakeWire()
        segments = list(region)
        if reverse:
            # Reverse the segment order and swap start/end of each segment
            segments = segments[::-1]

        for seg in segments:
            if isline(seg):
                p1 = seg[0]
                p2 = seg[1]
                if reverse:
                    p1, p2 = p2, p1  # Swap endpoints
                # Convert to XZ plane (Y=0, 2D Y -> 3D Z)
                start = gp_Pnt(p1[0], 0, p1[1] if len(p1) > 1 else 0)
                end = gp_Pnt(p2[0], 0, p2[1] if len(p2) > 1 else 0)
                edge = BRepBuilderAPI_MakeEdge(start, end).Edge()
                wire_builder.Add(edge)
            elif isarc(seg):
                center = seg[0]
                params = seg[1]
                radius = params[0]
                start_ang = math.radians(params[1])
                end_ang = math.radians(params[2])
                if reverse:
                    start_ang, end_ang = end_ang, start_ang  # Swap angles
                cx, cy = center[0], center[1] if len(center) > 1 else 0
                start_pt = gp_Pnt(cx + radius * math.cos(start_ang), 0, cy + radius * math.sin(start_ang))
                end_pt = gp_Pnt(cx + radius * math.cos(end_ang), 0, cy + radius * math.sin(end_ang))
                arc_center = gp_Pnt(cx, 0, cy)
                circ = gp_Circ(gp_Ax2(arc_center, gp_Dir(0, 1, 0)), radius)
                arc_maker = GC_MakeArcOfCircle(circ, start_pt, end_pt, True)
                if arc_maker.IsDone():
                    edge = BRepBuilderAPI_MakeEdge(arc_maker.Value()).Edge()
                    wire_builder.Add(edge)
        return wire_builder.Wire()

    # Build outer profile wire
    outer_wire = _build_wire_from_region2d(profile)

    # Build inner profile wire if provided (for hollow shapes)
    inner_wire = None
    if inner_profile is not None:
        inner_wire = _build_wire_from_region2d(inner_profile, reverse=True)

    # Build spine wire from path3d
    spine_wire = BRepBuilderAPI_MakeWire()
    segments = spine.get('segments', [])

    for seg in segments:
        seg_type = seg.get('type', 'line')
        if seg_type == 'line':
            start = seg['start']
            end = seg['end']
            p1 = gp_Pnt(start[0], start[1], start[2])
            p2 = gp_Pnt(end[0], end[1], end[2])
            edge = BRepBuilderAPI_MakeEdge(p1, p2).Edge()
            spine_wire.Add(edge)
        elif seg_type == 'full_circle':
            # Full circle segment - creates a single closed circular edge
            # This avoids the topology issues that arise with multi-arc paths
            center = seg['center']
            radius = seg['radius']
            normal = seg.get('normal', [0, 0, 1])

            circ = gp_Circ(gp_Ax2(gp_Pnt(center[0], center[1], center[2]),
                                   gp_Dir(normal[0], normal[1], normal[2])), radius)
            edge = BRepBuilderAPI_MakeEdge(circ).Edge()
            spine_wire.Add(edge)
        elif seg_type == 'tilted_circle':
            # Tilted circle segment - circle in a tilted plane
            center = seg['center']
            radius = seg['radius']
            normal = seg.get('normal', [0, 0, 1])

            circ = gp_Circ(gp_Ax2(gp_Pnt(center[0], center[1], center[2]),
                                   gp_Dir(normal[0], normal[1], normal[2])), radius)
            edge = BRepBuilderAPI_MakeEdge(circ).Edge()
            spine_wire.Add(edge)
        elif seg_type == 'arc':
            center = seg['center']
            start = seg['start']
            end = seg['end']
            normal = seg.get('normal', [0, 0, 1])
            radius = seg.get('radius')

            if radius is None:
                # Calculate radius from center to start
                dx = start[0] - center[0]
                dy = start[1] - center[1]
                dz = start[2] - center[2]
                radius = math.sqrt(dx*dx + dy*dy + dz*dz)

            arc_center = gp_Pnt(center[0], center[1], center[2])
            arc_start = gp_Pnt(start[0], start[1], start[2])
            arc_end = gp_Pnt(end[0], end[1], end[2])
            arc_normal = gp_Dir(normal[0], normal[1], normal[2])

            circ = gp_Circ(gp_Ax2(arc_center, arc_normal), radius)
            arc_maker = GC_MakeArcOfCircle(circ, arc_start, arc_end, True)
            if arc_maker.IsDone():
                edge = BRepBuilderAPI_MakeEdge(arc_maker.Value()).Edge()
                spine_wire.Add(edge)

    spine = spine_wire.Wire()

    # Get the spine start point and tangent to properly position and orient the profile
    # MakePipeShell needs the profile at the spine start, perpendicular to tangent
    from OCC.Core.BRepAdaptor import BRepAdaptor_CompCurve
    from OCC.Core.gp import gp_Trsf, gp_XYZ, gp_Ax1, gp_Ax3, gp_Pnt, gp_Dir
    from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_Transform

    spine_curve = BRepAdaptor_CompCurve(spine, True)
    param_start = spine_curve.FirstParameter()
    spine_start_pt = spine_curve.Value(param_start)

    # Get tangent direction at start
    tangent_vec = gp_Vec()
    spine_curve.D1(param_start, spine_start_pt, tangent_vec)
    tangent_vec.Normalize()
    tangent_dir = gp_Dir(tangent_vec)

    # Profile is built in XZ plane with normal pointing in Y direction
    # We need to transform it so the normal aligns with the spine tangent at start
    profile_normal = gp_Dir(0, 1, 0)  # Profile faces +Y

    # Check if transformation is needed
    dot = tangent_dir.X() * profile_normal.X() + \
          tangent_dir.Y() * profile_normal.Y() + \
          tangent_dir.Z() * profile_normal.Z()

    need_rotation = abs(dot - 1.0) > 1e-6  # Not already aligned
    need_translation = (abs(spine_start_pt.X()) > 1e-6 or
                        abs(spine_start_pt.Y()) > 1e-6 or
                        abs(spine_start_pt.Z()) > 1e-6)

    if need_rotation or need_translation:
        # Apply rotation first (if needed), then translation
        # Profile is in XZ plane with normal = +Y
        # We need normal to align with tangent

        trsf = gp_Trsf()

        if need_rotation:
            # Find rotation axis = profile_normal cross tangent
            # And rotation angle = angle between them
            cross = gp_Vec(profile_normal).Crossed(gp_Vec(tangent_dir))
            cross_mag = cross.Magnitude()

            if cross_mag > 1e-6:  # Not parallel/antiparallel
                cross.Normalize()
                rot_axis = gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(cross))
                # Angle from dot product (already computed above)
                angle = math.acos(max(-1.0, min(1.0, dot)))
                trsf.SetRotation(rot_axis, angle)

        # Apply translation to move profile to spine start
        trans = gp_Trsf()
        trans.SetTranslation(gp_Vec(spine_start_pt.X(), spine_start_pt.Y(), spine_start_pt.Z()))

        # Combine transformations: first rotate, then translate
        combined = trans.Multiplied(trsf)

        outer_transformer = BRepBuilderAPI_Transform(outer_wire, combined, True)
        outer_wire = outer_transformer.Shape()

        if inner_wire is not None:
            inner_transformer = BRepBuilderAPI_Transform(inner_wire, combined, True)
            inner_wire = inner_transformer.Shape()

    # Create pipe shell (sweep)
    # Using MakePipeShell instead of MakePipe for better handling of
    # closed paths (like circular rings) which otherwise produce invalid geometry.
    #
    # MakePipeShell provides better control over profile orientation and
    # produces valid manifold geometry for closed spine paths.
    pipe_shell = BRepOffsetAPI_MakePipeShell(spine)

    # Add the outer profile wire
    pipe_shell.Add(outer_wire)

    # Add inner profile wire if present (for hollow shapes)
    if inner_wire is not None:
        pipe_shell.Add(inner_wire)

    # Set mode to keep profile orientation consistent with Z axis
    # This helps avoid profile flipping/twisting along the path
    pipe_shell.SetMode(gp_Dir(0, 0, 1))

    pipe_shell.Build()

    if not pipe_shell.IsDone():
        raise RuntimeError("OCC pipe shell sweep failed")

    # Make it a solid (not just a shell)
    pipe_shell.MakeSolid()
    shape = pipe_shell.Shape()

    # Get volume for verification
    props = GProp_GProps()
    brepgprop.VolumeProperties(shape, props)
    volume = abs(props.Mass())

    # Create yapCAD solid with BREP attachment
    # For now, create a minimal solid structure and attach BREP
    from yapcad.geom3d import solid
    call = f"yapcad.geom3d_util.sweep_profile_along_path(profile, spine)"
    construction = ['procedure', call]

    if metadata is not None and isinstance(metadata, dict):
        sld = solid([], [], construction, metadata)
    else:
        sld = solid([], [], construction)

    # Attach BREP
    try:
        attach_brep_to_solid(sld, BrepSolid(shape))
    except Exception:
        pass

    return sld


def make_occ_helix(radius, pitch, height, left_hand=False):
    """
    Create a true helix wire using 2D parametric curve on cylindrical surface.

    This produces a mathematically exact helix, not a polyline approximation.
    Uses the standard OCC technique: 2D line on Geom_CylindricalSurface.

    :param radius: Radius of the helix cylinder
    :param pitch: Vertical rise per full turn (360 degrees)
    :param height: Total height of the helix
    :param left_hand: If True, create left-handed helix (default False = right-handed)
    :returns: An OCC TopoDS_Wire representing the helix

    .. note::
       This function requires pythonocc-core and will raise RuntimeError if not available.

    Example::

        >>> # Create a right-handed helix
        >>> helix = make_occ_helix(radius=10.0, pitch=5.0, height=20.0)
        >>> # Create a left-handed helix
        >>> left_helix = make_occ_helix(radius=10.0, pitch=5.0, height=20.0, left_hand=True)
    """
    if not occ_available():
        raise RuntimeError("make_occ_helix requires pythonocc-core")

    from OCC.Core.gp import gp_Ax3, gp_Pnt, gp_Dir, gp_Pnt2d, gp_Dir2d, gp_Lin2d
    from OCC.Core.Geom import Geom_CylindricalSurface
    from OCC.Core.GCE2d import GCE2d_MakeSegment
    from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire
    from OCC.Core.BRepLib import BRepLib

    # Create cylindrical surface centered at origin with axis along Z
    # gp_Ax3 defines a coordinate system: origin, Z direction, X direction
    # For XOY plane: origin at (0,0,0), Z along (0,0,1), X along (1,0,0)
    ax3_xoy = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0))
    cylinder = Geom_CylindricalSurface(ax3_xoy, radius)

    # Calculate the number of turns based on height and pitch
    turns = height / pitch

    # Direction in UV space: U is angle (in radians), V is height
    # For a right-handed helix going up, U increases as V increases
    # For a left-handed helix, U decreases as V increases
    sign = -1.0 if left_hand else 1.0

    # In UV space of cylinder: U = angle (radians), V = height along axis
    # A helix is a straight line in UV space: for each turn (2*pi in U), we rise by pitch in V
    # So direction is (sign * 2*pi, pitch) and we travel for 'turns' turns
    direction = gp_Dir2d(sign * 2.0 * math.pi, pitch)

    # Create 2D line from origin in the UV space
    line_2d = gp_Lin2d(gp_Pnt2d(0.0, 0.0), direction)

    # Parameter range: from 0 to 'turns' (each unit = one full turn)
    segment = GCE2d_MakeSegment(line_2d, 0.0, turns).Value()

    # Create 3D edge from 2D curve on cylindrical surface
    helix_edge = BRepBuilderAPI_MakeEdge(segment, cylinder, 0.0, turns).Edge()

    # Build the 3D curve representation
    BRepLib.BuildCurves3d_s(helix_edge)

    # Create wire from edge
    return BRepBuilderAPI_MakeWire(helix_edge).Wire()


def helical_extrude(profile, height, twist_angle_deg, *,
                    auxiliary_radius=10.0, segments=64, metadata=None):
    """
    Extrude a 2D profile along Z with smooth helical twist.

    This function creates a helical extrusion where the profile rotates around
    the Z axis as it extrudes upward. The implementation uses high-resolution
    lofting through many intermediate sections to produce smooth helical surfaces.

    For profiles centered at the origin, the twist rotates the entire profile
    around Z. This is ideal for helical gears, twisted columns, and spiral features.

    :param profile: A yapCAD region2d (list of 2D curve segments forming a closed loop)
                    in the XY plane, centered at or near the origin
    :param height: Total extrusion height along Z axis
    :param twist_angle_deg: Total twist angle in degrees over the full height.
                           Positive = counterclockwise when viewed from +Z.
                           Zero twist will produce a simple extrusion.
    :param auxiliary_radius: (Deprecated, kept for API compatibility) Previously used
                            for auxiliary helix. Now ignored.
    :param segments: Number of intermediate sections for lofting (default 64).
                     More segments = smoother helical surface. For helical gears,
                     use at least 64 segments to avoid visible stepping.
    :param metadata: Optional dict of metadata to attach
    :returns: A yapCAD solid with smooth helical surfaces

    .. note::
       This function requires pythonocc-core and will raise RuntimeError if not available.
       For small twist angles (<10 degrees), the result may be similar to a simple extrusion.
       The helical effect is most visible with larger twist angles. Using segments=64 or
       higher is recommended for smooth surfaces.

    Example::

        >>> # Create a twisted square prism with smooth surfaces
        >>> from yapcad.geom import line, point
        >>> square = [
        ...     line(point(-5, -5), point(5, -5)),
        ...     line(point(5, -5), point(5, 5)),
        ...     line(point(5, 5), point(-5, 5)),
        ...     line(point(-5, 5), point(-5, -5))
        ... ]
        >>> twisted = helical_extrude(square, height=20, twist_angle_deg=90)
    """
    if not occ_available():
        raise RuntimeError("helical_extrude requires pythonocc-core")

    from OCC.Core.gp import gp_Pnt, gp_Dir, gp_Ax2, gp_Circ, gp_Trsf, gp_Ax1
    from OCC.Core.BRepBuilderAPI import (
        BRepBuilderAPI_MakeWire, BRepBuilderAPI_MakeEdge,
        BRepBuilderAPI_MakeFace, BRepBuilderAPI_Transform
    )
    from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections
    from OCC.Core.GC import GC_MakeArcOfCircle
    from OCC.Core.GProp import GProp_GProps
    from OCC.Core.BRepGProp import brepgprop
    from OCC.Core.TopoDS import topods

    # Handle zero or very small twist as simple extrusion
    if abs(twist_angle_deg) < 1e-6:
        return extrude_region2d(profile, height, metadata=metadata)

    def _build_wire_from_region2d_xy(region):
        """Build an OCC wire from a yapCAD region2d in the XY plane (Z=0)."""
        wire_builder = BRepBuilderAPI_MakeWire()

        for seg in region:
            if isline(seg):
                p1 = seg[0]
                p2 = seg[1]
                # Keep in XY plane (Z=0)
                start = gp_Pnt(p1[0], p1[1] if len(p1) > 1 else 0, 0)
                end = gp_Pnt(p2[0], p2[1] if len(p2) > 1 else 0, 0)
                edge = BRepBuilderAPI_MakeEdge(start, end).Edge()
                wire_builder.Add(edge)
            elif isarc(seg):
                center = seg[0]
                params = seg[1]
                radius = params[0]
                start_ang = math.radians(params[1])
                end_ang = math.radians(params[2])
                cx, cy = center[0], center[1] if len(center) > 1 else 0
                start_pt = gp_Pnt(cx + radius * math.cos(start_ang),
                                  cy + radius * math.sin(start_ang), 0)
                end_pt = gp_Pnt(cx + radius * math.cos(end_ang),
                                cy + radius * math.sin(end_ang), 0)
                arc_center = gp_Pnt(cx, cy, 0)
                circ = gp_Circ(gp_Ax2(arc_center, gp_Dir(0, 0, 1)), radius)
                arc_maker = GC_MakeArcOfCircle(circ, start_pt, end_pt, True)
                if arc_maker.IsDone():
                    edge = BRepBuilderAPI_MakeEdge(arc_maker.Value()).Edge()
                    wire_builder.Add(edge)

        if not wire_builder.IsDone():
            raise RuntimeError("Failed to build wire from region2d")
        return wire_builder.Wire()

    def _transform_wire(wire, z_offset, rotation_deg):
        """Transform wire: translate along Z and rotate around Z axis."""
        trsf = gp_Trsf()

        # Combined transformation: first rotate around Z, then translate along Z
        # Create rotation around Z axis at origin
        rotation_rad = math.radians(rotation_deg)
        z_axis = gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1))
        trsf.SetRotation(z_axis, rotation_rad)

        # Apply rotation
        transformer = BRepBuilderAPI_Transform(wire, trsf, True)
        rotated_wire = topods.Wire(transformer.Shape())

        # Apply translation
        trsf_translate = gp_Trsf()
        trsf_translate.SetTranslation(gp_Pnt(0, 0, 0), gp_Pnt(0, 0, z_offset))
        transformer2 = BRepBuilderAPI_Transform(rotated_wire, trsf_translate, True)

        return topods.Wire(transformer2.Shape())

    # Build template wire at Z=0
    template_wire = _build_wire_from_region2d_xy(profile)

    # Create loft through rotated sections
    # Use ruled=True for better handling of profiles with straight edges
    loft = BRepOffsetAPI_ThruSections(True, True, 1.0e-6)

    # Add sections from bottom to top
    for i in range(segments + 1):
        t = i / segments  # 0 to 1
        z = t * height
        angle = t * twist_angle_deg

        transformed_wire = _transform_wire(template_wire, z, angle)
        loft.AddWire(transformed_wire)

    loft.Build()
    if not loft.IsDone():
        raise RuntimeError("Failed to create helical extrusion loft")

    shape = loft.Shape()

    # Normalize shape: if it's a Compound containing exactly one Solid, extract it
    from OCC.Core.TopExp import TopExp_Explorer
    from OCC.Core.TopAbs import TopAbs_SOLID

    if shape.ShapeType() == 0:  # Compound
        exp = TopExp_Explorer(shape, TopAbs_SOLID)
        solids = []
        while exp.More():
            solids.append(topods.Solid(exp.Current()))
            exp.Next()
        if len(solids) == 1:
            shape = solids[0]

    # Compute volume for verification
    props = GProp_GProps()
    brepgprop.VolumeProperties(shape, props)
    volume = abs(props.Mass())

    # Create BREP wrapper and tessellate to get yapCAD surface
    from yapcad.brep import BrepSolid, attach_brep_to_solid
    brep = BrepSolid(shape)
    surface = brep.tessellate()

    # Create solid from the tessellated surface
    construction = ['procedure',
                    f'helical_extrude(height={height}, twist={twist_angle_deg}deg)']
    sld = solid([surface], [], construction)

    # Attach BREP data to the solid
    attach_brep_to_solid(sld, brep)

    if metadata:
        from yapcad.metadata import get_solid_metadata
        meta = get_solid_metadata(sld, create=True)
        meta.update(metadata)

    return sld


def radial_pattern_solid(solid_geom, count, center=None, axis=None, angle=None):
    """
    Create a radial/circular pattern of solid copies.

    Creates multiple copies of a 3D solid arranged in a circular pattern
    around a center point, rotating about a specified axis.

    :param solid_geom: A yapCAD solid to pattern
    :param count: Number of copies (including original)
    :param center: Center point for rotation (default [0,0,0,1])
    :param axis: Rotation axis vector (default [0,0,1,0] for Z-axis)
    :param angle: Total angle to span in degrees (default 360 = full circle)
    :returns: List of solid copies, each rotated by angle/count increments

    Example::

        >>> # Create 6 holes around a circle
        >>> hole = conic(2.5, 2.5, 10)  # cylinder hole
        >>> holes = radial_pattern_solid(hole, count=6)  # 6 holes at 60 degree intervals
    """
    if count < 1:
        raise ValueError("count must be at least 1")

    if not issolid(solid_geom):
        raise ValueError("solid_geom must be a valid solid")

    # Set defaults
    if center is None:
        center = point(0, 0, 0)
    if axis is None:
        axis = point(0, 0, 1)
    if angle is None:
        angle = 360.0

    # Single item returns the original
    if count == 1:
        return [solid_geom]

    # Calculate angle increment between copies
    # For a full 360 pattern, don't duplicate at start/end
    if close(angle, 360.0):
        angle_step = angle / count
    else:
        angle_step = angle / (count - 1) if count > 1 else 0

    result = []
    for i in range(count):
        current_angle = i * angle_step
        if close(current_angle, 0.0):
            # No rotation needed for the first copy
            result.append(solid_geom)
        else:
            rotated = rotatesolid(solid_geom, current_angle, cent=center, axis=axis)
            result.append(rotated)

    return result


def linear_pattern_solid(solid_geom, count, spacing):
    """
    Create a linear pattern of solid copies.

    Creates multiple copies of a 3D solid arranged in a line with uniform spacing.

    :param solid_geom: A yapCAD solid to pattern
    :param count: Number of copies (including original)
    :param spacing: Vector defining direction and distance between copies.
                    Can be a 4-tuple [x, y, z, w] or 3-tuple/list [x, y, z]
    :returns: List of solid copies, each translated by spacing increments

    Example::

        >>> # Create 5 mounting holes in a row
        >>> hole = conic(3, 3, 10)  # cylinder hole
        >>> holes = linear_pattern_solid(hole, count=5, spacing=[20, 0, 0])  # 5 holes, 20mm apart
    """
    if count < 1:
        raise ValueError("count must be at least 1")

    if not issolid(solid_geom):
        raise ValueError("solid_geom must be a valid solid")

    # Normalize spacing to a proper vector
    if isinstance(spacing, (list, tuple)):
        if len(spacing) == 3:
            spacing = point(spacing[0], spacing[1], spacing[2])
        elif len(spacing) == 4:
            spacing = point(spacing[0], spacing[1], spacing[2])
        elif len(spacing) == 2:
            spacing = point(spacing[0], spacing[1], 0)
        else:
            raise ValueError("spacing must be a 2D, 3D, or 4D vector")
    else:
        raise ValueError("spacing must be a list or tuple")

    result = []
    for i in range(count):
        if i == 0:
            # No translation for the first copy
            result.append(solid_geom)
        else:
            # Calculate total offset for this copy
            delta = scale3(spacing, float(i))
            translated = translatesolid(solid_geom, delta)
            result.append(translated)

    return result


def radial_pattern_surface(surf, count, center=None, axis=None, angle=None):
    """
    Create a radial/circular pattern of surface copies.

    Creates multiple copies of a 3D surface arranged in a circular pattern
    around a center point, rotating about a specified axis.

    :param surf: A yapCAD surface to pattern
    :param count: Number of copies (including original)
    :param center: Center point for rotation (default [0,0,0,1])
    :param axis: Rotation axis vector (default [0,0,1,0] for Z-axis)
    :param angle: Total angle to span in degrees (default 360 = full circle)
    :returns: List of surface copies, each rotated by angle/count increments

    Example::

        >>> # Create 8 fins around a rocket body
        >>> fin_surface = triangulate_region(fin_profile)
        >>> fins = radial_pattern_surface(fin_surface, count=8)
    """
    if count < 1:
        raise ValueError("count must be at least 1")

    if not issurface(surf):
        raise ValueError("surf must be a valid surface")

    # Set defaults
    if center is None:
        center = point(0, 0, 0)
    if axis is None:
        axis = point(0, 0, 1)
    if angle is None:
        angle = 360.0

    # Single item returns the original
    if count == 1:
        return [surf]

    # Calculate angle increment between copies
    if close(angle, 360.0):
        angle_step = angle / count
    else:
        angle_step = angle / (count - 1) if count > 1 else 0

    result = []
    for i in range(count):
        current_angle = i * angle_step
        if close(current_angle, 0.0):
            result.append(surf)
        else:
            rotated = rotatesurface(surf, current_angle, cent=center, axis=axis)
            result.append(rotated)

    return result


def linear_pattern_surface(surf, count, spacing):
    """
    Create a linear pattern of surface copies.

    Creates multiple copies of a 3D surface arranged in a line with uniform spacing.

    :param surf: A yapCAD surface to pattern
    :param count: Number of copies (including original)
    :param spacing: Vector defining direction and distance between copies.
                    Can be a 2D, 3D, or 4D vector
    :returns: List of surface copies, each translated by spacing increments

    Example::

        >>> # Create a row of ribs along a structure
        >>> rib_surface = triangulate_region(rib_profile)
        >>> ribs = linear_pattern_surface(rib_surface, count=10, spacing=[5, 0, 0])
    """
    if count < 1:
        raise ValueError("count must be at least 1")

    if not issurface(surf):
        raise ValueError("surf must be a valid surface")

    # Normalize spacing to a proper vector
    if isinstance(spacing, (list, tuple)):
        if len(spacing) == 3:
            spacing = point(spacing[0], spacing[1], spacing[2])
        elif len(spacing) == 4:
            spacing = point(spacing[0], spacing[1], spacing[2])
        elif len(spacing) == 2:
            spacing = point(spacing[0], spacing[1], 0)
        else:
            raise ValueError("spacing must be a 2D, 3D, or 4D vector")
    else:
        raise ValueError("spacing must be a list or tuple")

    result = []
    for i in range(count):
        if i == 0:
            result.append(surf)
        else:
            delta = scale3(spacing, float(i))
            translated = translatesurface(surf, delta)
            result.append(translated)

    return result
