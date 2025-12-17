"""OCC ↔ Native BREP conversion utilities.

This module provides bidirectional conversion between OCC (Open CASCADE) BREP
representations and yapCAD native BREP representations.

The native BREP representation is the primary/canonical form for yapCAD geometry.
OCC representations are derived and used for:
- STEP file import/export
- Boolean operations (when OCC engine is selected)
- Complex surface operations

Conversion functions:
- occ_surface_to_native: Convert OCC Geom_Surface to native analytic surface
- occ_solid_to_native_brep: Convert OCC TopoDS_Solid to native BREP topology
- native_surface_to_occ: Convert native analytic surface to OCC Geom_Surface
- native_brep_to_occ: Convert native BREP topology to OCC TopoDS_Solid

Copyright (c) 2025 Richard DeVaul
MIT License
"""

from math import sqrt, pi, atan2, acos
from typing import Optional, Any, Tuple, List

from yapcad.geom import point

# OCC imports - gracefully handle missing dependency
try:
    from OCC.Core.TopoDS import (
        TopoDS_Shape, TopoDS_Solid, TopoDS_Shell, TopoDS_Face,
        TopoDS_Wire, TopoDS_Edge, TopoDS_Vertex, topods
    )
    from OCC.Core.TopExp import TopExp_Explorer, topexp
    from OCC.Core.TopAbs import (
        TopAbs_SOLID, TopAbs_SHELL, TopAbs_FACE, TopAbs_WIRE,
        TopAbs_EDGE, TopAbs_VERTEX, TopAbs_FORWARD, TopAbs_REVERSED
    )
    from OCC.Core.BRep import BRep_Tool
    from OCC.Core.BRepAdaptor import BRepAdaptor_Surface, BRepAdaptor_Curve
    from OCC.Core.GeomAbs import (
        GeomAbs_Plane, GeomAbs_Cylinder, GeomAbs_Cone,
        GeomAbs_Sphere, GeomAbs_Torus, GeomAbs_BSplineSurface,
        GeomAbs_BezierSurface, GeomAbs_SurfaceOfRevolution,
        GeomAbs_SurfaceOfExtrusion, GeomAbs_OffsetSurface,
        GeomAbs_OtherSurface,
        GeomAbs_Line, GeomAbs_Circle, GeomAbs_Ellipse,
        GeomAbs_BSplineCurve, GeomAbs_BezierCurve, GeomAbs_OtherCurve
    )
    from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir, gp_Ax1, gp_Ax2, gp_Ax3
    from OCC.Core.TopLoc import TopLoc_Location
    from OCC.Core.BRepMesh import BRepMesh_IncrementalMesh
    from OCC.Core.Geom import (
        Geom_Plane, Geom_CylindricalSurface, Geom_ConicalSurface,
        Geom_SphericalSurface, Geom_ToroidalSurface, Geom_BSplineSurface
    )
    from OCC.Core.BRepBuilderAPI import (
        BRepBuilderAPI_MakeVertex, BRepBuilderAPI_MakeEdge,
        BRepBuilderAPI_MakeWire, BRepBuilderAPI_MakeFace,
        BRepBuilderAPI_MakeShell, BRepBuilderAPI_MakeSolid,
        BRepBuilderAPI_Sewing
    )
    from OCC.Core.TColgp import TColgp_Array2OfPnt, TColgp_Array1OfPnt
    from OCC.Core.TColStd import TColStd_Array1OfReal, TColStd_Array1OfInteger
    from OCC.Core.ShapeFix import ShapeFix_Solid, ShapeFix_Shell
    from OCC.Core.Geom import Geom_Line, Geom_Circle, Geom_BSplineCurve
    from OCC.Core.GC import GC_MakeCircle, GC_MakeLine, GC_MakeArcOfCircle
    _OCC_AVAILABLE = True
except ImportError:
    _OCC_AVAILABLE = False
    # Type stubs for when OCC is not available
    TopoDS_Shape = TopoDS_Solid = TopoDS_Shell = TopoDS_Face = Any
    TopoDS_Wire = TopoDS_Edge = TopoDS_Vertex = Any
    BRepAdaptor_Surface = BRepAdaptor_Curve = Any

from yapcad.analytic_surfaces import (
    plane_surface, sphere_surface, cylinder_surface, cone_surface, torus_surface,
    bspline_surface, tessellated_surface, is_analytic_surface
)
from yapcad.native_brep import (
    brep_vertex, brep_edge, brep_trim, brep_loop, brep_face, brep_shell, brep_solid,
    line_edge, circle_edge, bspline_edge,
    vertex_location, edge_vertices,
    TopologyGraph
)


def occ_available() -> bool:
    """Return True if OCC is available."""
    return _OCC_AVAILABLE


def require_occ() -> None:
    """Raise error if OCC is not available."""
    if not _OCC_AVAILABLE:
        raise RuntimeError(
            "pythonocc-core is not available. Activate the yapcad-brep conda "
            "environment to use OCC conversion features."
        )


# -----------------------------------------------------------------------------
# OCC Surface → Native Surface Conversion
# -----------------------------------------------------------------------------

def _gp_pnt_to_point(pnt) -> list:
    """Convert gp_Pnt to yapCAD point."""
    return point(pnt.X(), pnt.Y(), pnt.Z())


def _gp_dir_to_vector(d) -> list:
    """Convert gp_Dir to yapCAD vector."""
    return [d.X(), d.Y(), d.Z(), 0.0]


def occ_surface_to_native(face: TopoDS_Face, tessellate_fallback: bool = True):
    """Convert an OCC face's surface to a native analytic surface.

    Parameters
    ----------
    face : TopoDS_Face
        The OCC face to convert.
    tessellate_fallback : bool
        If True, fall back to tessellation for unsupported surface types.

    Returns
    -------
    native_surface
        A native analytic surface representation, or None if conversion fails.
    """
    require_occ()

    adaptor = BRepAdaptor_Surface(face)
    surface_type = adaptor.GetType()

    # Get parameter bounds
    u_min, u_max, v_min, v_max = adaptor.FirstUParameter(), adaptor.LastUParameter(), \
                                  adaptor.FirstVParameter(), adaptor.LastVParameter()

    if surface_type == GeomAbs_Plane:
        return _convert_plane(adaptor, u_min, u_max, v_min, v_max)
    elif surface_type == GeomAbs_Sphere:
        return _convert_sphere(adaptor, u_min, u_max, v_min, v_max)
    elif surface_type == GeomAbs_Cylinder:
        return _convert_cylinder(adaptor, u_min, u_max, v_min, v_max)
    elif surface_type == GeomAbs_Cone:
        return _convert_cone(adaptor, u_min, u_max, v_min, v_max)
    elif surface_type == GeomAbs_Torus:
        return _convert_torus(adaptor, u_min, u_max, v_min, v_max)
    elif surface_type == GeomAbs_BSplineSurface:
        return _convert_bspline_surface(adaptor, u_min, u_max, v_min, v_max)
    elif surface_type == GeomAbs_BezierSurface:
        # Convert Bezier to BSpline representation
        return _convert_bezier_surface(adaptor, u_min, u_max, v_min, v_max)
    else:
        # Fallback: tessellate the face
        if tessellate_fallback:
            return _tessellate_occ_face(face, u_min, u_max, v_min, v_max)
        return None


def _convert_plane(adaptor, u_min, u_max, v_min, v_max):
    """Convert OCC plane to native plane surface."""
    pln = adaptor.Plane()
    loc = pln.Location()
    axis = pln.Axis()

    origin = _gp_pnt_to_point(loc)
    normal = _gp_dir_to_vector(axis.Direction())

    return plane_surface(origin, normal, u_range=(u_min, u_max), v_range=(v_min, v_max))


def _convert_sphere(adaptor, u_min, u_max, v_min, v_max):
    """Convert OCC sphere to native sphere surface."""
    sph = adaptor.Sphere()
    center = _gp_pnt_to_point(sph.Location())
    radius = sph.Radius()

    return sphere_surface(center, radius, u_range=(u_min, u_max), v_range=(v_min, v_max))


def _convert_cylinder(adaptor, u_min, u_max, v_min, v_max):
    """Convert OCC cylinder to native cylinder surface."""
    cyl = adaptor.Cylinder()
    axis = cyl.Axis()
    loc = axis.Location()
    direction = axis.Direction()
    radius = cyl.Radius()

    axis_point = _gp_pnt_to_point(loc)
    axis_dir = _gp_dir_to_vector(direction)

    return cylinder_surface(axis_point, axis_dir, radius,
                           u_range=(u_min, u_max), v_range=(v_min, v_max))


def _convert_cone(adaptor, u_min, u_max, v_min, v_max):
    """Convert OCC cone to native cone surface."""
    cone = adaptor.Cone()
    apex = cone.Apex()
    axis = cone.Axis()
    half_angle = cone.SemiAngle()

    apex_pt = _gp_pnt_to_point(apex)
    axis_dir = _gp_dir_to_vector(axis.Direction())

    # OCC uses negative half_angle for cones pointing downward - use absolute value
    half_angle = abs(half_angle)

    return cone_surface(apex_pt, axis_dir, half_angle,
                       u_range=(u_min, u_max), v_range=(v_min, v_max))


def _convert_torus(adaptor, u_min, u_max, v_min, v_max):
    """Convert OCC torus to native torus surface."""
    torus = adaptor.Torus()
    axis = torus.Axis()
    loc = axis.Location()
    direction = axis.Direction()
    major_radius = torus.MajorRadius()
    minor_radius = torus.MinorRadius()

    center = _gp_pnt_to_point(loc)
    axis_dir = _gp_dir_to_vector(direction)

    return torus_surface(center, axis_dir, major_radius, minor_radius,
                        u_range=(u_min, u_max), v_range=(v_min, v_max))


def _convert_bspline_surface(adaptor, u_min, u_max, v_min, v_max):
    """Convert OCC B-spline surface to native B-spline surface."""
    bspl = adaptor.BSpline()

    # Get degrees
    u_degree = bspl.UDegree()
    v_degree = bspl.VDegree()

    # Get pole counts
    n_u = bspl.NbUPoles()
    n_v = bspl.NbVPoles()

    # Get control points
    control_points = []
    for j in range(1, n_v + 1):
        row = []
        for i in range(1, n_u + 1):
            pnt = bspl.Pole(i, j)
            row.append([pnt.X(), pnt.Y(), pnt.Z()])
        control_points.append(row)

    # Get weights (if rational)
    weights = None
    if bspl.IsURational() or bspl.IsVRational():
        weights = []
        for j in range(1, n_v + 1):
            row = []
            for i in range(1, n_u + 1):
                row.append(bspl.Weight(i, j))
            weights.append(row)

    # Get knot vectors
    u_knots = []
    for i in range(1, bspl.NbUKnots() + 1):
        knot = bspl.UKnot(i)
        mult = bspl.UMultiplicity(i)
        u_knots.extend([knot] * mult)

    v_knots = []
    for i in range(1, bspl.NbVKnots() + 1):
        knot = bspl.VKnot(i)
        mult = bspl.VMultiplicity(i)
        v_knots.extend([knot] * mult)

    return bspline_surface(control_points, u_knots, v_knots, u_degree, v_degree,
                          weights=weights, u_range=(u_min, u_max), v_range=(v_min, v_max))


def _convert_bezier_surface(adaptor, u_min, u_max, v_min, v_max):
    """Convert OCC Bezier surface to native B-spline surface.

    Bezier surfaces are a special case of B-splines with specific knot vectors.
    """
    bez = adaptor.Bezier()

    # Get degrees
    u_degree = bez.UDegree()
    v_degree = bez.VDegree()

    # Get pole counts
    n_u = bez.NbUPoles()
    n_v = bez.NbVPoles()

    # Get control points
    control_points = []
    for j in range(1, n_v + 1):
        row = []
        for i in range(1, n_u + 1):
            pnt = bez.Pole(i, j)
            row.append([pnt.X(), pnt.Y(), pnt.Z()])
        control_points.append(row)

    # Get weights (if rational)
    weights = None
    if bez.IsURational() or bez.IsVRational():
        weights = []
        for j in range(1, n_v + 1):
            row = []
            for i in range(1, n_u + 1):
                row.append(bez.Weight(i, j))
            weights.append(row)

    # Bezier knot vectors: [0, 0, ..., 0, 1, 1, ..., 1]
    # with (degree + 1) zeros and (degree + 1) ones
    u_knots = [0.0] * (u_degree + 1) + [1.0] * (u_degree + 1)
    v_knots = [0.0] * (v_degree + 1) + [1.0] * (v_degree + 1)

    return bspline_surface(control_points, u_knots, v_knots, u_degree, v_degree,
                          weights=weights, u_range=(u_min, u_max), v_range=(v_min, v_max))


def _tessellate_occ_face(face: TopoDS_Face, u_min, u_max, v_min, v_max,
                         deflection: float = 0.1):
    """Tessellate an OCC face and return as tessellated_surface.

    Used as fallback for unsupported surface types.
    """
    require_occ()

    # Ensure mesh is generated
    BRepMesh_IncrementalMesh(face, deflection)

    loc = TopLoc_Location()
    triangulation = BRep_Tool.Triangulation(face, loc)

    if triangulation is None:
        return None

    trsf = loc.Transformation()
    orientation = face.Orientation()
    reverse = (orientation == TopAbs_REVERSED)

    vertices = []
    normals = []

    # Extract vertices
    for i in range(1, triangulation.NbNodes() + 1):
        pnt = triangulation.Node(i).Transformed(trsf)
        vertices.append([pnt.X(), pnt.Y(), pnt.Z()])

    # Extract normals if available
    if triangulation.HasNormals():
        for i in range(1, triangulation.NbNodes() + 1):
            n = triangulation.Normal(i)
            vec = gp_Vec(n.X(), n.Y(), n.Z())
            vec.Transform(trsf)
            if reverse:
                normals.append([-vec.X(), -vec.Y(), -vec.Z()])
            else:
                normals.append([vec.X(), vec.Y(), vec.Z()])
    else:
        # Compute normals from triangles
        accum = [[0.0, 0.0, 0.0] for _ in range(triangulation.NbNodes())]
        for i in range(1, triangulation.NbTriangles() + 1):
            t = triangulation.Triangle(i)
            n1, n2, n3 = t.Get()
            if reverse:
                n2, n3 = n3, n2
            v0 = vertices[n1 - 1]
            v1 = vertices[n2 - 1]
            v2 = vertices[n3 - 1]
            # Cross product
            nx = (v1[1] - v0[1]) * (v2[2] - v0[2]) - (v1[2] - v0[2]) * (v2[1] - v0[1])
            ny = (v1[2] - v0[2]) * (v2[0] - v0[0]) - (v1[0] - v0[0]) * (v2[2] - v0[2])
            nz = (v1[0] - v0[0]) * (v2[1] - v0[1]) - (v1[1] - v0[1]) * (v2[0] - v0[0])
            for idx in [n1 - 1, n2 - 1, n3 - 1]:
                accum[idx][0] += nx
                accum[idx][1] += ny
                accum[idx][2] += nz

        for vec in accum:
            length = sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)
            if length > 1e-12:
                normals.append([vec[0]/length, vec[1]/length, vec[2]/length])
            else:
                normals.append([0.0, 0.0, 1.0])

    # Extract faces
    faces = []
    for i in range(1, triangulation.NbTriangles() + 1):
        t = triangulation.Triangle(i)
        n1, n2, n3 = t.Get()
        if reverse:
            n2, n3 = n3, n2
        faces.append([n1 - 1, n2 - 1, n3 - 1])

    return tessellated_surface(vertices, normals, faces,
                              u_range=(u_min, u_max), v_range=(v_min, v_max))


# -----------------------------------------------------------------------------
# OCC Edge → Native Edge Conversion
# -----------------------------------------------------------------------------

def occ_edge_to_native(edge: TopoDS_Edge, vertex_map: dict = None):
    """Convert an OCC edge to a native BREP edge.

    Parameters
    ----------
    edge : TopoDS_Edge
        The OCC edge to convert.
    vertex_map : dict, optional
        Map from OCC vertex hash to native vertex, for topology consistency.

    Returns
    -------
    tuple
        (native_edge, start_vertex, end_vertex)
    """
    require_occ()

    if vertex_map is None:
        vertex_map = {}

    adaptor = BRepAdaptor_Curve(edge)
    curve_type = adaptor.GetType()

    # Get vertices
    v1 = topexp.FirstVertex(edge)
    v2 = topexp.LastVertex(edge)

    p1 = BRep_Tool.Pnt(v1)
    p2 = BRep_Tool.Pnt(v2)

    # Get or create native vertices
    v1_hash = hash((round(p1.X(), 8), round(p1.Y(), 8), round(p1.Z(), 8)))
    v2_hash = hash((round(p2.X(), 8), round(p2.Y(), 8), round(p2.Z(), 8)))

    if v1_hash in vertex_map:
        native_v1 = vertex_map[v1_hash]
    else:
        native_v1 = brep_vertex([p1.X(), p1.Y(), p1.Z()])
        vertex_map[v1_hash] = native_v1

    if v2_hash in vertex_map:
        native_v2 = vertex_map[v2_hash]
    else:
        native_v2 = brep_vertex([p2.X(), p2.Y(), p2.Z()])
        vertex_map[v2_hash] = native_v2

    # Get parameter bounds
    t_min, t_max = adaptor.FirstParameter(), adaptor.LastParameter()

    if curve_type == GeomAbs_Line:
        # Line edge
        native_edge = line_edge(native_v1, native_v2)
    elif curve_type == GeomAbs_Circle:
        # Circle edge
        circ = adaptor.Circle()
        center = _gp_pnt_to_point(circ.Location())
        axis = _gp_dir_to_vector(circ.Axis().Direction())
        radius = circ.Radius()
        native_edge = circle_edge(native_v1, native_v2, center, axis, radius,
                                  t_start=t_min, t_end=t_max)
    elif curve_type in (GeomAbs_BSplineCurve, GeomAbs_BezierCurve):
        # B-spline edge
        native_edge = _convert_bspline_edge(adaptor, native_v1, native_v2, t_min, t_max)
    else:
        # Fallback: sample the curve and create B-spline approximation
        native_edge = _approximate_edge_as_bspline(adaptor, native_v1, native_v2,
                                                   t_min, t_max, num_samples=20)

    return native_edge, native_v1, native_v2


def _convert_bspline_edge(adaptor, v1, v2, t_min, t_max):
    """Convert OCC B-spline curve to native B-spline edge."""
    bspl = adaptor.BSpline()

    degree = bspl.Degree()
    n_poles = bspl.NbPoles()

    # Get control points
    control_points = []
    for i in range(1, n_poles + 1):
        pnt = bspl.Pole(i)
        control_points.append([pnt.X(), pnt.Y(), pnt.Z()])

    # Get weights
    weights = None
    if bspl.IsRational():
        weights = [bspl.Weight(i) for i in range(1, n_poles + 1)]

    # Get knot vector
    knots = []
    for i in range(1, bspl.NbKnots() + 1):
        knot = bspl.Knot(i)
        mult = bspl.Multiplicity(i)
        knots.extend([knot] * mult)

    return bspline_edge(v1, v2, control_points, knots, degree,
                       weights=weights, t_start=t_min, t_end=t_max)


def _approximate_edge_as_bspline(adaptor, v1, v2, t_min, t_max, num_samples=20):
    """Approximate an arbitrary OCC curve as a B-spline edge.

    Samples the curve and fits a cubic B-spline through the points.
    """
    # Sample points along the curve
    points = []
    for i in range(num_samples):
        t = t_min + (t_max - t_min) * i / (num_samples - 1)
        pnt = adaptor.Value(t)
        points.append([pnt.X(), pnt.Y(), pnt.Z()])

    # Create cubic B-spline approximation
    # For simplicity, use clamped uniform knot vector
    degree = 3
    n = len(points)

    if n <= degree:
        # Too few points, use linear interpolation
        return line_edge(v1, v2)

    # Clamped uniform knot vector
    knots = [0.0] * (degree + 1)
    num_internal = n - degree - 1
    for i in range(1, num_internal + 1):
        knots.append(i / (num_internal + 1))
    knots.extend([1.0] * (degree + 1))

    # Use sampled points as control points (approximation)
    return bspline_edge(v1, v2, points, knots, degree, t_start=t_min, t_end=t_max)


# -----------------------------------------------------------------------------
# OCC Solid → Native BREP Conversion
# -----------------------------------------------------------------------------

def occ_solid_to_native_brep(shape: TopoDS_Shape, tessellate_fallback: bool = True):
    """Convert an OCC solid to native BREP representation.

    Parameters
    ----------
    shape : TopoDS_Shape
        The OCC shape (solid) to convert.
    tessellate_fallback : bool
        If True, use tessellation for unsupported surface types.

    Returns
    -------
    tuple
        (native_solid, topology_graph) where native_solid is the BREP solid
        and topology_graph is the TopologyGraph containing all entities.
    """
    require_occ()

    graph = TopologyGraph()
    vertex_map = {}  # OCC vertex hash -> native vertex
    edge_map = {}    # OCC edge hash -> native edge

    # Get the solid (handle both Solid and Shell inputs)
    if shape.ShapeType() == TopAbs_SOLID:
        solid = topods.Solid(shape)
    elif shape.ShapeType() == TopAbs_SHELL:
        # Wrap shell in a solid
        solid = shape
    else:
        raise ValueError(f"Expected SOLID or SHELL, got {shape.ShapeType()}")

    shells = []

    # Iterate over shells
    shell_explorer = TopExp_Explorer(shape, TopAbs_SHELL)
    while shell_explorer.More():
        occ_shell = topods.Shell(shell_explorer.Current())
        native_faces = []

        # Iterate over faces in shell
        face_explorer = TopExp_Explorer(occ_shell, TopAbs_FACE)
        while face_explorer.More():
            occ_face = topods.Face(face_explorer.Current())

            # Convert surface
            native_surface = occ_surface_to_native(occ_face, tessellate_fallback)
            if native_surface is None:
                face_explorer.Next()
                continue

            # Convert loops (wires)
            native_loops = []
            wire_explorer = TopExp_Explorer(occ_face, TopAbs_WIRE)
            while wire_explorer.More():
                occ_wire = topods.Wire(wire_explorer.Current())
                native_trims = []

                # Iterate over edges in wire
                edge_explorer = TopExp_Explorer(occ_wire, TopAbs_EDGE)
                while edge_explorer.More():
                    occ_edge = topods.Edge(edge_explorer.Current())

                    # Get or create native edge
                    edge_hash = occ_edge.__hash__()
                    if edge_hash in edge_map:
                        native_edge = edge_map[edge_hash]
                    else:
                        native_edge, v1, v2 = occ_edge_to_native(occ_edge, vertex_map)
                        edge_map[edge_hash] = native_edge
                        graph.add_vertex(v1)
                        graph.add_vertex(v2)
                        graph.add_edge(native_edge)

                    # Determine trim sense
                    sense = occ_edge.Orientation() != TopAbs_REVERSED
                    native_trim = brep_trim(native_edge, sense=sense)
                    graph.add_trim(native_trim)
                    native_trims.append(native_trim)

                    edge_explorer.Next()

                if native_trims:
                    # Determine loop type (outer vs inner)
                    # First wire is typically outer loop
                    loop_type = 'outer' if not native_loops else 'inner'
                    native_loop = brep_loop(native_trims, loop_type=loop_type)
                    graph.add_loop(native_loop)
                    native_loops.append(native_loop)

                wire_explorer.Next()

            # Create native face
            native_face = brep_face(native_surface, native_loops)
            graph.add_face(native_face)
            native_faces.append(native_face)

            face_explorer.Next()

        if native_faces:
            # Create shell (auto-detect closure)
            native_shell = brep_shell(native_faces, closed=None)
            graph.add_shell(native_shell, validate_closure=False)
            shells.append(native_shell)

        shell_explorer.Next()

    if not shells:
        raise ValueError("No shells found in OCC shape")

    # Create solid
    native_solid = brep_solid(shells)
    graph.add_solid(native_solid, validate=False)

    return native_solid, graph


# -----------------------------------------------------------------------------
# Native Surface → OCC Surface Conversion
# -----------------------------------------------------------------------------

def native_surface_to_occ(surf):
    """Convert a native analytic surface to an OCC Geom_Surface.

    Parameters
    ----------
    surf : native surface
        A native analytic surface.

    Returns
    -------
    Geom_Surface or None
        The OCC surface, or None if conversion is not supported.
    """
    require_occ()

    if not is_analytic_surface(surf):
        return None

    surface_type = surf[0]

    if surface_type == 'plane_surface':
        return _native_plane_to_occ(surf)
    elif surface_type == 'sphere_surface':
        return _native_sphere_to_occ(surf)
    elif surface_type == 'cylinder_surface':
        return _native_cylinder_to_occ(surf)
    elif surface_type == 'cone_surface':
        return _native_cone_to_occ(surf)
    elif surface_type == 'torus_surface':
        return _native_torus_to_occ(surf)
    elif surface_type == 'bspline_surface':
        return _native_bspline_to_occ(surf)
    else:
        return None


def _native_plane_to_occ(surf):
    """Convert native plane surface to OCC Geom_Plane."""
    origin = surf[1]
    meta = surf[2]
    normal = meta['normal']

    pnt = gp_Pnt(origin[0], origin[1], origin[2])
    direction = gp_Dir(normal[0], normal[1], normal[2])

    return Geom_Plane(pnt, direction)


def _native_sphere_to_occ(surf):
    """Convert native sphere surface to OCC Geom_SphericalSurface."""
    center = surf[1]
    meta = surf[2]
    radius = meta['radius']

    pnt = gp_Pnt(center[0], center[1], center[2])
    ax3 = gp_Ax3(pnt, gp_Dir(0, 0, 1))  # Default axis

    return Geom_SphericalSurface(ax3, radius)


def _native_cylinder_to_occ(surf):
    """Convert native cylinder surface to OCC Geom_CylindricalSurface."""
    origin = surf[1]
    meta = surf[2]
    axis = meta['axis']
    radius = meta['radius']

    pnt = gp_Pnt(origin[0], origin[1], origin[2])
    direction = gp_Dir(axis[0], axis[1], axis[2])
    ax3 = gp_Ax3(pnt, direction)

    return Geom_CylindricalSurface(ax3, radius)


def _native_cone_to_occ(surf):
    """Convert native cone surface to OCC Geom_ConicalSurface."""
    apex = surf[1]
    meta = surf[2]
    axis = meta['axis']
    half_angle = meta['half_angle']

    pnt = gp_Pnt(apex[0], apex[1], apex[2])
    direction = gp_Dir(axis[0], axis[1], axis[2])
    ax3 = gp_Ax3(pnt, direction)

    return Geom_ConicalSurface(ax3, half_angle, 0.0)  # Reference radius at apex


def _native_torus_to_occ(surf):
    """Convert native torus surface to OCC Geom_ToroidalSurface."""
    center = surf[1]
    meta = surf[2]
    axis = meta['axis']
    major_radius = meta['major_radius']
    minor_radius = meta['minor_radius']

    pnt = gp_Pnt(center[0], center[1], center[2])
    direction = gp_Dir(axis[0], axis[1], axis[2])
    ax3 = gp_Ax3(pnt, direction)

    return Geom_ToroidalSurface(ax3, major_radius, minor_radius)


def _native_bspline_to_occ(surf):
    """Convert native B-spline surface to OCC Geom_BSplineSurface."""
    cpts = surf[1]
    meta = surf[2]

    n_u = meta['n_u']
    n_v = meta['n_v']
    u_degree = meta['u_degree']
    v_degree = meta['v_degree']
    u_knots_flat = meta['u_knots']
    v_knots_flat = meta['v_knots']
    weights = meta['weights']

    # Create control points array
    poles = TColgp_Array2OfPnt(1, n_u, 1, n_v)
    for j in range(n_v):
        for i in range(n_u):
            cp = cpts[j][i]
            poles.SetValue(i + 1, j + 1, gp_Pnt(cp[0], cp[1], cp[2]))

    # Convert flat knot vector to knots + multiplicities
    def knots_to_occ(flat_knots):
        """Convert flat knot vector to unique knots and multiplicities."""
        unique_knots = []
        multiplicities = []
        prev = None
        for k in flat_knots:
            if prev is None or abs(k - prev) > 1e-10:
                unique_knots.append(k)
                multiplicities.append(1)
            else:
                multiplicities[-1] += 1
            prev = k
        return unique_knots, multiplicities

    u_knots, u_mults = knots_to_occ(u_knots_flat)
    v_knots, v_mults = knots_to_occ(v_knots_flat)

    # Create OCC arrays
    u_knots_arr = TColStd_Array1OfReal(1, len(u_knots))
    for i, k in enumerate(u_knots):
        u_knots_arr.SetValue(i + 1, k)

    v_knots_arr = TColStd_Array1OfReal(1, len(v_knots))
    for i, k in enumerate(v_knots):
        v_knots_arr.SetValue(i + 1, k)

    u_mults_arr = TColStd_Array1OfInteger(1, len(u_mults))
    for i, m in enumerate(u_mults):
        u_mults_arr.SetValue(i + 1, m)

    v_mults_arr = TColStd_Array1OfInteger(1, len(v_mults))
    for i, m in enumerate(v_mults):
        v_mults_arr.SetValue(i + 1, m)

    # Check if rational (non-uniform weights)
    is_rational = False
    for row in weights:
        for w in row:
            if abs(w - 1.0) > 1e-10:
                is_rational = True
                break

    if is_rational:
        # Create weights array
        from OCC.Core.TColStd import TColStd_Array2OfReal
        weights_arr = TColStd_Array2OfReal(1, n_u, 1, n_v)
        for j in range(n_v):
            for i in range(n_u):
                weights_arr.SetValue(i + 1, j + 1, weights[j][i])

        return Geom_BSplineSurface(poles, weights_arr,
                                   u_knots_arr, v_knots_arr,
                                   u_mults_arr, v_mults_arr,
                                   u_degree, v_degree)
    else:
        return Geom_BSplineSurface(poles,
                                   u_knots_arr, v_knots_arr,
                                   u_mults_arr, v_mults_arr,
                                   u_degree, v_degree)


# -----------------------------------------------------------------------------
# Native BREP → OCC Conversion (for boolean operations)
# -----------------------------------------------------------------------------

def native_vertex_to_occ(vertex):
    """Convert a native BREP vertex to an OCC TopoDS_Vertex.

    Parameters
    ----------
    vertex : brep_vertex
        The native vertex to convert.

    Returns
    -------
    TopoDS_Vertex
        The OCC vertex.
    """
    require_occ()
    from yapcad.native_brep import vertex_location

    loc = vertex_location(vertex)
    pnt = gp_Pnt(loc[0], loc[1], loc[2])
    return BRepBuilderAPI_MakeVertex(pnt).Vertex()


def native_edge_to_occ(edge, graph):
    """Convert a native BREP edge to an OCC TopoDS_Edge.

    Parameters
    ----------
    edge : brep_edge
        The native edge to convert.
    graph : TopologyGraph
        The topology graph containing the edge's vertices.

    Returns
    -------
    TopoDS_Edge
        The OCC edge.
    """
    require_occ()
    from yapcad.native_brep import (
        edge_curve_type, edge_vertices, edge_curve_params, vertex_location
    )

    curve_type = edge_curve_type(edge)
    start_vid, end_vid = edge_vertices(edge)
    params = edge_curve_params(edge)

    # Get vertex locations - handle missing vertices
    if start_vid is None or start_vid not in graph.vertices:
        raise ValueError(f"Edge missing start vertex: {start_vid}")
    if end_vid is None or end_vid not in graph.vertices:
        raise ValueError(f"Edge missing end vertex: {end_vid}")

    start_loc = vertex_location(graph.vertices[start_vid])
    end_loc = vertex_location(graph.vertices[end_vid])

    p1 = gp_Pnt(start_loc[0], start_loc[1], start_loc[2])
    p2 = gp_Pnt(end_loc[0], end_loc[1], end_loc[2])

    if curve_type == 'line':
        # Line edge
        edge_builder = BRepBuilderAPI_MakeEdge(p1, p2)

    elif curve_type == 'circle':
        # Circle/arc edge
        center = params.get('center', [0, 0, 0])
        axis = params.get('axis', [0, 0, 1])
        radius = params.get('radius', 1.0)
        t_start = params.get('t_start', 0.0)
        t_end = params.get('t_end', 2 * 3.141592653589793)

        center_pnt = gp_Pnt(center[0], center[1], center[2])
        axis_dir = gp_Dir(axis[0], axis[1], axis[2])
        ax2 = gp_Ax2(center_pnt, axis_dir)

        # Check if this is a full circle (start == end vertex)
        is_full_circle = (start_vid == end_vid) or (p1.Distance(p2) < 1e-9)

        try:
            from OCC.Core.Geom import Geom_Circle
            from OCC.Core.gp import gp_Circ

            # Create a circle
            circ = gp_Circ(ax2, radius)
            geom_circle = Geom_Circle(circ)

            if is_full_circle:
                # Full circle - use parameter range
                edge_builder = BRepBuilderAPI_MakeEdge(geom_circle, t_start, t_end)
            else:
                # Arc - use start/end points
                edge_builder = BRepBuilderAPI_MakeEdge(geom_circle, p1, p2)

            if not edge_builder.IsDone():
                # Fallback: try arc maker
                arc_maker = GC_MakeArcOfCircle(p1, p2, center_pnt)
                if arc_maker.IsDone():
                    edge_builder = BRepBuilderAPI_MakeEdge(arc_maker.Value())
                else:
                    edge_builder = BRepBuilderAPI_MakeEdge(p1, p2)
        except Exception:
            edge_builder = BRepBuilderAPI_MakeEdge(p1, p2)

    elif curve_type == 'arc':
        # Arc with bulge or center
        if 'center' in params:
            center = params['center']
            center_pnt = gp_Pnt(center[0], center[1], center[2])
            try:
                arc_maker = GC_MakeArcOfCircle(p1, p2, center_pnt)
                if arc_maker.IsDone():
                    edge_builder = BRepBuilderAPI_MakeEdge(arc_maker.Value())
                else:
                    edge_builder = BRepBuilderAPI_MakeEdge(p1, p2)
            except Exception:
                edge_builder = BRepBuilderAPI_MakeEdge(p1, p2)
        else:
            # Bulge-based arc - compute center
            edge_builder = BRepBuilderAPI_MakeEdge(p1, p2)

    elif curve_type == 'bspline':
        # B-spline edge
        control_points = params.get('control_points', [])
        knots = params.get('knots', [])
        degree = params.get('degree', 3)
        weights = params.get('weights', None)

        if len(control_points) < 2:
            edge_builder = BRepBuilderAPI_MakeEdge(p1, p2)
        else:
            try:
                occ_curve = _native_bspline_curve_to_occ(
                    control_points, knots, degree, weights
                )
                if occ_curve is not None:
                    edge_builder = BRepBuilderAPI_MakeEdge(occ_curve)
                else:
                    edge_builder = BRepBuilderAPI_MakeEdge(p1, p2)
            except Exception:
                edge_builder = BRepBuilderAPI_MakeEdge(p1, p2)

    else:
        # Unknown curve type - fallback to line
        edge_builder = BRepBuilderAPI_MakeEdge(p1, p2)

    if not edge_builder.IsDone():
        # Fallback to simple line
        edge_builder = BRepBuilderAPI_MakeEdge(p1, p2)

    try:
        return edge_builder.Edge()
    except RuntimeError:
        # Edge builder failed - return a minimal degenerate edge
        # Create a line edge as last resort
        edge_builder = BRepBuilderAPI_MakeEdge(p1, p2)
        if edge_builder.IsDone():
            return edge_builder.Edge()
        # If even that fails, create edge from a very small offset
        p2_offset = gp_Pnt(p1.X() + 1e-6, p1.Y(), p1.Z())
        return BRepBuilderAPI_MakeEdge(p1, p2_offset).Edge()


def _native_bspline_curve_to_occ(control_points, knots, degree, weights=None):
    """Convert native B-spline curve data to OCC Geom_BSplineCurve."""
    require_occ()

    n = len(control_points)
    if n < 2:
        return None

    # Create control points array
    poles = TColgp_Array1OfPnt(1, n)
    for i, cp in enumerate(control_points):
        poles.SetValue(i + 1, gp_Pnt(cp[0], cp[1], cp[2]))

    # Convert flat knot vector to unique knots + multiplicities
    unique_knots = []
    multiplicities = []
    prev = None
    for k in knots:
        if prev is None or abs(k - prev) > 1e-10:
            unique_knots.append(k)
            multiplicities.append(1)
        else:
            multiplicities[-1] += 1
        prev = k

    knots_arr = TColStd_Array1OfReal(1, len(unique_knots))
    for i, k in enumerate(unique_knots):
        knots_arr.SetValue(i + 1, k)

    mults_arr = TColStd_Array1OfInteger(1, len(multiplicities))
    for i, m in enumerate(multiplicities):
        mults_arr.SetValue(i + 1, m)

    # Check if rational
    is_rational = False
    if weights is not None:
        for w in weights:
            if abs(w - 1.0) > 1e-10:
                is_rational = True
                break

    if is_rational and weights is not None:
        weights_arr = TColStd_Array1OfReal(1, n)
        for i, w in enumerate(weights):
            weights_arr.SetValue(i + 1, w)
        return Geom_BSplineCurve(poles, weights_arr, knots_arr, mults_arr, degree)
    else:
        return Geom_BSplineCurve(poles, knots_arr, mults_arr, degree)


def native_loop_to_occ_wire(loop, graph, edge_cache=None):
    """Convert a native BREP loop to an OCC TopoDS_Wire.

    Parameters
    ----------
    loop : brep_loop
        The native loop to convert.
    graph : TopologyGraph
        The topology graph containing the loop's edges.
    edge_cache : dict, optional
        Cache of edge_id -> TopoDS_Edge for reuse.

    Returns
    -------
    TopoDS_Wire
        The OCC wire.
    """
    require_occ()
    from yapcad.native_brep import is_brep_loop, trim_edge_id

    if edge_cache is None:
        edge_cache = {}

    # Get trim IDs from loop
    trim_ids = loop[1]  # ['brep_loop', trim_refs, metadata]

    wire_builder = BRepBuilderAPI_MakeWire()

    for tid in trim_ids:
        trim = graph.trims.get(tid)
        if trim is None:
            continue

        edge_id = trim_edge_id(trim)
        edge = graph.edges.get(edge_id)
        if edge is None:
            continue

        # Get or create OCC edge
        if edge_id in edge_cache:
            occ_edge = edge_cache[edge_id]
        else:
            occ_edge = native_edge_to_occ(edge, graph)
            edge_cache[edge_id] = occ_edge

        # Add edge to wire (handle sense/orientation)
        sense = trim[2].get('sense', True)
        if sense:
            wire_builder.Add(occ_edge)
        else:
            # Reversed edge
            reversed_edge = occ_edge.Reversed()
            wire_builder.Add(topods.Edge(reversed_edge))

    if wire_builder.IsDone():
        return wire_builder.Wire()
    else:
        return None


def native_face_to_occ(face, graph, edge_cache=None):
    """Convert a native BREP face to an OCC TopoDS_Face.

    Parameters
    ----------
    face : brep_face
        The native face to convert.
    graph : TopologyGraph
        The topology graph containing the face's loops.
    edge_cache : dict, optional
        Cache of edge_id -> TopoDS_Edge for reuse.

    Returns
    -------
    TopoDS_Face
        The OCC face, or None if conversion fails.
    """
    require_occ()
    from math import pi

    if edge_cache is None:
        edge_cache = {}

    # Get surface and loops
    surface = face[1]  # ['brep_face', surface, metadata]
    meta = face[2]
    loop_ids = meta.get('loop_ids', [])
    face_sense = meta.get('sense', True)

    # Convert surface to OCC
    occ_surface = native_surface_to_occ(surface)
    if occ_surface is None:
        return None

    # Check if this is a complete analytic surface that should be created
    # without wire boundaries (e.g., full sphere, full torus)
    # These surfaces have singularities (poles) that make wire reconstruction complex
    surface_type = surface[0] if isinstance(surface, (list, tuple)) else None
    is_complete_surface = False

    if surface_type == 'sphere_surface':
        # Check if it's a full sphere (u: 0 to 2π, v: -π/2 to π/2)
        surf_meta = surface[2] if len(surface) > 2 else {}
        u_range = surf_meta.get('u_range', (0, 2 * pi))
        v_range = surf_meta.get('v_range', (-pi / 2, pi / 2))
        if (abs(u_range[1] - u_range[0] - 2 * pi) < 0.01 and
                abs(v_range[1] - v_range[0] - pi) < 0.01):
            is_complete_surface = True

    if is_complete_surface:
        # Create face from surface natural bounds, ignoring wires
        # (wires on surfaces with singularities are complex to reconstruct)
        try:
            face_builder = BRepBuilderAPI_MakeFace(occ_surface, 1e-6)
            if face_builder.IsDone():
                occ_face = face_builder.Face()
                if not face_sense:
                    occ_face = topods.Face(occ_face.Reversed())
                return occ_face
        except Exception:
            pass
        return None

    # Build wires from loops
    wires = []
    for lid in loop_ids:
        loop = graph.loops.get(lid)
        if loop is None:
            continue
        wire = native_loop_to_occ_wire(loop, graph, edge_cache)
        if wire is not None:
            wires.append(wire)

    if not wires:
        # No wires - create face from surface bounds
        try:
            face_builder = BRepBuilderAPI_MakeFace(occ_surface, 1e-6)
            if face_builder.IsDone():
                return face_builder.Face()
        except Exception:
            pass
        return None

    # First wire is outer boundary
    try:
        face_builder = BRepBuilderAPI_MakeFace(occ_surface, wires[0], True)
        if not face_builder.IsDone():
            return None

        # Add inner wires (holes)
        for inner_wire in wires[1:]:
            face_builder.Add(inner_wire)

        occ_face = face_builder.Face()

        # Handle face sense
        if not face_sense:
            occ_face = topods.Face(occ_face.Reversed())

        return occ_face

    except Exception:
        return None


def native_brep_to_occ(graph, fix_shape=True):
    """Convert a native BREP TopologyGraph to an OCC TopoDS_Solid.

    This function converts a complete native BREP representation back to
    an OCC solid, which can then be used for boolean operations or export.

    Parameters
    ----------
    graph : TopologyGraph
        The native BREP topology graph.
    fix_shape : bool, optional
        If True, apply shape fixing to ensure a valid solid. Default True.

    Returns
    -------
    TopoDS_Solid or TopoDS_Shell
        The OCC solid (or shell if conversion to solid fails).
    """
    require_occ()
    from yapcad.native_brep import TopologyGraph

    if not isinstance(graph, TopologyGraph):
        raise ValueError("Expected a TopologyGraph")

    edge_cache = {}

    # Convert all faces
    occ_faces = []
    for fid, face in graph.faces.items():
        occ_face = native_face_to_occ(face, graph, edge_cache)
        if occ_face is not None:
            occ_faces.append(occ_face)

    if not occ_faces:
        raise ValueError("No faces could be converted")

    # Sew faces into shell(s)
    sewing = BRepBuilderAPI_Sewing(1e-6)
    for occ_face in occ_faces:
        sewing.Add(occ_face)

    sewing.Perform()
    sewn_shape = sewing.SewedShape()

    # Try to create a solid from the sewn shell
    if sewn_shape.ShapeType() == TopAbs_SHELL:
        shell = topods.Shell(sewn_shape)

        if fix_shape:
            # Fix the shell (may fail for some geometries)
            try:
                shell_fixer = ShapeFix_Shell(shell)
                shell_fixer.Perform()
                shell = shell_fixer.Shell()
            except RuntimeError:
                # Shell fixing failed, use unfixed shell
                pass

        # Create solid from shell
        solid_builder = BRepBuilderAPI_MakeSolid(shell)
        if solid_builder.IsDone():
            solid = solid_builder.Solid()

            if fix_shape:
                # Fix the solid (may fail for some geometries)
                try:
                    solid_fixer = ShapeFix_Solid(solid)
                    solid_fixer.Perform()
                    solid = solid_fixer.Solid()
                except RuntimeError:
                    # Solid fixing failed, use unfixed solid
                    pass

            return solid
        else:
            # Return shell if solid creation fails
            return shell

    elif sewn_shape.ShapeType() == TopAbs_SOLID:
        solid = topods.Solid(sewn_shape)
        if fix_shape:
            try:
                solid_fixer = ShapeFix_Solid(solid)
                solid_fixer.Perform()
                solid = solid_fixer.Solid()
            except RuntimeError:
                # Solid fixing failed, use unfixed solid
                pass
        return solid

    else:
        # Return whatever we got
        return sewn_shape


def native_brep_to_occ_from_solid(yapcad_solid, fix_shape=True):
    """Convert a yapCAD solid's native BREP to OCC representation.

    Convenience function that extracts the native BREP from a yapCAD solid
    and converts it to OCC format.

    Parameters
    ----------
    yapcad_solid : solid
        A yapCAD solid with native BREP data attached.
    fix_shape : bool, optional
        If True, apply shape fixing to ensure a valid solid. Default True.

    Returns
    -------
    TopoDS_Solid or TopoDS_Shell or None
        The OCC shape, or None if no native BREP is attached.
    """
    require_occ()
    from yapcad.native_brep import native_brep_from_solid, has_native_brep

    if not has_native_brep(yapcad_solid):
        return None

    graph = native_brep_from_solid(yapcad_solid)
    if graph is None:
        return None

    return native_brep_to_occ(graph, fix_shape=fix_shape)


__all__ = [
    # OCC availability
    'occ_available',
    'require_occ',
    # OCC → Native conversion
    'occ_surface_to_native',
    'occ_edge_to_native',
    'occ_solid_to_native_brep',
    # Native → OCC conversion
    'native_surface_to_occ',
    'native_vertex_to_occ',
    'native_edge_to_occ',
    'native_loop_to_occ_wire',
    'native_face_to_occ',
    'native_brep_to_occ',
    'native_brep_to_occ_from_solid',
]
