"""Native BREP topology graph for yapCAD.

This module provides native (OCC-independent) BREP data structures for
representing boundary representations of 3D solids. These structures can be
used alongside or instead of OCC-based BREP data.

Topology hierarchy:
- BrepVertex: point in 3D space with tolerance
- BrepEdge: parametric curve segment bounded by vertices
- BrepTrim: oriented edge reference with UV parameterization on a face
- BrepLoop: closed sequence of trims forming a boundary
- BrepFace: surface bounded by loops
- BrepShell: connected set of faces
- BrepSolid: collection of shells (outer + inner voids)

Design principle: Vertices are the source of truth for endpoint positions.
Edges store curve *type* and *parameters*, not absolute endpoint coordinates.
This ensures transformations only need to update vertex locations.

Copyright (c) 2025 Richard DeVaul
MIT License
"""

from copy import deepcopy
from math import sqrt, atan2, sin, cos, pi
import uuid

from yapcad.geom import point, ispoint, vclose, epsilon, line, arc, add, sub, scale3


# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

def _generate_id():
    """Generate a unique ID for a BREP entity."""
    return str(uuid.uuid4())


def _normalize_point(p):
    """Convert various point representations to a standard yapCAD point."""
    if ispoint(p):
        return deepcopy(p)
    elif isinstance(p, (list, tuple)) and len(p) >= 3:
        return point(p[0], p[1], p[2])
    else:
        return point(p)


def _point_distance(p1, p2):
    """Compute distance between two points."""
    dx = p1[0] - p2[0]
    dy = p1[1] - p2[1]
    dz = p1[2] - p2[2]
    return sqrt(dx*dx + dy*dy + dz*dz)


def _midpoint(p1, p2):
    """Compute midpoint between two points."""
    return point(
        (p1[0] + p2[0]) / 2,
        (p1[1] + p2[1]) / 2,
        (p1[2] + p2[2]) / 2
    )


# -----------------------------------------------------------------------------
# BrepVertex
# -----------------------------------------------------------------------------

def brep_vertex(location, *, tolerance=epsilon, tags=None):
    """Create a BREP vertex.

    Parameters
    ----------
    location : point
        World-space position of the vertex.
    tolerance : float, optional
        Geometric tolerance for coincidence tests.
    tags : dict, optional
        User metadata.

    Returns
    -------
    list
        Vertex definition: ['brep_vertex', point, metadata_dict]
    """
    loc = _normalize_point(location)

    meta = {
        'id': _generate_id(),
        'tolerance': float(tolerance),
        'edges': [],  # IDs of incident edges
        'tags': tags or {},
    }

    return ['brep_vertex', loc, meta]


def is_brep_vertex(obj):
    """Return True if obj is a BREP vertex."""
    return (isinstance(obj, list) and len(obj) == 3 and obj[0] == 'brep_vertex'
            and ispoint(obj[1]) and isinstance(obj[2], dict))


def vertex_location(v):
    """Return the location point of a vertex."""
    if not is_brep_vertex(v):
        raise ValueError("Not a BREP vertex")
    return deepcopy(v[1])


def vertex_id(v):
    """Return the ID of a vertex."""
    if not is_brep_vertex(v):
        raise ValueError("Not a BREP vertex")
    return v[2].get('id')


def set_vertex_location(v, location):
    """Set the location of a vertex (for transformations)."""
    if not is_brep_vertex(v):
        raise ValueError("Not a BREP vertex")
    v[1] = _normalize_point(location)


def vertices_coincident(v1, v2, tolerance=None):
    """Return True if two vertices are coincident within tolerance."""
    if not (is_brep_vertex(v1) and is_brep_vertex(v2)):
        raise ValueError("Arguments must be BREP vertices")

    if tolerance is None:
        tolerance = max(v1[2].get('tolerance', epsilon),
                       v2[2].get('tolerance', epsilon))

    return _point_distance(v1[1], v2[1]) <= tolerance


# -----------------------------------------------------------------------------
# BrepEdge - Parametric Curve Design
# -----------------------------------------------------------------------------

# Supported curve types and their parameters:
#
# 'line':
#   No parameters needed - curve is defined by start/end vertices
#
# 'arc':
#   'bulge': float - signed distance from chord midpoint to arc midpoint
#                    positive = CCW arc, negative = CW arc
#   OR
#   'center': point - arc center (must be transformed with vertices)
#   'orientation': int - 1 for CCW, -1 for CW
#
# 'ellipse':
#   'center': point - ellipse center (must be transformed)
#   'semi_major': float
#   'semi_minor': float
#   'rotation': float - rotation angle in degrees
#   'start_angle': float - start angle in degrees
#   'end_angle': float - end angle in degrees
#
# 'nurbs':
#   'control_points': list of points (must be transformed)
#   'degree': int
#   'knots': list of floats
#   'weights': list of floats (optional)

def brep_edge(curve_type, start_vertex, end_vertex, *,
              curve_params=None, sense=True, tolerance=epsilon, tags=None):
    """Create a BREP edge with parametric curve definition.

    The edge stores curve type and parameters, NOT absolute endpoint
    coordinates. Endpoint positions are derived from the referenced vertices.

    Parameters
    ----------
    curve_type : str
        Type of curve: 'line', 'arc', 'ellipse', 'nurbs'.
    start_vertex : brep_vertex or str
        Start vertex or vertex ID.
    end_vertex : brep_vertex or str
        End vertex or vertex ID.
    curve_params : dict, optional
        Type-specific curve parameters (see module docstring).
    sense : bool, optional
        True if edge direction agrees with curve direction.
    tolerance : float, optional
        Geometric tolerance.
    tags : dict, optional
        User metadata.

    Returns
    -------
    list
        Edge definition: ['brep_edge', curve_type, metadata_dict]
    """
    # Get vertex IDs
    start_id = vertex_id(start_vertex) if is_brep_vertex(start_vertex) else start_vertex
    end_id = vertex_id(end_vertex) if is_brep_vertex(end_vertex) else end_vertex

    meta = {
        'id': _generate_id(),
        'start_vertex_id': start_id,
        'end_vertex_id': end_id,
        'curve_params': deepcopy(curve_params) if curve_params else {},
        'sense': bool(sense),
        'tolerance': float(tolerance),
        'faces': [],  # IDs of faces using this edge
        'tags': tags or {},
    }

    return ['brep_edge', curve_type, meta]


def is_brep_edge(obj):
    """Return True if obj is a BREP edge."""
    return (isinstance(obj, list) and len(obj) == 3 and obj[0] == 'brep_edge'
            and isinstance(obj[1], str) and isinstance(obj[2], dict))


def edge_curve_type(e):
    """Return the curve type of an edge."""
    if not is_brep_edge(e):
        raise ValueError("Not a BREP edge")
    return e[1]


def edge_curve_params(e):
    """Return the curve parameters of an edge."""
    if not is_brep_edge(e):
        raise ValueError("Not a BREP edge")
    return deepcopy(e[2].get('curve_params', {}))


def edge_id(e):
    """Return the ID of an edge."""
    if not is_brep_edge(e):
        raise ValueError("Not a BREP edge")
    return e[2].get('id')


def edge_vertices(e):
    """Return (start_vertex_id, end_vertex_id) for an edge."""
    if not is_brep_edge(e):
        raise ValueError("Not a BREP edge")
    return (e[2].get('start_vertex_id'), e[2].get('end_vertex_id'))


def evaluate_edge_curve(edge, start_point, end_point):
    """Construct a yapCAD curve from edge definition and vertex positions.

    Parameters
    ----------
    edge : brep_edge
        The edge definition.
    start_point : point
        Position of start vertex.
    end_point : point
        Position of end vertex.

    Returns
    -------
    curve
        A yapCAD curve (line, arc, ellipse, etc.).
    """
    if not is_brep_edge(edge):
        raise ValueError("Not a BREP edge")

    curve_type = edge[1]
    params = edge[2].get('curve_params', {})

    if curve_type == 'line':
        return line(start_point, end_point)

    elif curve_type == 'arc':
        if 'bulge' in params:
            # Compute arc from bulge factor
            bulge = params['bulge']
            if abs(bulge) < epsilon:
                return line(start_point, end_point)

            # Bulge = tan(theta/4) where theta is the included angle
            # Positive bulge = arc bulges to the left (CCW)
            mid = _midpoint(start_point, end_point)
            chord_len = _point_distance(start_point, end_point)

            # Perpendicular direction (in XY plane for now)
            dx = end_point[0] - start_point[0]
            dy = end_point[1] - start_point[1]
            perp_len = sqrt(dx*dx + dy*dy)
            if perp_len < epsilon:
                return line(start_point, end_point)

            perp = point(-dy/perp_len, dx/perp_len, 0)

            # Arc midpoint offset from chord midpoint
            arc_height = bulge * chord_len / 2
            arc_mid = add(mid, scale3(perp, arc_height))

            # Compute center and radius
            # Using the relationship: radius = (chord^2/4 + height^2) / (2*height)
            h = abs(arc_height)
            if h < epsilon:
                return line(start_point, end_point)
            radius = (chord_len*chord_len/4 + h*h) / (2*h)

            # Center is offset from arc midpoint
            center_offset = radius - h
            if bulge > 0:
                center = sub(arc_mid, scale3(perp, center_offset))
            else:
                center = add(arc_mid, scale3(perp, center_offset))

            # Compute start/end angles
            start_angle = atan2(start_point[1] - center[1],
                               start_point[0] - center[0]) * 180 / pi
            end_angle = atan2(end_point[1] - center[1],
                             end_point[0] - center[0]) * 180 / pi

            orientation = 1 if bulge > 0 else -1
            return arc(center, [radius, start_angle, end_angle, orientation])

        elif 'center' in params:
            # Arc defined by center point
            center = params['center']
            orientation = params.get('orientation', 1)
            radius = _point_distance(center, start_point)

            start_angle = atan2(start_point[1] - center[1],
                               start_point[0] - center[0]) * 180 / pi
            end_angle = atan2(end_point[1] - center[1],
                             end_point[0] - center[0]) * 180 / pi

            return arc(center, [radius, start_angle, end_angle, orientation])

        else:
            # Default to line if no arc params
            return line(start_point, end_point)

    elif curve_type == 'ellipse':
        from yapcad.geom import ellipse
        center = params.get('center', _midpoint(start_point, end_point))
        semi_major = params.get('semi_major', _point_distance(start_point, end_point) / 2)
        semi_minor = params.get('semi_minor', semi_major)
        rotation = params.get('rotation', 0.0)
        start_angle = params.get('start_angle', 0)
        end_angle = params.get('end_angle', 360)

        return ellipse(center, semi_major, semi_minor,
                      rotation=rotation, start=start_angle, end=end_angle)

    else:
        # Unknown curve type - default to line
        return line(start_point, end_point)


# Convenience constructors for common edge types

def line_edge(start_vertex, end_vertex, **kwargs):
    """Create a line edge between two vertices."""
    return brep_edge('line', start_vertex, end_vertex, **kwargs)


def arc_edge(start_vertex, end_vertex, *, bulge=None, center=None, orientation=1, **kwargs):
    """Create an arc edge between two vertices.

    Specify either bulge OR center, not both.

    Parameters
    ----------
    bulge : float, optional
        Signed arc height factor. Positive = CCW, negative = CW.
    center : point, optional
        Arc center point.
    orientation : int, optional
        Arc direction: 1 for CCW, -1 for CW (used with center).
    """
    if bulge is not None:
        params = {'bulge': float(bulge)}
    elif center is not None:
        params = {'center': _normalize_point(center), 'orientation': int(orientation)}
    else:
        raise ValueError("Must specify either bulge or center for arc edge")

    return brep_edge('arc', start_vertex, end_vertex, curve_params=params, **kwargs)


def circle_edge(start_vertex, end_vertex, center, axis, radius, *,
                t_start=0.0, t_end=None, **kwargs):
    """Create a circular arc edge between two vertices.

    Parameters
    ----------
    start_vertex, end_vertex : brep_vertex
        The endpoints of the edge.
    center : point
        Center of the circle.
    axis : vector
        Normal to the plane of the circle.
    radius : float
        Radius of the circle.
    t_start : float, optional
        Start parameter (angle in radians). Default 0.
    t_end : float, optional
        End parameter (angle in radians). Default computed from endpoints.
    """
    from math import pi

    params = {
        'center': _normalize_point(center),
        'axis': _normalize_vector(axis),
        'radius': float(radius),
        't_start': float(t_start),
        't_end': float(t_end) if t_end is not None else 2 * pi,
    }

    return brep_edge('circle', start_vertex, end_vertex, curve_params=params, **kwargs)


def bspline_edge(start_vertex, end_vertex, control_points, knots, degree, *,
                 weights=None, t_start=None, t_end=None, **kwargs):
    """Create a B-spline edge between two vertices.

    Parameters
    ----------
    start_vertex, end_vertex : brep_vertex
        The endpoints of the edge.
    control_points : list of points
        Control points for the B-spline.
    knots : list of float
        Knot vector.
    degree : int
        Degree of the B-spline.
    weights : list of float, optional
        Weights for rational B-splines. If None, uniform weights (1.0).
    t_start : float, optional
        Start parameter. Default is knots[degree].
    t_end : float, optional
        End parameter. Default is knots[-degree-1].
    """
    n = len(control_points)
    expected_knots = n + degree + 1
    if len(knots) != expected_knots:
        raise ValueError(f"Expected {expected_knots} knots, got {len(knots)}")

    cpts = [_normalize_point(cp) for cp in control_points]

    if weights is None:
        w = [1.0] * n
    else:
        w = [float(wi) for wi in weights]

    params = {
        'control_points': cpts,
        'knots': [float(k) for k in knots],
        'degree': int(degree),
        'weights': w,
        't_start': float(t_start) if t_start is not None else knots[degree],
        't_end': float(t_end) if t_end is not None else knots[n],
    }

    return brep_edge('bspline', start_vertex, end_vertex, curve_params=params, **kwargs)


def _normalize_vector(v):
    """Normalize a vector to unit length with 4 components."""
    from math import sqrt
    x, y, z = float(v[0]), float(v[1]), float(v[2])
    mag = sqrt(x*x + y*y + z*z)
    if mag < 1e-12:
        return [0.0, 0.0, 1.0, 0.0]
    return [x/mag, y/mag, z/mag, 0.0]


# -----------------------------------------------------------------------------
# Edge Transformation Support
# -----------------------------------------------------------------------------

def transform_edge_params(edge, transform_point_func):
    """Transform internal geometry points in edge curve parameters.

    For edges with internal geometry (like arc centers), this function
    applies a transformation to those points. Call this after transforming
    the vertices.

    Parameters
    ----------
    edge : brep_edge
        The edge to transform (modified in place).
    transform_point_func : callable
        Function that takes a point and returns transformed point.
    """
    if not is_brep_edge(edge):
        raise ValueError("Not a BREP edge")

    curve_type = edge[1]
    params = edge[2].get('curve_params', {})

    if curve_type == 'arc' and 'center' in params:
        params['center'] = transform_point_func(params['center'])

    elif curve_type == 'ellipse' and 'center' in params:
        params['center'] = transform_point_func(params['center'])

    elif curve_type == 'nurbs' and 'control_points' in params:
        params['control_points'] = [
            transform_point_func(pt) for pt in params['control_points']
        ]


# -----------------------------------------------------------------------------
# BrepTrim
# -----------------------------------------------------------------------------

def brep_trim(edge, *, sense=True, uv_curve=None, tags=None):
    """Create a BREP trim (oriented edge use on a face).

    A trim represents the use of an edge on a particular face, with
    optional UV parameterization for trimmed surface representation.

    Parameters
    ----------
    edge : brep_edge or str
        The edge or edge ID being trimmed.
    sense : bool, optional
        True if trim direction agrees with edge direction.
    uv_curve : curve, optional
        Curve in the (u, v) parameter space of the face.
    tags : dict, optional
        User metadata.

    Returns
    -------
    list
        Trim definition: ['brep_trim', edge_ref, metadata_dict]
    """
    edge_ref = edge_id(edge) if is_brep_edge(edge) else edge

    meta = {
        'id': _generate_id(),
        'edge_id': edge_ref,
        'sense': bool(sense),
        'uv_curve': deepcopy(uv_curve) if uv_curve else None,
        'tags': tags or {},
    }

    return ['brep_trim', edge_ref, meta]


def is_brep_trim(obj):
    """Return True if obj is a BREP trim."""
    return (isinstance(obj, list) and len(obj) == 3 and obj[0] == 'brep_trim'
            and isinstance(obj[2], dict))


def trim_id(t):
    """Return the ID of a trim."""
    if not is_brep_trim(t):
        raise ValueError("Not a BREP trim")
    return t[2].get('id')


def trim_edge_id(t):
    """Return the edge ID referenced by a trim."""
    if not is_brep_trim(t):
        raise ValueError("Not a BREP trim")
    return t[2].get('edge_id')


# -----------------------------------------------------------------------------
# BrepLoop
# -----------------------------------------------------------------------------

def brep_loop(trims, *, loop_type='outer', tags=None):
    """Create a BREP loop (closed boundary of a face).

    Parameters
    ----------
    trims : list
        Ordered list of trims forming the closed boundary.
    loop_type : str, optional
        'outer' for outer boundary, 'inner' for hole.
    tags : dict, optional
        User metadata.

    Returns
    -------
    list
        Loop definition: ['brep_loop', trim_refs, metadata_dict]
    """
    # Get trim IDs
    trim_refs = []
    for t in trims:
        if is_brep_trim(t):
            trim_refs.append(trim_id(t))
        else:
            trim_refs.append(t)  # Assume it's already an ID

    meta = {
        'id': _generate_id(),
        'loop_type': loop_type,
        'tags': tags or {},
    }

    return ['brep_loop', trim_refs, meta]


def is_brep_loop(obj):
    """Return True if obj is a BREP loop."""
    return (isinstance(obj, list) and len(obj) == 3 and obj[0] == 'brep_loop'
            and isinstance(obj[1], list) and isinstance(obj[2], dict))


def loop_id(loop):
    """Return the ID of a loop."""
    if not is_brep_loop(loop):
        raise ValueError("Not a BREP loop")
    return loop[2].get('id')


def loop_trims(loop):
    """Return the list of trim IDs in a loop."""
    if not is_brep_loop(loop):
        raise ValueError("Not a BREP loop")
    return deepcopy(loop[1])


def loop_type(loop):
    """Return the loop type ('outer' or 'inner')."""
    if not is_brep_loop(loop):
        raise ValueError("Not a BREP loop")
    return loop[2].get('loop_type', 'outer')


# -----------------------------------------------------------------------------
# BrepFace
# -----------------------------------------------------------------------------

def brep_face(surface, loops, *, sense=True, tags=None):
    """Create a BREP face.

    Parameters
    ----------
    surface : analytic_surface or tessellated surface
        The underlying surface geometry.
    loops : list
        List of loops bounding the face (first is outer, rest are holes).
    sense : bool, optional
        True if face normal agrees with surface normal.
    tags : dict, optional
        User metadata.

    Returns
    -------
    list
        Face definition: ['brep_face', surface, metadata_dict]
    """
    # Get loop IDs
    loop_refs = []
    for loop in loops:
        if is_brep_loop(loop):
            loop_refs.append(loop_id(loop))
        else:
            loop_refs.append(loop)  # Assume it's already an ID

    meta = {
        'id': _generate_id(),
        'loop_ids': loop_refs,
        'sense': bool(sense),
        'shell_id': None,  # Will be set when added to a shell
        'tags': tags or {},
    }

    return ['brep_face', deepcopy(surface), meta]


def is_brep_face(obj):
    """Return True if obj is a BREP face."""
    return (isinstance(obj, list) and len(obj) == 3 and obj[0] == 'brep_face'
            and isinstance(obj[2], dict))


def face_surface(f):
    """Return the surface geometry of a face."""
    if not is_brep_face(f):
        raise ValueError("Not a BREP face")
    return deepcopy(f[1])


def face_id(f):
    """Return the ID of a face."""
    if not is_brep_face(f):
        raise ValueError("Not a BREP face")
    return f[2].get('id')


def face_loops(f):
    """Return the list of loop IDs in a face."""
    if not is_brep_face(f):
        raise ValueError("Not a BREP face")
    return deepcopy(f[2].get('loop_ids', []))


# -----------------------------------------------------------------------------
# Shell Closure Validation Exception
# -----------------------------------------------------------------------------

class ShellClosureError(ValueError):
    """Exception raised when shell closure validation fails."""

    def __init__(self, message, details=None):
        super().__init__(message)
        self.details = details or {}


class SolidValidationError(ValueError):
    """Exception raised when solid validation fails."""

    def __init__(self, message, details=None):
        super().__init__(message)
        self.details = details or {}


# -----------------------------------------------------------------------------
# BrepShell
# -----------------------------------------------------------------------------

def brep_shell(faces, *, shell_type='outer', closed=None, tags=None):
    """Create a BREP shell.

    Parameters
    ----------
    faces : list
        List of faces forming the shell.
    shell_type : str, optional
        'outer' for outer shell, 'void' for inner void.
    closed : bool or None, optional
        True if shell should be closed (validated when added to TopologyGraph).
        False if shell is explicitly open.
        None to auto-detect closure when added to TopologyGraph.
    tags : dict, optional
        User metadata.

    Returns
    -------
    list
        Shell definition: ['brep_shell', face_refs, metadata_dict]

    Notes
    -----
    Closure validation is performed when the shell is added to a TopologyGraph.
    A closed shell satisfies:
    - Every edge is used by exactly two faces
    - The Euler characteristic V - E + F = 2 (for simple closed surfaces)
    """
    # Get face IDs
    face_refs = []
    for f in faces:
        if is_brep_face(f):
            face_refs.append(face_id(f))
        else:
            face_refs.append(f)  # Assume it's already an ID

    meta = {
        'id': _generate_id(),
        'shell_type': shell_type,
        'closed': closed,  # None means "to be determined"
        'solid_id': None,  # Will be set when added to a solid
        'tags': tags or {},
    }

    return ['brep_shell', face_refs, meta]


def is_brep_shell(obj):
    """Return True if obj is a BREP shell."""
    return (isinstance(obj, list) and len(obj) == 3 and obj[0] == 'brep_shell'
            and isinstance(obj[1], list) and isinstance(obj[2], dict))


def shell_id(s):
    """Return the ID of a shell."""
    if not is_brep_shell(s):
        raise ValueError("Not a BREP shell")
    return s[2].get('id')


def shell_faces(s):
    """Return the list of face IDs in a shell."""
    if not is_brep_shell(s):
        raise ValueError("Not a BREP shell")
    return deepcopy(s[1])


def shell_closed(s):
    """Return the closure status of a shell.

    Returns
    -------
    bool or None
        True if closed, False if open, None if not yet determined.
    """
    if not is_brep_shell(s):
        raise ValueError("Not a BREP shell")
    return s[2].get('closed')


def set_shell_closed(s, closed):
    """Set the closure status of a shell."""
    if not is_brep_shell(s):
        raise ValueError("Not a BREP shell")
    s[2]['closed'] = closed


# -----------------------------------------------------------------------------
# BrepSolid (Native)
# -----------------------------------------------------------------------------

def brep_solid(shells, *, tags=None):
    """Create a native BREP solid.

    Parameters
    ----------
    shells : list
        List of shells (first is outer, rest are voids).
    tags : dict, optional
        User metadata.

    Returns
    -------
    list
        Solid definition: ['brep_solid', shell_refs, metadata_dict]
    """
    # Get shell IDs
    shell_refs = []
    for s in shells:
        if is_brep_shell(s):
            shell_refs.append(shell_id(s))
        else:
            shell_refs.append(s)  # Assume it's already an ID

    meta = {
        'id': _generate_id(),
        'tags': tags or {},
    }

    return ['brep_solid', shell_refs, meta]


def is_brep_solid_native(obj):
    """Return True if obj is a native BREP solid."""
    return (isinstance(obj, list) and len(obj) == 3 and obj[0] == 'brep_solid'
            and isinstance(obj[1], list) and isinstance(obj[2], dict))


def solid_id(s):
    """Return the ID of a solid."""
    if not is_brep_solid_native(s):
        raise ValueError("Not a native BREP solid")
    return s[2].get('id')


def solid_shells(s):
    """Return the list of shell IDs in a solid."""
    if not is_brep_solid_native(s):
        raise ValueError("Not a native BREP solid")
    return deepcopy(s[1])


# -----------------------------------------------------------------------------
# Topology Container
# -----------------------------------------------------------------------------

class TopologyGraph:
    """Container for a complete BREP topology graph.

    This class manages the relationships between all topology entities
    (vertices, edges, trims, loops, faces, shells, solids) and provides
    lookup and traversal operations.
    """

    def __init__(self):
        """Initialize an empty topology graph."""
        self.vertices = {}  # id -> vertex
        self.edges = {}     # id -> edge
        self.trims = {}     # id -> trim
        self.loops = {}     # id -> loop
        self.faces = {}     # id -> face
        self.shells = {}    # id -> shell
        self.solids = {}    # id -> solid

    def add_vertex(self, vertex):
        """Add a vertex to the graph."""
        if not is_brep_vertex(vertex):
            raise ValueError("Not a BREP vertex")
        vid = vertex_id(vertex)
        self.vertices[vid] = vertex
        return vid

    def add_edge(self, edge):
        """Add an edge to the graph."""
        if not is_brep_edge(edge):
            raise ValueError("Not a BREP edge")
        eid = edge_id(edge)
        self.edges[eid] = edge

        # Register edge with incident vertices
        start_id, end_id = edge_vertices(edge)
        if start_id in self.vertices:
            self.vertices[start_id][2]['edges'].append(eid)
        if end_id in self.vertices and end_id != start_id:
            self.vertices[end_id][2]['edges'].append(eid)

        return eid

    def add_trim(self, trim):
        """Add a trim to the graph."""
        if not is_brep_trim(trim):
            raise ValueError("Not a BREP trim")
        tid = trim_id(trim)
        self.trims[tid] = trim
        return tid

    def add_loop(self, loop):
        """Add a loop to the graph."""
        if not is_brep_loop(loop):
            raise ValueError("Not a BREP loop")
        lid = loop_id(loop)
        self.loops[lid] = loop
        return lid

    def add_face(self, face):
        """Add a face to the graph."""
        if not is_brep_face(face):
            raise ValueError("Not a BREP face")
        fid = face_id(face)
        self.faces[fid] = face

        # Register face with its edges
        for lid in face_loops(face):
            if lid in self.loops:
                for tid in loop_trims(self.loops[lid]):
                    if tid in self.trims:
                        eid = trim_edge_id(self.trims[tid])
                        if eid in self.edges:
                            self.edges[eid][2]['faces'].append(fid)

        return fid

    def add_shell(self, shell, validate_closure=True):
        """Add a shell to the graph.

        Parameters
        ----------
        shell : brep_shell
            The shell to add.
        validate_closure : bool, optional
            If True (default), validate closure when shell.closed is True.
            If shell.closed is None, auto-detect and set the closure status.

        Raises
        ------
        ShellClosureError
            If shell.closed is True but the shell is not topologically closed.
        """
        if not is_brep_shell(shell):
            raise ValueError("Not a BREP shell")
        sid = shell_id(shell)
        self.shells[sid] = shell

        # Register shell with its faces
        for fid in shell_faces(shell):
            if fid in self.faces:
                self.faces[fid][2]['shell_id'] = sid

        # Handle closure validation
        if validate_closure:
            closed_flag = shell_closed(shell)

            if closed_flag is None:
                # Auto-detect closure
                is_closed, _ = self.compute_shell_closure(sid)
                set_shell_closed(shell, is_closed)

            elif closed_flag is True:
                # Validate that shell is actually closed
                self.validate_shell_closure(sid)

            # If closed_flag is False, no validation needed (explicitly open)

        return sid

    def add_solid(self, solid, validate=True, tolerance=1e-10):
        """Add a solid to the graph.

        Parameters
        ----------
        solid : brep_solid
            The solid to add.
        validate : bool, optional
            If True (default), validate the solid (shell closure, no
            intersections, proper containment hierarchy).
        tolerance : float, optional
            Geometric tolerance for validation tests.

        Raises
        ------
        SolidValidationError
            If validation is enabled and the solid is invalid.
        """
        if not is_brep_solid_native(solid):
            raise ValueError("Not a native BREP solid")
        sid = solid_id(solid)
        self.solids[sid] = solid

        # Register solid with its shells
        for shid in solid_shells(solid):
            if shid in self.shells:
                self.shells[shid][2]['solid_id'] = sid

        # Validate if requested
        if validate:
            self.validate_solid(sid, tolerance)

        return sid

    def get_vertex(self, vid):
        """Get a vertex by ID."""
        return self.vertices.get(vid)

    def get_edge(self, eid):
        """Get an edge by ID."""
        return self.edges.get(eid)

    def get_trim(self, tid):
        """Get a trim by ID."""
        return self.trims.get(tid)

    def get_loop(self, lid):
        """Get a loop by ID."""
        return self.loops.get(lid)

    def get_face(self, fid):
        """Get a face by ID."""
        return self.faces.get(fid)

    def get_shell(self, shid):
        """Get a shell by ID."""
        return self.shells.get(shid)

    def get_solid(self, sid):
        """Get a solid by ID."""
        return self.solids.get(sid)

    def vertex_edges(self, vid):
        """Return all edges incident to a vertex."""
        vertex = self.get_vertex(vid)
        if vertex is None:
            return []
        return [self.edges[eid] for eid in vertex[2].get('edges', [])
                if eid in self.edges]

    def edge_faces(self, eid):
        """Return all faces using an edge."""
        edge = self.get_edge(eid)
        if edge is None:
            return []
        return [self.faces[fid] for fid in edge[2].get('faces', [])
                if fid in self.faces]

    def face_edges(self, fid):
        """Return all edges bounding a face."""
        face = self.get_face(fid)
        if face is None:
            return []
        edges = []
        for lid in face_loops(face):
            loop = self.get_loop(lid)
            if loop:
                for tid in loop_trims(loop):
                    trim = self.get_trim(tid)
                    if trim:
                        eid = trim_edge_id(trim)
                        edge = self.get_edge(eid)
                        if edge and edge not in edges:
                            edges.append(edge)
        return edges

    def evaluate_edge(self, eid):
        """Evaluate an edge to get its yapCAD curve geometry.

        This constructs the actual curve from vertex positions and edge parameters.
        """
        edge = self.get_edge(eid)
        if edge is None:
            raise ValueError(f"Edge {eid} not found")

        start_id, end_id = edge_vertices(edge)
        start_vertex = self.get_vertex(start_id)
        end_vertex = self.get_vertex(end_id)

        if start_vertex is None or end_vertex is None:
            raise ValueError("Edge references missing vertex")

        return evaluate_edge_curve(edge, vertex_location(start_vertex),
                                   vertex_location(end_vertex))

    def transform(self, transform_point_func):
        """Transform the entire topology graph.

        Parameters
        ----------
        transform_point_func : callable
            Function that takes a point and returns the transformed point.

        This transforms:
        - All vertex locations
        - All internal geometry points in edges (arc centers, etc.)
        - Face surfaces would need separate handling
        """
        # Transform vertices
        for vertex in self.vertices.values():
            old_loc = vertex_location(vertex)
            new_loc = transform_point_func(old_loc)
            set_vertex_location(vertex, new_loc)

        # Transform internal edge geometry
        for edge in self.edges.values():
            transform_edge_params(edge, transform_point_func)

    def summary(self):
        """Return a summary of the topology graph contents."""
        return {
            'vertices': len(self.vertices),
            'edges': len(self.edges),
            'trims': len(self.trims),
            'loops': len(self.loops),
            'faces': len(self.faces),
            'shells': len(self.shells),
            'solids': len(self.solids),
        }

    # -------------------------------------------------------------------------
    # Shell Closure Validation
    # -------------------------------------------------------------------------

    def shell_edge_usage(self, shell_id):
        """Compute edge usage counts for a shell.

        Returns a dict mapping edge_id -> list of face_ids that use that edge
        within this shell.
        """
        shell = self.get_shell(shell_id)
        if shell is None:
            raise ValueError(f"Shell {shell_id} not found")

        shell_face_ids = set(shell_faces(shell))
        edge_usage = {}

        for fid in shell_face_ids:
            face = self.get_face(fid)
            if face is None:
                continue

            for lid in face_loops(face):
                loop = self.get_loop(lid)
                if loop is None:
                    continue

                for tid in loop_trims(loop):
                    trim = self.get_trim(tid)
                    if trim is None:
                        continue

                    eid = trim_edge_id(trim)
                    if eid not in edge_usage:
                        edge_usage[eid] = []
                    edge_usage[eid].append(fid)

        return edge_usage

    def shell_vertices(self, shell_id):
        """Get the set of vertex IDs used by a shell."""
        edge_usage = self.shell_edge_usage(shell_id)
        vertex_ids = set()

        for eid in edge_usage:
            edge = self.get_edge(eid)
            if edge:
                start_id, end_id = edge_vertices(edge)
                vertex_ids.add(start_id)
                vertex_ids.add(end_id)

        return vertex_ids

    def compute_shell_closure(self, shell_id):
        """Compute whether a shell is topologically closed.

        A shell is closed if:
        1. Every edge is used by exactly two faces (manifold condition)
        2. The Euler characteristic V - E + F = 2 (for simple closed surfaces)

        Parameters
        ----------
        shell_id : str
            ID of the shell to check.

        Returns
        -------
        tuple
            (is_closed: bool, details: dict)
            details contains:
            - 'edge_usage': dict mapping edge_id -> list of face_ids
            - 'boundary_edges': list of edges used by only 1 face
            - 'non_manifold_edges': list of edges used by >2 faces
            - 'euler_characteristic': computed V - E + F
            - 'vertex_count': number of vertices
            - 'edge_count': number of edges
            - 'face_count': number of faces
        """
        shell = self.get_shell(shell_id)
        if shell is None:
            raise ValueError(f"Shell {shell_id} not found")

        edge_usage = self.shell_edge_usage(shell_id)
        vertex_ids = self.shell_vertices(shell_id)
        face_ids = set(shell_faces(shell))

        # Classify edges
        boundary_edges = []      # Used by exactly 1 face (open boundary)
        non_manifold_edges = []  # Used by >2 faces (non-manifold)

        for eid, faces in edge_usage.items():
            usage_count = len(faces)
            if usage_count == 1:
                boundary_edges.append(eid)
            elif usage_count > 2:
                non_manifold_edges.append(eid)

        # Compute Euler characteristic
        V = len(vertex_ids)
        E = len(edge_usage)
        F = len(face_ids)
        euler = V - E + F

        details = {
            'edge_usage': edge_usage,
            'boundary_edges': boundary_edges,
            'non_manifold_edges': non_manifold_edges,
            'euler_characteristic': euler,
            'vertex_count': V,
            'edge_count': E,
            'face_count': F,
        }

        # Shell is closed if no boundary edges and no non-manifold edges
        # (Euler characteristic check is informational - it equals 2 for
        # a sphere, but can differ for other topologies like torus)
        is_closed = (len(boundary_edges) == 0 and len(non_manifold_edges) == 0)

        return is_closed, details

    def validate_shell_closure(self, shell_id):
        """Validate that a shell is topologically closed.

        Parameters
        ----------
        shell_id : str
            ID of the shell to validate.

        Raises
        ------
        ShellClosureError
            If the shell is not closed, with details about the failure.

        Returns
        -------
        dict
            Closure details if validation passes.
        """
        is_closed, details = self.compute_shell_closure(shell_id)

        if not is_closed:
            messages = []

            if details['boundary_edges']:
                n = len(details['boundary_edges'])
                messages.append(f"{n} boundary edge(s) (used by only 1 face)")

            if details['non_manifold_edges']:
                n = len(details['non_manifold_edges'])
                messages.append(f"{n} non-manifold edge(s) (used by >2 faces)")

            msg = f"Shell is not closed: {'; '.join(messages)}"
            raise ShellClosureError(msg, details)

        return details

    # -------------------------------------------------------------------------
    # Solid Validation - Shell Geometry
    # -------------------------------------------------------------------------

    def shell_bbox(self, shell_id):
        """Compute bounding box for a shell from its face surfaces.

        Returns
        -------
        tuple or None
            (min_point, max_point) bounding box, or None if shell is empty.
        """
        from yapcad.analytic_surfaces import is_analytic_surface, tessellate_surface

        shell = self.get_shell(shell_id)
        if shell is None:
            raise ValueError(f"Shell {shell_id} not found")

        min_pt = [float('inf'), float('inf'), float('inf')]
        max_pt = [float('-inf'), float('-inf'), float('-inf')]
        has_points = False

        for fid in shell_faces(shell):
            face = self.get_face(fid)
            if face is None:
                continue

            surface = face_surface(face)
            if surface is None:
                continue

            # Tessellate the surface to get vertex positions
            if is_analytic_surface(surface):
                mesh = tessellate_surface(surface, u_divisions=8, v_divisions=8)
                if mesh and len(mesh) > 1:
                    verts = mesh[1]
                    for v in verts:
                        has_points = True
                        min_pt[0] = min(min_pt[0], v[0])
                        min_pt[1] = min(min_pt[1], v[1])
                        min_pt[2] = min(min_pt[2], v[2])
                        max_pt[0] = max(max_pt[0], v[0])
                        max_pt[1] = max(max_pt[1], v[1])
                        max_pt[2] = max(max_pt[2], v[2])

        if not has_points:
            return None

        return (point(min_pt[0], min_pt[1], min_pt[2]),
                point(max_pt[0], max_pt[1], max_pt[2]))

    def tessellate_shell(self, shell_id, u_divisions=16, v_divisions=16):
        """Tessellate all faces in a shell into a single surface mesh.

        Returns
        -------
        list or None
            yapCAD surface ['surface', verts, normals, faces, ...] or None.
        """
        from yapcad.analytic_surfaces import (
            is_analytic_surface, tessellate_surface, surface_normal
        )

        shell = self.get_shell(shell_id)
        if shell is None:
            raise ValueError(f"Shell {shell_id} not found")

        all_verts = []
        all_normals = []
        all_faces = []

        for fid in shell_faces(shell):
            face = self.get_face(fid)
            if face is None:
                continue

            surface = face_surface(face)
            if surface is None or not is_analytic_surface(surface):
                continue

            mesh = tessellate_surface(surface, u_divisions=u_divisions,
                                      v_divisions=v_divisions)
            if not mesh or len(mesh) < 4:
                continue

            verts = mesh[1]
            normals = mesh[2]
            faces = mesh[3]

            # Offset face indices by current vertex count
            offset = len(all_verts)
            for f in faces:
                all_faces.append([f[0] + offset, f[1] + offset, f[2] + offset])

            all_verts.extend(verts)
            all_normals.extend(normals)

        if not all_faces:
            return None

        return ['surface', all_verts, all_normals, all_faces, [], {}]

    def shell_octree(self, shell_id, rebuild=False):
        """Get or build an octree for a shell's tessellated triangles.

        The octree is cached in shell metadata for performance.

        Returns
        -------
        NTree or None
            Octree of shell triangles, or None if shell is empty.
        """
        from yapcad.octtree import NTree

        shell = self.get_shell(shell_id)
        if shell is None:
            raise ValueError(f"Shell {shell_id} not found")

        # Check cache
        meta = shell[2]
        if not rebuild and '_octree' in meta and not meta.get('_octree_dirty'):
            return meta.get('_octree')

        # Build tessellation
        mesh = self.tessellate_shell(shell_id)
        if mesh is None or not mesh[3]:
            meta['_octree'] = None
            meta['_octree_dirty'] = False
            return None

        verts = mesh[1]
        faces = mesh[3]

        # Compute extent for mindim
        extent = 0.0
        for axis in range(3):
            vals = [v[axis] for v in verts]
            if vals:
                extent = max(extent, max(vals) - min(vals))
        mindim = max(extent / 32.0, epsilon)

        # Build octree
        tree = NTree(n=8, mindim=mindim)
        for face in faces:
            if len(face) >= 3:
                tri = [verts[face[0]], verts[face[1]], verts[face[2]]]
                tree.addElement(tri)
        tree.updateTree()

        meta['_octree'] = tree
        meta['_octree_dirty'] = False
        return tree

    def invalidate_shell_octree(self, shell_id):
        """Mark a shell's octree as needing rebuild."""
        shell = self.get_shell(shell_id)
        if shell is not None:
            shell[2]['_octree_dirty'] = True

    # -------------------------------------------------------------------------
    # Solid Validation - Shell Intersection Testing
    # -------------------------------------------------------------------------

    def shells_bboxes_overlap(self, shell_id1, shell_id2, tolerance=0.0):
        """Quick check if two shells' bounding boxes overlap.

        Returns False if boxes don't overlap (shells definitely don't intersect).
        Returns True if boxes overlap (shells might intersect).
        """
        from yapcad.octtree import boxoverlap2

        bbox1 = self.shell_bbox(shell_id1)
        bbox2 = self.shell_bbox(shell_id2)

        if bbox1 is None or bbox2 is None:
            return False

        # Expand boxes by tolerance
        if tolerance > 0:
            bbox1 = (
                point(bbox1[0][0] - tolerance, bbox1[0][1] - tolerance,
                      bbox1[0][2] - tolerance),
                point(bbox1[1][0] + tolerance, bbox1[1][1] + tolerance,
                      bbox1[1][2] + tolerance)
            )
            bbox2 = (
                point(bbox2[0][0] - tolerance, bbox2[0][1] - tolerance,
                      bbox2[0][2] - tolerance),
                point(bbox2[1][0] + tolerance, bbox2[1][1] + tolerance,
                      bbox2[1][2] + tolerance)
            )

        return boxoverlap2(list(bbox1), list(bbox2), dim3=True)

    def _iter_shell_triangles(self, shell_id):
        """Iterate over triangles in a shell's tessellation."""
        mesh = self.tessellate_shell(shell_id)
        if mesh is None:
            return
        verts = mesh[1]
        faces = mesh[3]
        for face in faces:
            if len(face) >= 3:
                yield [verts[face[0]], verts[face[1]], verts[face[2]]]

    def shells_intersect(self, shell_id1, shell_id2, tolerance=1e-10):
        """Test if two shells intersect using octree acceleration.

        Returns True if shells share interior volume or surface intersection.
        """
        from yapcad.octtree import boxoverlap2
        from yapcad.geom3d import triTriIntersect

        # Quick bbox rejection
        if not self.shells_bboxes_overlap(shell_id1, shell_id2, tolerance):
            return False

        # Get octrees
        tree1 = self.shell_octree(shell_id1)
        tree2 = self.shell_octree(shell_id2)

        if tree1 is None or tree2 is None:
            return False

        # For each triangle in shell1, query shell2's octree for candidates
        for tri1 in self._iter_shell_triangles(shell_id1):
            # Compute bbox for triangle
            tri_bbox = [
                point(min(tri1[0][0], tri1[1][0], tri1[2][0]) - tolerance,
                      min(tri1[0][1], tri1[1][1], tri1[2][1]) - tolerance,
                      min(tri1[0][2], tri1[1][2], tri1[2][2]) - tolerance),
                point(max(tri1[0][0], tri1[1][0], tri1[2][0]) + tolerance,
                      max(tri1[0][1], tri1[1][1], tri1[2][1]) + tolerance,
                      max(tri1[0][2], tri1[1][2], tri1[2][2]) + tolerance)
            ]

            # Query octree for candidate triangles
            candidates = tree2.getElements(tri_bbox)
            for elem in candidates:
                if isinstance(elem, tuple):
                    tri2 = elem[0]
                else:
                    tri2 = elem

                if triTriIntersect(tri1, tri2):
                    return True

        return False

    # -------------------------------------------------------------------------
    # Solid Validation - Containment Testing
    # -------------------------------------------------------------------------

    def shell_contains_point(self, shell_id, p, tolerance=1e-10):
        """Test if a point is inside a closed shell using ray casting.

        Returns True if point is inside, False otherwise.
        """
        from yapcad.geom import vect, dot, cross, sub, mag

        shell = self.get_shell(shell_id)
        if shell is None:
            raise ValueError(f"Shell {shell_id} not found")

        # Quick bbox check
        bbox = self.shell_bbox(shell_id)
        if bbox is None:
            return False

        expanded_bbox = [
            point(bbox[0][0] - tolerance, bbox[0][1] - tolerance,
                  bbox[0][2] - tolerance),
            point(bbox[1][0] + tolerance, bbox[1][1] + tolerance,
                  bbox[1][2] + tolerance)
        ]

        if not (expanded_bbox[0][0] <= p[0] <= expanded_bbox[1][0] and
                expanded_bbox[0][1] <= p[1] <= expanded_bbox[1][1] and
                expanded_bbox[0][2] <= p[2] <= expanded_bbox[1][2]):
            return False

        # Ray casting in +X direction
        extent = max(bbox[1][0] - bbox[0][0],
                     bbox[1][1] - bbox[0][1],
                     bbox[1][2] - bbox[0][2])
        ray_length = extent * 2.0

        direction = vect(1, 0, 0, 0)
        far_point = point(p[0] + ray_length, p[1], p[2])

        # Query octree for triangles along ray
        tree = self.shell_octree(shell_id)
        query_bbox = [
            point(min(p[0], far_point[0]) - tolerance,
                  p[1] - tolerance, p[2] - tolerance),
            point(max(p[0], far_point[0]) + tolerance,
                  p[1] + tolerance, p[2] + tolerance)
        ]

        if tree is None:
            candidates = list(self._iter_shell_triangles(shell_id))
        else:
            elems = tree.getElements(query_bbox)
            candidates = []
            for elem in elems:
                if isinstance(elem, tuple):
                    candidates.append(elem[0])
                else:
                    candidates.append(elem)

        # Count ray-triangle intersections
        hits = 0
        for tri in candidates:
            result = self._ray_triangle_intersect(p, direction, tri, tolerance)
            if result is not None:
                t, _ = result
                if t > tolerance:  # Hit is in front of point
                    hits += 1

        # Odd number of hits = inside
        return (hits % 2) == 1

    def _ray_triangle_intersect(self, origin, direction, tri, tolerance):
        """Ray-triangle intersection using Möller–Trumbore algorithm.

        Returns (t, hit_point) or None.
        """
        from yapcad.geom import sub, cross, dot

        v0, v1, v2 = tri
        e1 = sub(v1, v0)
        e2 = sub(v2, v0)

        h = cross(direction, e2)
        a = dot(e1, h)

        if abs(a) < tolerance:
            return None  # Ray parallel to triangle

        f = 1.0 / a
        s = sub(origin, v0)
        u = f * dot(s, h)

        if u < -tolerance or u > 1.0 + tolerance:
            return None

        q = cross(s, e1)
        v = f * dot(direction, q)

        if v < -tolerance or u + v > 1.0 + tolerance:
            return None

        t = f * dot(e2, q)

        if t < -tolerance:
            return None

        hit = point(origin[0] + direction[0] * t,
                    origin[1] + direction[1] * t,
                    origin[2] + direction[2] * t)
        return (t, hit)

    # -------------------------------------------------------------------------
    # Solid Validation - Full Validation
    # -------------------------------------------------------------------------

    def validate_solid(self, solid_id, tolerance=1e-10):
        """Validate a BREP solid.

        Validates:
        1. All shells are closed
        2. Shells do not intersect each other
        3. Void shells are contained within the outer shell

        Parameters
        ----------
        solid_id : str
            ID of the solid to validate.
        tolerance : float
            Geometric tolerance for intersection tests.

        Raises
        ------
        SolidValidationError
            If validation fails, with details about the failure.

        Returns
        -------
        dict
            Validation details if successful.
        """
        solid = self.get_solid(solid_id)
        if solid is None:
            raise ValueError(f"Solid {solid_id} not found")

        shell_ids = solid_shells(solid)
        if not shell_ids:
            raise SolidValidationError("Solid has no shells",
                                       {'solid_id': solid_id})

        details = {
            'solid_id': solid_id,
            'shell_count': len(shell_ids),
            'shell_closure': {},
            'shell_intersections': [],
            'containment_violations': [],
        }

        # 1. Validate all shells are closed
        for sid in shell_ids:
            shell = self.get_shell(sid)
            if shell is None:
                raise SolidValidationError(f"Shell {sid} not found in graph",
                                          details)

            closed_flag = shell_closed(shell)
            if closed_flag is False:
                raise SolidValidationError(
                    f"Solid contains explicitly open shell {sid}",
                    details)

            if closed_flag is None:
                # Auto-detect
                is_closed, closure_details = self.compute_shell_closure(sid)
                details['shell_closure'][sid] = closure_details
                if not is_closed:
                    raise SolidValidationError(
                        f"Shell {sid} is not closed: "
                        f"{len(closure_details['boundary_edges'])} boundary edges",
                        details)
            else:
                # Verify claimed closure
                is_closed, closure_details = self.compute_shell_closure(sid)
                details['shell_closure'][sid] = closure_details
                if not is_closed:
                    raise SolidValidationError(
                        f"Shell {sid} claims to be closed but is not",
                        details)

        # 2. Check for shell intersections
        if len(shell_ids) > 1:
            for i in range(len(shell_ids)):
                for j in range(i + 1, len(shell_ids)):
                    sid1, sid2 = shell_ids[i], shell_ids[j]
                    if self.shells_intersect(sid1, sid2, tolerance):
                        details['shell_intersections'].append((sid1, sid2))
                        raise SolidValidationError(
                            f"Shells {sid1} and {sid2} intersect",
                            details)

        # 3. Validate containment hierarchy (void shells inside outer)
        if len(shell_ids) > 1:
            # First shell is outer, rest are voids
            outer_shell_id = shell_ids[0]

            for void_shell_id in shell_ids[1:]:
                # Get a point on the void shell (center of bbox)
                bbox = self.shell_bbox(void_shell_id)
                if bbox is None:
                    continue

                center = point(
                    (bbox[0][0] + bbox[1][0]) / 2,
                    (bbox[0][1] + bbox[1][1]) / 2,
                    (bbox[0][2] + bbox[1][2]) / 2
                )

                # Void shell center should be inside outer shell
                if not self.shell_contains_point(outer_shell_id, center,
                                                  tolerance):
                    details['containment_violations'].append({
                        'void_shell': void_shell_id,
                        'outer_shell': outer_shell_id,
                        'test_point': center,
                    })
                    raise SolidValidationError(
                        f"Void shell {void_shell_id} is not contained "
                        f"within outer shell {outer_shell_id}",
                        details)

        return details


# -----------------------------------------------------------------------------
# Serialization and Deserialization
# -----------------------------------------------------------------------------

def _serialize_point(p):
    """Convert a point to a JSON-serializable list."""
    if p is None:
        return None
    return [float(p[0]), float(p[1]), float(p[2]), float(p[3]) if len(p) > 3 else 1.0]


def _deserialize_point(data):
    """Reconstruct a point from serialized data."""
    if data is None:
        return None
    return point(data[0], data[1], data[2])


def _serialize_vertex(v):
    """Serialize a BREP vertex to a dict."""
    if not is_brep_vertex(v):
        raise ValueError("Not a BREP vertex")
    return {
        'type': 'brep_vertex',
        'location': _serialize_point(v[1]),
        'id': v[2].get('id'),
        'tolerance': v[2].get('tolerance', epsilon),
        'tags': v[2].get('tags', {}),
    }


def _deserialize_vertex(data):
    """Reconstruct a BREP vertex from serialized data."""
    loc = _deserialize_point(data['location'])
    v = brep_vertex(loc, tolerance=data.get('tolerance', epsilon),
                    tags=data.get('tags', {}))
    # Restore original ID
    v[2]['id'] = data['id']
    return v


def _serialize_edge(e):
    """Serialize a BREP edge to a dict."""
    if not is_brep_edge(e):
        raise ValueError("Not a BREP edge")
    meta = e[2]
    result = {
        'type': 'brep_edge',
        'curve_type': e[1],
        'id': meta.get('id'),
        'start_vertex': meta.get('start_vertex_id'),  # Note: key is start_vertex_id in metadata
        'end_vertex': meta.get('end_vertex_id'),      # Note: key is end_vertex_id in metadata
        'tags': meta.get('tags', {}),
    }
    # Serialize curve parameters
    if 'curve_params' in meta:
        params = meta['curve_params']
        serialized_params = {}
        for key, value in params.items():
            if key in ('center', 'axis', 'control_points'):
                if isinstance(value, list) and value and isinstance(value[0], (list, tuple)):
                    # List of points
                    serialized_params[key] = [_serialize_point(p) for p in value]
                else:
                    serialized_params[key] = _serialize_point(value)
            elif key == 'knots':
                serialized_params[key] = list(value) if value else []
            elif key == 'weights':
                serialized_params[key] = list(value) if value else None
            else:
                serialized_params[key] = value
        result['curve_params'] = serialized_params
    return result


def _deserialize_edge(data, vertex_map=None):
    """Reconstruct a BREP edge from serialized data."""
    curve_type = data['curve_type']
    start_id = data['start_vertex']
    end_id = data['end_vertex']

    # Reconstruct curve parameters
    curve_params = {}
    if 'curve_params' in data:
        params = data['curve_params']
        for key, value in params.items():
            if key in ('center', 'axis'):
                curve_params[key] = _deserialize_point(value) if value else None
            elif key == 'control_points':
                if value and isinstance(value[0], list):
                    curve_params[key] = [_deserialize_point(p) for p in value]
                else:
                    curve_params[key] = _deserialize_point(value) if value else None
            elif key == 'knots':
                curve_params[key] = list(value) if value else []
            elif key == 'weights':
                curve_params[key] = list(value) if value else None
            else:
                curve_params[key] = value

    # Create edge with appropriate type
    start_v = vertex_map.get(start_id) if vertex_map else start_id
    end_v = vertex_map.get(end_id) if vertex_map else end_id

    if curve_type == 'line':
        e = line_edge(start_v, end_v, tags=data.get('tags', {}))
    elif curve_type == 'arc':
        # arc_edge uses bulge OR (center, orientation), not radius
        if curve_params.get('bulge') is not None:
            e = arc_edge(start_v, end_v, bulge=curve_params.get('bulge'),
                         tags=data.get('tags', {}))
        else:
            e = arc_edge(start_v, end_v, center=curve_params.get('center'),
                         orientation=curve_params.get('orientation', 1),
                         tags=data.get('tags', {}))
    elif curve_type == 'circle':
        e = circle_edge(start_v, end_v, curve_params.get('center'),
                        curve_params.get('axis'), curve_params.get('radius'),
                        t_start=curve_params.get('t_start', 0.0),
                        t_end=curve_params.get('t_end'),
                        tags=data.get('tags', {}))
    elif curve_type == 'bspline':
        e = bspline_edge(start_v, end_v, curve_params.get('control_points'),
                         curve_params.get('knots'), curve_params.get('degree'),
                         weights=curve_params.get('weights'),
                         t_start=curve_params.get('t_start'),
                         t_end=curve_params.get('t_end'),
                         tags=data.get('tags', {}))
    else:
        # Generic edge
        e = brep_edge(curve_type, start_v, end_v,
                      curve_params=curve_params, tags=data.get('tags', {}))

    # Restore original ID
    e[2]['id'] = data['id']
    return e


def _serialize_trim(t):
    """Serialize a BREP trim to a dict."""
    if not is_brep_trim(t):
        raise ValueError("Not a BREP trim")
    meta = t[2]
    return {
        'type': 'brep_trim',
        'edge_id': t[1],
        'id': meta.get('id'),
        'sense': meta.get('sense', True),
        'uv_curve': meta.get('uv_curve'),
        'face_id': meta.get('face_id'),
        'loop_id': meta.get('loop_id'),
        'tags': meta.get('tags', {}),
    }


def _deserialize_trim(data, edge_map=None):
    """Reconstruct a BREP trim from serialized data."""
    edge_ref = edge_map.get(data['edge_id']) if edge_map else data['edge_id']

    t = brep_trim(edge_ref, sense=data.get('sense', True),
                  uv_curve=data.get('uv_curve'),
                  tags=data.get('tags', {}))
    t[2]['id'] = data['id']
    if data.get('face_id'):
        t[2]['face_id'] = data['face_id']
    if data.get('loop_id'):
        t[2]['loop_id'] = data['loop_id']
    return t


def _serialize_loop(l):
    """Serialize a BREP loop to a dict."""
    if not is_brep_loop(l):
        raise ValueError("Not a BREP loop")
    meta = l[2]
    return {
        'type': 'brep_loop',
        'trim_ids': list(l[1]),
        'id': meta.get('id'),
        'loop_type': meta.get('loop_type', 'outer'),
        'face_id': meta.get('face_id'),
        'tags': meta.get('tags', {}),
    }


def _deserialize_loop(data, trim_map=None):
    """Reconstruct a BREP loop from serialized data."""
    trim_refs = [trim_map.get(tid) if trim_map else tid for tid in data['trim_ids']]

    l = brep_loop(trim_refs, loop_type=data.get('loop_type', 'outer'),
                  tags=data.get('tags', {}))
    l[2]['id'] = data['id']
    if data.get('face_id'):
        l[2]['face_id'] = data['face_id']
    return l


def _serialize_surface(surf):
    """Serialize an analytic surface to a dict.

    Analytic surfaces are already list-based structures that are mostly
    JSON-serializable, but we need to handle points properly.
    """
    if surf is None:
        return None

    surf_type = surf[0] if isinstance(surf, list) and surf else None
    if surf_type is None:
        return None

    # The surface is stored as ['type', data, metadata]
    # We need to serialize it properly
    result = {
        'surface_type': surf_type,
    }

    if surf_type == 'plane_surface':
        result['origin'] = _serialize_point(surf[1])
        result['metadata'] = _serialize_surface_metadata(surf[2])

    elif surf_type == 'sphere_surface':
        result['center'] = _serialize_point(surf[1])
        result['metadata'] = _serialize_surface_metadata(surf[2])

    elif surf_type == 'cylinder_surface':
        result['origin'] = _serialize_point(surf[1])
        result['metadata'] = _serialize_surface_metadata(surf[2])

    elif surf_type == 'cone_surface':
        result['apex'] = _serialize_point(surf[1])
        result['metadata'] = _serialize_surface_metadata(surf[2])

    elif surf_type == 'torus_surface':
        result['center'] = _serialize_point(surf[1])
        result['metadata'] = _serialize_surface_metadata(surf[2])

    elif surf_type == 'bspline_surface':
        # Control points are nested list of points
        control_points = surf[1]
        serialized_cpts = []
        for row in control_points:
            serialized_row = [_serialize_point(p) for p in row]
            serialized_cpts.append(serialized_row)
        result['control_points'] = serialized_cpts
        result['metadata'] = _serialize_surface_metadata(surf[2])

    elif surf_type == 'tessellated_surface':
        # Vertices and normals are lists of points
        result['vertices'] = [_serialize_point(v) for v in surf[1].get('vertices', [])]
        result['normals'] = [_serialize_point(n) for n in surf[1].get('normals', [])]
        result['faces'] = surf[1].get('faces', [])
        result['metadata'] = _serialize_surface_metadata(surf[2]) if len(surf) > 2 else {}

    else:
        # Unknown surface type - try generic serialization
        result['data'] = surf[1] if len(surf) > 1 else None
        result['metadata'] = surf[2] if len(surf) > 2 else {}

    return result


def _serialize_surface_metadata(meta):
    """Serialize surface metadata, handling vectors/points."""
    if meta is None:
        return {}
    result = {}
    for key, value in meta.items():
        if key in ('normal', 'axis', 'u_axis', 'v_axis'):
            result[key] = list(value) if value else None
        elif key == 'u_range' or key == 'v_range':
            result[key] = list(value) if value else None
        elif key == 'u_knots' or key == 'v_knots':
            result[key] = list(value) if value else []
        elif key == 'weights':
            if value and isinstance(value[0], list):
                result[key] = [list(row) for row in value]
            else:
                result[key] = list(value) if value else None
        else:
            result[key] = value
    return result


def _deserialize_surface(data):
    """Reconstruct an analytic surface from serialized data."""
    from yapcad.analytic_surfaces import (
        plane_surface, sphere_surface, cylinder_surface,
        cone_surface, torus_surface, bspline_surface, tessellated_surface
    )

    if data is None:
        return None

    surf_type = data.get('surface_type')
    if surf_type is None:
        return None

    meta = data.get('metadata', {})

    if surf_type == 'plane_surface':
        origin = _deserialize_point(data['origin'])
        normal = meta.get('normal', [0, 0, 1, 0])
        return plane_surface(origin, normal,
                            u_range=tuple(meta.get('u_range', (-1, 1))),
                            v_range=tuple(meta.get('v_range', (-1, 1))))

    elif surf_type == 'sphere_surface':
        center = _deserialize_point(data['center'])
        radius = meta.get('radius', 1.0)
        return sphere_surface(center, radius,
                             u_range=tuple(meta.get('u_range', (0, 2*pi))),
                             v_range=tuple(meta.get('v_range', (-pi/2, pi/2))))

    elif surf_type == 'cylinder_surface':
        origin = _deserialize_point(data['origin'])
        axis = meta.get('axis', [0, 0, 1, 0])
        radius = meta.get('radius', 1.0)
        return cylinder_surface(origin, axis, radius,
                               u_range=tuple(meta.get('u_range', (0, 2*pi))),
                               v_range=tuple(meta.get('v_range', (0, 1))))

    elif surf_type == 'cone_surface':
        apex = _deserialize_point(data['apex'])
        axis = meta.get('axis', [0, 0, 1, 0])
        half_angle = meta.get('half_angle', pi/4)
        return cone_surface(apex, axis, half_angle,
                           u_range=tuple(meta.get('u_range', (0, 2*pi))),
                           v_range=tuple(meta.get('v_range', (0, 1))))

    elif surf_type == 'torus_surface':
        center = _deserialize_point(data['center'])
        axis = meta.get('axis', [0, 0, 1, 0])
        major_radius = meta.get('major_radius', 2.0)
        minor_radius = meta.get('minor_radius', 0.5)
        return torus_surface(center, axis, major_radius, minor_radius,
                            u_range=tuple(meta.get('u_range', (0, 2*pi))),
                            v_range=tuple(meta.get('v_range', (0, 2*pi))))

    elif surf_type == 'bspline_surface':
        control_points = [[_deserialize_point(p) for p in row]
                         for row in data['control_points']]
        u_knots = meta.get('u_knots', [])
        v_knots = meta.get('v_knots', [])
        u_degree = meta.get('u_degree', 3)
        v_degree = meta.get('v_degree', 3)
        weights = meta.get('weights')
        return bspline_surface(control_points, u_knots, v_knots,
                              u_degree, v_degree, weights=weights,
                              u_range=tuple(meta.get('u_range', (0, 1))),
                              v_range=tuple(meta.get('v_range', (0, 1))))

    elif surf_type == 'tessellated_surface':
        vertices = [_deserialize_point(v) for v in data.get('vertices', [])]
        normals = [_deserialize_point(n) for n in data.get('normals', [])]
        faces = data.get('faces', [])
        return tessellated_surface(vertices, normals, faces,
                                  u_range=tuple(meta.get('u_range', (0, 1))),
                                  v_range=tuple(meta.get('v_range', (0, 1))))

    else:
        # Unknown type - return generic structure
        return [surf_type, data.get('data'), data.get('metadata', {})]


def _serialize_face(f):
    """Serialize a BREP face to a dict."""
    if not is_brep_face(f):
        raise ValueError("Not a BREP face")
    meta = f[2]
    return {
        'type': 'brep_face',
        'surface': _serialize_surface(f[1]),
        'id': meta.get('id'),
        'loop_ids': list(meta.get('loop_ids', [])),
        'sense': meta.get('sense', True),
        'shell_id': meta.get('shell_id'),
        'tags': meta.get('tags', {}),
    }


def _deserialize_face(data, loop_map=None):
    """Reconstruct a BREP face from serialized data."""
    surface = _deserialize_surface(data['surface'])
    loop_refs = [loop_map.get(lid) if loop_map else lid
                 for lid in data.get('loop_ids', [])]

    f = brep_face(surface, loop_refs, sense=data.get('sense', True),
                  tags=data.get('tags', {}))
    f[2]['id'] = data['id']
    if data.get('shell_id'):
        f[2]['shell_id'] = data['shell_id']
    return f


def _serialize_shell(s):
    """Serialize a BREP shell to a dict."""
    if not is_brep_shell(s):
        raise ValueError("Not a BREP shell")
    meta = s[2]
    return {
        'type': 'brep_shell',
        'face_ids': list(s[1]),
        'id': meta.get('id'),
        'closed': meta.get('closed'),
        'solid_id': meta.get('solid_id'),
        'tags': meta.get('tags', {}),
    }


def _deserialize_shell(data, face_map=None):
    """Reconstruct a BREP shell from serialized data."""
    face_refs = [face_map.get(fid) if face_map else fid
                 for fid in data.get('face_ids', [])]

    s = brep_shell(face_refs, closed=data.get('closed'),
                   tags=data.get('tags', {}))
    s[2]['id'] = data['id']
    if data.get('solid_id'):
        s[2]['solid_id'] = data['solid_id']
    return s


def _serialize_solid(s):
    """Serialize a native BREP solid to a dict."""
    if not is_brep_solid_native(s):
        raise ValueError("Not a native BREP solid")
    meta = s[2]
    return {
        'type': 'brep_solid',
        'shell_ids': list(s[1]),
        'id': meta.get('id'),
        'tags': meta.get('tags', {}),
    }


def _deserialize_solid(data, shell_map=None):
    """Reconstruct a native BREP solid from serialized data."""
    shell_refs = [shell_map.get(sid) if shell_map else sid
                  for sid in data.get('shell_ids', [])]

    s = brep_solid(shell_refs, tags=data.get('tags', {}))
    s[2]['id'] = data['id']
    return s


def serialize_topology_graph(graph):
    """Serialize a TopologyGraph to a JSON-serializable dict.

    Parameters
    ----------
    graph : TopologyGraph
        The topology graph to serialize.

    Returns
    -------
    dict
        JSON-serializable dictionary containing all topology entities.
    """
    return {
        'version': '1.0',
        'vertices': {vid: _serialize_vertex(v) for vid, v in graph.vertices.items()},
        'edges': {eid: _serialize_edge(e) for eid, e in graph.edges.items()},
        'trims': {tid: _serialize_trim(t) for tid, t in graph.trims.items()},
        'loops': {lid: _serialize_loop(l) for lid, l in graph.loops.items()},
        'faces': {fid: _serialize_face(f) for fid, f in graph.faces.items()},
        'shells': {sid: _serialize_shell(s) for sid, s in graph.shells.items()},
        'solids': {solid_id: _serialize_solid(s) for solid_id, s in graph.solids.items()},
    }


def deserialize_topology_graph(data):
    """Reconstruct a TopologyGraph from serialized data.

    Parameters
    ----------
    data : dict
        Serialized topology graph data.

    Returns
    -------
    TopologyGraph
        Reconstructed topology graph.
    """
    graph = TopologyGraph()

    # Deserialize vertices first (no dependencies)
    vertex_map = {}  # old_id -> vertex
    for vid, vdata in data.get('vertices', {}).items():
        vertex = _deserialize_vertex(vdata)
        vertex_map[vid] = vertex
        graph.vertices[vertex[2]['id']] = vertex

    # Deserialize edges (depend on vertices)
    edge_map = {}
    for eid, edata in data.get('edges', {}).items():
        edge = _deserialize_edge(edata, vertex_map=None)  # Use IDs directly
        edge_map[eid] = edge
        graph.edges[edge[2]['id']] = edge

    # Deserialize trims (depend on edges)
    trim_map = {}
    for tid, tdata in data.get('trims', {}).items():
        trim = _deserialize_trim(tdata)
        trim_map[tid] = trim
        graph.trims[trim[2]['id']] = trim

    # Deserialize loops (depend on trims)
    loop_map = {}
    for lid, ldata in data.get('loops', {}).items():
        loop = _deserialize_loop(ldata)
        loop_map[lid] = loop
        graph.loops[loop[2]['id']] = loop

    # Deserialize faces (depend on loops)
    face_map = {}
    for fid, fdata in data.get('faces', {}).items():
        face = _deserialize_face(fdata)
        face_map[fid] = face
        graph.faces[face[2]['id']] = face

    # Deserialize shells (depend on faces)
    shell_map = {}
    for sid, sdata in data.get('shells', {}).items():
        shell = _deserialize_shell(sdata)
        shell_map[sid] = shell
        graph.shells[shell[2]['id']] = shell

    # Deserialize solids (depend on shells)
    for solid_id, sdata in data.get('solids', {}).items():
        solid = _deserialize_solid(sdata)
        graph.solids[solid[2]['id']] = solid

    return graph


# -----------------------------------------------------------------------------
# Solid Metadata Integration
# -----------------------------------------------------------------------------

def attach_native_brep_to_solid(yapcad_solid, native_brep_graph):
    """Attach a native BREP topology graph to a yapCAD solid's metadata.

    Parameters
    ----------
    yapcad_solid : list
        A yapCAD solid (from geom3d).
    native_brep_graph : TopologyGraph
        The native BREP topology graph to attach.

    Notes
    -----
    This stores the serialized native BREP in the solid's metadata under
    the 'native_brep' key. The format mirrors the OCC BREP storage pattern
    used in brep.py.
    """
    from yapcad.metadata import get_solid_metadata, ensure_solid_id

    ensure_solid_id(yapcad_solid)
    meta = get_solid_metadata(yapcad_solid, create=True)

    serialized = serialize_topology_graph(native_brep_graph)
    meta['native_brep'] = {
        'encoding': 'topology-graph-json',
        'version': serialized.get('version', '1.0'),
        'data': serialized,
    }


def native_brep_from_solid(yapcad_solid):
    """Retrieve native BREP topology graph from a yapCAD solid's metadata.

    Parameters
    ----------
    yapcad_solid : list
        A yapCAD solid (from geom3d).

    Returns
    -------
    TopologyGraph or None
        The deserialized topology graph, or None if not present.
    """
    from yapcad.metadata import get_solid_metadata

    meta = get_solid_metadata(yapcad_solid, create=False)
    if not meta:
        return None

    native_brep_info = meta.get('native_brep')
    if not native_brep_info:
        return None

    data = native_brep_info.get('data')
    if not data:
        return None

    return deserialize_topology_graph(data)


def has_native_brep(yapcad_solid):
    """Check if a yapCAD solid has native BREP data attached.

    Parameters
    ----------
    yapcad_solid : list
        A yapCAD solid (from geom3d).

    Returns
    -------
    bool
        True if native BREP data is present.
    """
    from yapcad.metadata import get_solid_metadata

    meta = get_solid_metadata(yapcad_solid, create=False)
    return bool(meta and meta.get('native_brep'))


def clear_native_brep(yapcad_solid):
    """Remove native BREP data from a yapCAD solid's metadata.

    Parameters
    ----------
    yapcad_solid : list
        A yapCAD solid (from geom3d).

    Notes
    -----
    This is useful when the native BREP becomes stale (e.g., after
    modifications that don't update the BREP).
    """
    from yapcad.metadata import get_solid_metadata

    meta = get_solid_metadata(yapcad_solid, create=False)
    if meta and 'native_brep' in meta:
        del meta['native_brep']


# -----------------------------------------------------------------------------
# Native BREP Transformations
# -----------------------------------------------------------------------------

def _transform_point(p, transform_fn):
    """Apply a transformation function to a point."""
    if p is None:
        return None
    result = transform_fn(p)
    return point(result[0], result[1], result[2])


def _transform_vector(v, transform_fn):
    """Apply a transformation function to a vector (direction only)."""
    if v is None:
        return None
    # For vectors, transform from origin to get direction only
    origin = point(0, 0, 0)
    end = point(v[0], v[1], v[2])
    t_end = transform_fn(end)
    # Don't translate vectors, just rotate/scale them
    # For pure rotation/scale, we can compute by transforming the vector as a point
    # and subtracting transformed origin
    t_origin = transform_fn(origin)
    result = [t_end[0] - t_origin[0], t_end[1] - t_origin[1], t_end[2] - t_origin[2], 0.0]
    # Normalize
    mag = sqrt(result[0]**2 + result[1]**2 + result[2]**2)
    if mag > epsilon:
        result = [result[0]/mag, result[1]/mag, result[2]/mag, 0.0]
    return result


def transform_topology_graph(graph, transform_fn):
    """Apply a transformation to all geometry in a TopologyGraph.

    Parameters
    ----------
    graph : TopologyGraph
        The topology graph to transform (modified in-place).
    transform_fn : callable
        A function that takes a point and returns a transformed point.

    Notes
    -----
    This transforms:
    - All vertex locations
    - Edge curve parameters (centers, control points, axes)
    - Face surface parameters (origins, centers, normals, axes)
    - Invalidates cached octrees
    """
    # Transform all vertices
    for vid, vertex in graph.vertices.items():
        loc = vertex[1]
        vertex[1] = _transform_point(loc, transform_fn)

    # Transform edge curve parameters
    for eid, edge in graph.edges.items():
        meta = edge[2]
        if 'curve_params' in meta:
            params = meta['curve_params']
            if 'center' in params:
                params['center'] = _transform_point(params['center'], transform_fn)
            if 'axis' in params:
                params['axis'] = _transform_vector(params['axis'], transform_fn)
            if 'control_points' in params and params['control_points']:
                cpts = params['control_points']
                if isinstance(cpts[0], (list, tuple)) and len(cpts[0]) >= 3:
                    params['control_points'] = [_transform_point(p, transform_fn) for p in cpts]

    # Transform face surfaces
    for fid, face in graph.faces.items():
        surface = face[1]
        if surface is not None and isinstance(surface, list) and len(surface) >= 3:
            _transform_surface_inplace(surface, transform_fn)

    # Invalidate all shell octrees
    for sid, shell in graph.shells.items():
        shell[2]['_octree_dirty'] = True


def _transform_surface_inplace(surface, transform_fn):
    """Transform a surface's parameters in-place."""
    surf_type = surface[0]
    meta = surface[2] if len(surface) > 2 else {}

    if surf_type == 'plane_surface':
        surface[1] = _transform_point(surface[1], transform_fn)
        if 'normal' in meta:
            meta['normal'] = _transform_vector(meta['normal'], transform_fn)
        if 'u_axis' in meta:
            meta['u_axis'] = _transform_vector(meta['u_axis'], transform_fn)
        if 'v_axis' in meta:
            meta['v_axis'] = _transform_vector(meta['v_axis'], transform_fn)

    elif surf_type == 'sphere_surface':
        surface[1] = _transform_point(surface[1], transform_fn)

    elif surf_type == 'cylinder_surface':
        surface[1] = _transform_point(surface[1], transform_fn)
        if 'axis' in meta:
            meta['axis'] = _transform_vector(meta['axis'], transform_fn)

    elif surf_type == 'cone_surface':
        surface[1] = _transform_point(surface[1], transform_fn)
        if 'axis' in meta:
            meta['axis'] = _transform_vector(meta['axis'], transform_fn)

    elif surf_type == 'torus_surface':
        surface[1] = _transform_point(surface[1], transform_fn)
        if 'axis' in meta:
            meta['axis'] = _transform_vector(meta['axis'], transform_fn)

    elif surf_type == 'bspline_surface':
        # Control points are nested list
        control_points = surface[1]
        for i, row in enumerate(control_points):
            for j, pt in enumerate(row):
                control_points[i][j] = _transform_point(pt, transform_fn)

    elif surf_type == 'tessellated_surface':
        # Vertices and normals in dict at [1]
        data = surface[1]
        if 'vertices' in data:
            data['vertices'] = [_transform_point(v, transform_fn) for v in data['vertices']]
        if 'normals' in data:
            data['normals'] = [_transform_vector(n, transform_fn) for n in data['normals']]


def translate_native_brep(yapcad_solid, delta):
    """Apply a translation to native BREP data in a solid.

    Parameters
    ----------
    yapcad_solid : list
        A yapCAD solid with native BREP data.
    delta : point/list
        Translation vector [dx, dy, dz].
    """
    graph = native_brep_from_solid(yapcad_solid)
    if graph is None:
        return

    dx = float(delta[0]) if len(delta) > 0 else 0.0
    dy = float(delta[1]) if len(delta) > 1 else 0.0
    dz = float(delta[2]) if len(delta) > 2 else 0.0

    def translate_point(p):
        return point(p[0] + dx, p[1] + dy, p[2] + dz)

    transform_topology_graph(graph, translate_point)
    attach_native_brep_to_solid(yapcad_solid, graph)


def rotate_native_brep(yapcad_solid, ang, center, axis):
    """Apply a rotation to native BREP data in a solid.

    Parameters
    ----------
    yapcad_solid : list
        A yapCAD solid with native BREP data.
    ang : float
        Rotation angle in degrees.
    center : point
        Center of rotation.
    axis : point/vector
        Axis of rotation.
    """
    from yapcad.geom import rotate as geom_rotate

    graph = native_brep_from_solid(yapcad_solid)
    if graph is None:
        return

    def rotate_point(p):
        return geom_rotate(p, ang, cent=center, axis=axis)

    transform_topology_graph(graph, rotate_point)
    attach_native_brep_to_solid(yapcad_solid, graph)


def mirror_native_brep(yapcad_solid, plane):
    """Apply a mirror transformation to native BREP data in a solid.

    Parameters
    ----------
    yapcad_solid : list
        A yapCAD solid with native BREP data.
    plane : str
        Mirror plane: 'xy', 'xz', or 'yz'.
    """
    from yapcad.geom import mirror as geom_mirror

    graph = native_brep_from_solid(yapcad_solid)
    if graph is None:
        return

    def mirror_point(p):
        return geom_mirror(p, plane)

    transform_topology_graph(graph, mirror_point)
    attach_native_brep_to_solid(yapcad_solid, graph)


def scale_native_brep(yapcad_solid, factor, center=None):
    """Apply a uniform scale to native BREP data in a solid.

    Parameters
    ----------
    yapcad_solid : list
        A yapCAD solid with native BREP data.
    factor : float
        Scale factor.
    center : point, optional
        Center of scaling. Defaults to origin.
    """
    graph = native_brep_from_solid(yapcad_solid)
    if graph is None:
        return

    cx = float(center[0]) if center and len(center) > 0 else 0.0
    cy = float(center[1]) if center and len(center) > 1 else 0.0
    cz = float(center[2]) if center and len(center) > 2 else 0.0
    f = float(factor)

    def scale_point(p):
        # Scale relative to center
        return point(
            cx + (p[0] - cx) * f,
            cy + (p[1] - cy) * f,
            cz + (p[2] - cz) * f
        )

    transform_topology_graph(graph, scale_point)

    # Also update radii in surfaces
    for fid, face in graph.faces.items():
        surface = face[1]
        if surface is not None and isinstance(surface, list) and len(surface) >= 3:
            meta = surface[2] if len(surface) > 2 else {}
            if 'radius' in meta:
                meta['radius'] = meta['radius'] * f
            if 'major_radius' in meta:
                meta['major_radius'] = meta['major_radius'] * f
            if 'minor_radius' in meta:
                meta['minor_radius'] = meta['minor_radius'] * f

    attach_native_brep_to_solid(yapcad_solid, graph)
