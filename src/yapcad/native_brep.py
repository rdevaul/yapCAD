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
