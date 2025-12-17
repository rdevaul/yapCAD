"""Analytic surface definitions for yapCAD.

This module provides native analytic surface representations that store
parametric definitions rather than (or in addition to) tessellated meshes.
These surfaces can be evaluated at arbitrary (u, v) parameters and tessellated
on demand.

Surface types:
- PlaneSurface: infinite plane defined by point and normal
- SphereSurface: sphere defined by center and radius
- CylinderSurface: cylinder defined by axis, radius, and height bounds
- ConeSurface: cone defined by apex, axis, half-angle, and height bounds
- TorusSurface: torus defined by center, axis, major/minor radii

Each surface stores:
- Parameter domain (u, v bounds)
- Local coordinate system (origin, axis)
- Geometric parameters (radii, angles, etc.)
- Optional tessellation cache

Copyright (c) 2025 Richard DeVaul
MIT License
"""

from copy import deepcopy
from math import pi, sin, cos, sqrt, atan2, asin

from yapcad.geom import point, vect, ispoint, epsilon


# -----------------------------------------------------------------------------
# Plane Surface
# -----------------------------------------------------------------------------

def plane_surface(origin, normal, *, u_range=(-1.0, 1.0), v_range=(-1.0, 1.0)):
    """Create an analytic plane surface.

    Parameters
    ----------
    origin : point
        A point on the plane.
    normal : vector
        Normal vector to the plane (will be normalized).
    u_range : tuple, optional
        Parameter range in u direction (default (-1, 1)).
    v_range : tuple, optional
        Parameter range in v direction (default (-1, 1)).

    Returns
    -------
    list
        Plane surface definition: ['plane_surface', origin, metadata_dict]
    """
    if isinstance(origin, (list, tuple)) and len(origin) >= 3:
        orig = point(origin[0], origin[1], origin[2])
    else:
        orig = point(origin)

    # Normalize normal vector
    nx, ny, nz = float(normal[0]), float(normal[1]), float(normal[2])
    mag = sqrt(nx*nx + ny*ny + nz*nz)
    if mag < epsilon:
        raise ValueError("Normal vector cannot be zero")
    norm = [nx/mag, ny/mag, nz/mag, 0.0]

    # Compute local u, v axes (tangent vectors)
    # Choose u axis perpendicular to normal
    if abs(norm[2]) < 0.9:
        u_axis = [norm[1], -norm[0], 0.0, 0.0]
    else:
        u_axis = [0.0, norm[2], -norm[1], 0.0]
    u_mag = sqrt(u_axis[0]**2 + u_axis[1]**2 + u_axis[2]**2)
    u_axis = [u_axis[0]/u_mag, u_axis[1]/u_mag, u_axis[2]/u_mag, 0.0]

    # v axis = normal x u_axis
    v_axis = [
        norm[1]*u_axis[2] - norm[2]*u_axis[1],
        norm[2]*u_axis[0] - norm[0]*u_axis[2],
        norm[0]*u_axis[1] - norm[1]*u_axis[0],
        0.0
    ]

    meta = {
        'normal': norm,
        'u_axis': u_axis,
        'v_axis': v_axis,
        'u_range': tuple(u_range),
        'v_range': tuple(v_range),
    }

    return ['plane_surface', orig, meta]


def is_plane_surface(obj):
    """Return True if obj is a plane surface."""
    return (isinstance(obj, list) and len(obj) == 3 and obj[0] == 'plane_surface'
            and ispoint(obj[1]) and isinstance(obj[2], dict))


def evaluate_plane_surface(surf, u, v):
    """Evaluate a point on a plane surface at parameters (u, v).

    Parameters
    ----------
    surf : plane_surface
        The plane surface.
    u, v : float
        Parameter values.

    Returns
    -------
    point
        The point on the surface.
    """
    if not is_plane_surface(surf):
        raise ValueError("Not a plane surface")

    origin = surf[1]
    meta = surf[2]
    u_axis = meta['u_axis']
    v_axis = meta['v_axis']

    return point(
        origin[0] + u * u_axis[0] + v * v_axis[0],
        origin[1] + u * u_axis[1] + v * v_axis[1],
        origin[2] + u * u_axis[2] + v * v_axis[2]
    )


def plane_surface_normal(surf, u, v):
    """Return the normal vector at (u, v) on a plane surface.

    For a plane, the normal is constant everywhere.
    """
    if not is_plane_surface(surf):
        raise ValueError("Not a plane surface")
    return deepcopy(surf[2]['normal'])


# -----------------------------------------------------------------------------
# Sphere Surface
# -----------------------------------------------------------------------------

def sphere_surface(center, radius, *, u_range=(0.0, 2*pi), v_range=(-pi/2, pi/2)):
    """Create an analytic sphere surface.

    The sphere is parameterized as:
    - u: longitude angle (0 to 2*pi)
    - v: latitude angle (-pi/2 to pi/2)

    Parameters
    ----------
    center : point
        Center of the sphere.
    radius : float
        Radius of the sphere.
    u_range : tuple, optional
        Longitude parameter range (default (0, 2*pi)).
    v_range : tuple, optional
        Latitude parameter range (default (-pi/2, pi/2)).

    Returns
    -------
    list
        Sphere surface definition: ['sphere_surface', center, metadata_dict]
    """
    if isinstance(center, (list, tuple)) and len(center) >= 3:
        cen = point(center[0], center[1], center[2])
    else:
        cen = point(center)

    if radius <= 0:
        raise ValueError("Radius must be positive")

    meta = {
        'radius': float(radius),
        'u_range': tuple(u_range),
        'v_range': tuple(v_range),
    }

    return ['sphere_surface', cen, meta]


def is_sphere_surface(obj):
    """Return True if obj is a sphere surface."""
    return (isinstance(obj, list) and len(obj) == 3 and obj[0] == 'sphere_surface'
            and ispoint(obj[1]) and isinstance(obj[2], dict) and 'radius' in obj[2])


def evaluate_sphere_surface(surf, u, v):
    """Evaluate a point on a sphere surface at parameters (u, v).

    Parameters
    ----------
    surf : sphere_surface
        The sphere surface.
    u : float
        Longitude angle (radians).
    v : float
        Latitude angle (radians).

    Returns
    -------
    point
        The point on the surface.
    """
    if not is_sphere_surface(surf):
        raise ValueError("Not a sphere surface")

    center = surf[1]
    r = surf[2]['radius']

    # Spherical to Cartesian
    cos_v = cos(v)
    x = center[0] + r * cos_v * cos(u)
    y = center[1] + r * cos_v * sin(u)
    z = center[2] + r * sin(v)

    return point(x, y, z)


def sphere_surface_normal(surf, u, v):
    """Return the outward normal vector at (u, v) on a sphere surface."""
    if not is_sphere_surface(surf):
        raise ValueError("Not a sphere surface")

    cos_v = cos(v)
    nx = cos_v * cos(u)
    ny = cos_v * sin(u)
    nz = sin(v)

    return [nx, ny, nz, 0.0]


# -----------------------------------------------------------------------------
# Cylinder Surface
# -----------------------------------------------------------------------------

def cylinder_surface(axis_point, axis_direction, radius, *,
                     u_range=(0.0, 2*pi), v_range=(0.0, 1.0)):
    """Create an analytic cylinder surface.

    The cylinder is parameterized as:
    - u: angular position around the axis (0 to 2*pi)
    - v: position along the axis (0 to 1 maps to v_range)

    Parameters
    ----------
    axis_point : point
        A point on the cylinder axis.
    axis_direction : vector
        Direction of the cylinder axis (will be normalized).
    radius : float
        Radius of the cylinder.
    u_range : tuple, optional
        Angular parameter range (default (0, 2*pi)).
    v_range : tuple, optional
        Height parameter range (default (0, 1)).

    Returns
    -------
    list
        Cylinder surface definition.
    """
    if isinstance(axis_point, (list, tuple)) and len(axis_point) >= 3:
        origin = point(axis_point[0], axis_point[1], axis_point[2])
    else:
        origin = point(axis_point)

    # Normalize axis direction
    ax, ay, az = float(axis_direction[0]), float(axis_direction[1]), float(axis_direction[2])
    mag = sqrt(ax*ax + ay*ay + az*az)
    if mag < epsilon:
        raise ValueError("Axis direction cannot be zero")
    axis = [ax/mag, ay/mag, az/mag, 0.0]

    if radius <= 0:
        raise ValueError("Radius must be positive")

    # Compute radial reference direction (perpendicular to axis)
    if abs(axis[2]) < 0.9:
        ref = [axis[1], -axis[0], 0.0, 0.0]
    else:
        ref = [0.0, axis[2], -axis[1], 0.0]
    ref_mag = sqrt(ref[0]**2 + ref[1]**2 + ref[2]**2)
    ref = [ref[0]/ref_mag, ref[1]/ref_mag, ref[2]/ref_mag, 0.0]

    meta = {
        'axis': axis,
        'radius': float(radius),
        'ref_direction': ref,
        'u_range': tuple(u_range),
        'v_range': tuple(v_range),
    }

    return ['cylinder_surface', origin, meta]


def is_cylinder_surface(obj):
    """Return True if obj is a cylinder surface."""
    return (isinstance(obj, list) and len(obj) == 3 and obj[0] == 'cylinder_surface'
            and ispoint(obj[1]) and isinstance(obj[2], dict) and 'radius' in obj[2])


def evaluate_cylinder_surface(surf, u, v):
    """Evaluate a point on a cylinder surface at parameters (u, v).

    Parameters
    ----------
    surf : cylinder_surface
        The cylinder surface.
    u : float
        Angular position (radians).
    v : float
        Axial position.

    Returns
    -------
    point
        The point on the surface.
    """
    if not is_cylinder_surface(surf):
        raise ValueError("Not a cylinder surface")

    origin = surf[1]
    meta = surf[2]
    axis = meta['axis']
    r = meta['radius']
    ref = meta['ref_direction']

    # Compute perpendicular direction at angle u
    # perp = ref * cos(u) + (axis x ref) * sin(u)
    cross = [
        axis[1]*ref[2] - axis[2]*ref[1],
        axis[2]*ref[0] - axis[0]*ref[2],
        axis[0]*ref[1] - axis[1]*ref[0]
    ]

    cos_u = cos(u)
    sin_u = sin(u)

    # Point = origin + v*axis + r*(ref*cos(u) + cross*sin(u))
    x = origin[0] + v * axis[0] + r * (ref[0] * cos_u + cross[0] * sin_u)
    y = origin[1] + v * axis[1] + r * (ref[1] * cos_u + cross[1] * sin_u)
    z = origin[2] + v * axis[2] + r * (ref[2] * cos_u + cross[2] * sin_u)

    return point(x, y, z)


def cylinder_surface_normal(surf, u, v):
    """Return the outward normal vector at (u, v) on a cylinder surface."""
    if not is_cylinder_surface(surf):
        raise ValueError("Not a cylinder surface")

    meta = surf[2]
    ref = meta['ref_direction']
    axis = meta['axis']

    # Normal = ref * cos(u) + (axis x ref) * sin(u)
    cross = [
        axis[1]*ref[2] - axis[2]*ref[1],
        axis[2]*ref[0] - axis[0]*ref[2],
        axis[0]*ref[1] - axis[1]*ref[0]
    ]

    cos_u = cos(u)
    sin_u = sin(u)

    return [
        ref[0] * cos_u + cross[0] * sin_u,
        ref[1] * cos_u + cross[1] * sin_u,
        ref[2] * cos_u + cross[2] * sin_u,
        0.0
    ]


# -----------------------------------------------------------------------------
# Cone Surface
# -----------------------------------------------------------------------------

def cone_surface(apex, axis_direction, half_angle, *,
                 u_range=(0.0, 2*pi), v_range=(0.0, 1.0)):
    """Create an analytic cone surface.

    The cone is parameterized as:
    - u: angular position around the axis (0 to 2*pi)
    - v: distance from apex along the surface

    Parameters
    ----------
    apex : point
        Apex (tip) of the cone.
    axis_direction : vector
        Direction of the cone axis (from apex, will be normalized).
    half_angle : float
        Half-angle of the cone in radians.
    u_range : tuple, optional
        Angular parameter range (default (0, 2*pi)).
    v_range : tuple, optional
        Distance parameter range from apex (default (0, 1)).

    Returns
    -------
    list
        Cone surface definition.
    """
    if isinstance(apex, (list, tuple)) and len(apex) >= 3:
        ap = point(apex[0], apex[1], apex[2])
    else:
        ap = point(apex)

    # Normalize axis direction
    ax, ay, az = float(axis_direction[0]), float(axis_direction[1]), float(axis_direction[2])
    mag = sqrt(ax*ax + ay*ay + az*az)
    if mag < epsilon:
        raise ValueError("Axis direction cannot be zero")
    axis = [ax/mag, ay/mag, az/mag, 0.0]

    if half_angle <= 0 or half_angle >= pi/2:
        raise ValueError("Half angle must be in (0, pi/2)")

    # Compute radial reference direction
    if abs(axis[2]) < 0.9:
        ref = [axis[1], -axis[0], 0.0, 0.0]
    else:
        ref = [0.0, axis[2], -axis[1], 0.0]
    ref_mag = sqrt(ref[0]**2 + ref[1]**2 + ref[2]**2)
    ref = [ref[0]/ref_mag, ref[1]/ref_mag, ref[2]/ref_mag, 0.0]

    meta = {
        'axis': axis,
        'half_angle': float(half_angle),
        'ref_direction': ref,
        'u_range': tuple(u_range),
        'v_range': tuple(v_range),
    }

    return ['cone_surface', ap, meta]


def is_cone_surface(obj):
    """Return True if obj is a cone surface."""
    return (isinstance(obj, list) and len(obj) == 3 and obj[0] == 'cone_surface'
            and ispoint(obj[1]) and isinstance(obj[2], dict) and 'half_angle' in obj[2])


def evaluate_cone_surface(surf, u, v):
    """Evaluate a point on a cone surface at parameters (u, v).

    Parameters
    ----------
    surf : cone_surface
        The cone surface.
    u : float
        Angular position (radians).
    v : float
        Distance from apex along surface.

    Returns
    -------
    point
        The point on the surface.
    """
    if not is_cone_surface(surf):
        raise ValueError("Not a cone surface")

    apex = surf[1]
    meta = surf[2]
    axis = meta['axis']
    half_angle = meta['half_angle']
    ref = meta['ref_direction']

    # Compute perpendicular direction at angle u
    cross = [
        axis[1]*ref[2] - axis[2]*ref[1],
        axis[2]*ref[0] - axis[0]*ref[2],
        axis[0]*ref[1] - axis[1]*ref[0]
    ]

    cos_u = cos(u)
    sin_u = sin(u)
    sin_a = sin(half_angle)
    cos_a = cos(half_angle)

    # Point = apex + v*(cos_a*axis + sin_a*(ref*cos(u) + cross*sin(u)))
    x = apex[0] + v * (cos_a * axis[0] + sin_a * (ref[0] * cos_u + cross[0] * sin_u))
    y = apex[1] + v * (cos_a * axis[1] + sin_a * (ref[1] * cos_u + cross[1] * sin_u))
    z = apex[2] + v * (cos_a * axis[2] + sin_a * (ref[2] * cos_u + cross[2] * sin_u))

    return point(x, y, z)


def cone_surface_normal(surf, u, v):
    """Return the outward normal vector at (u, v) on a cone surface."""
    if not is_cone_surface(surf):
        raise ValueError("Not a cone surface")

    meta = surf[2]
    axis = meta['axis']
    half_angle = meta['half_angle']
    ref = meta['ref_direction']

    cross = [
        axis[1]*ref[2] - axis[2]*ref[1],
        axis[2]*ref[0] - axis[0]*ref[2],
        axis[0]*ref[1] - axis[1]*ref[0]
    ]

    cos_u = cos(u)
    sin_u = sin(u)
    sin_a = sin(half_angle)
    cos_a = cos(half_angle)

    # Normal = sin_a*axis - cos_a*(ref*cos(u) + cross*sin(u))
    # (outward normal for cone surface)
    return [
        sin_a * axis[0] - cos_a * (ref[0] * cos_u + cross[0] * sin_u),
        sin_a * axis[1] - cos_a * (ref[1] * cos_u + cross[1] * sin_u),
        sin_a * axis[2] - cos_a * (ref[2] * cos_u + cross[2] * sin_u),
        0.0
    ]


# -----------------------------------------------------------------------------
# Torus Surface
# -----------------------------------------------------------------------------

def torus_surface(center, axis_direction, major_radius, minor_radius, *,
                  u_range=(0.0, 2*pi), v_range=(0.0, 2*pi)):
    """Create an analytic torus surface.

    The torus is parameterized as:
    - u: angle around the major circle (0 to 2*pi)
    - v: angle around the minor circle (0 to 2*pi)

    Parameters
    ----------
    center : point
        Center of the torus.
    axis_direction : vector
        Direction of the torus axis (will be normalized).
    major_radius : float
        Major radius (distance from center to tube center).
    minor_radius : float
        Minor radius (tube radius).
    u_range : tuple, optional
        Major angle parameter range (default (0, 2*pi)).
    v_range : tuple, optional
        Minor angle parameter range (default (0, 2*pi)).

    Returns
    -------
    list
        Torus surface definition.
    """
    if isinstance(center, (list, tuple)) and len(center) >= 3:
        cen = point(center[0], center[1], center[2])
    else:
        cen = point(center)

    # Normalize axis direction
    ax, ay, az = float(axis_direction[0]), float(axis_direction[1]), float(axis_direction[2])
    mag = sqrt(ax*ax + ay*ay + az*az)
    if mag < epsilon:
        raise ValueError("Axis direction cannot be zero")
    axis = [ax/mag, ay/mag, az/mag, 0.0]

    if major_radius <= 0 or minor_radius <= 0:
        raise ValueError("Radii must be positive")
    if minor_radius >= major_radius:
        raise ValueError("Minor radius must be less than major radius")

    # Compute radial reference direction
    if abs(axis[2]) < 0.9:
        ref = [axis[1], -axis[0], 0.0, 0.0]
    else:
        ref = [0.0, axis[2], -axis[1], 0.0]
    ref_mag = sqrt(ref[0]**2 + ref[1]**2 + ref[2]**2)
    ref = [ref[0]/ref_mag, ref[1]/ref_mag, ref[2]/ref_mag, 0.0]

    meta = {
        'axis': axis,
        'major_radius': float(major_radius),
        'minor_radius': float(minor_radius),
        'ref_direction': ref,
        'u_range': tuple(u_range),
        'v_range': tuple(v_range),
    }

    return ['torus_surface', cen, meta]


def is_torus_surface(obj):
    """Return True if obj is a torus surface."""
    return (isinstance(obj, list) and len(obj) == 3 and obj[0] == 'torus_surface'
            and ispoint(obj[1]) and isinstance(obj[2], dict)
            and 'major_radius' in obj[2] and 'minor_radius' in obj[2])


def evaluate_torus_surface(surf, u, v):
    """Evaluate a point on a torus surface at parameters (u, v).

    Parameters
    ----------
    surf : torus_surface
        The torus surface.
    u : float
        Major angle (radians).
    v : float
        Minor angle (radians).

    Returns
    -------
    point
        The point on the surface.
    """
    if not is_torus_surface(surf):
        raise ValueError("Not a torus surface")

    center = surf[1]
    meta = surf[2]
    axis = meta['axis']
    R = meta['major_radius']
    r = meta['minor_radius']
    ref = meta['ref_direction']

    # Compute perpendicular directions
    cross = [
        axis[1]*ref[2] - axis[2]*ref[1],
        axis[2]*ref[0] - axis[0]*ref[2],
        axis[0]*ref[1] - axis[1]*ref[0]
    ]

    cos_u = cos(u)
    sin_u = sin(u)
    cos_v = cos(v)
    sin_v = sin(v)

    # Direction from center to tube center
    dir_x = ref[0] * cos_u + cross[0] * sin_u
    dir_y = ref[1] * cos_u + cross[1] * sin_u
    dir_z = ref[2] * cos_u + cross[2] * sin_u

    # Point = center + (R + r*cos(v))*dir + r*sin(v)*axis
    factor = R + r * cos_v
    x = center[0] + factor * dir_x + r * sin_v * axis[0]
    y = center[1] + factor * dir_y + r * sin_v * axis[1]
    z = center[2] + factor * dir_z + r * sin_v * axis[2]

    return point(x, y, z)


def torus_surface_normal(surf, u, v):
    """Return the outward normal vector at (u, v) on a torus surface."""
    if not is_torus_surface(surf):
        raise ValueError("Not a torus surface")

    meta = surf[2]
    axis = meta['axis']
    ref = meta['ref_direction']

    cross = [
        axis[1]*ref[2] - axis[2]*ref[1],
        axis[2]*ref[0] - axis[0]*ref[2],
        axis[0]*ref[1] - axis[1]*ref[0]
    ]

    cos_u = cos(u)
    sin_u = sin(u)
    cos_v = cos(v)
    sin_v = sin(v)

    # Direction from center to tube center
    dir_x = ref[0] * cos_u + cross[0] * sin_u
    dir_y = ref[1] * cos_u + cross[1] * sin_u
    dir_z = ref[2] * cos_u + cross[2] * sin_u

    # Normal = cos(v)*dir + sin(v)*axis
    return [
        cos_v * dir_x + sin_v * axis[0],
        cos_v * dir_y + sin_v * axis[1],
        cos_v * dir_z + sin_v * axis[2],
        0.0
    ]


# -----------------------------------------------------------------------------
# NURBS/B-Spline Surface
# -----------------------------------------------------------------------------

def bspline_surface(control_points, u_knots, v_knots, u_degree, v_degree, *,
                    weights=None, u_range=None, v_range=None):
    """Create a NURBS/B-spline surface.

    The surface is defined by a grid of control points, knot vectors in both
    parameter directions, and optional weights (for rational B-splines).

    Parameters
    ----------
    control_points : list of list of points
        2D grid of control points [v_rows][u_cols]. Each point is [x, y, z] or [x, y, z, w].
    u_knots : list of float
        Knot vector in u direction.
    v_knots : list of float
        Knot vector in v direction.
    u_degree : int
        Degree of the B-spline in u direction.
    v_degree : int
        Degree of the B-spline in v direction.
    weights : list of list of float, optional
        Weights for rational B-splines. If None, all weights are 1.0 (non-rational).
    u_range : tuple, optional
        Parameter range in u direction. Defaults to (min(u_knots), max(u_knots)).
    v_range : tuple, optional
        Parameter range in v direction. Defaults to (min(v_knots), max(v_knots)).

    Returns
    -------
    list
        B-spline surface definition: ['bspline_surface', control_points, metadata_dict]
    """
    # Validate dimensions
    n_v = len(control_points)
    if n_v < 2:
        raise ValueError("Need at least 2 rows of control points")
    n_u = len(control_points[0])
    if n_u < 2:
        raise ValueError("Need at least 2 columns of control points")
    for row in control_points:
        if len(row) != n_u:
            raise ValueError("All rows must have the same number of control points")

    # Convert control points to standard format
    cpts = []
    for row in control_points:
        cpts_row = []
        for cp in row:
            if len(cp) >= 3:
                cpts_row.append([float(cp[0]), float(cp[1]), float(cp[2]), 1.0])
            else:
                raise ValueError("Control points must have at least 3 coordinates")
        cpts.append(cpts_row)

    # Validate knot vectors
    u_knots = [float(k) for k in u_knots]
    v_knots = [float(k) for k in v_knots]

    expected_u_knots = n_u + u_degree + 1
    expected_v_knots = n_v + v_degree + 1
    if len(u_knots) != expected_u_knots:
        raise ValueError(f"u_knots should have {expected_u_knots} elements, got {len(u_knots)}")
    if len(v_knots) != expected_v_knots:
        raise ValueError(f"v_knots should have {expected_v_knots} elements, got {len(v_knots)}")

    # Handle weights
    if weights is None:
        w = [[1.0] * n_u for _ in range(n_v)]
    else:
        w = [[float(weights[j][i]) for i in range(n_u)] for j in range(n_v)]

    # Default parameter ranges from knot vectors
    if u_range is None:
        u_range = (u_knots[u_degree], u_knots[n_u])
    if v_range is None:
        v_range = (v_knots[v_degree], v_knots[n_v])

    meta = {
        'u_knots': u_knots,
        'v_knots': v_knots,
        'u_degree': int(u_degree),
        'v_degree': int(v_degree),
        'weights': w,
        'n_u': n_u,
        'n_v': n_v,
        'u_range': tuple(u_range),
        'v_range': tuple(v_range),
    }

    return ['bspline_surface', cpts, meta]


def is_bspline_surface(obj):
    """Return True if obj is a B-spline surface."""
    return (isinstance(obj, list) and len(obj) == 3 and obj[0] == 'bspline_surface'
            and isinstance(obj[1], list) and isinstance(obj[2], dict)
            and 'u_knots' in obj[2] and 'v_knots' in obj[2])


def _bspline_basis(knots, i, p, u):
    """Compute B-spline basis function N_{i,p}(u) using Cox-de Boor recursion.

    Parameters
    ----------
    knots : list of float
        Knot vector.
    i : int
        Basis function index.
    p : int
        Degree.
    u : float
        Parameter value.

    Returns
    -------
    float
        Basis function value.
    """
    if p == 0:
        # Base case: piecewise constant
        if knots[i] <= u < knots[i + 1]:
            return 1.0
        # Handle end of knot span
        if i + 1 == len(knots) - 1 and u == knots[i + 1]:
            return 1.0
        return 0.0

    # Recursive case
    N1 = 0.0
    denom1 = knots[i + p] - knots[i]
    if denom1 > 1e-12:
        N1 = ((u - knots[i]) / denom1) * _bspline_basis(knots, i, p - 1, u)

    N2 = 0.0
    denom2 = knots[i + p + 1] - knots[i + 1]
    if denom2 > 1e-12:
        N2 = ((knots[i + p + 1] - u) / denom2) * _bspline_basis(knots, i + 1, p - 1, u)

    return N1 + N2


def evaluate_bspline_surface(surf, u, v):
    """Evaluate a point on a B-spline surface at parameters (u, v).

    Parameters
    ----------
    surf : bspline_surface
        The B-spline surface.
    u, v : float
        Parameter values.

    Returns
    -------
    point
        The point on the surface.
    """
    if not is_bspline_surface(surf):
        raise ValueError("Not a B-spline surface")

    cpts = surf[1]
    meta = surf[2]
    u_knots = meta['u_knots']
    v_knots = meta['v_knots']
    p = meta['u_degree']
    q = meta['v_degree']
    weights = meta['weights']
    n_u = meta['n_u']
    n_v = meta['n_v']

    # Compute weighted sum
    x, y, z, w_sum = 0.0, 0.0, 0.0, 0.0

    for j in range(n_v):
        Nv = _bspline_basis(v_knots, j, q, v)
        if Nv == 0.0:
            continue
        for i in range(n_u):
            Nu = _bspline_basis(u_knots, i, p, u)
            if Nu == 0.0:
                continue
            wij = weights[j][i]
            basis = Nu * Nv * wij
            cp = cpts[j][i]
            x += basis * cp[0]
            y += basis * cp[1]
            z += basis * cp[2]
            w_sum += basis

    if w_sum < 1e-12:
        # Fallback: return first control point
        return point(cpts[0][0][0], cpts[0][0][1], cpts[0][0][2])

    return point(x / w_sum, y / w_sum, z / w_sum)


def bspline_surface_normal(surf, u, v, delta=1e-6):
    """Return the normal vector at (u, v) on a B-spline surface.

    Computed via finite differences of partial derivatives.
    """
    if not is_bspline_surface(surf):
        raise ValueError("Not a B-spline surface")

    meta = surf[2]
    u_range = meta['u_range']
    v_range = meta['v_range']

    # Clamp delta to stay in bounds
    du = min(delta, (u_range[1] - u_range[0]) * 0.01)
    dv = min(delta, (v_range[1] - v_range[0]) * 0.01)

    # Sample points for partial derivatives
    p0 = evaluate_bspline_surface(surf, u, v)
    pu = evaluate_bspline_surface(surf, min(u + du, u_range[1]), v)
    pv = evaluate_bspline_surface(surf, u, min(v + dv, v_range[1]))

    # Tangent vectors
    tu = [(pu[i] - p0[i]) / du for i in range(3)]
    tv = [(pv[i] - p0[i]) / dv for i in range(3)]

    # Normal = tu x tv
    nx = tu[1] * tv[2] - tu[2] * tv[1]
    ny = tu[2] * tv[0] - tu[0] * tv[2]
    nz = tu[0] * tv[1] - tu[1] * tv[0]

    # Normalize
    mag = sqrt(nx * nx + ny * ny + nz * nz)
    if mag < 1e-12:
        return [0.0, 0.0, 1.0, 0.0]

    return [nx / mag, ny / mag, nz / mag, 0.0]


# -----------------------------------------------------------------------------
# Tessellated Surface (fallback for non-analytic geometry)
# -----------------------------------------------------------------------------

def tessellated_surface(vertices, normals, faces, *, u_range=(0.0, 1.0), v_range=(0.0, 1.0)):
    """Create a pre-tessellated surface (fallback for non-analytic geometry).

    This surface type stores an already-tessellated mesh. It cannot be
    re-tessellated at different resolutions but provides a consistent
    interface with other analytic surfaces.

    Parameters
    ----------
    vertices : list of points
        List of vertex positions.
    normals : list of vectors
        List of normal vectors (one per vertex or one per face).
    faces : list of [i, j, k]
        Triangle indices into vertices.
    u_range : tuple, optional
        Nominal parameter range in u direction.
    v_range : tuple, optional
        Nominal parameter range in v direction.

    Returns
    -------
    list
        Tessellated surface definition.
    """
    # Store vertices as standard points
    verts = []
    for v in vertices:
        if len(v) >= 3:
            verts.append([float(v[0]), float(v[1]), float(v[2]), 1.0])
        else:
            raise ValueError("Vertices must have at least 3 coordinates")

    norms = []
    for n in normals:
        if len(n) >= 3:
            norms.append([float(n[0]), float(n[1]), float(n[2]), 0.0])
        else:
            raise ValueError("Normals must have at least 3 coordinates")

    face_list = [[int(f[0]), int(f[1]), int(f[2])] for f in faces]

    meta = {
        'u_range': tuple(u_range),
        'v_range': tuple(v_range),
        'is_tessellated': True,
    }

    return ['tessellated_surface', {'vertices': verts, 'normals': norms, 'faces': face_list}, meta]


def is_tessellated_surface(obj):
    """Return True if obj is a pre-tessellated surface."""
    return (isinstance(obj, list) and len(obj) == 3 and obj[0] == 'tessellated_surface'
            and isinstance(obj[1], dict) and isinstance(obj[2], dict))


def evaluate_tessellated_surface(surf, u, v):
    """Evaluate a point on a tessellated surface at parameters (u, v).

    Since the surface is pre-tessellated, this uses bilinear interpolation
    over the bounding box as an approximation.
    """
    if not is_tessellated_surface(surf):
        raise ValueError("Not a tessellated surface")

    data = surf[1]
    vertices = data['vertices']
    if not vertices:
        raise ValueError("Tessellated surface has no vertices")

    # Simple approximation: find bounds and interpolate
    xs = [v[0] for v in vertices]
    ys = [v[1] for v in vertices]
    zs = [v[2] for v in vertices]

    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    z_min, z_max = min(zs), max(zs)

    # Map u, v to position (crude approximation)
    meta = surf[2]
    u_range = meta['u_range']
    v_range = meta['v_range']

    t_u = (u - u_range[0]) / (u_range[1] - u_range[0]) if u_range[1] != u_range[0] else 0.5
    t_v = (v - v_range[0]) / (v_range[1] - v_range[0]) if v_range[1] != v_range[0] else 0.5

    x = x_min + t_u * (x_max - x_min)
    y = y_min + t_v * (y_max - y_min)
    z = (z_min + z_max) / 2  # Approximate

    return point(x, y, z)


def tessellated_surface_normal(surf, u, v):
    """Return a nominal normal for a tessellated surface.

    Returns average of all normals as an approximation.
    """
    if not is_tessellated_surface(surf):
        raise ValueError("Not a tessellated surface")

    data = surf[1]
    normals = data.get('normals', [])
    if not normals:
        return [0.0, 0.0, 1.0, 0.0]

    # Average all normals
    nx, ny, nz = 0.0, 0.0, 0.0
    for n in normals:
        nx += n[0]
        ny += n[1]
        nz += n[2]

    mag = sqrt(nx * nx + ny * ny + nz * nz)
    if mag < 1e-12:
        return [0.0, 0.0, 1.0, 0.0]

    return [nx / mag, ny / mag, nz / mag, 0.0]


# -----------------------------------------------------------------------------
# Generic Surface Operations
# -----------------------------------------------------------------------------

def is_analytic_surface(obj):
    """Return True if obj is any type of analytic surface."""
    return (is_plane_surface(obj) or is_sphere_surface(obj) or
            is_cylinder_surface(obj) or is_cone_surface(obj) or
            is_torus_surface(obj) or is_bspline_surface(obj) or
            is_tessellated_surface(obj))


def evaluate_surface(surf, u, v):
    """Evaluate a point on an analytic surface at parameters (u, v).

    This is a generic dispatcher for all analytic surface types.
    """
    if is_plane_surface(surf):
        return evaluate_plane_surface(surf, u, v)
    elif is_sphere_surface(surf):
        return evaluate_sphere_surface(surf, u, v)
    elif is_cylinder_surface(surf):
        return evaluate_cylinder_surface(surf, u, v)
    elif is_cone_surface(surf):
        return evaluate_cone_surface(surf, u, v)
    elif is_torus_surface(surf):
        return evaluate_torus_surface(surf, u, v)
    elif is_bspline_surface(surf):
        return evaluate_bspline_surface(surf, u, v)
    elif is_tessellated_surface(surf):
        return evaluate_tessellated_surface(surf, u, v)
    else:
        raise ValueError("Not an analytic surface")


def surface_normal(surf, u, v):
    """Return the normal vector at (u, v) on an analytic surface.

    This is a generic dispatcher for all analytic surface types.
    """
    if is_plane_surface(surf):
        return plane_surface_normal(surf, u, v)
    elif is_sphere_surface(surf):
        return sphere_surface_normal(surf, u, v)
    elif is_cylinder_surface(surf):
        return cylinder_surface_normal(surf, u, v)
    elif is_cone_surface(surf):
        return cone_surface_normal(surf, u, v)
    elif is_torus_surface(surf):
        return torus_surface_normal(surf, u, v)
    elif is_bspline_surface(surf):
        return bspline_surface_normal(surf, u, v)
    elif is_tessellated_surface(surf):
        return tessellated_surface_normal(surf, u, v)
    else:
        raise ValueError("Not an analytic surface")


def tessellate_surface(surf, *, u_divisions=16, v_divisions=16):
    """Convert an analytic surface to a tessellated mesh surface.

    Parameters
    ----------
    surf : analytic surface
        Any analytic surface type.
    u_divisions : int
        Number of divisions in u parameter.
    v_divisions : int
        Number of divisions in v parameter.

    Returns
    -------
    list
        Tessellated surface: ['surface', vertices, normals, faces, boundary, holes]
    """
    if not is_analytic_surface(surf):
        raise ValueError("Not an analytic surface")

    # Special case: pre-tessellated surfaces return their stored mesh
    if is_tessellated_surface(surf):
        data = surf[1]
        return ['surface', data['vertices'], data['normals'], data['faces'], [], []]

    meta = surf[2]
    u_range = meta.get('u_range', (0.0, 1.0))
    v_range = meta.get('v_range', (0.0, 1.0))

    u_min, u_max = u_range
    v_min, v_max = v_range

    vertices = []
    normals = []

    # Generate grid of vertices and normals
    for j in range(v_divisions + 1):
        v = v_min + (v_max - v_min) * j / v_divisions
        for i in range(u_divisions + 1):
            u = u_min + (u_max - u_min) * i / u_divisions
            pt = evaluate_surface(surf, u, v)
            nrm = surface_normal(surf, u, v)
            vertices.append(pt)
            normals.append(nrm)

    # Generate triangle faces
    faces = []
    for j in range(v_divisions):
        for i in range(u_divisions):
            # Indices of quad corners
            i00 = j * (u_divisions + 1) + i
            i10 = i00 + 1
            i01 = i00 + (u_divisions + 1)
            i11 = i01 + 1

            # Two triangles per quad
            faces.append([i00, i10, i11])
            faces.append([i00, i11, i01])

    # Compute boundary (outer edge vertices)
    boundary = []
    # Bottom edge
    for i in range(u_divisions + 1):
        boundary.append(i)
    # Right edge (excluding first)
    for j in range(1, v_divisions + 1):
        boundary.append(j * (u_divisions + 1) + u_divisions)
    # Top edge (excluding last, reversed)
    for i in range(u_divisions - 1, -1, -1):
        boundary.append(v_divisions * (u_divisions + 1) + i)
    # Left edge (excluding first and last, reversed)
    for j in range(v_divisions - 1, 0, -1):
        boundary.append(j * (u_divisions + 1))

    return ['surface', vertices, normals, faces, boundary, []]
