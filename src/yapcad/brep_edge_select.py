"""Edge selection helper functions for yapCAD BREP operations.

This module provides functions to select edges from BREP solids based on
various geometric criteria such as direction, length, position, and
association with specific faces.

These functions are particularly useful for applying selective fillet or
chamfer operations to specific edges rather than all edges of a solid.

Example usage:
    from yapcad.brep import brep_from_solid, fillet_edges
    from yapcad.brep_edge_select import select_vertical_edges
    from yapcad.geom3d_util import prism

    # Create a pocket (prism with a hole)
    solid = prism(20, 20, 10)
    brep = brep_from_solid(solid)

    # Select only vertical edges
    vertical_edges = select_vertical_edges(brep)

    # Apply fillet only to vertical edges
    filleted = fillet_edges(brep, vertical_edges, radius=1.0)

Copyright (c) 2025 Richard W. DeVaul
Copyright (c) 2025 yapCAD contributors
MIT License
"""

import math
from typing import List, Optional, Tuple, Union

from yapcad.brep import (
    BrepSolid,
    BrepEdge,
    require_occ,
    occ_available,
    _get_all_edges,
)
from yapcad.geom import point, vclose, epsilon

# Lazy imports for OCC modules
_BRepAdaptor_Curve = None
_GCPnts_AbscissaPoint = None
_GProp_GProps = None
_brepgprop = None
_gp_Vec = None
_gp_Pnt = None
_TopAbs_EDGE = None
_TopAbs_FACE = None
_TopExp_Explorer = None
_TopExp = None
_TopTools_IndexedDataMapOfShapeListOfShape = None
_topods = None
_BRep_Tool = None


def _ensure_occ_imports():
    """Lazily import OCC modules when needed."""
    global _BRepAdaptor_Curve, _GCPnts_AbscissaPoint, _GProp_GProps
    global _brepgprop, _gp_Vec, _gp_Pnt, _TopAbs_EDGE, _TopAbs_FACE
    global _TopExp_Explorer, _TopExp, _TopTools_IndexedDataMapOfShapeListOfShape
    global _topods, _BRep_Tool

    if _BRepAdaptor_Curve is not None:
        return  # Already imported

    require_occ()

    from OCC.Core.BRepAdaptor import BRepAdaptor_Curve
    from OCC.Core.GProp import GProp_GProps
    from OCC.Core.BRepGProp import brepgprop
    from OCC.Core.gp import gp_Vec, gp_Pnt
    from OCC.Core.TopAbs import TopAbs_EDGE, TopAbs_FACE
    from OCC.Core.TopExp import TopExp_Explorer
    from OCC.Core.TopoDS import topods
    from OCC.Core.BRep import BRep_Tool

    _BRepAdaptor_Curve = BRepAdaptor_Curve
    _GProp_GProps = GProp_GProps
    _brepgprop = brepgprop
    _gp_Vec = gp_Vec
    _gp_Pnt = gp_Pnt
    _TopAbs_EDGE = TopAbs_EDGE
    _TopAbs_FACE = TopAbs_FACE
    _TopExp_Explorer = TopExp_Explorer
    _topods = topods
    _BRep_Tool = BRep_Tool

    # These may not be available in all OCC versions
    try:
        from OCC.Core.TopExp import TopExp
        _TopExp = TopExp
    except ImportError:
        _TopExp = None

    try:
        from OCC.Core.TopTools import TopTools_IndexedDataMapOfShapeListOfShape
        _TopTools_IndexedDataMapOfShapeListOfShape = TopTools_IndexedDataMapOfShapeListOfShape
    except ImportError:
        _TopTools_IndexedDataMapOfShapeListOfShape = None

    try:
        from OCC.Core.GCPnts import GCPnts_AbscissaPoint
        _GCPnts_AbscissaPoint = GCPnts_AbscissaPoint
    except ImportError:
        _GCPnts_AbscissaPoint = None


def _cast_to_edge(edge):
    """Ensure edge is cast to TopoDS_Edge for use with BRepAdaptor."""
    _ensure_occ_imports()
    if _topods is not None:
        try:
            return _topods.Edge(edge)
        except Exception:
            pass
    return edge


def _get_edge_direction(edge) -> Optional[Tuple[float, float, float]]:
    """Get the direction vector of an edge.

    For linear edges, returns the normalized direction vector.
    For curved edges, returns None (curves don't have a single direction).

    Args:
        edge: TopoDS_Edge object

    Returns:
        Normalized direction vector as (x, y, z) tuple, or None for curved edges
    """
    _ensure_occ_imports()

    try:
        # Cast to TopoDS_Edge if needed
        edge = _cast_to_edge(edge)
        adaptor = _BRepAdaptor_Curve(edge)

        # Get start and end parameters
        first = adaptor.FirstParameter()
        last = adaptor.LastParameter()

        # Check if edge is linear by comparing midpoint to line
        p1 = adaptor.Value(first)
        p2 = adaptor.Value(last)
        pmid = adaptor.Value((first + last) / 2.0)

        # Expected midpoint for a line
        expected_mid_x = (p1.X() + p2.X()) / 2.0
        expected_mid_y = (p1.Y() + p2.Y()) / 2.0
        expected_mid_z = (p1.Z() + p2.Z()) / 2.0

        # Check if actual midpoint matches expected (indicating linear edge)
        tol = 1e-6
        if (abs(pmid.X() - expected_mid_x) > tol or
            abs(pmid.Y() - expected_mid_y) > tol or
            abs(pmid.Z() - expected_mid_z) > tol):
            # Non-linear edge
            return None

        # Calculate direction vector
        dx = p2.X() - p1.X()
        dy = p2.Y() - p1.Y()
        dz = p2.Z() - p1.Z()

        length = math.sqrt(dx*dx + dy*dy + dz*dz)
        if length < 1e-12:
            return None

        return (dx/length, dy/length, dz/length)

    except Exception:
        return None


def _get_edge_endpoints(edge) -> Optional[Tuple[Tuple[float, float, float], Tuple[float, float, float]]]:
    """Get the start and end points of an edge.

    Args:
        edge: TopoDS_Edge object

    Returns:
        Tuple of ((x1, y1, z1), (x2, y2, z2)) or None on error
    """
    _ensure_occ_imports()

    try:
        # Cast to TopoDS_Edge if needed
        edge = _cast_to_edge(edge)
        adaptor = _BRepAdaptor_Curve(edge)
        first = adaptor.FirstParameter()
        last = adaptor.LastParameter()

        p1 = adaptor.Value(first)
        p2 = adaptor.Value(last)

        return ((p1.X(), p1.Y(), p1.Z()), (p2.X(), p2.Y(), p2.Z()))

    except Exception:
        return None


def _get_edge_length(edge) -> float:
    """Calculate the length of an edge.

    Args:
        edge: TopoDS_Edge object

    Returns:
        Length of the edge, or 0.0 on error
    """
    _ensure_occ_imports()

    try:
        # Cast to TopoDS_Edge if needed
        edge = _cast_to_edge(edge)
        props = _GProp_GProps()
        _brepgprop.LinearProperties(edge, props)
        return props.Mass()  # For linear properties, Mass() returns length
    except Exception:
        return 0.0


def _get_edge_midpoint(edge) -> Optional[Tuple[float, float, float]]:
    """Get the midpoint of an edge.

    Args:
        edge: TopoDS_Edge object

    Returns:
        Midpoint as (x, y, z) tuple, or None on error
    """
    _ensure_occ_imports()

    try:
        # Cast to TopoDS_Edge if needed
        edge = _cast_to_edge(edge)
        adaptor = _BRepAdaptor_Curve(edge)
        first = adaptor.FirstParameter()
        last = adaptor.LastParameter()

        pmid = adaptor.Value((first + last) / 2.0)
        return (pmid.X(), pmid.Y(), pmid.Z())

    except Exception:
        return None


def _get_unique_edges(shape) -> list:
    """Get unique edges from a shape using TopTools_IndexedMapOfShape.

    This avoids the duplicate edges that TopExp_Explorer returns
    (which visits each edge once per adjacent face).

    Args:
        shape: TopoDS_Shape object

    Returns:
        List of unique TopoDS_Edge objects
    """
    _ensure_occ_imports()

    from OCC.Core.TopTools import TopTools_IndexedMapOfShape

    edge_map = TopTools_IndexedMapOfShape()

    # pythonocc-core exposes OCCT static methods differently across versions:
    #   * 7.7.x: lowercase singleton `topexp` with plain method names
    #   * newer (7.8+/OCP-style): PascalCase class with `_s`-suffixed statics
    # Try both so a silent fallback can't mask the dedup path (which caused
    # duplicate edges — one per adjacent face — to leak into edge selection).
    try:
        from OCC.Core.TopExp import topexp  # 7.7.x style
        topexp.MapShapes(shape, _TopAbs_EDGE, edge_map)
    except ImportError:
        from OCC.Core.TopExp import TopExp  # newer style
        TopExp.MapShapes_s(shape, _TopAbs_EDGE, edge_map)

    # Map size accessor also differs across versions (Extent vs Size).
    count = edge_map.Extent() if hasattr(edge_map, "Extent") else edge_map.Size()
    return [edge_map.FindKey(i) for i in range(1, count + 1)]


def get_all_edges(brep_solid: BrepSolid) -> List[BrepEdge]:
    """Get all unique edges from a BREP solid as BrepEdge objects.

    Args:
        brep_solid: A BrepSolid object

    Returns:
        List of BrepEdge objects (unique, no duplicates)
    """
    require_occ()

    raw_edges = _get_unique_edges(brep_solid.shape)
    return [BrepEdge(e) for e in raw_edges]


def select_vertical_edges(brep_solid: BrepSolid,
                          tolerance_deg: float = 1.0) -> List[BrepEdge]:
    """Select edges that are parallel to the Z axis (vertical).

    Args:
        brep_solid: A BrepSolid object
        tolerance_deg: Angular tolerance in degrees (default 1.0)

    Returns:
        List of BrepEdge objects that are vertical
    """
    return select_edges_by_direction(brep_solid, (0, 0, 1), tolerance_deg)


def select_horizontal_edges(brep_solid: BrepSolid,
                            tolerance_deg: float = 1.0) -> List[BrepEdge]:
    """Select edges that are perpendicular to the Z axis (horizontal).

    This includes edges in the XY plane at any Z height.

    Args:
        brep_solid: A BrepSolid object
        tolerance_deg: Angular tolerance in degrees (default 1.0)

    Returns:
        List of BrepEdge objects that are horizontal
    """
    require_occ()
    _ensure_occ_imports()

    tolerance_rad = math.radians(tolerance_deg)
    z_axis = (0.0, 0.0, 1.0)

    result = []
    raw_edges = _get_unique_edges(brep_solid.shape)

    for edge in raw_edges:
        direction = _get_edge_direction(edge)
        if direction is None:
            continue  # Skip curved edges

        # Calculate dot product with Z axis
        # For horizontal edges, this should be near 0
        dot = abs(direction[0] * z_axis[0] +
                  direction[1] * z_axis[1] +
                  direction[2] * z_axis[2])

        # dot = cos(angle), so for perpendicular edges, dot should be near 0
        # angle = 90 degrees means cos(90) = 0
        # We check if angle is within tolerance of 90 degrees
        # cos(90 - tol) = sin(tol) ~ tol for small angles
        if dot <= math.sin(tolerance_rad):
            result.append(BrepEdge(edge))

    return result


def select_edges_by_direction(brep_solid: BrepSolid,
                              direction: Union[Tuple[float, float, float], list],
                              tolerance_deg: float = 1.0) -> List[BrepEdge]:
    """Select linear edges parallel to a given direction.

    Args:
        brep_solid: A BrepSolid object
        direction: Direction vector as (x, y, z) tuple or list
        tolerance_deg: Angular tolerance in degrees (default 1.0)

    Returns:
        List of BrepEdge objects parallel to the direction
    """
    require_occ()
    _ensure_occ_imports()

    # Normalize the input direction
    dx, dy, dz = direction[0], direction[1], direction[2]
    length = math.sqrt(dx*dx + dy*dy + dz*dz)
    if length < 1e-12:
        return []
    target_dir = (dx/length, dy/length, dz/length)

    tolerance_rad = math.radians(tolerance_deg)
    cos_tol = math.cos(tolerance_rad)

    result = []
    raw_edges = _get_unique_edges(brep_solid.shape)

    for edge in raw_edges:
        edge_dir = _get_edge_direction(edge)
        if edge_dir is None:
            continue  # Skip curved edges

        # Calculate absolute dot product (parallel edges can be in either direction)
        dot = abs(edge_dir[0] * target_dir[0] +
                  edge_dir[1] * target_dir[1] +
                  edge_dir[2] * target_dir[2])

        if dot >= cos_tol:
            result.append(BrepEdge(edge))

    return result


def select_edges_by_length(brep_solid: BrepSolid,
                           min_length: Optional[float] = None,
                           max_length: Optional[float] = None) -> List[BrepEdge]:
    """Select edges within a length range.

    Args:
        brep_solid: A BrepSolid object
        min_length: Minimum edge length (inclusive), None for no minimum
        max_length: Maximum edge length (inclusive), None for no maximum

    Returns:
        List of BrepEdge objects within the length range
    """
    require_occ()
    _ensure_occ_imports()

    result = []
    raw_edges = _get_unique_edges(brep_solid.shape)

    for edge in raw_edges:
        length = _get_edge_length(edge)

        if min_length is not None and length < min_length:
            continue
        if max_length is not None and length > max_length:
            continue

        result.append(BrepEdge(edge))

    return result


def select_edges_at_z(brep_solid: BrepSolid,
                      z_value: float,
                      tolerance: float = 0.001) -> List[BrepEdge]:
    """Select edges that lie entirely at a specific Z height.

    This selects horizontal edges where both endpoints are at the given Z value.

    Args:
        brep_solid: A BrepSolid object
        z_value: The Z coordinate to match
        tolerance: Position tolerance (default 0.001)

    Returns:
        List of BrepEdge objects at the specified Z height
    """
    require_occ()
    _ensure_occ_imports()

    result = []
    raw_edges = _get_unique_edges(brep_solid.shape)

    for edge in raw_edges:
        endpoints = _get_edge_endpoints(edge)
        if endpoints is None:
            continue

        p1, p2 = endpoints

        # Check if both endpoints are at the target Z
        if (abs(p1[2] - z_value) <= tolerance and
            abs(p2[2] - z_value) <= tolerance):
            result.append(BrepEdge(edge))

    return result


def select_edges_in_z_range(brep_solid: BrepSolid,
                            z_min: float,
                            z_max: float,
                            tolerance: float = 0.001) -> List[BrepEdge]:
    """Select edges where both endpoints are within a Z range.

    Args:
        brep_solid: A BrepSolid object
        z_min: Minimum Z coordinate
        z_max: Maximum Z coordinate
        tolerance: Position tolerance (default 0.001)

    Returns:
        List of BrepEdge objects within the Z range
    """
    require_occ()
    _ensure_occ_imports()

    result = []
    raw_edges = _get_unique_edges(brep_solid.shape)

    for edge in raw_edges:
        endpoints = _get_edge_endpoints(edge)
        if endpoints is None:
            continue

        p1, p2 = endpoints

        # Check if both endpoints are within the Z range
        if (z_min - tolerance <= p1[2] <= z_max + tolerance and
            z_min - tolerance <= p2[2] <= z_max + tolerance):
            result.append(BrepEdge(edge))

    return result


def select_edges_crossing_z(brep_solid: BrepSolid,
                            z_value: float,
                            tolerance: float = 0.001) -> List[BrepEdge]:
    """Select edges that cross a specific Z height (vertical edges spanning Z).

    This is useful for selecting vertical edges of a pocket or boss.

    Args:
        brep_solid: A BrepSolid object
        z_value: The Z coordinate that edges must span
        tolerance: Position tolerance (default 0.001)

    Returns:
        List of BrepEdge objects that cross the Z value
    """
    require_occ()
    _ensure_occ_imports()

    result = []
    raw_edges = _get_unique_edges(brep_solid.shape)

    for edge in raw_edges:
        endpoints = _get_edge_endpoints(edge)
        if endpoints is None:
            continue

        p1, p2 = endpoints
        z1, z2 = p1[2], p2[2]

        # Check if the edge spans the Z value
        z_low = min(z1, z2) - tolerance
        z_high = max(z1, z2) + tolerance

        if z_low < z_value < z_high:
            result.append(BrepEdge(edge))

    return result


def select_top_edges(brep_solid: BrepSolid,
                     tolerance: float = 0.001) -> List[BrepEdge]:
    """Select edges at the maximum Z height of the solid.

    Args:
        brep_solid: A BrepSolid object
        tolerance: Position tolerance (default 0.001)

    Returns:
        List of BrepEdge objects at the top of the solid
    """
    require_occ()
    _ensure_occ_imports()

    # First, find the maximum Z coordinate
    max_z = None
    raw_edges = _get_unique_edges(brep_solid.shape)

    for edge in raw_edges:
        endpoints = _get_edge_endpoints(edge)
        if endpoints is None:
            continue
        p1, p2 = endpoints
        edge_max_z = max(p1[2], p2[2])
        if max_z is None or edge_max_z > max_z:
            max_z = edge_max_z

    if max_z is None:
        return []

    return select_edges_at_z(brep_solid, max_z, tolerance)


def select_bottom_edges(brep_solid: BrepSolid,
                        tolerance: float = 0.001) -> List[BrepEdge]:
    """Select edges at the minimum Z height of the solid.

    Args:
        brep_solid: A BrepSolid object
        tolerance: Position tolerance (default 0.001)

    Returns:
        List of BrepEdge objects at the bottom of the solid
    """
    require_occ()
    _ensure_occ_imports()

    # First, find the minimum Z coordinate
    min_z = None
    raw_edges = _get_unique_edges(brep_solid.shape)

    for edge in raw_edges:
        endpoints = _get_edge_endpoints(edge)
        if endpoints is None:
            continue
        p1, p2 = endpoints
        edge_min_z = min(p1[2], p2[2])
        if min_z is None or edge_min_z < min_z:
            min_z = edge_min_z

    if min_z is None:
        return []

    return select_edges_at_z(brep_solid, min_z, tolerance)


def select_edges_near_point(brep_solid: BrepSolid,
                            target_point: Union[Tuple[float, float, float], list],
                            max_distance: float) -> List[BrepEdge]:
    """Select edges whose midpoint is within a distance of a target point.

    Args:
        brep_solid: A BrepSolid object
        target_point: The reference point as (x, y, z)
        max_distance: Maximum distance from target point

    Returns:
        List of BrepEdge objects near the target point
    """
    require_occ()
    _ensure_occ_imports()

    tx, ty, tz = target_point[0], target_point[1], target_point[2]

    result = []
    raw_edges = _get_unique_edges(brep_solid.shape)

    for edge in raw_edges:
        midpoint = _get_edge_midpoint(edge)
        if midpoint is None:
            continue

        mx, my, mz = midpoint
        dist = math.sqrt((mx-tx)**2 + (my-ty)**2 + (mz-tz)**2)

        if dist <= max_distance:
            result.append(BrepEdge(edge))

    return result


def select_edges_in_cylinder(brep_solid: BrepSolid,
                             center: Union[Tuple[float, float, float], list],
                             radius: float,
                             axis: Union[Tuple[float, float, float], list] = (0, 0, 1)) -> List[BrepEdge]:
    """Select edges whose midpoint is within a cylindrical region.

    Useful for selecting edges around a hole or boss.

    Args:
        brep_solid: A BrepSolid object
        center: Center point of cylinder axis (x, y, z)
        radius: Radius of the cylindrical region
        axis: Cylinder axis direction (default Z axis)

    Returns:
        List of BrepEdge objects within the cylinder
    """
    require_occ()
    _ensure_occ_imports()

    cx, cy, cz = center[0], center[1], center[2]

    # Normalize axis
    ax, ay, az = axis[0], axis[1], axis[2]
    alen = math.sqrt(ax*ax + ay*ay + az*az)
    if alen < 1e-12:
        return []
    ax, ay, az = ax/alen, ay/alen, az/alen

    result = []
    raw_edges = _get_unique_edges(brep_solid.shape)

    for edge in raw_edges:
        midpoint = _get_edge_midpoint(edge)
        if midpoint is None:
            continue

        mx, my, mz = midpoint

        # Vector from center to midpoint
        vx, vy, vz = mx - cx, my - cy, mz - cz

        # Project onto axis
        proj = vx*ax + vy*ay + vz*az

        # Perpendicular component
        px = vx - proj * ax
        py = vy - proj * ay
        pz = vz - proj * az

        # Distance from axis
        dist = math.sqrt(px*px + py*py + pz*pz)

        if dist <= radius:
            result.append(BrepEdge(edge))

    return result


def filter_curved_edges(edges: List[BrepEdge]) -> List[BrepEdge]:
    """Filter to keep only curved (non-linear) edges.

    Args:
        edges: List of BrepEdge objects to filter

    Returns:
        List of curved BrepEdge objects
    """
    _ensure_occ_imports()

    result = []
    for brep_edge in edges:
        direction = _get_edge_direction(brep_edge.shape)
        if direction is None:  # Curved edge
            result.append(brep_edge)

    return result


def filter_linear_edges(edges: List[BrepEdge]) -> List[BrepEdge]:
    """Filter to keep only linear (straight) edges.

    Args:
        edges: List of BrepEdge objects to filter

    Returns:
        List of linear BrepEdge objects
    """
    _ensure_occ_imports()

    result = []
    for brep_edge in edges:
        direction = _get_edge_direction(brep_edge.shape)
        if direction is not None:  # Linear edge
            result.append(brep_edge)

    return result


def edge_info(edge: BrepEdge) -> dict:
    """Get detailed information about an edge.

    Args:
        edge: A BrepEdge object

    Returns:
        Dictionary with edge properties:
            - length: Edge length
            - endpoints: ((x1, y1, z1), (x2, y2, z2))
            - midpoint: (x, y, z)
            - direction: (dx, dy, dz) for linear edges, None for curved
            - is_linear: True if edge is straight
            - is_vertical: True if parallel to Z axis
            - is_horizontal: True if perpendicular to Z axis
    """
    _ensure_occ_imports()

    raw_edge = edge.shape if isinstance(edge, BrepEdge) else edge

    length = _get_edge_length(raw_edge)
    endpoints = _get_edge_endpoints(raw_edge)
    midpoint = _get_edge_midpoint(raw_edge)
    direction = _get_edge_direction(raw_edge)

    is_linear = direction is not None
    is_vertical = False
    is_horizontal = False

    if direction is not None:
        # Check if vertical (parallel to Z)
        z_dot = abs(direction[2])
        if z_dot >= 0.9998:  # ~1 degree tolerance
            is_vertical = True
        # Check if horizontal (perpendicular to Z)
        if z_dot <= 0.0175:  # ~1 degree tolerance
            is_horizontal = True

    return {
        'length': length,
        'endpoints': endpoints,
        'midpoint': midpoint,
        'direction': direction,
        'is_linear': is_linear,
        'is_vertical': is_vertical,
        'is_horizontal': is_horizontal,
    }


# Convenience functions for combining selections

def union_edges(*edge_lists: List[BrepEdge]) -> List[BrepEdge]:
    """Combine multiple edge lists, removing duplicates.

    Args:
        *edge_lists: Variable number of BrepEdge lists

    Returns:
        Combined list of unique BrepEdge objects
    """
    seen = set()
    result = []

    for edge_list in edge_lists:
        for edge in edge_list:
            # Use the shape's hash for uniqueness
            edge_hash = hash(edge.shape.HashCode(2147483647))
            if edge_hash not in seen:
                seen.add(edge_hash)
                result.append(edge)

    return result


def intersect_edges(*edge_lists: List[BrepEdge]) -> List[BrepEdge]:
    """Find edges common to all input lists.

    Args:
        *edge_lists: Variable number of BrepEdge lists

    Returns:
        List of BrepEdge objects present in all input lists
    """
    if not edge_lists:
        return []

    # Build hash sets for each list
    hash_sets = []
    edge_by_hash = {}

    for edge_list in edge_lists:
        hash_set = set()
        for edge in edge_list:
            edge_hash = hash(edge.shape.HashCode(2147483647))
            hash_set.add(edge_hash)
            edge_by_hash[edge_hash] = edge
        hash_sets.append(hash_set)

    # Find intersection
    common_hashes = hash_sets[0]
    for hash_set in hash_sets[1:]:
        common_hashes = common_hashes & hash_set

    return [edge_by_hash[h] for h in common_hashes]


def subtract_edges(base_edges: List[BrepEdge],
                   edges_to_remove: List[BrepEdge]) -> List[BrepEdge]:
    """Remove edges from a list.

    Args:
        base_edges: The original list of edges
        edges_to_remove: Edges to remove from the base list

    Returns:
        List of BrepEdge objects from base_edges not in edges_to_remove
    """
    remove_hashes = set()
    for edge in edges_to_remove:
        remove_hashes.add(hash(edge.shape.HashCode(2147483647)))

    return [e for e in base_edges
            if hash(e.shape.HashCode(2147483647)) not in remove_hashes]


__all__ = [
    # Core selection functions
    'get_all_edges',
    'select_vertical_edges',
    'select_horizontal_edges',
    'select_edges_by_direction',
    'select_edges_by_length',
    'select_edges_at_z',
    'select_edges_in_z_range',
    'select_edges_crossing_z',
    'select_top_edges',
    'select_bottom_edges',
    'select_edges_near_point',
    'select_edges_in_cylinder',
    # Filter functions
    'filter_curved_edges',
    'filter_linear_edges',
    # Utility functions
    'edge_info',
    'union_edges',
    'intersect_edges',
    'subtract_edges',
]
