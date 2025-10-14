"""Native boolean engine extracted from yapcad.geom3d."""

from __future__ import annotations

import math

from yapcad.geom import *
from yapcad.geom_util import *
from yapcad.xform import *
from yapcad.triangulator import triangulate_polygon
from yapcad.octtree import NTree


def _geom3d():
    from yapcad import geom3d as _g3
    return _g3

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
    box = _geom3d().solidbbox(target)
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
    p0, n = _geom3d().tri2p0n(tri)
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

    _, _, to_local, to_world = _geom3d().tri2p0n(tri, basis=True)

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
        # Don't call _orient_triangle! The winding order is already preserved through:
        # 1. Original triangle had correct outward normal
        # 2. Clipping preserves vertex order (just removes parts)
        # 3. Transform to 2D preserves order
        # 4. triangulate_polygon produces CCW triangles in 2D
        # 5. Transform back to 3D should maintain correct orientation
        # Check for degeneracy with better quality metric:
        v01 = sub(tri3d[1], tri3d[0])
        v02 = sub(tri3d[2], tri3d[0])
        cross_prod = cross(v01, v02)
        area = mag(cross_prod)
        if area >= tol:
            # Also filter out very thin slivers by checking aspect ratio
            # A degenerate sliver has tiny area relative to edge lengths
            edge_lengths = [mag(v01), mag(v02), dist(tri3d[1], tri3d[2])]
            max_edge = max(edge_lengths)
            # For a reasonable triangle, area should be at least 1% of max_edge^2
            # This filters out slivers where one edge is extremely small
            quality_threshold = 0.01 * max_edge * max_edge
            if area >= max(tol, quality_threshold):
                outside_triangles.append([point(tri3d[0]), point(tri3d[1]), point(tri3d[2])])

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

    if not _geom3d().issolid(sld, fast=False):
        raise ValueError('invalid solid passed to stitch_solid')

    triangles = list(_iter_triangles_from_solid(sld))
    stitched = stitch_open_edges(triangles, tol)
    if not stitched:
        return sld

    surface_result = _surface_from_triangles(stitched)
    if surface_result:
        return _geom3d().solid([surface_result], [], ['utility', 'stitch_open_edges'])
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
    _, n = _geom3d().tri2p0n(tri)
    return n
def solid_contains_point(sld, p, tol=_DEFAULT_RAY_TOL):
    if not _geom3d().issolid(sld, fast=False):
        raise ValueError('invalid solid passed to solid_contains_point')
    if not _geom3d().ispoint(p):
        raise ValueError('invalid point passed to solid_contains_point')

    surfaces = sld[1]
    if not surfaces:
        return False

    box = _geom3d().solidbbox(sld)
    if box:
        expanded = [
            point(box[0][0] - tol, box[0][1] - tol, box[0][2] - tol),
            point(box[1][0] + tol, box[1][1] + tol, box[1][2] + tol),
        ]
        if not _geom3d().isinsidebbox(expanded, p):
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
    if not _geom3d().issolid(a, fast=False) or not _geom3d().issolid(b, fast=False):
        raise ValueError('invalid solid passed to solids_intersect')
    if not a[1] or not b[1]:
        return False

    box_a = _geom3d().solidbbox(a)
    box_b = _geom3d().solidbbox(b)
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
            if _geom3d().triTriIntersect(tri_a, tri_b):
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

        # Filter out triangles in the interior of the overlap region
        # For union, if a triangle's center is inside both input solids,
        # it's in the interior and should not be on the surface
        filtered_tris = []
        for tri in result_tris:
            center = point((tri[0][0] + tri[1][0] + tri[2][0]) / 3.0,
                          (tri[0][1] + tri[1][1] + tri[2][1]) / 3.0,
                          (tri[0][2] + tri[1][2] + tri[2][2]) / 3.0)
            # Use a tighter tolerance for containment check to avoid false positives
            check_tol = tol * 100
            in_a = solid_contains_point(a, center, tol=check_tol)
            in_b = solid_contains_point(b, center, tol=check_tol)
            # Keep triangle only if it's not clearly inside both solids
            if not (in_a and in_b):
                filtered_tris.append(tri)
        result_tris = filtered_tris

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
        return _geom3d().solid([surface_result], [], ['boolean', operation])
    return _geom3d().solid([], [], ['boolean', operation])



