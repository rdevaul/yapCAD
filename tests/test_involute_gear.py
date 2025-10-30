import math

from examples.involute_gear_package.involute_gear import generate_involute_profile


def _segments(points):
    return list(zip(points[:-1], points[1:]))


def _near(pt_a, pt_b, tol=1e-9):
    return abs(pt_a[0] - pt_b[0]) <= tol and abs(pt_a[1] - pt_b[1]) <= tol


def _segments_intersect(seg_a, seg_b, tol=1e-9):
    (ax0, ay0), (ax1, ay1) = seg_a
    (bx0, by0), (bx1, by1) = seg_b

    def orient(px0, py0, px1, py1, px2, py2):
        return (px1 - px0) * (py2 - py0) - (py1 - py0) * (px2 - px0)

    def on_segment(px0, py0, px1, py1, qx, qy):
        return (min(px0, px1) - tol <= qx <= max(px0, px1) + tol and
                min(py0, py1) - tol <= qy <= max(py0, py1) + tol)

    o1 = orient(ax0, ay0, ax1, ay1, bx0, by0)
    o2 = orient(ax0, ay0, ax1, ay1, bx1, by1)
    o3 = orient(bx0, by0, bx1, by1, ax0, ay0)
    o4 = orient(bx0, by0, bx1, by1, ax1, ay1)

    if o1 * o2 < -tol and o3 * o4 < -tol:
        return True

    if abs(o1) <= tol and on_segment(ax0, ay0, ax1, ay1, bx0, by0):
        return True
    if abs(o2) <= tol and on_segment(ax0, ay0, ax1, ay1, bx1, by1):
        return True
    if abs(o3) <= tol and on_segment(bx0, by0, bx1, by1, ax0, ay0):
        return True
    if abs(o4) <= tol and on_segment(bx0, by0, bx1, by1, ax1, ay1):
        return True
    return False


def _has_real_self_intersections(points):
    segs = _segments(points)
    last = len(segs) - 1

    for i, seg_a in enumerate(segs):
        for j, seg_b in enumerate(segs):
            if j <= i:
                continue
            if abs(i - j) <= 1:
                continue
            if (i == 0 and j == last) or (j == 0 and i == last):
                continue
            if (_near(seg_a[0], seg_b[0]) or _near(seg_a[0], seg_b[1]) or
                    _near(seg_a[1], seg_b[0]) or _near(seg_a[1], seg_b[1])):
                continue
            if _segments_intersect(seg_a, seg_b):
                return True
    return False


def _signed_area(points):
    total = 0.0
    for (x0, y0), (x1, y1) in _segments(points):
        total += x0 * y1 - x1 * y0
    return total / 2.0


def test_involute_profile_orientation_and_simple_polygon():
    outline = generate_involute_profile(teeth=18, module_mm=2.0)
    assert outline[0] == outline[-1]
    area = _signed_area(outline)
    assert area > 0, "outer loop must follow right-hand rule"
    assert not _has_real_self_intersections(outline), "profile should be simple"


def test_tip_arc_traversal_is_short_path():
    outline = generate_involute_profile(teeth=18, module_mm=2.0)
    # Track the maximum angular difference between successive vertices; tip arcs
    # should not span almost an entire revolution.
    max_delta = 0.0
    prev_angle = None
    for x, y in outline:
        angle = math.atan2(y, x)
        if prev_angle is not None:
            diff = abs(angle - prev_angle)
            diff = min(diff, abs(diff - 2 * math.pi))
            max_delta = max(max_delta, diff)
        prev_angle = angle
    # With the corrected ordering, the outline never jumps across the circle.
    assert max_delta < math.pi / 2, "adjacent segments should be localised on the circle"
