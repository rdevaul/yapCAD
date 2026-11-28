"""``sdf.repeat``: copies of a field on a rectilinear lattice.

The operator folds a query point into the cell it falls in, so its cost does
not grow with the number of copies.  What has to be right is everything
folding alone gets wrong: copies that overlap their neighbours' cells, parts
that do not sit at the centre of their cell, and patterns that stop.  Each
case is checked against the thing ``repeat`` stands for -- the union of the
translated copies, built one by one.
"""

import itertools
import math

import numpy as np
import pytest

from yapcad import sdf
from yapcad.sdf.node import INF, SdfError


def brute(child, period, count, points, span=7):
    """The union of the translated copies, evaluated one copy at a time."""
    ranges = []
    for p, n in zip(period, count):
        if p == 0:
            ranges.append([0])
        elif n:
            ranges.append(range(n))
        else:
            ranges.append(range(-span, span + 1))
    out = np.full(len(points), np.inf)
    for k in itertools.product(*ranges):
        shift = [k[i] * period[i] for i in range(3)]
        out = np.minimum(out, sdf.evaluate(sdf.translate(child, shift),
                                           points))
    return out


CASES = {
    "fits its cell": (sdf.sphere(1.0), (3, 3, 0), (0, 0, 0)),
    "three axes": (sdf.sphere(1.0), (3, 4, 5), (0, 0, 0)),
    "overlaps its neighbours": (
        sdf.capsule((-2.5, -1, 0), (2.5, 1.5, 0), 0.6), (3, 2, 0), (0, 0, 0)),
    "off the cell centre": (
        sdf.translate(sdf.box((1, 1, 1)), (1.3, 0.2, 0)), (3, 3, 0),
        (0, 0, 0)),
    "a counted row": (sdf.sphere(1.0), (2.5, 0, 0), (4, 0, 0)),
    "a counted, overlapping sheet": (
        sdf.capsule((-2, 0, 0), (2, 1, 0), 0.5), (3, 2, 0), (3, 4, 0)),
    "several periods wide": (sdf.box((7, 1, 1)), (2, 3, 0), (0, 0, 0)),
}


@pytest.mark.parametrize("name", CASES)
def test_repeat_is_the_union_of_its_copies(name):
    child, period, count = CASES[name]
    node = sdf.repeat(child, period, count)
    rng = np.random.default_rng(1)
    points = rng.uniform(-9.0, 9.0, (20000, 3))
    if period[2] == 0:
        points[:, 2] *= 0.2
    got = sdf.evaluate(node, points)
    want = brute(child, period, count, points)
    # The same solid: inside and outside agree everywhere.
    assert (np.sign(got) == np.sign(want)).all()
    # Never further from the surface than it really is, so the mesher and
    # the ray tracer can step by it.
    assert (got - want).max() <= 1e-12
    # And the true distance close to the surface, where the mesh is placed.
    near = np.abs(want) < 0.3
    assert near.any()
    assert np.abs(got - want)[near].max() < 1e-9


def test_an_endless_axis_is_unbounded_and_a_counted_one_is_not():
    ball = sdf.sphere(1.0)
    sheet = sdf.repeat(ball, (3, 4, 0))
    assert sheet.bounds == ((-INF, -INF, -1.0), (INF, INF, 1.0))
    grid = sdf.repeat(ball, (3, 4, 0), (3, 2, 0))
    assert grid.bounds == ((-1.0, -1.0, -1.0), (7.0, 5.0, 1.0))
    row = sdf.repeat(ball, (3, 0, 0), 5)
    assert row.bounds == ((-1.0, -1.0, -1.0), (13.0, 1.0, 1.0))


def test_a_scalar_period_repeats_on_all_three_axes():
    node = sdf.repeat(sdf.sphere(0.5), 2.0)
    assert tuple(node.p["period"]) == (2.0, 2.0, 2.0)
    centres = [[2.0, 4.0, -6.0], [0.0, 0.0, 0.0], [-2.0, 2.0, 2.0]]
    assert sdf.evaluate(node, centres) == pytest.approx([-0.5] * 3)
    between = [[1.0, 1.0, 1.0]]
    assert sdf.evaluate(node, between)[0] > 0.0


def test_the_first_copy_is_the_part_where_it_stands():
    """A counted pattern grows along +x, +y, +z from the part, so moving the
    part moves the pattern and nothing appears behind it."""
    part = sdf.translate(sdf.sphere(0.5), (10.0, 0.0, 0.0))
    row = sdf.repeat(part, (2.0, 0, 0), 3)
    inside = sdf.evaluate(row, [[10, 0, 0], [12, 0, 0], [14, 0, 0]])
    outside = sdf.evaluate(row, [[8, 0, 0], [16, 0, 0]])
    assert (inside < 0.0).all() and (outside > 0.0).all()


def test_repeat_keeps_the_childs_lipschitz_bound():
    twisted = sdf.twist(sdf.box((2.0, 0.5, 2.0)), 0.6)
    assert sdf.repeat(twisted, (5, 5, 0)).lipschitz == twisted.lipschitz
    assert not sdf.repeat(sdf.sphere(1.0), 3.0).exact


def test_an_endless_pattern_meshes_once_it_is_cut_to_a_part():
    from yapcad.geom3d import issolidclosed, volumeof
    balls = sdf.intersect(sdf.box((9.0, 9.0, 9.0)),
                          sdf.repeat(sdf.sphere(1.0), 3.0))
    solid = sdf.to_solid(balls, cell_mm=0.1)
    assert issolidclosed(solid)
    assert volumeof(solid) == pytest.approx(27 * 4.0 / 3.0 * math.pi,
                                            rel=1e-2)


def test_a_repeat_round_trips_through_its_tree():
    node = sdf.repeat(sdf.sphere(1.0), (3, 0, 2), (0, 0, 4))
    again = sdf.tree_from_json(sdf.tree_to_json(node))
    assert again == node and again.digest == node.digest
    assert "repeat" in sdf.describe_blockers(node)


@pytest.mark.parametrize("period,count,message", [
    ((0, 0, 0), None, "at least one period"),
    ((-1, 0, 0), None, "not negative"),
    ((1, 2), None, "a number or three"),
    ((1, 1, 1), (-1, 0, 0), "counts must not be negative"),
    ("wide", None, "a number or three"),
])
def test_bad_arguments_are_refused(period, count, message):
    with pytest.raises(SdfError, match=message):
        sdf.repeat(sdf.sphere(1.0), period, count)


def test_a_part_unbounded_on_the_repeated_axis_is_refused():
    slab = sdf.half_space((1.0, 0.0, 0.0))
    with pytest.raises(SdfError, match="unbounded along x"):
        sdf.repeat(slab, (3, 0, 0)).bounds
    # Repeating it along an axis it is bounded on is fine.
    plate = sdf.intersect(sdf.half_space((0, 0, 1.0)),
                          sdf.half_space((0, 0, -1.0), 1.0))
    with pytest.raises(SdfError, match="unbounded along x"):
        sdf.repeat(plate, (3, 0, 0)).bounds


@pytest.mark.parametrize("name", CASES)
def test_repeat_is_continuous_where_the_copies_tested_change(name):
    """Crossing from one cell to the next changes which copies a point is
    tested against.  The field must not jump when it does -- at the ends of
    a counted pattern in particular, where there is no copy beyond to take
    over -- or it stops being Lipschitz and the ray tracer's step rule with
    it."""
    child, period, count = CASES[name]
    node = sdf.repeat(child, period, count)
    rng = np.random.default_rng(2)
    a = rng.uniform(-12.0, 12.0, (20000, 3))
    step = rng.normal(size=(20000, 3))
    step *= rng.uniform(1e-4, 0.4, (20000, 1)) / np.linalg.norm(
        step, axis=1, keepdims=True)
    fa, fb = sdf.evaluate(node, a), sdf.evaluate(node, a + step)
    moved = np.linalg.norm(step, axis=1)
    assert (np.abs(fa - fb) <= node.lipschitz * moved + 1e-9).all()


# -- repeat_polar ------------------------------------------------------------

def brute_polar(child, count, points):
    out = np.full(len(points), np.inf)
    for k in range(count):
        turned = sdf.rotate(child, (0.0, 0.0, 1.0), 360.0 * k / count)
        out = np.minimum(out, sdf.evaluate(turned, points))
    return out


POLAR = {
    "a spoke": (sdf.capsule((3, 0, 0), (9, 0, 0), 0.6), 7),
    "a fin wider than its sector": (
        sdf.translate(sdf.box((4, 9, 2)), (6, 1.5, 0)), 5),
    "a ball off the x axis": (
        sdf.translate(sdf.sphere(1.0), (2.5, 4.0, 0.3)), 12),
    "a part that reaches the axis": (sdf.box((8, 2, 2)), 6),
    "two copies": (sdf.translate(sdf.sphere(2.0), (5, 0, 0)), 2),
    "many, overlapping": (sdf.translate(sdf.sphere(1.5), (4, 0, 0)), 24),
}


@pytest.mark.parametrize("name", POLAR)
def test_repeat_polar_is_the_union_of_its_copies(name):
    child, count = POLAR[name]
    node = sdf.repeat_polar(child, count)
    rng = np.random.default_rng(3)
    points = rng.uniform(-11.0, 11.0, (20000, 3))
    points[:, 2] *= 0.3
    got = sdf.evaluate(node, points)
    want = brute_polar(child, count, points)
    assert (np.sign(got) == np.sign(want)).all()
    assert (got - want).max() <= 1e-9
    near = np.abs(want) < 0.3
    assert near.any()
    assert np.abs(got - want)[near].max() < 1e-9
    # Continuous as a point passes from one sector to the next.
    step = rng.normal(size=points.shape)
    step *= rng.uniform(1e-4, 0.3, (len(points), 1)) / np.linalg.norm(
        step, axis=1, keepdims=True)
    moved = np.linalg.norm(step, axis=1)
    jump = np.abs(got - sdf.evaluate(node, points + step))
    assert (jump <= node.lipschitz * moved + 1e-9).all()


def test_repeat_polar_is_bounded_by_the_circle_its_part_sweeps():
    fin = sdf.translate(sdf.box((4.0, 2.0, 6.0)), (8.0, 0.0, 1.0))
    ring = sdf.repeat_polar(fin, 3)
    r = math.hypot(10.0, 1.0)
    (x0, y0, z0), (x1, y1, z1) = ring.bounds
    assert (x0, y0, x1, y1) == pytest.approx((-r, -r, r, r))
    assert (z0, z1) == (-2.0, 4.0)


def test_one_copy_is_the_part_itself():
    ball = sdf.sphere(1.0)
    assert sdf.repeat_polar(ball, 1) is ball


@pytest.mark.parametrize("count", [0, -3, 2.5, "many"])
def test_repeat_polar_refuses_a_count_that_is_not_a_whole_number(count):
    with pytest.raises(SdfError, match="whole number"):
        sdf.repeat_polar(sdf.sphere(1.0), count)


def test_repeat_polar_round_trips_through_its_tree():
    node = sdf.repeat_polar(sdf.translate(sdf.sphere(1.0), (4, 0, 0)), 9)
    again = sdf.tree_from_json(sdf.tree_to_json(node))
    assert again == node and again.digest == node.digest


def test_repeat_polar_knows_the_axis_is_clear_of_a_ring_of_parts():
    """Every sector's wedge meets at the axis, so a bound from the wedges
    alone says nothing there.  The parts' own distance from the axis does:
    on the axis of a ring of rivets the field is the clearance, not zero."""
    ring = sdf.repeat_polar(sdf.translate(sdf.sphere(0.5), (10.0, 0.0, 0.0)),
                            24)
    on_axis = sdf.evaluate(ring, [[0.0, 0.0, 0.0], [0.0, 0.0, 3.0]])
    assert on_axis[0] == pytest.approx(9.5)
    assert on_axis[1] > 9.0
