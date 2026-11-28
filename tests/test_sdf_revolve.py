"""``sdf.revolve``: a planar profile swept round the z axis.

The profile is read with x as the distance from the axis and y as the
height, so the checks are against solids of revolution built another way.
"""

import math

import numpy as np
import pytest

from yapcad import sdf
from yapcad.geom3d import issolidclosed, volumeof
from yapcad.sdf.node import SdfError


def cloud(reach=7.0, count=20000, seed=5):
    return np.random.default_rng(seed).uniform(-reach, reach, (count, 3))


def test_a_rectangle_off_the_axis_revolves_to_a_tube():
    tube = sdf.revolve(sdf.polygon([(2, 0), (5, 0), (5, 3), (2, 3)]))
    same = sdf.translate(sdf.subtract(sdf.cylinder(5.0, 3.0),
                                      sdf.cylinder(2.0, 9.0)),
                         (0.0, 0.0, 1.5))
    points = cloud()
    got, want = sdf.evaluate(tube, points), sdf.evaluate(same, points)
    assert (np.sign(got) == np.sign(want)).all()
    assert tube.bounds == ((-5.0, -5.0, 0.0), (5.0, 5.0, 3.0))
    # A planar distance stays a distance when swept round an axis the
    # profile does not cross.
    assert tube.exact
    near = np.abs(got) < 0.5
    radial = np.hypot(points[:, 0], points[:, 1])
    box = np.maximum(np.abs(radial - 3.5) - 1.5, np.abs(points[:, 2] - 1.5) - 1.5)
    outside = np.hypot(np.maximum(np.abs(radial - 3.5) - 1.5, 0.0),
                       np.maximum(np.abs(points[:, 2] - 1.5) - 1.5, 0.0))
    exact = np.where(box > 0.0, outside, box)
    assert got[near] == pytest.approx(exact[near], abs=1e-9)


def test_a_profile_carried_across_the_axis_closes_on_it():
    """A frustum: its volume, and the axis inside it."""
    frustum = sdf.revolve(sdf.polygon([(-1, 0), (4, 0), (1, 6), (-1, 6)]))
    solid = sdf.to_solid(frustum, cell_mm=0.1)
    assert issolidclosed(solid)
    assert volumeof(solid) == pytest.approx(math.pi * 6 * (16 + 4 + 1) / 3,
                                            rel=5e-3)
    on_axis = sdf.evaluate(frustum, [[0.0, 0.0, 3.0]])
    assert on_axis[0] < -0.5
    assert not frustum.exact


def test_only_the_positive_half_of_the_profile_is_revolved():
    both = sdf.revolve(sdf.union(
        sdf.polygon([(2, 0), (3, 0), (3, 1), (2, 1)]),
        sdf.polygon([(-9, 0), (-8, 0), (-8, 1), (-9, 1)])))
    points = [[8.5, 0.0, 0.5], [2.5, 0.0, 0.5], [0.0, -2.5, 0.5]]
    inside = sdf.evaluate(both, points) < 0.0
    assert inside.tolist() == [False, True, True]


def test_a_profile_wholly_behind_the_axis_is_refused():
    behind = sdf.revolve(sdf.polygon([(-5, 0), (-2, 0), (-2, 3)]))
    with pytest.raises(SdfError, match="positive x"):
        behind.bounds


def test_revolve_keeps_the_profiles_lipschitz_bound_and_round_trips():
    profile = sdf.polygon([(1, 0), (4, 0), (4, 2), (1, 2)])
    node = sdf.revolve(profile)
    assert node.lipschitz == profile.lipschitz
    again = sdf.tree_from_json(sdf.tree_to_json(node))
    assert again == node and again.digest == node.digest
