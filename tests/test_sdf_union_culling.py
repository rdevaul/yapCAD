"""A union with many operands skips children that cannot change the result.

Evaluating every child at every point is what made field-built engraved
lettering slow: a union of a couple of hundred stroke boxes paid for all of
them at every sample, while any one point is near only a few.  The saving must
not change a single value, so these tests compare the culled path against the
plain fold rather than against recorded numbers.
"""
import numpy as np
import pytest

from yapcad import sdf
from yapcad.sdf import ops
from yapcad.sdf.evaluate import evaluate


@pytest.fixture
def no_culling(monkeypatch):
    """Force the plain fold over every child."""
    monkeypatch.setattr(ops, "_UNION_CULL_MIN_CHILDREN", 10 ** 9)


def many_boxes(count=40, pitch=3.0):
    return sdf.union(*[
        sdf.translate(sdf.box((1.0, 2.0, 0.5)), (i * pitch, 0.0, 0.0))
        for i in range(count)
    ])


def sample(seed=0, n=4000, span=((-5.0, 125.0), (-8.0, 8.0), (-4.0, 4.0))):
    rng = np.random.default_rng(seed)
    return np.column_stack([rng.uniform(lo, hi, n) for lo, hi in span])


def test_culling_changes_no_value(monkeypatch):
    node, pts = many_boxes(), sample()
    culled = evaluate(node, pts)
    monkeypatch.setattr(ops, "_UNION_CULL_MIN_CHILDREN", 10 ** 9)
    plain = evaluate(node, pts)
    # min() is exact and order-independent, so this is equality, not closeness.
    assert np.array_equal(culled, plain)


def test_culling_changes_no_value_under_a_subtraction(monkeypatch):
    # The shape engraving takes: a slab with many pockets cut out of it.
    node = sdf.subtract(sdf.box((140.0, 16.0, 4.0)), many_boxes())
    pts = sample(seed=1)
    culled = evaluate(node, pts)
    monkeypatch.setattr(ops, "_UNION_CULL_MIN_CHILDREN", 10 ** 9)
    assert np.array_equal(culled, evaluate(node, pts))


@pytest.mark.parametrize("split", (10 ** 9, 1024, 64))
def test_the_split_threshold_changes_no_value(monkeypatch, split):
    node, pts = many_boxes(), sample(seed=2)
    reference = evaluate(node, pts)
    monkeypatch.setattr(ops, "_UNION_SPLIT_MIN_POINTS", split)
    assert np.array_equal(evaluate(node, pts), reference)


def test_a_small_union_still_takes_the_plain_path():
    # Below the threshold the bookkeeping costs more than it saves, so the
    # fold must stay in place -- and still agree.
    node = sdf.union(sdf.box((1.0, 1.0, 1.0)),
                     sdf.translate(sdf.box((1.0, 1.0, 1.0)), (4.0, 0.0, 0.0)))
    assert len(node.children) < ops._UNION_CULL_MIN_CHILDREN
    pts = sample(seed=3, n=500, span=((-3.0, 7.0), (-3.0, 3.0), (-3.0, 3.0)))
    expected = np.minimum(evaluate(node.children[0], pts),
                          evaluate(node.children[1], pts))
    assert np.array_equal(evaluate(node, pts), expected)


def test_points_inside_a_child_are_never_culled():
    # A point inside a child reads a negative value that no bound on distance
    # constrains, so overlapping children must always be evaluated.
    node = many_boxes(count=20)
    centres = np.array([[i * 3.0, 0.0, 0.0] for i in range(20)], dtype=float)
    values = evaluate(node, centres)
    assert np.all(values < 0.0)


def test_an_empty_point_batch_is_handled():
    node = many_boxes()
    out = evaluate(node, np.zeros((0, 3)))
    assert out.shape == (0,)
