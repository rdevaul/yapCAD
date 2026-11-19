"""Threaded field evaluation and the simplifier's early rejections.

Both are pure speed: neither may change a single value.  yapCAD signs
packages, so "the same mesh, sooner" has to mean byte for byte.
"""

import multiprocessing
import sys

import numpy as np
import pytest

from yapcad import sdf
from yapcad.sdf.contour import dual_contour
from yapcad.sdf.evaluate import (
    THREAD_CHUNK,
    evaluate,
    evaluate_batches,
    gradient,
    thread_count,
)

simplify = sys.modules["yapcad.sdf.simplify"]


def bracket():
    """A part with every kind of node the evaluator treats specially: a
    many-child union (culled and split by batch), blends and transforms."""
    pegs = sdf.union(*[
        sdf.translate(sdf.cylinder(1.0, 6.0), (3.0 * i - 9.0, 0.0, 3.0))
        for i in range(7)])
    body = sdf.smooth_union(sdf.box((24.0, 8.0, 4.0)), pegs, 0.6)
    return sdf.subtract(body, sdf.rotate(sdf.capsule(
        (-8.0, 0.0, -4.0), (8.0, 0.0, 4.0), 0.9), (0.0, 0.0, 1.0), 20.0))


def cloud(count):
    rng = np.random.default_rng(7)
    return rng.uniform((-14.0, -6.0, -4.0), (14.0, 6.0, 8.0), (count, 3))


@pytest.fixture
def threads(monkeypatch):
    def use(count):
        monkeypatch.setenv("YAPCAD_SDF_THREADS", str(count))
    return use


def test_a_split_batch_gives_the_values_of_the_whole_one(threads):
    node, points = bracket(), cloud(4 * THREAD_CHUNK + 123)
    threads(1)
    values, slopes = evaluate(node, points), gradient(node, points)
    threads(4)
    assert np.array_equal(evaluate(node, points), values)
    assert np.array_equal(gradient(node, points), slopes)


def test_batches_come_back_in_order(threads):
    node = bracket()
    batches = [cloud(50 + 10 * n) + n for n in range(9)]
    threads(1)
    expected = evaluate_batches(node, batches)
    threads(4)
    got = evaluate_batches(node, batches)
    assert len(got) == len(batches)
    for a, b in zip(got, expected):
        assert np.array_equal(a, b)


def test_the_mesh_is_the_same_on_any_number_of_threads(threads):
    node = bracket()
    threads(1)
    one = dual_contour(node, resolution=72, chunk=16)
    threads(4)
    four = dual_contour(node, resolution=72, chunk=16)
    for a, b in zip(one, four):
        assert np.array_equal(a, b)


def _in_child(queue):
    queue.put(float(evaluate(bracket(), cloud(3 * THREAD_CHUNK)).sum()))


@pytest.mark.skipif("fork" not in multiprocessing.get_all_start_methods(),
                    reason="needs fork")
def test_the_pool_is_not_carried_across_a_fork(threads):
    """Parts are meshed in worker processes.  A forked child inherits the
    parent's pool object without its threads, and would wait on it for
    ever."""
    threads(4)
    expected = float(evaluate(bracket(), cloud(3 * THREAD_CHUNK)).sum())
    context = multiprocessing.get_context("fork")
    queue = context.Queue()
    child = context.Process(target=_in_child, args=(queue,))
    child.start()
    try:
        assert queue.get(timeout=120) == expected
    finally:
        child.join(timeout=10)
        if child.is_alive():
            child.kill()


def test_thread_count_setting(monkeypatch):
    monkeypatch.setenv("YAPCAD_SDF_THREADS", "3")
    assert thread_count() == 3
    monkeypatch.setenv("YAPCAD_SDF_THREADS", "0")
    assert thread_count() == 1
    monkeypatch.setenv("YAPCAD_SDF_THREADS", "many")
    with pytest.raises(sdf.SdfError):
        thread_count()
    monkeypatch.delenv("YAPCAD_SDF_THREADS")
    assert 1 <= thread_count() <= 8


# -- the simplifier's quick look ---------------------------------------------

def _simplified(node):
    mesh = dual_contour(node, resolution=56)
    spacing = 28.0 / 56
    return simplify.simplify_mesh(node, mesh.vertices, mesh.triangles,
                                  0.05 * spacing)


@pytest.mark.parametrize("subset", [
    lambda table: table,            # the quick look is the full look
    lambda table: table[:1],        # one point: it rejects almost nothing
    lambda table: table[::3],
])
def test_the_quick_look_changes_no_collapse(monkeypatch, subset):
    """Rejecting a collapse early on a few sample points is sound because
    those points are among the full set.  So the result cannot depend on
    which few they are -- including all of them, which is the simplifier
    with no quick look at all."""
    node = bracket()
    vertices, triangles = _simplified(node)
    monkeypatch.setattr(simplify, "_COARSE", {
        n: subset(table) for n, table in simplify._LATTICES.items()})
    again_v, again_t = _simplified(node)
    assert np.array_equal(again_v, vertices)
    assert np.array_equal(again_t, triangles)
    assert len(triangles) < 0.5 * len(
        dual_contour(node, resolution=56).triangles)


def test_coarse_points_are_among_the_full_ones():
    for n, table in simplify._LATTICES.items():
        full = {tuple(row) for row in table.tolist()}
        assert {tuple(row) for row in simplify._COARSE[n].tolist()} <= full
