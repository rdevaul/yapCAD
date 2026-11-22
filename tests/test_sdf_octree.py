"""Octree sampling of the corner lattice (stage A of
docs/sdf-adaptive-contouring.md): fewer field evaluations, and the same mesh
byte for byte."""

import numpy as np
import pytest

from yapcad import sdf
from yapcad.sdf import contour, octree, tape
from yapcad.sdf.node import SdfError


def _box_off_lattice():
    return sdf.translate(sdf.box([10.3, 7.1, 4.9]), [0.37, -0.21, 0.13])


def _plate():
    cutters = [sdf.translate(sdf.box([0.4, 3.0, 1.0]), [x, 0.0, 1.6])
               for x in np.arange(-8.0, 8.5, 1.1)]
    return sdf.subtract(sdf.box([20.0, 6.0, 4.0]), sdf.union(*cutters))


FIXTURES = {
    "box on lattice": (sdf.box([8.0, 8.0, 8.0]), dict(cell_mm=0.5)),
    "box off lattice": (_box_off_lattice(), dict(resolution=40)),
    "sphere": (sdf.sphere(6.0), dict(resolution=48)),
    "torus": (sdf.torus(8.0, 2.0), dict(resolution=56)),
    "gyroid slab": (sdf.intersect(sdf.gyroid(6.0, 0.8),
                                  sdf.box([20.0, 20.0, 4.0])),
                    dict(resolution=60)),
    "wall 1.5 cells": (sdf.box([10.0, 10.0, 0.375]), dict(cell_mm=0.25)),
    "pocket": (sdf.subtract(sdf.box([20.0, 20.0, 6.0]),
                            sdf.translate(sdf.box([0.5, 8.0, 3.0]),
                                          [0.0, 0.0, 1.6])),
               dict(cell_mm=0.2)),
    "shared face": (sdf.union(sdf.box([5.0, 5.0, 5.0]),
                              sdf.translate(sdf.box([5.0, 5.0, 5.0]),
                                            [5.0, 0.0, 0.0])),
                    dict(resolution=20)),
    "nested shells": (sdf.union(sdf.shell(sdf.sphere(8.0), 0.6),
                                sdf.sphere(4.0)),
                      dict(resolution=64)),
    "engraved plate": (_plate(), dict(cell_mm=0.1)),
}


@pytest.mark.parametrize("name", sorted(FIXTURES))
def test_octree_sampling_gives_the_same_mesh_byte_for_byte(name):
    node, kw = FIXTURES[name]
    chunks = sdf.dual_contour(node, **kw)
    tree = sdf.dual_contour(node, sampling="octree", **kw)
    assert len(chunks.triangles) > 0
    for a, b in zip(chunks, tree):
        assert a.dtype == b.dtype and a.shape == b.shape
        assert a.tobytes() == b.tobytes()


@pytest.mark.skipif(not tape.numba_available(), reason="numba not installed")
@pytest.mark.parametrize("name", ["engraved plate", "gyroid slab"])
def test_the_compiled_backend_gives_the_same_mesh_too(name):
    node, kw = FIXTURES[name]
    chunks = sdf.dual_contour(node, backend="numba", **kw)
    tree = sdf.dual_contour(node, backend="numba", sampling="octree", **kw)
    for a, b in zip(chunks, tree):
        assert a.tobytes() == b.tobytes()


@pytest.mark.parametrize("leaf", [1, 2, 4, 8])
def test_every_leaf_size_gives_true_signs_and_exact_crossings(leaf):
    node, kw = FIXTURES["torus"]
    origin, spacing, counts = contour._resolve_region(
        node, None, kw["resolution"], None)
    dense = contour.sample_corners(node, origin, spacing, counts,
                                   prune=False)
    sample = octree.sample_corners_octree(node, origin, spacing, counts,
                                          leaf=leaf)
    assert np.array_equal(sample.values < 0.0, dense < 0.0)
    for axis in range(3):
        crossing = contour._edge_crossings(dense, axis)
        lower = [slice(None)] * 3
        upper = [slice(None)] * 3
        lower[axis] = slice(0, -1)
        upper[axis] = slice(1, None)
        for part in (lower, upper):
            assert np.array_equal(sample.values[tuple(part)][crossing],
                                  dense[tuple(part)][crossing])


def test_it_samples_fewer_points_than_chunk_pruning(monkeypatch):
    node = sdf.sphere(20.0)
    origin, spacing, counts = contour._resolve_region(node, None, 256, None)
    seen = []
    real = contour.evaluate_batches

    def counting(node_, batches, backend=None):
        seen.extend(len(b) for b in batches)
        return real(node_, batches, backend=backend)

    monkeypatch.setattr(contour, "evaluate_batches", counting)
    contour.sample_corners(node, origin, spacing, counts)
    monkeypatch.undo()
    sample = octree.sample_corners_octree(node, origin, spacing, counts)
    total = sample.centre_evaluations + sample.corner_evaluations
    assert total < 0.3 * sum(seen)


def test_a_lattice_that_is_not_cubic_is_covered():
    node = sdf.box([30.0, 3.0, 1.0])
    origin, spacing, counts = contour._resolve_region(node, None, 97, None)
    assert len(set(counts)) == 3
    sample = octree.sample_corners_octree(node, origin, spacing, counts)
    assert sample.values.shape == tuple(c + 1 for c in counts)
    assert np.isfinite(sample.values).all()
    assert octree.root_size(counts) >= max(counts)


def test_to_solid_passes_sampling_through():
    node, kw = FIXTURES["sphere"]
    a = sdf.to_solid(node, sampling="chunks", **kw)
    b = sdf.to_solid(node, sampling="octree", **kw)
    assert a[1][0][1] == b[1][0][1]
    assert a[1][0][3] == b[1][0][3]


def test_it_is_deterministic():
    node, kw = FIXTURES["gyroid slab"]
    first = sdf.dual_contour(node, sampling="octree", **kw)
    second = sdf.dual_contour(node, sampling="octree", **kw)
    for a, b in zip(first, second):
        assert a.tobytes() == b.tobytes()


def test_bad_arguments_are_refused():
    with pytest.raises(SdfError, match="sampling"):
        sdf.dual_contour(sdf.sphere(1.0), sampling="adaptive")
    origin, spacing, counts = contour._resolve_region(
        sdf.sphere(1.0), None, 8, None)
    with pytest.raises(SdfError, match="power of two"):
        octree.sample_corners_octree(sdf.sphere(1.0), origin, spacing,
                                     counts, leaf=3)
