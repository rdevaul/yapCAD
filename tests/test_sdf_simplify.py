"""Field-checked simplification of meshed fields.

Dual contouring meshes a part at the cell size its finest feature needs,
everywhere.  Simplifying against the field recovers most of that, with the
field itself proving the result has not strayed from the part.
"""

import numpy as np
import pytest

from yapcad import sdf
from yapcad.construction import get_construction
from yapcad.geom3d import issolidclosed, volumeof
from yapcad.geom3d_util import prism
from yapcad.sdf.simplify import simplify_mesh, simplify_solid

PLATE = sdf.subtract(sdf.rounded_box((40, 30, 6), 1.0),
                     sdf.translate(sdf.cylinder(4, 10), (8, 0, 0)),
                     sdf.translate(sdf.cylinder(4, 10), (-8, 0, 0)))


def arrays(solid):
    surf = solid[1][0]
    return (np.array([p[:3] for p in surf[1]], dtype=float),
            np.array(surf[3], dtype=np.int64))


def area_error(node, V, F, n=200_000, seed=0):
    """|f| at points spread uniformly by area -- not per triangle, which
    would over-weight the many small triangles left at sharp edges."""
    rng = np.random.default_rng(seed)
    A, B, C = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    area = 0.5 * np.linalg.norm(np.cross(B - A, C - A), axis=1)
    pick = rng.choice(len(F), size=n, p=area / area.sum())
    u, v = rng.random((n, 1)), rng.random((n, 1))
    flip = (u + v) > 1
    u, v = np.where(flip, 1 - u, u), np.where(flip, 1 - v, v)
    P = A[pick] + u * (B - A)[pick] + v * (C - A)[pick]
    return np.abs(sdf.evaluate(node, P))


@pytest.fixture(scope="module")
def plate():
    return sdf.to_solid(PLATE, resolution=64)


@pytest.fixture(scope="module")
def simplified(plate):
    return simplify_solid(plate, 0.01)


def test_simplification_keeps_a_closed_solid(simplified):
    assert issolidclosed(simplified)


def test_it_removes_most_of_the_triangles(plate, simplified):
    assert len(simplified[1][0][3]) < 0.3 * len(plate[1][0][3])


def test_it_is_no_less_accurate_than_the_uniform_mesh(plate, simplified):
    before = area_error(PLATE, *arrays(plate))
    after = area_error(PLATE, *arrays(simplified))
    tol = 0.01
    assert np.percentile(after, 99) <= max(tol, np.percentile(before, 99))
    assert (after > tol).mean() <= (before > tol).mean() + 0.01
    assert after.max() <= before.max() * 1.06


def test_volume_is_preserved(plate, simplified):
    assert volumeof(simplified) == pytest.approx(volumeof(plate), rel=2e-3)


def test_the_record_keeps_the_field_and_the_tolerance(simplified):
    record = get_construction(simplified)
    assert sdf.from_construction(record) == PLATE
    meshing = sdf.meshing_from_construction(record)
    assert meshing["simplify"] == {"method": "field-qem", "tolerance": 0.01}
    assert meshing["resolution"] == 64


def test_simplification_is_deterministic(plate):
    V, F = arrays(plate)
    first = simplify_mesh(PLATE, V, F, 0.01)
    second = simplify_mesh(PLATE, V, F, 0.01)
    assert np.array_equal(first[0], second[0])
    assert np.array_equal(first[1], second[1])


def test_a_sphere_simplifies_within_tolerance():
    node = sdf.sphere(10.0)
    V, F = arrays(sdf.to_solid(node, resolution=48))
    V2, F2 = simplify_mesh(node, V, F, 0.02)
    assert len(F2) < 0.5 * len(F)
    # An exact field, so |f| is the distance.  Merged vertices are moved onto
    # the sphere; the rest keep dual contouring's placement, so no vertex
    # ends up further off than one already was.
    before = np.abs(np.linalg.norm(V, axis=1) - 10.0).max()
    assert np.abs(np.linalg.norm(V2, axis=1) - 10.0).max() <= before
    assert np.percentile(area_error(node, V2, F2), 99) <= 0.02


def test_bad_arguments_are_refused(plate):
    V, F = arrays(plate)
    with pytest.raises(sdf.SdfError):
        simplify_mesh(PLATE, V, F, 0.0)
    with pytest.raises(sdf.SdfError):
        simplify_mesh(PLATE, V, F, "fine")
    with pytest.raises(sdf.SdfError):
        simplify_solid(prism(1, 1, 1), 0.01)
