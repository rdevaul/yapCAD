"""Degeneracy is a question of shape, asked in one place (#85), and dual
contouring does not leave the slivers and near-coincident vertices that
break float32 consumers (#121)."""

import collections
import io
import struct
import sys
import warnings

import numpy as np
import pytest

from yapcad import sdf
from yapcad.geom import point
from yapcad.geom3d import issolidclosed, solid, surface, translatesolid
from yapcad.geometry_utils import (
    degenerate_faces,
    triangle_is_degenerate,
    triangle_normal,
)
from yapcad.io.stl import write_stl


def _equilateral(scale, offset=(0.0, 0.0, 0.0)):
    ox, oy, oz = offset
    return [(ox, oy, oz), (ox + scale, oy, oz),
            (ox + 0.5 * scale, oy + 0.8660254 * scale, oz)]


@pytest.mark.parametrize("scale", [1e-6, 1e-3, 1.0, 1e3, 1e6])
def test_a_healthy_triangle_has_a_normal_at_any_scale(scale):
    assert triangle_normal(*_equilateral(scale)) is not None
    assert not triangle_is_degenerate(*_equilateral(scale))


@pytest.mark.parametrize("scale", [1e-3, 1.0, 1e3])
def test_a_collinear_triangle_is_degenerate_at_any_scale(scale):
    tri = [(0.0, 0.0, 0.0), (scale, 0.0, 0.0), (0.3 * scale, 1e-12 * scale, 0.0)]
    assert triangle_normal(*tri) is None
    assert triangle_is_degenerate(*tri)


def test_the_floor_follows_the_coordinates_not_the_size():
    # A 0.01 mm triangle that is exactly collinear in exact arithmetic,
    # stored 100 m from the origin: only the representation floor catches
    # the noise its cross product carries.
    rng = np.random.default_rng(3)
    caught = 0
    for _ in range(50):
        d = rng.normal(size=3)
        d /= np.linalg.norm(d)
        base = np.array([1e5, -7e4, 3e4])
        tri = [tuple(base), tuple(base + 0.01 * d), tuple(base + 0.004 * d)]
        caught += triangle_normal(*tri) is None
    assert caught == 50
    # A healthy triangle there still has its normal.
    assert triangle_normal(*_equilateral(0.01, (1e5, -7e4, 3e4))) is not None


def test_vectorised_and_scalar_tests_agree():
    rng = np.random.default_rng(7)
    V = rng.normal(size=(40, 3)) * 10.0
    V[5] = V[4] + 1e-13 * (V[6] - V[4])           # nearly on edge 4-6
    V[9] = 0.5 * (V[8] + V[10])                    # exactly collinear
    faces = [[i, i + 1, i + 2] for i in range(0, 37)]
    mask = degenerate_faces([tuple(v) for v in V], faces)
    expected = [triangle_normal(*(tuple(V[k]) for k in f)) is None for f in faces]
    assert mask.tolist() == expected
    assert sum(expected) >= 2


def test_a_solid_closed_only_by_zero_area_faces_is_not_closed():
    # A tetrahedron with one face split at a point on an edge, giving a
    # zero-area face that keeps the index topology closed.
    pts = [point(0, 0, 0), point(1, 0, 0), point(0, 1, 0), point(0, 0, 1),
           point(0.5, 0, 0)]
    good = [[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]]
    assert issolidclosed(solid([surface(pts[:4], [[0, 0, 1, 0]] * 4, good)], [], []))
    # Replace face [0, 1, 3] by [0, 4, 3] + [4, 1, 3], and [0, 2, 1] by
    # [0, 2, 4] + [4, 2, 1]: still closed, nothing degenerate.
    split = [[0, 2, 4], [4, 2, 1], [0, 4, 3], [4, 1, 3], [1, 2, 3], [0, 3, 2]]
    assert issolidclosed(solid([surface(pts, [[0, 0, 1, 0]] * 5, split)], [], []))
    # Instead close the split edge with a zero-area sliver [0, 4, 1]: the
    # index topology is closed, but no consumer can use the sliver.
    sliver = [[0, 2, 1], [0, 4, 3], [4, 1, 3], [1, 4, 0], [1, 2, 3], [0, 3, 2]]
    faces_closed = collections.Counter(
        tuple(sorted((f[k], f[(k + 1) % 3]))) for f in sliver for k in range(3))
    assert all(c == 2 for c in faces_closed.values())
    assert not issolidclosed(solid([surface(pts, [[0, 0, 1, 0]] * 5, sliver)], [], []))


# ---------------------------------------------------------------------------
# dual contouring
# ---------------------------------------------------------------------------

BLOCK = sdf.subtract(sdf.box([60, 40, 20]),
                     sdf.translate(sdf.cylinder(6, 30), [10, 5, 0]),
                     sdf.translate(sdf.cylinder(4, 30), [-17, -9, 0]))


def _stl_defects(obj):
    buf = io.BytesIO()
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        write_stl(obj, buf)
    data = buf.getvalue()
    n = struct.unpack("<I", data[80:84])[0]
    index, faces = {}, []
    for i in range(n):
        v = struct.unpack("<12f", data[84 + 50 * i: 84 + 50 * i + 48])[3:]
        faces.append([index.setdefault(v[k:k + 3], len(index)) for k in (0, 3, 6)])
    edges = collections.Counter(
        tuple(sorted((f[k], f[(k + 1) % 3]))) for f in faces for k in range(3))
    return (sum(len(set(f)) < 3 for f in faces),
            sum(c != 2 for c in edges.values()))


def test_no_sub_micron_edges_in_a_lattice_clamped_mesh():
    mesh = sdf.dual_contour(BLOCK, cell_mm=0.5)
    V, F = np.asarray(mesh.vertices), np.asarray(mesh.triangles)
    T = V[F]
    shortest = min(np.linalg.norm(T[:, i] - T[:, j], axis=1).min()
                   for i, j in ((0, 1), (1, 2), (2, 0)))
    assert shortest > 1e-3
    assert not any(sdf.manifold_defects(mesh).values())


@pytest.mark.parametrize("offset", [0.0, 600.0, 5000.0])
def test_a_field_mesh_stays_closed_in_float32_stl_away_from_the_origin(offset):
    solid_ = sdf.to_solid(BLOCK, cell_mm=0.5, refine=False)
    assert _stl_defects(translatesolid(solid_, [offset, 0.0, 0.0])) == (0, 0)


def test_write_stl_warns_when_float32_collapses_an_edge():
    far = 1000.0
    pts = [point(far, 0, 0), point(far + 1e-5, 0, 0), point(far, 1, 0)]
    tri = solid([surface(pts, [[0, 0, 1, 0]] * 3, [[0, 1, 2]])], [], [])
    with pytest.warns(RuntimeWarning, match="float32"):
        write_stl(tri, io.BytesIO())


def test_the_field_miter_gear_has_no_degenerate_faces():
    from yapcad.gears.bevel import (StraightBevelGearSpec,
                                    make_straight_bevel_gear_sdf)
    spec = StraightBevelGearSpec(teeth=24, mate_teeth=24, outer_module_mm=1.5,
                                 face_width_mm=8.0, backlash_mm=0.25,
                                 bore_diameter_mm=8.4)
    for simplify in (False, True):
        gear = make_straight_bevel_gear_sdf(spec, simplify=simplify)
        surf = gear[1][0]
        assert not degenerate_faces(surf[1], surf[3]).any()
        assert issolidclosed(gear)


def test_strokes_thinner_than_a_few_cells_mesh_without_coincident_vertices():
    # Engraved block text at the resolutions #94 reported refused.
    from yapcad.text3d import _text_to_polygons_block
    cutters = []
    for row, text in enumerate(["YAPCAD FIELD TEST", "ENGRAVED PLATE 42"]):
        polys = _text_to_polygons_block(text, 3.0, 1.0)
        width = max(p[1][0] for p in polys)
        for p in polys:
            x0, y0, x1, y1 = p[0][0], p[0][1], p[2][0], p[2][1]
            cutters.append(sdf.translate(
                sdf.box([x1 - x0, y1 - y0, 1.0]),
                [(x0 + x1) / 2 - width / 2, (y0 + y1) / 2 + 3 - row * 6, 2.0]))
    plate = sdf.subtract(sdf.box([70, 20, 4]), sdf.union(*cutters))
    for resolution in (350, 490):
        mesh = sdf.dual_contour(plate, resolution=resolution)
        assert sdf.manifold_defects(mesh) == {
            "boundary": 0, "nonmanifold": 0, "coincident": 0, "degenerate": 0,
        }, resolution
