"""Solid mass-property queries: area, volume centroid and the DSL builtins.

Two families of defect are pinned here (issues #107 and #109).

The DSL ``volume()`` and ``surface_area()`` builtins imported ``geom3d``
names that had never been written, and ``is_empty()`` type-checked but had
no runtime, so all three compiled and then crashed.

``centroid()`` averaged mesh vertices.  That is biased toward densely
tessellated regions, and for a solid with no surfaces -- a deferred field,
or an empty CSG result -- it reported the origin wherever the part was.
It now computes the volume centroid and refuses what it cannot measure;
``vertex_centroid()`` keeps the old estimate.
"""

import math

import pytest

from yapcad import sdf
from yapcad.dsl import compile_and_run
from yapcad.dsl.runtime.builtins import call_builtin
from yapcad.dsl.runtime.values import solid_val
from yapcad.geom import point
from yapcad.geom3d import (
    centroidof,
    reversesurface,
    solid,
    surface_area,
    surface_center,
    surface_normal,
    translatesolid,
    vertex_centroid,
    volumeof,
)
from yapcad.geom3d_util import conic, prism, sphere
from yapcad.sdf.deferred import DeferredMeshError


def xyz(p):
    return [p[0], p[1], p[2]]


def two_boxes():
    """A 2x2x2 box at the origin (V=8) and a 4x2x2 box at x=10 (V=16).

    Both have 24 mesh vertices, so the vertex mean sits halfway between
    them while the volume centroid is pulled toward the bigger one.
    """
    a = prism(2, 2, 2)
    b = prism(4, 2, 2, center=point(10, 0, 0))
    return solid(a[1] + b[1])


def hollow_box():
    """10-unit cube with a 2-unit void centred at x=3 (inward-wound)."""
    outer = prism(10, 10, 10)
    inner = prism(2, 2, 2, center=point(3, 0, 0))
    return solid(outer[1] + [reversesurface(s) for s in inner[1]])


def deferred_block(at=(100.0, 0.0, 0.0)):
    with sdf.deferring():
        part = sdf.to_solid(sdf.box((20.0, 16.0, 10.0)))
    return translatesolid(part, point(*at))


class TestSurfaceArea:
    def test_cube(self):
        assert surface_area(prism(2, 3, 4)) == pytest.approx(52.0)

    def test_accepts_a_single_surface(self):
        assert surface_area(prism(2, 3, 4)[1][0]) > 0.0

    def test_void_shells_count(self):
        assert surface_area(hollow_box()) == pytest.approx(600.0 + 24.0)

    def test_cylinder_converges_on_the_analytic_area(self):
        r, h = 3.0, 5.0
        cyl = conic(r, r, h, max_chord=0.05)
        exact = 2 * math.pi * r * h + 2 * math.pi * r * r
        assert surface_area(cyl) == pytest.approx(exact, rel=2e-3)

    def test_empty_solid_has_zero_area(self):
        assert surface_area(solid([])) == 0.0

    def test_deferred_solid_is_refused(self):
        with pytest.raises(DeferredMeshError):
            surface_area(deferred_block())


class TestSurfaceCenterAndNormal:
    @pytest.mark.parametrize("surf", prism(2, 4, 6, center=point(1, 2, 3))[1])
    def test_each_prism_face(self, surf):
        c = surface_center(surf)
        n = surface_normal(surf)
        # A face centre lies on one bounding plane, offset from the body
        # centre along exactly the outward normal.
        offset = [c[k] - (1, 2, 3)[k] for k in range(3)]
        half = {0: 1.0, 1: 2.0, 2: 3.0}
        axis = max(range(3), key=lambda k: abs(n[k]))
        assert abs(n[axis]) == pytest.approx(1.0)
        assert offset[axis] == pytest.approx(math.copysign(half[axis], n[axis]))
        for k in range(3):
            if k != axis:
                assert offset[k] == pytest.approx(0.0, abs=1e-12)

    def test_closed_shell_has_no_mean_normal(self):
        shell = sphere(4.0)[1][0]
        with pytest.raises(ValueError, match="no mean normal"):
            surface_normal(shell)

    def test_center_is_area_weighted_not_vertex_weighted(self):
        # The sphere's skin is symmetric about its centre however unevenly
        # the subdivision spaces the vertices.
        c = surface_center(sphere(4.0, center=point(1, -2, 3))[1][0])
        assert xyz(c) == pytest.approx([1, -2, 3], abs=1e-9)


class TestCentroid:
    def test_translated_cube(self):
        c = centroidof(translatesolid(prism(2, 2, 2), point(5, 1, -2)))
        assert xyz(c) == pytest.approx([5, 1, -2])

    def test_volume_weighted_across_bodies(self):
        part = two_boxes()
        assert centroidof(part)[0] == pytest.approx((0 * 8 + 10 * 16) / 24)
        # The legacy estimate weights the bodies by vertex count instead.
        assert vertex_centroid(part)[0] == pytest.approx(5.0)

    def test_cone_is_a_quarter_up_the_axis(self):
        h = 8.0
        cone = conic(3.0, 0.0, h, max_chord=0.1)
        c = centroidof(cone)
        assert c[2] == pytest.approx(h / 4, rel=1e-3)
        # The apex is one vertex and the base has many, so the vertex mean
        # sits far below the true centroid.
        assert vertex_centroid(cone)[2] < h / 8

    def test_independent_of_tessellation(self):
        coarse = centroidof(conic(3.0, 1.0, 6.0, max_chord=1.0))
        fine = centroidof(conic(3.0, 1.0, 6.0, max_chord=0.05))
        assert coarse[2] == pytest.approx(fine[2], rel=1e-2)

    def test_void_shell_subtracts(self):
        c = centroidof(hollow_box())
        assert c[0] == pytest.approx(-8 * 3 / (1000 - 8))
        assert volumeof(hollow_box()) == pytest.approx(992.0)

    def test_deferred_solid_is_refused_not_placed_at_the_origin(self):
        part = deferred_block()
        with pytest.raises(DeferredMeshError):
            centroidof(part)
        with pytest.raises(DeferredMeshError):
            vertex_centroid(part)

    def test_empty_solid_has_no_centroid(self):
        with pytest.raises(ValueError, match="empty"):
            centroidof(solid([]))
        with pytest.raises(ValueError, match="no vertices"):
            vertex_centroid(solid([]))

    def test_open_solid_is_refused(self):
        open_box = solid(prism(2, 2, 2)[1][:5])
        with pytest.raises(ValueError, match="closed"):
            centroidof(open_box)


class TestDslQueries:
    def test_volume_and_area(self):
        box = solid_val(prism(2, 3, 4))
        assert call_builtin("volume", [box]).data == pytest.approx(24.0)
        assert call_builtin("surface_area", [box]).data == pytest.approx(52.0)

    def test_centroid_is_the_volume_centroid(self):
        part = solid_val(two_boxes())
        assert call_builtin("centroid", [part]).data[0] == pytest.approx(20 / 3)
        assert call_builtin("vertex_centroid", [part]).data[0] == pytest.approx(5.0)

    def test_centroid_of_a_deferred_part_is_refused(self):
        part = solid_val(deferred_block())
        with pytest.raises(DeferredMeshError):
            call_builtin("centroid", [part])

    def test_centroid_of_a_realized_part_is_where_the_part_is(self):
        part = deferred_block()
        sdf.realize(part)
        c = call_builtin("centroid", [solid_val(part)]).data
        assert xyz(c) == pytest.approx([100, 0, 0], abs=1e-6)

    def test_is_empty(self):
        assert call_builtin("is_empty", [solid_val(solid([]))]).data is True
        assert call_builtin("is_empty", [solid_val(prism(1, 1, 1))]).data is False
        with pytest.raises(DeferredMeshError):
            call_builtin("is_empty", [solid_val(deferred_block())])

    def test_queries_compile_and_run(self):
        src = """module q

command MAKE() -> solid:
    b: solid = translate(box(2.0, 2.0, 2.0), 5.0, 0.0, 0.0)
    v: float = volume(b)
    a: float = surface_area(b)
    c: point3d = centroid(b)
    vc: point3d = vertex_centroid(b)
    e: bool = is_empty(b)
    e2: bool = is_empty(empty_solid())
    emit b
"""
        result = compile_and_run(src, "MAKE", {})
        assert result.success, result.error_message
