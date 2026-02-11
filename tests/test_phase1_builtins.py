"""Tests for Phase 1 DSL builtins: the 9 functions that were in symbols.py but missing from builtins.py."""

import math
import pytest
from yapcad.dsl.runtime.builtins import get_builtin_registry, call_builtin
from yapcad.dsl.runtime.values import (
    Value, int_val, float_val, bool_val, string_val, list_val,
    point_val, vector_val, transform_val, solid_val, region2d_val,
)
from yapcad.dsl.types import (
    FLOAT, INT, BOOL, STRING, POINT, POINT3D, VECTOR3D, VECTOR2D,
    TRANSFORM, SOLID, REGION2D, BEZIER, ListType,
)
from yapcad.geom import point, vect


@pytest.fixture
def registry():
    return get_builtin_registry()


class TestScaleUniform:
    def test_exists(self, registry):
        func = registry.get_function("scale_uniform")
        assert func is not None, "scale_uniform should be registered"

    def test_returns_transform(self):
        result = call_builtin("scale_uniform", [float_val(2.0)])
        assert result.type == TRANSFORM

    def test_uniform_scaling(self):
        """Verify the matrix scales uniformly."""
        result = call_builtin("scale_uniform", [float_val(3.0)])
        mat = result.data
        # Apply to a point
        p = point(1, 1, 1)
        scaled = mat.mul(p)
        assert abs(scaled[0] - 3.0) < 1e-6
        assert abs(scaled[1] - 3.0) < 1e-6
        assert abs(scaled[2] - 3.0) < 1e-6


class TestRotate2d:
    def test_exists(self, registry):
        func = registry.get_function("rotate_2d")
        assert func is not None, "rotate_2d should be registered"

    def test_returns_transform(self):
        result = call_builtin("rotate_2d", [float_val(90.0)])
        assert result.type == TRANSFORM

    def test_90_degree_rotation(self):
        """Rotate (1,0,0) by 90° around Z should give (0,1,0)."""
        result = call_builtin("rotate_2d", [float_val(90.0)])
        mat = result.data
        p = point(1, 0, 0)
        rotated = mat.mul(p)
        assert abs(rotated[0] - 0.0) < 1e-6
        assert abs(rotated[1] - 1.0) < 1e-6


class TestMirrorY:
    def test_exists(self, registry):
        func = registry.get_function("mirror_y")
        assert func is not None, "mirror_y should be registered"

    def test_returns_transform(self):
        result = call_builtin("mirror_y", [])
        assert result.type == TRANSFORM

    def test_negates_x(self):
        """mirror_y should negate X coordinate."""
        result = call_builtin("mirror_y", [])
        mat = result.data
        p = point(5, 3, 1)
        mirrored = mat.mul(p)
        assert abs(mirrored[0] - (-5.0)) < 1e-6
        assert abs(mirrored[1] - 3.0) < 1e-6
        assert abs(mirrored[2] - 1.0) < 1e-6


class TestMirror2d:
    def test_exists(self, registry):
        func = registry.get_function("mirror_2d")
        assert func is not None, "mirror_2d should be registered"

    def test_returns_transform(self):
        axis = vector_val(vect(1, 0), is_2d=True)
        result = call_builtin("mirror_2d", [axis])
        assert result.type == TRANSFORM

    def test_mirror_across_x_axis(self):
        """Mirror across X axis should negate Y."""
        axis = vector_val(vect(1, 0), is_2d=True)
        result = call_builtin("mirror_2d", [axis])
        mat = result.data
        p = point(3, 5, 0)
        mirrored = mat.mul(p)
        assert abs(mirrored[0] - 3.0) < 1e-6
        assert abs(mirrored[1] - (-5.0)) < 1e-6


class TestMirrorSolid:
    def test_exists(self, registry):
        func = registry.get_function("mirror")
        assert func is not None, "mirror should be registered"

    def test_mirror_box(self):
        """Mirror a box across YZ plane (normal = 1,0,0)."""
        from yapcad.geom3d_util import prism
        box = prism(10, 10, 10)
        box_val = solid_val(box)
        normal = vector_val(vect(1, 0, 0), is_2d=False)
        result = call_builtin("mirror", [box_val, normal])
        assert result.type == SOLID
        assert result.data is not None


class TestBezier:
    def test_exists(self, registry):
        func = registry.get_function("bezier")
        assert func is not None, "bezier should be registered"

    def test_creates_curve(self):
        """Create a bezier from 4 control points."""
        pts = [
            point_val(point(0, 0, 0)),
            point_val(point(1, 2, 0)),
            point_val(point(3, 2, 0)),
            point_val(point(4, 0, 0)),
        ]
        pts_list = list_val(pts, POINT)
        result = call_builtin("bezier", [pts_list])
        # Should return some curve type
        assert result.data is not None


class TestCentroid:
    def test_exists(self, registry):
        func = registry.get_function("centroid")
        assert func is not None, "centroid should be registered"

    def test_box_centroid(self):
        """Centroid of a centered box should be near origin."""
        from yapcad.geom3d_util import prism
        # prism creates a box centered at origin (-w/2..w/2, etc.)
        box = prism(10, 10, 10)
        result = call_builtin("centroid", [solid_val(box)])
        c = result.data
        assert abs(c[0]) < 1.0  # Should be near origin
        assert abs(c[1]) < 1.0
        assert abs(c[2]) < 1.0


class TestDistance:
    def test_exists(self, registry):
        func = registry.get_function("distance")
        assert func is not None, "distance should be registered"

    def test_simple_distance(self):
        a = point_val(point(0, 0, 0))
        b = point_val(point(3, 4, 0))
        tol = float_val(0.001)
        result = call_builtin("distance", [a, b, tol])
        assert abs(result.data - 5.0) < 1e-6

    def test_3d_distance(self):
        a = point_val(point(1, 2, 3))
        b = point_val(point(4, 6, 3))
        tol = float_val(0.001)
        result = call_builtin("distance", [a, b, tol])
        expected = math.sqrt(9 + 16)  # 5.0
        assert abs(result.data - expected) < 1e-6


class TestLoft:
    def test_exists(self, registry):
        func = registry.get_function("loft")
        assert func is not None, "loft should be registered"

    def test_loft_two_rectangles(self):
        """Loft between two rectangles at different heights."""
        from yapcad.geom import point, line

        # Bottom rectangle at z=0
        r1 = [
            line(point(-5, -5, 0), point(5, -5, 0)),
            line(point(5, -5, 0), point(5, 5, 0)),
            line(point(5, 5, 0), point(-5, 5, 0)),
            line(point(-5, 5, 0), point(-5, -5, 0)),
        ]

        # Top rectangle (smaller) at z=10
        r2 = [
            line(point(-3, -3, 10), point(3, -3, 10)),
            line(point(3, -3, 10), point(3, 3, 10)),
            line(point(3, 3, 10), point(-3, 3, 10)),
            line(point(-3, 3, 10), point(-3, -3, 10)),
        ]

        profiles = list_val([region2d_val(r1), region2d_val(r2)], REGION2D)
        result = call_builtin("loft", [profiles])
        assert result.type == SOLID
        assert result.data is not None
