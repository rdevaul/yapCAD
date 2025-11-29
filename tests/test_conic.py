"""Tests for the parabola and hyperbola (conic) primitives."""

import pytest
from math import cosh, sinh

from yapcad.geom import (
    parabola, isparabola, parabola_sample, parabola_tangent,
    parabola_length, parabola_bbox,
    hyperbola, ishyperbola, hyperbola_sample, hyperbola_tangent,
    hyperbola_length, hyperbola_bbox,
    sample, length, bbox, center, unsample, point
)


# ============================================================================
# Parabola Tests
# ============================================================================


class TestParabolaConstruction:
    """Test parabola construction and validation."""

    def test_basic_parabola(self):
        """Test basic parabola creation."""
        p = parabola([0, 0, 0], 1.0)
        assert isparabola(p)
        assert p[0] == 'parabola'
        assert p[2]['focal_length'] == 1.0

    def test_parabola_with_vertex(self):
        """Test parabola with non-origin vertex."""
        p = parabola([1, 2, 3], 2.0)
        assert isparabola(p)
        assert p[1][0] == 1
        assert p[1][1] == 2
        assert p[1][2] == 3

    def test_parabola_with_rotation(self):
        """Test rotated parabola."""
        p = parabola([0, 0, 0], 1.0, rotation=45)
        assert isparabola(p)
        assert p[2]['rotation'] == 45.0

    def test_parabola_with_range(self):
        """Test parabola with custom parameter range."""
        p = parabola([0, 0, 0], 1.0, start=-5, end=5)
        assert isparabola(p)
        assert p[2]['start'] == -5.0
        assert p[2]['end'] == 5.0

    def test_parabola_with_normal(self):
        """Test parabola with custom normal."""
        p = parabola([0, 0, 0], 1.0, normal=[1, 0, 0])
        assert isparabola(p)
        norm = p[2]['normal']
        assert abs(norm[0] - 1.0) < 1e-10
        assert abs(norm[1]) < 1e-10
        assert abs(norm[2]) < 1e-10

    def test_parabola_invalid_focal_length(self):
        """Test that invalid focal_length raises error."""
        with pytest.raises(ValueError):
            parabola([0, 0, 0], -1.0)
        with pytest.raises(ValueError):
            parabola([0, 0, 0], 0)


class TestParabolaSampling:
    """Test parabola sampling and evaluation."""

    def test_sample_at_vertex(self):
        """Test that u=0.5 is at vertex for symmetric range."""
        p = parabola([0, 0, 0], 1.0, start=-10, end=10)
        # At t=0 (which is u=0.5 for this range), we should be at the vertex
        pt = parabola_sample(p, 0.5)
        assert abs(pt[0]) < 1e-10
        assert abs(pt[1]) < 1e-10

    def test_sample_at_start(self):
        """Test sampling at start of range."""
        p = parabola([0, 0, 0], 1.0, start=-2, end=2)
        pt = parabola_sample(p, 0)
        # At t=-2: x = 4/(4*1) = 1, y = -2
        assert abs(pt[0] - 1.0) < 1e-10
        assert abs(pt[1] - (-2.0)) < 1e-10

    def test_sample_at_end(self):
        """Test sampling at end of range."""
        p = parabola([0, 0, 0], 1.0, start=-2, end=2)
        pt = parabola_sample(p, 1)
        # At t=2: x = 4/(4*1) = 1, y = 2
        assert abs(pt[0] - 1.0) < 1e-10
        assert abs(pt[1] - 2.0) < 1e-10

    def test_generic_sample_function(self):
        """Test that generic sample() works with parabola."""
        p = parabola([0, 0, 0], 1.0, start=-2, end=2)
        pt = sample(p, 0.5)
        assert abs(pt[0]) < 1e-10  # At vertex

    def test_parabola_with_offset_vertex(self):
        """Test parabola with offset vertex."""
        p = parabola([5, 5, 0], 1.0, start=-2, end=2)
        pt = parabola_sample(p, 0.5)
        assert abs(pt[0] - 5.0) < 1e-10
        assert abs(pt[1] - 5.0) < 1e-10


class TestParabolaTangent:
    """Test parabola tangent computation."""

    def test_tangent_at_vertex(self):
        """Test tangent at vertex (t=0)."""
        p = parabola([0, 0, 0], 1.0, start=-10, end=10)
        # At t=0, dx/dt = 0, dy/dt = 1, so tangent is (0, 1)
        t = parabola_tangent(p, 0.5)
        assert abs(t[0]) < 1e-10
        assert abs(t[1] - 1.0) < 1e-10

    def test_tangent_is_unit_vector(self):
        """Test that tangent is normalized."""
        p = parabola([0, 0, 0], 1.0, start=-5, end=5)
        for u in [0, 0.25, 0.5, 0.75, 1.0]:
            t = parabola_tangent(p, u)
            mag = (t[0]**2 + t[1]**2 + t[2]**2)**0.5
            assert abs(mag - 1.0) < 1e-10


class TestParabolaLength:
    """Test parabola length computation."""

    def test_length_positive(self):
        """Test that length is positive."""
        p = parabola([0, 0, 0], 1.0, start=-5, end=5)
        l = parabola_length(p)
        assert l > 0

    def test_generic_length_function(self):
        """Test that generic length() works with parabola."""
        p = parabola([0, 0, 0], 1.0, start=-2, end=2)
        l = length(p)
        assert l > 0


class TestParabolaBbox:
    """Test parabola bounding box computation."""

    def test_symmetric_parabola_bbox(self):
        """Test bbox of symmetric parabola."""
        p = parabola([0, 0, 0], 1.0, start=-2, end=2)
        bb = parabola_bbox(p)
        # x ranges from 0 (vertex) to 1 (at t=±2)
        # y ranges from -2 to 2
        assert bb[0][0] >= -0.1  # min x near 0
        assert bb[1][0] <= 1.1   # max x near 1
        assert bb[0][1] >= -2.1  # min y near -2
        assert bb[1][1] <= 2.1   # max y near 2

    def test_generic_bbox_function(self):
        """Test that generic bbox() works with parabola."""
        p = parabola([0, 0, 0], 1.0, start=-2, end=2)
        bb = bbox(p)
        assert bb[0][0] < bb[1][0]  # min x < max x


class TestParabolaCenter:
    """Test parabola center (vertex) computation."""

    def test_center_function(self):
        """Test that generic center() returns parabola vertex."""
        p = parabola([1, 2, 3], 1.0)
        c = center(p)
        assert abs(c[0] - 1) < 1e-10
        assert abs(c[1] - 2) < 1e-10
        assert abs(c[2] - 3) < 1e-10


class TestParabolaUnsample:
    """Test parabola unsample (inverse sampling)."""

    def test_unsample_roundtrip(self):
        """Test that unsample inverts sample."""
        p = parabola([0, 0, 0], 1.0, start=-5, end=5)
        for u in [0.0, 0.25, 0.5, 0.75, 1.0]:
            pt = sample(p, u)
            u_back = unsample(p, pt)
            assert u_back is not False
            assert abs(u_back - u) < 0.05


# ============================================================================
# Hyperbola Tests
# ============================================================================


class TestHyperbolaConstruction:
    """Test hyperbola construction and validation."""

    def test_basic_hyperbola(self):
        """Test basic hyperbola creation."""
        h = hyperbola([0, 0, 0], 2.0, 1.0)
        assert ishyperbola(h)
        assert h[0] == 'hyperbola'
        assert h[2]['semi_major'] == 2.0
        assert h[2]['semi_minor'] == 1.0

    def test_hyperbola_with_center(self):
        """Test hyperbola with non-origin center."""
        h = hyperbola([1, 2, 3], 2.0, 1.0)
        assert ishyperbola(h)
        assert h[1][0] == 1
        assert h[1][1] == 2
        assert h[1][2] == 3

    def test_hyperbola_with_rotation(self):
        """Test rotated hyperbola."""
        h = hyperbola([0, 0, 0], 2.0, 1.0, rotation=45)
        assert ishyperbola(h)
        assert h[2]['rotation'] == 45.0

    def test_hyperbola_with_range(self):
        """Test hyperbola with custom parameter range."""
        h = hyperbola([0, 0, 0], 2.0, 1.0, start=-1, end=1)
        assert ishyperbola(h)
        assert h[2]['start'] == -1.0
        assert h[2]['end'] == 1.0

    def test_hyperbola_left_branch(self):
        """Test hyperbola on left branch."""
        h = hyperbola([0, 0, 0], 2.0, 1.0, branch=-1)
        assert ishyperbola(h)
        assert h[2]['branch'] == -1

    def test_hyperbola_with_normal(self):
        """Test hyperbola with custom normal."""
        h = hyperbola([0, 0, 0], 2.0, 1.0, normal=[1, 0, 0])
        assert ishyperbola(h)
        norm = h[2]['normal']
        assert abs(norm[0] - 1.0) < 1e-10
        assert abs(norm[1]) < 1e-10
        assert abs(norm[2]) < 1e-10

    def test_hyperbola_invalid_axes(self):
        """Test that invalid axes raise errors."""
        with pytest.raises(ValueError):
            hyperbola([0, 0, 0], -1.0, 1.0)
        with pytest.raises(ValueError):
            hyperbola([0, 0, 0], 1.0, 0)

    def test_hyperbola_invalid_branch(self):
        """Test that invalid branch raises error."""
        with pytest.raises(ValueError):
            hyperbola([0, 0, 0], 2.0, 1.0, branch=0)


class TestHyperbolaSampling:
    """Test hyperbola sampling and evaluation."""

    def test_sample_at_vertex(self):
        """Test that u=0.5 is at vertex (t=0) for symmetric range."""
        h = hyperbola([0, 0, 0], 2.0, 1.0, start=-2, end=2)
        # At t=0: x = a*cosh(0) = a = 2, y = b*sinh(0) = 0
        pt = hyperbola_sample(h, 0.5)
        assert abs(pt[0] - 2.0) < 1e-10
        assert abs(pt[1]) < 1e-10

    def test_sample_right_branch(self):
        """Test sampling on right branch (positive x)."""
        h = hyperbola([0, 0, 0], 2.0, 1.0, start=-1, end=1, branch=1)
        # All points on right branch should have x > 0
        for u in [0, 0.25, 0.5, 0.75, 1.0]:
            pt = hyperbola_sample(h, u)
            assert pt[0] > 0

    def test_sample_left_branch(self):
        """Test sampling on left branch (negative x)."""
        h = hyperbola([0, 0, 0], 2.0, 1.0, start=-1, end=1, branch=-1)
        # All points on left branch should have x < 0
        for u in [0, 0.25, 0.5, 0.75, 1.0]:
            pt = hyperbola_sample(h, u)
            assert pt[0] < 0

    def test_generic_sample_function(self):
        """Test that generic sample() works with hyperbola."""
        h = hyperbola([0, 0, 0], 2.0, 1.0, start=-1, end=1)
        pt = sample(h, 0.5)
        assert abs(pt[0] - 2.0) < 1e-10  # At vertex

    def test_hyperbola_equation_satisfied(self):
        """Test that sampled points satisfy hyperbola equation."""
        a, b = 3.0, 2.0
        h = hyperbola([0, 0, 0], a, b, start=-1.5, end=1.5)
        for u in [0.1, 0.3, 0.5, 0.7, 0.9]:
            pt = hyperbola_sample(h, u)
            x, y = pt[0], pt[1]
            # x²/a² - y²/b² should equal 1
            result = (x**2 / a**2) - (y**2 / b**2)
            assert abs(result - 1.0) < 1e-9


class TestHyperbolaTangent:
    """Test hyperbola tangent computation."""

    def test_tangent_at_vertex(self):
        """Test tangent at vertex (t=0)."""
        h = hyperbola([0, 0, 0], 2.0, 1.0, start=-2, end=2)
        # At t=0: dx/dt = a*sinh(0) = 0, dy/dt = b*cosh(0) = b
        # So tangent should be (0, 1) after normalization
        t = hyperbola_tangent(h, 0.5)
        assert abs(t[0]) < 1e-10
        assert abs(t[1] - 1.0) < 1e-10

    def test_tangent_is_unit_vector(self):
        """Test that tangent is normalized."""
        h = hyperbola([0, 0, 0], 2.0, 1.0, start=-1, end=1)
        for u in [0, 0.25, 0.5, 0.75, 1.0]:
            t = hyperbola_tangent(h, u)
            mag = (t[0]**2 + t[1]**2 + t[2]**2)**0.5
            assert abs(mag - 1.0) < 1e-10


class TestHyperbolaLength:
    """Test hyperbola length computation."""

    def test_length_positive(self):
        """Test that length is positive."""
        h = hyperbola([0, 0, 0], 2.0, 1.0, start=-1, end=1)
        l = hyperbola_length(h)
        assert l > 0

    def test_generic_length_function(self):
        """Test that generic length() works with hyperbola."""
        h = hyperbola([0, 0, 0], 2.0, 1.0, start=-1, end=1)
        l = length(h)
        assert l > 0


class TestHyperbolaBbox:
    """Test hyperbola bounding box computation."""

    def test_right_branch_bbox(self):
        """Test bbox of right branch hyperbola."""
        h = hyperbola([0, 0, 0], 2.0, 1.0, start=-1, end=1, branch=1)
        bb = hyperbola_bbox(h)
        # x should be positive (right branch)
        assert bb[0][0] > 0  # min x positive

    def test_generic_bbox_function(self):
        """Test that generic bbox() works with hyperbola."""
        h = hyperbola([0, 0, 0], 2.0, 1.0, start=-1, end=1)
        bb = bbox(h)
        assert bb[0][0] < bb[1][0]  # min x < max x


class TestHyperbolaCenter:
    """Test hyperbola center computation."""

    def test_center_function(self):
        """Test that generic center() returns hyperbola center."""
        h = hyperbola([1, 2, 3], 2.0, 1.0)
        c = center(h)
        assert abs(c[0] - 1) < 1e-10
        assert abs(c[1] - 2) < 1e-10
        assert abs(c[2] - 3) < 1e-10


class TestHyperbolaUnsample:
    """Test hyperbola unsample (inverse sampling)."""

    def test_unsample_roundtrip(self):
        """Test that unsample inverts sample."""
        h = hyperbola([0, 0, 0], 2.0, 1.0, start=-1, end=1)
        for u in [0.0, 0.25, 0.5, 0.75, 1.0]:
            pt = sample(h, u)
            u_back = unsample(h, pt)
            assert u_back is not False
            assert abs(u_back - u) < 0.05
