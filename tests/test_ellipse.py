"""Tests for the ellipse primitive."""

import pytest
from math import pi

from yapcad.geom import (
    ellipse, isellipse, isfullellipse,
    ellipse_sample, ellipse_tangent, ellipse_length, ellipse_bbox,
    sample, length, bbox, center, unsample, point
)


class TestEllipseConstruction:
    """Test ellipse construction and validation."""

    def test_basic_ellipse(self):
        """Test basic ellipse creation."""
        e = ellipse([0, 0, 0], 2.0, 1.0)
        assert isellipse(e)
        assert isfullellipse(e)
        assert e[0] == 'ellipse'
        assert e[2]['semi_major'] == 2.0
        assert e[2]['semi_minor'] == 1.0

    def test_ellipse_with_center(self):
        """Test ellipse with non-origin center."""
        e = ellipse([1, 2, 3], 3.0, 2.0)
        assert isellipse(e)
        assert e[1][0] == 1
        assert e[1][1] == 2
        assert e[1][2] == 3

    def test_ellipse_with_rotation(self):
        """Test rotated ellipse."""
        e = ellipse([0, 0, 0], 2.0, 1.0, rotation=45)
        assert isellipse(e)
        assert e[2]['rotation'] == 45.0

    def test_ellipse_arc(self):
        """Test elliptical arc (not full ellipse)."""
        e = ellipse([0, 0, 0], 2.0, 1.0, start=0, end=90)
        assert isellipse(e)
        assert not isfullellipse(e)
        assert e[2]['start'] == 0
        assert e[2]['end'] == 90.0

    def test_ellipse_swaps_axes(self):
        """Test that semi_minor > semi_major swaps axes."""
        e = ellipse([0, 0, 0], 1.0, 2.0)  # minor > major
        # Should swap so semi_major is always larger
        assert e[2]['semi_major'] == 2.0
        assert e[2]['semi_minor'] == 1.0

    def test_ellipse_with_normal(self):
        """Test ellipse with custom normal."""
        e = ellipse([0, 0, 0], 2.0, 1.0, normal=[1, 0, 0])
        assert isellipse(e)
        norm = e[2]['normal']
        assert abs(norm[0] - 1.0) < 1e-10
        assert abs(norm[1]) < 1e-10
        assert abs(norm[2]) < 1e-10


class TestEllipseSampling:
    """Test ellipse sampling and evaluation."""

    def test_sample_at_cardinal_points(self):
        """Test sampling at 0, 0.25, 0.5, 0.75."""
        e = ellipse([0, 0, 0], 2.0, 1.0)

        # u=0 should be at (semi_major, 0)
        p0 = ellipse_sample(e, 0)
        assert abs(p0[0] - 2.0) < 1e-10
        assert abs(p0[1]) < 1e-10

        # u=0.25 should be at (0, semi_minor)
        p25 = ellipse_sample(e, 0.25)
        assert abs(p25[0]) < 1e-10
        assert abs(p25[1] - 1.0) < 1e-10

    def test_generic_sample_function(self):
        """Test that generic sample() works with ellipse."""
        e = ellipse([0, 0, 0], 2.0, 1.0)
        pt = sample(e, 0)
        assert abs(pt[0] - 2.0) < 1e-10

    def test_ellipse_arc_sampling(self):
        """Test sampling on ellipse arc."""
        e = ellipse([0, 0, 0], 2.0, 1.0, start=0, end=90)

        # u=0 should be at start of arc
        p0 = ellipse_sample(e, 0)
        assert abs(p0[0] - 2.0) < 1e-10

        # u=1 should be at end of arc (90 degrees)
        p1 = ellipse_sample(e, 1)
        assert abs(p1[0]) < 1e-10
        assert abs(p1[1] - 1.0) < 1e-10


class TestEllipseTangent:
    """Test ellipse tangent computation."""

    def test_tangent_at_cardinal_points(self):
        """Test tangent vectors at cardinal points."""
        e = ellipse([0, 0, 0], 2.0, 1.0)

        # At u=0 (right side), tangent should point up (0, 1)
        t0 = ellipse_tangent(e, 0)
        assert abs(t0[0]) < 1e-10
        assert abs(t0[1] - 1.0) < 1e-10

    def test_tangent_is_unit_vector(self):
        """Test that tangent is normalized."""
        e = ellipse([0, 0, 0], 2.0, 1.0)
        for u in [0, 0.1, 0.25, 0.5, 0.75]:
            t = ellipse_tangent(e, u)
            mag = (t[0]**2 + t[1]**2 + t[2]**2)**0.5
            assert abs(mag - 1.0) < 1e-10


class TestEllipseLength:
    """Test ellipse length computation."""

    def test_circle_length(self):
        """Test that a circle (a=b) has correct perimeter."""
        e = ellipse([0, 0, 0], 1.0, 1.0)  # Actually a circle
        l = ellipse_length(e)
        expected = 2 * pi * 1.0  # 2*pi*r
        assert abs(l - expected) < 0.01  # Numerical approximation

    def test_generic_length_function(self):
        """Test that generic length() works with ellipse."""
        e = ellipse([0, 0, 0], 2.0, 1.0)
        l = length(e)
        assert l > 0


class TestEllipseBbox:
    """Test ellipse bounding box computation."""

    def test_unrotated_ellipse_bbox(self):
        """Test bbox of axis-aligned ellipse."""
        e = ellipse([0, 0, 0], 2.0, 1.0)
        bb = ellipse_bbox(e)
        assert abs(bb[0][0] - (-2.0)) < 0.1
        assert abs(bb[1][0] - 2.0) < 0.1
        assert abs(bb[0][1] - (-1.0)) < 0.1
        assert abs(bb[1][1] - 1.0) < 0.1

    def test_generic_bbox_function(self):
        """Test that generic bbox() works with ellipse."""
        e = ellipse([0, 0, 0], 2.0, 1.0)
        bb = bbox(e)
        assert bb[0][0] < 0  # min x negative
        assert bb[1][0] > 0  # max x positive


class TestEllipseCenter:
    """Test ellipse center computation."""

    def test_center_function(self):
        """Test that generic center() returns ellipse center."""
        e = ellipse([1, 2, 3], 2.0, 1.0)
        c = center(e)
        assert abs(c[0] - 1) < 1e-10
        assert abs(c[1] - 2) < 1e-10
        assert abs(c[2] - 3) < 1e-10


class TestEllipseUnsample:
    """Test ellipse unsample (inverse sampling)."""

    def test_unsample_roundtrip(self):
        """Test that unsample inverts sample."""
        e = ellipse([0, 0, 0], 2.0, 1.0)
        for u in [0.0, 0.25, 0.5, 0.75]:
            pt = sample(e, u)
            u_back = unsample(e, pt)
            assert u_back is not False
            assert abs(u_back - u) < 0.01
