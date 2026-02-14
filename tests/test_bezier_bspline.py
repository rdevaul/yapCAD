"""Test suite for Bezier and B-spline curve support in yapCAD."""

import math
import pytest

from yapcad.geom import (
    bezier,
    bspline,
    isbezier,
    isbspline,
    dist,
    length,
    center,
    bbox,
    point,
    sample,
)
from yapcad.spline import (
    is_bezier,
    is_bspline,
    bezier_point,
    bezier_curve,
    bezier_tangent,
    sample_bezier,
    evaluate_bezier,
    bspline_point,
    bspline_curve,
    bspline_tangent,
    sample_bspline,
    evaluate_bspline,
)


def _close(a, b, tol=1e-6):
    """Assert two points are within tolerance."""
    assert dist(point(a), point(b)) <= tol


class TestBezierConstruction:
    """Tests for Bezier curve construction."""

    def test_bezier_cubic(self):
        """Test creating a cubic Bezier curve."""
        cp = [point(0, 0), point(10, 20), point(30, 20), point(40, 0)]
        curve = bezier(cp)
        assert isbezier(curve)
        assert is_bezier(curve)
        assert curve[0] == 'bezier'
        assert len(curve[1]) == 4
        assert curve[2]['degree'] == 3

    def test_bezier_quadratic(self):
        """Test creating a quadratic Bezier curve."""
        cp = [point(0, 0), point(15, 30), point(30, 0)]
        curve = bezier(cp)
        assert isbezier(curve)
        assert curve[2]['degree'] == 2

    def test_bezier_linear(self):
        """Test creating a linear Bezier curve (straight line)."""
        cp = [point(0, 0), point(10, 10)]
        curve = bezier(cp)
        assert isbezier(curve)
        assert curve[2]['degree'] == 1

    def test_bezier_requires_two_points(self):
        """Test that Bezier requires at least 2 control points."""
        with pytest.raises(ValueError):
            bezier([point(0, 0)])


class TestBezierEvaluation:
    """Tests for Bezier curve evaluation."""

    def test_bezier_endpoints(self):
        """Test that Bezier passes through first and last control points."""
        cp = [point(0, 0), point(10, 20), point(30, 20), point(40, 0)]
        curve = bezier(cp)

        # Start point
        p0 = sample(curve, 0.0)
        _close(p0, cp[0])

        # End point
        p1 = sample(curve, 1.0)
        _close(p1, cp[-1])

    def test_bezier_midpoint(self):
        """Test Bezier curve midpoint for known case."""
        # For a symmetric cubic Bezier with control points forming a parabola,
        # the midpoint should be at the peak
        cp = [point(0, 0), point(10, 20), point(30, 20), point(40, 0)]
        curve = bezier(cp)
        mid = sample(curve, 0.5)
        # Expected midpoint: ((0+3*10+3*30+40)/8, (0+3*20+3*20+0)/8) = (20, 15)
        _close(mid, point(20, 15), tol=1e-4)

    def test_bezier_linear_is_line(self):
        """Test that linear Bezier is equivalent to a line."""
        cp = [point(0, 0), point(10, 10)]
        curve = bezier(cp)

        for t in [0.0, 0.25, 0.5, 0.75, 1.0]:
            p = sample(curve, t)
            expected = point(t * 10, t * 10)
            _close(p, expected)


class TestBezierCurveSampling:
    """Tests for sampling Bezier curves as polylines."""

    def test_bezier_curve_segments(self):
        """Test that bezier_curve returns correct number of points."""
        cp = [point(0, 0), point(10, 20), point(30, 20), point(40, 0)]
        polyline = bezier_curve(cp, segments=32)
        assert len(polyline) == 33  # segments + 1 points

    def test_sample_bezier(self):
        """Test sample_bezier on curve definition."""
        cp = [point(0, 0), point(10, 20), point(30, 20), point(40, 0)]
        curve = bezier(cp)
        polyline = sample_bezier(curve, segments=16)
        assert len(polyline) == 17

        # First and last points match endpoints
        _close(polyline[0], cp[0])
        _close(polyline[-1], cp[-1])


class TestBezierTangent:
    """Tests for Bezier curve tangent calculation."""

    def test_bezier_tangent_direction(self):
        """Test that tangent at t=0 points toward second control point."""
        cp = [point(0, 0), point(10, 20), point(30, 20), point(40, 0)]
        tang = bezier_tangent(cp, 0.0)

        # Tangent at t=0 should be proportional to P1 - P0
        # Direction: (10, 20, 0)
        mag_t = (tang[0]**2 + tang[1]**2 + tang[2]**2)**0.5
        if mag_t > 0:
            normalized = [tang[0]/mag_t, tang[1]/mag_t, tang[2]/mag_t]
            expected_dir = [10/math.sqrt(500), 20/math.sqrt(500), 0]
            assert abs(normalized[0] - expected_dir[0]) < 1e-6
            assert abs(normalized[1] - expected_dir[1]) < 1e-6


class TestBSplineConstruction:
    """Tests for B-spline curve construction."""

    def test_bspline_cubic(self):
        """Test creating a cubic B-spline."""
        cp = [point(0, 0), point(10, 20), point(20, -10), point(30, 15), point(40, 0)]
        curve = bspline(cp, degree=3)
        assert isbspline(curve)
        assert is_bspline(curve)
        assert curve[0] == 'bspline'
        assert curve[2]['degree'] == 3
        assert curve[2]['closed'] is False

    def test_bspline_closed(self):
        """Test creating a closed B-spline."""
        cp = [point(0, 0), point(10, 20), point(20, 0), point(10, -20)]
        curve = bspline(cp, degree=3, closed=True)
        assert isbspline(curve)
        assert curve[2]['closed'] is True

    def test_bspline_requires_minimum_points(self):
        """Test that B-spline requires degree+1 control points."""
        with pytest.raises(ValueError):
            bspline([point(0, 0), point(10, 10)], degree=3)


class TestBSplineEvaluation:
    """Tests for B-spline curve evaluation."""

    def test_bspline_clamped_endpoints(self):
        """Test that clamped B-spline passes through endpoints."""
        cp = [point(0, 0), point(10, 20), point(20, -10), point(30, 15), point(40, 0)]
        curve = bspline(cp, degree=3, closed=False)

        # Start point
        p0 = sample(curve, 0.0)
        _close(p0, cp[0], tol=1e-4)

        # End point
        p1 = sample(curve, 1.0)
        _close(p1, cp[-1], tol=1e-4)

    def test_bspline_curve_sampling(self):
        """Test sampling B-spline as polyline."""
        cp = [point(0, 0), point(10, 20), point(20, -10), point(30, 15), point(40, 0)]
        polyline = bspline_curve(cp, degree=3, segments=32)
        assert len(polyline) == 33


class TestGeometryFunctions:
    """Tests for yapCAD geometry functions with Bezier/B-spline."""

    def test_bezier_length(self):
        """Test length calculation for Bezier curve."""
        cp = [point(0, 0), point(10, 0)]  # Straight line
        curve = bezier(cp)
        l = length(curve)
        assert abs(l - 10.0) < 0.1  # Should be close to 10

    def test_bezier_center(self):
        """Test center calculation for Bezier curve."""
        cp = [point(0, 0), point(10, 20), point(30, 20), point(40, 0)]
        curve = bezier(cp)
        c = center(curve)
        # Center should be roughly in the middle of the curve's extent
        assert c[0] >= 0 and c[0] <= 40
        assert c[1] >= 0

    def test_bezier_bbox(self):
        """Test bounding box calculation for Bezier curve."""
        cp = [point(0, 0), point(10, 20), point(30, 20), point(40, 0)]
        curve = bezier(cp)
        box = bbox(curve)
        # Bounding box should contain all control points
        assert box[0][0] <= 0 and box[1][0] >= 40
        # Y extent: curve doesn't exceed control polygon
        assert box[0][1] >= -1 and box[1][1] <= 21

    def test_bspline_length(self):
        """Test length calculation for B-spline curve."""
        cp = [point(0, 0), point(0, 10), point(10, 10), point(10, 0)]
        curve = bspline(cp, degree=2)
        l = length(curve)
        assert l > 0  # Should have positive length


class TestBezierAndBSplineIntegration:
    """Integration tests combining Bezier and B-spline with other yapCAD features."""

    def test_bezier_sample_evaluate_consistency(self):
        """Test that sample() and evaluate_bezier() give same results."""
        cp = [point(0, 0), point(10, 20), point(30, 20), point(40, 0)]
        curve = bezier(cp)

        for t in [0.0, 0.25, 0.5, 0.75, 1.0]:
            p1 = sample(curve, t)
            p2 = evaluate_bezier(curve, t)
            _close(p1, p2)

    def test_bspline_sample_evaluate_consistency(self):
        """Test that sample() and evaluate_bspline() give same results."""
        cp = [point(0, 0), point(10, 20), point(20, -10), point(30, 15), point(40, 0)]
        curve = bspline(cp, degree=3)

        for t in [0.0, 0.25, 0.5, 0.75, 1.0]:
            p1 = sample(curve, t)
            p2 = evaluate_bspline(curve, t)
            _close(p1, p2)


class TestBezier3D:
    """Tests for 3D Bezier curves."""

    def test_bezier_3d(self):
        """Test Bezier curve in 3D space."""
        cp = [
            point(0, 0, 0),
            point(10, 20, 5),
            point(30, 20, 10),
            point(40, 0, 15)
        ]
        curve = bezier(cp)

        # Endpoints
        p0 = sample(curve, 0.0)
        _close(p0, cp[0])

        p1 = sample(curve, 1.0)
        _close(p1, cp[-1])

        # Midpoint should have Z between 0 and 15
        mid = sample(curve, 0.5)
        assert mid[2] >= 0 and mid[2] <= 15


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
