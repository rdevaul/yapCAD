"""
Tests for radial_pattern() and linear_pattern() functions.
"""

import pytest
import math
from yapcad.geom import point, arc, line, dist, close


class TestRadialPattern:
    """Tests for 2D radial pattern function."""

    def test_radial_pattern_6_copies(self):
        """Test creating 6 holes in a radial pattern."""
        from yapcad.geom_util import radial_pattern

        # Create a circle at radius 10
        hole = arc(point(10, 0), 2.5)
        holes = radial_pattern(hole, count=6)

        assert len(holes) == 6

        # Check that holes are evenly distributed (60 degree intervals)
        expected_angles = [0, 60, 120, 180, 240, 300]
        for i, h in enumerate(holes):
            center = h[0]
            angle = math.degrees(math.atan2(center[1], center[0]))
            if angle < 0:
                angle += 360
            # Normalize expected angle
            expected = expected_angles[i]
            # Allow some floating point tolerance - check within 1 degree
            angle_diff = abs(angle - expected)
            if angle_diff > 180:
                angle_diff = 360 - angle_diff
            assert angle_diff < 1.0, f"Hole {i}: expected {expected}, got {angle}"
            # Check radius is maintained
            r = math.sqrt(center[0]**2 + center[1]**2)
            assert close(r, 10.0)

    def test_radial_pattern_count_1(self):
        """Test that count=1 returns just the original."""
        from yapcad.geom_util import radial_pattern

        hole = arc(point(10, 0), 2.5)
        holes = radial_pattern(hole, count=1)

        assert len(holes) == 1
        assert holes[0] is hole  # Should be the same object

    def test_radial_pattern_custom_center(self):
        """Test radial pattern with custom center point."""
        from yapcad.geom_util import radial_pattern

        hole = arc(point(15, 5), 2.5)  # hole at (15,5)
        center = point(5, 5)  # rotate around (5,5)
        holes = radial_pattern(hole, count=4, center=center)

        assert len(holes) == 4

        # All holes should be 10mm from center (5,5)
        for h in holes:
            hcenter = h[0]
            r = dist(hcenter, center)
            assert close(r, 10.0)

    def test_radial_pattern_partial_angle(self):
        """Test radial pattern with less than 360 degrees."""
        from yapcad.geom_util import radial_pattern

        hole = arc(point(10, 0), 2.5)
        holes = radial_pattern(hole, count=3, angle=90)  # 3 holes over 90 degrees

        assert len(holes) == 3

        # First should be at 0 degrees, last at 90 degrees
        first_center = holes[0][0]
        last_center = holes[2][0]

        assert close(first_center[0], 10.0)
        assert close(first_center[1], 0.0)
        assert close(last_center[0], 0.0)
        assert close(last_center[1], 10.0)

    def test_radial_pattern_invalid_count(self):
        """Test that count < 1 raises error."""
        from yapcad.geom_util import radial_pattern

        hole = arc(point(10, 0), 2.5)

        with pytest.raises(ValueError):
            radial_pattern(hole, count=0)


class TestLinearPattern:
    """Tests for 2D linear pattern function."""

    def test_linear_pattern_4_copies(self):
        """Test creating 4 holes in a linear pattern."""
        from yapcad.geom_util import linear_pattern

        hole = arc(point(0, 0), 3)
        holes = linear_pattern(hole, count=4, spacing=[20, 0, 0])

        assert len(holes) == 4

        # Check positions
        for i, h in enumerate(holes):
            center = h[0]
            assert close(center[0], i * 20.0)
            assert close(center[1], 0.0)

    def test_linear_pattern_count_1(self):
        """Test that count=1 returns just the original."""
        from yapcad.geom_util import linear_pattern

        hole = arc(point(0, 0), 3)
        holes = linear_pattern(hole, count=1, spacing=[20, 0, 0])

        assert len(holes) == 1
        assert holes[0] is hole

    def test_linear_pattern_diagonal(self):
        """Test linear pattern with diagonal spacing."""
        from yapcad.geom_util import linear_pattern

        hole = arc(point(0, 0), 3)
        holes = linear_pattern(hole, count=3, spacing=[10, 10, 0])

        assert len(holes) == 3

        for i, h in enumerate(holes):
            center = h[0]
            assert close(center[0], i * 10.0)
            assert close(center[1], i * 10.0)

    def test_linear_pattern_with_line(self):
        """Test linear pattern on a line segment."""
        from yapcad.geom_util import linear_pattern

        l = line(point(0, 0), point(5, 0))
        lines = linear_pattern(l, count=3, spacing=[10, 0, 0])

        assert len(lines) == 3

        # Check start points of each line
        assert close(lines[0][0][0], 0.0)
        assert close(lines[1][0][0], 10.0)
        assert close(lines[2][0][0], 20.0)

    def test_linear_pattern_invalid_count(self):
        """Test that count < 1 raises error."""
        from yapcad.geom_util import linear_pattern

        hole = arc(point(0, 0), 3)

        with pytest.raises(ValueError):
            linear_pattern(hole, count=0, spacing=[20, 0, 0])


class TestSolidPatterns:
    """Tests for 3D solid pattern functions."""

    def test_radial_pattern_solid(self):
        """Test radial pattern with 3D solid."""
        from yapcad.geom3d_util import conic, radial_pattern_solid
        from yapcad.geom3d import translatesolid, issolid

        cylinder = conic(2.5, 2.5, 10)
        cylinder = translatesolid(cylinder, point(10, 0, 0))

        copies = radial_pattern_solid(cylinder, count=4)

        assert len(copies) == 4
        for s in copies:
            assert issolid(s)

    def test_linear_pattern_solid(self):
        """Test linear pattern with 3D solid."""
        from yapcad.geom3d_util import conic, linear_pattern_solid
        from yapcad.geom3d import issolid

        cylinder = conic(3, 3, 10)
        copies = linear_pattern_solid(cylinder, count=3, spacing=[20, 0, 0])

        assert len(copies) == 3
        for s in copies:
            assert issolid(s)

    def test_pattern_solid_invalid(self):
        """Test that non-solid raises error."""
        from yapcad.geom3d_util import radial_pattern_solid

        not_solid = [1, 2, 3]

        with pytest.raises(ValueError):
            radial_pattern_solid(not_solid, count=4)
