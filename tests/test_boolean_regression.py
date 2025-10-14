#!/usr/bin/env python3
"""
Test for Boolean operation regression found in example12.py

This test documents a regression where Boolean union operations were crashing
due to issurface() and issolid() not properly checking list lengths before
accessing elements.

The crash has been fixed, but there remains a logic bug in combineglist where
union operations don't correctly include geometry from both operands.
"""

import pytest
from yapcad.geom import *
from yapcad.poly import *
from yapcad.combine import *


class TestBooleanRegression:
    """Tests for Boolean operation bugs found in example12.py"""

    def test_boolean_no_crash_on_circles(self):
        """Boolean operations on circles should not crash"""
        circle1 = Circle(point(0, 0), 10)
        circle2 = Circle(point(15, 0), 10)

        # This used to crash with IndexError in issolid()
        nb = Boolean('union', [circle1, circle2])
        result = nb.geom

        # Should return something (even if logic is wrong)
        assert isinstance(result, list)

    def test_boolean_no_crash_empty_result(self):
        """Boolean operations should handle empty results without crashing"""
        # Create two non-overlapping, non-contained shapes
        circle1 = Circle(point(0, 0), 5)
        circle2 = Circle(point(100, 0), 5)

        # Should not crash even if result handling is imperfect
        nb = Boolean('union', [circle1, circle2])
        result = nb.geom

        assert isinstance(result, list)

    def test_boolean_union_includes_both_circles(self):
        """Boolean union should include geometry from both operands."""
        circle1 = Circle(point(0, 0), 10)
        circle2 = Circle(point(15, 0), 10)

        nb = Boolean('union', [circle1, circle2])
        result = nb.geom

        # Boolean operations convert to line segments, not arcs
        # Check that result has geometry spanning both circles
        assert len(result) > 0, "Result should not be empty"

        # Check that some line segments are near circle1's extent
        has_circle1_geom = any(
            isline(elem) and (elem[0][0] < 0 or elem[1][0] < 0)  # Left of origin
            for elem in result
        )

        # Check that some line segments are near circle2's extent
        has_circle2_geom = any(
            isline(elem) and (elem[0][0] > 20 or elem[1][0] > 20)  # Right of circle2 center
            for elem in result
        )

        assert has_circle1_geom, "Missing geometry from first circle"
        assert has_circle2_geom, "Missing geometry from second circle"

    def test_boolean_union_bbox_spans_both(self):
        """Boolean union bbox should span both operands."""
        circle1 = Circle(point(0, 0), 10)
        circle2 = Circle(point(15, 0), 10)

        nb = Boolean('union', [circle1, circle2])
        result = nb.geom

        result_bbox = bbox(result)

        # Union should span from x=-10 (left of circle1) to x=25 (right of circle2)
        assert result_bbox[0][0] <= -9.9, f"Left bound {result_bbox[0][0]} should be near -10"
        assert result_bbox[1][0] >= 24.9, f"Right bound {result_bbox[1][0]} should be near 25"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
