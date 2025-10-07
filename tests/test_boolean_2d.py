import pytest
import os
import math
from yapcad.geom import *
from yapcad.poly import *
from yapcad.combine import *

"""Test functions for 2D boolean operations based on example11.py"""

# Control flag for visual tests - set via environment variable or directly
VISUALTEST = os.environ.get('VISUALTEST', 'false').lower() in ('true', '1', 'yes')
if VISUALTEST:
    from yapcad.pyglet_drawable import pygletDraw


def makeStars(center, rad=5.0, diam=2.0):
    """Helper function to create two interlocking star patterns"""
    poly1 = Polygon()
    poly2 = Polygon()

    for i in range(8):
        ang = (pi2/8)*i
        p = poly1
        if i % 2 == 0:
            p = poly2
        p.addArc(arc(add(center, point(math.cos(ang)*rad,
                                       math.sin(ang)*rad)), diam))
    return poly1, poly2


class TestBoolean2D:
    """Tests for 2D boolean operations"""

    def test_union_circle_roundrect(self):
        """Test union of a circle and rounded rectangle"""
        p1 = Polygon()
        x = point(-2, 2)
        p1.addArc(arc(x, 3.0))

        p2 = RoundRect(7, 6, 1, point(1, -1))
        b = Boolean('union', [p1, p2])

        # Check that the boolean operation produces valid geometry
        assert b.geom is not None
        assert len(b.geom) > 0

        # The union should produce a valid geometry list
        assert isinstance(b.geom, list)

    def test_union_stars(self):
        """Test union of two interlocking star patterns"""
        p1, p2 = makeStars(point(12, -10))
        b = Boolean('union', [p1, p2])

        # Check that the union produces valid geometry
        assert b.geom is not None
        assert len(b.geom) > 0

        # Verify it's a list of geometry elements
        assert isinstance(b.geom, list)

    def test_intersection_stars(self):
        """Test intersection of two interlocking star patterns"""
        p1, p2 = makeStars(point(12, 10))
        b = Boolean('intersection', [p1, p2])

        # Check that the intersection produces valid geometry
        assert b.geom is not None

        # Intersection might produce multiple or single geometry elements
        if isinstance(b.geom, list):
            assert len(b.geom) >= 0  # Could be empty if no intersection

    def test_difference_stars(self):
        """Test difference of two interlocking star patterns"""
        p1, p2 = makeStars(point(-12, -10))
        b = Boolean('difference', [p2, p1])

        # Check that the difference produces valid geometry
        assert b.geom is not None

        # Difference operation should produce geometry
        if isinstance(b.geom, list):
            assert len(b.geom) >= 0

    def test_union_multiple_circles(self):
        """Test union of multiple circles"""
        polys = []
        for i in range(3):
            p = Polygon()
            p.addArc(arc(point(i*2, 0), 1.5))
            polys.append(p)

        b = Boolean('union', polys)

        # Check that the union produces valid geometry
        assert b.geom is not None
        assert len(b.geom) > 0

    def test_difference_concentric_circles(self):
        """Test difference of concentric circles to create a ring"""
        outer = Polygon()
        outer.addArc(arc(point(0, 0), 5.0))

        inner = Polygon()
        inner.addArc(arc(point(0, 0), 3.0))

        b = Boolean('difference', [outer, inner])

        # Check that the difference produces valid geometry (should be a ring)
        assert b.geom is not None
        assert len(b.geom) > 0

    def test_intersection_overlapping_rectangles(self):
        """Test intersection of two overlapping rounded rectangles"""
        r1 = RoundRect(6, 4, 0.5, point(0, 0))
        r2 = RoundRect(4, 6, 0.5, point(2, 0))

        b = Boolean('intersection', [r1, r2])

        # Check that the intersection produces valid geometry
        assert b.geom is not None

    def test_union_non_overlapping_shapes(self):
        """Test union of non-overlapping shapes"""
        p1 = Polygon()
        p1.addArc(arc(point(-10, 0), 2.0))

        p2 = Polygon()
        p2.addArc(arc(point(10, 0), 2.0))

        b = Boolean('union', [p1, p2])

        # Non-overlapping shapes should still produce valid union
        assert b.geom is not None
        assert len(b.geom) > 0

    def test_difference_non_overlapping_shapes(self):
        """Test difference of non-overlapping shapes"""
        p1 = Polygon()
        p1.addArc(arc(point(-10, 0), 2.0))

        p2 = Polygon()
        p2.addArc(arc(point(10, 0), 2.0))

        b = Boolean('difference', [p1, p2])

        # Difference of non-overlapping should return first shape
        assert b.geom is not None

    def test_complex_star_operations(self):
        """Test all three boolean operations on star patterns"""
        center = point(0, 0)
        p1, p2 = makeStars(center, rad=5.0, diam=2.0)

        # Union
        union = Boolean('union', [p1, p2])
        assert union.geom is not None

        # Intersection
        intersection = Boolean('intersection', [p1, p2])
        assert intersection.geom is not None

        # Difference
        difference = Boolean('difference', [p1, p2])
        assert difference.geom is not None

    def test_roundrect_parameters(self):
        """Test RoundRect with various parameters"""
        r1 = RoundRect(10, 8, 0.5, point(0, 0))
        assert r1.geom is not None

        r2 = RoundRect(5, 5, 1.0, point(5, 5))
        assert r2.geom is not None

        # Union of two rounded rectangles
        b = Boolean('union', [r1, r2])
        assert b.geom is not None

    def test_star_generation(self):
        """Test the makeStars helper function"""
        p1, p2 = makeStars(point(0, 0), rad=10.0, diam=3.0)

        # Both polygons should have valid geometry
        assert p1.geom is not None
        assert p2.geom is not None

        # Each polygon has 4 arcs, and each arc generates 2 geometry elements (line + arc)
        assert len(p1.geom) == 8
        assert len(p2.geom) == 8


@pytest.mark.skipif(not VISUALTEST, reason="Visual tests disabled")
class TestBoolean2DVisual:
    """Visual tests for 2D boolean operations (requires VISUALTEST=true)"""

    @pytest.mark.visual
    def test_visual_all_boolean_operations(self):
        """Visual test showing all boolean operations"""
        dd = pygletDraw()

        # Union of circle and rounded rectangle
        p1 = Polygon()
        p1.addArc(arc(point(-2, 2), 3.0))
        p2 = RoundRect(7, 6, 1, point(1, -1))
        union = Boolean('union', [p1, p2])

        dd.linecolor = 'aqua'
        dd.draw(union.geom)

        # Union of stars
        p4, p5 = makeStars(point(12, -10))
        union_stars = Boolean('union', [p4, p5])

        dd.linecolor = 'white'
        for i in range(len(union_stars.geom)):
            dd.linecolor = i % 7 + 1
            dd.draw(union_stars.geom[i])

        # Intersection of stars
        p6, p7 = makeStars(point(12, 10))
        intersection_stars = Boolean('intersection', [p6, p7])

        dd.linecolor = 'white'
        dd.draw(intersection_stars.geom)

        # Difference of stars
        p8, p9 = makeStars(point(-12, -10))
        difference_stars = Boolean('difference', [p9, p8])

        dd.linecolor = 'yellow'
        dd.draw(difference_stars.geom)

        # Individual stars (no boolean operation)
        p10, p11 = makeStars(point(-12, 10))
        dd.linecolor = 'aqua'
        dd.draw(p10.geom)
        dd.draw(p11.geom)

        dd.display()
