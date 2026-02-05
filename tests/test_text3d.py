"""Tests for the text3d module."""

import pytest
from yapcad.text3d import (
    text_to_polygons, text_solid, text_width,
    engrave_text, get_supported_characters, BLOCK_FONT
)
from yapcad.geom import point, vect, dist, epsilon
from yapcad.geom3d import issolid, solidbbox
from yapcad.geom3d_util import prism


class TestBlockFont:
    """Test the BLOCK_FONT data structure."""

    def test_font_has_letters(self):
        """Font should contain A-Z."""
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            assert c in BLOCK_FONT, f"Missing letter {c}"

    def test_font_has_digits(self):
        """Font should contain 0-9."""
        for c in '0123456789':
            assert c in BLOCK_FONT, f"Missing digit {c}"

    def test_font_has_punctuation(self):
        """Font should contain basic punctuation."""
        for c in '-_.:/()':
            assert c in BLOCK_FONT, f"Missing punctuation {c}"

    def test_font_has_space(self):
        """Font should contain space (empty list)."""
        assert ' ' in BLOCK_FONT
        assert BLOCK_FONT[' '] == []

    def test_rectangles_are_valid(self):
        """Each rectangle should be (x, y, w, h) with non-negative values."""
        for char, rects in BLOCK_FONT.items():
            for rect in rects:
                assert len(rect) == 4, f"Bad rectangle in {char}: {rect}"
                x, y, w, h = rect
                assert x >= 0, f"Negative x in {char}"
                assert y >= 0, f"Negative y in {char}"
                assert w > 0, f"Non-positive width in {char}"
                assert h > 0, f"Non-positive height in {char}"


class TestTextToPolygons:
    """Test the text_to_polygons function."""

    def test_single_character(self):
        """Single character should produce polygons."""
        polys = text_to_polygons('A', height=5.0)
        assert len(polys) > 0, "Should produce at least one polygon"

    def test_empty_string(self):
        """Empty string should produce no polygons."""
        polys = text_to_polygons('', height=5.0)
        assert polys == []

    def test_space_only(self):
        """Space-only string should produce no polygons."""
        polys = text_to_polygons(' ', height=5.0)
        assert polys == []

    def test_polygon_closed(self):
        """Each polygon should be closed (first point == last point)."""
        polys = text_to_polygons('X', height=5.0)
        for poly in polys:
            assert dist(poly[0], poly[-1]) < epsilon, "Polygon not closed"

    def test_height_scaling(self):
        """Polygon coordinates should scale with height parameter."""
        polys_5 = text_to_polygons('A', height=5.0)
        polys_10 = text_to_polygons('A', height=10.0)

        # Max Y should be roughly proportional to height
        max_y_5 = max(p[1] for poly in polys_5 for p in poly)
        max_y_10 = max(p[1] for poly in polys_10 for p in poly)
        ratio = max_y_10 / max_y_5
        assert 1.9 < ratio < 2.1, f"Height scaling incorrect: {ratio}"


class TestTextWidth:
    """Test the text_width function."""

    def test_empty_string(self):
        """Empty string should have zero width."""
        assert text_width('') == 0.0

    def test_single_character(self):
        """Single character width should be positive."""
        w = text_width('A', height=5.0)
        assert w > 0

    def test_longer_text(self):
        """Longer text should have greater width."""
        w1 = text_width('A', height=5.0)
        w3 = text_width('ABC', height=5.0)
        assert w3 > w1

    def test_height_scaling(self):
        """Width should scale with height."""
        w5 = text_width('HELLO', height=5.0)
        w10 = text_width('HELLO', height=10.0)
        ratio = w10 / w5
        assert 1.9 < ratio < 2.1, f"Width scaling incorrect: {ratio}"


class TestTextSolid:
    """Test the text_solid function."""

    def test_produces_valid_solid(self):
        """text_solid should produce a valid yapCAD solid."""
        solid = text_solid('A', height=5.0, depth=1.0)
        assert issolid(solid), "Result should be a valid solid"

    def test_has_surfaces(self):
        """Result solid should have surfaces."""
        solid = text_solid('HI', height=5.0, depth=1.0)
        assert len(solid[1]) > 0, "Solid should have surfaces"

    def test_empty_string(self):
        """Empty string should produce empty solid."""
        solid = text_solid('', height=5.0, depth=1.0)
        assert issolid(solid)
        assert len(solid[1]) == 0

    def test_depth_parameter(self):
        """Depth should affect Z extent of bounding box."""
        solid_thin = text_solid('X', height=5.0, depth=0.5)
        solid_thick = text_solid('X', height=5.0, depth=2.0)

        bbox_thin = solidbbox(solid_thin)
        bbox_thick = solidbbox(solid_thick)

        z_thin = bbox_thin[1][2] - bbox_thin[0][2]
        z_thick = bbox_thick[1][2] - bbox_thick[0][2]

        assert z_thick > z_thin, "Thicker depth should produce larger Z extent"


class TestEngraveText:
    """Test the engrave_text function."""

    def test_produces_valid_solid(self):
        """engrave_text should produce a valid solid."""
        plate = prism(30, 20, 3)
        result = engrave_text(
            plate, 'X',
            position=point(0, 0, 1.5),
            normal=vect(0, 0, 1, 0),
            height=4.0,
            depth=0.5
        )
        assert issolid(result)

    def test_invalid_target_raises(self):
        """Should raise if target is not a solid."""
        with pytest.raises(ValueError):
            engrave_text(
                [1, 2, 3],  # Not a solid
                'X',
                position=point(0, 0, 0),
                normal=vect(0, 0, 1, 0)
            )

    def test_empty_text_returns_copy(self):
        """Empty text should return copy of original."""
        plate = prism(30, 20, 3)
        result = engrave_text(
            plate, '',
            position=point(0, 0, 1.5),
            normal=vect(0, 0, 1, 0)
        )
        # Should have same number of surfaces as original
        assert len(result[1]) == len(plate[1])


class TestGetSupportedCharacters:
    """Test the get_supported_characters function."""

    def test_returns_string(self):
        """Should return a string."""
        chars = get_supported_characters()
        assert isinstance(chars, str)

    def test_contains_expected_characters(self):
        """Should contain letters, digits, and punctuation."""
        chars = get_supported_characters()
        for c in 'ABCXYZ0189-_.':
            assert c in chars


class TestUnknownCharacter:
    """Test handling of unknown characters."""

    def test_unknown_char_produces_placeholder(self):
        """Unknown characters should produce a placeholder polygon."""
        # Use a character not in the font
        polys = text_to_polygons('@', height=5.0)
        assert len(polys) > 0, "Should produce placeholder"

    def test_mixed_known_unknown(self):
        """Mixed known and unknown characters should work."""
        polys = text_to_polygons('A@B', height=5.0)
        # Should have polygons for all three characters
        assert len(polys) > 0
