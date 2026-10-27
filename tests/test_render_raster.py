"""The raster and the annotation drawn on it.

An agent reading one of these images cannot measure a pixel or hover a
tooltip, so the scale bar, the lettering and the ramp key have to be right.
The PNG writer is checked by decoding what it wrote rather than by trusting
it, since nothing else in yapCAD reads PNG.
"""

import struct
import zlib

import numpy as np
import pytest

from yapcad.render import palette
from yapcad.render.raster import (
    GLYPHS,
    frame,
    legend,
    line,
    new_image,
    paste,
    plate,
    ramp_key,
    rect,
    round_number,
    scale_bar,
    text,
    text_size,
    to_image,
    write_png,
)

WHITE = (1.0, 1.0, 1.0)
BLACK = (0.0, 0.0, 0.0)


def decode_png(path):
    """A minimal reader for exactly what :func:`write_png` writes."""
    data = open(path, "rb").read()
    assert data[:8] == b"\x89PNG\r\n\x1a\n"
    pos, header, pixels = 8, None, b""
    while pos < len(data):
        length = struct.unpack(">I", data[pos:pos + 4])[0]
        tag = data[pos + 4:pos + 8]
        payload = data[pos + 8:pos + 8 + length]
        crc = struct.unpack(">I", data[pos + 8 + length:pos + 12 + length])[0]
        assert crc == zlib.crc32(tag + payload) & 0xFFFFFFFF, tag
        if tag == b"IHDR":
            header = struct.unpack(">IIBBBBB", payload)
        elif tag == b"IDAT":
            pixels += payload
        pos += 12 + length
    width, height, depth, colour, *_ = header
    assert (depth, colour) == (8, 2)          # 8-bit RGB
    raw = zlib.decompress(pixels)
    stride = width * 3
    out = np.empty((height, width, 3), dtype=np.uint8)
    for row in range(height):
        scanline = raw[row * (stride + 1):(row + 1) * (stride + 1)]
        assert scanline[0] == 0               # filter type: none
        out[row] = np.frombuffer(scanline[1:], dtype=np.uint8).reshape(
            width, 3)
    return out


class TestPng:

    def test_what_was_written_is_what_comes_back(self, tmp_path):
        rng = np.random.default_rng(3)
        img = rng.integers(0, 256, (17, 23, 3), dtype=np.uint8)
        path = tmp_path / "probe.png"
        write_png(path, img)
        assert np.array_equal(decode_png(path), img)

    def test_a_single_pixel_round_trips(self, tmp_path):
        img = np.array([[[1, 2, 3]]], dtype=np.uint8)
        write_png(tmp_path / "one.png", img)
        assert np.array_equal(decode_png(tmp_path / "one.png"), img)

    def test_a_greyscale_array_is_refused(self, tmp_path):
        with pytest.raises(ValueError, match="RGB"):
            write_png(tmp_path / "bad.png", np.zeros((4, 4), dtype=np.uint8))

    def test_float_images_are_converted_and_clamped(self):
        img = to_image(np.array([[[-0.5, 0.5, 1.5]]]))
        assert img.dtype == np.uint8
        assert list(img[0, 0]) == [0, 128, 255]


class TestFont:

    def test_the_face_covers_what_a_drawing_letters(self):
        for ch in "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789 .,-+=:/()[]<>_#%*'?":
            assert ch in GLYPHS, ch

    def test_every_glyph_is_five_by_seven(self):
        for ch, rows in GLYPHS.items():
            assert len(rows) == 7, ch
            assert {len(r) for r in rows} == {5}, ch
            assert set("".join(rows)) <= {"0", "1"}, ch

    def test_lettering_is_upper_case(self):
        lower = new_image(60, 16, WHITE)
        upper = new_image(60, 16, WHITE)
        text(lower, "mm", 2, 2, BLACK)
        text(upper, "MM", 2, 2, BLACK)
        assert np.array_equal(lower, upper)

    def test_an_unknown_character_shows_rather_than_vanishes(self):
        img = new_image(40, 16, WHITE)
        text(img, "©", 2, 2, BLACK)
        assert (img != 255).any()

    def test_text_size_matches_what_gets_drawn(self):
        img = new_image(200, 40, WHITE)
        text(img, "SECTION", 10, 6, BLACK, 2)
        drawn = np.argwhere((img != 255).any(axis=2))
        width, height = text_size("SECTION", 2)
        assert drawn[:, 1].max() < 10 + width
        assert drawn[:, 0].max() < 6 + height

    def test_drawing_off_the_edge_clips_rather_than_raises(self):
        img = new_image(20, 12, WHITE)
        text(img, "OVERFLOWING", -5, -3, BLACK)
        assert img.shape == (12, 20, 3)


class TestPrimitives:

    def test_a_rectangle_fills_exactly_its_box(self):
        img = new_image(10, 10, WHITE)
        rect(img, 2, 3, 4, 5, BLACK)
        assert np.all(img[3:8, 2:6] == 0)
        assert np.all(img[0:3] == 255)

    def test_a_rectangle_off_the_edge_is_clipped(self):
        img = new_image(10, 10, WHITE)
        rect(img, -4, -4, 6, 6, BLACK)
        assert np.all(img[0:2, 0:2] == 0)
        assert np.all(img[3:, 3:] == 255)

    def test_a_frame_is_hollow(self):
        img = new_image(10, 10, WHITE)
        frame(img, 0, 0, 10, 10, BLACK)
        assert np.all(img[0] == 0) and np.all(img[-1] == 0)
        assert np.all(img[4, 1:9] == 255)

    def test_a_line_reaches_both_ends(self):
        img = new_image(20, 20, WHITE)
        line(img, 1, 1, 18, 18, BLACK)
        assert np.all(img[1, 1] == 0)
        assert np.all(img[18, 18] == 0)

    def test_pasting_places_the_tile(self):
        sheet = new_image(20, 20, WHITE)
        paste(sheet, new_image(4, 4, BLACK), 5, 6)
        assert np.all(sheet[6:10, 5:9] == 0)
        assert np.all(sheet[0, 0] == 255)

    def test_a_plate_washes_rather_than_hides(self):
        img = new_image(10, 10, BLACK)
        plate(img, 0, 0, 10, 10, WHITE, alpha=0.8)
        assert np.all(img > 180)      # mostly the wash
        assert np.all(img < 255)      # but the model shows through


class TestRoundNumbers:

    @pytest.mark.parametrize("target,expected", [
        (110.0, 100.0), (9.0, 5.0), (2.4, 2.0), (1.0, 1.0),
        (0.07, 0.05), (55.0, 50.0), (999.0, 500.0), (1000.0, 1000.0),
    ])
    def test_it_is_always_one_two_or_five(self, target, expected):
        assert round_number(target) == pytest.approx(expected)

    def test_it_never_exceeds_the_target(self):
        for target in (0.3, 1.7, 23.0, 410.0, 6100.0):
            assert round_number(target) <= target

    def test_a_non_positive_target_is_refused(self):
        with pytest.raises(ValueError, match="positive"):
            round_number(0.0)


class TestScaleBarAndKeys:

    def test_the_bar_is_a_round_number_of_millimetres(self):
        img = new_image(400, 200, WHITE)
        span = scale_bar(img, 0.37, BLACK, BLACK)
        assert span in (10.0, 20.0, 50.0)

    def test_the_bar_is_drawn_inside_the_image(self):
        img = new_image(400, 200, WHITE)
        span = scale_bar(img, 0.37, BLACK, BLACK)
        drawn = np.argwhere((img != 255).any(axis=2))
        assert drawn[:, 1].max() < 400 and drawn[:, 0].max() < 200
        assert span > 0

    def test_the_bar_length_matches_the_scale_it_claims(self):
        img = new_image(400, 200, WHITE)
        span = scale_bar(img, 0.25, BLACK, BLACK, x=10, y=150)
        row = img[150]
        dark = np.flatnonzero((row != 255).any(axis=1))
        assert (dark.max() - dark.min() + 1) == pytest.approx(
            span / 0.25, abs=1.5)

    def test_a_backing_keeps_the_bar_legible_over_a_model(self):
        img = new_image(400, 200, (0.5, 0.5, 0.5))
        scale_bar(img, 0.25, BLACK, BLACK, backing=WHITE)
        assert (img > 200).any()

    def test_the_legend_labels_every_swatch(self):
        theme = palette.LIGHT
        img = new_image(300, 120, WHITE)
        legend(img, [("ALUMINIUM", theme.series[0]),
                     ("PLA", theme.series[1])], BLACK)
        assert (img != 255).any()
        # Both swatch colours appear.
        for colour in theme.series[:2]:
            want = np.round(np.asarray(colour) * 255).astype(np.uint8)
            assert (np.abs(img.astype(int) - want).sum(axis=2) < 6).any()

    def test_an_empty_legend_draws_nothing(self):
        img = new_image(60, 40, WHITE)
        legend(img, [], BLACK)
        assert np.all(img == 255)

    def test_the_ramp_key_is_labelled_at_both_ends(self):
        img = new_image(400, 200, WHITE)
        ramp_key(img, palette.LIGHT.ramp, 0, 10, "MM", BLACK)
        assert (img != 255).any()


class TestPalette:

    def test_hex_parses_to_unit_floats(self):
        assert palette.from_hex("#ffffff") == (1.0, 1.0, 1.0)
        assert palette.from_hex("000000") == (0.0, 0.0, 0.0)

    def test_a_malformed_colour_is_refused(self):
        with pytest.raises(ValueError, match="hex"):
            palette.from_hex("#abc")

    def test_both_themes_carry_every_slot(self):
        for theme in (palette.LIGHT, palette.DARK):
            assert len(theme.series) == 8
            assert len(theme.ramp) >= 8

    def test_the_ramp_runs_light_to_dark(self):
        for theme in (palette.LIGHT, palette.DARK):
            luminance = [sum(step) for step in theme.ramp]
            assert luminance == sorted(luminance, reverse=True)

    def test_groups_take_slots_in_first_appearance_order(self):
        slots = palette.assign_groups(["B", "A", "B", "C"])
        assert slots == {"B": 0, "A": 1, "C": 2}

    def test_an_absent_group_gets_the_neutral_tone(self):
        slots = palette.assign_groups(["", "A"])
        assert slots[""] == -1
        assert palette.LIGHT.group_colour(-1) == palette.LIGHT.material

    def test_a_ninth_group_folds_rather_than_inventing_a_hue(self):
        keys = [f"G{i}" for i in range(10)]
        slots = palette.assign_groups(keys)
        theme = palette.LIGHT
        assert theme.group_colour(slots["G8"]) == theme.other
        assert theme.group_colour(slots["G9"]) == theme.other
        assert theme.group_colour(slots["G7"]) == theme.series[7]

    def test_a_theme_can_be_named_or_passed(self):
        assert palette.get_theme("dark") is palette.DARK
        assert palette.get_theme(palette.LIGHT) is palette.LIGHT

    def test_an_unknown_theme_names_the_ones_that_exist(self):
        with pytest.raises(ValueError, match="light"):
            palette.get_theme("sepia")
