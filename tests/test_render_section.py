"""Sections, which are one evaluation of the field over a plane.

No marching, no meshing, no boolean: the slice *is* the field restricted to a
plane.  That makes it both the cheapest view and the only one that carries
information a mesh renderer cannot -- outside the material the value is the
clearance to the nearest surface, and inside it is the local half-thickness.
"""

import math

import numpy as np
import pytest

from yapcad import sdf
from yapcad.render import palette
from yapcad.render.section import PLANES, section, section_colours, section_image

BALL = sdf.sphere(20.0)
BLOCK = sdf.box((40.0, 30.0, 20.0))


class TestSampling:

    def test_the_centre_of_a_sphere_is_one_radius_inside(self):
        cut = section(BALL, axis="z", where=0.0, pixels=(129, 129))
        assert cut.distance[64, 64] == pytest.approx(-20.0, abs=1e-9)

    def test_the_centre_of_a_box_is_its_nearest_face_away(self):
        cut = section(BLOCK, axis="z", where=0.0, pixels=(129, 129))
        assert cut.distance[64, 64] == pytest.approx(-10.0, abs=1e-9)

    def test_the_cut_area_is_the_material_the_plane_passes_through(self):
        cut = section(BALL, axis="z", where=0.0, pixels=(257, 257))
        assert cut.area_mm2 == pytest.approx(math.pi * 400.0, rel=0.01)

    def test_an_off_centre_plane_cuts_a_smaller_circle(self):
        cut = section(BALL, axis="z", where=16.0, pixels=(257, 257))
        assert cut.area_mm2 == pytest.approx(math.pi * (400.0 - 256.0),
                                             rel=0.02)

    def test_the_default_plane_is_the_middle_of_the_model(self):
        shifted = sdf.translate(BALL, (0.0, 0.0, 50.0))
        assert section(shifted, axis="z").where == pytest.approx(50.0)

    def test_a_plane_that_misses_the_model_warns_and_is_empty(self):
        with pytest.warns(RuntimeWarning, match="outside the model"):
            cut = section(BALL, axis="z", where=90.0, pixels=(33, 33))
        assert cut.is_empty()
        assert cut.area_mm2 == 0.0

    def test_outside_the_material_the_value_is_the_clearance(self):
        cut = section(BALL, axis="z", where=0.0, pixels=(129, 129),
                      bounds=((-40.0, -40.0, -40.0), (40.0, 40.0, 40.0)))
        corner = cut.distance[0, 0]
        assert corner > 0.0
        # The far corner of the frame, measured from the sphere.
        (u0, v0), _ = cut.span
        reach = math.hypot(u0 + 0.5 * cut.mm_per_px, v0 + 0.5 * cut.mm_per_px)
        assert corner == pytest.approx(reach - 20.0, abs=1e-6)

    def test_each_axis_gets_the_screen_orientation_of_its_view(self):
        assert section(BLOCK, axis="z").labels == ("+X", "+Y")
        assert section(BLOCK, axis="y").labels == ("+X", "+Z")
        assert section(BLOCK, axis="x").labels == ("+Y", "+Z")

    def test_a_section_normal_to_x_cuts_the_right_rectangle(self):
        cut = section(BLOCK, axis="x", where=0.0, pixels=(257, 257))
        assert cut.area_mm2 == pytest.approx(30.0 * 20.0, rel=0.02)

    def test_the_slice_fills_the_frame_it_was_asked_for(self):
        cut = section(BLOCK, axis="z", pixels=(80, 40))
        assert cut.shape == (40, 80)
        assert cut.mm_per_px == pytest.approx(
            (cut.span[1][1] - cut.span[0][1]) / 40)

    def test_a_bad_axis_names_the_ones_that_exist(self):
        with pytest.raises(ValueError, match="x, y, z"):
            section(BALL, axis="w")

    def test_an_unbounded_field_is_refused_with_a_reason(self):
        with pytest.raises(ValueError, match="unbounded"):
            section(sdf.half_space((0.0, 0.0, 1.0), 0.0))

    def test_every_plane_is_declared(self):
        assert set(PLANES) == {"x", "y", "z"}


class TestColouring:

    def cut(self):
        return section(BALL, axis="z", where=0.0, pixels=(129, 129),
                       bounds=((-40.0, -40.0, -40.0), (40.0, 40.0, 40.0)))

    def test_material_is_flat_and_clearance_is_graded(self):
        theme = palette.LIGHT
        rgb, limit = section_colours(self.cut(), theme=theme)
        inside = np.allclose(rgb[64, 64], theme.material)
        assert inside
        assert limit > 0.0
        # The ramp actually ramps: the far corner is lighter than the wall.
        wall = rgb[64, 20]
        corner = rgb[0, 0]
        assert corner.sum() > wall.sum()

    def test_the_zero_contour_is_drawn_in_ink(self):
        theme = palette.LIGHT
        rgb, _ = section_colours(self.cut(), theme=theme)
        assert np.isclose(rgb, theme.line).all(axis=2).any()

    def test_the_limit_is_a_round_number(self):
        _, limit = section_colours(self.cut())
        assert limit in (0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0)

    def test_an_explicit_limit_is_honoured(self):
        _, limit = section_colours(self.cut(), limit=3.0)
        assert limit == 3.0

    def test_an_empty_plane_still_colours(self):
        with pytest.warns(RuntimeWarning):
            cut = section(BALL, axis="z", where=90.0, pixels=(17, 17))
        rgb, limit = section_colours(cut)
        assert rgb.shape == (17, 17, 3)
        assert limit > 0.0

    def test_the_image_is_eight_bit_rgb(self):
        img = section_image(self.cut())
        assert img.dtype == np.uint8
        assert img.shape == (129, 129, 3)

    def test_both_themes_colour_a_section(self):
        for theme in ("light", "dark"):
            rgb, _ = section_colours(self.cut(), theme=theme)
            assert rgb.shape == (129, 129, 3)
