"""Shading and -- the part that does the work -- the line pass.

A shaded surface says roughly which way a face points.  A drawn crease says
exactly where one face stops and the next begins, which is what makes a
chamfer, a gusset or a bolt hole legible at a few hundred pixels.  Every edge
is found on the composited buffer, not on any one part, so a silhouette
between two parts costs the same as one against the background -- the
property that lets a second backend join without reworking this file.
"""

import numpy as np
import pytest

from yapcad import sdf
from yapcad.render import palette
from yapcad.render.buffers import PartBuffer, composite
from yapcad.render.camera import Camera
from yapcad.render.shade import edges, legend_entries, base_colours, shade
from yapcad.render.trace import render_part

BLOCK = sdf.box((20.0, 20.0, 20.0))


def flat(value, normal=(0.0, -1.0, 0.0), shape=(16, 16), name="", group=""):
    depth = np.full(shape, value, dtype=float)
    normals = np.zeros(shape + (3,))
    normals[:] = normal
    return PartBuffer(depth=depth, normal=normals, name=name, group=group)


def camera(node=BLOCK, view="iso", pixels=(65, 65)):
    return Camera.fit(node.bounds, view=view, pixels=pixels, margin=1.3)


class TestEdges:

    def test_a_silhouette_is_drawn_against_the_background(self):
        buffer = flat(10.0)
        buffer.depth[:, 8:] = np.inf
        mask = edges(composite([buffer]), 0.5)
        assert mask[:, 7].all()
        assert not mask[:, 0].any()

    def test_a_crease_is_drawn_where_the_normal_turns(self):
        buffer = flat(10.0)
        buffer.normal[:, 8:] = (1.0, 0.0, 0.0)
        mask = edges(composite([buffer]), 0.5)
        assert mask[:, 7].all()

    def test_a_gentle_curve_is_not_a_crease(self):
        buffer = flat(10.0)
        angles = np.linspace(0.0, 0.25, 16)
        buffer.normal[:] = np.stack(
            [np.sin(angles), -np.cos(angles), np.zeros(16)], axis=1)
        assert not edges(composite([buffer]), 0.5).any()

    def test_a_depth_step_within_one_part_is_drawn(self):
        buffer = flat(10.0)
        buffer.depth[:, 8:] = 14.0
        mask = edges(composite([buffer]), 0.5)
        assert mask[:, 7].all()

    def test_a_shallow_slope_is_not_a_step(self):
        buffer = flat(10.0)
        buffer.depth[:] = 10.0 + 0.1 * np.arange(16)[None, :]
        assert not edges(composite([buffer]), 0.5).any()

    def test_the_boundary_between_two_parts_is_drawn(self):
        left, right = flat(10.0, name="L"), flat(10.0, name="R")
        left.depth[:, 8:] = np.inf
        right.depth[:, :8] = np.inf
        mask = edges(composite([left, right]), 0.5)
        # Two parts flush against each other still read as two parts.
        assert mask[:, 7].all() or mask[:, 8].all()

    def test_two_parts_meeting_flush_still_read_as_two(self):
        # Same normal, and a depth step far too small to be a step: only
        # knowing which part owns the pixel draws this line, which is why
        # the composited buffer carries a part index at all.
        back = flat(10.0, name="BACK")
        front_part = flat(9.9, name="FRONT")
        front_part.depth[:, :8] = np.inf
        gbuf = composite([back, front_part], tie=0.01)
        assert len(np.unique(gbuf.part)) == 2
        mask = edges(gbuf, 0.5)
        assert mask[:, 7].all()

    def test_a_flat_face_has_no_lines_in_it(self):
        assert not edges(composite([flat(10.0)]), 0.5).any()

    def test_a_real_box_gets_its_three_visible_creases(self):
        cam = camera()
        gbuf = composite([render_part(BLOCK, cam)])
        mask = edges(gbuf, cam.mm_per_px)
        # The silhouette is a hexagon and the three faces meet in a Y, so
        # there is appreciably more line than silhouette alone.
        assert mask.sum() > 150


class TestColourAssignment:

    def test_one_ungrouped_part_takes_the_neutral_material(self):
        gbuf = composite([flat(10.0, name="PLATE")])
        colours = base_colours(gbuf, palette.LIGHT)
        assert np.allclose(colours[0, 0], palette.LIGHT.material)

    def test_distinct_groups_take_distinct_slots(self):
        left = flat(10.0, name="A", group="ALUMINIUM")
        right = flat(12.0, name="B", group="PLA")
        right.depth[:, :8] = np.inf
        left.depth[:, 8:] = np.inf
        colours = base_colours(composite([left, right]), palette.LIGHT)
        assert np.allclose(colours[0, 0], palette.LIGHT.series[0])
        assert np.allclose(colours[0, 15], palette.LIGHT.series[1])

    def test_background_is_left_black_for_the_surface_to_cover(self):
        gbuf = composite([flat(np.inf)])
        assert np.all(base_colours(gbuf, palette.LIGHT) == 0.0)


class TestLegend:

    def test_a_single_material_needs_no_legend(self):
        gbuf = composite([flat(10.0, name="PLATE")])
        assert legend_entries(gbuf, "light") == []

    def test_each_group_gets_a_labelled_row(self):
        gbuf = composite([flat(10.0, group="ALUMINIUM"),
                          flat(12.0, group="PLA")])
        entries = legend_entries(gbuf, "light")
        assert [label for label, _ in entries] == ["ALUMINIUM", "PLA"]

    def test_a_repeated_group_appears_once(self):
        gbuf = composite([flat(10.0, group="PLA"), flat(12.0, group="PLA")])
        assert len(legend_entries(gbuf, "light")) == 1

    def test_past_eight_groups_the_rest_are_named_other(self):
        gbuf = composite([flat(10.0 + i, group=f"G{i}") for i in range(11)])
        entries = legend_entries(gbuf, "light")
        assert len(entries) == 9
        assert entries[-1][0] == "OTHER"
        assert entries[-1][1] == palette.LIGHT.other


class TestShade:

    def test_background_is_the_theme_surface(self):
        gbuf = composite([flat(np.inf)])
        rgb = shade(gbuf, camera(), theme="light")
        assert np.allclose(rgb, palette.LIGHT.surface)

    def test_a_lit_face_is_not_the_surface_colour(self):
        rgb = shade(composite([flat(10.0)]), camera(), theme="light")
        assert not np.allclose(rgb[8, 8], palette.LIGHT.surface)

    def test_a_face_turned_from_the_light_is_darker_than_one_toward_it(self):
        cam = camera()
        _, right, up = cam.basis
        toward = flat(10.0, normal=tuple(-0.9 * np.asarray(cam.direction)
                                         + 0.45 * up - 0.35 * right))
        away = flat(10.0, normal=tuple(np.asarray(cam.direction)))
        bright = shade(composite([toward]), cam, draw_edges=False)
        dim = shade(composite([away]), cam, draw_edges=False)
        assert bright[8, 8].sum() > dim[8, 8].sum()

    def test_nothing_is_ever_pitch_black(self):
        cam = camera()
        gbuf = composite([flat(10.0, normal=tuple(cam.direction))])
        rgb = shade(gbuf, cam, draw_edges=False)
        assert rgb[8, 8].min() > 0.0

    def test_edges_are_drawn_in_the_line_colour(self):
        buffer = flat(10.0)
        buffer.depth[:, 8:] = np.inf
        rgb = shade(composite([buffer]), camera(), theme="light")
        assert np.allclose(rgb[8, 7], palette.LIGHT.line)

    def test_edges_can_be_turned_off(self):
        buffer = flat(10.0)
        buffer.depth[:, 8:] = np.inf
        rgb = shade(composite([buffer]), camera(), theme="light",
                    draw_edges=False)
        assert not np.allclose(rgb[8, 7], palette.LIGHT.line)

    def test_the_dark_theme_is_darker(self):
        gbuf = composite([flat(np.inf)])
        light = shade(gbuf, camera(), theme="light")
        dark = shade(gbuf, camera(), theme="dark")
        assert dark.sum() < light.sum()

    def test_lighting_is_fixed_in_camera_space(self):
        # The same face, seen from two views, is lit the same way: that is
        # what lets a reader compare a feature between tiles of a sheet.
        front = camera(view="front")
        top = camera(view="top")
        facing = flat(10.0, normal=tuple(-np.asarray(front.direction)))
        up_facing = flat(10.0, normal=tuple(-np.asarray(top.direction)))
        a = shade(composite([facing]), front, draw_edges=False)[8, 8]
        b = shade(composite([up_facing]), top, draw_edges=False)[8, 8]
        assert np.allclose(a, b)

    def test_the_image_is_in_unit_range(self):
        cam = camera()
        rgb = shade(composite([render_part(BLOCK, cam)]), cam)
        assert rgb.min() >= 0.0 and rgb.max() <= 1.0
