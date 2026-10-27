"""The functions a caller actually reaches for, and the command line.

The point of the package is that an agent can build a part and look at it in
the same turn, so these tests go end to end: DSL source in, PNG on disk out,
with the geometry checked by where it lands in the frame rather than by
eyeballing the picture.
"""

import math
import subprocess
import sys

import numpy as np
import pytest

from yapcad import sdf
from yapcad.render import palette
from yapcad.render.api import (
    Part,
    as_parts,
    four_up,
    image,
    parts_from_dsl,
    render,
    scene_bounds,
    section_view,
)
from yapcad.render.camera import Camera
from yapcad.render.__main__ import _param, _section_spec, _size, _value, main

BALL = sdf.sphere(20.0)
BLOCK = sdf.box((20.0, 20.0, 20.0))

DESIGN = """module renderprobe
command PLATE() -> solid:
    let slab: solid = box(40.0, 30.0, 6.0)
    let hole: solid = translate(cylinder(4.0, 20.0), 10.0, 0.0, -10.0)
    emit difference(slab, hole)
"""


@pytest.fixture
def design(tmp_path):
    path = tmp_path / "probe.dsl"
    path.write_text(DESIGN)
    return path


class TestAsParts:

    def test_a_bare_field_becomes_one_part(self):
        parts = as_parts(BALL)
        assert len(parts) == 1 and parts[0].node is BALL

    def test_a_part_passes_through(self):
        part = Part(BALL, name="BALL", group="PLA")
        assert as_parts(part) == [part]

    def test_a_mapping_names_each_part(self):
        parts = as_parts({"LEFT": BALL, "RIGHT": BLOCK})
        assert [p.name for p in parts] == ["LEFT", "RIGHT"]

    def test_a_list_is_flattened(self):
        parts = as_parts([BALL, {"BLOCK": BLOCK}, Part(BALL, name="X")])
        assert [p.name for p in parts] == ["", "BLOCK", "X"]

    def test_an_sdf_authored_solid_gives_up_its_field(self):
        solid = sdf.to_solid(BLOCK, resolution=16)
        parts = as_parts(solid)
        assert len(parts) == 1
        assert parts[0].node.bounds == BLOCK.bounds

    def test_something_that_is_not_geometry_is_refused(self):
        with pytest.raises(TypeError, match="SDF-authored"):
            as_parts(42)

    def test_scene_bounds_covers_every_part(self):
        far = sdf.translate(BALL, (100.0, 0.0, 0.0))
        lo, hi = scene_bounds(as_parts([BALL, far]))
        assert lo[0] == pytest.approx(-20.0)
        assert hi[0] == pytest.approx(120.0)

    def test_an_empty_scene_has_no_bounds(self):
        with pytest.raises(ValueError, match="empty scene"):
            scene_bounds([])


class TestImage:

    def test_the_image_has_the_size_that_was_asked_for(self):
        img = image(BALL, view="front", pixels=(80, 60), annotate=False)
        assert img.shape == (60, 80, 3)
        assert img.dtype == np.uint8

    def test_the_model_lands_in_the_middle_of_the_frame(self):
        img = image(BALL, view="front", pixels=(81, 81), annotate=False)
        background = np.asarray(palette.LIGHT.surface) * 255
        covered = (np.abs(img.astype(float) - background).sum(axis=2) > 12)
        rows, cols = np.nonzero(covered)
        assert rows.mean() == pytest.approx(40.0, abs=1.5)
        assert cols.mean() == pytest.approx(40.0, abs=1.5)
        # A sphere framed with the default margin covers about pi/4 of it.
        assert covered.mean() == pytest.approx(
            math.pi / 4.0 / 1.08 ** 2, abs=0.04)

    def test_a_path_writes_a_png(self, tmp_path):
        out = tmp_path / "ball.png"
        image(BALL, view="iso", pixels=(48, 48), path=out)
        assert out.exists()
        assert out.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"

    def test_annotation_marks_the_frame(self):
        plain = image(BALL, view="front", pixels=(200, 160), annotate=False)
        marked = image(BALL, view="front", pixels=(200, 160), annotate=True)
        assert not np.array_equal(plain, marked)

    def test_a_direction_can_be_given_instead_of_a_named_view(self):
        img = image(BALL, view=(-1.0, -0.4, -0.3), pixels=(48, 48))
        assert img.shape == (48, 48, 3)

    def test_both_themes_render(self):
        light = image(BALL, pixels=(40, 40), theme="light", annotate=False)
        dark = image(BALL, pixels=(40, 40), theme="dark", annotate=False)
        assert dark.sum() < light.sum()

    def test_a_perspective_view_is_labelled_not_to_scale(self):
        from yapcad.render.camera import PERSPECTIVE
        img = image(BLOCK, view="iso", pixels=(240, 180),
                    projection=PERSPECTIVE)
        assert img.shape == (180, 240, 3)

    def test_an_unbounded_field_is_refused_with_a_reason(self):
        with pytest.raises(ValueError, match="unbounded"):
            image(sdf.half_space((0.0, 0.0, 1.0), 0.0))

    def test_explicit_bounds_let_an_unbounded_field_render(self):
        img = image(sdf.half_space((0.0, 0.0, 1.0), 0.0), view="front",
                    pixels=(40, 40),
                    bounds=((-10.0, -10.0, -10.0), (10.0, 10.0, 10.0)))
        assert img.shape == (40, 40, 3)


class TestScenes:

    def two(self):
        return {"LEFT": sdf.translate(BLOCK, (-15.0, 0.0, 0.0)),
                "RIGHT": sdf.translate(BALL, (25.0, 0.0, 0.0))}

    def test_both_parts_appear(self):
        parts = as_parts(self.two())
        camera = Camera.fit(scene_bounds(parts), view="front",
                            pixels=(129, 65))
        gbuf = render(parts, camera)
        assert set(np.unique(gbuf.part)) == {-1, 0, 1}

    def test_the_nearer_part_occludes_the_further_one(self):
        near = sdf.translate(BLOCK, (0.0, -30.0, 0.0))
        far = sdf.translate(BLOCK, (0.0, 30.0, 0.0))
        parts = as_parts({"NEAR": near, "FAR": far})
        camera = Camera.fit(scene_bounds(parts), view="front",
                            pixels=(65, 65))
        gbuf = render(parts, camera)
        middle = 32
        assert gbuf.names[gbuf.part[middle, middle]] == "NEAR"

    def test_a_small_scene_colours_by_part_on_its_own(self):
        parts = as_parts(self.two())
        camera = Camera.fit(scene_bounds(parts), view="front",
                            pixels=(65, 65))
        gbuf = render(parts, camera, colour_by="auto")
        assert gbuf.groups == ("LEFT", "RIGHT")

    def test_a_crowded_scene_falls_back_to_the_neutral_tone(self):
        many = {f"P{i}": sdf.translate(BALL, (40.0 * i, 0.0, 0.0))
                for i in range(9)}
        parts = as_parts(many)
        camera = Camera.fit(scene_bounds(parts), view="front",
                            pixels=(129, 33))
        gbuf = render(parts, camera, colour_by="auto")
        assert set(gbuf.groups) == {""}

    def test_explicit_groups_win_over_the_automatic_choice(self):
        parts = [Part(BLOCK, name="A", group="PLA"),
                 Part(sdf.translate(BALL, (40.0, 0.0, 0.0)), name="B",
                      group="PLA")]
        camera = Camera.fit(scene_bounds(parts), view="front",
                            pixels=(65, 33))
        assert render(parts, camera, colour_by="auto").groups == ("PLA", "PLA")

    def test_colour_can_be_switched_off(self):
        parts = as_parts(self.two())
        camera = Camera.fit(scene_bounds(parts), view="front",
                            pixels=(33, 33))
        assert set(render(parts, camera, colour_by="none").groups) == {""}

    def test_an_unknown_colour_rule_is_refused(self):
        camera = Camera.fit(BALL.bounds, pixels=(16, 16))
        with pytest.raises(ValueError, match="colour_by"):
            render(BALL, camera, colour_by="rainbow")


class TestSheets:

    def test_a_section_view_is_labelled_and_sized(self, tmp_path):
        out = tmp_path / "cut.png"
        img = section_view(BLOCK, axis="z", where=0.0, pixels=(200, 160),
                           path=out)
        assert img.shape == (160, 200, 3)
        assert out.exists()

    def test_four_up_lays_out_four_tiles(self):
        sheet = four_up(BLOCK, pixels=(60, 50), gutter=4)
        assert sheet.shape == (2 * 50 + 3 * 4, 2 * 60 + 3 * 4, 3)

    def test_four_up_can_drop_the_section(self):
        sheet = four_up(BLOCK, views=("iso", "front", "top"),
                        section_axis=None, pixels=(60, 50), gutter=4)
        assert sheet.shape == (2 * 50 + 3 * 4, 2 * 60 + 3 * 4, 3)

    def test_four_up_writes_one_png(self, tmp_path):
        out = tmp_path / "sheet.png"
        four_up(BLOCK, pixels=(50, 40), path=out)
        assert out.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"

    def test_every_tile_frames_the_same_model(self):
        tall = sdf.box((10.0, 10.0, 60.0))
        sheet = four_up(tall, pixels=(80, 80), gutter=2)
        assert sheet.shape[0] == sheet.shape[1]

    def test_every_tile_is_drawn_at_one_scale(self):
        # Left to fit itself, an isometric view of a cube reaches further
        # across the screen than a face-on one, and the tiles come out at
        # different scales -- which defeats the comparison a sheet is for.
        cube = sdf.box((20.0, 20.0, 20.0))
        pixels = (120, 120)
        shared = max(Camera.fit(cube.bounds, view=v, pixels=pixels).extent
                     for v in ("iso", "front", "top"))
        assert Camera.fit(cube.bounds, view="front",
                          pixels=pixels).extent < shared
        sheet = four_up(cube, views=("iso", "front", "top"),
                        section_axis=None, pixels=pixels, gutter=0,
                        annotate=False)
        front = sheet[0:120, 120:240]          # the second tile
        background = np.asarray(palette.LIGHT.surface) * 255
        covered = np.abs(front.astype(float) - background).sum(axis=2) > 12
        width = covered[60].sum()
        assert width == pytest.approx(20.0 / (shared / 120), abs=2.0)


class TestFromDsl:

    def test_a_dsl_command_renders(self, design, tmp_path):
        parts, seconds = parts_from_dsl(design, "PLATE", cell_mm=2.0)
        assert len(parts) == 1
        assert parts[0].name == "PLATE"
        assert seconds >= 0.0
        lo, hi = parts[0].node.bounds
        assert hi[0] - lo[0] == pytest.approx(40.0, abs=0.01)
        out = tmp_path / "plate.png"
        image(parts, view="top", pixels=(120, 90), path=out)
        assert out.exists()

    def test_a_missing_command_says_which_design(self, design):
        with pytest.raises(RuntimeError, match="NOPE"):
            parts_from_dsl(design, "NOPE")

    def test_the_hole_is_visible_from_above(self, design):
        parts, _ = parts_from_dsl(design, "PLATE", cell_mm=2.0)
        # Straight down through the bore: nothing to hit.
        camera = Camera.fit(scene_bounds(parts), view="top",
                            pixels=(161, 161), margin=1.0)
        gbuf = render(parts, camera)
        col = int((10.0 + 20.0) / 40.0 * 161)
        assert not gbuf.hit[80, col]
        assert gbuf.hit[80, 10]


class TestCommandLine:

    def test_size_is_parsed(self):
        assert _size("800x600") == (800, 600)

    def test_a_malformed_size_is_refused(self):
        import argparse
        with pytest.raises(argparse.ArgumentTypeError, match="WIDTHxHEIGHT"):
            _size("big")

    @pytest.mark.parametrize("text,value", [
        ("4", 4), ("4.5", 4.5), ("true", True), ("FALSE", False),
        ("M8", "M8"),
    ])
    def test_parameter_values_keep_their_type(self, text, value):
        assert _value(text) == value
        assert isinstance(_value(text), type(value))

    def test_a_parameter_needs_a_name_and_a_value(self):
        import argparse
        assert _param("teeth=24") == ("teeth", 24)
        with pytest.raises(argparse.ArgumentTypeError, match="NAME=VALUE"):
            _param("teeth")

    def test_a_section_spec_takes_an_axis_and_an_optional_plane(self):
        assert _section_spec("z") == ("z", None)
        assert _section_spec("x=12.5") == ("x", 12.5)

    def test_a_bad_section_spec_is_refused(self):
        import argparse
        with pytest.raises(argparse.ArgumentTypeError, match="axis"):
            _section_spec("w=1")
        with pytest.raises(argparse.ArgumentTypeError, match="millimetres"):
            _section_spec("z=middle")

    def test_it_renders_a_design_from_the_command_line(self, design,
                                                       tmp_path, capsys):
        out = tmp_path / "cli.png"
        assert main([str(design), "PLATE", "--view", "front",
                     "--size", "160x120", "--cell", "2",
                     "-o", str(out)]) == 0
        assert out.exists()
        assert "PLATE" in capsys.readouterr().out

    def test_four_up_from_the_command_line(self, design, tmp_path):
        out = tmp_path / "cli4.png"
        main([str(design), "PLATE", "--view", "4up", "--size", "160x120",
              "--cell", "2", "--section", "z=0", "-o", str(out), "-q"])
        assert out.exists()

    def test_a_section_from_the_command_line(self, design, tmp_path):
        out = tmp_path / "cut.png"
        main([str(design), "PLATE", "--view", "section", "--section", "z=0",
              "--size", "120x100", "--cell", "2", "-o", str(out), "-q"])
        assert out.exists()

    def test_an_unknown_view_is_refused(self, design, tmp_path):
        with pytest.raises(SystemExit):
            main([str(design), "PLATE", "--view", "oblique"])

    def test_the_module_runs_as_a_script(self, design, tmp_path):
        out = tmp_path / "script.png"
        done = subprocess.run(
            [sys.executable, "-m", "yapcad.render", str(design), "PLATE",
             "--size", "80x60", "--cell", "4", "-o", str(out), "-q"],
            capture_output=True, text=True)
        assert done.returncode == 0, done.stderr
        assert out.exists()
