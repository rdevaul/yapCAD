"""The camera, which both render backends have to agree about.

Everything downstream -- the scale bar's honesty, where the axis triad points,
whether a section lines up with the plan view it was cut from -- rests on this
module getting the screen basis right, so the conventions are pinned here
rather than left to be noticed in a picture.
"""

import math

import numpy as np
import pytest

from yapcad import sdf
from yapcad.render.camera import (
    ORTHOGRAPHIC,
    PERSPECTIVE,
    VIEWS,
    Camera,
)

BALL = sdf.sphere(20.0)


def fit(view="front", pixels=(65, 65), **kwargs):
    return Camera.fit(BALL.bounds, view=view, pixels=pixels, margin=1.0,
                      **kwargs)


class TestScreenBasis:
    """Drawing conventions: which way is right, and which way is up."""

    def test_the_front_view_puts_x_right_and_z_up(self):
        _, right, up = fit("front").basis
        assert np.allclose(right, (1.0, 0.0, 0.0))
        assert np.allclose(up, (0.0, 0.0, 1.0))

    def test_the_top_view_puts_x_right_and_y_up(self):
        _, right, up = fit("top").basis
        assert np.allclose(right, (1.0, 0.0, 0.0))
        assert np.allclose(up, (0.0, 1.0, 0.0))

    def test_the_basis_is_right_handed(self):
        for view in VIEWS:
            direction, right, up = fit(view).basis
            assert np.allclose(np.cross(right, up), -direction), view

    def test_every_named_view_is_orthonormal(self):
        for view in VIEWS:
            for vector in fit(view).basis:
                assert math.isclose(float(np.linalg.norm(vector)), 1.0,
                                    rel_tol=1e-12)

    def test_a_skewed_up_vector_is_re_orthogonalised(self):
        camera = Camera(target=(0.0, 0.0, 0.0), direction=(0.0, 1.0, 0.0),
                        up=(0.0, 0.4, 1.0), extent=10.0, pixels=(8, 8),
                        distance=20.0, far=40.0)
        direction, right, up = camera.basis
        assert abs(float(np.dot(direction, up))) < 1e-12
        assert abs(float(np.dot(right, up))) < 1e-12

    def test_an_up_hint_along_the_view_direction_falls_back(self):
        # Looking straight down with "up" also straight down has no answer
        # in it; the camera has to pick one rather than divide by zero.
        camera = Camera.fit(BALL.bounds, direction=(0.0, 0.0, -1.0),
                            up=(0.0, 0.0, 1.0), pixels=(8, 8))
        for vector in camera.basis:
            assert np.all(np.isfinite(vector))


class TestFraming:
    """``fit`` frames the box itself, not a sphere around it."""

    def test_a_tight_fit_is_exactly_the_model(self):
        assert fit("front").extent == pytest.approx(40.0)

    def test_the_margin_is_a_fraction_of_the_model(self):
        loose = Camera.fit(BALL.bounds, view="front", pixels=(65, 65),
                           margin=1.25)
        assert loose.extent == pytest.approx(50.0)

    def test_a_wide_frame_still_fits_the_whole_model(self):
        # Twice as wide as tall: the height is what binds for a cube.
        wide = Camera.fit(BALL.bounds, view="front", pixels=(128, 64),
                          margin=1.0)
        assert wide.extent == pytest.approx(40.0)
        assert wide.extent * wide.aspect == pytest.approx(80.0)

    def test_a_tall_frame_widens_the_extent_to_fit(self):
        tall = Camera.fit(BALL.bounds, view="front", pixels=(64, 128),
                          margin=1.0)
        assert tall.extent * tall.aspect == pytest.approx(40.0)

    def test_millimetres_per_pixel_follows_the_extent(self):
        camera = fit("front", pixels=(80, 40))
        assert camera.mm_per_px == pytest.approx(camera.extent / 40)

    def test_the_rays_start_clear_of_the_model(self):
        for view in VIEWS:
            camera = fit(view)
            assert camera.distance > 20.0, view

    def test_a_flat_box_still_frames(self):
        flat = Camera.fit(((0.0, 0.0, 0.0), (10.0, 10.0, 0.0)), view="top",
                          pixels=(16, 16))
        assert flat.extent > 0.0

    def test_an_unbounded_field_is_refused_with_a_reason(self):
        with pytest.raises(ValueError, match="unbounded"):
            Camera.fit(sdf.half_space((0.0, 0.0, 1.0), 0.0).bounds)

    def test_inverted_bounds_are_refused(self):
        with pytest.raises(ValueError, match="inverted"):
            Camera.fit(((10.0, 0.0, 0.0), (0.0, 1.0, 1.0)))

    def test_an_unknown_view_names_the_ones_that_exist(self):
        with pytest.raises(ValueError, match="iso"):
            Camera.fit(BALL.bounds, view="oblique")


class TestRays:

    def test_an_orthographic_view_has_one_direction(self):
        _, dirs = fit("front").rays()
        assert np.allclose(dirs, (0.0, 1.0, 0.0))

    def test_orthographic_origins_span_the_extent(self):
        camera = fit("front", pixels=(65, 65))
        origins, _ = camera.rays()
        half = 0.5 * camera.extent - 0.5 * camera.mm_per_px
        assert origins[:, 0].max() == pytest.approx(half)
        assert origins[:, 0].min() == pytest.approx(-half)

    def test_rays_run_row_major_from_the_top_left(self):
        camera = fit("front", pixels=(9, 9))
        origins, _ = camera.rays()
        assert origins[0][0] < 0.0 and origins[0][2] > 0.0     # top left
        assert origins[-1][0] > 0.0 and origins[-1][2] < 0.0   # bottom right

    def test_a_perspective_view_shares_one_eye(self):
        camera = fit("front", projection=PERSPECTIVE)
        origins, dirs = camera.rays()
        assert np.allclose(origins, origins[0])
        assert np.allclose(np.linalg.norm(dirs, axis=1), 1.0)

    def test_the_perspective_frame_holds_the_model(self):
        camera = Camera.fit(BALL.bounds, view="front", pixels=(65, 65),
                            projection=PERSPECTIVE, fov=32.0)
        assert camera.extent >= 40.0

    def test_one_ray_per_pixel(self):
        camera = fit("iso", pixels=(11, 7))
        origins, dirs = camera.rays()
        assert origins.shape == (77, 3) == dirs.shape


class TestProjection:

    def test_the_target_lands_in_the_middle(self):
        camera = fit("iso", pixels=(64, 48))
        col, row = camera.project([camera.target])[0]
        assert col == pytest.approx(31.5)
        assert row == pytest.approx(23.5)

    def test_moving_right_moves_right(self):
        camera = fit("front", pixels=(64, 64))
        _, right, up = camera.basis
        here = camera.project([camera.target])[0]
        there = camera.project([np.asarray(camera.target) + right * 5.0])[0]
        assert there[0] > here[0]
        above = camera.project([np.asarray(camera.target) + up * 5.0])[0]
        assert above[1] < here[1]


class TestValidation:

    def test_a_bad_projection_is_named(self):
        with pytest.raises(ValueError, match="orthographic"):
            Camera(target=(0, 0, 0), direction=(0, 1, 0), up=(0, 0, 1),
                   extent=1.0, projection="fisheye")

    def test_an_infinite_far_is_refused(self):
        with pytest.raises(ValueError, match="unbounded field"):
            Camera(target=(0, 0, 0), direction=(0, 1, 0), up=(0, 0, 1),
                   extent=1.0, far=float("inf"))

    def test_a_zero_extent_is_refused(self):
        with pytest.raises(ValueError, match="extent"):
            Camera(target=(0, 0, 0), direction=(0, 1, 0), up=(0, 0, 1),
                   extent=0.0)

    def test_a_degenerate_frame_is_refused(self):
        with pytest.raises(ValueError, match="2x2"):
            Camera(target=(0, 0, 0), direction=(0, 1, 0), up=(0, 0, 1),
                   extent=1.0, pixels=(1, 10))

    def test_orthographic_is_the_default(self):
        assert fit("iso").projection == ORTHOGRAPHIC
        assert fit("iso").is_orthographic
