"""Sphere tracing, checked against geometry whose answer is known in closed
form.

A renderer is easy to get wrong in ways that still look plausible, so these
tests ask for numbers rather than pictures: the depth of a sphere's pole, the
exact plane a box's front face sits on, the direction of a normal.  The
Lipschitz step rule is the load-bearing claim -- that a march guided by a
*lower bound* on distance cannot step through a surface -- and it is checked
on a hard CSG tree, where the bound is not tight.
"""

import math
import warnings

import numpy as np
import pytest

from yapcad import sdf
from yapcad.render.buffers import PartBuffer, composite
from yapcad.render.camera import Camera
from yapcad.render.trace import march, render_part, slab

BALL = sdf.sphere(20.0)
BLOCK = sdf.box((20.0, 20.0, 20.0))


def front(node=BALL, pixels=(65, 65), margin=1.0):
    return Camera.fit(node.bounds, view="front", pixels=pixels, margin=margin)


class TestSlab:
    """Ray/box clipping, including the cases a naive division turns to nan."""

    def box(self):
        return ((-1.0, -1.0, -1.0), (1.0, 1.0, 1.0))

    def test_a_ray_through_the_middle_enters_and_leaves(self):
        enter, leave = slab(np.array([[0.0, -5.0, 0.0]]),
                            np.array([[0.0, 1.0, 0.0]]), self.box())
        assert enter[0] == pytest.approx(4.0)
        assert leave[0] == pytest.approx(6.0)

    def test_a_ray_that_starts_inside_enters_at_zero(self):
        enter, leave = slab(np.array([[0.0, 0.0, 0.0]]),
                            np.array([[0.0, 1.0, 0.0]]), self.box())
        assert enter[0] == 0.0
        assert leave[0] == pytest.approx(1.0)

    def test_an_axis_parallel_ray_outside_the_slab_misses(self):
        # Travelling along +y at z = 9: no division can tell you this, the
        # z component of the direction is exactly zero.
        enter, leave = slab(np.array([[0.0, -5.0, 9.0]]),
                            np.array([[0.0, 1.0, 0.0]]), self.box())
        assert enter[0] >= leave[0]

    def test_an_axis_parallel_ray_inside_the_slab_hits(self):
        enter, leave = slab(np.array([[0.0, -5.0, 0.5]]),
                            np.array([[0.0, 1.0, 0.0]]), self.box())
        assert enter[0] < leave[0]

    def test_a_ray_pointing_away_is_clipped_to_the_present(self):
        enter, _ = slab(np.array([[0.0, 5.0, 0.0]]),
                        np.array([[0.0, 1.0, 0.0]]), self.box())
        assert enter[0] == 0.0

    def test_infinite_bounds_never_end(self):
        unbounded = ((-math.inf,) * 3, (math.inf,) * 3)
        enter, leave = slab(np.array([[0.0, -5.0, 0.0]]),
                            np.array([[0.0, 1.0, 0.0]]), unbounded)
        assert enter[0] == 0.0 and leave[0] == math.inf


class TestAgainstClosedForm:

    def test_a_sphere_pole_is_exactly_one_radius_from_the_frame(self):
        camera = front()
        buffer = render_part(BALL, camera)
        middle = camera.height // 2
        assert buffer.depth[middle, middle] == pytest.approx(
            camera.distance - 20.0, abs=1e-9)

    def test_every_hit_lands_on_the_surface(self):
        # The invariant the march actually promises: it stops where the
        # field is within the hit tolerance of zero, so every point it
        # reports is on the sphere to that tolerance.
        camera = front(pixels=(97, 97))
        buffer = render_part(BALL, camera)
        origins, dirs = camera.rays()
        depth = buffer.depth.ravel()
        hit = np.isfinite(depth)
        points = origins[hit] + depth[hit, None] * dirs[hit]
        radius = np.linalg.norm(points, axis=1)
        assert np.all(np.abs(radius - 20.0) <= buffer.stats["tolerance"])

    def test_a_grazing_ray_stops_short_rather_than_overshooting(self):
        # Near the limb the field changes slowly along the ray, so a hit is
        # declared well before the true depth -- by much more than the
        # tolerance, which bounds the *field*, not the depth.  It errs
        # toward the camera, never through the surface.
        camera = front(pixels=(97, 97))
        buffer = render_part(BALL, camera)
        origins, _ = camera.rays()
        depth = buffer.depth.ravel()
        hit = np.isfinite(depth)
        radial = origins[hit][:, 0] ** 2 + origins[hit][:, 2] ** 2
        inside = radial < 400.0
        truth = camera.distance - np.sqrt(400.0 - radial[inside])
        shortfall = truth - depth[hit][inside]
        assert np.all(shortfall >= -1e-9)
        assert shortfall.max() > buffer.stats["tolerance"]
        # Away from the limb it stays within a pixel or so.
        core = radial[inside] < 100.0
        assert np.all(shortfall[core] <= 2.0 * buffer.stats["tolerance"])

    def test_a_box_face_lands_exactly_on_its_plane(self):
        camera = front(BLOCK)
        buffer = render_part(BLOCK, camera)
        inner = buffer.depth[8:-8, 8:-8]
        assert np.all(np.isfinite(inner))
        assert np.allclose(inner, camera.distance - 10.0, atol=1e-9)

    def test_a_missed_ray_reports_infinite_depth(self):
        camera = Camera.fit(BALL.bounds, view="front", pixels=(65, 65),
                            margin=2.0)
        buffer = render_part(BALL, camera)
        assert buffer.depth[0, 0] == math.inf
        assert not buffer.hit[0, 0]

    def test_the_covered_area_is_the_silhouette(self):
        camera = front(pixels=(129, 129))
        buffer = render_part(BALL, camera)
        # pi r^2 over the framed square, within one ring of edge pixels.
        assert buffer.coverage == pytest.approx(math.pi / 4.0, abs=0.02)

    def test_normals_face_the_camera_at_the_pole(self):
        camera = front()
        buffer = render_part(BALL, camera)
        middle = camera.height // 2
        assert np.allclose(buffer.normal[middle, middle], (0.0, -1.0, 0.0),
                           atol=1e-6)

    def test_normals_are_radial_on_a_sphere(self):
        camera = front(pixels=(65, 65))
        buffer = render_part(BALL, camera)
        origins, dirs = camera.rays()
        depth = buffer.depth.ravel()
        hit = np.isfinite(depth)
        points = origins[hit] + depth[hit, None] * dirs[hit]
        truth = points / np.linalg.norm(points, axis=1)[:, None]
        got = buffer.normal.reshape(-1, 3)[hit]
        # Away from the limb, where the finite difference straddles the edge.
        steep = np.abs(got[:, 1]) > 0.3
        assert np.allclose(got[steep], truth[steep], atol=5e-3)


class TestHardCsg:
    """The conservative-bound claim, on a tree whose field is not exact."""

    def test_a_subtraction_does_not_leak_through_its_wall(self):
        # A 2 mm wall left by a subtraction: thin enough that a march which
        # overshot would sail straight through it and report the far face.
        shell = sdf.subtract(sdf.box((40.0, 40.0, 40.0)),
                             sdf.box((36.0, 36.0, 36.0)))
        assert not shell.exact
        camera = Camera.fit(shell.bounds, view="front", pixels=(65, 65),
                            margin=1.0)
        buffer = render_part(shell, camera)
        wall = buffer.depth[2:-2, 2:-2]
        assert np.all(np.isfinite(wall))
        assert np.allclose(wall, camera.distance - 20.0, atol=1e-9)

    def test_the_far_wall_is_found_through_the_opening(self):
        tube = sdf.subtract(sdf.box((40.0, 40.0, 40.0)),
                            sdf.cylinder(8.0, 60.0))
        camera = Camera.fit(tube.bounds, view="top", pixels=(65, 65),
                            margin=1.0)
        buffer = render_part(tube, camera)
        middle = camera.height // 2
        assert buffer.depth[middle, middle] == math.inf


class TestMarchMechanics:

    def test_clipping_skips_rays_that_miss_the_box(self):
        camera = Camera.fit(BALL.bounds, view="front", pixels=(65, 65),
                            margin=3.0)
        buffer = render_part(BALL, camera)
        assert buffer.stats["clipped"] > 0.5 * buffer.stats["rays"]

    def test_a_frame_with_nothing_in_it_marches_nothing(self):
        # Looking at empty space beside the model: every ray is clipped and
        # the field is never evaluated.
        camera = Camera.fit(((100.0, 100.0, 100.0), (140.0, 140.0, 140.0)),
                            view="front", pixels=(33, 33))
        buffer = render_part(BALL, camera)
        assert buffer.stats["clipped"] == buffer.stats["rays"]
        assert buffer.stats["steps"] == 0
        assert buffer.stats["evaluations"] == 0
        assert not buffer.hit.any()

    def test_chunking_does_not_change_the_answer(self):
        camera = front(pixels=(33, 33))
        whole = render_part(BALL, camera)
        split = render_part(BALL, camera, chunk=64)
        assert np.array_equal(whole.depth, split.depth)
        assert np.allclose(whole.normal, split.normal)

    def test_stats_account_for_every_ray(self):
        camera = front(pixels=(33, 33))
        stats = render_part(BALL, camera).stats
        assert stats["rays"] == 33 * 33
        assert stats["evaluations"] > stats["rays"]
        assert stats["unconverged"] == 0

    def test_running_out_of_steps_warns_rather_than_lying_quietly(self):
        camera = front(pixels=(33, 33))
        with pytest.warns(RuntimeWarning, match="did not converge"):
            render_part(BALL, camera, max_steps=1, name="BALL")

    def test_a_field_with_no_usable_lipschitz_bound_says_so(self, monkeypatch):
        monkeypatch.setattr(type(BALL), "lipschitz",
                            property(lambda self: math.inf))
        camera = front(pixels=(17, 17))
        with pytest.warns(RuntimeWarning, match="Lipschitz"):
            march(BALL, *camera.rays(), BALL.bounds, camera.far, 0.1)

    def test_an_unbounded_field_is_clipped_to_the_cameras_reach(self):
        ground = sdf.half_space((0.0, 0.0, 1.0), 0.0)
        camera = Camera.fit(((-20.0, -20.0, -20.0), (20.0, 20.0, 20.0)),
                            view="front", pixels=(17, 17))
        buffer = render_part(ground, camera)
        assert buffer.hit.any()
        assert np.all(buffer.depth[buffer.hit] <= camera.far)


class TestComposite:

    def flat(self, value, normal=(0.0, -1.0, 0.0), name="", group=""):
        depth = np.full((8, 8), value, dtype=float)
        normals = np.zeros((8, 8, 3))
        normals[:] = normal
        return PartBuffer(depth=depth, normal=normals, name=name, group=group)

    def test_the_nearest_part_wins(self):
        gbuf = composite([self.flat(10.0, name="FAR"),
                          self.flat(4.0, name="NEAR")])
        assert np.all(gbuf.part == 1)
        assert np.all(gbuf.depth == 4.0)

    def test_background_is_marked_and_not_shaded(self):
        gbuf = composite([self.flat(math.inf), self.flat(math.inf)])
        assert np.all(gbuf.part == -1)
        assert np.all(np.isinf(gbuf.depth))
        assert np.all(gbuf.normal == 0.0)

    def test_a_shared_face_resolves_to_one_part_rather_than_speckle(self):
        # Two datum-mated parts meeting on a plane, with the float noise a
        # real transform chain leaves behind.
        rng = np.random.default_rng(7)
        first = self.flat(25.0, name="A")
        second = self.flat(25.0, name="B")
        second.depth = second.depth + rng.normal(0.0, 1e-9, second.shape)
        speckled = composite([first, second], tie=0.0)
        assert len(np.unique(speckled.part)) == 2      # the failure mode
        clean = composite([first, second], tie=0.05)
        assert np.all(clean.part == 0)

    def test_quantising_does_not_reorder_parts_that_really_differ(self):
        gbuf = composite([self.flat(25.0), self.flat(24.0)], tie=0.05)
        assert np.all(gbuf.part == 1)

    def test_normals_come_from_the_winner(self):
        gbuf = composite([self.flat(10.0, normal=(1.0, 0.0, 0.0)),
                          self.flat(4.0, normal=(0.0, 0.0, 1.0))])
        assert np.allclose(gbuf.normal, (0.0, 0.0, 1.0))

    def test_names_and_groups_survive(self):
        gbuf = composite([self.flat(4.0, name="PLATE", group="ALUMINIUM")])
        assert gbuf.names == ("PLATE",)
        assert gbuf.groups == ("ALUMINIUM",)

    def test_depth_range_reports_the_covered_pixels(self):
        near = self.flat(4.0)
        near.depth[:4] = math.inf
        assert composite([near]).depth_range() == (4.0, 4.0)

    def test_an_empty_scene_is_refused(self):
        with pytest.raises(ValueError, match="at least one"):
            composite([])

    def test_mismatched_buffers_are_refused(self):
        odd = PartBuffer(depth=np.zeros((4, 4)), normal=np.zeros((4, 4, 3)))
        with pytest.raises(ValueError, match="image size"):
            composite([self.flat(1.0), odd])

    def test_a_buffer_must_be_self_consistent(self):
        with pytest.raises(ValueError, match="disagree"):
            PartBuffer(depth=np.zeros((4, 4)), normal=np.zeros((8, 8, 3)))
