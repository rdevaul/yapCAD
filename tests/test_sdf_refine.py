"""Meshing a field finer when the cell asked for cannot find the surface.

A boolean can leave a feature thinner than anything either operand had: the
yapRover differential cradle cuts a 12 mm gusset flat 2 mm above its base,
and the 2 mm remnant is real material that a 2 mm cell steps straight over.
Dual contouring then returned an empty mesh, and the empty mesh reached
``geom3d.surface`` as an empty vertex list, four frames down, as
``ValueError: bad arguments to surface`` -- a part silently lost, reported
as a type error.

``to_solid`` now halves the cell and tries again, says in a warning what it
had to do and which cell worked, and takes that cell back as ``cell_mm`` so
the next build pays for one pass.  A region the Lipschitz bound proves holds
nothing is not refined at all: that is an empty CSG result, and it meshes to
the empty solid yapCAD already allows for one.
"""

import math
import warnings

import pytest

from yapcad import sdf
from yapcad.construction import get_construction
from yapcad.geom3d import issolidclosed, volumeof
from yapcad.sdf.contour import proves_empty


def meshing_of(solid):
    return sdf.meshing_from_construction(get_construction(solid))


def resolution_of(solid):
    return meshing_of(solid)["resolution"]


# ---------------------------------------------------------------------------
# The case this is for: the yapRover differential cradle, as a field
# ---------------------------------------------------------------------------


def wedge_y(width, h, cx, y0, z0, direction):
    """A 45 degree gusset: the design's own WEDGE_Y, as a field."""
    bounds = sdf.translate(sdf.box((width, h, h)),
                           (cx, y0 + direction * h / 2.0, z0 + h / 2.0))
    cube = sdf.rotate(sdf.box((width + 2.0, 4.0 * h, 4.0 * h)),
                      (1, 0, 0), 45.0)
    keep = sdf.translate(cube, (cx, y0 - direction * 1.41421356 * h,
                                z0 + h - 1.41421356 * h))
    return sdf.intersect(bounds, keep)


#: Two gussets cut flat 2 mm above the base they stand on.  The remnant is
#: a 20 x 2 mm trapezoidal sliver of each: 440 mm**3 apiece, exactly.
CRADLE_GUSSETS = sdf.subtract(
    sdf.union(sdf.rotate(wedge_y(20.0, 12.0, 0.0, -36.0, -11.0, -1.0),
                         (0, 0, 1), -90.0),
              sdf.rotate(wedge_y(20.0, 12.0, 0.0, -36.0, -11.0, -1.0),
                         (0, 0, 1), 90.0)),
    sdf.translate(sdf.box((120.0, 30.0, 20.0)), (0.0, 0.0, 1.0)))

CRADLE_VOLUME = 880.0


class TestTheCradle:

    def test_a_2_mm_remnant_used_to_be_lost_and_is_not(self):
        with pytest.warns(RuntimeWarning, match="found no surface"):
            out = sdf.to_solid(CRADLE_GUSSETS, resolution=48)
        assert issolidclosed(out)
        assert volumeof(out) == pytest.approx(CRADLE_VOLUME, rel=1e-6)
        assert resolution_of(out) == 96

    def test_the_warning_hands_back_a_cell_that_works_first_time(self):
        with pytest.warns(RuntimeWarning) as caught:
            sdf.to_solid(CRADLE_GUSSETS, resolution=48)
        message = str(caught[0].message)
        assert "cell 2 mm (resolution 48)" in message
        assert "refined to cell 1 mm (resolution 96)" in message
        hint = float(message.split("cell_mm=")[1].split()[0])
        with warnings.catch_warnings():
            warnings.simplefilter("error")      # one pass, or this raises
            out = sdf.to_solid(CRADLE_GUSSETS, cell_mm=hint)
        assert issolidclosed(out)
        assert volumeof(out) == pytest.approx(CRADLE_VOLUME, rel=1e-6)

    def test_the_recorded_resolution_is_the_one_it_meshed_at(self):
        """So a later boolean on this solid starts from the cell that
        worked, instead of rediscovering it."""
        with pytest.warns(RuntimeWarning):
            out = sdf.to_solid(CRADLE_GUSSETS, resolution=48)
        lo, hi = CRADLE_GUSSETS.bounds
        longest = max(hi[i] - lo[i] for i in range(3))
        assert longest / resolution_of(out) == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# Empty, or too thin? -- the test that tells them apart
# ---------------------------------------------------------------------------


class TestProvesEmpty:

    GONE = sdf.subtract(sdf.box(10.0), sdf.box(20.0))

    def test_a_result_with_nothing_left_is_proved_at_once(self):
        """Including at the coarsest lattice, so an empty boolean costs one
        pass rather than an escalation up to the limit."""
        assert proves_empty(self.GONE, 8)
        assert proves_empty(self.GONE, 48)

    def test_a_part_that_is_there_is_never_proved_empty(self):
        for node in (sdf.box(10.0), sdf.sphere(5.0), CRADLE_GUSSETS):
            assert not proves_empty(node, 48)

    def test_a_thin_part_is_not_proved_empty_either(self):
        """The whole point: an empty lattice is not evidence of an empty
        part, and this is the case where they differ."""
        thin = sdf.subtract(sdf.box(10.0),
                            sdf.translate(sdf.box(20.0), (0.0, 0.0, 5.5)))
        assert len(sdf.dual_contour(thin, resolution=16).triangles) == 0
        assert not proves_empty(thin, 16)

    def test_an_unusable_lipschitz_constant_proves_nothing(self,
                                                           monkeypatch):
        assert proves_empty(self.GONE, 16)
        monkeypatch.setattr(type(self.GONE), "lipschitz",
                            property(lambda self: math.inf))
        assert not proves_empty(self.GONE, 16)

    def test_an_empty_result_meshes_to_the_empty_solid(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error")      # nothing to warn about
            out = sdf.to_solid(self.GONE, resolution=16)
        assert out[1] == []
        assert sdf.is_sdf_solid(out)            # it still carries its tree
        assert sdf.from_construction(get_construction(out)) == self.GONE


# ---------------------------------------------------------------------------
# The refinement itself
# ---------------------------------------------------------------------------


THIN = sdf.subtract(sdf.box(10.0),
                    sdf.translate(sdf.box(20.0), (0.0, 0.0, 5.5)))
THIN_VOLUME = 50.0

HAIRLINE = sdf.subtract(sdf.box(10.0),
                        sdf.translate(sdf.box(20.0), (0.0, 0.0, 4.99)))


class TestRefinement:

    def test_it_halves_the_cell_until_the_surface_shows(self):
        with pytest.warns(RuntimeWarning, match="found no surface"):
            out = sdf.to_solid(THIN, resolution=16)
        assert resolution_of(out) == 32
        assert volumeof(out) == pytest.approx(THIN_VOLUME, rel=1e-6)

    def test_nothing_extra_happens_on_the_happy_path(self, monkeypatch):
        """Refinement costs a mesh that was already going to work nothing:
        one contouring pass, and no emptiness proof."""
        import yapcad.sdf.convert as convert
        calls = {"contour": 0, "proof": 0}
        real = convert.dual_contour

        def counted(*args, **kwargs):
            calls["contour"] += 1
            return real(*args, **kwargs)

        monkeypatch.setattr(convert, "dual_contour", counted)
        monkeypatch.setattr(convert, "proves_empty",
                            lambda *a, **k: calls.__setitem__(
                                "proof", calls["proof"] + 1) or False)
        sdf.to_solid(sdf.box(10.0), resolution=16)
        assert calls == {"contour": 1, "proof": 0}

    def test_refine_false_refuses_instead(self):
        with pytest.raises(sdf.SdfError) as excinfo:
            sdf.to_solid(THIN, resolution=16, refine=False)
        message = str(excinfo.value)
        assert "found no surface" in message
        assert "refine=False" in message
        assert "cell_mm=" in message

    def test_the_limit_is_respected_and_reported(self):
        with pytest.raises(sdf.SdfError) as excinfo:
            sdf.to_solid(HAIRLINE, resolution=16, refine_limit=32)
        message = str(excinfo.value)
        assert "cell 0.625 mm down to 0.312 mm" in message
        assert "resolution 16 to 32" in message
        assert "cell_mm=" in message
        assert "refine_limit above 32" in message

    def test_a_limit_below_the_request_does_not_coarsen_anything(self):
        """refine_limit under the resolution asked for means no refinement,
        not a coarser mesh than was asked for."""
        with pytest.raises(sdf.SdfError, match="resolution 16"):
            sdf.to_solid(THIN, resolution=16, refine_limit=4)

    def test_an_unclosed_mesh_is_refined_too_and_says_so(self):
        """The other way a cell can be too coarse: the surface crosses one
        cell twice.  Same cure, different sentence."""
        node = sdf.intersect(sdf.box(20.0), sdf.gyroid(8.0, 0.5))
        with pytest.warns(RuntimeWarning,
                          match="could not close the mesh") as caught:
            out = sdf.to_solid(node, resolution=32)
        assert issolidclosed(out)
        assert resolution_of(out) == 64
        assert "cell_mm=" in str(caught[0].message)

    def test_a_compound_refines_as_a_whole(self):
        """Its bodies share one cell by design, so the cell that finds the
        thin body is the cell they all get, and the record stays a single
        resolution that reproduces every one of them."""
        pair = sdf.compound(THIN, sdf.translate(sdf.box(6.0), (20.0, 0, 0)))
        with pytest.warns(RuntimeWarning, match="found no surface"):
            out = sdf.to_solid(pair, resolution=16)
        assert len(out[1]) == 2
        # Refinement stops at *found*, not at well resolved: here the slab
        # ends up 1.2 cells thick, which is enough to mesh and not enough to
        # measure to a thousandth.  The warning names the cell for anyone
        # who needs better than that.
        assert volumeof(out) == pytest.approx(THIN_VOLUME + 216.0, rel=0.02)

    def test_a_compound_body_that_is_empty_is_dropped_not_refined_for(self):
        pair = sdf.compound(sdf.box(6.0),
                            sdf.translate(
                                sdf.subtract(sdf.box(10.0), sdf.box(20.0)),
                                (20.0, 0.0, 0.0)))
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            out = sdf.to_solid(pair, resolution=16)
        assert len(out[1]) == 1
        assert volumeof(out) == pytest.approx(216.0, rel=1e-3)


# ---------------------------------------------------------------------------
# cell_mm, the hint
# ---------------------------------------------------------------------------


class TestCellHint:

    def test_it_means_the_same_as_the_matching_resolution(self):
        by_cell = sdf.to_solid(sdf.box(10.0), cell_mm=0.5)
        by_resolution = sdf.to_solid(sdf.box(10.0), resolution=20)
        assert resolution_of(by_cell) == resolution_of(by_resolution)
        assert by_cell[1][0][1] == by_resolution[1][0][1]

    def test_it_rounds_down_so_the_cell_is_never_coarser_than_asked(self):
        out = sdf.to_solid(sdf.box(10.0), cell_mm=0.7)
        assert resolution_of(out) == 15           # ceil(10 / 0.7)
        assert 10.0 / resolution_of(out) <= 0.7

    def test_it_overrides_resolution(self):
        out = sdf.to_solid(sdf.box(10.0), resolution=8, cell_mm=0.25)
        assert resolution_of(out) == 40

    @pytest.mark.parametrize("cell", [0.0, -1.0, math.inf, math.nan])
    def test_a_cell_that_is_not_a_size_is_refused(self, cell):
        with pytest.raises(sdf.SdfError, match="positive and finite"):
            sdf.to_solid(sdf.box(10.0), cell_mm=cell)

    def test_an_unbounded_field_needs_bounds_to_measure_a_cell_against(self):
        with pytest.raises(sdf.SdfError, match="finite extent"):
            sdf.to_solid(sdf.gyroid(5.0, 0.3), cell_mm=0.5)
        region = ((-10.0, -10.0, -10.0), (10.0, 10.0, 10.0))
        # check=False because a lattice clipped by its bounds is open at
        # them by construction, which is not what is under test here.
        out = sdf.to_solid(sdf.gyroid(5.0, 0.3), cell_mm=0.5, bounds=region,
                           check=False)
        assert resolution_of(out) == 40

    def test_the_kernel_takes_it_too(self):
        mesh = sdf.dual_contour(sdf.box(10.0), cell_mm=0.5)
        assert len(mesh.triangles) > 0
        assert mesh.vertices.shape == sdf.dual_contour(
            sdf.box(10.0), resolution=20).vertices.shape

    def test_the_kernel_refuses_a_cell_that_is_not_a_size(self):
        with pytest.raises(sdf.SdfError, match="positive and finite"):
            sdf.dual_contour(sdf.box(10.0), cell_mm=0.0)
