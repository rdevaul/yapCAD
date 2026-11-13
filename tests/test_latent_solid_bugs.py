"""Three latent bugs in the solid machinery, found while building something else.

None of these were reported by a user, because each is the kind that returns a
plausible answer rather than an error: an empty file that looks like a file, a
centre that is a centre of *something*, an import that nothing happens to
reach.  Each test below fails on the code as it was.
"""

import math

import pytest

from yapcad import sdf
from yapcad.geom3d import center, issolid, solid, solidbbox
from yapcad.geometry import Geometry
from yapcad.io.stl import write_stl


def two_bodies():
    """A solid whose first surface is nowhere near its centre."""
    tree = sdf.compound(sdf.box((10.0, 10.0, 10.0)),
                        sdf.translate(sdf.box((10.0, 10.0, 10.0)),
                                      (100.0, 0.0, 0.0)))
    built = sdf.to_solid(tree, resolution=16)
    assert len(built[1]) == 2, "expected one surface per body"
    return built


class TestEmptySolidConstructor:
    """``solid()`` documents an empty solid and raised instead.

    ``args`` is a tuple, so the ``args == []`` test was always False and the
    no-argument call fell through to the error at the bottom of the function.
    """

    def test_it_returns_an_empty_solid(self):
        assert solid() == ['solid', [], [], []]

    def test_the_documented_spelling_still_works(self):
        assert solid([]) == ['solid', [], [], []]

    def test_the_result_is_a_solid(self):
        assert issolid(solid())
        assert issolid(solid(), fast=False)

    def test_a_genuinely_bad_argument_is_still_refused(self):
        with pytest.raises(ValueError, match="bad arguments"):
            solid("not geometry")


class TestSolidCentre:
    """``Geometry.center`` reported the centre of a solid's *first surface*.

    For a one-surface solid that is the same answer, which is why it went
    unnoticed; for anything else it is the centre of an arbitrary piece.
    """

    def test_the_centre_is_the_centre_of_the_whole_solid(self):
        built = two_bodies()
        got = Geometry(built).center
        assert got[0] == pytest.approx(50.0, abs=0.5)

    def test_it_agrees_with_geom3d(self):
        built = two_bodies()
        assert Geometry(built).center[:3] == \
            pytest.approx(center(built)[:3], abs=1e-9)

    def test_it_is_the_midpoint_of_the_bounding_box(self):
        built = two_bodies()
        box = solidbbox(built)
        want = [(box[0][i] + box[1][i]) / 2.0 for i in range(3)]
        assert Geometry(built).center[:3] == pytest.approx(want, abs=1e-9)

    def test_it_is_not_the_first_surfaces_centre(self):
        # The old answer, pinned so a regression is unambiguous.
        built = two_bodies()
        first = center(built[1][0])
        assert Geometry(built).center[0] != pytest.approx(first[0], abs=1.0)

    def test_a_single_surface_solid_is_unchanged(self):
        one = sdf.to_solid(sdf.box((10.0, 10.0, 10.0)), resolution=16)
        assert Geometry(one).center[:3] == pytest.approx([0.0, 0.0, 0.0],
                                                           abs=1e-6)


class TestEmptyStl:
    """An STL export that produced no triangles wrote the file anyway.

    A zero-facet STL is a valid file describing nothing, which is worse than
    a refusal.  ``write_step`` already refuses; this makes the pair agree.
    """

    def test_an_empty_solid_is_refused(self, tmp_path):
        with pytest.raises(ValueError, match="no triangles"):
            write_stl(solid([]), tmp_path / "empty.stl")

    def test_no_file_is_left_behind(self, tmp_path):
        target = tmp_path / "empty.stl"
        with pytest.raises(ValueError):
            write_stl(solid([]), target)
        assert not target.exists()

    def test_a_real_solid_still_writes(self, tmp_path):
        target = tmp_path / "box.stl"
        write_stl(sdf.to_solid(sdf.box((10.0, 10.0, 10.0)), resolution=16),
                  target)
        assert target.stat().st_size > 0

    def test_it_says_the_same_thing_step_says(self, tmp_path):
        # Both refusals should read alike; only the format differs.
        with pytest.raises(ValueError) as caught:
            write_stl(solid([]), tmp_path / "empty.stl")
        assert "STL export" in str(caught.value)
