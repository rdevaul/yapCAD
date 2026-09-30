"""Transforming an SDF-authored solid must carry its tree along.

Before this, ``geom3d.translatesolid`` moved an SDF solid's mesh and left the
tree that *defines* the solid where it was: the document then described a
sphere at the origin carrying a mesh a hundred units away, and anything
regenerating from the tree silently undid the move.
"""

import json

import numpy as np
import pytest

from yapcad import sdf
from yapcad.brep import occ_available
from yapcad.construction import get_construction
from yapcad.geom import point
from yapcad.geom3d import (
    mirrorsolid,
    rotatesolid,
    scalesolid,
    solidbbox,
    translatesolid,
)
from yapcad.geom3d_util import prism
from yapcad.io.geometry_json import geometry_from_json, geometry_to_json

needs_occ = pytest.mark.skipif(not occ_available(),
                               reason="pythonocc-core is not available")


def tree_of(solid):
    return sdf.from_construction(get_construction(solid))


def mesh_points(solid):
    return np.array([v[:3] for v in solid[1][0][1]])


def assert_tree_matches_mesh(solid, tolerance):
    """Every preview vertex lies on the (updated) tree's zero set."""
    residual = np.abs(sdf.evaluate(tree_of(solid), mesh_points(solid)))
    assert residual.max() < tolerance


BOX = sdf.box(10.0)


@pytest.fixture
def base():
    return sdf.to_solid(BOX, resolution=12)


class TestTreeFollowsTheMesh:

    def test_translate(self, base):
        moved = translatesolid(base, point(100.0, 0.0, 0.0))
        assert sdf.evaluate(tree_of(moved), [100.0, 0.0, 0.0]) == \
            pytest.approx(-5.0)
        assert_tree_matches_mesh(moved, 1e-6)

    def test_rotate_about_an_off_origin_centre(self, base):
        turned = rotatesolid(translatesolid(base, point(10.0, 0.0, 0.0)),
                             90.0, cent=point(50.0, 0.0, 0.0))
        assert_tree_matches_mesh(turned, 1e-6)
        lo, hi = solidbbox(turned)[:2]
        assert tree_of(turned).bounds[0][0] == pytest.approx(lo[0])
        assert tree_of(turned).bounds[1][1] == pytest.approx(hi[1])

    def test_rotate_with_an_explicit_matrix(self, base):
        from yapcad.xform import Rotation
        turned = rotatesolid(base, 30.0, mat=Rotation(point(0, 0, 1), 30.0))
        assert_tree_matches_mesh(turned, 1e-6)

    def test_uniform_scale(self, base):
        grown = scalesolid(base, 3.0)
        assert tree_of(grown).exact
        assert_tree_matches_mesh(grown, 1e-6)

    def test_nonuniform_scale(self, base):
        stretched = scalesolid(base, 3.0, 1.0, 1.0)
        assert not tree_of(stretched).exact
        assert_tree_matches_mesh(stretched, 1e-6)

    def test_scale_about_a_centre(self, base):
        grown = scalesolid(base, 2.0, cent=point(20.0, 0.0, 0.0))
        assert_tree_matches_mesh(grown, 1e-6)

    @pytest.mark.parametrize("plane", ["yz", "xz", "xy"])
    def test_mirror(self, base, plane):
        moved = translatesolid(base, point(30.0, 20.0, 10.0))
        assert_tree_matches_mesh(mirrorsolid(moved, plane), 1e-6)

    def test_successive_transforms_compose(self, base):
        solid = translatesolid(base, point(5.0, 0.0, 0.0))
        solid = rotatesolid(solid, 45.0, axis=point(1.0, 0.0, 0.0))
        solid = scalesolid(solid, 1.5)
        assert_tree_matches_mesh(solid, 1e-6)

    def test_the_original_is_left_untouched(self, base):
        translatesolid(base, point(100.0, 0.0, 0.0))
        assert tree_of(base) == BOX

    def test_meshing_parameters_are_dropped(self, base):
        """A mesh moved after meshing is not what those parameters would
        regenerate, so the record stops claiming it."""
        assert sdf.meshing_from_construction(get_construction(base))
        moved = translatesolid(base, point(1.0, 0.0, 0.0))
        assert sdf.meshing_from_construction(get_construction(moved)) is None

    def test_solids_without_a_tree_are_unaffected(self):
        cube = prism(2.0, 2.0, 2.0)
        before = get_construction(cube)
        after = get_construction(translatesolid(cube, point(5.0, 0.0, 0.0)))
        assert after == before

    def test_the_document_describes_the_moved_solid(self, base):
        moved = translatesolid(base, point(100.0, 0.0, 0.0))
        document = json.loads(json.dumps(geometry_to_json([moved],
                                                          units="mm")))
        entry = next(e for e in document["entities"] if e["type"] == "solid")
        assert entry["representations"]["sdf"]["bounds"][0] == \
            pytest.approx(entry["boundingBox"][0])
        restored = geometry_from_json(document)
        assert tree_of(restored[0]) == tree_of(moved)


@needs_occ
class TestDerivedBrepFollows:

    @pytest.fixture
    def with_brep(self):
        return sdf.to_solid(BOX, resolution=12, brep=True)

    def _brep_record(self, solid):
        from yapcad.metadata import get_solid_metadata
        return get_solid_metadata(solid).get("brep")

    def test_a_similarity_retags_the_brep_to_the_new_tree(self, with_brep):
        moved = translatesolid(with_brep, point(100.0, 0.0, 0.0))
        assert self._brep_record(moved)["derivedFrom"] == \
            {"sdf": tree_of(moved).digest}

    def test_the_retagged_brep_is_the_replay_of_the_new_tree(self,
                                                             with_brep):
        """The re-tag is only honest if it is literally true."""
        from yapcad.brep import brep_from_solid
        from yapcad.sdf import occ
        turned = rotatesolid(with_brep, 30.0, cent=point(5.0, 5.0, 0.0))
        carried = occ.bounding_box(brep_from_solid(turned).shape)
        replayed = occ.bounding_box(occ.to_occ_shape(tree_of(turned)))
        np.testing.assert_allclose(carried, replayed, atol=1e-6)

    def test_it_still_serialises_after_a_similarity(self, with_brep):
        solid = scalesolid(mirrorsolid(
            translatesolid(with_brep, point(3.0, 0.0, 0.0)), "xz"), 2.0)
        document = json.loads(json.dumps(geometry_to_json([solid],
                                                          units="mm")))
        entry = next(e for e in document["entities"] if e["type"] == "solid")
        assert entry["representations"]["brep"]["role"] == "derived"
        geometry_from_json(document)

    def test_a_nonuniform_scale_drops_the_stale_brep(self, with_brep):
        from yapcad.brep import has_brep_data
        stretched = scalesolid(with_brep, 3.0, 1.0, 1.0)
        assert not has_brep_data(stretched)
        json.dumps(geometry_to_json([stretched], units="mm"))


def test_transform_honours_an_xform_matrix_transpose_flag():
    """sdf.transform used to read a Matrix's raw rows and ignore its
    transpose flag, silently applying M^T."""
    from yapcad.xform import Rotation
    forward = Rotation(point(0, 0, 1), 30.0)
    flagged = Rotation(point(0, 0, 1), 30.0)
    flagged.trans = True                  # now denotes the -30 degree turn
    reverse = Rotation(point(0, 0, 1), -30.0)
    base = sdf.box((10.0, 2.0, 2.0))
    pts = np.random.default_rng(7).uniform(-8, 8, size=(200, 3))
    assert not np.allclose(sdf.evaluate(sdf.transform(base, forward), pts),
                           sdf.evaluate(sdf.transform(base, flagged), pts))
    np.testing.assert_allclose(
        sdf.evaluate(sdf.transform(base, flagged), pts),
        sdf.evaluate(sdf.transform(base, reverse), pts), atol=1e-12)
