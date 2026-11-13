"""transformsolid, and text_on_surface's construction record (issue #108).

``text_on_surface`` rotated its text by rewriting vertex lists and copied
the old construction record onto the result, so the record described the
text before the rotation.  Text carries a ``procedure`` record today, which
nothing replays, so that could not be observed through text alone -- the
last test here substitutes a field-authored solid for the glyphs to make
the record observable on the real code path.
"""

import math

import pytest

from yapcad import sdf, text3d, xform
from yapcad.construction import get_construction
from yapcad.geom import point
from yapcad.geom3d import (
    rotatesolid,
    solidbbox,
    transformsolid,
    translatesolid,
)
from yapcad.geom3d_util import prism
from yapcad.sdf.node import from_construction


def placement(axis=(1, 1, 1), angle=40.0, delta=(3, -2, 5)):
    """An off-axis rotation followed by a translation."""
    return xform.Translation(point(*delta)).mul(
        xform.Rotation(point(*axis), angle))


def extents(box):
    return [box[1][k] - box[0][k] for k in range(3)]


def tree_bounds(sld):
    lo, hi = from_construction(get_construction(sld)).bounds
    return [list(lo), list(hi)]


class TestTransformsolid:
    def test_matches_rotate_then_translate(self):
        part = prism(2, 3, 4)
        moved = transformsolid(part, placement())
        expect = translatesolid(
            rotatesolid(part, 40.0, axis=point(1, 1, 1)), point(3, -2, 5))
        for s, e in zip(moved[1], expect[1]):
            for v, w in zip(s[1], e[1]):
                assert list(v) == pytest.approx(list(w), abs=1e-12)
            for v, w in zip(s[2], e[2]):
                assert list(v) == pytest.approx(list(w), abs=1e-12)

    def test_keeps_material_and_metadata(self):
        part = prism(1, 1, 1)
        part[2] = ["steel"]
        moved = transformsolid(part, placement())
        assert moved[2] == ["steel"]
        assert moved[3] == part[3]   # a procedure record is not a field

    @pytest.mark.parametrize("mat, match", [
        (xform.Scale(2.0), "rigid"),
        (xform.Scale(-1.0, 1.0, 1.0), "mirrors"),
    ])
    def test_refuses_non_rigid(self, mat, match):
        with pytest.raises(ValueError, match=match):
            transformsolid(prism(1, 1, 1), mat)

    def test_refuses_a_plain_list(self):
        with pytest.raises(ValueError, match="Matrix"):
            transformsolid(prism(1, 1, 1), [[1, 0, 0, 0]] * 4)

    def test_field_follows_the_mesh(self):
        node = sdf.box((6.0, 4.0, 2.0))
        mat = placement()
        meshed = transformsolid(sdf.to_solid(node, resolution=24), mat)
        with sdf.deferring():
            deferred = sdf.to_solid(node, resolution=24)
        moved = transformsolid(deferred, mat)
        # The tree alone, before any meshing, already knows where it went.
        assert solidbbox(moved)[0][:3] == pytest.approx(
            solidbbox(meshed)[0][:3], abs=0.5)
        sdf.realize(moved)
        a, b = solidbbox(moved), solidbbox(meshed)
        for k in range(3):
            assert a[0][k] == pytest.approx(b[0][k], abs=0.5)
            assert a[1][k] == pytest.approx(b[1][k], abs=0.5)


def test_text_on_surface_record_follows_the_rotation(monkeypatch):
    # A 20 x 4 x 1 field-authored "glyph block" in the text's local frame:
    # width along x, height along y, depth along z.
    block = sdf.to_solid(sdf.box((20.0, 4.0, 1.0)), resolution=24)
    monkeypatch.setattr(text3d, "text_solid_fitted",
                        lambda *a, **k: (block, 4.0))

    # Facing +x with z up: width runs along y, height along z, depth along x.
    placed = text3d.text_on_surface("ignored", (10, 0, 0), (1, 0, 0),
                                    (0, 0, 1), 20.0, 1.0)

    mesh_box = solidbbox(placed)
    assert extents(mesh_box) == pytest.approx([1, 20, 4], abs=0.5)
    tree = tree_bounds(placed)
    # Before the fix the record was the unrotated block -- 20 x 4 x 1.
    assert extents(tree) == pytest.approx(extents(mesh_box), abs=0.5)
    for k in range(3):
        assert tree[0][k] == pytest.approx(mesh_box[0][k], abs=0.5)
