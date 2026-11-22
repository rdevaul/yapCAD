"""Engraved and embossed block text is watertight, cuts to the depth asked
for, and engraves field-authored parts as fields (#82)."""

import collections
import io

import pytest

from yapcad import sdf
from yapcad.geom import point, vect
from yapcad.geom3d import issolidclosed, volumeof
from yapcad.geom3d_util import prism
from yapcad.io.stl import write_stl
from yapcad.text3d import (
    BLOCK_FONT,
    CHAR_HEIGHT,
    _glyph_cells,
    engrave_text,
    text_field,
    text_solid,
)


def _stl_edge_defects(obj):
    """(boundary, non-manifold) edges of ``obj`` welded by position, as an
    STL reader sees it."""
    trimesh = pytest.importorskip("trimesh")
    buf = io.BytesIO()
    write_stl(obj, buf)
    buf.seek(0)
    mesh = trimesh.load(buf, file_type="stl", process=True)
    edges = collections.Counter(map(tuple, mesh.edges_sorted.tolist()))
    return (sum(c == 1 for c in edges.values()),
            sum(c > 2 for c in edges.values()))


@pytest.mark.parametrize("char", sorted(BLOCK_FONT))
def test_no_glyph_has_strokes_meeting_only_at_a_corner(char):
    cells = _glyph_cells(BLOCK_FONT[char])
    for (i, j) in cells:
        for di, dj in ((1, 1), (1, -1)):
            if (i + di, j + dj) in cells:
                assert (i + di, j) in cells or (i, j + dj) in cells, (char, i, j)


@pytest.mark.parametrize("char", sorted(c for c in BLOCK_FONT if c != " "))
def test_every_block_glyph_extrudes_to_a_closed_solid(char):
    glyph = text_solid(char, 3.0, 0.6, 1.0, font="block")
    assert issolidclosed(glyph)
    cells = len(_glyph_cells(BLOCK_FONT[char]))
    assert volumeof(glyph) == pytest.approx(cells * (3.0 / CHAR_HEIGHT) ** 2 * 0.6)


def _stroke_area(text):
    scale = 3.0 / CHAR_HEIGHT
    return sum(len(_glyph_cells(BLOCK_FONT[c])) for c in text) * scale ** 2


def test_engraving_a_mesh_is_watertight_and_cuts_full_depth():
    box = prism(40.0, 20.0, 4.0)
    cut = engrave_text(box, "ABCD", point(0, 0, 2), vect(0, 0, 1, 0),
                       3.0, 0.5, 1.0, font="block")
    assert _stl_edge_defects(cut) == (0, 0)
    assert volumeof(cut) == pytest.approx(3200.0 - _stroke_area("ABCD") * 0.5,
                                          rel=1e-6)


def test_engraving_a_field_stays_a_field():
    plate = sdf.to_solid(sdf.box((40.0, 20.0, 4.0)), cell_mm=1.0, defer=True)
    cut = engrave_text(plate, "ABCD", point(0, 0, 2), vect(0, 0, 1, 0),
                       3.0, 0.5, 1.0, font="block")
    assert sdf.is_sdf_solid(cut)
    meshed = sdf.realize(cut)
    assert issolidclosed(meshed)
    assert volumeof(meshed) == pytest.approx(
        3200.0 - _stroke_area("ABCD") * 0.5, rel=2e-3)


def test_the_text_field_matches_the_text_solid():
    field = text_field("AB8", 3.0, 0.6, 1.0, font="block")
    meshed = sdf.to_solid(field, cell_mm=0.1)
    solid_ = text_solid("AB8", 3.0, 0.6, 1.0, font="block")
    # Dual contouring rounds the strokes' corners at this cell; the field
    # itself is exact, which the samples below check directly.
    assert volumeof(meshed) == pytest.approx(volumeof(solid_), rel=1.5e-2)
    import numpy as np
    from yapcad.text3d import _glyph_cells as cells_of
    scale = 3.0 / CHAR_HEIGHT
    inside = [((i + 0.5) * scale, (j + 0.5) * scale, 0.3)
              for (i, j) in cells_of(BLOCK_FONT["A"])]
    assert (sdf.evaluate(field, np.array(inside)) < 0).all()
    hollow = [(2.5 * scale, 4.5 * scale, 0.3), (2.5 * scale, 1.5 * scale, 0.3)]
    assert (sdf.evaluate(field, np.array(hollow)) > 0).all()


def test_a_failed_engraving_raises_instead_of_returning_the_target(monkeypatch):
    import yapcad.text3d as text3d

    def broken(*_args, **_kwargs):
        raise ValueError("no engine")
    monkeypatch.setattr(text3d, "solid_boolean", broken)
    with pytest.raises(RuntimeError, match="engrave_text"):
        engrave_text(prism(10.0, 10.0, 4.0), "A", point(0, 0, 2),
                     vect(0, 0, 1, 0), font="block")
