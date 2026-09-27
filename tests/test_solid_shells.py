"""Tests for derived outer/void shell partitioning.

yapCAD mesh solids carry cavities as inward-wound shells inside the surface
list; there is no voids slot. ``geom3d.solid_shells`` recovers the partition
on demand, and geometry JSON uses it at the interchange boundary.
"""

from __future__ import annotations

import json

import pytest

from yapcad.geom3d import (
    issolidclosed,
    reversesurface,
    solid,
    solid_shells,
    volumeof,
)
from yapcad.geom3d_util import prism
from yapcad.io.geometry_json import geometry_from_json, geometry_to_json


def _cube_surfaces(size, center=None):
    if center is None:
        return prism(size, size, size)[1]
    return prism(size, size, size, center=center)[1]


def _hollow_cube():
    """4x4x4 cube with a 2x2x2 cavity: volume 64 - 8 = 56."""
    outer = _cube_surfaces(4)
    inner = [reversesurface(s) for s in _cube_surfaces(2)]
    return solid(outer + inner)


def _roundtrip(entities):
    document = json.loads(json.dumps(geometry_to_json(entities, units="mm")))
    return document, geometry_from_json(document)


def _solid_entry(document):
    return next(e for e in document["entities"] if e["type"] == "solid")


# ---------------------------------------------------------------------------
# The premise: cavities already work, encoded as normal orientation
# ---------------------------------------------------------------------------

def test_hollow_solid_is_closed_and_has_correct_volume():
    hollow = _hollow_cube()
    assert issolidclosed(hollow)
    assert volumeof(hollow) == pytest.approx(56.0)


def test_unreversed_inner_shell_adds_volume_instead():
    """Orientation is what carries the meaning, not a slot."""
    both_outward = solid(_cube_surfaces(4) + _cube_surfaces(2))
    assert volumeof(both_outward) == pytest.approx(72.0)


# ---------------------------------------------------------------------------
# solid_shells
# ---------------------------------------------------------------------------

def test_simple_solid_is_one_outer_shell():
    outer, voids = solid_shells(solid(_cube_surfaces(4)))
    assert len(outer) == 1
    assert len(outer[0]) == 6
    assert voids == []


def test_hollow_solid_splits_into_outer_and_void():
    outer, voids = solid_shells(_hollow_cube())
    assert len(outer) == 1 and len(outer[0]) == 6
    assert len(voids) == 1 and len(voids[0]) == 6


def test_disjoint_shells_are_both_outer():
    """A compound of two separate bodies is not a solid with a cavity."""
    compound = solid(_cube_surfaces(4) + _cube_surfaces(2, center=[20, 0, 0, 1]))
    outer, voids = solid_shells(compound)
    assert len(outer) == 2
    assert voids == []


def test_empty_solid():
    assert solid_shells(solid([])) == ([], [])


def test_open_shell_is_never_classified_as_a_void():
    """An open surface cannot bound a cavity, whatever its winding."""
    single_face = reversesurface(_cube_surfaces(2)[0])
    outer, voids = solid_shells(solid([single_face]))
    assert voids == []
    assert len(outer) == 1


def test_partition_is_deterministic():
    hollow = _hollow_cube()
    first = solid_shells(hollow)
    for _ in range(3):
        assert solid_shells(hollow) == first


def test_rejects_non_solids():
    with pytest.raises(ValueError):
        solid_shells(["not", "a", "solid"])


def test_surface_granularity_limit_is_documented_behaviour():
    """A one-surface boundary (the shape of an OCC tessellation) cannot be
    split, so it is conservatively reported as a single outer shell."""
    merged = _cube_surfaces(4)[0]
    outer, voids = solid_shells(solid([merged]))
    assert len(outer) == 1
    assert voids == []


# ---------------------------------------------------------------------------
# Interchange
# ---------------------------------------------------------------------------

def test_ordinary_solid_emits_no_voids():
    document, _ = _roundtrip([prism(2, 2, 2)])
    entry = _solid_entry(document)
    assert entry["voids"] == []
    assert len(entry["shell"]) == 6


def test_hollow_solid_emits_the_derived_partition():
    document, _ = _roundtrip([_hollow_cube()])
    entry = _solid_entry(document)
    assert len(entry["shell"]) == 6
    assert len(entry["voids"]) == 1
    assert len(entry["voids"][0]) == 6
    # Shell and void surfaces are distinct entities.
    assert not set(entry["shell"]) & set(entry["voids"][0])


def test_void_surfaces_rejoin_the_shell_list_on_read():
    """A mesh solid has no voids slot, so they come back as shell surfaces."""
    _, restored = _roundtrip([_hollow_cube()])
    sld = restored[0]
    assert len(sld[1]) == 12
    assert sld[2] == []          # material, not voids


def test_hollow_solid_survives_roundtrip():
    original = _hollow_cube()
    _, restored = _roundtrip([original])
    sld = restored[0]
    assert issolidclosed(sld)
    assert volumeof(sld) == pytest.approx(volumeof(original))
    outer, voids = solid_shells(sld)
    assert len(outer) == 1 and len(voids) == 1


def test_roundtrip_is_idempotent():
    """Two passes must produce the same partition, for stable package hashes."""
    once_doc, once = _roundtrip([_hollow_cube()])
    twice_doc, _ = _roundtrip(once)
    assert len(_solid_entry(twice_doc)["shell"]) == len(_solid_entry(once_doc)["shell"])
    assert len(_solid_entry(twice_doc)["voids"]) == len(_solid_entry(once_doc)["voids"])


def test_document_voids_are_honoured_on_read():
    """A document from any producer that emits voids must load correctly."""
    document, _ = _roundtrip([_hollow_cube()])
    restored = geometry_from_json(document)[0]
    assert volumeof(restored) == pytest.approx(56.0)
