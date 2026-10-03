"""Helical threads and hex nuts as fields.

The swept-mesh nut from ``metric_hex_nut`` is not a closed solid, so it
cannot be measured or combined reliably.  The field nut is the same part --
same thread profile, from ``threadgen``'s own definition -- as a closed,
SDF-authoritative solid.
"""

import json
import math

import numpy as np
import pytest

from yapcad import sdf
from yapcad.fasteners import metric_hex_nut
from yapcad.fasteners.builders import _make_profile_from_catalog
from yapcad.fasteners.catalog import get_nut_data
from yapcad.geom3d import issolidclosed, volumeof
from yapcad.metadata import get_solid_metadata
from yapcad.threadgen import ThreadProfile

M8 = get_nut_data("metric_coarse", "M8", "6H")
M8_THREAD = _make_profile_from_catalog(M8["thread"], internal=True)
WF, T = M8["body"]["width_across_flats"], M8["body"]["thickness"]


def in_nut(q, spec, wf, hand=1.0):
    """Independent point classification: inside the hexagon, outside the
    helical hole r < R(u)."""
    prof = np.array(sdf.thread_profile(spec))
    pitch = spec.P_pitch
    lead = pitch * spec.starts
    r = np.hypot(q[:, 0], q[:, 1])
    theta = np.arctan2(q[:, 1], q[:, 0])
    u = np.mod(q[:, 2] - hand * lead * theta / (2 * math.pi), pitch)
    flats = np.max([q[:, 0] * math.cos(k * math.pi / 3)
                    + q[:, 1] * math.sin(k * math.pi / 3)
                    for k in range(6)], axis=0)
    return (flats <= wf / 2) & (r > np.interp(u, prof[:, 0], prof[:, 1]))


def analytic_volume(spec, wf, t):
    """Hexagon prism less the hole, whose section has area pi mean(R^2)."""
    prof = np.array(sdf.thread_profile(spec))
    u = np.linspace(0, spec.P_pitch, 200001)
    hole = math.pi * t * np.mean(np.interp(u, prof[:, 0], prof[:, 1]) ** 2)
    return 2 * math.sqrt(3) * (wf / 2) ** 2 * t - hole


def sample(n, seed=0):
    rng = np.random.default_rng(seed)
    return rng.uniform([-WF * 0.6, -WF * 0.6, 0.0], [WF * 0.6, WF * 0.6, T],
                       (n, 3))


def test_thread_profile_matches_threadgen():
    from yapcad.threadgen import _radius_at
    prof = np.array(sdf.thread_profile(M8_THREAD))
    for u in np.linspace(0, M8_THREAD.P_pitch, 97):
        assert np.interp(u, prof[:, 0], prof[:, 1]) == pytest.approx(
            _radius_at(M8_THREAD, u, u), abs=1e-12)


@pytest.mark.parametrize("hand,starts", [("right", 1), ("left", 1),
                                         ("right", 2)])
def test_nut_field_classifies_every_point_correctly(hand, starts):
    spec = ThreadProfile(D_nominal=M8_THREAD.D_nominal,
                         P_pitch=M8_THREAD.P_pitch,
                         crest_flat_ratio=M8_THREAD.crest_flat_ratio,
                         root_flat_ratio=M8_THREAD.root_flat_ratio,
                         thread_depth_ratio=M8_THREAD.thread_depth_ratio,
                         handedness=hand, starts=starts, internal=True)
    node = sdf.hex_nut(spec, WF, T)
    q = sample(100000)
    expected = in_nut(q, spec, WF, 1.0 if hand == "right" else -1.0)
    assert np.array_equal(sdf.evaluate(node, q) < 0, expected)


def test_the_axis_is_inside_the_hole():
    """The guard-radius blend once had its sign flipped, putting a thin
    spurious cylinder of material along the axis."""
    node = sdf.hex_nut(M8_THREAD, WF, T)
    axis = np.column_stack([np.zeros(50), np.zeros(50),
                            np.linspace(0.1, T - 0.1, 50)])
    assert np.all(sdf.evaluate(node, axis) > 3.0)       # well outside


def test_the_lipschitz_bound_holds():
    node = sdf.hex_nut(M8_THREAD, WF, T)
    a = sample(20000, seed=1)
    b = a + np.random.default_rng(2).normal(scale=0.05, size=a.shape)
    ratio = np.abs(sdf.evaluate(node, a) - sdf.evaluate(node, b)) / \
        np.linalg.norm(a - b, axis=1)
    assert ratio.max() <= node.lipschitz


def test_sdf_nut_is_closed_and_has_the_analytic_volume():
    nut = metric_hex_nut("M8", representation="sdf")
    assert issolidclosed(nut)
    assert sdf.is_sdf_solid(nut)
    assert volumeof(nut) == pytest.approx(analytic_volume(M8_THREAD, WF, T),
                                          rel=1e-3)
    meta = get_solid_metadata(nut)
    assert "hex_nut" in meta.get("tags", []) and meta["hex_nut"]["pitch"]


def test_the_nut_tree_is_small_and_round_trips():
    node = sdf.hex_nut(M8_THREAD, WF, T)
    doc = sdf.tree_to_json(node)
    assert len(json.dumps(doc)) < 1000
    assert sdf.tree_from_json(doc) == node


def test_mesh_remains_the_default_and_bad_choices_are_refused():
    assert not sdf.is_sdf_solid(metric_hex_nut("M8"))
    with pytest.raises(ValueError, match="representation"):
        metric_hex_nut("M8", representation="voxels")


def test_thread_validation():
    good = sdf.thread_profile(M8_THREAD)
    with pytest.raises(sdf.SdfError, match="pitch"):
        sdf.thread(good, 0.0)
    with pytest.raises(sdf.SdfError, match="u = pitch"):
        sdf.thread(good, 2.0)
    with pytest.raises(sdf.SdfError, match="hand"):
        sdf.thread(good, M8_THREAD.P_pitch, hand="up")
    with pytest.raises(sdf.SdfError, match="wider"):
        sdf.hex_nut(M8_THREAD, 7.0, T)
