import math
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from yapcad.threadgen import metric_profile, sample_thread_profile, unified_profile

def test_thread_profile_sampling_consistent_counts():
    profile = metric_profile(10.0, 1.5)
    pts_theta0, wrap0 = sample_thread_profile(profile, 0.0, 4.5, 0.0)
    pts_theta90, wrap90 = sample_thread_profile(profile, 0.0, 4.5, 90.0)
    assert len(pts_theta0) == len(pts_theta90) >= 2
    assert wrap0 == wrap90 > 0
    radii = [pt[1] for pt in pts_theta0]
    assert max(radii) > min(radii) >= 0

def test_thread_profile_handles_internal_unified():
    profile = unified_profile(0.5, 13, internal=True)
    pts, wrap = sample_thread_profile(profile, 0.0, 50.0, 30.0)
    radii = [pt[1] for pt in pts]
    assert all(r > 0 for r in radii)
    assert 0.0 <= pts[0][0] <= 50.0
    assert 0.0 <= pts[-1][0] <= 50.0
    assert wrap > 0


def test_thread_profile_theta_offsets_translate_points_and_clamp():
    profile = metric_profile(12.0, 2.0)
    span = 30.0
    pts0, _ = sample_thread_profile(profile, 0.0, span, 0.0, samples_per_pitch=4)
    pts90, _ = sample_thread_profile(profile, 0.0, span, 90.0, samples_per_pitch=4)
    assert len(pts0) == len(pts90) > 2

    offset = profile.lead() * 90.0 / 360.0
    mid_idx = len(pts0) // 2
    base_mid = pts0[mid_idx][0]
    expected_mid = base_mid + offset
    assert math.isclose(pts90[mid_idx][0], expected_mid, rel_tol=0.0, abs_tol=1e-6)

    near_end_idx = len(pts0) - 2
    base_end = pts0[near_end_idx][0]
    expected_end = min(base_end + offset, span)
    assert math.isclose(pts90[near_end_idx][0], expected_end, rel_tol=0.0, abs_tol=1e-6)
