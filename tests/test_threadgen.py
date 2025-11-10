import math
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from yapcad.threadgen import metric_profile, sample_thread_profile, unified_profile

def test_thread_profile_sampling_consistent_counts():
    profile = metric_profile(10.0, 1.5)
    pts_theta0 = sample_thread_profile(profile, 0.0, 4.5, 0.0)
    pts_theta90 = sample_thread_profile(profile, 0.0, 4.5, 90.0)
    assert len(pts_theta0) == len(pts_theta90) >= 2
    radii = [pt[1] for pt in pts_theta0]
    assert max(radii) > min(radii) >= 0

def test_thread_profile_handles_internal_unified():
    profile = unified_profile(0.5, 13, internal=True)
    pts = sample_thread_profile(profile, 0.0, 50.0, 30.0)
    radii = [pt[1] for pt in pts]
    assert all(r > 0 for r in radii)
    assert math.isclose(pts[0][0], 0.0)
    assert math.isclose(pts[-1][0], 50.0)
