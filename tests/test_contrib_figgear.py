from math import isclose

import pytest

from yapcad.contrib.figgear import make_gear_figure


def polygon_area(points):
    area = 0.0
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        area += x0 * y1 - x1 * y0
    return area / 2.0


@pytest.mark.parametrize("bottom_type", ["line", "spline"])
def test_make_gear_figure_properties(bottom_type):
    points, info = make_gear_figure(
        m=2.0,
        z=24,
        alpha_deg=20.0,
        bottom_type=bottom_type,
        involute_step=0.3,
        spline_division_num=16,
    )

    assert len(points) > 0
    assert points[0] == points[-1]

    expected_keys = {
        "diameter_addendum",
        "diameter_pitch",
        "diameter_base",
        "diameter_dedendum",
        "radius_addendum",
        "radius_pitch",
        "radius_base",
        "radius_dedendum",
    }
    assert expected_keys == set(info.keys())

    assert info["radius_addendum"] > info["radius_pitch"] > info["radius_base"]
    assert info["radius_dedendum"] < info["radius_base"]

    area = polygon_area(points)
    assert area > 0.0

    # The addendum radius should match module-based expectations.
    expected_addendum = (24 * 2.0 + 2 * 2.0) / 2.0
    assert isclose(info["radius_addendum"], expected_addendum, rel_tol=1e-9)
