"""OCC-free extrusion of 2D regions through manifold3d."""
import math
import textwrap

import pytest

pytest.importorskip("manifold3d")

from yapcad.geom import line, point  # noqa: E402
from yapcad.geom3d import solidbbox, volumeof  # noqa: E402
from yapcad.geom3d_util import extrude_region2d_mesh  # noqa: E402


def square(cx, cy, w, h):
    a = point(cx - w / 2, cy - h / 2)
    b = point(cx + w / 2, cy - h / 2)
    c = point(cx + w / 2, cy + h / 2)
    d = point(cx - w / 2, cy + h / 2)
    return [line(a, b), line(b, c), line(c, d), line(d, a)]


def test_prism_volume_and_extent():
    solid = extrude_region2d_mesh(square(0, 0, 4, 2), 3.0)
    assert volumeof(solid) == pytest.approx(24.0, rel=1e-6)
    lo, hi = solidbbox(solid)
    assert (lo[2], hi[2]) == pytest.approx((0.0, 3.0))


def test_oblique_direction_shears():
    solid = extrude_region2d_mesh(square(0, 0, 2, 2), 2.0, direction=[1, 0, 1])
    lo, hi = solidbbox(solid)
    assert hi[0] == pytest.approx(1.0 + math.sqrt(2.0))
    assert volumeof(solid) == pytest.approx(4.0 * math.sqrt(2.0), rel=1e-6)


def test_twist_is_counter_clockwise_and_volume_preserving():
    solid = extrude_region2d_mesh(square(6, 0, 4, 1), 20.0, twist_angle_deg=90.0)
    assert volumeof(solid) == pytest.approx(80.0, rel=0.015)
    top = [v for s in solid[1] for v in s[1] if v[2] > 19.999]
    cx = sum(v[0] for v in top) / len(top)
    cy = sum(v[1] for v in top) / len(top)
    assert (cx, cy) == pytest.approx((0.0, 6.0), abs=0.05)


def test_dsl_regions_with_holes_and_arcs():
    from yapcad.brep import occ_available
    from yapcad.dsl import compile_and_run

    if occ_available():
        pytest.skip("checks the OCC-free path")
    source = textwrap.dedent('''
        module ex
        command PLATE() -> solid:
            emit extrude(difference2d(rectangle(20.0, 10.0), disk(point(0.0, 0.0), 3.0)), 2.0)
        command D() -> solid:
            emit extrude(difference2d(disk(point(0.0, 0.0), 10.0), rectangle(30.0, 30.0, point2d(-15.0, 0.0))), 1.0)
        command HELIX() -> solid:
            emit helical_extrude(rectangle(4.0, 1.0, point2d(6.0, 0.0)), 20.0, 360.0)
    ''')
    expect = {"PLATE": (200 - math.pi * 9) * 2, "D": math.pi * 50, "HELIX": 80.0}
    for command, volume in expect.items():
        result = compile_and_run(source, command, {})
        assert result.success, result.error_message
        assert volumeof(result.geometry) == pytest.approx(volume, rel=0.015)
