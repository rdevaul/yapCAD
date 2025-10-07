import math
from pathlib import Path

import pytest

from yapcad.geom import (
    catmullrom,
    dist,
    isnurbs,
    iscatmullrom,
    nurbs,
    point,
)
from yapcad.geom_util import geomlist2poly
from yapcad.spline import (
    evaluate_catmullrom,
    evaluate_nurbs,
    sample_catmullrom,
    sample_nurbs,
)
from yapcad.ezdxf_drawable import ezdxfDraw


def _close(a, b, tol=1e-6):
    assert dist(point(a), point(b)) <= tol


def test_catmullrom_samples_match_evaluator():
    ctrl = [
        point(0, 0),
        point(1, 1),
        point(2, 0),
        point(3, 1),
    ]
    curve = catmullrom(ctrl, closed=False)
    assert iscatmullrom(curve)

    segments = len(ctrl) - 1
    samples = sample_catmullrom(curve, segments_per_span=8)

    _close(samples[0], point(0, 0))
    _close(samples[-1], point(3, 1))

    idx = 0
    for segment in range(segments):
        for step in range(1, 8 + 1):
            idx += 1
            u = (segment + step / 8.0) / segments
            _close(samples[idx], evaluate_catmullrom(curve, u), tol=5e-4)


def test_catmullrom_closed_geomlist_support():
    ctrl = [
        point(0, 0),
        point(2, 0),
        point(2, 2),
        point(0, 2),
    ]
    curve = catmullrom(ctrl, closed=True)
    poly = geomlist2poly([curve], minlen=0.1)
    assert poly
    assert dist(poly[0], poly[-1]) <= 1e-6
    assert len(poly) > len(ctrl)


def _make_simple_nurbs():
    ctrl = [point(0, 0), point(1, 2), point(2, 0)]
    degree = 2
    return nurbs(ctrl, degree=degree)


def test_default_open_uniform_knot_vector():
    curve = _make_simple_nurbs()
    knots = curve[2]['knots']
    assert knots == [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]

    higher = nurbs([
        point(0, 0),
        point(1, 1),
        point(2, -1),
        point(3, 0),
        point(4, 2),
    ], degree=3)
    assert higher[2]['knots'] == [0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0]


def test_nurbs_samples_match_evaluator():
    curve = _make_simple_nurbs()
    assert isnurbs(curve)

    samples = 40
    sampled = sample_nurbs(curve, samples=samples)
    _close(sampled[0], evaluate_nurbs(curve, 0.0))
    _close(sampled[-1], evaluate_nurbs(curve, 1.0))

    for i in range(1, samples - 1):
        u = i / (samples - 1)
        _close(sampled[i], evaluate_nurbs(curve, u), tol=5e-4)


@pytest.mark.slow
def test_spline_dxf_output(tmp_path):
    cat = catmullrom([
        point(0, 0),
        point(4, 3),
        point(8, -1),
        point(12, 2),
    ])
    nurb = nurbs([
        point(0, 0),
        point(2, 4),
        point(4, 0),
        point(6, 2),
    ], degree=3)

    drawer = ezdxfDraw()
    drawer.filename = str(tmp_path / 'splines')
    drawer.pointstyle = 'xo'
    drawer.pointsize = 0.2
    drawer.linecolor='coral'

    for ctrl in cat[1]:
        drawer.draw_point(ctrl)

    drawer.linecolor='red'
    drawer.draw(cat)

    drawer.linecolor='aqua'
    for ctrl in nurb[1]:
        drawer.draw_point(ctrl)

    drawer.linecolor='blue'
    drawer.draw(nurb)

    drawer.display()
    output = tmp_path / 'splines.dxf'
    assert output.exists()
    assert output.stat().st_size > 0
