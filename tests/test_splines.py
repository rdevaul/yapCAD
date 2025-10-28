import math
import subprocess
import sys
from pathlib import Path

import ezdxf

import pytest

from yapcad.geom import (
    catmullrom,
    dist,
    isnurbs,
    iscatmullrom,
    nurbs,
    point,
    sample,
)
from yapcad.geom3d import solidbbox
from yapcad.geom3d_util import extrude, poly2surfaceXY
from yapcad.geom_util import geomlist2poly
from yapcad.package import create_package_from_entities
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


def _sample_closed_curve(curve, segments):
    pts = [point(sample(curve, i / segments)) for i in range(segments)]
    pts.append(point(sample(curve, 0.0)))
    return pts


def test_spline_extrusion_package(tmp_path):
    outer_ctrl = [
        point(-40, 0),
        point(-30, 30),
        point(0, 40),
        point(30, 30),
        point(40, 0),
        point(30, -30),
        point(0, -40),
        point(-30, -30),
    ]
    outer_curve = catmullrom(outer_ctrl, closed=True)

    hole1_ctrl = [
        point(-10, 0),
        point(-6, 8),
        point(0, 10),
        point(6, 8),
        point(10, 0),
        point(6, -8),
        point(0, -10),
        point(-6, -8),
        point(-10, 0),
    ]
    hole1_curve = catmullrom(hole1_ctrl, closed=True)

    hole2_ctrl = [
        point(15, 5),
        point(18, 8),
        point(22, 5),
        point(18, 2),
        point(15, 5),
    ]
    hole2_curve = nurbs(hole2_ctrl, degree=3)

    outer_poly = _sample_closed_curve(outer_curve, segments=256)
    hole_polys = [
        _sample_closed_curve(hole1_curve, segments=160),
        _sample_closed_curve(hole2_curve, segments=120),
    ]
    surface, _ = poly2surfaceXY(outer_poly, hole_polys, minlen=0.1)
    solid = extrude(surface, 6.0)
    bbox = solidbbox(solid)
    assert pytest.approx(bbox[1][2] - bbox[0][2], abs=1e-6) == 6.0

    sketch_geom = [outer_curve, hole1_curve, hole2_curve]

    pkg_root = tmp_path / "spline_extrusion.ycpkg"
    manifest = create_package_from_entities(
        [
            {"geometry": solid, "metadata": {"layer": "solid", "tags": ["spline"]}},
            {"geometry": sketch_geom, "metadata": {"layer": "sketch"}},
        ],
        pkg_root,
        name="Spline Extrusion",
        version="0.1",
        units="mm",
        overwrite=True,
    )
    manifest.recompute_hashes()
    manifest.save()

    primary_json = pkg_root / "geometry" / "primary.json"
    import json
    with primary_json.open() as fh:
        doc = json.load(fh)
    sketch_entries = [entry for entry in doc["entities"] if entry["type"] == "sketch"]
    assert sketch_entries
    primitive_kinds = {prim["kind"] for entry in sketch_entries for prim in entry.get("primitives", [])}
    assert {"catmullrom", "nurbs"}.issubset(primitive_kinds)

    export_dir = tmp_path / "exports"
    subprocess.check_call(
        [
            sys.executable,
            "tools/ycpkg_export.py",
            str(pkg_root),
            "--format",
            "dxf",
            "--format",
            "step",
            "--output",
            str(export_dir),
            "--overwrite",
        ],
        cwd=Path(__file__).resolve().parents[1],
    )
    dxf_path = export_dir / "sketches.dxf"
    step_path = export_dir / "solid_01.step"
    assert dxf_path.exists() and dxf_path.stat().st_size > 0
    assert step_path.exists() and step_path.stat().st_size > 0
    assert "LWPOLYLINE" in dxf_path.read_text()

    doc = ezdxf.readfile(dxf_path)
    polylines = list(doc.modelspace().query("LWPOLYLINE"))
    assert polylines
    # Pick the smallest loop (nurbs hole)
    def loop_span(pl):
        xs = [v[0] for v in pl.get_points()]
        ys = [v[1] for v in pl.get_points()]
        return max(xs) - min(xs) + max(ys) - min(ys)

    small_loop = min(polylines, key=loop_span)
    assert small_loop.closed
