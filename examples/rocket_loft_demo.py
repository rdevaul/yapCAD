"""Lofted rocket body demo with BREP support.

Builds a rocket airframe by lofting spline-based cross-sections, unions a
revolved nose cone, hollows the interior with a BREP difference, and either
visualizes the result or exports a .ycpkg.
"""

import argparse
import math
from pathlib import Path

import ezdxf

from yapcad.geom import point, mag, sub
from yapcad.geom import polybbox


def _circle_loop(radius: float, samples: int = 64):
    pts = []
    for i in range(samples):
        ang = 2 * math.pi * i / samples
        pts.append(point(radius * math.cos(ang), radius * math.sin(ang), 0.0))
    return pts
def _clean_loop(loop):
    if loop and mag(sub(loop[0], loop[-1])) < 1e-9:
        return loop[:-1]
    return loop


def _perimeter(loop):
    if len(loop) < 2:
        return 0.0
    total = 0.0
    for i in range(len(loop)):
        p0 = loop[i]
        p1 = loop[(i + 1) % len(loop)]
        total += mag(sub(p1, p0))
    return total


def _resample_loop(loop, target):
    """Resample a closed loop to have exactly target points."""
    loop = list(loop)
    n = len(loop)
    if n == target:
        return loop
    if n < 2:
        raise ValueError("loop too small to resample")
    cum = [0.0]
    for i in range(n):
        seg = mag(sub(loop[(i + 1) % n], loop[i]))
        cum.append(cum[-1] + seg)
    total = cum[-1]
    if total < 1e-12:
        # Fallback: rebuild as a circle at the same Z using the loop's radius
        zval = loop[0][2] if loop else 0.0
        radius = max(_loop_radius(loop), 0.5)
        circ = _clean_loop(_circle_loop(radius, samples=target))
        for p in circ:
            p[2] = zval
        return circ
    step = total / target
    res = []
    seg_idx = 0
    for k in range(target):
        dist = k * step
        while seg_idx < n and cum[seg_idx + 1] < dist - 1e-12:
            seg_idx += 1
        seg_idx = min(seg_idx, n - 1)
        p0 = loop[seg_idx % n]
        p1 = loop[(seg_idx + 1) % n]
        seg_len = cum[seg_idx + 1] - cum[seg_idx]
        t = 0.0 if seg_len < 1e-12 else (dist - cum[seg_idx]) / seg_len
        res.append(point(
            p0[0] + t * (p1[0] - p0[0]),
            p0[1] + t * (p1[1] - p0[1]),
            p0[2] + t * (p1[2] - p0[2]),
        ))
    return res


def _superellipse_loop(a: float, b: float, exp: float, samples: int, z: float = 0.0):
    """Generate a 4-lobed superellipse loop (p < 1 gives cross-like fins)."""
    pts = []
    for i in range(samples):
        ang = 2 * math.pi * i / samples
        ca = math.cos(ang)
        sa = math.sin(ang)
        x = a * math.copysign(abs(ca) ** exp, ca)
        y = b * math.copysign(abs(sa) ** exp, sa)
        pts.append(point(x, y, z))
    return _clean_loop(pts)


def _cross_loop(radius: float, fin_amp: float, fin_exp: float, samples: int, z: float = 0.0):
    """Polar loop with four fins along axes; fin_amp controls protrusion."""
    pts = []
    for i in range(samples):
        ang = 2 * math.pi * i / samples
        r = radius * (1.0 + fin_amp * (abs(math.cos(2 * ang)) ** fin_exp))
        x = r * math.cos(ang)
        y = r * math.sin(ang)
        pts.append(point(x, y, z))
    return _clean_loop(pts)


def _scale_loop(loop, scale: float):
    return [point(p[0] * scale, p[1] * scale, p[2]) for p in loop]


def _loop_radius(loop):
    if not loop:
        return 0.0
    bbox = polybbox([[p[0], p[1], 0.0] for p in loop])
    dx = bbox[1][0] - bbox[0][0]
    dy = bbox[1][1] - bbox[0][1]
    return 0.5 * max(dx, dy)


def _build_body_sections():
    samples = 96
    base = _cross_loop(1.25, fin_amp=0.30, fin_exp=3.0, samples=samples, z=0.0)
    mid = _cross_loop(1.6, fin_amp=0.50, fin_exp=2.5, samples=samples, z=2.0)
    upper = _cross_loop(1.0, fin_amp=0.35, fin_exp=3.0, samples=samples, z=4.0)
    top = _clean_loop(_circle_loop(0.55, samples=samples))

    zs = [0.0, 2.0, 4.0, 6.0]
    for loop, z in zip([base, mid, upper, top], zs):
        for p in loop:
            p[2] = z
    target = len(base)
    loops = [
        _resample_loop(base, target),
        _resample_loop(mid, target),
        _resample_loop(upper, target),
        _resample_loop(top, target),
    ]
    return loops


def _loft_chain(loops, engine="occ"):
    from yapcad.geom3d import solid_boolean
    from yapcad.geom3d_util import makeLoftSolid
    from yapcad.brep import has_brep_data

    solids = []
    for a, b in zip(loops[:-1], loops[1:]):
        solids.append(makeLoftSolid(a, b))
    result = solids[0]
    for nxt in solids[1:]:
        eng = engine if (has_brep_data(result) and has_brep_data(nxt)) else None
        result = solid_boolean(result, nxt, "union", engine=eng)
    return result


def _build_nose(z_start, radius, length):
    from yapcad.geom3d_util import makeRevolutionSolid
    def contour(z):
        t = max(0.0, min(1.0, (z - z_start) / length))
        return radius * math.sqrt(1 - t)

    return makeRevolutionSolid(contour, z_start, z_start + length, steps=48, arcSamples=96)


def _hollow_body(outer_solid, loops, wall, nose_length, engine="occ"):
    from yapcad.geom3d import solid_boolean
    from yapcad.geom3d_util import tube
    from yapcad.brep import has_brep_data

    inner_loops = []
    for loop in loops:
        r = _loop_radius(loop)
        scale = max(0.05, (r - wall) / r)
        inner_loops.append(_scale_loop(loop, scale))
    inner_body = _loft_chain(inner_loops, engine=engine)
    z_top = max(p[2] for p in inner_loops[-1])
    inner_nose = _build_nose(z_top, _loop_radius(inner_loops[-1]), nose_length * 0.95)
    eng = engine if (has_brep_data(inner_body) and has_brep_data(inner_nose)) else None
    inner = solid_boolean(inner_body, inner_nose, "union", engine=eng)
    eng_diff = engine if (has_brep_data(outer_solid) and has_brep_data(inner)) else None
    return solid_boolean(outer_solid, inner, "difference", engine=eng_diff)


def build_rocket(wall=0.15, nose_length=1.0, engine="occ"):
    from yapcad.geom3d import solid_boolean
    from yapcad.brep import has_brep_data

    loops = _build_body_sections()
    body = _loft_chain(loops, engine=engine)
    z_top = max(p[2] for p in loops[-1])
    nose = _build_nose(z_top, _loop_radius(loops[-1]), nose_length)
    eng = engine if (has_brep_data(body) and has_brep_data(nose)) else None
    outer = solid_boolean(body, nose, "union", engine=eng)
    hollow = _hollow_body(outer, loops, wall, nose_length, engine=engine)
    return hollow


def _view_solid(sld):
    from yapcad.pyglet_drawable import pygletDraw

    viewer = pygletDraw()
    viewer.make_object("rocket", lighting=True, linecolor="white", material="pearl")
    viewer.draw_solid(sld, name="rocket")
    viewer.display()


def _write_package(sld, output):
    from yapcad.package import create_package_from_entities

    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    manifest = create_package_from_entities([sld], out,
                                            name="rocket_loft_demo",
                                            version="0.1")
    print(f"Wrote package to {out} with manifest {manifest}")


def _write_profiles_dxf(output, loops):
    path = Path(output)
    if path.suffix.lower() != ".dxf":
        path = path.with_suffix(".dxf")
    doc = ezdxf.new(dxfversion="R2010", setup=True)
    msp = doc.modelspace()
    colors = [1, 2, 3, 4]
    for idx, loop in enumerate(loops):
        pts2d = [(p[0], p[1]) for p in loop]
        msp.add_lwpolyline(pts2d, dxfattribs={"layer": "profiles", "color": colors[idx % len(colors)]}, close=True)
    doc.saveas(path)
    print(f"Wrote profile DXF to {path}")


def main():
    parser = argparse.ArgumentParser(description="Lofted rocket body demo (BREP-aware).")
    parser.add_argument("--mode", choices=["view", "package", "profiles"], default="view")
    parser.add_argument("--output", default="build/rocket_loft")
    parser.add_argument("--wall", type=float, default=0.15, help="Wall thickness")
    parser.add_argument("--nose-length", type=float, default=1.0)
    parser.add_argument("--engine", default="occ", help="Boolean engine (occ|native|trimesh:manifold)")
    args = parser.parse_args()

    if args.mode == "profiles":
        loops = _build_body_sections()
        _write_profiles_dxf(args.output, loops)
        return

    rocket = build_rocket(wall=args.wall, nose_length=args.nose_length, engine=args.engine)

    if args.mode == "package":
        _write_package(rocket, args.output)
    else:
        _view_solid(rocket)


if __name__ == "__main__":
    main()
