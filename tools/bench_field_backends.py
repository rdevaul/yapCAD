#!/usr/bin/env python3
"""Time every field backend this machine has, on the same work.

    python tools/bench_field_backends.py
    python tools/bench_field_backends.py --links 400 --size 1280x960

The part is a double helix of blended capsules, the shape that makes the
numpy backend slow: hundreds of nodes and no hard union to cull.  For each
backend that is available it reports

* raw field evaluation, in nanoseconds a point;
* one preview frame, and how many pixels differ from numpy's;
* for the double-precision backends, a dual-contour mesh, and how far its
  vertices sit from numpy's.

numpy is the reference.  ``numba`` needs ``pip install "yapcad[fast]"`` and
``wgpu`` needs ``pip install "yapcad[gpu]"`` and a GPU.  The GPU backend is
single precision and is timed for previews only.
"""

import argparse
import math
import platform
import sys
import time
import warnings

import numpy as np

from yapcad import render, sdf
from yapcad.render.api import Part
from yapcad.sdf import tape
from yapcad.sdf.contour import dual_contour
from yapcad.sdf.evaluate import evaluate


def helix(links, radius=14.0, pitch=18.0, per_turn=28):
    """Two strands of ``links`` blended capsules each, with rungs."""
    def at(i, phase):
        a = 2.0 * math.pi * i / per_turn + phase
        return (radius * math.cos(a), radius * math.sin(a),
                pitch * i / per_turn)
    body = None
    for phase in (0.0, math.pi):
        for i in range(links):
            link = sdf.capsule(at(i, phase), at(i + 1, phase), 1.6)
            body = link if body is None else sdf.smooth_union(body, link, 1.0)
    for i in range(0, links + 1, 4):
        body = sdf.smooth_union(
            body, sdf.capsule(at(i, 0.0), at(i, math.pi), 0.9), 1.0)
    base = sdf.translate(sdf.cylinder(radius + 4.0, 3.0), (0.0, 0.0, -1.5))
    return sdf.smooth_union(body, base, 1.0)


def backends():
    """``(name, note)`` for each backend that works here, numpy first."""
    found = [("numpy", "reference")]
    if tape.numba_available():
        import numba
        found.append(("numba", f"numba {numba.__version__}, "
                               f"{numba.get_num_threads()} threads"))
    else:
        print('numba: not installed (pip install "yapcad[fast]")')
    try:
        from yapcad.sdf import tape_wgpu
        if tape_wgpu.available():
            found.append(("wgpu", tape_wgpu.adapter_summary()))
        else:
            print("wgpu: no usable GPU adapter")
    except Exception as exc:                # an optional backend
        print(f"wgpu: unavailable ({exc})")
    return found


def timed(fn, repeat):
    best, value = math.inf, None
    for _ in range(repeat):
        started = time.perf_counter()
        value = fn()
        best = min(best, time.perf_counter() - started)
    return best, value


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--links", type=int, default=112,
                        help="capsules per strand (default: %(default)s, "
                             "about 500 nodes)")
    parser.add_argument("--size", default="640x480", metavar="WxH")
    parser.add_argument("--points", type=int, default=400000)
    parser.add_argument("--cells", type=int, default=160,
                        help="mesh resolution (default: %(default)s)")
    parser.add_argument("--repeat", type=int, default=2)
    args = parser.parse_args(argv)
    width, height = (int(v) for v in args.size.lower().split("x"))
    sys.setrecursionlimit(max(sys.getrecursionlimit(), 20 * args.links))
    warnings.simplefilter("ignore")

    node = helix(args.links)
    flat = tape.flatten(node)
    print(f"{platform.platform()}, Python {platform.python_version()}, "
          f"numpy {np.__version__}")
    print(f"part: {len(flat.ops)} tape instructions, "
          f"{flat.registers} registers")
    found = backends()
    for name, note in found:
        print(f"  {name}: {note}")

    lo, hi = (np.asarray(v) for v in node.bounds)
    points = lo + np.random.default_rng(0).random((args.points, 3)) * (hi - lo)
    scene = [Part(node, name="helix")]
    reference = {}
    print(f"\n{'backend':8} {'evaluate':>14} {'vs numpy':>9}   "
          f"{'frame ' + args.size:>16} {'pixels off':>11}   "
          f"{'mesh':>8} {'vertex shift':>13}")
    for name, _note in found:
        evaluate(node, points[:64], backend=name)        # load the kernel
        t_eval, values = timed(
            lambda: evaluate(node, points, backend=name), args.repeat)
        t_frame, frame = timed(
            lambda: render.image(scene, view="iso", pixels=(width, height),
                                 backend=name), args.repeat)
        if name == "numpy":
            reference.update(values=values, frame=frame)
        worst = float(np.abs(values - reference["values"]).max())
        off = int((frame != reference["frame"]).any(axis=2).sum())
        row = (f"{name:8} {1e9 * t_eval / args.points:9.0f} ns/pt "
               f"{worst:9.1e}   {t_frame:14.2f} s {off:11d}   ")
        if name == "wgpu":
            row += f"{'--':>8} {'previews only':>13}"
        else:
            t_mesh, mesh = timed(
                lambda: dual_contour(node, resolution=args.cells,
                                     backend=name), 1)
            if name == "numpy":
                reference["mesh"] = mesh
            same = (mesh.triangles.shape == reference["mesh"].triangles.shape
                    and np.array_equal(mesh.triangles,
                                       reference["mesh"].triangles))
            shift = (f"{np.abs(mesh.vertices - reference['mesh'].vertices).max():.1e}"
                     if same else "topology differs")
            row += f"{t_mesh:6.2f} s {shift:>13}"
        print(row)
    return 0


if __name__ == "__main__":
    sys.exit(main())
