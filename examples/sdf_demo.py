#!/usr/bin/env python3
"""SDF gallery: mesh a set of signed-distance fields to STL and PNG.

Each model is defined purely as a field -- primitives combined with booleans,
blends, offsets and transforms -- and then dual-contoured into an ordinary
yapCAD solid.  Export goes through ``yapcad.io.stl.write_stl``, the same path
any other yapCAD solid takes, which is the point: downstream code needs no
knowledge of fields.

The gallery doubles as a visual regression check.  Every model reports its
volume, whether the mesh is a closed manifold, whether the field is an exact
distance or only bounds one, its Lipschitz constant, and whether the tree is
replayable as exact CSG.  ``lattice-in-a-part`` is deliberately non-manifold
and is kept as an exhibit of dual contouring's one-vertex-per-cell limit.

See ``docs/SDF-DESIGN.md`` for the representation and its phases.

Usage:
    # Write STL meshes and PNG previews to build/sdf-demo
    python examples/sdf_demo.py

    # Somewhere else, at higher resolution
    python examples/sdf_demo.py --output /tmp/gallery --scale 1.5

    # Just one model, meshes only, no rendering
    python examples/sdf_demo.py --only cube-frame --no-render

    # Also write analytic STEP for every model that replays as exact CSG
    python examples/sdf_demo.py --step

    # List what is available
    python examples/sdf_demo.py --list
"""

import argparse
import contextlib
import json
import math
import os
import struct
import sys
import time

from yapcad import sdf
from yapcad.brep import has_brep_data
from yapcad.geom3d import issolidclosed, volumeof
from yapcad.io.step import write_step_analytic
from yapcad.io.stl import write_stl

try:                                    # run as a script, or as a module
    from sdf_preview import render, write_png
except ImportError:                     # pragma: no cover
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sdf_preview import render, write_png


def _bracket_plate():
    """A plate with one bore and two mounting holes: ordinary hard CSG."""
    return sdf.subtract(
        sdf.rounded_box((40.0, 40.0, 10.0), 3.0),
        sdf.translate(sdf.cylinder(6.0, 20.0), (10.0, 0.0, 0.0)),
        sdf.translate(sdf.cylinder(3.0, 20.0), (-12.0, -12.0, 0.0)),
        sdf.translate(sdf.cylinder(3.0, 20.0), (-12.0, 12.0, 0.0)),
    )


def _lattice_in_a_part():
    """A solid skin over a gyroid core, opened at one corner.

    Kept despite being non-manifold: the three-plane corner cut leaves knife
    edges in the shell, a few cells there see the surface twice, and dual
    contouring cannot represent that.  It is the documented limit rather
    than a defect, and it is worth being able to look at.
    """
    skin = sdf.subtract(sdf.rounded_box((26.0, 26.0, 26.0), 3.0),
                        sdf.box(21.0))
    opening = sdf.translate(sdf.box(22.0), (13.0, 13.0, 13.0))
    core = sdf.intersect(sdf.box(22.0), sdf.gyroid(8.0, 0.6))
    return sdf.union(sdf.subtract(skin, opening), core)


#: name, field, resolution, what it demonstrates, camera overrides
MODELS = [
    ("cube-frame",
     sdf.subtract(sdf.box(10.0), sdf.sphere(6.0)), 128,
     "box minus a sphere that breaches all six faces: a genus-5 frame",
     {}),
    ("box",
     sdf.box(8.0), 32,
     "sharp edges and corners preserved exactly by the QEF solve",
     {}),
    ("sphere",
     sdf.sphere(5.0), 72,
     "the smoothness baseline; normals are analytic, not estimated",
     {}),
    ("torus",
     sdf.torus(9.0, 2.5), 96,
     "genus 1, so Euler characteristic 0",
     {"elevation": 42.0}),
    ("bracket-plate",
     _bracket_plate(), 96,
     "hard CSG: rounded plate with three bores, csgExact and STEP-ready",
     {"elevation": 40.0}),
    ("smooth-union",
     sdf.smooth_union(sdf.translate(sdf.sphere(4.5), (-3.5, 0.0, 0.0)),
                      sdf.translate(sdf.sphere(4.5), (3.5, 0.0, 0.0)),
                      2.5), 96,
     "field-native fillet; exactness is lost but the Lipschitz bound is not",
     {}),
    ("shell-cutaway",
     sdf.intersect(
         sdf.shell(sdf.rounded_box((18.0, 18.0, 18.0), 3.0), 2.0),
         sdf.half_space((1.0, 1.0, 0.4), 1.0)), 96,
     "a 2mm wall cut open by a half space -- both are field operations",
     {"azimuth": 60.0}),
    ("gyroid-lattice",
     sdf.intersect(sdf.box(24.0), sdf.gyroid(9.0, 0.55)), 112,
     "infinite lattice bounded by a solid; the node with Lipschitz > 1",
     {}),
    ("lattice-in-a-part",
     _lattice_in_a_part(), 112,
     "solid skin over a lattice core, corner opened; a few knife-edge cells "
     "make this one non-manifold -- the documented dual-contouring limit",
     {"azimuth": 50.0}),
]


def stl_facet_count(path):
    """Read the facet count from a binary STL header."""
    with open(path, "rb") as handle:
        return struct.unpack("<I", handle.read(84)[80:84])[0]


@contextlib.contextmanager
def _quiet_stdout():
    """Silence OCC's STEP transfer banner, which it writes to fd 1 directly."""
    sys.stdout.flush()
    saved = os.dup(1)
    with open(os.devnull, "w") as sink:
        os.dup2(sink.fileno(), 1)
        try:
            yield
        finally:
            sys.stdout.flush()
            os.dup2(saved, 1)
            os.close(saved)


def build(model, output, scale, want_stl, want_render, size,
          want_step=False):
    """Mesh one model, write its artefacts and return a manifest entry."""
    name, node, resolution, note, view = model
    resolution = max(2, int(round(resolution * scale)))

    started = time.time()
    mesh = sdf.dual_contour(node, resolution=resolution)
    elapsed = time.time() - started
    defects = sdf.manifold_defects(mesh)
    sound = not any(defects.values())

    entry = {
        "name": name,
        "note": note,
        "resolution": resolution,
        "vertices": int(len(mesh.vertices)),
        "triangles": int(len(mesh.triangles)),
        "manifold": sound,
        "defects": defects,
        "exact": bool(node.exact),
        "lipschitz": round(float(node.lipschitz), 6),
        "csgExact": bool(node.csg_exact),
        "seconds": round(elapsed, 3),
    }

    # check=False so the non-manifold exhibit still builds; the manifest
    # records the truth either way.  Wrapping the mesh is cheap next to
    # meshing it, so the solid's own measurements are always reported.
    solid = sdf.to_solid(node, resolution=resolution, check=False,
                         brep="auto" if want_step else False)
    entry["closed"] = bool(issolidclosed(solid))
    entry["volume"] = (round(float(volumeof(solid)), 4)
                       if entry["closed"] else None)

    if want_stl:
        path = os.path.join(output, f"{name}.stl")
        write_stl(solid, path)
        facets = stl_facet_count(path)
        entry["stlFacets"] = facets
        # yapcad.mesh skips faces with area <= geom.epsilon, so a mesh with
        # legitimate slivers can lose a few facets on the way out.
        entry["stlDropped"] = entry["triangles"] - facets
        entry["stlBytes"] = os.path.getsize(path)

    # Which nodes, if any, keep this part from replaying as exact CSG --
    # the authoring-time answer to "will this round-trip to STEP?".
    entry["csgBlockers"] = [n.kind for n in sdf.csg_blockers(node)]
    if want_step:
        entry["step"] = None
        if has_brep_data(solid):
            path = os.path.join(output, f"{name}.step")
            with _quiet_stdout():
                analytic = write_step_analytic(solid, path,
                                               fallback_to_faceted=False)
            entry["step"] = "analytic" if analytic else None

    if want_render:
        image = render(mesh.vertices, mesh.normals, mesh.triangles,
                       size=size, **view)
        write_png(os.path.join(output, f"{name}.png"), image)

    return entry


def _step_cell(entry):
    """What the README says about a model's STEP export."""
    if entry.get("step") == "analytic":
        return "analytic"
    if entry["csgBlockers"]:
        return "no: " + ", ".join(sorted(set(entry["csgBlockers"])))
    return "--" if "step" not in entry else "no"


def write_readme(output, manifest):
    """Write a short index describing what landed in the output directory."""
    lines = [
        "# SDF gallery",
        "",
        "Generated by `examples/sdf_demo.py`. Each model is defined as a",
        "signed distance field, dual-contoured to a yapCAD solid, and",
        "exported through `yapcad.io.stl.write_stl` -- the same path any",
        "other yapCAD solid takes.",
        "",
        "| model | res | triangles | volume | exact | L | csgExact "
        "| manifold | facets | dropped | STEP |",
        "|---|--:|--:|--:|:-:|--:|:-:|:-:|--:|--:|---|",
    ]
    for e in manifest:
        volume = e.get("volume")
        lines.append(
            f"| `{e['name']}` | {e['resolution']} | {e['triangles']:,} "
            f"| {'--' if volume is None else format(volume, '.3f')} "
            f"| {'yes' if e['exact'] else 'no'} | {e['lipschitz']:g} "
            f"| {'yes' if e['csgExact'] else 'no'} "
            f"| {'yes' if e['manifold'] else 'NO'} "
            f"| {e.get('stlFacets', 0):,} | {e.get('stlDropped', 0)} "
            f"| {_step_cell(e)} |"
        )
    lines += [
        "",
        "`dropped` counts triangles the STL writer discarded: `yapcad.mesh`",
        "skips any face with area <= `geom.epsilon`, and dual contouring",
        "makes a few slivers that small on curved geometry. Cosmetically",
        "irrelevant, but a strict watertight checker will notice.",
        "",
        "`lattice-in-a-part` is deliberately non-manifold; see the module",
        "docstring of `examples/sdf_demo.py`.",
        "",
    ]
    for e in manifest:
        lines.append(f"- `{e['name']}` -- {e['note']}")
    lines.append("")
    with open(os.path.join(output, "README.md"), "w") as handle:
        handle.write("\n".join(lines))


def main(argv=None):
    names = [m[0] for m in MODELS]
    parser = argparse.ArgumentParser(
        description="Mesh a gallery of SDF models to STL and PNG.")
    parser.add_argument("--output", default=os.path.join("build", "sdf-demo"),
                        help="output directory (default: build/sdf-demo)")
    parser.add_argument("--only", action="append", choices=names,
                        metavar="NAME",
                        help="build just this model; repeatable")
    parser.add_argument("--scale", type=float, default=1.0,
                        help="multiply every model's resolution")
    parser.add_argument("--size", type=int, default=900,
                        help="rendered image edge in pixels")
    parser.add_argument("--no-stl", action="store_true",
                        help="skip STL export")
    parser.add_argument("--no-render", action="store_true",
                        help="skip PNG rendering")
    parser.add_argument("--step", action="store_true",
                        help="also write analytic STEP for models that "
                             "replay as exact CSG (needs pythonocc-core)")
    parser.add_argument("--list", action="store_true",
                        help="list the models and exit")
    args = parser.parse_args(argv)

    if args.list:
        for name, _, resolution, note, _ in MODELS:
            print(f"{name:20s} res={resolution:4d}  {note}")
        return 0
    if args.scale <= 0 or not math.isfinite(args.scale):
        parser.error("--scale must be positive")

    selected = [m for m in MODELS if not args.only or m[0] in args.only]
    os.makedirs(args.output, exist_ok=True)

    manifest = []
    for model in selected:
        entry = build(model, args.output, args.scale,
                      not args.no_stl, not args.no_render, args.size,
                      want_step=args.step)
        manifest.append(entry)
        volume = entry.get("volume")
        shown = "(open)" if volume is None else format(volume, ".3f")
        print(f"{entry['name']:20s} res={entry['resolution']:4d} "
              f"tris={entry['triangles']:7d} vol={shown:>11s} "
              f"csgExact={str(entry['csgExact']):5s} "
              f"manifold={str(entry['manifold']):5s} "
              f"{entry['seconds']:5.2f}s")

    with open(os.path.join(args.output, "manifest.json"), "w") as handle:
        json.dump(manifest, handle, indent=2)
    write_readme(args.output, manifest)
    print(f"\nwrote {len(manifest)} models to {args.output}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
