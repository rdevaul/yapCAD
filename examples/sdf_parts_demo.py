#!/usr/bin/env python3
"""Real parts built as signed distance fields, with no OpenCASCADE.

``sdf_demo.py`` shows the field primitives and operators.  This shows what
they are for: the parts a mechanism needs, built the ways a yapCAD user
would build them.

``yaprover-wheel``
    The yapRover wheel, written in the DSL exactly as the rover's design
    writes it, and built with ``representation="sdf"``.  Its filleted rim is
    an exact field, so no OCC is needed for the fillet.
``gyroid-wheel``
    The same wheel with a gyroid lattice in place of the spokes, blended
    into the rim and hub: a part that is easy as a field and hard as a BREP.
``m8-nut``, ``m8x30-bolt``
    ISO metric fasteners with helical threads, from ``yapcad.fasteners``.
``miter-gear``, ``herringbone-gear``
    A straight bevel gear and a double-helical involute gear, from
    ``yapcad.gears``.

Every part is meshed by dual contouring and then simplified against its
own field.  yapCAD's builders simplify finished parts by default; this demo
asks them for the uniform mesh and simplifies it itself
(``sdf.simplify_solid``), so the report can give both triangle counts.

Usage:
    # STL meshes, PNG previews and a summary in build/sdf-parts
    python examples/sdf_parts_demo.py

    # One part, no simplification, somewhere else
    python examples/sdf_parts_demo.py --only m8-nut --tolerance 0 \
        --output /tmp/parts

    # What is available
    python examples/sdf_parts_demo.py --list
"""

import argparse
import json
import os
import sys
import time

from yapcad import sdf
from yapcad.dsl import compile_and_run
from yapcad.fasteners import metric_hex_bolt, metric_hex_nut
from yapcad.gears import make_involute_gear_sdf, make_straight_bevel_gear_sdf
from yapcad.gears.bevel import StraightBevelGearSpec
from yapcad.geom3d import issolidclosed, volumeof
from yapcad.io.stl import write_stl

try:                                    # run as a script, or as a module
    from sdf_preview import render, write_png
except ImportError:                     # pragma: no cover
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sdf_preview import render, write_png

#: PETG, the yapRover print material, in g/mm^3.
PETG_DENSITY = 1.27e-3

#: The yapRover wheel, as yaprover_suspension_detailed.dsl writes it.
WHEEL_DSL = """
module wheel_demo

command WHEEL_HUB() -> solid:
    let outer: solid = fillet(cylinder(65.0, 36.0), 2.5)
    let rim_cutter: solid = translate(cylinder(55.0, 38.0), 0.0, 0.0, -1.0)
    let rim: solid = difference(outer, rim_cutter)
    let hub: solid = cylinder(18.0, 36.0)
    let spoke: solid = translate(box(112.0, 8.0, 8.0), 0.0, 0.0, 14.0)
    let angles: list<float> = [0.0, 45.0, 90.0, 135.0]
    let spokes: list<solid> = [rotate(spoke, 0.0, 0.0, a) for a in angles]
    let blank: solid = union_all([rim, hub] + spokes)
    let left_seat: solid = translate(cylinder(11.075, 7.4), 0.0, 0.0, -0.1)
    let right_seat: solid = translate(cylinder(11.075, 7.4),
                                      0.0, 0.0, 28.7)
    let bore: solid = translate(cylinder(6.25, 36.4), 0.0, 0.0, -0.2)
    emit difference_all(blank, [left_seat, right_seat, bore])
"""


def yaprover_wheel():
    result = compile_and_run(WHEEL_DSL, "WHEEL_HUB", {},
                             representation="sdf", sdf_cell_mm=0.5,
                             sdf_simplify=False)
    if not result.success:
        raise RuntimeError(result.error_message)
    geometry = result.geometry
    if isinstance(geometry, list) and geometry and geometry[0] != "solid":
        geometry = geometry[0]
    return geometry


def _zcyl(radius, height, z0):
    """A cylinder standing on ``z = z0``, as the DSL's ``cylinder`` does."""
    return sdf.translate(sdf.cylinder(radius, height), (0, 0, z0 + height / 2))


def gyroid_wheel():
    """The yapRover wheel with a gyroid web in place of its four spokes."""
    outer = sdf.translate(sdf.rounded_cylinder(65.0, 36.0, 2.5), (0, 0, 18.0))
    rim = sdf.subtract(outer, _zcyl(55.0, 38.0, -1.0))
    hub = _zcyl(18.0, 36.0, 0.0)
    band = sdf.subtract(_zcyl(56.0, 28.0, 4.0), _zcyl(17.0, 36.0, 0.0))
    web = sdf.intersect(sdf.gyroid(16.0, 0.6), band)
    body = sdf.smooth_union(sdf.smooth_union(rim, web, 1.5), hub, 1.5)
    seats = [_zcyl(11.075, 7.4, -0.1), _zcyl(11.075, 7.4, 28.7),
             _zcyl(6.25, 36.4, -0.2)]
    return sdf.to_solid(sdf.subtract(body, *seats), resolution=260)


def miter_gear():
    return make_straight_bevel_gear_sdf(StraightBevelGearSpec(
        teeth=24, mate_teeth=24, outer_module_mm=1.5, face_width_mm=8.0,
        backlash_mm=0.25, bore_diameter_mm=8.4), simplify=False)


def herringbone_gear():
    return make_involute_gear_sdf(20, 1.5, 8.0, helix_angle_deg=25.0,
                                  herringbone=True, simplify=False)


#: name, builder, what it shows, camera overrides
PARTS = [
    ("yaprover-wheel", yaprover_wheel,
     "the rover's own DSL, built with representation='sdf'; exact fillets, "
     "no OCC", {"elevation": 40.0}),
    ("gyroid-wheel", gyroid_wheel,
     "a gyroid web blended into rim and hub in place of the spokes",
     {"elevation": 40.0}),
    ("m8-nut", lambda: metric_hex_nut("M8", representation="sdf",
                                      simplify=False),
     "ISO 4032 M8 nut with a helical thread field",
     {"elevation": 50.0, "zoom": 3.4}),
    ("m8x30-bolt", lambda: metric_hex_bolt("M8", 30.0, representation="sdf",
                                           simplify=False),
     "M8x30 hex bolt: external thread, shank, washer face, head",
     {"zoom": 3.6}),
    ("miter-gear", miter_gear,
     "24-tooth straight bevel gear: the outer tooth section carried to the "
     "pitch apex", {"elevation": 45.0}),
    ("herringbone-gear", herringbone_gear,
     "20-tooth involute gear twisted 25 degrees each way from mid face",
     {"elevation": 30.0}),
]


def _arrays(solid):
    import numpy as np
    vertices, normals, triangles, offset = [], [], [], 0
    for surface in solid[1]:
        vertices.append(np.array([p[:3] for p in surface[1]], dtype=float))
        normals.append(np.array([n[:3] for n in surface[2]], dtype=float))
        triangles.append(np.array(surface[3], dtype=np.int64) + offset)
        offset += len(surface[1])
    return np.vstack(vertices), np.vstack(normals), np.vstack(triangles)


def _triangles(solid):
    return sum(len(surface[3]) for surface in solid[1])


def build(part, output, tolerance, want_render, size):
    name, make, note, view = part
    started = time.time()
    solid = make()
    meshed = time.time() - started
    entry = {"name": name, "note": note, "meshSeconds": round(meshed, 2),
             "uniformTriangles": _triangles(solid)}
    if tolerance > 0.0:
        started = time.time()
        solid = sdf.simplify_solid(solid, tolerance)
        entry["simplifySeconds"] = round(time.time() - started, 2)
    entry["triangles"] = _triangles(solid)
    entry["closed"] = bool(issolidclosed(solid))
    if entry["closed"]:
        volume = float(volumeof(solid))
        entry["volume_mm3"] = round(volume, 1)
        entry["mass_g_petg"] = round(volume * PETG_DENSITY, 1)
    write_stl(solid, os.path.join(output, f"{name}.stl"))
    if want_render:
        vertices, normals, triangles = _arrays(solid)
        write_png(os.path.join(output, f"{name}.png"),
                  render(vertices, normals, triangles, size=size, **view))
    return entry


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--output", default=os.path.join("build", "sdf-parts"))
    parser.add_argument("--only", action="append", metavar="NAME",
                        help="build only these parts (repeatable)")
    parser.add_argument("--tolerance", type=float, default=0.01,
                        help="simplification tolerance in mm; 0 to skip "
                             "(default 0.01)")
    parser.add_argument("--no-render", action="store_true",
                        help="skip the PNG previews")
    parser.add_argument("--size", type=int, default=720,
                        help="rendered image edge in pixels")
    parser.add_argument("--list", action="store_true",
                        help="list the parts and exit")
    args = parser.parse_args(argv)

    if args.list:
        for name, _, note, _ in PARTS:
            print(f"{name:18s} {note}")
        return 0
    chosen = [p for p in PARTS if not args.only or p[0] in args.only]
    unknown = set(args.only or ()) - {p[0] for p in PARTS}
    if unknown:
        parser.error(f"unknown part(s): {', '.join(sorted(unknown))}")

    os.makedirs(args.output, exist_ok=True)
    report = []
    for part in chosen:
        entry = build(part, args.output, args.tolerance, not args.no_render,
                      args.size)
        report.append(entry)
        mass = (f"{entry['mass_g_petg']:7.1f} g PETG"
                if "mass_g_petg" in entry else "   (open)")
        print(f"{entry['name']:18s} {entry['uniformTriangles']:8,d} -> "
              f"{entry['triangles']:7,d} triangles  {mass}  "
              f"closed={entry['closed']}", flush=True)
    with open(os.path.join(args.output, "parts.json"), "w") as handle:
        json.dump(report, handle, indent=2)
    return 0 if all(e["closed"] for e in report) else 1


if __name__ == "__main__":
    sys.exit(main())
