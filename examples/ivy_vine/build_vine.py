"""Build, check and export Ivy's Vine.

    python build_vine.py                 # everything -> out/ivy_vine.stl + reports
    python build_vine.py --quick         # no STL: pictures, clearances, sweeps
    python build_vine.py --links 2 --max-level 0

The design itself is all in ivy_vine.dsl.  This script only drives it:

  1. runs IVY_VINE with representation="sdf" and no meshing, so every part
     is a field and the build takes a fraction of a second;
  2. measures the as-printed clearance between every pair of parts on their
     fields (yapcad.rom) -- a print-in-place part must not fuse;
  3. sweeps every joint across its declared limits and reports the closest
     approach of each jointed pair;
  4. sphere-traces pictures of the assembly straight from the fields;
  5. meshes and simplifies the parts in parallel and writes the STL.

Run it from this directory.
"""

import argparse
import json
import math
import pathlib
import time

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent


def build(params):
    """Run IVY_VINE as fields, without meshing anything.

    Clearances, sweeps and pictures all read the fields, so the design is
    built deferred (yapCAD's ``sdf_mesh=False``): the parts come back as a
    field and a meshing plan each, in well under a second.  Only writing the
    STL needs triangles, and finish_parts() builds those.
    """
    from yapcad.dsl import compile_and_run
    path = HERE / "ivy_vine.dsl"
    started = time.time()
    result = compile_and_run(
        path.read_text(), "IVY_VINE", params, representation="sdf",
        sdf_mesh=False, source_path=str(path))
    if not result.success:
        raise SystemExit(f"build failed: {result.error_message}")
    print(f"built IVY_VINE{params} in {time.time() - started:.1f}s (deferred)")
    return result


def _finish_one(solid):
    """Mesh one deferred part and simplify the mesh against its field."""
    from yapcad import sdf
    sdf.realize(solid, simplify=False)
    return sdf.simplify_solid(solid)


def finish_parts(assembly, workers):
    """Mesh and simplify every part of the assembly, in parallel.

    Simplifying a mesh against its field is the slow step -- minutes a part
    -- and the parts are independent, so they go to a process pool.  Links
    built from the same parameters have the same field; each distinct field
    is finished once and shared.
    """
    import multiprocessing
    from yapcad import sdf
    from yapcad.construction import get_construction

    started = time.time()
    by_field = {}
    for name, solid in assembly.geometry.items():
        key = sdf.from_construction(get_construction(solid)).digest
        by_field.setdefault(key, []).append(name)
    firsts = [names[0] for names in by_field.values()]
    print(f"\nfinishing {len(assembly.geometry)} parts "
          f"({len(firsts)} distinct) on {workers} processes ...", flush=True)
    with multiprocessing.get_context("spawn").Pool(workers) as pool:
        done = pool.map(_finish_one,
                        [assembly.geometry[name] for name in firsts],
                        chunksize=1)
    for names, finished in zip(by_field.values(), done):
        for name in names:
            # In place: the assembly and its compound hold this same list.
            assembly.geometry[name][:] = finished
    print(f"finished in {time.time() - started:.0f}s")
    return assembly.compound_geometry()


def world_fields(model, poses):
    from yapcad import sdf
    return {name: sdf.transform(model.field(name), poses[name].tolist())
            for name in model.parts}


def rest_clearances(model, poses, pitch):
    """Clearance of every pair of parts in the printed pose."""
    from yapcad.rom import bounds_clearance, clearance_sdf
    rows = []
    names = list(model.parts)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            fa, fb = model.field(a), model.field(b)
            if bounds_clearance(fa, fb, poses[a], poses[b]) > 2.0:
                continue
            # Sample the smaller part, evaluate the larger.
            if _extent(fa) < _extent(fb):
                a_, b_, fa, fb = b, a, fb, fa
            else:
                a_, b_ = a, b
            c = clearance_sdf(fa, fb, poses[a_], poses[b_], pitch_mm=pitch,
                              max_samples=4_000_000)
            rows.append((c.min_distance, c.tolerance_mm, a, b))
    return sorted(rows)


def joint_sweeps(model, rest, step_deg=10.0, pitch=0.25):
    """Each joint through its declared limits, the others at rest.

    Only the jointed pair is measured: this is the check that a hinge
    turns freely over the range the design declares for it.
    """
    from yapcad.rom import SweepConfig, sweep_rom
    rows = []
    print("\nJoint sweeps (each joint alone, across its declared limits):")
    for joint in model.movable_joints():
        lo, hi = math.degrees(joint.lower), math.degrees(joint.upper)
        fixed = {k: v for k, v in rest.items() if k != joint.name}
        # Evenly spaced, both limits included, held a hair inside them so
        # that rounding cannot put a sample outside the declared range.
        count = max(2, int(round((hi - lo) / step_deg)) + 1)
        values = [math.radians(v)
                  for v in np.linspace(lo + 1e-6, hi - 1e-6, count)]
        sweep = sweep_rom(
            model, {joint.name: values}, fixed=fixed,
            pairs=[(joint.parent, joint.child)],
            config=SweepConfig(pitch_mm=pitch, clearance_threshold=0.3,
                               max_samples=2_000_000))
        trace = sweep.trace(joint.parent, joint.child)
        clear = [math.degrees(s.config[joint.name]) for s in trace
                 if s.clearance.min_distance > 0.0]
        worst = min(trace, key=lambda s: s.clearance.min_distance)
        row = {"joint": joint.name, "limits_deg": [round(lo, 1), round(hi, 1)],
               "min_clearance_mm": round(worst.clearance.min_distance, 3),
               "at_deg": round(math.degrees(worst.config[joint.name]), 1),
               "clear_from_deg": round(min(clear), 1) if clear else None,
               "clear_to_deg": round(max(clear), 1) if clear else None,
               "steps": len(trace)}
        rows.append(row)
        print(f"  {joint.name:8s} {lo:7.1f} .. {hi:6.1f} deg   "
              f"closest {row['min_clearance_mm']:+.3f} mm at "
              f"{row['at_deg']:.0f} deg   clear at "
              f"{len(clear)}/{len(trace)} steps")
    return rows


def export_stl(solid, path, min_volume=1.0):
    """Write the assembly as one STL of separate bodies, one per part.

    Meshing a thin feature can leave a speck of surface a fraction of a
    millimetre across, detached from its part.  A slicer would ignore it,
    but it has no business in the file: bodies under ``min_volume`` mm^3 are
    dropped, and the report says how many.
    """
    import trimesh
    from check_print import body_labels

    vertices, faces, base = [], [], 0
    for surface in solid[1]:
        vertices.extend(v[:3] for v in surface[1])
        faces.extend([a + base, b + base, c + base] for a, b, c in surface[3])
        base += len(surface[1])
    mesh = trimesh.Trimesh(np.asarray(vertices, float),
                           np.asarray(faces, np.int64), process=True)
    labels = body_labels(mesh)
    tris = mesh.vertices[mesh.faces]
    signed = np.einsum("ij,ij->i", tris[:, 0],
                       np.cross(tris[:, 1], tris[:, 2])) / 6.0
    volumes = np.bincount(labels, weights=signed)
    keep = np.abs(volumes) >= min_volume
    mesh.update_faces(keep[labels])
    mesh.remove_unreferenced_vertices()
    mesh.export(str(path))
    lo, hi = mesh.bounds
    info = {"triangles": int(len(mesh.faces)), "bodies": int(keep.sum()),
            "specks_dropped": int((~keep).sum()),
            "size_mm": [round(float(v), 2) for v in hi - lo],
            "z_min": round(float(lo[2]), 4),
            "volume_cm3": round(float(np.abs(volumes[keep]).sum()) / 1000.0, 2)}
    print(f"\nwrote {path}  ({info['triangles']} triangles in "
          f"{info['bodies']} bodies, {info['size_mm'][0]} x "
          f"{info['size_mm'][1]} x {info['size_mm'][2]} mm, "
          f"{info['volume_cm3']} cm^3; dropped {info['specks_dropped']} "
          f"speck(s))")
    return info


def _extent(field):
    lo, hi = field.bounds
    return float(np.prod([hi[i] - lo[i] for i in range(3)]))


def pictures(fields, out, tag):
    from yapcad import render
    started = time.time()
    render.image(fields, view="top", pixels=(1200, 900),
                 path=str(out / f"{tag}_top.png"), colour_by="none")
    render.image(fields, view="iso", pixels=(1200, 800),
                 path=str(out / f"{tag}_iso.png"), colour_by="part")
    print(f"pictures in {time.time() - started:.0f}s")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--quick", action="store_true",
                    help="no STL: pictures, clearances and sweeps only")
    ap.add_argument("--workers", type=int, default=6,
                    help="processes used to mesh and simplify the parts")
    ap.add_argument("--links", type=int, default=None)
    ap.add_argument("--max-level", type=int, default=None)
    ap.add_argument("--swing", type=float, default=None)
    ap.add_argument("--wander", type=float, default=None)
    ap.add_argument("--bud", type=float, default=None)
    ap.add_argument("--step", type=float, default=10.0,
                    help="joint sweep step in degrees")
    ap.add_argument("--no-sweep", action="store_true")
    ap.add_argument("--tag", default="ivy_vine")
    args = ap.parse_args()

    params = {}
    for key, value in (("links", args.links), ("max_level", args.max_level),
                       ("swing_deg", args.swing), ("wander_deg", args.wander),
                       ("bud_deg", args.bud)):
        if value is not None:
            params[key] = value

    out = HERE / "out"
    out.mkdir(exist_ok=True)
    result = build(params)

    from yapcad.rom import RomModel, SweepConfig, sweep_rom
    model = RomModel.from_result(result)
    print(model.describe())
    rest = model.joint_values()
    poses = model.pose(rest)

    report = {"params": params, "parts": len(model.parts)}

    rows = rest_clearances(model, poses, pitch=0.2)
    print("\nAs-printed clearance, closest pairs (mm):")
    for dist, tol, a, b in rows[:12]:
        print(f"  {dist:+.3f} (+/- {tol:.3f})  {a} - {b}")
    report["rest_clearance"] = [
        {"a": a, "b": b, "mm": round(d, 4), "tolerance_mm": round(t, 4)}
        for d, t, a, b in rows]

    pictures(world_fields(model, poses), out, args.tag)

    if not args.no_sweep:
        report["sweep"] = joint_sweeps(model, rest, step_deg=args.step)
        model.pose(rest)

    if not args.quick:
        stl = out / f"{args.tag}.stl"
        report["stl"] = export_stl(
            finish_parts(model.assembly, args.workers), stl)
    (out / f"{args.tag}_report.json").write_text(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
