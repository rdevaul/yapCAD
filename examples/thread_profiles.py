"""Thread profile demo: preview threads or export them (STL/STEP/.ycpkg)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from yapcad.geom import point
from yapcad.geom3d import solid
from yapcad.geom3d_util import makeRevolutionThetaSamplingSurface
from yapcad.io.step import write_step
from yapcad.io.stl import write_stl
from yapcad.metadata import add_tags, get_solid_metadata, set_layer
from yapcad.package import create_package_from_entities
from yapcad.threadgen import ThreadProfile, metric_profile, sample_thread_profile, unified_profile


def _profile_from_args(args: argparse.Namespace) -> ThreadProfile:
    if args.standard == "metric":
        return metric_profile(args.diameter, args.pitch, internal=args.internal)
    return unified_profile(args.diameter, args.tpi, internal=args.internal)


def _build_surface(profile: ThreadProfile, length: float, arc_samples: int):
    def contour(z0: float, z1: float, theta: float):
        pts, wrap = sample_thread_profile(profile, z0, z1, theta)
        return (pts, wrap)

    return makeRevolutionThetaSamplingSurface(
        contour,
        0.0,
        length,
        arcSamples=arc_samples,
        endcaps=False,
        degrees=True,
    )


def build_thread_solid(profile: ThreadProfile, length: float, arc_samples: int = 180):
    surface = _build_surface(profile, length, arc_samples)
    body = solid([surface])
    meta = get_solid_metadata(body, create=True)
    set_layer(meta, "thread")
    add_tags(meta, ["thread", profile.handedness, profile.starts])
    return body


def export(thread, args):
    if args.export == "stl":
        write_stl(thread, args.output.with_suffix(".stl"))
        print(f"Wrote {args.output.with_suffix('.stl')}")
    elif args.export == "step":
        write_step(thread, args.output.with_suffix(".step"))
        print(f"Wrote {args.output.with_suffix('.step')}")
    elif args.export == "ycpkg":
        pkg_root = args.output.with_suffix(".ycpkg")
        if pkg_root.exists():
            import shutil

            shutil.rmtree(pkg_root)
        manifest = create_package_from_entities(
            [thread],
            pkg_root,
            name="thread_profile",
            version="0.1.0",
            units="mm",
            description="Generated thread profile",
        )
        manifest.save()
        print(f"Wrote package {pkg_root}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Thread profile demo")
    parser.add_argument("--standard", choices=["metric", "unified"], default="metric")
    parser.add_argument("--diameter", type=float, default=10.0,
                        help="Nominal diameter (mm for metric, inches for unified)")
    parser.add_argument("--pitch", type=float, default=1.5, help="Pitch in mm (metric only)")
    parser.add_argument("--tpi", type=float, default=20.0, help="Threads per inch (unified only)")
    parser.add_argument("--length", type=float, default=15.0, help="Threaded length in mm")
    parser.add_argument("--internal", default=False,  action="store_true", help="Generate internal/tapped profile")
    parser.add_argument("--arc-samples", type=int, default=180, help="Angular samples for revolution")
    parser.add_argument("--mode", choices=["view", "stl", "step", "ycpkg"], default="view")
    parser.add_argument("--output", type=Path, default=Path("thread_profile"),
                        help="Output path (for export modes)")
    return parser.parse_args()


def main():
    args = parse_args()
    profile = _profile_from_args(args)
    thread = build_thread_solid(profile, args.length, args.arc_samples)
    args.export = args.mode
    export(thread, args)


if __name__ == "__main__":
    main()
