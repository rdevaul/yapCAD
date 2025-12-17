"""Generate a .ycpkg for threaded fasteners (screws, nuts, or washers)."""

from __future__ import annotations

from dataclasses import  replace
import argparse
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from yapcad.fasteners import (
    HexCapScrewSpec,
    HexNutSpec,
    build_hex_cap_screw,
    build_hex_nut,
    metric_hex_cap_catalog,
    metric_hex_cap_screw,
    metric_hex_nut,
    metric_hex_nut_catalog,
    unified_hex_cap_catalog,
    unified_hex_cap_screw,
    unified_hex_nut,
    unified_hex_nut_catalog,
)
from yapcad.geom import point
from yapcad.geom3d import translate
from yapcad.metadata import get_solid_metadata
from yapcad.geom3d_util import conic
from yapcad.package import create_package_from_entities
from yapcad.threadgen import metric_profile, unified_profile


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create threaded fastener packages (.ycpkg).")
    parser.add_argument(
        "--fastener",
        choices=["screw", "nut", "washer"],
        default="screw",
        help="Which fastener to generate.",
    )
    parser.add_argument(
        "--standard",
        choices=["metric", "unified"],
        default="metric",
        help="Which catalog to use for defaults.",
    )
    parser.add_argument("--size", default="M8", help="Standard designation (e.g. M8 or 1/4-20).")
    parser.add_argument("--length", type=float, help="Screw length (mm for metric, inches for unified).")
    parser.add_argument("--thread-length", type=float, help="Optional screw thread length override.")
    parser.add_argument("--width-flat", type=float, help="Override width across flats for nuts.")
    parser.add_argument("--thickness", type=float, help="Override nut/washer thickness.")
    parser.add_argument("--diameter", type=float, help="Override nominal diameter (mm).")
    parser.add_argument("--pitch", type=float, help="Override thread pitch (mm).")
    parser.add_argument("--tpi", type=float, help="Threads per inch (unified).")
    parser.add_argument("--starts", type=int, default=1, help="Thread starts.")
    parser.add_argument("--handedness", choices=["right", "left"], default="right", help="Thread handedness.")
    parser.add_argument("--output", type=Path, default=Path("threaded_fastener"), help="Base output path.")
    parser.add_argument("--with-nut", action="store_true", help="When generating a screw, include a mated nut positioned mid-thread.")
    parser.add_argument("--name", default=None, help="Package name override.")
    parser.add_argument("--version", default="0.1.0", help="Package version.")
    parser.add_argument("--description", default="Threaded fastener generated via yapCAD", help="Package description.")
    return parser.parse_args()


def _select_catalog_entry(args: argparse.Namespace, kind: str | None = None):
    fastener_kind = kind or args.fastener
    if fastener_kind == "screw":
        catalog = metric_hex_cap_catalog() if args.standard == "metric" else unified_hex_cap_catalog()
    else:
        catalog = metric_hex_nut_catalog() if args.standard == "metric" else unified_hex_nut_catalog()
    if args.size.upper() == "CUSTOM":
        return None
    key = args.size.upper() if args.standard == "metric" else args.size.lower()
    return catalog.get(key)


def _build_screw(args: argparse.Namespace, catalog_entry: dict | None):
    length = args.length if args.length is not None else (25.0 if args.standard == "metric" else 1.0)
    thread_length_value = args.thread_length if args.thread_length is not None else length * 0.7

    if args.standard == "metric":
        if catalog_entry and not any(args.__dict__.get(k) for k in ("diameter", "pitch")):
            return metric_hex_cap_screw(
                args.size,
                length=length,
                thread_length=thread_length_value,
                starts=args.starts,
            )
        diameter = args.diameter or (catalog_entry["diameter"] if catalog_entry else 8.0)
        pitch = args.pitch or (catalog_entry["pitch"] if catalog_entry else 1.25)
        profile = metric_profile(diameter, pitch)
        profile = replace(profile, starts=args.starts)
        length_mm = length
        thread_length_mm = thread_length_value
    else:
        if catalog_entry and not any(args.__dict__.get(k) for k in ("diameter", "pitch", "tpi")):
            return unified_hex_cap_screw(
                args.size,
                length_in=length,
                thread_length_in=thread_length_value,
                starts=args.starts,
            )
        diameter = args.diameter or (catalog_entry["diameter"] if catalog_entry else 6.35)
        tpi = args.tpi or (catalog_entry["tpi"] if catalog_entry else 20.0)
        pitch = 25.4 / tpi
        profile = unified_profile(diameter / 25.4, tpi)
        profile = replace(profile, starts=args.starts)
        length_mm = length * 25.4
        thread_length_mm = thread_length_value * 25.4
    spec = HexCapScrewSpec(
        diameter=diameter,
        thread_length=thread_length_mm,
        shank_length=length_mm,
        head_height=catalog_entry["head_height"] if catalog_entry else diameter * 0.8,
        head_flat_diameter=catalog_entry["head_flat"] if catalog_entry else diameter * 1.6,
        washer_thickness=catalog_entry["washer_thickness"] if catalog_entry else diameter * 0.05,
        washer_diameter=None,
        thread_arc_samples=180,
        thread_samples_per_pitch=6,
    )
    return build_hex_cap_screw(profile, spec)


def _build_nut(args: argparse.Namespace, catalog_entry: dict | None):
    if args.standard == "metric":
        if catalog_entry and not any(args.__dict__.get(k) for k in ("diameter", "pitch", "width_flat", "thickness")):
            return metric_hex_nut(
                args.size,
                starts=args.starts,
                handedness=args.handedness,
            )
        diameter = args.diameter or (catalog_entry["diameter"] if catalog_entry else 8.0)
        pitch = args.pitch or (catalog_entry["pitch"] if catalog_entry else 1.25)
        profile = metric_profile(diameter, pitch, internal=True)
        profile = replace(profile, starts=args.starts)
    else:
        if catalog_entry and not any(args.__dict__.get(k) for k in ("diameter", "tpi", "width_flat", "thickness")):
            return unified_hex_nut(
                args.size,
                starts=args.starts,
                handedness=args.handedness,
            )
        diameter = args.diameter or (catalog_entry["diameter"] if catalog_entry else 6.35)
        tpi = args.tpi or (catalog_entry["tpi"] if catalog_entry else 20.0)
        pitch = 25.4 / tpi
        profile = unified_profile(diameter / 25.4, tpi, internal=True)
        profile = replace(profile, starts=args.starts)
    spec = HexNutSpec(
        diameter=diameter,
        pitch=pitch,
        width_flat=args.width_flat or (catalog_entry["width_flat"] if catalog_entry else diameter * 1.6),
        thickness=args.thickness or (catalog_entry["thickness"] if catalog_entry else diameter * 0.8),
        handedness=args.handedness,
        starts=args.starts
    )
    return build_hex_nut(profile, spec)


def _build_washer(args: argparse.Namespace, catalog_entry: dict | None):
    default_diameter = 8.0 if args.standard == "metric" else 6.35
    diameter = args.diameter or (catalog_entry["diameter"] if catalog_entry else default_diameter)
    width_flat = args.width_flat or (catalog_entry.get("width_flat") if catalog_entry else diameter * 1.5)
    if width_flat is None:
        width_flat = diameter * 1.5
    thickness = args.thickness or (catalog_entry.get("thickness") if catalog_entry else diameter * 0.5)
    if thickness is None:
        thickness = diameter * 0.5
    return conic(width_flat / 2.0, width_flat / 2.0, thickness)


def _build_fastener(args: argparse.Namespace):
    catalog_entry = _select_catalog_entry(args)
    if args.fastener == "screw":
        screw = _build_screw(args, catalog_entry)
        if args.with_nut:
            nut = _build_nut(args, _select_catalog_entry(args, kind="nut"))
            nut_meta = get_solid_metadata(nut, create=False)
            screw_meta = get_solid_metadata(screw, create=False)
            thread_length = screw_meta["hex_cap_screw"]["thread_length"]
            nut_thickness = nut_meta["hex_nut"]["thickness"]
            offset = thread_length / 2.0 - nut_thickness / 2.0
            nut = translate(nut, point(0, 0, offset))
            return [screw, nut]
        return [screw]
    if args.fastener == "nut":
        return [_build_nut(args, catalog_entry)]
    return [_build_washer(args, catalog_entry)]


def main():
    args = _parse_args()
    if args.with_nut and args.fastener != "screw":
        raise ValueError("--with-nut can only be used when --fastener screw")
    fasteners = _build_fastener(args)
    pkg_root = args.output.with_suffix(".ycpkg")
    if pkg_root.exists():
        shutil.rmtree(pkg_root)
    kind = args.fastener
    name = args.name or f"{args.standard}-{args.size}-{kind}".replace("/", "_")
    manifest = create_package_from_entities(
        fasteners,
        pkg_root,
        name=name,
        version=args.version,
        units="mm",
        description=args.description,
    )
    manifest.save()
    print(f"Wrote package {pkg_root}")


if __name__ == "__main__":
    main()
