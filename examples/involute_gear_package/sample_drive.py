"""Create a simple two-gear drive assembly using the involute gear helper."""

from __future__ import annotations

import argparse
from pathlib import Path
from copy import deepcopy
from typing import Optional, Tuple

from yapcad.metadata import get_solid_metadata, set_layer, add_tags
from yapcad.package import create_package_from_entities, PackageManifest, load_geometry
from yapcad.geom3d import issolid

from .involute_gear import generate_involute_spur, position_gear


def _load_solid_from_package(pkg_path: Path) -> Tuple[list, dict]:
    manifest = PackageManifest.load(pkg_path)
    entities = load_geometry(manifest)
    solids = [geom for geom in entities if issolid(geom)]
    if not solids:
        raise ValueError(f"package {pkg_path} did not contain any solids")
    solid = deepcopy(solids[0])
    meta = get_solid_metadata(solid, create=False) or {}
    return solid, meta


def _extract_invocation(meta: dict) -> dict:
    invocation = meta.get("invocation") if meta else None
    if isinstance(invocation, dict):
        params = invocation.get("parameters") or {}
        if isinstance(params, dict):
            return params
    return {}


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a simple two-gear drive package.")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--driver-teeth", type=int, default=18)
    parser.add_argument("--driven-teeth", type=int, default=36)
    parser.add_argument("--module-mm", type=float, default=1.5)
    parser.add_argument("--face-width-mm", type=float, default=6.0)
    parser.add_argument(
        "--driver-package",
        type=Path,
        help="Optional .ycpkg providing a canonical driver gear to reuse.",
    )
    parser.add_argument(
        "--driven-package",
        type=Path,
        help="Optional .ycpkg providing a canonical driven gear to reuse.",
    )
    args = parser.parse_args()

    driver_meta: Optional[dict] = None
    driven_meta: Optional[dict] = None

    if args.driver_package:
        driver_solid, driver_meta = _load_solid_from_package(args.driver_package)
    else:
        _, driver_solid = generate_involute_spur(
            teeth=args.driver_teeth,
            module_mm=args.module_mm,
            face_width_mm=args.face_width_mm,
        )

    if args.driven_package:
        driven_solid, driven_meta = _load_solid_from_package(args.driven_package)
    else:
        _, driven_solid = generate_involute_spur(
            teeth=args.driven_teeth,
            module_mm=args.module_mm,
            face_width_mm=args.face_width_mm,
        )

    driver_params = _extract_invocation(driver_meta) if driver_meta else {}
    driven_params = _extract_invocation(driven_meta) if driven_meta else {}
    driver_teeth = driver_params.get("teeth", args.driver_teeth)
    driven_teeth = driven_params.get("teeth", args.driven_teeth)
    driver_module = driver_params.get("module_mm", args.module_mm)
    driven_module = driven_params.get("module_mm", args.module_mm)

    if abs(driver_module - driven_module) > 1e-9:
        raise ValueError(
            f"driver module ({driver_module}) does not match driven module ({driven_module}); "
            "gears must share the same module to mesh."
        )
    module_mm = driver_module

    center_distance = module_mm * (driver_teeth + driven_teeth) / 2.0
    placed_driver = position_gear(driver_solid, centre=(0.0, 0.0, 0.0))
    placed_driven = position_gear(driven_solid, centre=(center_distance, 0.0, 0.0))

    driver_meta = get_solid_metadata(placed_driver, create=True)
    set_layer(driver_meta, "gear-driver")
    add_tags(driver_meta, ["gear", "driver"])

    driven_meta = get_solid_metadata(placed_driven, create=True)
    set_layer(driven_meta, "gear-driven")
    add_tags(driven_meta, ["gear", "driven"])

    manifest = create_package_from_entities(
        [
            {"geometry": placed_driver, "metadata": driver_meta},
            {"geometry": placed_driven, "metadata": driven_meta},
        ],
        args.output,
        name="Simple Gear Drive",
        version="0.1.0",
        description="Two-gear drive referencing involute gear helpers",
        units="mm",
        overwrite=True,
    )

    manifest.save()
    print(f"Assembly package written to {manifest.manifest_path}")


if __name__ == "__main__":
    main()
