"""Build a .ycpkg package containing a canonical involute spur gear."""

from __future__ import annotations

import argparse
from pathlib import Path
import shutil
from typing import Any, Dict

from yapcad.metadata import (
    add_tags,
    get_solid_metadata,
    get_surface_metadata,
    set_layer,
)
from yapcad.package import create_package_from_entities
from yapcad.package.core import _compute_hash

from .involute_gear import generate_involute_spur


def _add_invocation(meta: Dict[str, Any], command: str, params: Dict[str, Any]) -> None:
    meta["invocation"] = {
        "package": "involute_gear",
        "command": command,
        "version": "0.1.0",
        "parameters": params,
    }


def _copy_into_package(src: Path, dest_rel: str, pkg_root: Path) -> Path:
    """Copy ``src`` into the package under ``dest_rel`` and return the destination path."""
    dest = pkg_root / dest_rel
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dest)
    return dest


def main() -> None:
    parser = argparse.ArgumentParser(description="Build an involute gear package (.ycpkg).")
    parser.add_argument("--output", type=Path, required=True, help="Target directory for the package.")
    parser.add_argument("--teeth", type=int, default=24)
    parser.add_argument("--module-mm", type=float, default=2.0)
    parser.add_argument("--face-width-mm", type=float, default=8.0)
    parser.add_argument("--pressure-angle-deg", type=float, default=20.0)
    args = parser.parse_args()

    profile_surface, gear_solid = generate_involute_spur(
        teeth=args.teeth,
        module_mm=args.module_mm,
        face_width_mm=args.face_width_mm,
        pressure_angle_deg=args.pressure_angle_deg,
    )

    params = {
        "teeth": args.teeth,
        "module_mm": args.module_mm,
        "face_width_mm": args.face_width_mm,
        "pressure_angle_deg": args.pressure_angle_deg,
    }

    profile_meta = get_surface_metadata(profile_surface, create=True)
    set_layer(profile_meta, "gear-profile")
    add_tags(profile_meta, ["gear-profile"])
    _add_invocation(profile_meta, "INVOLUTE_SPUR2D", params)

    solid_meta = get_solid_metadata(gear_solid, create=True)
    set_layer(solid_meta, "gear")
    add_tags(solid_meta, ["gear", f"{args.teeth}-teeth"])
    _add_invocation(solid_meta, "INVOLUTE_SPUR", params)

    pkg_root = args.output
    manifest = create_package_from_entities(
        [
            {"geometry": profile_surface, "metadata": profile_meta},
            {"geometry": gear_solid, "metadata": solid_meta},
        ],
        pkg_root,
        name="Involute Spur Gear",
        version="0.1.0",
        description="Parametric involute spur gear",
        author="involute-example",
        units="mm",
        overwrite=True,
    )

    # Bundle DSL + helper scripts inside the package for reproducibility.
    example_root = Path(__file__).resolve().parent
    dsl_src = example_root / "src" / "involute_gear.dsl"
    helper_src = example_root / "scripts" / "involute_helpers.py"
    dsl_dest = _copy_into_package(dsl_src, "src/involute_gear.dsl", pkg_root)
    helper_dest = _copy_into_package(helper_src, "scripts/involute_helpers.py", pkg_root)

    source_block = manifest.data.setdefault("source", {})
    modules = source_block.setdefault("modules", [])
    modules = [entry for entry in modules if entry.get("id") != "involute_gear"]
    modules.append(
        {
            "id": "involute_gear",
            "path": "src/involute_gear.dsl",
            "language": "yapdsl-v0.1",
            "exports": ["INVOLUTE_SPUR", "INVOLUTE_SPUR2D"],
            "hash": _compute_hash(dsl_dest),
        }
    )
    source_block["modules"] = modules

    runtime = source_block.setdefault("runtime", {})
    python_runtime = runtime.setdefault("python", {})
    helpers = python_runtime.setdefault("helpers", [])
    helpers = [entry for entry in helpers if entry.get("path") != "scripts/involute_helpers.py"]
    helpers.append(
        {
            "path": "scripts/involute_helpers.py",
            "hash": _compute_hash(helper_dest),
        }
    )
    python_runtime["helpers"] = helpers
    runtime["python"] = python_runtime
    source_block["runtime"] = runtime

    manifest.save()
    print(f"Package written to {manifest.manifest_path}")


if __name__ == "__main__":
    main()
