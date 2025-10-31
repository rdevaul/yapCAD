"""Generate a parametric rocket bulkhead and package it as a .ycpkg."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable, List, Sequence

from yapcad.geom import arc, point, samplearc
from yapcad.geom3d import poly2surfaceXY, rotatesolid, solid_boolean, translatesolid
from yapcad.geom3d_util import extrude, prism
from yapcad.metadata import (
    add_design_history_entry,
    add_tags,
    ensure_solid_id,
    get_solid_metadata,
    set_layer,
    set_material,
)
from yapcad.package import create_package_from_entities

INCH_TO_MM = 25.4


def _circle_arc(radius_mm: float) -> list:
    """Return a full 360° arc in the XY plane."""
    return arc(point(0.0, 0.0), radius_mm, 0.0, 360.0)


def _arc_to_poly(arc_entity: list, segments: int) -> List[List[float]]:
    """Sample an arc into a polygon loop (closing the loop explicitly)."""
    pts = [samplearc(arc_entity, i / segments) for i in range(segments)]
    pts.append(pts[0])
    return pts


def make_bulkhead_solid(thickness_mm: float) -> list:
    disk_radius = 6.0 * INCH_TO_MM
    engine_bore_radius = 2.0 * INCH_TO_MM
    stringer_width = 1.0 * INCH_TO_MM
    stringer_height = 0.5 * INCH_TO_MM

    outer_loop = _arc_to_poly(_circle_arc(disk_radius), segments=256)
    inner_loop = _arc_to_poly(_circle_arc(engine_bore_radius), segments=192)

    base_surface, _ = poly2surfaceXY(outer_loop, holepolys=[inner_loop])
    bulkhead = extrude(base_surface, thickness_mm)

    # Prepare cutter prism (tangential width x radial height x bulkhead thickness).
    cutter_height = thickness_mm + 2.0
    cutter_core = prism(stringer_width, stringer_height, cutter_height)
    radial_center = disk_radius - stringer_height / 2.0
    cutter_core = translatesolid(
        cutter_core,
        point(0.0, radial_center, thickness_mm / 2.0),
    )

    for i in range(3):
        angle = i * 120.0
        cutter = rotatesolid(cutter_core, angle, axis=point(0.0, 0.0, 1.0))
        bulkhead = solid_boolean(bulkhead, cutter, "difference")

    meta = get_solid_metadata(bulkhead, create=True)
    ensure_solid_id(bulkhead)
    set_layer(meta, "bulkhead")
    add_tags(meta, ["bulkhead", "rocket", "thrust-structure"])
    set_material(meta, name="6061 Aluminum", grade="6061-T6", density_kg_m3=2700.0)
    add_design_history_entry(
        meta,
        author="automation",
        source="bulkhead.py",
        notes="Parametric disk bulkhead with engine bore and stringer cutouts",
    )
    return bulkhead


def create_bulkhead_package(target_dir: Path, thickness_mm: float) -> None:
    bulkhead = make_bulkhead_solid(thickness_mm)
    manifest = create_package_from_entities(
        [bulkhead],
        target_dir,
        name="rocket_bulkhead",
        version="0.6.1",
        description="Circular thrust bulkhead with engine mount bore and stringer notches",
        units="mm",
    )
    manifest.save()


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Generate a parametric rocket bulkhead package")
    parser.add_argument("--thickness-mm", type=float, default=12.0,
                        help="Bulkhead thickness in millimetres (default: 12 mm)")
    parser.add_argument("--output", type=Path, default=Path("build/rocket_bulkhead.ycpkg"),
                        help="Output .ycpkg directory")
    args = parser.parse_args()

    if args.output.exists():
        import shutil

        shutil.rmtree(args.output)

    create_bulkhead_package(args.output, args.thickness_mm)
    print(f"Bulkhead package written to {args.output}")


if __name__ == "__main__":
    main()
