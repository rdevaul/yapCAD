#!/usr/bin/env python3
"""Export geometry from a yapCAD `.ycpkg` package to interchange formats."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Iterable, List, Sequence

from yapcad.geom import (
    isgeomlist,
    isline,
    isarc,
    iscircle,
    ispoint,
    iscatmullrom,
    isnurbs,
    sample,
)
from yapcad.geom3d import issolid, issurface
from yapcad.geom_util import geomlist2poly_components, geomlist2poly_with_holes
from yapcad.io import write_step, write_stl
from yapcad.io.geometry_json import geometry_from_json, SCHEMA_ID
from yapcad.package import PackageManifest, load_geometry


def _export_step(solids: Iterable[Sequence], out_dir: Path) -> List[Path]:
    paths: List[Path] = []
    for idx, solid in enumerate(solids, 1):
        target = out_dir / f"solid_{idx:02d}.step"
        write_step(solid, target)
        paths.append(target)
    return paths


def _export_stl(solids: Iterable[Sequence], out_dir: Path) -> List[Path]:
    paths: List[Path] = []
    for idx, solid in enumerate(solids, 1):
        target = out_dir / f"solid_{idx:02d}.stl"
        write_stl(solid, target)
        paths.append(target)
    return paths


def _export_surfaces_step(surfaces: Iterable[Sequence], out_dir: Path) -> List[Path]:
    paths: List[Path] = []
    for idx, surface in enumerate(surfaces, 1):
        target = out_dir / f"surface_{idx:02d}.step"
        write_step(surface, target)
        paths.append(target)
    return paths


def _export_sketches_dxf(sketches: Iterable[Sequence], out_path: Path) -> Path:
    import ezdxf  # imported lazily to avoid dependency when DXF export not requested

    doc = ezdxf.new("R2010")
    msp = doc.modelspace()

    def _as_xy(pt: Sequence[float]) -> tuple[float, float]:
        return float(pt[0]), float(pt[1])

    def _close_pts(a: Sequence[float], b: Sequence[float], tol: float = 1e-4) -> bool:
        dx = float(a[0]) - float(b[0])
        dy = float(a[1]) - float(b[1])
        return (dx * dx + dy * dy) ** 0.5 <= tol

    def _handle_element(element: Sequence, layer: str) -> int:
        if isline(element):
            start = _as_xy(element[0])
            end = _as_xy(element[1])
            msp.add_line(start, end, dxfattribs={"layer": layer})
            return 1
        if iscircle(element):
            center = _as_xy(element[0])
            radius = float(element[1][0])
            msp.add_circle(center, radius, dxfattribs={"layer": layer})
            return 1
        if isarc(element):
            center = _as_xy(element[0])
            radius = float(element[1][0])
            start_angle_raw = float(element[1][1])
            end_angle_raw = float(element[1][2])
            start_angle = start_angle_raw % 360.0
            end_angle = end_angle_raw % 360.0
            if abs(end_angle_raw - start_angle_raw) >= 360.0:
                msp.add_circle(center, radius, dxfattribs={"layer": layer})
                return 1
            orient = element[1][3]
            if orient == -2:
                start_angle, end_angle = end_angle, start_angle
            msp.add_arc(
                center,
                radius,
                start_angle=start_angle,
                end_angle=end_angle,
                dxfattribs={"layer": layer},
            )
            return 1
        if iscatmullrom(element) or isnurbs(element):
            segs = 64
            coords = [_as_xy(sample(element, i / segs)) for i in range(segs + 1)]
            closed = False
            if iscatmullrom(element):
                closed = bool(element[2].get("closed", False))
            elif isnurbs(element):
                ctrl = element[1]
                if ctrl:
                    closed = _close_pts(ctrl[0], ctrl[-1])
            if len(coords) >= 2:
                close = closed or _close_pts(coords[0], coords[-1])
                if close:
                    coords[-1] = coords[0]
                    coords = coords[:-1]
                msp.add_lwpolyline(coords, format="xy", close=close, dxfattribs={"layer": layer})
                return 1
            return 0
        if isinstance(element, list) and element and all(isinstance(pt, list) and len(pt) >= 2 for pt in element):
            coords = [_as_xy(pt) for pt in element]
            if len(coords) < 2:
                return 0
            close = _close_pts(coords[0], coords[-1])
            if close:
                coords = coords[:-1]
            if len(coords) < 2:
                return 0
            msp.add_lwpolyline(coords, format="xy", close=close, dxfattribs={"layer": layer})
            return 1
        return 0

    for idx, geom in enumerate(sketches, 1):
        layer_name = f"SKETCH_{idx:02d}"
        if layer_name not in doc.layers:
            doc.layers.add(name=layer_name)
        written = 0
        for element in geom:
            written += _handle_element(element, layer_name)
        if written == 0:
            raise ValueError("failed to derive DXF entities from sketch geometry")

    doc.saveas(out_path)
    return out_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Export geometry from a yapCAD .ycpkg package.",
    )
    parser.add_argument("package", type=Path, help="Path to the .ycpkg directory.")
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Destination directory for exported files (default: PACKAGE/exports).",
    )
    parser.add_argument(
        "--format",
        choices=["step", "stl", "dxf", "all"],
        action="append",
        help="Export format(s) to generate. Repeat this flag for multiple outputs.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Allow replacing existing files in the output directory.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    manifest = PackageManifest.load(args.package)
    geometry = load_geometry(manifest)

    solids = []
    surfaces = []
    sketches: List[Sequence] = []

    for entity in geometry:
        if issolid(entity):
            solids.append(entity)
        elif issurface(entity):
            surfaces.append(entity)
    doc_path = manifest.geometry_primary_path()
    with doc_path.open("r", encoding="utf-8") as fh:
        doc = json.load(fh)
    sketch_entities = [entry for entry in doc.get("entities", []) if entry.get("type") == "sketch"]
    for entry in sketch_entities:
        subdoc = {"schema": SCHEMA_ID, "entities": [entry]}
        sketches.extend(geometry_from_json(subdoc))
    formats = args.format or []
    requested = set()
    for fmt in formats:
        if fmt == "all":
            requested.update({"step", "stl", "dxf"})
        else:
            requested.add(fmt)
    if not requested:
        requested = {"step"}

    out_dir = args.output or (manifest.root / "exports")
    out_dir.mkdir(parents=True, exist_ok=True)

    written: List[Path] = []

    def _check_target(path: Path) -> None:
        if path.exists() and not args.overwrite:
            raise FileExistsError(f"export target already exists: {path}")

    if solids:
        if "step" in requested:
            targets = [out_dir / f"solid_{idx:02d}.step" for idx in range(1, len(solids) + 1)]
            for target in targets:
                _check_target(target)
            written.extend(_export_step(solids, out_dir))
        if "stl" in requested:
            targets = [out_dir / f"solid_{idx:02d}.stl" for idx in range(1, len(solids) + 1)]
            for target in targets:
                _check_target(target)
            written.extend(_export_stl(solids, out_dir))
    if surfaces and "step" in requested:
        targets = [out_dir / f"surface_{idx:02d}.step" for idx in range(1, len(surfaces) + 1)]
        for target in targets:
            _check_target(target)
        written.extend(_export_surfaces_step(surfaces, out_dir))
    if sketches and "dxf" in requested:
        dxf_path = out_dir / "sketches.dxf"
        _check_target(dxf_path)
        written.append(_export_sketches_dxf(sketches, dxf_path))

    if not written:
        print("No geometry exported (check formats and package contents).", file=sys.stderr)
        return 1

    for path in written:
        print(path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
