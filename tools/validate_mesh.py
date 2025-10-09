#!/usr/bin/env python3
"""Run CAM-oriented linting on triangular meshes.

This script chains together three pragmatic checks that mirror the
downstream tools typically used in fabrication workflows:

1. ``admesh`` for statistics about the STL (holes, connected components,
   facet orientation, volume, bounding box, etc.).
2. ``meshfix`` for geometry healing (fills small holes, removes dangling
   triangles, produces a cleaned STL for further work).
3. A slicer front-end such as ``prusa-slicer``/``prusa-slicer-console``/
   ``slic3r``/``curaengine`` to prove the model slices successfully.

Each step is optional – if a tool is not installed the script will emit a
warning and carry on.  The goal is to provide a quick report that correlates
better with real-world CAM robustness than the strict, topological
``issolidclosed`` test.

Example
=======

.. code-block:: bash

    PYTHONPATH=src python examples/solid_boolean_demo.py \
        --mode stl --operation difference --shapes box_hole \
        --output box_hole_difference

    python tools/validate_mesh.py box_hole_difference.stl \
        --workdir build/mesh_checks

The command above keeps all intermediate files (repaired STL, generated
G-code, log snippets) under ``build/mesh_checks`` so they can be inspected
or archived alongside the STEP export.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Iterable, Optional


def _run_pymeshfix(stl: Path, output: Path) -> Optional[dict[str, object]]:
    """Attempt to repair ``stl`` using the PyMeshFix Python bindings.

    Returns a summary dictionary if successful, otherwise ``None`` when the
    dependency is missing or the mesh could not be loaded.
    """

    try:
        import pymeshfix  # type: ignore[import-untyped]
        import trimesh  # type: ignore[import-untyped]
    except ImportError:
        return None

    try:
        mesh = trimesh.load(str(stl), force='mesh')
    except Exception as exc:  # pragma: no cover - depends on trimesh internals
        return {
            "returncode": 1,
            "stdout": "",
            "stderr": f"Failed to load STL with trimesh: {exc}",
        }

    if not isinstance(mesh, trimesh.Trimesh) or mesh.vertices.size == 0:
        return {
            "returncode": 1,
            "stdout": "",
            "stderr": "Input did not contain a valid triangular mesh",
        }

    fixer = pymeshfix.MeshFix(mesh.vertices, mesh.faces)
    fixer.repair(verbose=False)
    try:
        fixer.write_stl(str(output))
    except AttributeError:
        # Older PyMeshFix releases expose the repaired data via ``v``/``f``.
        repaired = trimesh.Trimesh(fixer.v, fixer.f, process=False)
        repaired.export(str(output))
    return {
        "returncode": 0,
        "stdout": "PyMeshFix repair completed",
        "stderr": "",
        "output": str(output),
    }


ToolName = str


def which(names: Iterable[ToolName]) -> Optional[Path]:
    """Return the first executable in *names* that exists on PATH."""

    for name in names:
        path = shutil.which(name)
        if path:
            return Path(path)
    return None


def run(cmd: Iterable[str], *, cwd: Optional[Path] = None) -> subprocess.CompletedProcess:
    """Execute *cmd*, capturing stdout/stderr for later inspection."""

    result = subprocess.run(
        list(cmd),
        cwd=str(cwd) if cwd else None,
        check=False,
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
    )
    return result


def parse_admesh_report(output: str) -> dict[str, str]:
    """Extract a few useful numbers from ``admesh`` output."""

    summary: dict[str, str] = {}
    for line in output.splitlines():
        line = line.strip()
        if not line:
            continue
        for prefix in (
            "Number of facets:",
            "Volume:",
            "Number of shells:",
            "Number of pieces of geometry:",
            "Number of edges (with 2 facets):",
            "Number of edges (with 1 facet):",
        ):
            if line.startswith(prefix):
                summary[prefix.rstrip(':')] = line.split(':', 1)[1].strip()
    return summary


def slicer_command(output_dir: Path, stl: Path) -> Optional[list[str]]:
    """Return a slicer command if a supported engine is installed."""

    slicers: list[tuple[list[str], list[str]]] = [
        (["prusa-slicer-console", "prusa-slicer", "PrusaSlicer"], ["--help"]),
        (["slic3r"], ["--help"]),
        (["curaengine"], ["--help"]),
    ]

    for names, probe_args in slicers:
        binary = which(names)
        if not binary:
            continue

        # Make sure the executable is functional by running the probe command.
        probe = run([str(binary), *probe_args])
        if probe.returncode != 0:
            continue

        gcode = output_dir / (stl.stem + ".gcode")

        if binary.name.startswith("curaengine"):
            profile = os.environ.get("CURA_PROFILE")
            if not profile:
                continue
            return [str(binary), "slice", "-v", "-o", str(gcode), "-j", profile, "-l", str(stl)]

        if binary.name.startswith("slic3r"):
            return [
                str(binary),
                str(stl),
                "--output",
                str(gcode),
                "--loglevel",
                "0",
            ]

        # Default for PrusaSlicer variants.
        return [
            str(binary),
            "--export-gcode",
            str(stl),
            "--output",
            str(gcode),
        ]

    return None


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run mesh linting/healing/slicing checks in one pass.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("stl", type=Path, help="Path to the STL file to validate")
    parser.add_argument(
        "--workdir",
        type=Path,
        help="Directory to store intermediate files (default: temporary directory).",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Emit a machine-readable JSON summary instead of human text.",
    )

    args = parser.parse_args(argv)
    stl_path = args.stl

    if not stl_path.exists():
        parser.error(f"STL file not found: {stl_path}")

    if not stl_path.suffix.lower() == ".stl":
        parser.error("The input file must be an STL")

    if args.workdir:
        workdir = args.workdir
        workdir.mkdir(parents=True, exist_ok=True)
        cleanup = False
    else:
        workdir = Path(tempfile.mkdtemp(prefix="yapcad-meshcheck-"))
        cleanup = True

    summary: dict[str, object] = {
        "input": str(stl_path.resolve()),
        "workdir": str(workdir.resolve()),
        "admesh": None,
        "meshfix": None,
        "slicer": None,
        "notes": [],
    }

    # 1. admesh ------------------------------------------
    admesh_bin = which(["admesh"])
    if admesh_bin:
        repaired_stl = workdir / f"{stl_path.stem}-admesh.stl"
        cmd = [str(admesh_bin), str(stl_path)]
        result = run(cmd)
        admesh_info = {
            "returncode": result.returncode,
            "stdout": result.stdout,
            "stderr": result.stderr,
        }
        admesh_info.update(parse_admesh_report(result.stdout))
        summary["admesh"] = admesh_info
    else:
        summary["notes"].append("admesh not found on PATH; skipping mesh statistics")

    # 2. meshfix -----------------------------------------
    meshfix_bin = which(["meshfix"])
    if meshfix_bin:
        meshfix_out = workdir / f"{stl_path.stem}-meshfix.stl"
        cmd = [str(meshfix_bin), str(stl_path), "-o", str(meshfix_out)]
        result = run(cmd)
        summary["meshfix"] = {
            "returncode": result.returncode,
            "stdout": result.stdout,
            "stderr": result.stderr,
            "output": str(meshfix_out),
        }
        if result.returncode != 0:
            summary["notes"].append("meshfix reported an error; see stderr for details")
    else:
        meshfix_out = workdir / f"{stl_path.stem}-pymeshfix.stl"
        pymeshfix_info = _run_pymeshfix(stl_path, meshfix_out)
        if pymeshfix_info is not None:
            summary["meshfix"] = pymeshfix_info
            if pymeshfix_info.get("returncode") != 0:
                summary["notes"].append("PyMeshFix reported an issue; see stderr for details")
        else:
            summary["notes"].append("meshfix not found and PyMeshFix unavailable; skipping repair step")

    # 3. slicer ------------------------------------------
    slicer_cmd = slicer_command(workdir, stl_path)
    if slicer_cmd:
        result = run(slicer_cmd)
        summary["slicer"] = {
            "command": slicer_cmd,
            "returncode": result.returncode,
            "stdout": result.stdout,
            "stderr": result.stderr,
        }
        if result.returncode != 0:
            summary["notes"].append("Slicer failed – inspect stdout/stderr for diagnostics")
    else:
        summary["notes"].append(
            "No supported slicer (prusa-slicer, slic3r, curaengine) found; skipping G-code export"
        )

    if args.json:
        print(json.dumps(summary, indent=2))
    else:
        lines = [f"Mesh validation report for {stl_path}"]
        lines.append("Output directory: " + str(workdir))
        lines.append("")

        if summary["admesh"]:
            lines.append("admesh summary:")
            admesh_info = summary["admesh"]  # type: ignore[assignment]
            for key in (
                "Number of facets",
                "Number of facets connected to:",
                "Number of edges (with 1 facet)",
                "Number of edges (with 2 facets)",
                "Number of shells",
                "Volume",
            ):
                value = admesh_info.get(key)
                if value is not None:
                    lines.append(f"  {key}: {value}")
            lines.append("")
        else:
            lines.append("admesh: skipped")
            lines.append("")

        if summary["meshfix"]:
            meshfix_info = summary["meshfix"]  # type: ignore[assignment]
            lines.append("meshfix: return code %s" % meshfix_info["returncode"])
            lines.append("  output: %s" % meshfix_info["output"])
            lines.append("")
        else:
            lines.append("meshfix: skipped")
            lines.append("")

        if summary["slicer"]:
            slicer_info = summary["slicer"]  # type: ignore[assignment]
            cmd_str = " ".join(str(x) for x in slicer_info.get("command", []))
            lines.append("slicer command: " + cmd_str)
            lines.append("  return code: %s" % slicer_info["returncode"])
            lines.append("")
        else:
            lines.append("slicer: skipped")
            lines.append("")

        if summary["notes"]:
            lines.append("Notes:")
            for note in summary["notes"]:
                lines.append("  - " + note)
        else:
            lines.append("No warnings.")

        print("\n".join(lines))

    if cleanup:
        shutil.rmtree(workdir, ignore_errors=True)

    return 0


if __name__ == "__main__":
    sys.exit(main())
