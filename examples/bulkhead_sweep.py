"""Parameter sweep for the rocket bulkhead thickness.

Generates a `.ycpkg` for each thickness value, writes a CalculiX plan, runs
the analysis, and prints the maximum axial deflection reported by `ccx`.

Usage::

    python examples/bulkhead_sweep.py --start 5 --stop 12 --step 1

Requires the CalculiX executable (`ccx`) to be on the PATH.
"""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path
from typing import Iterable, List, Tuple

from bulkhead import create_bulkhead_package
from yapcad.package.analysis.cli import analyze_package

INCH_TO_MM = 25.4
DISK_RADIUS_MM = 6.0 * INCH_TO_MM
ENGINE_BORE_RADIUS_MM = 2.0 * INCH_TO_MM

PLAN_TEMPLATE = """id: bulkhead-fea
kind: fea
backend: calculix
name: Bulkhead thickness sweep
geometry:
  inner_radius_mm: {inner_radius}
  outer_radius_mm: {outer_radius}
backendOptions:
  inner_radius_mm: {inner_radius}
  outer_radius_mm: {outer_radius}
  thickness_mm: {thickness}
  thrust_n: 2224.0
  radial_divisions: 32
  thickness_divisions: 3
  youngs_modulus_mpa: 68900.0
  poisson_ratio: 0.33
  density_tonemm3: 2.70e-6
loads:
  - id: engine-thrust
    type: axial
    magnitude_n: 2224.0
    application: inner-bore
boundaryConditions:
  - id: stringer-constraint
    type: fixed
    surfaces: [stringer]
acceptance:
  displacement.max:
    limit_mm: 10.0
    comparison: "<="
execution:
  mode: local
  command: ccx
  env:
    OMP_NUM_THREADS: "8"
"""


def _thickness_range(start: float, stop: float, step: float) -> Iterable[float]:
    value = start
    while value <= stop + 1e-9:
        yield round(value, 6)
        value += step


def _write_plan(package_root: Path, thickness_mm: float) -> Path:
    plan_dir = package_root / "validation" / "plans"
    plan_dir.mkdir(parents=True, exist_ok=True)
    plan_path = plan_dir / "bulkhead_fea.yaml"
    plan_text = PLAN_TEMPLATE.format(
        inner_radius=f"{ENGINE_BORE_RADIUS_MM:.6f}",
        outer_radius=f"{DISK_RADIUS_MM:.6f}",
        thickness=f"{thickness_mm:.6f}",
    )
    plan_path.write_text(plan_text, encoding="utf-8")
    return plan_path.relative_to(package_root)


def sweep_thicknesses(
    start_mm: float,
    stop_mm: float,
    step_mm: float,
    output_dir: Path,
) -> List[Tuple[float, float | None, str]]:
    output_dir.mkdir(parents=True, exist_ok=True)

    results: List[Tuple[float, float | None, str]] = []
    for thickness in _thickness_range(start_mm, stop_mm, step_mm):
        package_dir = output_dir / f"bulkhead_{thickness:.1f}mm.ycpkg"
        if package_dir.exists():
            shutil.rmtree(package_dir)

        create_bulkhead_package(package_dir, thickness)
        plan_rel = _write_plan(package_dir, thickness)

        results_root = package_dir / "validation" / "results" / "bulkhead-fea"
        if results_root.exists():
            shutil.rmtree(results_root)

        summary_path = analyze_package(package_dir, plan_rel)
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        disp_mm = summary.get("metrics", {}).get("displacement.max_mm")
        status = summary.get("status", "pending")
        results.append((thickness, disp_mm, status))

    return results


def _print_results(results: List[Tuple[float, float | None, str]]) -> None:
    header = f"{'Thickness (mm)':>15}  {'Max Deflection (mm)':>22}  {'Status':>10}"
    print(header)
    print("-" * len(header))
    for thickness, disp, status in results:
        disp_text = f"{disp:.6f}" if disp is not None else "--"
        print(f"{thickness:15.3f}  {disp_text:>22}  {status:>10}")

    passed = [t for t, _, status in results if status == "passed"]
    if passed:
        print(f"\nMinimum passing thickness: {min(passed):.3f} mm")


def main() -> None:
    parser = argparse.ArgumentParser(description="Sweep bulkhead thickness and run CalculiX analysis")
    parser.add_argument("--start", type=float, default=5.0, help="Starting thickness in mm (default: 5)")
    parser.add_argument("--stop", type=float, default=12.0, help="Ending thickness in mm (default: 12)")
    parser.add_argument("--step", type=float, default=1.0, help="Thickness increment in mm (default: 1)")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("build/bulkhead_sweep"),
        help="Directory to store generated packages (default: build/bulkhead_sweep)",
    )
    args = parser.parse_args()

    results = sweep_thicknesses(args.start, args.stop, args.step, args.output_dir)
    _print_results(results)


if __name__ == "__main__":
    main()
