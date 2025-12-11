#!/usr/bin/env python3
"""
Optimization driver for thrust structure design.

This script compares different hole configurations to find
the minimum-mass design that meets deflection requirements.

Uses actual FEA (FEniCSx) for deflection calculations when available,
falling back to analytical approximation otherwise.

Usage:
    python optimize.py [--output DIR] [--use-fea] [--skip-fea]
"""

import argparse
import json
import math
import shutil
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# Add src to path for development
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))


@dataclass
class DesignConfig:
    """Configuration for a single design evaluation."""
    name: str
    command: str
    num_holes_total: int
    hole_radius_mm: float
    hole_radial_fraction: float

    def config_id(self) -> str:
        return self.name.lower().replace(" ", "_").replace("(", "").replace(")", "").replace(",", "")


@dataclass
class EvaluationResult:
    """Result of evaluating a single design."""
    config: DesignConfig
    mass_kg: float
    volume_mm3: float
    max_displacement_mm: float
    max_stress_mpa: float
    passed: bool
    method: str  # "fea" or "analytical"
    notes: str = ""


def get_configurations() -> List[DesignConfig]:
    """Get the available design configurations."""
    return [
        DesignConfig(
            name="Baseline (no holes)",
            command="MAKE_BASELINE_PLATE",
            num_holes_total=0,
            hole_radius_mm=0.0,
            hole_radial_fraction=0.0,
        ),
        DesignConfig(
            name="Light (9 holes, 20mm)",
            command="MAKE_LIGHTENED_PLATE_3",
            num_holes_total=9,
            hole_radius_mm=20.0,
            hole_radial_fraction=0.6,
        ),
        DesignConfig(
            name="Medium (12 holes, 20mm)",
            command="MAKE_LIGHTENED_PLATE_4",
            num_holes_total=12,
            hole_radius_mm=20.0,
            hole_radial_fraction=0.6,
        ),
        DesignConfig(
            name="Optimized (12 holes, 25mm)",
            command="MAKE_OPTIMIZED_PLATE",
            num_holes_total=12,
            hole_radius_mm=25.0,
            hole_radial_fraction=0.65,
        ),
    ]


def check_fea_available() -> bool:
    """Check if FEA dependencies are available."""
    try:
        from yapcad.package.analysis.gmsh_mesher import gmsh_available
        from yapcad.package.analysis.fenics import fenics_available
        return gmsh_available() and fenics_available()
    except ImportError:
        return False


def generate_geometry(config: DesignConfig) -> Tuple[Any, float]:
    """Generate geometry for a configuration and return (solid, volume_mm3)."""
    from yapcad.dsl import compile_and_run

    # Read DSL source
    dsl_path = Path(__file__).parent / "thrust_structure.dsl"
    source = dsl_path.read_text()

    # Run DSL to generate geometry (no parameters needed, using defaults)
    result = compile_and_run(source, config.command, {})

    if not result.success:
        raise RuntimeError(f"DSL execution failed: {result.error_message}")

    solid = result.geometry

    # Get volume from the solid
    from yapcad.geom3d import volumeof
    volume = volumeof(solid)

    return solid, volume


def run_fea_analysis(
    solid: Any,
    config: DesignConfig,
    work_dir: Path
) -> Tuple[float, float]:
    """
    Run actual FEA analysis on the solid.

    Returns (max_displacement_mm, max_stress_mpa)
    """
    from yapcad.package import create_package_from_entities
    from yapcad.package.core import PackageManifest
    from yapcad.package.analysis import load_plan
    from yapcad.package.analysis.fenics import FenicsxAdapter

    # Create a temporary package for this design
    pkg_path = work_dir / f"{config.config_id()}.ycpkg"
    if pkg_path.exists():
        shutil.rmtree(pkg_path)

    create_package_from_entities(
        [solid],
        pkg_path,
        name=config.name,
        version="1.0.0",
        units="mm"
    )

    # Copy the analysis plan
    plan_src = Path(__file__).parent / "validation" / "plans" / "thrust-fea.yaml"
    plan_dst = pkg_path / "validation" / "plans"
    plan_dst.mkdir(parents=True, exist_ok=True)
    shutil.copy(plan_src, plan_dst / "thrust-fea.yaml")

    # Load manifest and plan
    manifest = PackageManifest.load(pkg_path)
    plan = load_plan(plan_dst / "thrust-fea.yaml")

    # Create workspace for results
    workspace = pkg_path / "validation" / "results" / "thrust-fea"
    workspace.mkdir(parents=True, exist_ok=True)

    # Run FEA
    adapter = FenicsxAdapter()
    result = adapter.run(manifest, plan, workspace)

    if result.status == "error":
        error_msg = result.summary.get("error", "Unknown error")
        raise RuntimeError(f"FEA failed: {error_msg}")

    # Extract results
    max_disp = result.metrics.get("displacement.max_mm", 0.0)
    max_stress = result.metrics.get("stress.von_mises.max_mpa", 0.0)

    return max_disp, max_stress


def analytical_deflection(
    volume_mm3: float,
    num_holes: int,
    hole_radius_mm: float
) -> Tuple[float, float]:
    """
    Analytical approximation for deflection of circular plate.

    Returns (max_displacement_mm, max_stress_mpa)
    """
    # Material properties
    E = 68.9e9  # Pa (6061-T6)
    nu = 0.33

    # Geometry (convert to meters)
    outer_diameter_mm = 304.8
    motor_mount_diameter_mm = 101.6
    thickness_mm = 8.0

    a = outer_diameter_mm / 2 / 1000  # outer radius in m
    b = motor_mount_diameter_mm / 2 / 1000  # inner radius in m
    t = thickness_mm / 1000  # thickness in m

    # Load
    F = 6672  # N (1500 lbf)
    q = F / (2 * math.pi * b)  # N/m line load at inner edge

    # Plate stiffness
    D = E * t**3 / (12 * (1 - nu**2))

    # Account for holes reducing stiffness
    if num_holes > 0:
        hole_area = num_holes * math.pi * (hole_radius_mm/1000)**2
        plate_area = math.pi * (a**2 - b**2)
        stiffness_factor = max(0.3, 1 - 1.5 * hole_area / plate_area)
    else:
        stiffness_factor = 1.0

    D_eff = D * stiffness_factor

    # Deflection calculation (simplified annular plate)
    beta = b / a
    deflection_factor = (1 - beta**2)**2 / 16
    w_max = q * b * deflection_factor * a**2 / D_eff

    # Convert to mm
    max_disp_mm = abs(w_max) * 1000

    # Estimate stress (bending stress at inner edge)
    M_max = q * b * (1 - beta**2) / 4
    sigma_max = 6 * M_max / t**2
    max_stress_mpa = sigma_max / 1e6

    return max_disp_mm, max_stress_mpa


def evaluate_config(
    config: DesignConfig,
    deflection_limit_mm: float = 10.0,
    use_fea: bool = True,
    work_dir: Optional[Path] = None
) -> EvaluationResult:
    """Evaluate a single configuration."""
    print(f"\n  Evaluating: {config.name}")

    try:
        # Generate geometry
        solid, volume_mm3 = generate_geometry(config)

        # Calculate mass (density = 2700 kg/m³)
        volume_m3 = volume_mm3 / 1e9
        mass_kg = volume_m3 * 2700

        print(f"    Volume: {volume_mm3:.0f} mm³")
        print(f"    Mass: {mass_kg:.3f} kg ({mass_kg*2.205:.3f} lbs)")

        # Get deflection - FEA or analytical
        if use_fea and work_dir is not None:
            try:
                print(f"    Running FEA analysis...")
                max_disp, max_stress = run_fea_analysis(solid, config, work_dir)
                method = "fea"
            except Exception as e:
                print(f"    FEA failed ({e}), falling back to analytical")
                max_disp, max_stress = analytical_deflection(
                    volume_mm3, config.num_holes_total, config.hole_radius_mm
                )
                method = "analytical"
        else:
            max_disp, max_stress = analytical_deflection(
                volume_mm3, config.num_holes_total, config.hole_radius_mm
            )
            method = "analytical"

        passed = max_disp <= deflection_limit_mm
        status = "PASS" if passed else "FAIL"
        print(f"    Deflection: {max_disp:.3f} mm [{status}] ({method})")
        print(f"    Max Stress: {max_stress:.1f} MPa")

        return EvaluationResult(
            config=config,
            mass_kg=mass_kg,
            volume_mm3=volume_mm3,
            max_displacement_mm=max_disp,
            max_stress_mpa=max_stress,
            passed=passed,
            method=method,
        )

    except Exception as e:
        print(f"    Error: {e}")
        import traceback
        traceback.print_exc()
        return EvaluationResult(
            config=config,
            mass_kg=0.0,
            volume_mm3=0.0,
            max_displacement_mm=float('inf'),
            max_stress_mpa=0.0,
            passed=False,
            method="error",
            notes=str(e),
        )


def run_optimization(
    output_dir: Path,
    deflection_limit_mm: float = 10.0,
    use_fea: bool = True
) -> Dict[str, Any]:
    """Run the optimization comparison."""
    print("=" * 60)
    print("Thrust Structure Design Optimization")
    print("=" * 60)
    print(f"\nDeflection limit: {deflection_limit_mm} mm")
    print(f"Material: 6061-T6 Aluminum")
    print(f"Thrust load: 1500 lbf (6672 N)")

    # Check FEA availability
    fea_available = check_fea_available()
    if use_fea and not fea_available:
        print("\nWARNING: FEA not available (missing gmsh or dolfinx)")
        print("         Using analytical approximation instead")
        use_fea = False
    elif use_fea:
        print(f"\nUsing FEA (FEniCSx) for deflection analysis")
    else:
        print(f"\nUsing analytical approximation for deflection")

    configs = get_configurations()
    results: List[EvaluationResult] = []

    # Create work directory for FEA packages
    work_dir = output_dir / "fea_work"
    work_dir.mkdir(parents=True, exist_ok=True)

    for config in configs:
        result = evaluate_config(
            config,
            deflection_limit_mm,
            use_fea=use_fea,
            work_dir=work_dir
        )
        results.append(result)

    # Find best passing design
    passing = [r for r in results if r.passed]
    if passing:
        best = min(passing, key=lambda r: r.mass_kg)
    else:
        best = None

    # Baseline for comparison
    baseline = results[0]

    print("\n" + "=" * 60)
    print("RESULTS SUMMARY")
    print("=" * 60)

    print(f"\n{'Design':<28} {'Mass (kg)':<10} {'Deflection':<12} {'Stress':<10} {'Status'}")
    print("-" * 75)
    for r in results:
        status = "PASS" if r.passed else "FAIL"
        method_tag = f"[{r.method[:3]}]" if r.method else ""
        print(f"{r.config.name:<28} {r.mass_kg:>7.3f}    {r.max_displacement_mm:>7.3f} mm   {r.max_stress_mpa:>6.1f} MPa  {status} {method_tag}")

    mass_reduction = 0.0
    if best:
        mass_reduction = (baseline.mass_kg - best.mass_kg) / baseline.mass_kg * 100
        print(f"\nBest design: {best.config.name}")
        print(f"Mass reduction: {mass_reduction:.1f}% vs baseline")
    else:
        print("\nNo design met the deflection requirement!")

    # Prepare results dict
    result_data = {
        "timestamp": datetime.now().isoformat(),
        "deflection_limit_mm": deflection_limit_mm,
        "analysis_method": "fea" if use_fea else "analytical",
        "fea_available": fea_available,
        "designs": [
            {
                "name": r.config.name,
                "command": r.config.command,
                "num_holes": r.config.num_holes_total,
                "hole_radius_mm": r.config.hole_radius_mm,
                "mass_kg": r.mass_kg,
                "volume_mm3": r.volume_mm3,
                "max_displacement_mm": r.max_displacement_mm,
                "max_stress_mpa": r.max_stress_mpa,
                "passed": r.passed,
                "method": r.method,
            }
            for r in results
        ],
        "best_design": best.config.name if best else None,
        "mass_reduction_pct": mass_reduction,
    }

    return result_data


def export_best_design(results: Dict[str, Any], output_dir: Path) -> None:
    """Export the best design to STEP."""
    best_name = results.get("best_design")
    if not best_name:
        print("\nNo best design to export")
        return

    # Find the command for the best design
    for design in results["designs"]:
        if design["name"] == best_name:
            command = design["command"]
            break
    else:
        return

    print(f"\nExporting best design: {best_name}")

    from yapcad.dsl import compile_and_run

    dsl_path = Path(__file__).parent / "thrust_structure.dsl"
    source = dsl_path.read_text()

    result = compile_and_run(source, command, {})

    if result.success and result.geometry:
        # Export STEP
        step_path = output_dir / "best_thrust_plate.step"
        try:
            from yapcad.io import write_step
            write_step(result.geometry, str(step_path))
            print(f"  STEP: {step_path}")
        except Exception as e:
            print(f"  STEP export failed: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Compare thrust structure designs for minimum mass using FEA"
    )
    parser.add_argument(
        "--output", type=Path, default=Path("./output/optimization"),
        help="Output directory for results"
    )
    parser.add_argument(
        "--deflection-limit", type=float, default=10.0,
        help="Maximum allowable deflection in mm (default: 10)"
    )
    parser.add_argument(
        "--use-fea", action="store_true", default=True,
        help="Use FEA for deflection analysis (default: True)"
    )
    parser.add_argument(
        "--skip-fea", action="store_true",
        help="Skip FEA, use analytical approximation only"
    )

    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    use_fea = args.use_fea and not args.skip_fea

    # Run comparison
    results = run_optimization(args.output, args.deflection_limit, use_fea=use_fea)

    # Save results
    results_path = args.output / "comparison_results.json"
    with open(results_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to: {results_path}")

    # Export best design
    export_best_design(results, args.output)

    return 0


if __name__ == "__main__":
    sys.exit(main())
