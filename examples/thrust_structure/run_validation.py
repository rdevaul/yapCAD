#!/usr/bin/env python3
"""
Run FEA validation on the thrust structure.

This script:
1. Creates a yapCAD package with the optimized thrust plate geometry
2. Runs the FEA analysis using the thrust-fea.yaml plan
3. Outputs results including displacement and stress fields

Prerequisites:
    conda install -c conda-forge fenics-dolfinx gmsh meshio pyvista

Usage:
    python run_validation.py [--design baseline|light|optimized]

Outputs:
    thrust_structure.ycpkg/
    └── validation/
        └── results/
            └── thrust-fea/
                ├── summary.json      # Analysis results
                ├── mesh.msh          # Gmsh mesh
                ├── mesh.xdmf         # XDMF mesh for DOLFINx
                ├── displacement.vtu  # Displacement field (ParaView)
                └── stress.vtu        # Von Mises stress field (ParaView)

Visualization:
    paraview thrust_structure.ycpkg/validation/results/thrust-fea/stress.vtu
"""

import argparse
import json
import shutil
import sys
from pathlib import Path

# Add src to path for development
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))


DESIGNS = {
    "baseline": "MAKE_BASELINE_PLATE",
    "light": "MAKE_LIGHTENED_PLATE_3",
    "optimized": "MAKE_OPTIMIZED_PLATE",
}


def check_dependencies():
    """Check if required dependencies are available."""
    missing = []

    try:
        from yapcad.package.analysis.gmsh_mesher import gmsh_available
        if not gmsh_available():
            missing.append("gmsh")
    except ImportError:
        missing.append("gmsh")

    try:
        from yapcad.package.analysis.fenics import fenics_available
        if not fenics_available():
            missing.append("fenics-dolfinx")
    except ImportError:
        missing.append("fenics-dolfinx")

    if missing:
        print("Missing dependencies:", ", ".join(missing))
        print("\nInstall with:")
        print("  conda install -c conda-forge fenics-dolfinx gmsh meshio")
        return False

    return True


def create_package(design: str, output_path: Path) -> bool:
    """Create a yapCAD package with the thrust plate geometry."""
    from yapcad.dsl import compile_and_run
    from yapcad.package import create_package_from_entities

    print(f"\n1. Generating {design} geometry...")

    # Read DSL source
    dsl_path = Path(__file__).parent / "thrust_structure.dsl"
    source = dsl_path.read_text()

    # Generate geometry
    command = DESIGNS[design]
    result = compile_and_run(source, command, {})

    if not result.success:
        print(f"   Error: {result.error_message}")
        return False

    solid = result.geometry
    print(f"   Generated solid with volume {result.volume:.0f} mm³")

    # Create package
    print(f"\n2. Creating package at {output_path}...")
    if output_path.exists():
        shutil.rmtree(output_path)

    create_package_from_entities(
        [solid],
        output_path,
        name=f"Thrust Structure ({design})",
        version="1.0.0",
        description="Liquid rocket motor thrust plate for FEA validation",
        units="mm",
        materials={
            "6061-T6": {
                "name": "6061-T6 Aluminum",
                "youngs_modulus_pa": 68.9e9,
                "poisson_ratio": 0.33,
                "density_kgm3": 2700,
                "yield_strength_pa": 276e6,
            }
        },
    )

    # Copy validation plans
    plans_src = Path(__file__).parent / "validation" / "plans"
    plans_dst = output_path / "validation" / "plans"
    plans_dst.mkdir(parents=True, exist_ok=True)
    for plan_file in plans_src.glob("*.yaml"):
        shutil.copy(plan_file, plans_dst)
        print(f"   Copied {plan_file.name}")

    return True


def run_fea(package_path: Path) -> dict:
    """Run FEA analysis on the package."""
    from yapcad.package import PackageManifest
    from yapcad.package.analysis import load_analysis_plan
    from yapcad.package.analysis.fenics import FenicsxAdapter

    print("\n3. Running FEA analysis...")

    manifest = PackageManifest.load(package_path)
    plan_path = package_path / "validation" / "plans" / "thrust-fea.yaml"
    plan = load_analysis_plan(plan_path)

    workspace = package_path / "validation" / "results" / "thrust-fea"
    workspace.mkdir(parents=True, exist_ok=True)

    adapter = FenicsxAdapter()
    result = adapter.run(manifest, plan, workspace)

    return {
        "status": result.status,
        "metrics": result.metrics,
        "artifacts": result.artifacts,
    }


def print_results(results: dict):
    """Print analysis results."""
    print("\n" + "=" * 60)
    print("ANALYSIS RESULTS")
    print("=" * 60)

    status = results["status"]
    metrics = results["metrics"]

    status_str = "✓ PASSED" if status == "passed" else "✗ FAILED" if status == "failed" else f"⚠ {status.upper()}"
    print(f"\nStatus: {status_str}")

    print("\nMetrics:")
    if "displacement.max_mm" in metrics:
        print(f"  Max displacement:  {metrics['displacement.max_mm']:.3f} mm")
    if "stress.von_mises.max_mpa" in metrics:
        print(f"  Max von Mises:     {metrics['stress.von_mises.max_mpa']:.1f} MPa")

    print("\nArtifacts:")
    for artifact in results.get("artifacts", []):
        desc = artifact.get("description", artifact.get("kind", ""))
        print(f"  - {artifact['path']}: {desc}")

    print("\nVisualization:")
    print("  paraview validation/results/thrust-fea/displacement.vtu")
    print("  paraview validation/results/thrust-fea/stress.vtu")


def main():
    parser = argparse.ArgumentParser(
        description="Run FEA validation on thrust structure"
    )
    parser.add_argument(
        "--design", choices=["baseline", "light", "optimized"],
        default="optimized",
        help="Design variant to analyze (default: optimized)"
    )
    parser.add_argument(
        "--output", type=Path, default=Path("thrust_structure.ycpkg"),
        help="Output package path"
    )
    parser.add_argument(
        "--skip-fea", action="store_true",
        help="Skip FEA (just create package)"
    )

    args = parser.parse_args()

    # Check dependencies
    if not args.skip_fea and not check_dependencies():
        print("\nUse --skip-fea to create package without running FEA")
        return 1

    # Create package
    if not create_package(args.design, args.output):
        return 1

    print(f"\n   Package created: {args.output}")

    # Run FEA
    if not args.skip_fea:
        try:
            results = run_fea(args.output)
            print_results(results)

            # Save summary
            summary_path = args.output / "validation" / "results" / "thrust-fea" / "summary.json"
            with open(summary_path, "w") as f:
                json.dump(results, f, indent=2)

        except Exception as e:
            print(f"\nFEA Error: {e}")
            import traceback
            traceback.print_exc()
            return 1
    else:
        print("\n   Skipping FEA (--skip-fea)")
        print(f"\n   To run FEA later:")
        print(f"   python run_validation.py --design {args.design}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
