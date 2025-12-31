#!/usr/bin/env python3
"""
Run validation checks on the thrust structure.

This script:
1. Creates a yapCAD package with the optimized thrust plate geometry
2. Runs the mass budget check using mass-check.yaml (no external deps)
3. Runs the FEA analysis using thrust-fea.yaml (requires fenics)
4. Outputs results including displacement and stress fields

Prerequisites:
    For FEA analysis only:
    conda install -c conda-forge fenics-dolfinx gmsh meshio pyvista

Usage:
    python run_validation.py [--design baseline|light|optimized]
    python run_validation.py --skip-fea   # Run only mass check
    python run_validation.py --mass-only  # Run only mass check (alias)

Outputs:
    thrust_structure.ycpkg/
    └── validation/
        └── results/
            ├── mass-check/
            │   └── summary.json      # Mass budget results
            └── thrust-fea/
                ├── summary.json      # FEA analysis results
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

    # Direct import test first
    try:
        import gmsh
        gmsh_direct = True
    except ImportError as e:
        gmsh_direct = False
        print(f"DEBUG: Direct gmsh import failed: {e}")

    try:
        import dolfinx
        dolfinx_direct = True
    except ImportError as e:
        dolfinx_direct = False
        print(f"DEBUG: Direct dolfinx import failed: {e}")

    # Now test via yapCAD wrappers
    try:
        from yapcad.package.analysis.gmsh_mesher import gmsh_available
        gmsh_wrapper = gmsh_available()
        if not gmsh_wrapper:
            missing.append("gmsh")
            print(f"DEBUG: gmsh_available() returned False (direct import: {gmsh_direct})")
    except ImportError as e:
        missing.append("gmsh")
        print(f"DEBUG: Failed to import gmsh_mesher: {e}")

    try:
        from yapcad.package.analysis.fenics import fenics_available
        fenics_wrapper = fenics_available()
        if not fenics_wrapper:
            missing.append("fenics-dolfinx")
            print(f"DEBUG: fenics_available() returned False (direct import: {dolfinx_direct})")
    except ImportError as e:
        missing.append("fenics-dolfinx")
        print(f"DEBUG: Failed to import fenics: {e}")

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

    # Calculate volume
    from yapcad.geom3d import volumeof
    volume = volumeof(solid)
    print(f"   Generated solid with volume {volume:.0f} mm³")

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


def run_mass_check(package_path: Path) -> dict:
    """Run mass budget check on the package (no external deps required)."""
    from yapcad.package import PackageManifest
    from yapcad.package.analysis import load_plan
    from yapcad.package.analysis.yapcad_native import YapCADNativeAdapter

    print("\n3. Running mass budget check...")

    manifest = PackageManifest.load(package_path)
    plan_path = package_path / "validation" / "plans" / "mass-check.yaml"
    plan = load_plan(plan_path)

    workspace = package_path / "validation" / "results" / "mass-check"
    workspace.mkdir(parents=True, exist_ok=True)

    adapter = YapCADNativeAdapter()
    result = adapter.run(manifest, plan, workspace)

    # Check for error details
    if result.status == "error":
        error_msg = result.summary.get("error", "Unknown error")
        status_detail = result.summary.get("statusDetail", "")
        print(f"\n   Mass check error: {error_msg}")
        if status_detail:
            print(f"   Details: {status_detail}")

    return {
        "status": result.status,
        "metrics": result.metrics,
        "artifacts": result.artifacts,
        "error": result.summary.get("error"),
    }


def run_fea(package_path: Path) -> dict:
    """Run FEA analysis on the package."""
    from yapcad.package import PackageManifest
    from yapcad.package.analysis import load_plan
    from yapcad.package.analysis.fenics import FenicsxAdapter

    print("\n4. Running FEA analysis...")

    manifest = PackageManifest.load(package_path)
    plan_path = package_path / "validation" / "plans" / "thrust-fea.yaml"
    plan = load_plan(plan_path)

    workspace = package_path / "validation" / "results" / "thrust-fea"
    workspace.mkdir(parents=True, exist_ok=True)

    adapter = FenicsxAdapter()
    result = adapter.run(manifest, plan, workspace)

    # Check for error details
    if result.status == "error":
        error_msg = result.summary.get("error", "Unknown error")
        status_detail = result.summary.get("statusDetail", "")
        print(f"\n   Analysis error: {error_msg}")
        if status_detail:
            print(f"   Details: {status_detail}")

    return {
        "status": result.status,
        "metrics": result.metrics,
        "artifacts": result.artifacts,
        "error": result.summary.get("error"),
    }


def print_mass_results(results: dict):
    """Print mass check results."""
    print("\n" + "=" * 60)
    print("MASS BUDGET CHECK")
    print("=" * 60)

    status = results["status"]
    metrics = results["metrics"]

    status_str = "PASSED" if status == "passed" else "FAILED" if status == "failed" else status.upper()
    print(f"\nStatus: {status_str}")

    print("\nMetrics:")
    if "mass_kg" in metrics:
        print(f"  Mass:     {metrics['mass_kg']:.3f} kg ({metrics.get('mass_g', 0):.1f} g)")
    if "volume_mm3" in metrics:
        print(f"  Volume:   {metrics['volume_mm3']:.0f} mm^3")
    if "density_kgm3" in metrics:
        print(f"  Density:  {metrics['density_kgm3']:.0f} kg/m^3")

    # Show limit comparison
    if "mass_kg" in metrics:
        limit = 1.0  # From mass-check.yaml
        margin = limit - metrics["mass_kg"]
        print(f"\n  Budget:   < {limit:.1f} kg")
        print(f"  Margin:   {margin:.3f} kg ({margin/limit*100:.1f}%)")


def print_fea_results(results: dict, output_path: Path):
    """Print FEA analysis results."""
    print("\n" + "=" * 60)
    print("FEA ANALYSIS RESULTS")
    print("=" * 60)

    status = results["status"]
    metrics = results["metrics"]

    status_str = "PASSED" if status == "passed" else "FAILED" if status == "failed" else status.upper()
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

    print("\nVisualization (open in ParaView):")
    print(f"  paraview {output_path}/validation/results/thrust-fea/displacement.xdmf")
    print(f"  paraview {output_path}/validation/results/thrust-fea/stress.xdmf")


def print_summary(mass_results: dict, fea_results: dict = None):
    """Print combined validation summary."""
    print("\n" + "=" * 60)
    print("VALIDATION SUMMARY")
    print("=" * 60)

    all_passed = True

    # Mass check
    mass_status = mass_results["status"]
    mass_str = "PASS" if mass_status == "passed" else "FAIL"
    if mass_status != "passed":
        all_passed = False
    mass_kg = mass_results["metrics"].get("mass_kg", 0)
    print(f"\n  Mass Budget:  [{mass_str}]  {mass_kg:.3f} kg < 1.0 kg")

    # FEA check
    if fea_results:
        fea_status = fea_results["status"]
        fea_str = "PASS" if fea_status == "passed" else "FAIL"
        if fea_status != "passed":
            all_passed = False
        disp = fea_results["metrics"].get("displacement.max_mm", 0)
        print(f"  FEA Analysis: [{fea_str}]  {disp:.3f} mm < 1.0 mm")
    else:
        print(f"  FEA Analysis: [SKIP]  (use without --skip-fea to run)")

    # Overall
    overall = "ALL CHECKS PASSED" if all_passed else "SOME CHECKS FAILED"
    print(f"\n  Overall: {overall}")


def main():
    parser = argparse.ArgumentParser(
        description="Run validation checks on thrust structure"
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
        help="Skip FEA analysis (run mass check only)"
    )
    parser.add_argument(
        "--mass-only", action="store_true",
        help="Run mass check only (alias for --skip-fea)"
    )

    args = parser.parse_args()

    # --mass-only is alias for --skip-fea
    skip_fea = args.skip_fea or args.mass_only

    # Check dependencies for FEA
    if not skip_fea and not check_dependencies():
        print("\nUse --skip-fea to run mass check without FEA")
        return 1

    # Create package
    if not create_package(args.design, args.output):
        return 1

    print(f"\n   Package created: {args.output}")

    mass_results = None
    fea_results = None

    # Run mass check (always runs - no external deps)
    try:
        mass_results = run_mass_check(args.output)
        print_mass_results(mass_results)
    except Exception as e:
        print(f"\nMass check error: {e}")
        import traceback
        traceback.print_exc()
        return 1

    # Run FEA (if not skipped)
    if not skip_fea:
        try:
            fea_results = run_fea(args.output)
            print_fea_results(fea_results, args.output)
        except Exception as e:
            print(f"\nFEA Error: {e}")
            import traceback
            traceback.print_exc()
            return 1

    # Print combined summary
    print_summary(mass_results, fea_results)

    # Return non-zero if any check failed
    if mass_results["status"] != "passed":
        return 1
    if fea_results and fea_results["status"] != "passed":
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
