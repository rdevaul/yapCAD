#!/usr/bin/env python3
"""
Globe Stand FEA Example - DSL-to-FEA workflow demonstration.

This example demonstrates:
1. Generating 3D geometry from yapCAD DSL
2. Exporting geometry as STEP files
3. Running FEA analysis with FEniCSx for multiple materials
4. Creating ycpkg packages with analysis results

The globe stand is designed to hold a 30.48cm (12 inch) diameter Mars globe.

Prerequisites for FEA analysis (optional):
    conda install -c conda-forge fenics-dolfinx gmsh meshio

Usage:
    # Generate geometry and STEP export only
    python globe_stand_fea.py

    # Run FEA for all materials (requires FEniCSx)
    python globe_stand_fea.py --fea --material all

    # Create packages with analysis results
    python globe_stand_fea.py --fea --package
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

MATERIALS = {
    "pla": {
        "name": "Generic PLA",
        "youngs_modulus_pa": 3.5e9,
        "poisson_ratio": 0.36,
        "density_kgm3": 1300,
        "yield_strength_pa": 50e6,
        "color": [0.2, 0.6, 0.9],
    },
    "petg": {
        "name": "Generic PETG",
        "youngs_modulus_pa": 2.2e9,
        "poisson_ratio": 0.38,
        "density_kgm3": 1290,
        "yield_strength_pa": 53e6,
        "color": [0.2, 0.9, 0.4],
    },
    "abs": {
        "name": "Generic ABS",
        "youngs_modulus_pa": 2.0e9,
        "poisson_ratio": 0.35,
        "density_kgm3": 1100,
        "yield_strength_pa": 41e6,
        "color": [0.9, 0.7, 0.2],
    },
}


def check_fea_dependencies():
    """Check if FEA dependencies are available."""
    missing = []
    try:
        import gmsh  # noqa: F401
    except ImportError:
        missing.append("gmsh")
    try:
        import dolfinx  # noqa: F401
    except ImportError:
        missing.append("fenics-dolfinx")
    try:
        import meshio  # noqa: F401
    except ImportError:
        missing.append("meshio")

    if missing:
        print("Missing FEA dependencies:", ", ".join(missing))
        print("\nInstall with:")
        print("  conda install -c conda-forge fenics-dolfinx gmsh meshio")
        return False
    return True


#def generate_geometry_from_dsl(dsl_file: Path, command: str = "CENTERED_GLOBE_STAND_ONLY"):
def generate_geometry_from_dsl(dsl_file: Path, command: str = "CENTERED_GLOBE_STAND_HIRES"):
    """Generate geometry from DSL file."""
    from yapcad.dsl import compile_and_run

    print(f"\n1. Generating geometry from DSL")
    print(f"   File: {dsl_file}")
    print(f"   Command: {command}")

    with open(dsl_file) as f:
        source = f.read()

    result = compile_and_run(source, command, {})
    if not result.success:
        raise RuntimeError(f"DSL execution failed: {result.error_message}")

    print("   Geometry generated successfully")
    return result.emit_result.data


def export_step(solid, output_path: Path):
    """Export solid to STEP format."""
    from yapcad.io.step import write_step_analytic

    print(f"\n2. Exporting STEP file")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    write_step_analytic(solid, str(output_path))
    print(f"   Wrote: {output_path}")
    return output_path


def mesh_solid(solid, element_size: float = 5.0, output_dir: Path = None):
    """Mesh the solid geometry using Gmsh."""
    from yapcad.brep import brep_from_solid
    from OCC.Core.BRepTools import breptools
    import gmsh

    if output_dir is None:
        output_dir = Path("/tmp/globe_stand_fea")
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n3. Meshing geometry (element size={element_size}mm)...")

    brep = brep_from_solid(solid)
    if brep is None or brep.shape is None:
        print("   ERROR: Could not extract BREP from solid")
        return None

    brep_file = output_dir / "stand.brep"
    breptools.Write(brep.shape, str(brep_file))

    gmsh.initialize()
    gmsh.model.add("stand")
    gmsh.option.setNumber("General.Terminal", 1)

    try:
        gmsh.model.occ.importShapes(str(brep_file))
        gmsh.model.occ.synchronize()

        volumes = gmsh.model.getEntities(dim=3)
        print(f"   Found {len(volumes)} volumes")
        if not volumes:
            return None

        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", element_size * 0.3)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", element_size)
        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)
        gmsh.option.setNumber("Mesh.Optimize", 1)
        gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)

        gmsh.model.mesh.generate(3)

        element_types, element_tags, _ = gmsh.model.mesh.getElements()
        tetra_count = sum(
            len(tags) for et, tags in zip(element_types, element_tags)
            if "Tetrahedron" in gmsh.model.mesh.getElementProperties(et)[0]
        )
        print(f"   Generated {tetra_count} tetrahedra")

        mesh_file = output_dir / "stand.msh"
        gmsh.write(str(mesh_file))
        print(f"   Wrote mesh: {mesh_file}")
        return mesh_file
    finally:
        gmsh.finalize()


def run_fenics_analysis(mesh_file: Path, material: dict, output_dir: Path,
                        globe_mass_kg: float = 1.0):
    """Run FEniCSx linear elastic analysis."""
    import numpy as np
    import meshio
    import dolfinx
    from dolfinx import fem, io
    from dolfinx.fem.petsc import LinearProblem
    import ufl
    from mpi4py import MPI

    print(f"\n4. Running FEniCSx analysis")
    print(f"   Material: {material['name']}")
    print(f"   E = {material['youngs_modulus_pa']/1e9:.1f} GPa")
    print(f"   Globe mass: {globe_mass_kg} kg")

    output_dir.mkdir(parents=True, exist_ok=True)

    msh = meshio.read(str(mesh_file))
    tetra_cells = next((c.data for c in msh.cells if c.type == "tetra"), None)
    if tetra_cells is None:
        print("   ERROR: No tetrahedral cells found")
        return None

    xdmf_file = output_dir / "mesh.xdmf"
    meshio.write(xdmf_file, meshio.Mesh(points=msh.points, cells=[("tetra", tetra_cells)]),
                 file_format="xdmf")

    with io.XDMFFile(MPI.COMM_WORLD, str(xdmf_file), "r") as xdmf:
        domain = xdmf.read_mesh(name="Grid")

    print(f"   Mesh cells: {domain.topology.index_map(3).size_global}")

    E = material["youngs_modulus_pa"] / 1e6
    nu = material["poisson_ratio"]
    mu = E / (2 * (1 + nu))
    lmbda = E * nu / ((1 + nu) * (1 - 2 * nu))

    V = fem.functionspace(domain, ("Lagrange", 1, (3,)))

    def epsilon(u):
        return ufl.sym(ufl.grad(u))

    def sigma(u):
        return lmbda * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2 * mu * epsilon(u)

    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)

    # Body force: rho [kg/m³] * g [m/s²] / 1e9 = N/mm³
    body_force = material["density_kgm3"] * 9.81 / 1e9
    f_body = fem.Constant(domain, (0.0, 0.0, -body_force))

    # Globe load on cradle (z > 200mm)
    globe_force = globe_mass_kg * 9.81 / 100000.0
    x = ufl.SpatialCoordinate(domain)
    load_factor = ufl.conditional(ufl.gt(x[2], 200.0), 1.0, 0.0)
    f_globe = ufl.as_vector([0.0, 0.0, -globe_force * load_factor])

    a = ufl.inner(sigma(u), epsilon(v)) * ufl.dx
    L = ufl.dot(f_body, v) * ufl.dx + ufl.dot(f_globe, v) * ufl.dx

    def bottom_boundary(x):
        return x[2] < -4.0

    bottom_dofs = fem.locate_dofs_geometrical(V, bottom_boundary)
    bc = fem.dirichletbc(np.array([0.0, 0.0, 0.0], dtype=dolfinx.default_scalar_type),
                         bottom_dofs, V)

    print(f"   Fixed DOFs: {len(bottom_dofs)}")
    print("   Solving...")

    problem = LinearProblem(a, L, bcs=[bc], petsc_options_prefix="fea_",
                            petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
    uh = problem.solve()

    u_values = uh.x.array.reshape(-1, 3)
    max_disp = float(np.max(np.linalg.norm(u_values, axis=1)))
    print(f"   Max displacement: {max_disp:.4f} mm")

    results = {
        "material": material["name"],
        "max_displacement_mm": max_disp,
        "mesh_cells": domain.topology.index_map(3).size_global,
        "fixed_dofs": len(bottom_dofs),
        "globe_mass_kg": globe_mass_kg,
        "status": "PASS" if max_disp <= 5.0 else "FAIL",
    }

    with io.XDMFFile(MPI.COMM_WORLD, str(output_dir / "displacement.xdmf"), "w") as xdmf:
        xdmf.write_mesh(domain)
        uh.name = "displacement"
        xdmf.write_function(uh)

    return results


def print_fea_results(results: list):
    """Print FEA results summary."""
    print("\n" + "=" * 60)
    print("FEA RESULTS SUMMARY")
    print("=" * 60)
    print(f"\n{'Material':<15} {'Max Disp (mm)':<15} {'Status':<10}")
    print("-" * 40)
    for r in results:
        if r:
            print(f"{r['material']:<15} {r['max_displacement_mm']:<15.4f} {r['status']:<10}")
    print("-" * 40)


def create_package(solid, output_dir: Path, material_key: str, fea_results: dict = None):
    """Create a ycpkg package with geometry and FEA results."""
    from yapcad.package import create_package_from_entities, add_geometry_file

    material = MATERIALS[material_key]
    package_dir = output_dir / f"globe_stand_{material_key}.ycpkg"

    print(f"\n5. Creating package: {package_dir}")

    # Create package with solid geometry
    # Note: pass solid directly (not wrapped in extra list)
    manifest = create_package_from_entities(
        [solid],  # List of entities to serialize
        package_dir,
        name=f"Mars Globe Stand ({material['name']})",
        version="1.0.0",
        description=f"Globe stand for 12-inch Mars globe, designed for {material['name']}",
        units="mm",
        materials={
            material_key: {
                "source": {"type": "custom", "custom": {
                    "name": material["name"],
                    "youngs_modulus_pa": material["youngs_modulus_pa"],
                    "poisson_ratio": material["poisson_ratio"],
                    "density_kgm3": material["density_kgm3"],
                }},
                "visual": {"color": material["color"], "metallic": 0.0, "roughness": 0.5},
            }
        },
        overwrite=True,
    )

    # Add STEP export as derived geometry
    step_file = output_dir / f"globe_stand_{material_key}.step"
    if step_file.exists():
        add_geometry_file(manifest, step_file, purpose="STEP export for CAM", overwrite=True)

    # Add FEA results if available
    if fea_results:
        analysis_dir = package_dir / "analysis"
        analysis_dir.mkdir(parents=True, exist_ok=True)
        with open(analysis_dir / "fea_results.json", "w") as f:
            json.dump(fea_results, f, indent=2)

        manifest.data.setdefault("analysis", []).append({
            "type": "fea",
            "method": "linear_elastic",
            "material": material_key,
            "results": {
                "max_displacement_mm": fea_results["max_displacement_mm"],
                "status": fea_results["status"],
            },
            "path": "analysis/fea_results.json",
        })
        manifest.save()

    print(f"   Created: {manifest.manifest_path}")
    return manifest


def main(argv: list[str] = None) -> int:
    parser = argparse.ArgumentParser(description="Globe Stand FEA Example")
    parser.add_argument("--output", type=Path, default=Path("../output/globe_stand_fea"))
    parser.add_argument("--fea", action="store_true", help="Run FEA analysis")
    parser.add_argument("--material", choices=["pla", "petg", "abs", "all"], default="all")
    parser.add_argument("--globe-mass", type=float, default=1.0)
    parser.add_argument("--element-size", type=float, default=5.0)
    parser.add_argument("--package", action="store_true", help="Create ycpkg packages")
    parser.add_argument("--dsl-file", type=Path, default=None)

    args = parser.parse_args(argv)

    if args.dsl_file:
        dsl_file = args.dsl_file
    else:
        dsl_file = Path(__file__).parent / "globe_stand_v4.dsl"
        if not dsl_file.exists():
            print(f"ERROR: DSL file not found: {dsl_file}")
            return 1

    if args.fea and not check_fea_dependencies():
        return 1

    output_dir = args.output
    if not output_dir.is_absolute():
        output_dir = Path(__file__).parent / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("GLOBE STAND FEA EXAMPLE")
    print("=" * 60)

    solid = generate_geometry_from_dsl(dsl_file)
    step_file = output_dir / "globe_stand.step"
    export_step(solid, step_file)

    materials_to_run = list(MATERIALS.keys()) if args.material == "all" else [args.material]
    all_results = {}

    if args.fea:
        mesh_file = mesh_solid(solid, args.element_size, output_dir / "mesh")
        if mesh_file is None:
            print("\nMeshing failed!")
            return 1

        for mat_key in materials_to_run:
            try:
                result = run_fenics_analysis(mesh_file, MATERIALS[mat_key],
                                             output_dir / mat_key, args.globe_mass)
                all_results[mat_key] = result
                if result:
                    with open(output_dir / mat_key / "results.json", "w") as f:
                        json.dump(result, f, indent=2)
            except Exception as e:
                print(f"\nFEA Error for {mat_key}: {e}")
                all_results[mat_key] = None

        print_fea_results(list(all_results.values()))

    if args.package:
        for mat_key in materials_to_run:
            mat_step = output_dir / f"globe_stand_{mat_key}.step"
            if not mat_step.exists():
                shutil.copy(step_file, mat_step)
            create_package(solid, output_dir, mat_key, all_results.get(mat_key))

    print(f"\nOutputs saved to: {output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
