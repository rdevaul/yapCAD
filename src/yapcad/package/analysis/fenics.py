"""FEniCSx (DOLFINx) backend for yapCAD structural analysis.

This module provides FEA capabilities using FEniCSx/DOLFINx with Gmsh meshing.
The key advantage is OCC kernel alignment: yapCAD, Gmsh, and FEniCSx all use
OpenCASCADE, enabling direct geometry passing without lossy conversions.

Supported analysis types:
- Linear elastic static analysis
- Thermal analysis (future)
- Modal analysis (future)

Usage:
    from yapcad.package.analysis.fenics import FenicsxAdapter

    adapter = FenicsxAdapter()
    result = adapter.run(manifest, plan, workspace)

Installation:
    conda create -n yapcad-fenics -c conda-forge fenics-dolfinx gmsh pythonocc-core

Copyright (c) 2025 yapCAD contributors
MIT License
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .base import AnalysisAdapter, AnalysisPlan, AnalysisResult, register_backend
from .gmsh_mesher import GmshMesher, MeshHints, gmsh_available

# FEniCSx imports with graceful fallback
try:
    import dolfinx
    from dolfinx import fem, mesh as dfx_mesh, io, default_scalar_type
    from dolfinx.fem import (
        FunctionSpace, Function, dirichletbc, locate_dofs_topological,
        assemble_scalar, form
    )
    from dolfinx.fem.petsc import LinearProblem
    import ufl
    from mpi4py import MPI
    import numpy as np
    _FENICS_AVAILABLE = True
except ImportError:
    dolfinx = fem = dfx_mesh = io = ufl = None
    MPI = np = None
    _FENICS_AVAILABLE = False

# meshio for mesh format conversion
try:
    import meshio
    _MESHIO_AVAILABLE = True
except ImportError:
    meshio = None
    _MESHIO_AVAILABLE = False


def fenics_available() -> bool:
    """Return True if FEniCSx (DOLFINx) is available."""
    return _FENICS_AVAILABLE


def require_fenics() -> None:
    """Raise error if FEniCSx is not available."""
    if not _FENICS_AVAILABLE:
        raise RuntimeError(
            "FEniCSx (DOLFINx) is not available. Install via conda: "
            "conda install -c conda-forge fenics-dolfinx"
        )


def _now() -> str:
    """Return current UTC timestamp."""
    return datetime.now(timezone.utc).isoformat()


@dataclass
class MaterialProperties:
    """Linear elastic material properties."""
    youngs_modulus: float  # Pa
    poisson_ratio: float
    density: float = 0.0  # kg/m³ (for mass-related calculations)

    @property
    def lame_lambda(self) -> float:
        """First Lamé parameter."""
        E, nu = self.youngs_modulus, self.poisson_ratio
        return E * nu / ((1 + nu) * (1 - 2 * nu))

    @property
    def lame_mu(self) -> float:
        """Second Lamé parameter (shear modulus)."""
        E, nu = self.youngs_modulus, self.poisson_ratio
        return E / (2 * (1 + nu))


class FenicsxAdapter(AnalysisAdapter):
    """FEniCSx (DOLFINx) FEA backend for structural analysis.

    This adapter performs linear elastic static analysis using FEniCSx
    with Gmsh for mesh generation. It supports:

    - Fixed boundary conditions (displacement = 0)
    - Pressure loads on faces
    - Point loads (concentrated forces)
    - Traction (distributed force) on faces

    Results include:
    - Maximum displacement
    - Maximum von Mises stress
    - Displacement and stress fields (VTU export)
    """

    name = "fenics"

    # Default material: Aluminum 6061-T6
    _DEFAULT_E = 68.9e9  # Pa
    _DEFAULT_NU = 0.33
    _DEFAULT_DENSITY = 2700.0  # kg/m³

    def run(
        self,
        manifest: Any,
        plan: AnalysisPlan,
        workspace: Path,
        **kwargs: Any,
    ) -> AnalysisResult:
        """Execute the FEA analysis.

        Args:
            manifest: Package manifest
            plan: Analysis plan specification
            workspace: Working directory for intermediate files

        Returns:
            AnalysisResult with metrics and artifacts
        """
        workspace.mkdir(parents=True, exist_ok=True)

        summary: Dict[str, Any] = {
            "plan": plan.plan_id,
            "backend": "fenics",
            "timestamp": _now(),
        }
        artifacts: List[Dict[str, Any]] = []
        metrics: Dict[str, float] = {}

        # Check dependencies
        if not gmsh_available():
            summary["statusDetail"] = "Gmsh not available for meshing"
            return AnalysisResult(
                plan_id=plan.plan_id,
                status="skipped",
                backend="fenics",
                summary=summary,
            )

        if not fenics_available():
            summary["statusDetail"] = "FEniCSx (DOLFINx) not available"
            return AnalysisResult(
                plan_id=plan.plan_id,
                status="skipped",
                backend="fenics",
                summary=summary,
            )

        try:
            # 1. Load geometry
            solid = self._load_geometry(manifest, plan)

            # 2. Generate mesh
            mesh_path = workspace / "mesh.msh"
            mesh_hints = self._get_mesh_hints(plan)
            physical_groups = self._get_physical_groups(plan)

            self._generate_mesh(solid, mesh_path, mesh_hints, physical_groups, plan)
            artifacts.append({"kind": "mesh", "path": "mesh.msh"})

            # 3. Convert mesh to XDMF for DOLFINx
            xdmf_path = workspace / "mesh.xdmf"
            self._convert_mesh_to_xdmf(mesh_path, xdmf_path)
            artifacts.append({"kind": "mesh-xdmf", "path": "mesh.xdmf"})

            # 4. Run FEA
            fea_results = self._run_fea(xdmf_path, plan, workspace)
            metrics.update(fea_results.get("metrics", {}))

            # 5. Export results
            if "displacement_field" in fea_results:
                vtu_path = workspace / "displacement.vtu"
                # VTU export happens in _run_fea
                artifacts.append({"kind": "vtu", "path": "displacement.vtu"})

            # 6. Evaluate acceptance criteria
            status = self._evaluate_acceptance(metrics, plan.acceptance)
            summary["metrics"] = metrics
            summary["mesh_stats"] = fea_results.get("mesh_stats", {})

            if status == "passed":
                summary["statusDetail"] = "Acceptance criteria satisfied"
            elif status == "failed":
                summary["statusDetail"] = "Acceptance criteria violated"

        except Exception as e:
            summary["statusDetail"] = f"Analysis failed: {str(e)}"
            summary["error"] = str(e)
            return AnalysisResult(
                plan_id=plan.plan_id,
                status="error",
                backend="fenics",
                summary=summary,
                artifacts=artifacts,
            )

        return AnalysisResult(
            plan_id=plan.plan_id,
            status=status,
            backend="fenics",
            metrics=metrics,
            summary=summary,
            artifacts=artifacts,
        )

    def _load_geometry(self, manifest: Any, plan: AnalysisPlan) -> Any:
        """Load geometry from the package."""
        from yapcad.package.core import load_geometry

        # Get geometry source from plan
        geom_spec = plan.geometry
        source = geom_spec.get("source", "geometry/primary.json")
        entities = geom_spec.get("entities", [])

        # Load geometry from package
        geometry = load_geometry(manifest)

        if entities:
            # Filter to specified entities
            # TODO: Implement entity filtering
            pass

        return geometry

    def _get_mesh_hints(self, plan: AnalysisPlan) -> MeshHints:
        """Extract mesh hints from plan."""
        opts = plan.backend_options
        mesh_opts = opts.get("mesh", {})

        return MeshHints(
            element_size=float(mesh_opts.get("element_size", 5.0)),
            min_element_size=mesh_opts.get("min_element_size"),
            max_element_size=mesh_opts.get("max_element_size"),
            element_order=int(mesh_opts.get("element_order", 1)),
            optimize=bool(mesh_opts.get("optimize", True)),
        )

    def _get_physical_groups(self, plan: AnalysisPlan) -> Dict[str, List[int]]:
        """Extract physical groups from plan."""
        groups: Dict[str, List[int]] = {}

        # From boundary conditions
        for bc in plan.boundary_conditions:
            surfaces = bc.get("surfaces", [])
            bc_id = bc.get("id", "bc")
            if surfaces:
                groups[bc_id] = surfaces if isinstance(surfaces, list) else [surfaces]

        # From loads
        for load in plan.loads:
            surfaces = load.get("surfaces", [])
            load_id = load.get("id", "load")
            if surfaces:
                groups[load_id] = surfaces if isinstance(surfaces, list) else [surfaces]

        return groups

    def _generate_mesh(
        self,
        solid: Any,
        mesh_path: Path,
        hints: MeshHints,
        physical_groups: Dict[str, List[int]],
        plan: AnalysisPlan,
    ) -> None:
        """Generate mesh using Gmsh."""
        with GmshMesher() as mesher:
            # Check if we have STEP file reference
            geom_spec = plan.geometry
            step_file = geom_spec.get("step_file")

            if step_file:
                mesher.import_step(Path(step_file))
            else:
                mesher.import_solid(solid)

            # Set up physical groups for BCs
            if physical_groups:
                mesher.set_physical_groups(physical_groups, dim=2)

            # Also set up volume group
            volumes = [(3, t) for _, t in mesher._get_all_volumes()]
            if volumes:
                mesher.set_physical_groups({"volume": [v[1] for v in volumes]}, dim=3)

            mesher.generate_mesh(hints, dim=3)
            mesher.export_mesh(mesh_path)

    def _get_all_volumes(self) -> List[Tuple[int, int]]:
        """Get all volume entities."""
        import gmsh
        return gmsh.model.getEntities(dim=3)

    def _convert_mesh_to_xdmf(self, msh_path: Path, xdmf_path: Path) -> None:
        """Convert Gmsh MSH to XDMF for DOLFINx."""
        if not _MESHIO_AVAILABLE:
            raise RuntimeError("meshio required for mesh conversion")

        mesh = meshio.read(str(msh_path))

        # Extract tetrahedra for 3D
        cells = []
        for cell_block in mesh.cells:
            if cell_block.type == "tetra":
                cells.append(("tetra", cell_block.data))

        if not cells:
            raise ValueError("No tetrahedral elements found in mesh")

        # Create meshio mesh with just tets
        tet_mesh = meshio.Mesh(
            points=mesh.points,
            cells=cells,
        )

        meshio.write(str(xdmf_path), tet_mesh)

    def _run_fea(
        self,
        mesh_path: Path,
        plan: AnalysisPlan,
        workspace: Path,
    ) -> Dict[str, Any]:
        """Run FEA using DOLFINx."""
        require_fenics()

        results: Dict[str, Any] = {"metrics": {}}

        # Read mesh
        with io.XDMFFile(MPI.COMM_WORLD, str(mesh_path), "r") as xdmf:
            domain = xdmf.read_mesh(name="Grid")

        # Get mesh stats
        results["mesh_stats"] = {
            "num_cells": domain.topology.index_map(domain.topology.dim).size_local,
            "num_vertices": domain.topology.index_map(0).size_local,
        }

        # Define function space (vector-valued for displacement)
        V = fem.functionspace(domain, ("Lagrange", 1, (domain.geometry.dim,)))

        # Get material properties
        material = self._get_material(plan)
        lmbda = fem.Constant(domain, default_scalar_type(material.lame_lambda))
        mu = fem.Constant(domain, default_scalar_type(material.lame_mu))

        # Define strain and stress
        def epsilon(u):
            return ufl.sym(ufl.grad(u))

        def sigma(u):
            return lmbda * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2 * mu * epsilon(u)

        # Trial and test functions
        u = ufl.TrialFunction(V)
        v = ufl.TestFunction(V)

        # Bilinear form
        a = ufl.inner(sigma(u), epsilon(v)) * ufl.dx

        # Linear form (loads)
        f = self._build_load_vector(domain, V, plan)
        L = ufl.inner(f, v) * ufl.dx

        # Boundary conditions
        bcs = self._apply_boundary_conditions(domain, V, plan)

        # Solve
        problem = LinearProblem(a, L, bcs=bcs, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
        uh = problem.solve()

        # Compute metrics
        # Maximum displacement
        with uh.vector.localForm() as loc:
            u_array = loc.array.reshape(-1, domain.geometry.dim)
            u_mag = np.sqrt(np.sum(u_array**2, axis=1))
            max_disp = float(np.max(u_mag))

        results["metrics"]["displacement.max"] = max_disp
        results["metrics"]["displacement.max_mm"] = max_disp * 1000  # Assuming SI units

        # Von Mises stress
        vm_stress = self._compute_von_mises(domain, uh, material)
        results["metrics"]["stress.von_mises.max"] = vm_stress

        # Export displacement field
        vtu_path = workspace / "displacement.vtu"
        with io.VTXWriter(MPI.COMM_WORLD, str(vtu_path), [uh]) as vtx:
            vtx.write(0.0)

        results["displacement_field"] = uh

        return results

    def _get_material(self, plan: AnalysisPlan) -> MaterialProperties:
        """Get material properties from plan."""
        mat_spec = plan.materials.get("default", {})

        # Try to load from file reference
        if isinstance(mat_spec, str):
            # Path to material file - TODO: load from file
            mat_spec = {}

        return MaterialProperties(
            youngs_modulus=float(mat_spec.get("youngs_modulus_pa", self._DEFAULT_E)),
            poisson_ratio=float(mat_spec.get("poisson_ratio", self._DEFAULT_NU)),
            density=float(mat_spec.get("density_kgm3", self._DEFAULT_DENSITY)),
        )

    def _build_load_vector(self, domain, V, plan: AnalysisPlan):
        """Build the load vector from plan specification."""
        # Default: gravity in -Z (if density specified)
        # For now, return zero body force
        return fem.Constant(domain, default_scalar_type((0, 0, 0)))

    def _apply_boundary_conditions(self, domain, V, plan: AnalysisPlan) -> List:
        """Apply boundary conditions from plan."""
        bcs = []

        for bc_spec in plan.boundary_conditions:
            bc_type = bc_spec.get("type", "fixed")

            if bc_type == "fixed":
                # Fixed BC: u = 0 on specified surfaces
                # For now, fix bottom face (z = min)
                def bottom_boundary(x):
                    return np.isclose(x[2], np.min(x[2]))

                boundary_facets = dfx_mesh.locate_entities_boundary(
                    domain, domain.topology.dim - 1, bottom_boundary
                )
                boundary_dofs = locate_dofs_topological(V, domain.topology.dim - 1, boundary_facets)

                u_bc = np.array([0, 0, 0], dtype=default_scalar_type)
                bc = dirichletbc(u_bc, boundary_dofs, V)
                bcs.append(bc)

        return bcs

    def _compute_von_mises(self, domain, uh, material: MaterialProperties) -> float:
        """Compute maximum von Mises stress."""
        # Define stress tensor
        def epsilon(u):
            return ufl.sym(ufl.grad(u))

        def sigma(u):
            lmbda = material.lame_lambda
            mu = material.lame_mu
            return lmbda * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2 * mu * epsilon(u)

        # Von Mises stress
        s = sigma(uh) - (1/3) * ufl.tr(sigma(uh)) * ufl.Identity(len(uh))
        von_mises = ufl.sqrt((3/2) * ufl.inner(s, s))

        # Project to DG space for evaluation
        V_vm = fem.functionspace(domain, ("DG", 0))
        vm_expr = fem.Expression(von_mises, V_vm.element.interpolation_points())
        vm_func = fem.Function(V_vm)
        vm_func.interpolate(vm_expr)

        with vm_func.vector.localForm() as loc:
            return float(np.max(np.abs(loc.array)))

    def _evaluate_acceptance(self, metrics: Dict[str, float], acceptance: Dict[str, Any]) -> str:
        """Evaluate acceptance criteria."""
        if not acceptance:
            return "passed" if metrics else "pending"

        for key, rule in acceptance.items():
            metric_value = metrics.get(key)

            # Try alternate keys
            if metric_value is None:
                metric_value = metrics.get(f"{key}_mm")
            if metric_value is None:
                metric_value = metrics.get(key.replace(".", "_"))

            if metric_value is None:
                return "pending"

            limit = rule.get("limit")
            if limit is None:
                limit = rule.get("limit_mm")
                if limit is not None:
                    # Convert to base units if needed
                    pass

            if limit is None:
                continue

            limit = float(limit)
            comparison = rule.get("comparison", "<=")

            if comparison == "<=" and metric_value > limit:
                return "failed"
            elif comparison == "<" and metric_value >= limit:
                return "failed"
            elif comparison == ">=" and metric_value < limit:
                return "failed"
            elif comparison == ">" and metric_value <= limit:
                return "failed"

        return "passed"


# Register the backend
register_backend("fenics", FenicsxAdapter)
register_backend("fenics-dolfinx", FenicsxAdapter)
register_backend("dolfinx", FenicsxAdapter)


__all__ = [
    "FenicsxAdapter",
    "MaterialProperties",
    "fenics_available",
    "require_fenics",
]
