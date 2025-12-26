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
                # XDMF export happens in _run_fea
                artifacts.append({"kind": "xdmf", "path": "displacement.xdmf", "description": "Displacement field"})
            if "stress_field" in fea_results:
                artifacts.append({"kind": "xdmf", "path": "stress.xdmf", "description": "Von Mises stress field"})

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
        """Load geometry from the package.

        Returns the first solid from the package geometry.
        """
        from yapcad.package.core import load_geometry
        from yapcad.geom3d import issolid

        # Get geometry source from plan
        geom_spec = plan.geometry
        source = geom_spec.get("source", "geometry/primary.json")
        entities = geom_spec.get("entities", [])

        # Load geometry from package - returns List[list] of entities
        geometry = load_geometry(manifest)

        if entities:
            # Filter to specified entities
            # TODO: Implement entity filtering
            pass

        # Extract the first solid from the geometry list
        # load_geometry returns a list of entities (solids, surfaces, etc.)
        if isinstance(geometry, list):
            for entity in geometry:
                if issolid(entity):
                    return entity
            # If no solid found, maybe the geometry itself is a solid
            if geometry and issolid(geometry[0] if len(geometry) == 1 else geometry):
                return geometry[0] if len(geometry) == 1 else geometry
            raise ValueError("No solid geometry found in package")

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
            algorithm_2d=int(mesh_opts.get("algorithm_2d", 6)),  # Frontal-Delaunay
            algorithm_3d=int(mesh_opts.get("algorithm_3d", 1)),  # Delaunay
            optimize=bool(mesh_opts.get("optimize", True)),
            optimize_netgen=bool(mesh_opts.get("optimize_netgen", False)),
            geometry_tolerance=float(mesh_opts.get("geometry_tolerance", 1e-4)),
            use_stl=bool(mesh_opts.get("use_stl", False)),
            scale_factor=float(mesh_opts.get("scale_factor", 1.0)),
        )

    def _get_physical_groups(self, plan: AnalysisPlan) -> Dict[str, List[int]]:
        """Extract physical groups from plan.

        Note: Currently only supports integer face indices, not named face selectors.
        Named selectors (like 'motor_mount_inner') are skipped with a warning.
        """
        groups: Dict[str, List[int]] = {}

        # From boundary conditions
        for bc in plan.boundary_conditions:
            surfaces = bc.get("surfaces", [])
            bc_id = bc.get("id", "bc")
            if surfaces:
                # Filter to only integer indices, skip string names
                int_surfaces = []
                for s in (surfaces if isinstance(surfaces, list) else [surfaces]):
                    if isinstance(s, int):
                        int_surfaces.append(s)
                    # Skip string names - face naming via selectors not yet implemented
                if int_surfaces:
                    groups[bc_id] = int_surfaces

        # From loads
        for load in plan.loads:
            surfaces = load.get("surfaces", [])
            load_id = load.get("id", "load")
            if surfaces:
                # Filter to only integer indices, skip string names
                int_surfaces = []
                for s in (surfaces if isinstance(surfaces, list) else [surfaces]):
                    if isinstance(s, int):
                        int_surfaces.append(s)
                    # Skip string names - face naming via selectors not yet implemented
                if int_surfaces:
                    groups[load_id] = int_surfaces

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
        mesher = GmshMesher()
        mesher._geometry_tolerance = hints.geometry_tolerance
        mesher._use_stl = hints.use_stl
        mesher._scale_factor = hints.scale_factor
        with mesher:
            # Check if we have STEP file reference
            geom_spec = plan.geometry
            step_file = geom_spec.get("step_file")

            if step_file:
                mesher.import_step(Path(step_file))
            else:
                mesher.import_solid(solid, use_stl=hints.use_stl)

            # Set up physical groups for BCs
            if physical_groups:
                mesher.set_physical_groups(physical_groups, dim=2)

            # Also set up volume group for the entire solid
            import gmsh
            volume_entities = gmsh.model.getEntities(dim=3)
            if volume_entities:
                volume_tags = [tag for dim, tag in volume_entities]
                mesher.set_physical_groups({"volume": volume_tags}, dim=3)

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
        # Note: Geometry is in mm, so we use mm-N-MPa unit system
        # E in MPa = E in Pa / 1e6, pressure in MPa = pressure in Pa / 1e6
        # This gives displacement in mm directly
        material = self._get_material(plan)

        # Scale from Pa to MPa for mm-based geometry
        SCALE_PA_TO_MPA = 1e-6
        lmbda_mpa = material.lame_lambda * SCALE_PA_TO_MPA
        mu_mpa = material.lame_mu * SCALE_PA_TO_MPA

        lmbda = fem.Constant(domain, default_scalar_type(lmbda_mpa))
        mu = fem.Constant(domain, default_scalar_type(mu_mpa))

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

        # Add surface loads (pressure on top surface)
        surface_load = self._get_surface_load_form(domain, v, plan)
        if surface_load is not None:
            L = L + surface_load

        # Boundary conditions
        bcs = self._apply_boundary_conditions(domain, V, plan)

        # Solve
        problem = LinearProblem(
            a, L, bcs=bcs,
            petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
            petsc_options_prefix="yapCAD_"
        )
        uh = problem.solve()

        # Compute metrics
        # Maximum displacement
        # In newer DOLFINx, use .x.array instead of .vector.localForm()
        u_array = uh.x.array.reshape(-1, domain.geometry.dim)
        u_mag = np.sqrt(np.sum(u_array**2, axis=1))
        max_disp = float(np.max(u_mag))

        # With mm-based geometry and MPa material properties:
        # - displacement is already in mm
        # - stress is already in MPa
        results["metrics"]["displacement.max"] = max_disp
        results["metrics"]["displacement.max_mm"] = max_disp  # Already in mm

        # Von Mises stress - compute field and max value
        # Note: Using scaled material properties (MPa), stress is directly in MPa
        vm_func, vm_max = self._compute_von_mises_field(domain, uh, material, scale_to_mpa=True)
        results["metrics"]["stress.von_mises.max"] = vm_max * 1e6  # Store in Pa for consistency
        results["metrics"]["stress.von_mises.max_mpa"] = vm_max  # Already in MPa

        # Export displacement field using XDMF (more compatible with ParaView)
        disp_xdmf_path = workspace / "displacement.xdmf"
        with io.XDMFFile(MPI.COMM_WORLD, str(disp_xdmf_path), "w") as xdmf:
            xdmf.write_mesh(domain)
            uh.name = "displacement"
            xdmf.write_function(uh)

        # Export von Mises stress field for visualization
        stress_xdmf_path = workspace / "stress.xdmf"
        with io.XDMFFile(MPI.COMM_WORLD, str(stress_xdmf_path), "w") as xdmf:
            xdmf.write_mesh(domain)
            vm_func.name = "von_mises_stress"
            xdmf.write_function(vm_func)

        results["displacement_field"] = uh
        results["stress_field"] = vm_func

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
        """Build the body force vector from plan specification."""
        # Default: zero body force
        # Surface loads (pressure, traction) are handled in the weak form
        return fem.Constant(domain, default_scalar_type((0, 0, 0)))

    def _get_surface_load_form(self, domain, v, plan: AnalysisPlan):
        """Build surface load contribution to weak form.

        Supports multiple load application strategies based on plan configuration:
        - 'z_max': Apply load on top surface (highest Z)
        - 'z_min': Apply load on bottom surface (lowest Z)
        - 'z_max_ring': Apply load on ring-shaped region at top (for cradle rings)
        - 'motor_mount': Apply load on central cylindrical region (thrust structure)
        """
        coords = domain.geometry.x
        z_min = float(np.min(coords[:, 2]))
        z_max = float(np.max(coords[:, 2]))
        z_range = z_max - z_min
        z_tol = z_range * 0.05  # 5% tolerance

        # Get geometry parameters from plan metadata
        metadata = plan.metadata or {}
        design_params = metadata.get("design_parameters", {})

        tdim = domain.topology.dim
        fdim = tdim - 1
        domain.topology.create_connectivity(fdim, tdim)

        # Get pressure magnitude and load strategy from plan loads
        pressure_pa = 0.0
        load_strategy = "z_max"  # Default: apply to top surface
        load_direction = [0, 0, -1]  # Default: downward (-Z)

        for load in plan.loads:
            if load.get("type") == "pressure":
                pressure_pa = float(load.get("magnitude_pa", 0.0))
                load_strategy = load.get("strategy", "z_max")
                load_direction = load.get("direction", [0, 0, -1])
                break

        if pressure_pa == 0.0:
            return None

        # Build facet selector based on strategy
        if load_strategy == "motor_mount":
            # Original thrust structure logic
            motor_mount_diameter = float(design_params.get("motor_mount_diameter_mm", 101.6))
            motor_mount_radius = motor_mount_diameter / 2.0

            def load_region(x):
                r = np.sqrt(x[0]**2 + x[1]**2)
                on_top = np.isclose(x[2], z_max, atol=z_tol)
                in_motor_region = r <= motor_mount_radius * 1.1
                return on_top & in_motor_region

        elif load_strategy == "z_max_ring":
            # Ring-shaped load region at top (e.g., cradle ring of globe stand)
            # Load is applied to upper portion of geometry
            z_threshold = z_max - z_range * 0.15  # Top 15% of height

            def load_region(x):
                return x[2] >= z_threshold

        elif load_strategy == "z_min":
            # Bottom surface
            def load_region(x):
                return np.isclose(x[2], z_min, atol=z_tol)

        else:  # "z_max" - default
            # Top surface
            def load_region(x):
                return np.isclose(x[2], z_max, atol=z_tol)

        load_facet_indices = dfx_mesh.locate_entities_boundary(domain, fdim, load_region)

        if len(load_facet_indices) == 0:
            return None

        # Convert pressure from Pa to MPa for mm-based geometry
        SCALE_PA_TO_MPA = 1e-6
        pressure_mpa = pressure_pa * SCALE_PA_TO_MPA

        # Create facet tags for the load surface
        num_facets = domain.topology.index_map(fdim).size_local
        facet_values = np.zeros(num_facets, dtype=np.int32)
        facet_values[load_facet_indices] = 1

        facet_tags = dfx_mesh.meshtags(domain, fdim, np.arange(num_facets, dtype=np.int32), facet_values)
        ds = ufl.Measure("ds", domain=domain, subdomain_data=facet_tags)

        # Apply load in specified direction
        traction = fem.Constant(domain, default_scalar_type((
            load_direction[0] * pressure_mpa,
            load_direction[1] * pressure_mpa,
            load_direction[2] * pressure_mpa
        )))

        return ufl.dot(traction, v) * ds(1)

    def _apply_boundary_conditions(self, domain, V, plan: AnalysisPlan) -> List:
        """Apply boundary conditions from plan.

        Supports multiple BC strategies based on plan configuration:
        - 'z_min': Fix bottom surface (lowest Z) - default for stand-like structures
        - 'z_max': Fix top surface (highest Z)
        - 'stringer_notches': Fix at stringer notch locations (thrust structure)
        - 'outer_edge': Fix at outer radial boundary
        """
        bcs = []

        coords = domain.geometry.x
        z_min = float(np.min(coords[:, 2]))
        z_max = float(np.max(coords[:, 2]))
        z_range = z_max - z_min
        z_tol = z_range * 0.02  # 2% tolerance for surface detection

        # Get geometry parameters from plan metadata
        metadata = plan.metadata or {}
        design_params = metadata.get("design_parameters", {})

        tdim = domain.topology.dim
        fdim = tdim - 1

        for bc_spec in plan.boundary_conditions:
            bc_type = bc_spec.get("type", "fixed")
            bc_strategy = bc_spec.get("strategy", "z_min")  # Default: fix bottom

            if bc_type == "fixed":
                if bc_strategy == "stringer_notches":
                    # Original thrust structure logic
                    outer_diameter = float(design_params.get("outer_diameter_mm", 304.8))
                    outer_radius = outer_diameter / 2.0
                    stringer_width = 25.4  # mm
                    stringer_depth = 12.7  # mm
                    notch_angles_deg = [0.0, 120.0, 240.0]

                    def bc_region(x):
                        r = np.sqrt(x[0]**2 + x[1]**2)
                        theta = np.arctan2(x[1], x[0])
                        near_outer = r >= (outer_radius - stringer_depth * 1.5)
                        half_angle = (stringer_width / 2.0) / outer_radius * 1.5
                        in_notch = np.zeros_like(r, dtype=bool)
                        for angle_deg in notch_angles_deg:
                            angle_rad = np.radians(angle_deg)
                            angle_diff = np.abs(np.arctan2(np.sin(theta - angle_rad), np.cos(theta - angle_rad)))
                            in_notch = in_notch | (angle_diff <= half_angle)
                        return near_outer & in_notch

                elif bc_strategy == "z_max":
                    # Fix top surface
                    def bc_region(x):
                        return np.isclose(x[2], z_max, atol=z_tol)

                elif bc_strategy == "outer_edge":
                    # Fix outer radial boundary
                    x_min, x_max = float(np.min(coords[:, 0])), float(np.max(coords[:, 0]))
                    y_min, y_max = float(np.min(coords[:, 1])), float(np.max(coords[:, 1]))
                    r_max = max(abs(x_min), abs(x_max), abs(y_min), abs(y_max))
                    r_tol = r_max * 0.05

                    def bc_region(x):
                        r = np.sqrt(x[0]**2 + x[1]**2)
                        return r >= (r_max - r_tol)

                else:  # "z_min" - default
                    # Fix bottom surface (most common for stand structures)
                    def bc_region(x):
                        return np.isclose(x[2], z_min, atol=z_tol)

                boundary_facets = dfx_mesh.locate_entities_boundary(
                    domain, fdim, bc_region
                )

                if len(boundary_facets) > 0:
                    boundary_dofs = locate_dofs_topological(V, fdim, boundary_facets)
                    u_bc = np.array([0, 0, 0], dtype=default_scalar_type)
                    bc = dirichletbc(u_bc, boundary_dofs, V)
                    bcs.append(bc)

        return bcs

    def _compute_von_mises_field(self, domain, uh, material: MaterialProperties, scale_to_mpa: bool = False) -> Tuple[Any, float]:
        """Compute von Mises stress field and maximum value.

        Args:
            domain: The mesh domain
            uh: Displacement solution
            material: Material properties
            scale_to_mpa: If True, use MPa-scaled material properties (for mm geometry)

        Returns:
            Tuple of (stress_function, max_stress_value)
        """
        # Define stress tensor
        def epsilon(u):
            return ufl.sym(ufl.grad(u))

        # Use scaled properties if geometry is in mm
        if scale_to_mpa:
            SCALE_PA_TO_MPA = 1e-6
            lmbda = material.lame_lambda * SCALE_PA_TO_MPA
            mu = material.lame_mu * SCALE_PA_TO_MPA
        else:
            lmbda = material.lame_lambda
            mu = material.lame_mu

        def sigma(u):
            return lmbda * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2 * mu * epsilon(u)

        # Von Mises stress: sqrt(3/2 * s:s) where s is deviatoric stress
        s = sigma(uh) - (1/3) * ufl.tr(sigma(uh)) * ufl.Identity(len(uh))
        von_mises = ufl.sqrt((3/2) * ufl.inner(s, s))

        # Project to DG space for evaluation and export
        V_vm = fem.functionspace(domain, ("DG", 0))
        # In newer DOLFINx, interpolation_points is a property, not a method
        vm_expr = fem.Expression(von_mises, V_vm.element.interpolation_points)
        vm_func = fem.Function(V_vm, name="von_mises_stress")
        vm_func.interpolate(vm_expr)

        # In newer DOLFINx, use .x.array instead of .vector.localForm()
        max_stress = float(np.max(np.abs(vm_func.x.array)))

        return vm_func, max_stress

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
