"""Gmsh meshing integration for yapCAD analysis.

This module provides meshing capabilities using Gmsh's OCC integration,
enabling direct geometry transfer from yapCAD's OCC-based BREP representation.

The key advantage is that both yapCAD and Gmsh use the OpenCASCADE kernel,
so geometry can be passed directly without lossy STEP/IGES conversions.

Usage:
    from yapcad.package.analysis.gmsh_mesher import GmshMesher, MeshHints

    mesher = GmshMesher()
    mesh = mesher.mesh_from_solid(solid, hints=MeshHints(element_size=2.0))
    mesher.export_mesh(workspace / "model.msh")

Copyright (c) 2025 yapCAD contributors
MIT License
"""

from __future__ import annotations

import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# Gmsh import with graceful fallback
try:
    import gmsh
    _GMSH_AVAILABLE = True
except ImportError:
    gmsh = None  # type: ignore
    _GMSH_AVAILABLE = False

# OCC imports for geometry conversion
try:
    from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Solid
    from OCC.Core.BRepTools import breptools
    _OCC_AVAILABLE = True
except ImportError:
    TopoDS_Shape = TopoDS_Solid = Any
    breptools = None
    _OCC_AVAILABLE = False


def gmsh_available() -> bool:
    """Return True if Gmsh Python API is available."""
    return _GMSH_AVAILABLE


def require_gmsh() -> None:
    """Raise error if Gmsh is not available."""
    if not _GMSH_AVAILABLE:
        raise RuntimeError(
            "Gmsh Python API is not available. Install via conda: "
            "conda install -c conda-forge gmsh"
        )


@dataclass
class MeshHints:
    """Mesh generation hints for Gmsh.

    Attributes:
        element_size: Target element size (mesh density)
        min_element_size: Minimum element size
        max_element_size: Maximum element size
        algorithm_2d: 2D meshing algorithm (1=MeshAdapt, 2=Auto, 5=Delaunay, 6=Frontal-Delaunay)
        algorithm_3d: 3D meshing algorithm (1=Delaunay, 4=Frontal, 10=HXT)
        element_order: Element polynomial order (1=linear, 2=quadratic)
        optimize: Whether to optimize mesh quality
        optimize_netgen: Use Netgen optimizer for 3D meshes
        refinement_fields: List of refinement field specifications
    """
    element_size: float = 5.0
    min_element_size: Optional[float] = None
    max_element_size: Optional[float] = None
    algorithm_2d: int = 6  # Frontal-Delaunay
    algorithm_3d: int = 1  # Delaunay
    element_order: int = 1
    optimize: bool = True
    optimize_netgen: bool = False
    refinement_fields: List[Dict[str, Any]] = field(default_factory=list)


@dataclass
class PhysicalGroup:
    """A named group of mesh entities (for boundary conditions).

    Attributes:
        name: Human-readable name for the group
        dim: Dimension (0=point, 1=edge, 2=face, 3=volume)
        tags: Gmsh entity tags in this group
    """
    name: str
    dim: int
    tags: List[int]


class GmshMesher:
    """Primary meshing interface using Gmsh's OCC integration.

    This class provides meshing capabilities for yapCAD solids using Gmsh.
    It leverages the shared OCC kernel between yapCAD and Gmsh for direct
    geometry transfer without intermediate file formats.

    Example:
        mesher = GmshMesher()
        mesher.initialize()
        mesher.import_solid(solid)
        mesher.set_physical_groups({"fixed_face": [1, 2], "load_face": [3]})
        mesher.generate_mesh(hints)
        mesher.export_mesh(Path("output.msh"))
        mesher.finalize()
    """

    def __init__(self, model_name: str = "yapCAD_model"):
        """Initialize the mesher.

        Args:
            model_name: Name for the Gmsh model
        """
        require_gmsh()
        self._model_name = model_name
        self._initialized = False
        self._physical_groups: Dict[str, PhysicalGroup] = {}
        self._face_map: Dict[int, str] = {}  # Gmsh tag -> name

    def initialize(self) -> None:
        """Initialize Gmsh (must be called before other operations)."""
        if self._initialized:
            return
        gmsh.initialize()
        gmsh.model.add(self._model_name)
        gmsh.option.setNumber("General.Terminal", 0)  # Suppress terminal output
        self._initialized = True

    def finalize(self) -> None:
        """Finalize Gmsh and release resources."""
        if self._initialized:
            gmsh.finalize()
            self._initialized = False

    def __enter__(self) -> "GmshMesher":
        """Context manager entry."""
        self.initialize()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Context manager exit."""
        self.finalize()

    def import_solid(self, solid: Any, face_names: Optional[Dict[str, List[int]]] = None) -> List[Tuple[int, int]]:
        """Import a yapCAD solid into Gmsh.

        This uses Gmsh's OCC integration to import the geometry directly
        from the OCC representation, avoiding STEP/IGES conversion losses.

        Args:
            solid: yapCAD solid (must have OCC BREP representation)
            face_names: Optional mapping of face names to face indices

        Returns:
            List of (dim, tag) tuples for imported entities
        """
        if not self._initialized:
            raise RuntimeError("GmshMesher not initialized. Call initialize() first.")

        # Get OCC shape from yapCAD solid
        occ_shape = self._get_occ_shape(solid)

        # Import via temporary BREP file (Gmsh's importShapes needs a file)
        # TODO: Investigate direct OCC shape passing when gmsh Python API supports it
        with tempfile.NamedTemporaryFile(suffix=".brep", delete=False) as tmp:
            tmp_path = Path(tmp.name)

        try:
            # Write OCC shape to BREP file
            breptools.Write_s(occ_shape, str(tmp_path))

            # Import into Gmsh via OCC kernel
            entities = gmsh.model.occ.importShapes(str(tmp_path))
            gmsh.model.occ.synchronize()

        finally:
            tmp_path.unlink(missing_ok=True)

        # Store face mapping if provided
        if face_names:
            # Get all faces from the model
            faces = gmsh.model.getEntities(dim=2)
            for name, indices in face_names.items():
                for idx in indices:
                    if idx < len(faces):
                        self._face_map[faces[idx][1]] = name

        return entities

    def import_step(self, step_path: Path) -> List[Tuple[int, int]]:
        """Import geometry from a STEP file.

        Args:
            step_path: Path to the STEP file

        Returns:
            List of (dim, tag) tuples for imported entities
        """
        if not self._initialized:
            raise RuntimeError("GmshMesher not initialized.")

        entities = gmsh.model.occ.importShapes(str(step_path))
        gmsh.model.occ.synchronize()
        return entities

    def set_physical_groups(self, groups: Dict[str, List[int]], dim: int = 2) -> None:
        """Define physical groups for boundary conditions.

        Physical groups associate mesh entities with names that can be
        used to apply boundary conditions in the solver.

        Args:
            groups: Mapping of group names to entity tags
            dim: Dimension of entities (2 for faces, 3 for volumes)
        """
        for name, tags in groups.items():
            if tags:
                pg_tag = gmsh.model.addPhysicalGroup(dim, tags)
                gmsh.model.setPhysicalName(dim, pg_tag, name)
                self._physical_groups[name] = PhysicalGroup(name=name, dim=dim, tags=tags)

    def set_physical_groups_by_normal(
        self,
        groups: Dict[str, Tuple[float, float, float]],
        tolerance_deg: float = 5.0
    ) -> None:
        """Define physical groups by face normal direction.

        This is useful for automatically identifying faces like "top", "bottom",
        "front", etc. based on their orientation.

        Args:
            groups: Mapping of group names to normal vectors (x, y, z)
            tolerance_deg: Angular tolerance in degrees
        """
        import math

        faces = gmsh.model.getEntities(dim=2)
        tol_rad = math.radians(tolerance_deg)

        for name, target_normal in groups.items():
            # Normalize target
            mag = math.sqrt(sum(n*n for n in target_normal))
            if mag < 1e-10:
                continue
            target = tuple(n/mag for n in target_normal)

            matching_tags = []
            for dim, tag in faces:
                # Get face center and normal via Gmsh
                try:
                    # Get parametric center
                    bounds = gmsh.model.getParametrizationBounds(dim, tag)
                    u_mid = (bounds[0][0] + bounds[1][0]) / 2
                    v_mid = (bounds[0][1] + bounds[1][1]) / 2

                    # Get normal at center
                    normal = gmsh.model.getNormal(tag, [u_mid, v_mid])

                    # Check angle
                    dot = sum(a*b for a, b in zip(target, normal))
                    angle = math.acos(max(-1, min(1, dot)))

                    if angle < tol_rad:
                        matching_tags.append(tag)
                except Exception:
                    continue

            if matching_tags:
                pg_tag = gmsh.model.addPhysicalGroup(2, matching_tags)
                gmsh.model.setPhysicalName(2, pg_tag, name)
                self._physical_groups[name] = PhysicalGroup(name=name, dim=2, tags=matching_tags)

    def generate_mesh(self, hints: Optional[MeshHints] = None, dim: int = 3) -> None:
        """Generate the mesh.

        Args:
            hints: Mesh generation hints
            dim: Mesh dimension (2 for surface, 3 for volume)
        """
        if not self._initialized:
            raise RuntimeError("GmshMesher not initialized.")

        hints = hints or MeshHints()

        # Apply mesh size options
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin",
                             hints.min_element_size or hints.element_size * 0.1)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax",
                             hints.max_element_size or hints.element_size * 10)
        gmsh.option.setNumber("Mesh.CharacteristicLengthFactor", 1.0)

        # Set default element size
        gmsh.option.setNumber("Mesh.MeshSizeMin",
                             hints.min_element_size or hints.element_size * 0.1)
        gmsh.option.setNumber("Mesh.MeshSizeMax",
                             hints.max_element_size or hints.element_size * 10)

        # Set meshing algorithms
        gmsh.option.setNumber("Mesh.Algorithm", hints.algorithm_2d)
        gmsh.option.setNumber("Mesh.Algorithm3D", hints.algorithm_3d)

        # Element order
        gmsh.option.setNumber("Mesh.ElementOrder", hints.element_order)

        # Optimization
        gmsh.option.setNumber("Mesh.Optimize", 1 if hints.optimize else 0)
        gmsh.option.setNumber("Mesh.OptimizeNetgen", 1 if hints.optimize_netgen else 0)

        # Apply refinement fields if specified
        for i, field_spec in enumerate(hints.refinement_fields):
            self._apply_refinement_field(i + 1, field_spec)

        # Set uniform mesh size on all entities
        entities = gmsh.model.getEntities(0)  # vertices
        gmsh.model.mesh.setSize(entities, hints.element_size)

        # Generate mesh
        gmsh.model.mesh.generate(dim)

        # Optimize if requested
        if hints.optimize:
            gmsh.model.mesh.optimize("Laplace2D")
            if dim == 3:
                gmsh.model.mesh.optimize("Netgen" if hints.optimize_netgen else "Laplace3D")

    def _apply_refinement_field(self, field_id: int, spec: Dict[str, Any]) -> None:
        """Apply a mesh refinement field."""
        field_type = spec.get("type", "Box")

        gmsh.model.mesh.field.add(field_type, field_id)

        for key, value in spec.items():
            if key == "type":
                continue
            if isinstance(value, (int, float)):
                gmsh.model.mesh.field.setNumber(field_id, key, value)
            elif isinstance(value, list):
                gmsh.model.mesh.field.setNumbers(field_id, key, value)
            elif isinstance(value, str):
                gmsh.model.mesh.field.setString(field_id, key, value)

    def export_mesh(self, path: Path, format: Optional[str] = None) -> Path:
        """Export the mesh to a file.

        Args:
            path: Output path
            format: Optional format override (msh, vtk, xdmf, su2)

        Returns:
            Path to the exported file
        """
        if not self._initialized:
            raise RuntimeError("GmshMesher not initialized.")

        path = Path(path)

        # Determine format from extension if not specified
        if format is None:
            format = path.suffix.lstrip(".").lower()

        # Handle special formats
        if format == "xdmf":
            # XDMF for FEniCSx - export as MSH2 then convert
            msh_path = path.with_suffix(".msh")
            gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
            gmsh.write(str(msh_path))
            # Note: Actual XDMF conversion requires meshio or dolfinx
            return msh_path
        elif format == "su2":
            gmsh.write(str(path))
        else:
            gmsh.write(str(path))

        return path

    def get_mesh_stats(self) -> Dict[str, Any]:
        """Get mesh statistics.

        Returns:
            Dictionary with node count, element counts by type, quality metrics
        """
        if not self._initialized:
            raise RuntimeError("GmshMesher not initialized.")

        stats = {
            "nodes": 0,
            "elements": {},
            "physical_groups": list(self._physical_groups.keys()),
        }

        # Count nodes
        node_tags, coords, _ = gmsh.model.mesh.getNodes()
        stats["nodes"] = len(node_tags)

        # Count elements by type
        element_types, element_tags, _ = gmsh.model.mesh.getElements()
        for elem_type, tags in zip(element_types, element_tags):
            elem_name = gmsh.model.mesh.getElementProperties(elem_type)[0]
            stats["elements"][elem_name] = len(tags)

        return stats

    def _get_occ_shape(self, solid: Any) -> "TopoDS_Shape":
        """Extract OCC shape from yapCAD solid.

        Args:
            solid: yapCAD solid representation

        Returns:
            OCC TopoDS_Shape
        """
        # Check if solid already has OCC representation
        if hasattr(solid, '_occ_shape'):
            return solid._occ_shape

        # Try to get from BREP cache
        from yapcad.brep import BrepSolid, occ_available
        if occ_available():
            if isinstance(solid, BrepSolid):
                return solid.shape

        # Convert via native BREP
        from yapcad.occ_native_convert import native_brep_to_occ, occ_available as occ_convert_available
        if occ_convert_available():
            # yapCAD solid is typically a list of surfaces
            # Need to build OCC solid from it
            occ_shape = native_brep_to_occ(solid)
            if occ_shape is not None:
                return occ_shape

        raise ValueError(
            "Cannot extract OCC shape from solid. Ensure the solid has "
            "a valid BREP representation (created with OCC-enabled yapCAD)."
        )


def mesh_solid(
    solid: Any,
    output_path: Path,
    hints: Optional[MeshHints] = None,
    physical_groups: Optional[Dict[str, List[int]]] = None,
    dim: int = 3
) -> Dict[str, Any]:
    """Convenience function to mesh a solid and export.

    Args:
        solid: yapCAD solid to mesh
        output_path: Path for output mesh file
        hints: Mesh generation hints
        physical_groups: Optional face groups for BCs
        dim: Mesh dimension

    Returns:
        Mesh statistics dictionary
    """
    with GmshMesher() as mesher:
        mesher.import_solid(solid)
        if physical_groups:
            mesher.set_physical_groups(physical_groups)
        mesher.generate_mesh(hints, dim)
        mesher.export_mesh(output_path)
        return mesher.get_mesh_stats()


__all__ = [
    "GmshMesher",
    "MeshHints",
    "PhysicalGroup",
    "gmsh_available",
    "require_gmsh",
    "mesh_solid",
]
