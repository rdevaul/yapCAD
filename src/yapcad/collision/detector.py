"""Main collision detection class with BREP/mesh/AABB methods.

This module defines the CollisionDetector class and GeometryProvider protocol
for performing collision detection on assemblies.

Copyright (c) 2026 yapCAD contributors
License: MIT
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any, Protocol, Union, runtime_checkable

from .result import CollisionResult, CollisionMethod
from .registry import InterfaceRegistry

# ============================================================================
# Optional dependency detection
# ============================================================================

HAVE_NUMPY = False
try:
    import numpy as np
    HAVE_NUMPY = True
except ImportError:
    np = None

HAVE_OCC = False
try:
    from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common
    from OCC.Core.GProp import GProp_GProps
    from OCC.Core.BRepGProp import brepgprop
    from OCC.Core.STEPControl import STEPControl_Reader
    from OCC.Core.StlAPI import StlAPI_Reader
    from OCC.Core.TopoDS import TopoDS_Shape
    from OCC.Core.IFSelect import IFSelect_RetDone
    from OCC.Core.gp import gp_Trsf
    from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_Transform
    HAVE_OCC = True
except ImportError:
    pass

HAVE_TRIMESH = False
try:
    import trimesh
    HAVE_TRIMESH = True
except ImportError:
    trimesh = None


# ============================================================================
# Geometry Provider Protocol
# ============================================================================

@runtime_checkable
class GeometryProvider(Protocol):
    """Protocol for providing geometry data for collision detection.

    Implementations must provide a way to retrieve geometry for parts
    by name. The geometry can be returned as:
        - Path to STEP file (preferred for BREP detection)
        - Path to STL file (for mesh detection)
        - trimesh.Trimesh object (pre-loaded mesh)
        - TopoDS_Shape object (pre-loaded BREP)

    Example:
        >>> class FileBasedProvider:
        ...     def __init__(self, base_dir: Path):
        ...         self.base_dir = base_dir
        ...
        ...     def get_geometry(self, part_name: str) -> Optional[str]:
        ...         step_path = self.base_dir / f"{part_name}.step"
        ...         if step_path.exists():
        ...             return str(step_path)
        ...         stl_path = self.base_dir / f"{part_name}.stl"
        ...         if stl_path.exists():
        ...             return str(stl_path)
        ...         return None
    """

    def get_geometry(self, part_name: str) -> Optional[Any]:
        """Get geometry for a part by name.

        Args:
            part_name: Name/identifier of the part

        Returns:
            Path string, mesh object, or BREP shape for the part,
            or None if geometry not available
        """
        ...


# ============================================================================
# Collision Detector
# ============================================================================

class CollisionDetector:
    """Main collision detection class with multiple detection methods.

    CollisionDetector provides collision detection between assembly parts
    using the best available method:

    1. **BREP detection** (if pythonocc available):
       Uses boolean intersection for exact collision volume

    2. **Mesh detection** (if trimesh available):
       Uses point sampling and containment checks

    3. **AABB detection** (always available):
       Fast bounding box checks (used as pre-filter)

    The detector integrates with InterfaceRegistry to distinguish
    expected overlaps (gear mesh, threads) from actual collisions.

    Attributes:
        geometry_provider: Provider for part geometry
        interface_registry: Registry of allowed overlap interfaces
        verbose: Enable detailed logging
        min_collision_volume: Minimum volume to report (mm^3)
        mesh_sample_points: Number of sample points for mesh detection
        penetration_threshold: Minimum depth to count as collision (mm)

    Example:
        >>> provider = MyGeometryProvider()
        >>> detector = CollisionDetector(provider)
        >>>
        >>> # Optional: register allowed overlaps
        >>> registry = InterfaceRegistry()
        >>> registry.register(GearMeshInterface(...))
        >>> detector.set_interface_registry(registry)
        >>>
        >>> # Check single pair
        >>> result = detector.check_collision(
        ...     "PART_A", transform_a,
        ...     "PART_B", transform_b
        ... )
        >>>
        >>> # Check full assembly
        >>> transforms = {"PART_A": tf_a, "PART_B": tf_b, "PART_C": tf_c}
        >>> results = detector.check_assembly(transforms)
        >>> collisions = [r for r in results if r.is_error]
    """

    def __init__(
        self,
        geometry_provider: GeometryProvider,
        interface_registry: Optional[InterfaceRegistry] = None,
        verbose: bool = False,
        min_collision_volume: float = 0.1,
        mesh_sample_points: int = 1500,
        penetration_threshold: float = 0.1
    ):
        """Initialize collision detector.

        Args:
            geometry_provider: Provider for part geometry
            interface_registry: Registry of allowed overlaps (optional)
            verbose: Enable detailed logging
            min_collision_volume: Minimum collision volume to report (mm^3)
            mesh_sample_points: Number of sample points for mesh detection
            penetration_threshold: Minimum penetration to count as collision (mm)
        """
        self.geometry_provider = geometry_provider
        self.interface_registry = interface_registry or InterfaceRegistry()
        self.verbose = verbose
        self.min_collision_volume = min_collision_volume
        self.mesh_sample_points = mesh_sample_points
        self.penetration_threshold = penetration_threshold

        # Geometry cache
        self._geometry_cache: Dict[str, Any] = {}

        # Detection method availability
        self._has_brep = HAVE_OCC
        self._has_mesh = HAVE_TRIMESH and HAVE_NUMPY

    def set_interface_registry(self, registry: InterfaceRegistry) -> None:
        """Set the interface registry for allowed overlaps.

        Args:
            registry: InterfaceRegistry with interface definitions
        """
        self.interface_registry = registry

    def clear_cache(self) -> None:
        """Clear the geometry cache."""
        self._geometry_cache.clear()

    @property
    def preferred_method(self) -> CollisionMethod:
        """Get the preferred detection method based on available libraries."""
        if self._has_brep:
            return CollisionMethod.BREP
        if self._has_mesh:
            return CollisionMethod.MESH
        return CollisionMethod.AABB

    def check_collision(
        self,
        part_a: str,
        transform_a: Any,
        part_b: str,
        transform_b: Any,
        method: Optional[CollisionMethod] = None
    ) -> CollisionResult:
        """Check for collision between two parts.

        Args:
            part_a: Name of first part
            transform_a: 4x4 transformation matrix for part A
            part_b: Name of second part
            transform_b: 4x4 transformation matrix for part B
            method: Specific detection method (default: best available)

        Returns:
            CollisionResult with collision details
        """
        # Get geometry for both parts
        geom_a = self._get_geometry(part_a)
        geom_b = self._get_geometry(part_b)

        if geom_a is None:
            return CollisionResult.error(
                part_a, part_b,
                f"Geometry not available for part: {part_a}"
            )

        if geom_b is None:
            return CollisionResult.error(
                part_a, part_b,
                f"Geometry not available for part: {part_b}"
            )

        # Convert transforms to numpy arrays if needed
        tf_a = self._to_numpy_transform(transform_a)
        tf_b = self._to_numpy_transform(transform_b)

        if tf_a is None or tf_b is None:
            return CollisionResult.error(
                part_a, part_b,
                "Invalid transformation matrix"
            )

        # Select detection method
        if method is None:
            method = self.preferred_method

        # Perform collision detection
        if method == CollisionMethod.BREP and self._has_brep:
            result = self._check_brep_collision(
                part_a, geom_a, tf_a,
                part_b, geom_b, tf_b
            )
        elif method == CollisionMethod.MESH and self._has_mesh:
            result = self._check_mesh_collision(
                part_a, geom_a, tf_a,
                part_b, geom_b, tf_b
            )
        else:
            result = self._check_aabb_collision(
                part_a, geom_a, tf_a,
                part_b, geom_b, tf_b
            )

        # Check interface compatibility if collision detected
        if result.collides:
            is_compatible, compat_results = self.interface_registry.check_overlap_compatibility(
                part_a, part_b
            )
            if is_compatible:
                result.compatible_interface = True
                result.interface_names = self.interface_registry.get_compatible_interface_names(
                    part_a, part_b
                )

        return result

    def check_assembly(
        self,
        world_transforms: Dict[str, Any],
        skip_pairs: Optional[List[Tuple[str, str]]] = None
    ) -> List[CollisionResult]:
        """Check all pairs in an assembly for collisions.

        Args:
            world_transforms: Dictionary mapping part names to 4x4 transforms
            skip_pairs: Optional list of (part_a, part_b) pairs to skip

        Returns:
            List of CollisionResult for all pairs checked
        """
        skip_set = set()
        if skip_pairs:
            for a, b in skip_pairs:
                skip_set.add((a, b))
                skip_set.add((b, a))

        results = []
        part_names = list(world_transforms.keys())

        for i, part_a in enumerate(part_names):
            for part_b in part_names[i+1:]:
                # Skip if in skip list
                if (part_a, part_b) in skip_set:
                    continue

                # Get transforms
                tf_a = world_transforms[part_a]
                tf_b = world_transforms[part_b]

                # Check collision
                result = self.check_collision(part_a, tf_a, part_b, tf_b)
                results.append(result)

                if self.verbose and result.collides:
                    print(f"  {result}")

        return results

    def get_collision_summary(
        self,
        results: List[CollisionResult]
    ) -> Dict[str, Any]:
        """Generate summary statistics from collision results.

        Args:
            results: List of CollisionResult from check_assembly

        Returns:
            Dictionary with summary statistics
        """
        total = len(results)
        collisions = [r for r in results if r.collides]
        errors = [r for r in collisions if r.is_error]
        interfaces = [r for r in collisions if r.is_interface_overlap]
        failed = [r for r in results if r.error_message]

        return {
            "total_pairs": total,
            "collisions_detected": len(collisions),
            "actual_errors": len(errors),
            "expected_interfaces": len(interfaces),
            "detection_failures": len(failed),
            "error_pairs": [(r.part_a, r.part_b) for r in errors],
            "interface_pairs": [(r.part_a, r.part_b) for r in interfaces],
            "method": self.preferred_method.name,
        }

    # ========================================================================
    # Internal methods
    # ========================================================================

    def _get_geometry(self, part_name: str) -> Optional[Any]:
        """Get geometry for a part, using cache if available."""
        if part_name in self._geometry_cache:
            return self._geometry_cache[part_name]

        geom = self.geometry_provider.get_geometry(part_name)
        if geom is not None:
            self._geometry_cache[part_name] = geom

        return geom

    def _to_numpy_transform(self, transform: Any) -> Optional[Any]:
        """Convert transform to numpy 4x4 array."""
        if not HAVE_NUMPY:
            return transform  # Return as-is if no numpy

        if transform is None:
            return np.eye(4)

        # Handle Transform objects with .matrix attribute
        if hasattr(transform, 'matrix'):
            transform = transform.matrix

        try:
            arr = np.array(transform, dtype=np.float64)
            if arr.shape == (4, 4):
                return arr
            elif arr.shape == (16,):
                return arr.reshape(4, 4)
            else:
                return None
        except Exception:
            return None

    def _check_brep_collision(
        self,
        part_a: str, geom_a: Any, tf_a: Any,
        part_b: str, geom_b: Any, tf_b: Any
    ) -> CollisionResult:
        """Check collision using BREP boolean intersection."""
        if not HAVE_OCC:
            return CollisionResult.error(
                part_a, part_b,
                "pythonocc not available for BREP collision"
            )

        try:
            # Load/convert shapes
            shape_a = self._load_brep_shape(geom_a, tf_a)
            shape_b = self._load_brep_shape(geom_b, tf_b)

            if shape_a is None or shape_b is None:
                return CollisionResult.error(
                    part_a, part_b,
                    "Failed to load BREP shape"
                )

            # Perform boolean intersection
            common = BRepAlgoAPI_Common(shape_a, shape_b)
            common.Build()

            if not common.IsDone():
                return CollisionResult.error(
                    part_a, part_b,
                    "BREP intersection failed"
                )

            # Calculate intersection volume
            intersection_shape = common.Shape()
            props = GProp_GProps()
            brepgprop.VolumeProperties(intersection_shape, props)
            volume = abs(props.Mass())

            collides = volume > self.min_collision_volume

            return CollisionResult(
                part_a=part_a,
                part_b=part_b,
                collides=collides,
                method=CollisionMethod.BREP,
                intersection_volume=volume,
            )

        except Exception as e:
            return CollisionResult.error(
                part_a, part_b,
                f"BREP collision check failed: {e}"
            )

    def _load_brep_shape(self, geom: Any, transform: Any) -> Optional[Any]:
        """Load geometry as BREP shape with transform applied."""
        if not HAVE_OCC:
            return None

        shape = None

        # Handle different geometry types
        if isinstance(geom, str):
            path = Path(geom)
            if path.suffix.lower() == '.step' or path.suffix.lower() == '.stp':
                reader = STEPControl_Reader()
                status = reader.ReadFile(str(path))
                if status == IFSelect_RetDone:
                    reader.TransferRoots()
                    shape = reader.OneShape()
            elif path.suffix.lower() == '.stl':
                reader = StlAPI_Reader()
                shape = TopoDS_Shape()
                if not reader.Read(shape, str(path)):
                    shape = None
        elif hasattr(geom, 'ShapeType'):
            # Already a TopoDS_Shape
            shape = geom

        if shape is None:
            return None

        # Apply transformation
        if HAVE_NUMPY:
            trsf = gp_Trsf()
            tf = np.array(transform)
            trsf.SetValues(
                tf[0, 0], tf[0, 1], tf[0, 2], tf[0, 3],
                tf[1, 0], tf[1, 1], tf[1, 2], tf[1, 3],
                tf[2, 0], tf[2, 1], tf[2, 2], tf[2, 3]
            )
            transformer = BRepBuilderAPI_Transform(shape, trsf, True)
            if transformer.IsDone():
                return transformer.Shape()

        return shape

    def _check_mesh_collision(
        self,
        part_a: str, geom_a: Any, tf_a: Any,
        part_b: str, geom_b: Any, tf_b: Any
    ) -> CollisionResult:
        """Check collision using mesh sampling and containment."""
        if not HAVE_TRIMESH or not HAVE_NUMPY:
            return CollisionResult.error(
                part_a, part_b,
                "trimesh/numpy not available for mesh collision"
            )

        try:
            # Load meshes
            mesh_a = self._load_mesh(geom_a, tf_a)
            mesh_b = self._load_mesh(geom_b, tf_b)

            if mesh_a is None or mesh_b is None:
                return CollisionResult.error(
                    part_a, part_b,
                    "Failed to load mesh"
                )

            # Quick AABB check
            if not self._aabb_intersects(mesh_a.bounds, mesh_b.bounds):
                return CollisionResult.no_collision(
                    part_a, part_b, CollisionMethod.MESH
                )

            # Sample points and check containment
            n_samples = self.mesh_sample_points

            # Sample from mesh B, check in mesh A
            samples_b = mesh_b.sample(n_samples)
            contains_a = mesh_a.contains(samples_b)
            points_b_in_a = np.sum(contains_a)

            # Sample from mesh A, check in mesh B
            samples_a = mesh_a.sample(n_samples)
            contains_b = mesh_b.contains(samples_a)
            points_a_in_b = np.sum(contains_b)

            total_inside = points_b_in_a + points_a_in_b

            # Calculate penetration depth if overlap detected
            penetration_depth = 0.0
            contact_points = []

            if total_inside > 0:
                # Estimate penetration using signed distance
                try:
                    if points_b_in_a > 0:
                        prox_a = trimesh.proximity.ProximityQuery(mesh_a)
                        dists = prox_a.signed_distance(samples_b[contains_a])
                        if len(dists) > 0:
                            penetration_depth = max(penetration_depth, abs(np.min(dists)))

                    if points_a_in_b > 0:
                        prox_b = trimesh.proximity.ProximityQuery(mesh_b)
                        dists = prox_b.signed_distance(samples_a[contains_b])
                        if len(dists) > 0:
                            penetration_depth = max(penetration_depth, abs(np.min(dists)))
                except Exception:
                    pass

                # Get contact points
                if points_b_in_a > 0:
                    contact_points.extend([tuple(p) for p in samples_b[contains_a][:10]])
                if points_a_in_b > 0:
                    contact_points.extend([tuple(p) for p in samples_a[contains_b][:10]])

            # Determine if collision
            collides = (
                total_inside > 10 or
                (penetration_depth > self.penetration_threshold and total_inside > 0)
            )

            # Estimate intersection volume if collision
            intersection_volume = None
            if collides and total_inside > 10:
                fraction = total_inside / (2 * n_samples)
                smaller_vol = min(mesh_a.volume, mesh_b.volume)
                intersection_volume = fraction * smaller_vol

            return CollisionResult(
                part_a=part_a,
                part_b=part_b,
                collides=collides,
                method=CollisionMethod.MESH,
                intersection_volume=intersection_volume,
                penetration_depth=penetration_depth,
                contact_points=contact_points[:20],
                metadata={"points_inside": int(total_inside)}
            )

        except Exception as e:
            return CollisionResult.error(
                part_a, part_b,
                f"Mesh collision check failed: {e}"
            )

    def _load_mesh(self, geom: Any, transform: Any) -> Optional[Any]:
        """Load geometry as trimesh with transform applied."""
        if not HAVE_TRIMESH:
            return None

        mesh = None

        # Handle different geometry types
        if isinstance(geom, str):
            try:
                mesh = trimesh.load(geom, force='mesh')
            except Exception:
                return None
        elif isinstance(geom, trimesh.Trimesh):
            mesh = geom.copy()
        elif hasattr(geom, 'geometry'):
            # Scene object
            try:
                mesh = trimesh.util.concatenate(list(geom.geometry.values()))
            except Exception:
                return None

        if mesh is None:
            return None

        # Apply transformation
        if HAVE_NUMPY:
            tf = np.array(transform)
            mesh.apply_transform(tf)

        return mesh

    def _check_aabb_collision(
        self,
        part_a: str, geom_a: Any, tf_a: Any,
        part_b: str, geom_b: Any, tf_b: Any
    ) -> CollisionResult:
        """Check collision using axis-aligned bounding boxes only."""
        if not HAVE_TRIMESH or not HAVE_NUMPY:
            # Cannot compute AABB without trimesh/numpy
            return CollisionResult(
                part_a=part_a,
                part_b=part_b,
                collides=True,  # Conservative: assume collision
                method=CollisionMethod.AABB,
                metadata={"warning": "AABB check unavailable, assuming collision"}
            )

        try:
            mesh_a = self._load_mesh(geom_a, tf_a)
            mesh_b = self._load_mesh(geom_b, tf_b)

            if mesh_a is None or mesh_b is None:
                return CollisionResult.error(
                    part_a, part_b,
                    "Failed to load mesh for AABB check"
                )

            collides = self._aabb_intersects(mesh_a.bounds, mesh_b.bounds)

            return CollisionResult(
                part_a=part_a,
                part_b=part_b,
                collides=collides,
                method=CollisionMethod.AABB,
            )

        except Exception as e:
            return CollisionResult.error(
                part_a, part_b,
                f"AABB collision check failed: {e}"
            )

    def _aabb_intersects(self, bounds_a: Any, bounds_b: Any) -> bool:
        """Check if two axis-aligned bounding boxes intersect."""
        if not HAVE_NUMPY:
            return True  # Conservative

        # bounds format: [[min_x, min_y, min_z], [max_x, max_y, max_z]]
        for i in range(3):
            if bounds_a[1, i] < bounds_b[0, i] or bounds_b[1, i] < bounds_a[0, i]:
                return False
        return True
