"""Native yapCAD backend for geometric and measurement validation tests.

This module provides validation capabilities using yapCAD's built-in geometry
functions. No external solvers are required.

Supported test kinds:
- geometric: volume, area, bbox checks
- measurement: mass, centroid calculations

Usage:
    from yapcad.package.analysis.yapcad_native import YapCADNativeAdapter

    adapter = YapCADNativeAdapter()
    result = adapter.run(manifest, plan, workspace)

Copyright (c) 2025 yapCAD contributors
MIT License
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

from .base import AnalysisAdapter, AnalysisPlan, AnalysisResult, register_backend


def _now() -> str:
    """Return current UTC timestamp."""
    return datetime.now(timezone.utc).isoformat()


class YapCADNativeAdapter(AnalysisAdapter):
    """Native yapCAD backend for geometric validation tests.

    This adapter performs geometric and measurement checks using yapCAD's
    built-in functions. It supports:

    - Volume checks (solid volume against limits)
    - Area checks (surface area or 2D region area)
    - Bounding box checks (dimensions, diagonal)
    - Mass checks (volume * density)
    - Centroid checks (center of mass location)

    Results include computed metrics and pass/fail status based on
    acceptance criteria.
    """

    name = "yapcad"

    def run(
        self,
        manifest: Any,
        plan: AnalysisPlan,
        workspace: Path,
        **kwargs: Any,
    ) -> AnalysisResult:
        """Execute the validation check.

        Args:
            manifest: Package manifest
            plan: Analysis plan specification
            workspace: Working directory for intermediate files

        Returns:
            AnalysisResult with metrics and pass/fail status
        """
        workspace.mkdir(parents=True, exist_ok=True)

        summary: Dict[str, Any] = {
            "plan": plan.plan_id,
            "backend": "yapcad",
            "timestamp": _now(),
        }
        artifacts: List[Dict[str, Any]] = []
        metrics: Dict[str, Any] = {}

        try:
            # Load geometry from package
            solid = self._load_geometry(manifest, plan)

            # Get the check specification
            check_spec = plan.raw.get("check", {})
            check_property = check_spec.get("property", "volume")

            # Dispatch to appropriate check method
            if check_property == "volume":
                metrics.update(self._check_volume(solid, check_spec))
            elif check_property == "area":
                metrics.update(self._check_area(solid, check_spec))
            elif check_property == "bbox":
                metrics.update(self._check_bbox(solid, check_spec))
            elif check_property == "mass":
                metrics.update(self._check_mass(solid, check_spec))
            elif check_property == "centroid":
                metrics.update(self._check_centroid(solid, check_spec))
            elif check_property == "clearance":
                # Clearance requires multiple entities
                metrics.update(self._check_clearance(manifest, plan, check_spec))
            else:
                raise ValueError(f"Unknown check property: {check_property}")

            # Evaluate acceptance criteria
            status = self._evaluate_acceptance(metrics, plan.acceptance)
            summary["metrics"] = metrics

            if status == "passed":
                summary["statusDetail"] = "Acceptance criteria satisfied"
            elif status == "failed":
                summary["statusDetail"] = "Acceptance criteria violated"
                summary["failed_criteria"] = self._get_failed_criteria(metrics, plan.acceptance)

        except Exception as e:
            summary["statusDetail"] = f"Check failed: {str(e)}"
            summary["error"] = str(e)
            return AnalysisResult(
                plan_id=plan.plan_id,
                status="error",
                backend="yapcad",
                summary=summary,
                artifacts=artifacts,
            )

        # Write summary to workspace
        summary_path = workspace / "summary.json"
        with open(summary_path, "w") as f:
            json.dump({
                "plan_id": plan.plan_id,
                "status": status,
                "timestamp": summary["timestamp"],
                "backend": "yapcad",
                "metrics": metrics,
                "acceptance_results": self._build_acceptance_results(metrics, plan.acceptance),
            }, f, indent=2)
        artifacts.append({"kind": "summary", "path": "summary.json"})

        return AnalysisResult(
            plan_id=plan.plan_id,
            status=status,
            backend="yapcad",
            metrics=metrics,
            summary=summary,
            artifacts=artifacts,
            summary_path=summary_path,
        )

    def _load_geometry(self, manifest: Any, plan: AnalysisPlan) -> Any:
        """Load geometry from the package.

        Returns either a single solid or a list of solids (assembly).
        For aggregate checks like total mass or volume, all solids are used.
        """
        from yapcad.package.core import load_geometry
        from yapcad.geom3d import issolid

        # Get geometry source from plan
        geom_spec = plan.geometry
        entities = geom_spec.get("entities", [])
        aggregate = geom_spec.get("aggregate", False)  # Return all solids for assembly checks

        # Load geometry from package
        geometry = load_geometry(manifest)

        if entities:
            # Filter to specified entities (future: implement entity lookup)
            pass

        # Extract solids from geometry list
        if isinstance(geometry, list):
            solids = [entity for entity in geometry if issolid(entity)]
            if not solids:
                # Check if geometry itself is a solid
                if geometry and issolid(geometry[0] if len(geometry) == 1 else geometry):
                    return geometry[0] if len(geometry) == 1 else geometry
                raise ValueError("No solid geometry found in package")

            # Return all solids for aggregate checks, or first solid otherwise
            if aggregate or len(solids) > 1:
                return solids
            return solids[0]

        return geometry

    def _normalize_geometry(self, geometry: Any) -> List[Any]:
        """Normalize geometry to a list of solids for aggregate operations."""
        from yapcad.geom3d import issolid

        if isinstance(geometry, list):
            # Could be list of solids or a single solid structure
            if geometry and issolid(geometry):
                return [geometry]
            return [g for g in geometry if issolid(g)]
        elif issolid(geometry):
            return [geometry]
        return []

    def _check_volume(self, geometry: Any, check_spec: Dict[str, Any]) -> Dict[str, float]:
        """Compute volume of a solid or assembly.

        For assemblies (list of solids), computes total volume.
        """
        from yapcad.geom3d import volumeof

        solids = self._normalize_geometry(geometry)
        if not solids:
            raise ValueError("No solid geometry to check")

        total_volume = 0.0
        part_volumes = []

        for solid in solids:
            vol = volumeof(solid)
            if vol is None:
                raise ValueError("Could not compute volume - solid may lack BREP data")
            total_volume += vol
            part_volumes.append(vol)

        units = check_spec.get("units", "mm3")

        metrics = {
            "volume": total_volume,
            f"volume_{units}": total_volume,
            "part_count": len(solids),
        }

        # Include individual part volumes for assemblies
        if len(solids) > 1:
            for i, vol in enumerate(part_volumes):
                metrics[f"volume.part_{i}"] = vol

        return metrics

    def _check_area(self, geometry: Any, check_spec: Dict[str, Any]) -> Dict[str, float]:
        """Compute surface area of a solid or assembly.

        Uses BREP if available, falls back to tessellation.
        For assemblies, computes total surface area.
        """
        solids = self._normalize_geometry(geometry)
        if not solids:
            raise ValueError("No solid geometry to check")

        total_area = 0.0
        part_areas = []

        for solid in solids:
            area = self._compute_solid_area(solid)
            if area is None:
                raise ValueError("Could not compute surface area for one or more parts")
            total_area += area
            part_areas.append(area)

        units = check_spec.get("units", "mm2")

        metrics = {
            "area": total_area,
            f"area_{units}": total_area,
            "part_count": len(solids),
        }

        if len(solids) > 1:
            for i, a in enumerate(part_areas):
                metrics[f"area.part_{i}"] = a

        return metrics

    def _compute_solid_area(self, solid: Any) -> Optional[float]:
        """Compute surface area of a single solid."""
        area = None

        # Try BREP first
        try:
            from yapcad.brep import brep_from_solid
            from OCC.Core.GProp import GProp_GProps
            from OCC.Core.BRepGProp import brepgprop

            brep = brep_from_solid(solid)
            if brep is not None and brep.shape is not None:
                props = GProp_GProps()
                brepgprop.SurfaceProperties(brep.shape, props)
                area = props.Mass()  # For surfaces, Mass() returns area
        except Exception:
            pass

        # Fall back to tessellation if BREP failed
        if area is None:
            area = self._mesh_surface_area(solid)

        return area

    def _mesh_surface_area(self, solid: Any) -> Optional[float]:
        """Compute surface area from triangle mesh."""
        from yapcad.geom3d import surfacearea as surface_area_tess

        surfaces = solid[1] if len(solid) > 1 else []
        if not surfaces:
            return None

        return sum(surface_area_tess(surf) for surf in surfaces)

    def _check_bbox(self, geometry: Any, check_spec: Dict[str, Any]) -> Dict[str, float]:
        """Compute bounding box dimensions.

        For assemblies, computes combined bounding box of all parts.
        """
        from yapcad.geom3d import solidbbox
        import math

        solids = self._normalize_geometry(geometry)
        if not solids:
            raise ValueError("No solid geometry to check")

        # Compute combined bounding box
        combined_bbox = None
        for solid in solids:
            bbox = solidbbox(solid)
            if bbox:
                if combined_bbox is None:
                    combined_bbox = bbox
                else:
                    # Merge bounding boxes
                    combined_bbox = [
                        [
                            min(combined_bbox[0][0], bbox[0][0]),  # xmin
                            min(combined_bbox[0][1], bbox[0][1]),  # ymin
                            min(combined_bbox[0][2], bbox[0][2]),  # zmin
                            1,
                        ],
                        [
                            max(combined_bbox[1][0], bbox[1][0]),  # xmax
                            max(combined_bbox[1][1], bbox[1][1]),  # ymax
                            max(combined_bbox[1][2], bbox[1][2]),  # zmax
                            1,
                        ],
                    ]

        if not combined_bbox:
            raise ValueError("Could not compute bounding box")

        # bbox format: [[xmin, ymin, zmin, 1], [xmax, ymax, zmax, 1]]
        xmin, ymin, zmin = combined_bbox[0][0], combined_bbox[0][1], combined_bbox[0][2]
        xmax, ymax, zmax = combined_bbox[1][0], combined_bbox[1][1], combined_bbox[1][2]

        width = xmax - xmin
        depth = ymax - ymin
        height = zmax - zmin

        diagonal = math.sqrt(width**2 + depth**2 + height**2)

        axis = check_spec.get("axis", "all")
        units = check_spec.get("units", "mm")

        metrics: Dict[str, float] = {
            "bbox.xmin": xmin,
            "bbox.xmax": xmax,
            "bbox.ymin": ymin,
            "bbox.ymax": ymax,
            "bbox.zmin": zmin,
            "bbox.zmax": zmax,
            "bbox.width": width,
            "bbox.depth": depth,
            "bbox.height": height,
            "bbox.diagonal": diagonal,
            "part_count": len(solids),
        }

        # Add axis-specific keys
        if axis == "x" or axis == "all":
            metrics[f"bbox.x_{units}"] = width
        if axis == "y" or axis == "all":
            metrics[f"bbox.y_{units}"] = depth
        if axis == "z" or axis == "all":
            metrics[f"bbox.z_{units}"] = height

        return metrics

    def _check_mass(self, geometry: Any, check_spec: Dict[str, Any]) -> Dict[str, float]:
        """Compute mass given density.

        For assemblies, computes total mass.
        """
        from yapcad.geom3d import volumeof

        solids = self._normalize_geometry(geometry)
        if not solids:
            raise ValueError("No solid geometry to check")

        total_volume = 0.0
        part_volumes = []

        for solid in solids:
            vol = volumeof(solid)
            if vol is None:
                raise ValueError("Could not compute volume for mass calculation")
            total_volume += vol
            part_volumes.append(vol)

        # Density in kg/m^3, volume in mm^3
        density_kgm3 = check_spec.get("density_kgm3", 2700)  # Default: aluminum

        # Convert mm^3 to m^3: 1 mm^3 = 1e-9 m^3
        volume_m3 = total_volume * 1e-9
        mass_kg = volume_m3 * density_kgm3

        metrics = {
            "volume_mm3": total_volume,
            "volume_m3": volume_m3,
            "mass_kg": mass_kg,
            "mass_g": mass_kg * 1000,
            "density_kgm3": density_kgm3,
            "part_count": len(solids),
        }

        # Include per-part masses for assemblies
        if len(solids) > 1:
            for i, vol in enumerate(part_volumes):
                vol_m3 = vol * 1e-9
                metrics[f"mass_kg.part_{i}"] = vol_m3 * density_kgm3

        return metrics

    def _check_centroid(self, geometry: Any, check_spec: Dict[str, Any]) -> Dict[str, float]:
        """Compute center of mass (centroid) of a solid or assembly.

        Uses BREP if available, falls back to mesh-based calculation.
        For assemblies, computes volume-weighted centroid.
        """
        from yapcad.geom3d import volumeof

        solids = self._normalize_geometry(geometry)
        if not solids:
            raise ValueError("No solid geometry to check")

        # For assembly, compute volume-weighted centroid
        total_volume = 0.0
        weighted_cx, weighted_cy, weighted_cz = 0.0, 0.0, 0.0
        part_centroids = []

        for solid in solids:
            cent = self._compute_solid_centroid(solid)
            if cent is None:
                raise ValueError("Could not compute centroid for one or more parts")

            vol = volumeof(solid)
            if vol is None:
                vol = 1.0  # Fallback to equal weighting

            total_volume += vol
            weighted_cx += vol * cent[0]
            weighted_cy += vol * cent[1]
            weighted_cz += vol * cent[2]
            part_centroids.append((cent, vol))

        if abs(total_volume) < 1e-10:
            raise ValueError("Total volume is zero")

        final_centroid = [
            weighted_cx / total_volume,
            weighted_cy / total_volume,
            weighted_cz / total_volume,
        ]

        metrics = {
            "centroid.x": final_centroid[0],
            "centroid.y": final_centroid[1],
            "centroid.z": final_centroid[2],
            "part_count": len(solids),
        }

        # Include per-part centroids for assemblies
        if len(solids) > 1:
            for i, (cent, _) in enumerate(part_centroids):
                metrics[f"centroid.part_{i}.x"] = cent[0]
                metrics[f"centroid.part_{i}.y"] = cent[1]
                metrics[f"centroid.part_{i}.z"] = cent[2]

        return metrics

    def _compute_solid_centroid(self, solid: Any) -> Optional[List[float]]:
        """Compute centroid of a single solid."""
        cent = None

        # Try BREP first
        try:
            from yapcad.brep import brep_from_solid
            from OCC.Core.GProp import GProp_GProps
            from OCC.Core.BRepGProp import brepgprop

            brep = brep_from_solid(solid)
            if brep is not None and brep.shape is not None:
                props = GProp_GProps()
                brepgprop.VolumeProperties(brep.shape, props)
                cog = props.CentreOfMass()
                cent = [cog.X(), cog.Y(), cog.Z()]
        except Exception:
            pass

        # Fall back to mesh-based centroid calculation
        if cent is None:
            cent = self._mesh_centroid(solid)

        return cent

    def _mesh_centroid(self, solid: Any) -> Optional[List[float]]:
        """Compute centroid from triangle mesh surfaces.

        Uses the formula for centroid of a closed triangular mesh:
        weighted average of triangle centroids by signed volume contribution.
        """
        surfaces = solid[1] if len(solid) > 1 else []
        if not surfaces:
            return None

        total_volume = 0.0
        cx, cy, cz = 0.0, 0.0, 0.0

        for surf in surfaces:
            vertices = surf[1]
            faces = surf[3]

            for face in faces:
                # Get triangle vertices
                v0 = vertices[face[0]]
                v1 = vertices[face[1]]
                v2 = vertices[face[2]]

                # Signed volume of tetrahedron formed with origin
                # V = (1/6) * (v0 . (v1 x v2))
                cross = [
                    v1[1] * v2[2] - v1[2] * v2[1],
                    v1[2] * v2[0] - v1[0] * v2[2],
                    v1[0] * v2[1] - v1[1] * v2[0],
                ]
                signed_vol = (v0[0] * cross[0] + v0[1] * cross[1] + v0[2] * cross[2]) / 6.0

                # Centroid of tetrahedron is at (v0 + v1 + v2) / 4
                # (origin is at 0,0,0 so contributes nothing)
                tet_cx = (v0[0] + v1[0] + v2[0]) / 4.0
                tet_cy = (v0[1] + v1[1] + v2[1]) / 4.0
                tet_cz = (v0[2] + v1[2] + v2[2]) / 4.0

                cx += signed_vol * tet_cx
                cy += signed_vol * tet_cy
                cz += signed_vol * tet_cz
                total_volume += signed_vol

        if abs(total_volume) < 1e-10:
            return None

        return [cx / total_volume, cy / total_volume, cz / total_volume]

    def _check_clearance(
        self, manifest: Any, plan: AnalysisPlan, check_spec: Dict[str, Any]
    ) -> Dict[str, float]:
        """Check clearance between two entities.

        Note: This is a placeholder - full implementation would require
        entity lookup and distance computation between surfaces.
        """
        raise NotImplementedError(
            "Clearance check requires entity lookup - not yet implemented"
        )

    def _evaluate_acceptance(
        self, metrics: Dict[str, float], acceptance: Dict[str, Any]
    ) -> str:
        """Evaluate acceptance criteria against computed metrics."""
        if not acceptance:
            return "passed" if metrics else "pending"

        for key, rule in acceptance.items():
            # Try exact key match first
            metric_value = metrics.get(key)

            # Try alternate key formats
            if metric_value is None:
                # Try without dots (e.g., "bbox_x" for "bbox.x")
                metric_value = metrics.get(key.replace(".", "_"))
            if metric_value is None:
                # Try with underscores as dots
                metric_value = metrics.get(key.replace("_", "."))

            if metric_value is None:
                # Metric not found - can't evaluate
                continue

            limit = rule.get("limit")
            if limit is None:
                continue

            limit = float(limit)
            comparison = rule.get("comparison", "<=")
            tolerance = rule.get("tolerance", 0.0)

            if comparison == "<=":
                if metric_value > limit:
                    return "failed"
            elif comparison == "<":
                if metric_value >= limit:
                    return "failed"
            elif comparison == ">=":
                if metric_value < limit:
                    return "failed"
            elif comparison == ">":
                if metric_value <= limit:
                    return "failed"
            elif comparison == "==":
                if abs(metric_value - limit) > tolerance:
                    return "failed"
            elif comparison == "~=":
                # Approximate equality with tolerance
                if abs(metric_value - limit) > tolerance:
                    return "failed"

        return "passed"

    def _get_failed_criteria(
        self, metrics: Dict[str, float], acceptance: Dict[str, Any]
    ) -> List[Dict[str, Any]]:
        """Get list of failed acceptance criteria."""
        failed = []

        for key, rule in acceptance.items():
            metric_value = metrics.get(key)
            if metric_value is None:
                metric_value = metrics.get(key.replace(".", "_"))
            if metric_value is None:
                continue

            limit = rule.get("limit")
            if limit is None:
                continue

            limit = float(limit)
            comparison = rule.get("comparison", "<=")

            violated = False
            if comparison == "<=" and metric_value > limit:
                violated = True
            elif comparison == "<" and metric_value >= limit:
                violated = True
            elif comparison == ">=" and metric_value < limit:
                violated = True
            elif comparison == ">" and metric_value <= limit:
                violated = True

            if violated:
                failed.append({
                    "criterion": key,
                    "value": metric_value,
                    "limit": limit,
                    "comparison": comparison,
                })

        return failed

    def _build_acceptance_results(
        self, metrics: Dict[str, float], acceptance: Dict[str, Any]
    ) -> Dict[str, Dict[str, Any]]:
        """Build detailed acceptance results for each criterion."""
        results = {}

        for key, rule in acceptance.items():
            metric_value = metrics.get(key)
            if metric_value is None:
                metric_value = metrics.get(key.replace(".", "_"))

            limit = rule.get("limit")
            comparison = rule.get("comparison", "<=")

            result = {
                "value": metric_value,
                "limit": limit,
                "comparison": comparison,
                "passed": None,
            }

            if metric_value is not None and limit is not None:
                limit = float(limit)
                if comparison == "<=":
                    result["passed"] = metric_value <= limit
                elif comparison == "<":
                    result["passed"] = metric_value < limit
                elif comparison == ">=":
                    result["passed"] = metric_value >= limit
                elif comparison == ">":
                    result["passed"] = metric_value > limit
                elif comparison == "==":
                    tolerance = rule.get("tolerance", 0.0)
                    result["passed"] = abs(metric_value - limit) <= tolerance
                elif comparison == "~=":
                    tolerance = rule.get("tolerance", 0.0)
                    result["passed"] = abs(metric_value - limit) <= tolerance

            results[key] = result

        return results


# Register the backend
register_backend("yapcad", YapCADNativeAdapter)
register_backend("yapcad-native", YapCADNativeAdapter)


__all__ = [
    "YapCADNativeAdapter",
]
