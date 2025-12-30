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

        Returns the first solid from the package geometry.
        """
        from yapcad.package.core import load_geometry
        from yapcad.geom3d import issolid

        # Get geometry source from plan
        geom_spec = plan.geometry
        entities = geom_spec.get("entities", [])

        # Load geometry from package
        geometry = load_geometry(manifest)

        if entities:
            # Filter to specified entities (future: implement entity lookup)
            pass

        # Extract the first solid from the geometry list
        if isinstance(geometry, list):
            for entity in geometry:
                if issolid(entity):
                    return entity
            # Check if geometry itself is a solid
            if geometry and issolid(geometry[0] if len(geometry) == 1 else geometry):
                return geometry[0] if len(geometry) == 1 else geometry
            raise ValueError("No solid geometry found in package")

        return geometry

    def _check_volume(self, solid: Any, check_spec: Dict[str, Any]) -> Dict[str, float]:
        """Compute volume of a solid."""
        from yapcad.geom3d import volumeof

        volume = volumeof(solid)
        if volume is None:
            raise ValueError("Could not compute volume - solid may lack BREP data")

        units = check_spec.get("units", "mm3")

        # Store volume with unit suffix for clarity
        return {
            "volume": volume,
            f"volume_{units}": volume,
        }

    def _check_area(self, solid: Any, check_spec: Dict[str, Any]) -> Dict[str, float]:
        """Compute surface area of a solid.

        Uses BREP if available, falls back to tessellation.
        """
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
            from yapcad.geom3d import surfacearea as surface_area_tess

            # For solids, we need to compute area of all surfaces
            surfaces = solid[1] if len(solid) > 1 else []
            if surfaces:
                area = sum(surface_area_tess(surf) for surf in surfaces)

        if area is None:
            raise ValueError("Could not compute surface area")

        units = check_spec.get("units", "mm2")

        return {
            "area": area,
            f"area_{units}": area,
        }

    def _check_bbox(self, solid: Any, check_spec: Dict[str, Any]) -> Dict[str, float]:
        """Compute bounding box dimensions."""
        from yapcad.geom3d import solidbbox

        bbox = solidbbox(solid)
        if not bbox:
            raise ValueError("Could not compute bounding box")

        # bbox format: [[xmin, ymin, zmin, 1], [xmax, ymax, zmax, 1]]
        xmin, ymin, zmin = bbox[0][0], bbox[0][1], bbox[0][2]
        xmax, ymax, zmax = bbox[1][0], bbox[1][1], bbox[1][2]

        width = xmax - xmin
        depth = ymax - ymin
        height = zmax - zmin

        import math
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
        }

        # Add axis-specific keys
        if axis == "x" or axis == "all":
            metrics[f"bbox.x_{units}"] = width
        if axis == "y" or axis == "all":
            metrics[f"bbox.y_{units}"] = depth
        if axis == "z" or axis == "all":
            metrics[f"bbox.z_{units}"] = height

        return metrics

    def _check_mass(self, solid: Any, check_spec: Dict[str, Any]) -> Dict[str, float]:
        """Compute mass given density."""
        from yapcad.geom3d import volumeof

        volume = volumeof(solid)
        if volume is None:
            raise ValueError("Could not compute volume for mass calculation")

        # Density in kg/m^3, volume in mm^3
        density_kgm3 = check_spec.get("density_kgm3", 2700)  # Default: aluminum

        # Convert mm^3 to m^3: 1 mm^3 = 1e-9 m^3
        volume_m3 = volume * 1e-9
        mass_kg = volume_m3 * density_kgm3

        return {
            "volume_mm3": volume,
            "volume_m3": volume_m3,
            "mass_kg": mass_kg,
            "mass_g": mass_kg * 1000,
            "density_kgm3": density_kgm3,
        }

    def _check_centroid(self, solid: Any, check_spec: Dict[str, Any]) -> Dict[str, float]:
        """Compute center of mass (centroid) of a solid.

        Uses BREP if available.
        """
        cent = None

        # Try BREP
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

        if cent is None:
            raise ValueError("Could not compute centroid - BREP data required")

        return {
            "centroid.x": cent[0],
            "centroid.y": cent[1],
            "centroid.z": cent[2],
        }

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
