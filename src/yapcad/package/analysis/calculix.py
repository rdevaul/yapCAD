"""CalculiX backend for yapCAD analysis plans."""

from __future__ import annotations

import math
import shutil
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from ..core import PackageManifest
from .base import AnalysisAdapter, AnalysisPlan, AnalysisResult, register_backend


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


@dataclass
class _MeshConfig:
    inner_radius: float
    outer_radius: float
    thickness: float
    radial_divisions: int
    thickness_divisions: int
    youngs_modulus: float
    poisson_ratio: float
    density: float
    thrust_force: float


class CalculixAdapter(AnalysisAdapter):
    """Create a simplified axisymmetric disk model and execute CalculiX when available."""

    name = "calculix"

    _DEFAULT_E = 68.9e3  # MPa for 6061-T6
    _DEFAULT_NU = 0.33
    _DEFAULT_DENSITY = 2.70e-6  # tonne/mm^3 (~2700 kg/m^3)

    def run(
        self,
        manifest: PackageManifest,
        plan: AnalysisPlan,
        workspace: Path,
        **_: Any,
    ) -> AnalysisResult:
        workspace.mkdir(parents=True, exist_ok=True)

        mesh = self._mesh_config(plan)
        inp_basename = plan.plan_id or "calculix_job"
        inp_path = workspace / f"{inp_basename}.inp"
        self._write_input_file(inp_path, mesh)

        command = plan.execution.command or "ccx"
        executable = shutil.which(command)

        summary: Dict[str, Any] = {
            "plan": plan.plan_id,
            "backend": plan.backend,
            "timestamp": _now(),
            "execution": {
                "mode": plan.execution.mode,
                "command": command,
            },
            "mesh": {
                "inner_radius_mm": mesh.inner_radius,
                "outer_radius_mm": mesh.outer_radius,
                "thickness_mm": mesh.thickness,
                "radial_divisions": mesh.radial_divisions,
                "thickness_divisions": mesh.thickness_divisions,
            },
        }

        artifacts = [
            {"kind": "solver-input", "path": inp_path.name},
        ]

        if executable is None:
            summary["statusDetail"] = f"Executable '{command}' not found on PATH"
            return AnalysisResult(
                plan_id=plan.plan_id,
                status="skipped",
                backend=plan.backend,
                summary=summary,
                artifacts=artifacts,
            )

        run = subprocess.run(
            [executable, inp_basename],
            cwd=workspace,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        summary.setdefault("solver", {})
        summary["solver"].update(
            {
                "returncode": run.returncode,
                "stdout_tail": run.stdout[-2000:],
                "stderr_tail": run.stderr[-2000:],
            }
        )

        if run.returncode != 0:
            summary["statusDetail"] = "CalculiX execution failed"
            return AnalysisResult(
                plan_id=plan.plan_id,
                status="error",
                backend=plan.backend,
                summary=summary,
                artifacts=artifacts,
            )

        dat_path = workspace / f"{inp_basename}.dat"
        max_disp = self._parse_max_displacement(dat_path)
        metrics: Dict[str, float] = {}
        if max_disp is not None:
            metrics["displacement.max"] = max_disp
            metrics["displacement.max_mm"] = max_disp

        status = self._evaluate_acceptance(metrics, plan.acceptance)
        summary["metrics"] = metrics
        if status == "passed":
            summary["statusDetail"] = "Acceptance criteria satisfied"
        elif status == "failed":
            summary["statusDetail"] = "Acceptance criteria violated"
        else:
            summary.setdefault("statusDetail", "Acceptance evaluation incomplete")

        for suffix in (".dat", ".frd", ".sta", ".log"):
            candidate = workspace / f"{inp_basename}{suffix}"
            if candidate.exists():
                artifacts.append({"kind": f"ccx{suffix}", "path": candidate.name})

        return AnalysisResult(
            plan_id=plan.plan_id,
            status=status,
            backend=plan.backend,
            metrics=metrics,
            summary=summary,
            artifacts=artifacts,
        )

    # ------------------------------------------------------------------
    # Mesh helpers

    def _mesh_config(self, plan: AnalysisPlan) -> _MeshConfig:
        opts = dict(plan.backend_options)
        geom = dict(plan.geometry)

        inner = float(opts.get("inner_radius_mm") or geom.get("inner_radius_mm") or 0.0)
        outer = float(opts.get("outer_radius_mm") or geom.get("outer_radius_mm") or (inner + 1.0))
        thickness = float(opts.get("thickness_mm") or geom.get("thickness_mm") or 10.0)
        thrust = float(opts.get("thrust_n") or opts.get("thrust_force_n") or 2224.0)

        radial = int(opts.get("radial_divisions", 24))
        axial = int(opts.get("thickness_divisions", 2))

        youngs = float(opts.get("youngs_modulus_mpa", self._DEFAULT_E))
        poisson = float(opts.get("poisson_ratio", self._DEFAULT_NU))
        density = float(opts.get("density_tonemm3", self._DEFAULT_DENSITY))

        return _MeshConfig(
            inner_radius=inner,
            outer_radius=outer,
            thickness=thickness,
            radial_divisions=radial,
            thickness_divisions=axial,
            youngs_modulus=youngs,
            poisson_ratio=poisson,
            density=density,
            thrust_force=thrust,
        )

    def _write_input_file(self, path: Path, mesh: _MeshConfig) -> None:
        nodes, elements, inner_nodes = self._generate_mesh(mesh)

        def fmt(num: float) -> str:
            return f"{num:.6g}"

        with path.open("w", encoding="utf-8") as fp:
            fp.write("*HEADING\n")
            fp.write("yapCAD bulkhead axisymmetric approximation\n")
            fp.write("*NODE\n")
            for nid, (r, z) in enumerate(nodes, start=1):
                fp.write(f"{nid}, {fmt(r)}, {fmt(z)}\n")

            fp.write("*ELEMENT, TYPE=CAX4\n")
            for eid, conn in enumerate(elements, start=1):
                fp.write(f"{eid}, {', '.join(str(nid) for nid in conn)}\n")

            fp.write("*ELSET, ELSET=EALL, GENERATE\n1, {0}\n".format(len(elements)))
            fp.write("*SOLID SECTION, ELSET=EALL, MATERIAL=AL6061\n1.0\n")

            fp.write("*MATERIAL, NAME=AL6061\n")
            fp.write("*ELASTIC\n")
            fp.write(f"{mesh.youngs_modulus}, {mesh.poisson_ratio}\n")
            fp.write("*DENSITY\n")
            fp.write(f"{mesh.density}\n")

            outer_nodes = [nid for nid, (r, _) in enumerate(nodes, start=1) if math.isclose(r, mesh.outer_radius, rel_tol=1e-6)]
            fp.write("*NSET, NSET=N_OUTER\n")
            fp.write(", ".join(str(n) for n in outer_nodes) + "\n")
            fp.write("*NSET, NSET=N_INNER\n")
            fp.write(", ".join(str(n) for n in inner_nodes) + "\n")

            fp.write("*BOUNDARY\n")
            fp.write("N_OUTER, 1, 2\n")

            fp.write("*STEP\n")
            fp.write("*STATIC\n")
            if inner_nodes:
                load_per_node = -(mesh.thrust_force / len(inner_nodes))
                fp.write("*CLOAD\n")
                for nid in inner_nodes:
                    fp.write(f"{nid}, 2, {fmt(load_per_node)}\n")
            fp.write("*NODE PRINT, NSET=N_INNER\n")
            fp.write("U\n")
            fp.write("*NODE PRINT, NSET=N_OUTER\n")
            fp.write("U\n")
            fp.write("*END STEP\n")

    def _generate_mesh(
        self, mesh: _MeshConfig
    ) -> Tuple[List[Tuple[float, float]], List[Tuple[int, int, int, int]], List[int]]:
        nodes: List[Tuple[float, float]] = []
        elements: List[Tuple[int, int, int, int]] = []

        r_vals = [mesh.inner_radius + (mesh.outer_radius - mesh.inner_radius) * i / mesh.radial_divisions for i in range(mesh.radial_divisions + 1)]
        z_vals = [mesh.thickness * j / mesh.thickness_divisions for j in range(mesh.thickness_divisions + 1)]

        node_index: Dict[Tuple[int, int], int] = {}
        nid = 1
        for j, z in enumerate(z_vals):
            for i, r in enumerate(r_vals):
                nodes.append((r, z))
                node_index[(i, j)] = nid
                nid += 1

        for j in range(mesh.thickness_divisions):
            for i in range(mesh.radial_divisions):
                n1 = node_index[(i, j)]
                n2 = node_index[(i + 1, j)]
                n3 = node_index[(i + 1, j + 1)]
                n4 = node_index[(i, j + 1)]
                elements.append((n1, n2, n3, n4))

        inner_nodes = [node_index[(0, j)] for j in range(mesh.thickness_divisions + 1)]
        return nodes, elements, inner_nodes

    def _parse_max_displacement(self, dat_path: Path) -> Optional[float]:
        if not dat_path.exists():
            return None
        max_disp: Optional[float] = None
        with dat_path.open("r", encoding="utf-8", errors="ignore") as fp:
            for line in fp:
                parts = line.strip().replace("=", " ").split()
                if not parts:
                    continue
                try:
                    int(parts[0])
                except ValueError:
                    continue
                numeric_parts: List[float] = []
                for token in parts[1:]:
                    try:
                        numeric_parts.append(float(token))
                    except ValueError:
                        continue
                if not numeric_parts:
                    continue
                uz = numeric_parts[1] if len(numeric_parts) > 1 else numeric_parts[0]
                uz_abs = abs(uz)
                if max_disp is None or uz_abs > max_disp:
                    max_disp = uz_abs
        return max_disp

    def _evaluate_acceptance(self, metrics: Dict[str, float], acceptance: Dict[str, Any]) -> str:
        if not acceptance:
            return "pending" if not metrics else "passed"

        status = "passed"
        for key, rule in acceptance.items():
            metric_value = metrics.get(key)
            if metric_value is None:
                metric_value = metrics.get(f"{key}_mm")
            if metric_value is None:
                status = "pending"
                continue
            limit = rule.get("limit")
            if limit is None and "limit_mm" in rule:
                limit = rule["limit_mm"]
            if limit is None:
                continue
            limit = float(limit)
            comparison = rule.get("comparison", "<=")
            if comparison == "<=" and not (metric_value <= limit):
                return "failed"
            if comparison == "<" and not (metric_value < limit):
                return "failed"
            if comparison == ">=" and not (metric_value >= limit):
                return "failed"
            if comparison == ">" and not (metric_value > limit):
                return "failed"
        return status


register_backend("calculix", CalculixAdapter)

__all__ = ["CalculixAdapter"]
