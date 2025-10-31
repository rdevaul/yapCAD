import json
from pathlib import Path
from typing import Any

import pytest

from yapcad.geom import point
from yapcad.geom3d import poly2surfaceXY, solid
from yapcad.geom3d_util import extrude
from yapcad.metadata import get_solid_metadata, set_layer
from yapcad.package import PackageManifest, create_package_from_entities, load_analysis_plan
from yapcad.package.analysis import available_backends
from yapcad.package.analysis.base import AnalysisAdapter, register_backend
from yapcad.package.analysis.cli import analyze_package


def _make_solid():
    poly = [
        point(0, 0),
        point(100, 0),
        point(100, 100),
        point(0, 100),
        point(0, 0),
    ]
    surf, _ = poly2surfaceXY(poly)
    sld = extrude(surf, distance=3.0)
    set_layer(get_solid_metadata(sld, create=True), "structure")
    return sld


def test_load_plan_remote_fields(tmp_path: Path):
    plan_path = tmp_path / "remote_plan.yaml"
    plan_path.write_text(
        """
id: remote-test
kind: multiphysics
backend: comsol
description: Remote COMSOL validation
geometry:
  source: geometry/primary.json
execution:
  mode: remote
  transport: ssh
  host: comsol-cluster.internal
  options:
    queue: premium
  env:
    COMSOL_LICENSE_FILE: 1717@license-host
loads:
  - id: pressure
    type: pressure
    magnitude_pa: 5.0e5
boundaryConditions:
  - id: fixed-base
    type: fixed
acceptance:
  stress.von_mises.max:
    limit: 2.0e8
    comparison: "<="
""",
        encoding="utf-8",
    )

    plan = load_analysis_plan(plan_path)
    assert plan.plan_id == "remote-test"
    assert plan.execution.mode == "remote"
    assert plan.execution.transport == "ssh"
    assert plan.execution.host == "comsol-cluster.internal"
    assert float(plan.loads[0]["magnitude_pa"]) == pytest.approx(5.0e5)
    assert float(plan.acceptance["stress.von_mises.max"]["limit"]) == pytest.approx(2.0e8)


def test_analyze_package_records_summary(tmp_path: Path):
    sld = _make_solid()
    pkg_root = tmp_path / "demo.ycpkg"
    manifest = create_package_from_entities(
        [sld],
        pkg_root,
        name="Demo",
        version="0.1.0",
    )
    manifest.save()

    plan_dir = pkg_root / "validation" / "plans"
    plan_dir.mkdir(parents=True, exist_ok=True)
    plan_path = plan_dir / "bulkhead_fea.yaml"
    plan_path.write_text(
        """
id: bulkhead-fea
kind: fea
backend: calculix
geometry:
  source: geometry/primary.json
loads:
  - id: axial
    type: pressure
    magnitude_pa: 1.0e5
acceptance:
  stress.von_mises.max:
    limit: 1.5e8
""",
        encoding="utf-8",
    )

    summary_path = analyze_package(pkg_root, plan_path.relative_to(pkg_root))
    assert summary_path.exists()

    summary = json.loads(summary_path.read_text())
    assert summary["plan"] == "bulkhead-fea"
    assert summary["backend"] == "calculix"
    assert "command" in summary.get("execution", {})

    reloaded = load_analysis_plan(plan_path)
    assert reloaded.plan_id == "bulkhead-fea"
    manifest = PackageManifest.load(pkg_root)
    validation = manifest.data.get("validation", {})
    result_entry = next((item for item in validation.get("results", []) if item.get("plan") == "bulkhead-fea"), None)
    assert result_entry is not None
    assert result_entry.get("backend") == "calculix"
    assert "summary" in result_entry
    assert "statusDetail" in result_entry["summary"]
    assert any(item.get("id") == "bulkhead-fea" for item in validation.get("plans", []))


def test_register_backend_records_adapter():
    class DummyAdapter(AnalysisAdapter):
        name = "dummy"

        def run(self, manifest, plan, workspace: Path, **kwargs: Any):  # type: ignore[override]
            raise NotImplementedError

    register_backend("dummy", DummyAdapter)
    assert "dummy" in available_backends()
