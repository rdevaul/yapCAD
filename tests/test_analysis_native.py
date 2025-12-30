"""Tests for the native yapCAD validation backend.

Tests geometric and measurement validation checks using yapCAD's built-in
functions without external solvers.
"""

import json
from pathlib import Path

import pytest

from yapcad.geom import point
from yapcad.geom3d import issolid, poly2surfaceXY
from yapcad.geom3d_util import extrude
from yapcad.metadata import get_solid_metadata, set_layer
from yapcad.package import create_package_from_entities
from yapcad.package.analysis import available_backends, get_backend
from yapcad.package.analysis.base import load_plan
from yapcad.package.analysis.yapcad_native import YapCADNativeAdapter


def _make_test_solid():
    """Create a simple test solid (100x100x10 box)."""
    poly = [
        point(0, 0),
        point(100, 0),
        point(100, 100),
        point(0, 100),
        point(0, 0),
    ]
    surf, _ = poly2surfaceXY(poly)
    sld = extrude(surf, distance=10.0)
    set_layer(get_solid_metadata(sld, create=True), "test")
    return sld


def _create_test_package(tmp_path: Path, solid):
    """Create a test package with the given solid."""
    pkg_root = tmp_path / "test.ycpkg"
    manifest = create_package_from_entities(
        [solid],
        pkg_root,
        name="Test Package",
        version="0.1.0",
    )
    manifest.save()
    return pkg_root, manifest


class TestYapCADNativeBackend:
    """Tests for the native yapCAD validation backend."""

    def test_backend_registered(self):
        """Test that yapcad backend is registered."""
        backends = available_backends()
        assert "yapcad" in backends
        assert "yapcad-native" in backends

    def test_get_backend(self):
        """Test getting the backend by name."""
        adapter_cls = get_backend("yapcad")
        assert adapter_cls is YapCADNativeAdapter

    def test_volume_check_passes(self, tmp_path: Path):
        """Test volume check that should pass."""
        solid = _make_test_solid()
        pkg_root, manifest = _create_test_package(tmp_path, solid)

        # Create validation plan
        plan_dir = pkg_root / "validation" / "plans"
        plan_dir.mkdir(parents=True, exist_ok=True)
        plan_path = plan_dir / "volume-check.yaml"
        plan_path.write_text(
            """
id: volume-check
kind: geometric
backend: yapcad
check:
  property: volume
  units: mm3
acceptance:
  volume:
    limit: 50000.0
    comparison: ">="
    description: "Minimum volume requirement"
""",
            encoding="utf-8",
        )

        plan = load_plan(plan_path)
        adapter = YapCADNativeAdapter()
        workspace = pkg_root / "validation" / "results" / plan.plan_id
        result = adapter.run(manifest, plan, workspace)

        assert result.status == "passed"
        assert "volume" in result.metrics
        # 100 x 100 x 10 = 100000 mm^3
        assert result.metrics["volume"] == pytest.approx(100000.0, rel=0.01)

    def test_volume_check_fails(self, tmp_path: Path):
        """Test volume check that should fail."""
        solid = _make_test_solid()
        pkg_root, manifest = _create_test_package(tmp_path, solid)

        plan_dir = pkg_root / "validation" / "plans"
        plan_dir.mkdir(parents=True, exist_ok=True)
        plan_path = plan_dir / "volume-check.yaml"
        plan_path.write_text(
            """
id: volume-check-fail
kind: geometric
backend: yapcad
check:
  property: volume
  units: mm3
acceptance:
  volume:
    limit: 200000.0
    comparison: ">="
    description: "Volume too high"
""",
            encoding="utf-8",
        )

        plan = load_plan(plan_path)
        adapter = YapCADNativeAdapter()
        workspace = pkg_root / "validation" / "results" / plan.plan_id
        result = adapter.run(manifest, plan, workspace)

        assert result.status == "failed"
        assert result.summary.get("failed_criteria")

    def test_bbox_check(self, tmp_path: Path):
        """Test bounding box check."""
        solid = _make_test_solid()
        pkg_root, manifest = _create_test_package(tmp_path, solid)

        plan_dir = pkg_root / "validation" / "plans"
        plan_dir.mkdir(parents=True, exist_ok=True)
        plan_path = plan_dir / "bbox-check.yaml"
        plan_path.write_text(
            """
id: bbox-check
kind: geometric
backend: yapcad
check:
  property: bbox
  axis: all
acceptance:
  bbox.width:
    limit: 150.0
    comparison: "<="
  bbox.height:
    limit: 20.0
    comparison: "<="
""",
            encoding="utf-8",
        )

        plan = load_plan(plan_path)
        adapter = YapCADNativeAdapter()
        workspace = pkg_root / "validation" / "results" / plan.plan_id
        result = adapter.run(manifest, plan, workspace)

        assert result.status == "passed"
        assert "bbox.width" in result.metrics
        assert "bbox.height" in result.metrics
        assert result.metrics["bbox.width"] == pytest.approx(100.0, rel=0.01)
        assert result.metrics["bbox.height"] == pytest.approx(10.0, rel=0.01)

    def test_mass_check(self, tmp_path: Path):
        """Test mass calculation with density."""
        solid = _make_test_solid()
        pkg_root, manifest = _create_test_package(tmp_path, solid)

        plan_dir = pkg_root / "validation" / "plans"
        plan_dir.mkdir(parents=True, exist_ok=True)
        plan_path = plan_dir / "mass-check.yaml"
        plan_path.write_text(
            """
id: mass-check
kind: measurement
backend: yapcad
check:
  property: mass
  density_kgm3: 2700
acceptance:
  mass_kg:
    limit: 1.0
    comparison: "<="
    description: "Max mass budget"
""",
            encoding="utf-8",
        )

        plan = load_plan(plan_path)
        adapter = YapCADNativeAdapter()
        workspace = pkg_root / "validation" / "results" / plan.plan_id
        result = adapter.run(manifest, plan, workspace)

        assert result.status == "passed"
        assert "mass_kg" in result.metrics
        # Volume = 100000 mm^3 = 1e-4 m^3
        # Mass = 1e-4 * 2700 = 0.27 kg
        assert result.metrics["mass_kg"] == pytest.approx(0.27, rel=0.01)

    def test_summary_json_written(self, tmp_path: Path):
        """Test that summary.json is written to workspace."""
        solid = _make_test_solid()
        pkg_root, manifest = _create_test_package(tmp_path, solid)

        plan_dir = pkg_root / "validation" / "plans"
        plan_dir.mkdir(parents=True, exist_ok=True)
        plan_path = plan_dir / "volume-check.yaml"
        plan_path.write_text(
            """
id: summary-test
kind: geometric
backend: yapcad
check:
  property: volume
acceptance:
  volume:
    limit: 50000.0
    comparison: ">="
""",
            encoding="utf-8",
        )

        plan = load_plan(plan_path)
        adapter = YapCADNativeAdapter()
        workspace = pkg_root / "validation" / "results" / plan.plan_id
        result = adapter.run(manifest, plan, workspace)

        summary_path = workspace / "summary.json"
        assert summary_path.exists()

        summary = json.loads(summary_path.read_text())
        assert summary["plan_id"] == "summary-test"
        assert summary["status"] == "passed"
        assert "metrics" in summary
        assert "acceptance_results" in summary

    def test_acceptance_results_detail(self, tmp_path: Path):
        """Test detailed acceptance results in summary."""
        solid = _make_test_solid()
        pkg_root, manifest = _create_test_package(tmp_path, solid)

        plan_dir = pkg_root / "validation" / "plans"
        plan_dir.mkdir(parents=True, exist_ok=True)
        plan_path = plan_dir / "detail-check.yaml"
        plan_path.write_text(
            """
id: detail-test
kind: geometric
backend: yapcad
check:
  property: bbox
acceptance:
  bbox.width:
    limit: 150.0
    comparison: "<="
  bbox.diagonal:
    limit: 200.0
    comparison: "<="
""",
            encoding="utf-8",
        )

        plan = load_plan(plan_path)
        adapter = YapCADNativeAdapter()
        workspace = pkg_root / "validation" / "results" / plan.plan_id
        result = adapter.run(manifest, plan, workspace)

        summary_path = workspace / "summary.json"
        summary = json.loads(summary_path.read_text())

        acceptance = summary["acceptance_results"]
        assert "bbox.width" in acceptance
        assert acceptance["bbox.width"]["passed"] is True
        assert acceptance["bbox.width"]["value"] == pytest.approx(100.0, rel=0.01)
        assert acceptance["bbox.width"]["limit"] == 150.0


class TestAcceptanceCriteria:
    """Tests for acceptance criteria evaluation."""

    def test_less_than_or_equal(self, tmp_path: Path):
        """Test <= comparison."""
        solid = _make_test_solid()
        pkg_root, manifest = _create_test_package(tmp_path, solid)

        plan_dir = pkg_root / "validation" / "plans"
        plan_dir.mkdir(parents=True, exist_ok=True)

        # Test passing case - use generous limit since volume may vary slightly
        plan_path = plan_dir / "test.yaml"
        plan_path.write_text(
            """
id: le-pass
kind: geometric
backend: yapcad
check:
  property: volume
acceptance:
  volume:
    limit: 101000.0
    comparison: "<="
""",
            encoding="utf-8",
        )

        plan = load_plan(plan_path)
        adapter = YapCADNativeAdapter()
        workspace = pkg_root / "validation" / "results" / plan.plan_id
        result = adapter.run(manifest, plan, workspace)
        assert result.status == "passed"

    def test_greater_than_or_equal(self, tmp_path: Path):
        """Test >= comparison."""
        solid = _make_test_solid()
        pkg_root, manifest = _create_test_package(tmp_path, solid)

        plan_dir = pkg_root / "validation" / "plans"
        plan_dir.mkdir(parents=True, exist_ok=True)

        plan_path = plan_dir / "test.yaml"
        plan_path.write_text(
            """
id: ge-pass
kind: geometric
backend: yapcad
check:
  property: volume
acceptance:
  volume:
    limit: 100000.0
    comparison: ">="
""",
            encoding="utf-8",
        )

        plan = load_plan(plan_path)
        adapter = YapCADNativeAdapter()
        workspace = pkg_root / "validation" / "results" / plan.plan_id
        result = adapter.run(manifest, plan, workspace)
        assert result.status == "passed"

    def test_approximate_equality(self, tmp_path: Path):
        """Test ~= approximate comparison."""
        solid = _make_test_solid()
        pkg_root, manifest = _create_test_package(tmp_path, solid)

        plan_dir = pkg_root / "validation" / "plans"
        plan_dir.mkdir(parents=True, exist_ok=True)

        plan_path = plan_dir / "test.yaml"
        plan_path.write_text(
            """
id: approx-pass
kind: geometric
backend: yapcad
check:
  property: volume
acceptance:
  volume:
    limit: 100500.0
    comparison: "~="
    tolerance: 1000.0
""",
            encoding="utf-8",
        )

        plan = load_plan(plan_path)
        adapter = YapCADNativeAdapter()
        workspace = pkg_root / "validation" / "results" / plan.plan_id
        result = adapter.run(manifest, plan, workspace)
        assert result.status == "passed"
