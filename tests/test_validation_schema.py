"""Tests for validation test schema implementation."""

import json
import tempfile
from pathlib import Path

import pytest
import yaml

from yapcad.package.analysis.schema import (
    ComparisonOp,
    ResultStatus,
    ValidationKind,
    ValidationReport,
    validate_plan,
    validate_plan_file,
    validate_result,
    validate_result_file,
)


class TestValidationKind:
    """Test ValidationKind enum."""

    def test_valid_kinds(self):
        assert ValidationKind.GEOMETRIC.value == "geometric"
        assert ValidationKind.MEASUREMENT.value == "measurement"
        assert ValidationKind.STRUCTURAL.value == "structural"
        assert ValidationKind.THERMAL.value == "thermal"
        assert ValidationKind.CFD.value == "cfd"
        assert ValidationKind.MULTIPHYSICS.value == "multiphysics"
        assert ValidationKind.ASSEMBLY.value == "assembly"


class TestResultStatus:
    """Test ResultStatus enum."""

    def test_valid_statuses(self):
        assert ResultStatus.PASSED.value == "passed"
        assert ResultStatus.FAILED.value == "failed"
        assert ResultStatus.ERROR.value == "error"
        assert ResultStatus.SKIPPED.value == "skipped"
        assert ResultStatus.PENDING.value == "pending"


class TestComparisonOp:
    """Test ComparisonOp enum."""

    def test_valid_operators(self):
        assert ComparisonOp.LE.value == "<="
        assert ComparisonOp.GE.value == ">="
        assert ComparisonOp.LT.value == "<"
        assert ComparisonOp.GT.value == ">"
        assert ComparisonOp.EQ.value == "=="
        assert ComparisonOp.APPROX.value == "~="


class TestValidatePlan:
    """Test validate_plan function."""

    def test_minimal_geometric_plan(self):
        """Test minimal valid geometric plan."""
        data = {
            "id": "test-001",
            "kind": "geometric",
            "backend": "yapcad",
            "check": {
                "property": "volume",
            },
        }
        report = validate_plan(data)
        assert report.valid

    def test_minimal_measurement_plan(self):
        """Test minimal valid measurement plan."""
        data = {
            "id": "mass-check",
            "kind": "measurement",
            "backend": "yapcad",
            "check": {
                "property": "mass",
                "density_kgm3": 2700,
            },
        }
        report = validate_plan(data)
        assert report.valid

    def test_minimal_structural_plan(self):
        """Test minimal valid structural plan."""
        data = {
            "id": "fea-001",
            "kind": "structural",
            "backend": "fenics",
        }
        report = validate_plan(data)
        assert report.valid
        # Should have warning about missing materials
        assert any("materials" in w.message for w in report.warnings)

    def test_missing_required_fields(self):
        """Test that missing required fields cause errors."""
        # Missing id
        data = {"kind": "geometric", "backend": "yapcad"}
        report = validate_plan(data)
        assert not report.valid
        assert any("id" in e.message for e in report.errors)

        # Missing kind
        data = {"id": "test", "backend": "yapcad"}
        report = validate_plan(data)
        assert not report.valid
        assert any("kind" in e.message for e in report.errors)

        # Missing backend
        data = {"id": "test", "kind": "geometric"}
        report = validate_plan(data)
        assert not report.valid
        assert any("backend" in e.message for e in report.errors)

    def test_invalid_kind(self):
        """Test that invalid kind values cause errors."""
        data = {
            "id": "test",
            "kind": "invalid_kind",
            "backend": "yapcad",
        }
        report = validate_plan(data)
        assert not report.valid
        assert any("kind" in e.path for e in report.errors)

    def test_geometric_plan_missing_check(self):
        """Test that geometric plans require check section."""
        data = {
            "id": "test",
            "kind": "geometric",
            "backend": "yapcad",
        }
        report = validate_plan(data)
        assert not report.valid
        assert any("check" in e.message for e in report.errors)

    def test_invalid_check_property(self):
        """Test that invalid check properties cause errors."""
        data = {
            "id": "test",
            "kind": "geometric",
            "backend": "yapcad",
            "check": {
                "property": "invalid_property",
            },
        }
        report = validate_plan(data)
        assert not report.valid
        assert any("property" in e.path for e in report.errors)

    def test_valid_acceptance_criteria(self):
        """Test valid acceptance criteria."""
        data = {
            "id": "test",
            "kind": "geometric",
            "backend": "yapcad",
            "check": {"property": "volume"},
            "acceptance": {
                "volume": {
                    "limit": 1000.0,
                    "comparison": ">=",
                    "description": "Minimum volume",
                },
            },
        }
        report = validate_plan(data)
        assert report.valid

    def test_invalid_comparison_operator(self):
        """Test that invalid comparison operators cause errors."""
        data = {
            "id": "test",
            "kind": "geometric",
            "backend": "yapcad",
            "check": {"property": "volume"},
            "acceptance": {
                "volume": {
                    "limit": 1000.0,
                    "comparison": "invalid",
                },
            },
        }
        report = validate_plan(data)
        assert not report.valid
        assert any("comparison" in e.path for e in report.errors)

    def test_approx_comparison_without_tolerance_warning(self):
        """Test that ~= without tolerance produces warning."""
        data = {
            "id": "test",
            "kind": "geometric",
            "backend": "yapcad",
            "check": {"property": "volume"},
            "acceptance": {
                "volume": {
                    "limit": 1000.0,
                    "comparison": "~=",
                },
            },
        }
        report = validate_plan(data)
        assert report.valid  # Still valid, just warning
        assert any("tolerance" in w.message for w in report.warnings)

    def test_structural_plan_with_loads(self):
        """Test structural plan with loads and BCs."""
        data = {
            "id": "fea",
            "kind": "structural",
            "backend": "fenics",
            "materials": {
                "default": {
                    "name": "Steel",
                    "youngs_modulus_pa": 200e9,
                    "poisson_ratio": 0.3,
                },
            },
            "loads": [
                {
                    "id": "pressure_1",
                    "type": "pressure",
                    "surfaces": ["top"],
                    "magnitude_pa": 1e6,
                },
            ],
            "boundaryConditions": [
                {
                    "id": "fixed_base",
                    "type": "fixed",
                    "strategy": "z_min",
                },
            ],
        }
        report = validate_plan(data)
        assert report.valid

    def test_invalid_load_type(self):
        """Test that invalid load types cause errors."""
        data = {
            "id": "fea",
            "kind": "structural",
            "backend": "fenics",
            "loads": [
                {
                    "type": "invalid_load_type",
                },
            ],
        }
        report = validate_plan(data)
        assert not report.valid
        assert any("load type" in e.message.lower() for e in report.errors)

    def test_invalid_bc_type(self):
        """Test that invalid BC types cause errors."""
        data = {
            "id": "fea",
            "kind": "structural",
            "backend": "fenics",
            "boundaryConditions": [
                {
                    "type": "invalid_bc",
                },
            ],
        }
        report = validate_plan(data)
        assert not report.valid
        assert any("BC type" in e.message for e in report.errors)

    def test_load_without_id_warning(self):
        """Test that loads without id produce warning."""
        data = {
            "id": "fea",
            "kind": "structural",
            "backend": "fenics",
            "loads": [
                {
                    "type": "pressure",
                    "surfaces": ["top"],
                },
            ],
        }
        report = validate_plan(data)
        assert report.valid
        assert any("id" in w.message.lower() for w in report.warnings)

    def test_complete_plan(self):
        """Test a complete structural plan."""
        data = {
            "id": "thrust-fea",
            "kind": "structural",
            "backend": "fenics",
            "name": "Thrust Plate Analysis",
            "description": "FEA analysis of thrust plate",
            "geometry": {
                "source": "geometry/primary.json",
                "entities": [],
            },
            "materials": {
                "default": {
                    "name": "6061-T6 Aluminum",
                    "youngs_modulus_pa": 68.9e9,
                    "poisson_ratio": 0.33,
                },
            },
            "loads": [
                {
                    "id": "thrust",
                    "type": "pressure",
                    "surfaces": ["motor_mount"],
                    "magnitude_pa": 2.61e6,
                },
            ],
            "boundaryConditions": [
                {
                    "id": "fixed",
                    "type": "fixed",
                    "strategy": "z_min",
                },
            ],
            "acceptance": {
                "displacement.max_mm": {
                    "limit": 1.0,
                    "comparison": "<=",
                },
            },
            "backendOptions": {
                "mesh": {"element_size": 5.0},
                "solver": {"type": "linear"},
            },
            "execution": {
                "mode": "local",
            },
            "metadata": {
                "author": "Test",
                "created": "2025-12-30",
            },
        }
        report = validate_plan(data)
        assert report.valid


class TestValidateResult:
    """Test validate_result function."""

    def test_minimal_result(self):
        """Test minimal valid result."""
        data = {
            "plan_id": "test-001",
            "status": "passed",
        }
        report = validate_result(data)
        assert report.valid

    def test_missing_required_fields(self):
        """Test that missing required fields cause errors."""
        # Missing plan_id
        data = {"status": "passed"}
        report = validate_result(data)
        assert not report.valid

        # Missing status
        data = {"plan_id": "test"}
        report = validate_result(data)
        assert not report.valid

    def test_invalid_status(self):
        """Test that invalid status values cause errors."""
        data = {
            "plan_id": "test",
            "status": "invalid_status",
        }
        report = validate_result(data)
        assert not report.valid
        assert any("status" in e.path for e in report.errors)

    def test_complete_result(self):
        """Test a complete result."""
        data = {
            "plan_id": "thrust-fea",
            "status": "passed",
            "timestamp": "2025-12-30T15:30:00Z",
            "backend": "fenics",
            "backend_version": "0.8.0",
            "execution_time_s": 45.2,
            "metrics": {
                "displacement.max_mm": 0.73,
                "stress.von_mises.max_pa": 145e6,
            },
            "acceptance_results": {
                "displacement.max_mm": {
                    "value": 0.73,
                    "limit": 1.0,
                    "comparison": "<=",
                    "passed": True,
                    "margin": 0.27,
                },
            },
            "artifacts": [
                {"kind": "mesh", "path": "mesh.msh"},
                {"kind": "solution", "path": "displacement.pvd"},
            ],
            "errors": [],
            "warnings": [],
            "notes": "Analysis completed successfully.",
        }
        report = validate_result(data)
        assert report.valid

    def test_acceptance_result_missing_value_warning(self):
        """Test that missing value in acceptance_results produces warning."""
        data = {
            "plan_id": "test",
            "status": "passed",
            "acceptance_results": {
                "metric": {
                    "limit": 1.0,
                    "passed": True,
                },
            },
        }
        report = validate_result(data)
        assert report.valid
        assert any("value" in w.message for w in report.warnings)

    def test_artifact_missing_path_warning(self):
        """Test that artifacts without path produce warning."""
        data = {
            "plan_id": "test",
            "status": "passed",
            "artifacts": [
                {"kind": "mesh"},
            ],
        }
        report = validate_result(data)
        assert report.valid
        assert any("path" in w.message for w in report.warnings)


class TestValidateFiles:
    """Test file validation functions."""

    def test_validate_plan_file(self, tmp_path):
        """Test validating a plan from a YAML file."""
        plan = {
            "id": "test",
            "kind": "geometric",
            "backend": "yapcad",
            "check": {"property": "volume"},
        }
        plan_file = tmp_path / "plan.yaml"
        with plan_file.open("w") as f:
            yaml.dump(plan, f)

        report = validate_plan_file(plan_file)
        assert report.valid

    def test_validate_plan_file_not_found(self):
        """Test validating a non-existent file."""
        report = validate_plan_file("/nonexistent/path/plan.yaml")
        assert not report.valid
        assert any("not found" in e.message.lower() for e in report.errors)

    def test_validate_plan_file_invalid_yaml(self, tmp_path):
        """Test validating invalid YAML."""
        plan_file = tmp_path / "bad.yaml"
        plan_file.write_text("{ invalid yaml")

        report = validate_plan_file(plan_file)
        assert not report.valid
        assert any("parse" in e.message.lower() for e in report.errors)

    def test_validate_result_file(self, tmp_path):
        """Test validating a result from a JSON file."""
        result = {
            "plan_id": "test",
            "status": "passed",
        }
        result_file = tmp_path / "result.json"
        with result_file.open("w") as f:
            json.dump(result, f)

        report = validate_result_file(result_file)
        assert report.valid

    def test_validate_result_file_not_found(self):
        """Test validating a non-existent result file."""
        report = validate_result_file("/nonexistent/path/result.json")
        assert not report.valid

    def test_validate_result_file_invalid_json(self, tmp_path):
        """Test validating invalid JSON."""
        result_file = tmp_path / "bad.json"
        result_file.write_text("{ invalid json")

        report = validate_result_file(result_file)
        assert not report.valid
        assert any("parse" in e.message.lower() for e in report.errors)


class TestValidationReport:
    """Test ValidationReport formatting."""

    def test_str_valid(self):
        """Test string representation of valid report."""
        report = ValidationReport(valid=True)
        assert "passed" in str(report).lower()

    def test_str_invalid(self):
        """Test string representation of invalid report."""
        from yapcad.package.analysis.schema import SchemaError

        report = ValidationReport(
            valid=False,
            errors=[SchemaError("id", "Missing required field")],
        )
        s = str(report)
        assert "failed" in s.lower()
        assert "Missing required field" in s


class TestRealPlanFiles:
    """Test validation against real plan files in the repository."""

    @pytest.mark.parametrize(
        "plan_path",
        [
            "examples/thrust_structure/validation/plans/thrust-fea.yaml",
            "examples/thrust_structure/validation/plans/mass-check.yaml",
        ],
    )
    def test_existing_plans_valid(self, plan_path):
        """Test that existing plan files in the repo are valid."""
        path = Path(plan_path)
        if not path.exists():
            pytest.skip(f"Plan file not found: {plan_path}")
        report = validate_plan_file(path)
        assert report.valid, f"Plan {plan_path} failed validation: {report}"
