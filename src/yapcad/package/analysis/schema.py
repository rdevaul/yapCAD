"""Validation test schema implementation.

This module provides schema validation for yapCAD validation plans and results
as specified in ``docs/validation_schema.rst``.

Schema Version: validation-schema-v0.1
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Union


class ValidationKind(Enum):
    """Supported validation test kinds."""

    GEOMETRIC = "geometric"
    MEASUREMENT = "measurement"
    STRUCTURAL = "structural"
    THERMAL = "thermal"
    CFD = "cfd"
    MULTIPHYSICS = "multiphysics"
    ASSEMBLY = "assembly"


class ResultStatus(Enum):
    """Result status values."""

    PASSED = "passed"
    FAILED = "failed"
    ERROR = "error"
    SKIPPED = "skipped"
    PENDING = "pending"


class ComparisonOp(Enum):
    """Acceptance criteria comparison operators."""

    LE = "<="
    GE = ">="
    LT = "<"
    GT = ">"
    EQ = "=="
    APPROX = "~="


@dataclass
class SchemaError:
    """Represents a schema validation error."""

    path: str
    message: str
    severity: str = "error"  # "error" or "warning"

    def __str__(self) -> str:
        prefix = "WARNING" if self.severity == "warning" else "ERROR"
        return f"{prefix} at {self.path}: {self.message}"


@dataclass
class ValidationReport:
    """Result of schema validation."""

    valid: bool
    errors: List[SchemaError] = field(default_factory=list)
    warnings: List[SchemaError] = field(default_factory=list)
    schema_version: str = "validation-schema-v0.1"

    def __str__(self) -> str:
        if self.valid:
            return "Schema validation passed"
        lines = ["Schema validation failed:"]
        for err in self.errors:
            lines.append(f"  - {err}")
        if self.warnings:
            lines.append("Warnings:")
            for warn in self.warnings:
                lines.append(f"  - {warn}")
        return "\n".join(lines)


# -----------------------------------------------------------------------------
# Schema Definitions
# -----------------------------------------------------------------------------

# Required fields for all validation plans
PLAN_REQUIRED_FIELDS = {"id", "kind", "backend"}

# Optional common fields
PLAN_COMMON_FIELDS = {
    "name",
    "description",
    "geometry",
    "acceptance",
    "execution",
    "metadata",
}

# Kind-specific fields
KIND_SPECIFIC_FIELDS: Dict[str, set] = {
    "geometric": {"check"},
    "measurement": {"check"},
    "structural": {"materials", "loads", "boundaryConditions", "boundary_conditions", "backendOptions", "faceNaming"},
    "thermal": {"materials", "loads", "boundaryConditions", "boundary_conditions", "backendOptions", "faceNaming"},
    "cfd": {"materials", "boundaryConditions", "boundary_conditions", "backendOptions", "faceNaming"},
    "multiphysics": {"materials", "loads", "boundaryConditions", "boundary_conditions", "backendOptions", "faceNaming"},
    "assembly": {"check"},
}

# Valid check properties by kind
VALID_CHECK_PROPERTIES: Dict[str, set] = {
    "geometric": {"volume", "area", "bbox", "centroid", "distance", "clearance"},
    "measurement": {"dimension", "mass", "centroid", "moment"},
    "assembly": {"interference", "fit", "alignment", "clearance"},
}

# Valid load types
VALID_LOAD_TYPES = {"pressure", "force", "moment", "gravity", "thermal", "displacement"}

# Valid boundary condition types
VALID_BC_TYPES = {"fixed", "pinned", "roller", "spring", "symmetry", "thermal"}

# Result required fields
RESULT_REQUIRED_FIELDS = {"plan_id", "status"}

# Result optional fields
RESULT_OPTIONAL_FIELDS = {
    "timestamp",
    "backend",
    "backend_version",
    "execution_time_s",
    "metrics",
    "acceptance_results",
    "artifacts",
    "errors",
    "warnings",
    "notes",
}


# -----------------------------------------------------------------------------
# Validation Functions
# -----------------------------------------------------------------------------


def _check_required_fields(
    data: Dict[str, Any],
    required: set,
    path: str,
    errors: List[SchemaError],
) -> bool:
    """Check that all required fields are present."""
    missing = required - set(data.keys())
    if missing:
        errors.append(SchemaError(
            path=path,
            message=f"Missing required fields: {', '.join(sorted(missing))}",
        ))
        return False
    return True


def _check_field_type(
    data: Dict[str, Any],
    field: str,
    expected_type: type,
    path: str,
    errors: List[SchemaError],
) -> bool:
    """Check that a field has the expected type."""
    if field not in data:
        return True  # Missing fields are handled elsewhere
    value = data[field]
    if not isinstance(value, expected_type):
        errors.append(SchemaError(
            path=f"{path}.{field}",
            message=f"Expected {expected_type.__name__}, got {type(value).__name__}",
        ))
        return False
    return True


def _check_enum_value(
    value: Any,
    valid_values: Sequence[str],
    path: str,
    errors: List[SchemaError],
) -> bool:
    """Check that a value is one of the valid options."""
    if value not in valid_values:
        errors.append(SchemaError(
            path=path,
            message=f"Invalid value '{value}', expected one of: {', '.join(sorted(valid_values))}",
        ))
        return False
    return True


def validate_acceptance_criterion(
    criterion: Dict[str, Any],
    path: str,
    errors: List[SchemaError],
    warnings: List[SchemaError],
) -> bool:
    """Validate a single acceptance criterion."""
    valid = True

    # Check for limit
    if "limit" not in criterion:
        errors.append(SchemaError(
            path=path,
            message="Missing required field 'limit'",
        ))
        valid = False

    # Check comparison operator
    if "comparison" in criterion:
        comparison = criterion["comparison"]
        valid_ops = {op.value for op in ComparisonOp}
        if comparison not in valid_ops:
            errors.append(SchemaError(
                path=f"{path}.comparison",
                message=f"Invalid comparison '{comparison}', expected one of: {', '.join(sorted(valid_ops))}",
            ))
            valid = False

        # For approximate comparison, tolerance is recommended
        if comparison == "~=" and "tolerance" not in criterion:
            warnings.append(SchemaError(
                path=path,
                message="Approximate comparison (~=) used without tolerance",
                severity="warning",
            ))

    return valid


def validate_acceptance_criteria(
    acceptance: Dict[str, Any],
    path: str,
    errors: List[SchemaError],
    warnings: List[SchemaError],
) -> bool:
    """Validate acceptance criteria section."""
    if not isinstance(acceptance, dict):
        errors.append(SchemaError(
            path=path,
            message=f"Expected mapping, got {type(acceptance).__name__}",
        ))
        return False

    valid = True
    for metric_name, criterion in acceptance.items():
        if not isinstance(criterion, dict):
            errors.append(SchemaError(
                path=f"{path}.{metric_name}",
                message=f"Expected mapping, got {type(criterion).__name__}",
            ))
            valid = False
            continue
        valid = validate_acceptance_criterion(
            criterion, f"{path}.{metric_name}", errors, warnings
        ) and valid

    return valid


def validate_load(
    load: Dict[str, Any],
    index: int,
    path: str,
    errors: List[SchemaError],
    warnings: List[SchemaError],
) -> bool:
    """Validate a load definition."""
    valid = True
    load_path = f"{path}[{index}]"

    # Check required fields
    if "type" not in load:
        errors.append(SchemaError(
            path=load_path,
            message="Missing required field 'type'",
        ))
        valid = False
    else:
        load_type = load["type"]
        if load_type not in VALID_LOAD_TYPES:
            errors.append(SchemaError(
                path=f"{load_path}.type",
                message=f"Invalid load type '{load_type}', expected one of: {', '.join(sorted(VALID_LOAD_TYPES))}",
            ))
            valid = False

    # id is recommended
    if "id" not in load:
        warnings.append(SchemaError(
            path=load_path,
            message="Load missing 'id' field (recommended for clarity)",
            severity="warning",
        ))

    # Check for surfaces or strategy
    if "surfaces" not in load and "strategy" not in load:
        warnings.append(SchemaError(
            path=load_path,
            message="Load has no 'surfaces' or 'strategy' field (application region unclear)",
            severity="warning",
        ))

    return valid


def validate_boundary_condition(
    bc: Dict[str, Any],
    index: int,
    path: str,
    errors: List[SchemaError],
    warnings: List[SchemaError],
) -> bool:
    """Validate a boundary condition definition."""
    valid = True
    bc_path = f"{path}[{index}]"

    # Check required fields
    if "type" not in bc:
        errors.append(SchemaError(
            path=bc_path,
            message="Missing required field 'type'",
        ))
        valid = False
    else:
        bc_type = bc["type"]
        if bc_type not in VALID_BC_TYPES:
            errors.append(SchemaError(
                path=f"{bc_path}.type",
                message=f"Invalid BC type '{bc_type}', expected one of: {', '.join(sorted(VALID_BC_TYPES))}",
            ))
            valid = False

    # id is recommended
    if "id" not in bc:
        warnings.append(SchemaError(
            path=bc_path,
            message="Boundary condition missing 'id' field (recommended for clarity)",
            severity="warning",
        ))

    # Check for surfaces or strategy
    if "surfaces" not in bc and "strategy" not in bc:
        warnings.append(SchemaError(
            path=bc_path,
            message="BC has no 'surfaces' or 'strategy' field (application region unclear)",
            severity="warning",
        ))

    return valid


def validate_geometry_section(
    geometry: Dict[str, Any],
    path: str,
    errors: List[SchemaError],
    warnings: List[SchemaError],
) -> bool:
    """Validate the geometry section of a plan."""
    valid = True

    if not isinstance(geometry, dict):
        errors.append(SchemaError(
            path=path,
            message=f"Expected mapping, got {type(geometry).__name__}",
        ))
        return False

    # source is effectively required for most plans
    if "source" not in geometry:
        warnings.append(SchemaError(
            path=path,
            message="No 'source' field - geometry reference unclear",
            severity="warning",
        ))

    # entities should be a list
    if "entities" in geometry:
        entities = geometry["entities"]
        if not isinstance(entities, list):
            errors.append(SchemaError(
                path=f"{path}.entities",
                message=f"Expected list, got {type(entities).__name__}",
            ))
            valid = False

    return valid


def validate_check_section(
    check: Dict[str, Any],
    kind: str,
    path: str,
    errors: List[SchemaError],
    warnings: List[SchemaError],
) -> bool:
    """Validate the check section for geometric/measurement/assembly tests."""
    valid = True

    if not isinstance(check, dict):
        errors.append(SchemaError(
            path=path,
            message=f"Expected mapping, got {type(check).__name__}",
        ))
        return False

    # property is required
    if "property" not in check:
        errors.append(SchemaError(
            path=path,
            message="Missing required field 'property'",
        ))
        valid = False
    else:
        prop = check["property"]
        valid_props = VALID_CHECK_PROPERTIES.get(kind, set())
        if prop not in valid_props:
            errors.append(SchemaError(
                path=f"{path}.property",
                message=f"Invalid property '{prop}' for kind '{kind}', expected one of: {', '.join(sorted(valid_props))}",
            ))
            valid = False

    return valid


def validate_plan(data: Dict[str, Any]) -> ValidationReport:
    """Validate a validation plan against the schema.

    Args:
        data: The plan data as a dictionary (loaded from YAML)

    Returns:
        ValidationReport with validation results
    """
    errors: List[SchemaError] = []
    warnings: List[SchemaError] = []

    # Check required fields
    _check_required_fields(data, PLAN_REQUIRED_FIELDS, "", errors)

    if errors:
        # Can't continue without required fields
        return ValidationReport(valid=False, errors=errors, warnings=warnings)

    # Validate kind
    kind = data["kind"]
    valid_kinds = {k.value for k in ValidationKind}
    if kind not in valid_kinds:
        errors.append(SchemaError(
            path="kind",
            message=f"Invalid kind '{kind}', expected one of: {', '.join(sorted(valid_kinds))}",
        ))
        return ValidationReport(valid=False, errors=errors, warnings=warnings)

    # Validate id is a non-empty string
    plan_id = data["id"]
    if not isinstance(plan_id, str) or not plan_id.strip():
        errors.append(SchemaError(
            path="id",
            message="Plan 'id' must be a non-empty string",
        ))

    # Validate backend is a non-empty string
    backend = data["backend"]
    if not isinstance(backend, str) or not backend.strip():
        errors.append(SchemaError(
            path="backend",
            message="Plan 'backend' must be a non-empty string",
        ))

    # Validate geometry section if present
    if "geometry" in data:
        validate_geometry_section(data["geometry"], "geometry", errors, warnings)

    # Validate acceptance criteria if present
    if "acceptance" in data:
        validate_acceptance_criteria(data["acceptance"], "acceptance", errors, warnings)

    # Kind-specific validation
    if kind in ("geometric", "measurement", "assembly"):
        # These kinds require a 'check' section
        if "check" not in data:
            errors.append(SchemaError(
                path="",
                message=f"Kind '{kind}' requires a 'check' section",
            ))
        else:
            validate_check_section(data["check"], kind, "check", errors, warnings)

    elif kind in ("structural", "thermal", "cfd", "multiphysics"):
        # Simulation kinds - validate materials, loads, BCs
        if "materials" in data:
            materials = data["materials"]
            if not isinstance(materials, dict):
                errors.append(SchemaError(
                    path="materials",
                    message=f"Expected mapping, got {type(materials).__name__}",
                ))

        # Validate loads
        if "loads" in data:
            loads = data["loads"]
            if not isinstance(loads, list):
                errors.append(SchemaError(
                    path="loads",
                    message=f"Expected list, got {type(loads).__name__}",
                ))
            else:
                for i, load in enumerate(loads):
                    if isinstance(load, dict):
                        validate_load(load, i, "loads", errors, warnings)
                    else:
                        errors.append(SchemaError(
                            path=f"loads[{i}]",
                            message=f"Expected mapping, got {type(load).__name__}",
                        ))

        # Validate boundary conditions
        bc_key = "boundaryConditions"
        if bc_key not in data:
            bc_key = "boundary_conditions"

        if bc_key in data:
            bcs = data[bc_key]
            if not isinstance(bcs, list):
                errors.append(SchemaError(
                    path=bc_key,
                    message=f"Expected list, got {type(bcs).__name__}",
                ))
            else:
                for i, bc in enumerate(bcs):
                    if isinstance(bc, dict):
                        validate_boundary_condition(bc, i, bc_key, errors, warnings)
                    else:
                        errors.append(SchemaError(
                            path=f"{bc_key}[{i}]",
                            message=f"Expected mapping, got {type(bc).__name__}",
                        ))

        # Warn if no materials defined for structural/thermal
        if kind in ("structural", "thermal") and "materials" not in data:
            warnings.append(SchemaError(
                path="",
                message=f"Kind '{kind}' typically requires 'materials' section",
                severity="warning",
            ))

    # Check for unknown top-level fields
    known_fields = PLAN_REQUIRED_FIELDS | PLAN_COMMON_FIELDS
    if kind in KIND_SPECIFIC_FIELDS:
        known_fields = known_fields | KIND_SPECIFIC_FIELDS[kind]
    known_fields.add("attachments")  # Always allowed

    unknown = set(data.keys()) - known_fields
    # Don't warn about metadata subfields or common extras
    for unk in sorted(unknown):
        warnings.append(SchemaError(
            path=unk,
            message=f"Unknown field '{unk}' (may be ignored)",
            severity="warning",
        ))

    return ValidationReport(
        valid=len(errors) == 0,
        errors=errors,
        warnings=warnings,
    )


def validate_result(data: Dict[str, Any]) -> ValidationReport:
    """Validate a validation result against the schema.

    Args:
        data: The result data as a dictionary (loaded from JSON)

    Returns:
        ValidationReport with validation results
    """
    errors: List[SchemaError] = []
    warnings: List[SchemaError] = []

    # Check required fields
    _check_required_fields(data, RESULT_REQUIRED_FIELDS, "", errors)

    if errors:
        return ValidationReport(valid=False, errors=errors, warnings=warnings)

    # Validate status
    status = data["status"]
    valid_statuses = {s.value for s in ResultStatus}
    if status not in valid_statuses:
        errors.append(SchemaError(
            path="status",
            message=f"Invalid status '{status}', expected one of: {', '.join(sorted(valid_statuses))}",
        ))

    # Validate plan_id
    plan_id = data["plan_id"]
    if not isinstance(plan_id, str) or not plan_id.strip():
        errors.append(SchemaError(
            path="plan_id",
            message="Result 'plan_id' must be a non-empty string",
        ))

    # Validate metrics if present
    if "metrics" in data:
        metrics = data["metrics"]
        if not isinstance(metrics, dict):
            errors.append(SchemaError(
                path="metrics",
                message=f"Expected mapping, got {type(metrics).__name__}",
            ))

    # Validate acceptance_results if present
    if "acceptance_results" in data:
        ar = data["acceptance_results"]
        if not isinstance(ar, dict):
            errors.append(SchemaError(
                path="acceptance_results",
                message=f"Expected mapping, got {type(ar).__name__}",
            ))
        else:
            for metric_name, result in ar.items():
                if not isinstance(result, dict):
                    errors.append(SchemaError(
                        path=f"acceptance_results.{metric_name}",
                        message=f"Expected mapping, got {type(result).__name__}",
                    ))
                    continue

                # Check required fields in each acceptance result
                if "value" not in result:
                    warnings.append(SchemaError(
                        path=f"acceptance_results.{metric_name}",
                        message="Missing 'value' field",
                        severity="warning",
                    ))
                if "passed" not in result:
                    warnings.append(SchemaError(
                        path=f"acceptance_results.{metric_name}",
                        message="Missing 'passed' field",
                        severity="warning",
                    ))

    # Validate artifacts if present
    if "artifacts" in data:
        artifacts = data["artifacts"]
        if not isinstance(artifacts, list):
            errors.append(SchemaError(
                path="artifacts",
                message=f"Expected list, got {type(artifacts).__name__}",
            ))
        else:
            for i, artifact in enumerate(artifacts):
                if not isinstance(artifact, dict):
                    errors.append(SchemaError(
                        path=f"artifacts[{i}]",
                        message=f"Expected mapping, got {type(artifact).__name__}",
                    ))
                elif "path" not in artifact:
                    warnings.append(SchemaError(
                        path=f"artifacts[{i}]",
                        message="Artifact missing 'path' field",
                        severity="warning",
                    ))

    # Check for unknown fields
    known_fields = RESULT_REQUIRED_FIELDS | RESULT_OPTIONAL_FIELDS
    unknown = set(data.keys()) - known_fields
    for unk in sorted(unknown):
        warnings.append(SchemaError(
            path=unk,
            message=f"Unknown field '{unk}' (may be ignored)",
            severity="warning",
        ))

    return ValidationReport(
        valid=len(errors) == 0,
        errors=errors,
        warnings=warnings,
    )


def validate_plan_file(path: Union[str, Path]) -> ValidationReport:
    """Validate a validation plan YAML file.

    Args:
        path: Path to the YAML file

    Returns:
        ValidationReport with validation results
    """
    import yaml

    plan_path = Path(path)
    if not plan_path.exists():
        return ValidationReport(
            valid=False,
            errors=[SchemaError("", f"File not found: {plan_path}")],
        )

    try:
        with plan_path.open("r", encoding="utf-8") as fp:
            data = yaml.safe_load(fp)
    except yaml.YAMLError as e:
        return ValidationReport(
            valid=False,
            errors=[SchemaError("", f"YAML parse error: {e}")],
        )

    if data is None:
        return ValidationReport(
            valid=False,
            errors=[SchemaError("", "Empty YAML file")],
        )

    if not isinstance(data, dict):
        return ValidationReport(
            valid=False,
            errors=[SchemaError("", f"Expected mapping at root, got {type(data).__name__}")],
        )

    return validate_plan(data)


def validate_result_file(path: Union[str, Path]) -> ValidationReport:
    """Validate a validation result JSON file.

    Args:
        path: Path to the JSON file

    Returns:
        ValidationReport with validation results
    """
    import json

    result_path = Path(path)
    if not result_path.exists():
        return ValidationReport(
            valid=False,
            errors=[SchemaError("", f"File not found: {result_path}")],
        )

    try:
        with result_path.open("r", encoding="utf-8") as fp:
            data = json.load(fp)
    except json.JSONDecodeError as e:
        return ValidationReport(
            valid=False,
            errors=[SchemaError("", f"JSON parse error: {e}")],
        )

    if not isinstance(data, dict):
        return ValidationReport(
            valid=False,
            errors=[SchemaError("", f"Expected mapping at root, got {type(data).__name__}")],
        )

    return validate_result(data)


__all__ = [
    "ValidationKind",
    "ResultStatus",
    "ComparisonOp",
    "SchemaError",
    "ValidationReport",
    "validate_plan",
    "validate_result",
    "validate_plan_file",
    "validate_result_file",
]
