"""Public API for yapCAD `.ycpkg` package workflows."""

from __future__ import annotations

from .core import (
    MANIFEST_FILENAME,
    PACKAGE_SCHEMA,
    PackageManifest,
    create_package_from_entities,
    add_geometry_file,
    load_geometry,
)
from .validator import validate_package
from .analysis import (
    AnalysisAdapter,
    AnalysisPlan,
    AnalysisResult,
    ExecutionConfig,
    available_backends,
    load_plan as load_analysis_plan,
    register_backend,
)
from .analysis.cli import analyze_package
from .signing import (
    sign_package,
    verify_package,
    list_signatures,
    SignatureMethod,
    VerificationStatus,
    SignatureInfo,
    VerificationResult,
    SigningError,
)


def view_package(package_path, *, strict: bool = False):
    """Import and invoke the interactive viewer lazily."""

    from .viewer import view_package as _view  # local import to avoid pyglet init during tests

    return _view(package_path, strict=strict)

__all__ = [
    "PACKAGE_SCHEMA",
    "MANIFEST_FILENAME",
    "PackageManifest",
    "create_package_from_entities",
    "add_geometry_file",
    "load_geometry",
    "validate_package",
    "AnalysisAdapter",
    "AnalysisPlan",
    "AnalysisResult",
    "ExecutionConfig",
    "available_backends",
    "load_analysis_plan",
    "register_backend",
    "analyze_package",
    "view_package",
    # Signing
    "sign_package",
    "verify_package",
    "list_signatures",
    "SignatureMethod",
    "VerificationStatus",
    "SignatureInfo",
    "VerificationResult",
    "SigningError",
]
