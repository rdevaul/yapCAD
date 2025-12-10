"""DSL-to-Package integration.

Provides functions to compile DSL source, execute commands, and package
the resulting geometry with full provenance tracking.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

from .runtime import compile_and_run, ExecutionResult


def package_from_dsl(
    source: str,
    command_name: str,
    parameters: Dict[str, Any],
    target_dir: Path | str,
    *,
    name: str,
    version: str,
    description: Optional[str] = None,
    author: Optional[str] = None,
    units: Optional[str] = None,
    materials: Optional[Dict[str, Dict[str, Any]]] = None,
    overwrite: bool = False,
) -> "PackageResult":
    """Compile DSL source, execute a command, and package the result.

    This is a high-level function that combines DSL compilation, execution,
    and packaging into a single workflow. The resulting package includes
    full provenance metadata linking it to the original DSL source.

    Args:
        source: DSL source code containing the module and command definitions.
        command_name: Name of the command to execute (e.g., "MAKE_GEAR").
        parameters: Dictionary of parameter values for the command.
        target_dir: Directory where the package will be created.
        name: Name for the package (used in manifest).
        version: Version string for the package.
        description: Optional description for the package.
        author: Optional author name.
        units: Unit system (default "mm").
        materials: Optional materials dictionary for the package.
        overwrite: If True, overwrite existing package directory.

    Returns:
        PackageResult with success status, manifest, and any error info.

    Example:
        >>> source = '''
        ... module gear_design;
        ...
        ... command MAKE_GEAR(teeth: int, module_mm: float) -> solid {
        ...     let pitch_diameter: float = teeth * module_mm;
        ...     let gear: solid = cylinder(pitch_diameter / 2.0, 10.0);
        ...     emit gear;
        ... }
        ... '''
        >>> result = package_from_dsl(
        ...     source,
        ...     "MAKE_GEAR",
        ...     {"teeth": 24, "module_mm": 2.0},
        ...     "output/gear_pkg",
        ...     name="gear_24t",
        ...     version="1.0.0",
        ... )
        >>> if result.success:
        ...     print(f"Package created at {result.manifest.root}")
    """
    from yapcad.package import create_package_from_entities, PackageManifest
    from yapcad.geom3d import issolid, issurface

    # Step 1: Compile and execute the DSL
    exec_result = compile_and_run(source, command_name, parameters)

    if not exec_result.success:
        return PackageResult(
            success=False,
            error_message=f"DSL execution failed: {exec_result.error_message}",
            execution_result=exec_result,
        )

    # Step 2: Validate that we got geometry
    geometry = exec_result.geometry
    if geometry is None:
        return PackageResult(
            success=False,
            error_message="DSL execution produced no geometry output",
            execution_result=exec_result,
        )

    # Wrap single entity in list if needed
    if issolid(geometry) or issurface(geometry):
        entities = [geometry]
    elif isinstance(geometry, list):
        entities = geometry
    else:
        return PackageResult(
            success=False,
            error_message=f"Unexpected geometry type: {type(geometry).__name__}",
            execution_result=exec_result,
        )

    # Step 3: Build generator metadata from provenance
    generator: Dict[str, Any] = {
        "tool": "yapCAD-DSL",
        "version": _get_yapcad_version(),
    }

    if exec_result.provenance:
        prov = exec_result.provenance
        generator["dsl"] = {
            "module": prov.module_name,
            "command": prov.command_name,
            "parameters": prov.parameters,
            "source_signature": prov.source_signature,
        }
        if prov.version:
            generator["dsl"]["version"] = prov.version

    # Step 4: Create the package
    try:
        target_path = Path(target_dir)
        manifest = create_package_from_entities(
            entities,
            target_path,
            name=name,
            version=version,
            description=description,
            author=author,
            units=units,
            materials=materials,
            generator=generator,
            overwrite=overwrite,
        )

        # Step 5: Add DSL source as an attachment
        _add_dsl_source_attachment(manifest, source, command_name)

        return PackageResult(
            success=True,
            manifest=manifest,
            execution_result=exec_result,
        )

    except Exception as e:
        return PackageResult(
            success=False,
            error_message=f"Package creation failed: {e}",
            execution_result=exec_result,
        )


def _get_yapcad_version() -> str:
    """Get yapCAD version string."""
    try:
        from yapcad import __version__
        return __version__
    except ImportError:
        return "unknown"


def _add_dsl_source_attachment(
    manifest: "PackageManifest",
    source: str,
    command_name: str,
) -> None:
    """Add DSL source code as a package attachment."""
    import hashlib

    # Write source to attachments directory
    attachments_dir = manifest.root / "attachments"
    attachments_dir.mkdir(parents=True, exist_ok=True)

    source_file = attachments_dir / "source.dsl"
    source_file.write_text(source, encoding="utf-8")

    # Compute hash
    source_hash = hashlib.sha256(source.encode("utf-8")).hexdigest()

    # Add to manifest
    attachments = manifest.data.setdefault("attachments", [])
    attachments.append({
        "id": "dsl-source",
        "path": "attachments/source.dsl",
        "hash": f"sha256:{source_hash}",
        "format": "dsl",
        "purpose": f"DSL source code for {command_name} command",
    })
    manifest.save()


class PackageResult:
    """Result of DSL-to-package operation."""

    def __init__(
        self,
        *,
        success: bool,
        manifest: Optional["PackageManifest"] = None,
        execution_result: Optional[ExecutionResult] = None,
        error_message: Optional[str] = None,
    ):
        self.success = success
        self.manifest = manifest
        self.execution_result = execution_result
        self.error_message = error_message

    def __repr__(self) -> str:
        if self.success:
            return f"PackageResult(success=True, manifest={self.manifest})"
        return f"PackageResult(success=False, error_message={self.error_message!r})"


__all__ = [
    "package_from_dsl",
    "PackageResult",
]
