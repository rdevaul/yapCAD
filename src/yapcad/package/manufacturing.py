"""Manufacturing artifact generation for assembly-aware packages."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence

from yapcad.io import write_stl_brep
from yapcad.io.step import write_step_analytic

from .core import ASSEMBLY_PACKAGE_SCHEMA, PackageManifest, _compute_hash


SUPPORTED_COMPONENT_EXPORTS = frozenset({"stl", "step"})


def export_component_artifacts(
    assembly: Any,
    manifest: PackageManifest,
    *,
    formats: Sequence[str] = ("stl", "step"),
    dispositions: Iterable[str] = ("make",),
    strict_step: bool = True,
    strict_stl: bool = True,
    overwrite: bool = False,
) -> List[Dict[str, Any]]:
    """Export one canonical, component-local manufacturing model per component.

    The live ``Assembly`` is required because package geometry JSON contains a
    portable mesh representation, not the OCC BREP needed for analytic STEP.
    Repeated instances are deliberately exported once using their shared
    component definition.
    """
    if manifest.data.get("schema") != ASSEMBLY_PACKAGE_SCHEMA:
        raise ValueError("component exports require a ycpkg-spec-v0.2 package")

    normalized_formats = tuple(dict.fromkeys(item.lower() for item in formats))
    unsupported = sorted(set(normalized_formats) - SUPPORTED_COMPONENT_EXPORTS)
    if unsupported:
        raise ValueError(
            "unsupported component export format(s): " + ", ".join(unsupported)
        )
    selected_dispositions = set(dispositions)
    component_by_id = {
        component["id"]: component
        for component in manifest.data.get("components", [])
    }

    source_instances: Dict[str, str] = {}
    for instance_name, definition in assembly.parts.items():
        component_id = definition.component_id or definition.name
        source_instances.setdefault(component_id, instance_name)

    jobs = []
    for component_id, component in component_by_id.items():
        if component.get("disposition") not in selected_dispositions:
            continue
        instance_name = source_instances.get(component_id)
        if instance_name is None:
            raise ValueError(
                f"component {component_id!r} is not present in the source assembly"
            )
        if instance_name not in assembly.geometry:
            raise ValueError(f"component {component_id!r} has no exportable geometry")
        for file_format in normalized_formats:
            relative_path = (
                Path("exports") / "components" / component_id
                / f"{component_id}.{file_format}"
            )
            target = manifest.root / relative_path
            if target.exists() and not overwrite:
                raise FileExistsError(f"target file already exists: {target}")
            jobs.append((component, instance_name, file_format, relative_path, target))

    new_entries: List[Dict[str, Any]] = []
    for component, instance_name, file_format, relative_path, target in jobs:
        target.parent.mkdir(parents=True, exist_ok=True)
        geometry = assembly.get_part_geometry(instance_name, positioned=False)
        try:
            if file_format == "stl":
                brep_mesh = write_stl_brep(
                    geometry,
                    str(target),
                    fallback_to_mesh=not strict_stl,
                    validate_watertight=strict_stl,
                )
                analytic = None
            else:
                analytic = write_step_analytic(
                    geometry,
                    str(target),
                    name=component["id"],
                    fallback_to_faceted=not strict_step,
                )
        except Exception as exc:
            if target.exists():
                target.unlink()
            raise RuntimeError(
                f"{file_format.upper()} export failed for component "
                f"{component['id']!r}: {exc}"
            ) from exc

        entry: Dict[str, Any] = {
            "id": f"component-{component['id']}-{file_format}",
            "kind": "manufacturing-model",
            "format": file_format,
            "path": relative_path.as_posix(),
            "hash": _compute_hash(target),
            "sourceComponent": component["id"],
            "coordinateFrame": "component-local",
        }
        if file_format == "step":
            entry["analytic"] = bool(analytic)
        elif file_format == "stl":
            entry["brepTessellation"] = bool(brep_mesh)
            if strict_stl:
                entry["watertight"] = True
        new_entries.append(entry)

    new_ids = {entry["id"] for entry in new_entries}
    exports = [
        entry for entry in manifest.data.get("exports", []) or []
        if entry.get("id") not in new_ids
    ]
    exports.extend(new_entries)
    if exports:
        manifest.data["exports"] = exports

    exports_by_component: Dict[str, List[str]] = {}
    for entry in exports:
        component_id = entry.get("sourceComponent")
        if component_id:
            exports_by_component.setdefault(component_id, []).append(entry["id"])
    for component in component_by_id.values():
        component_exports = exports_by_component.get(component["id"])
        if component_exports:
            component["exports"] = component_exports
        else:
            component.pop("exports", None)
    manifest.save()
    return new_entries


__all__ = ["SUPPORTED_COMPONENT_EXPORTS", "export_component_artifacts"]
