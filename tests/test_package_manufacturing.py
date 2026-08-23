"""Manufacturing-export tests for assembly-aware yapCAD packages."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from yapcad.assembly import Assembly, PartDefinition
from yapcad.brep import _clear_brep_data, occ_available
from yapcad.dsl.packaging import package_from_dsl
from yapcad.geom3d_util import prism
from yapcad.package import (
    create_package_from_assembly,
    export_component_artifacts,
    validate_package,
)


def _definition(component_id: str, disposition: str) -> PartDefinition:
    definition = PartDefinition(component_id)
    definition.component_id = component_id
    definition.component_name = component_id
    definition.disposition = disposition
    definition.part_number = component_id.upper()
    definition.revision = "A"
    if disposition == "make":
        definition.manufacturing = {"process": "FDM"}
    else:
        definition.procurement = {"specification": "catalogue item"}
    return definition


def _assembly() -> Assembly:
    assembly = Assembly("manufacturing_test")
    bracket = _definition("printed-bracket", "make")
    bearing = _definition("bearing-608", "buy")
    shifted = np.eye(4)
    shifted[0, 3] = 100.0
    assembly.add_part(bracket, "left_bracket", geometry=prism(10, 8, 4))
    assembly.add_part(bracket, "right_bracket", shifted, prism(10, 8, 4))
    assembly.add_part(bearing, "bearing", geometry=prism(4, 4, 2))
    return assembly


def test_component_stl_exports_are_local_deduplicated_and_manifested(tmp_path: Path):
    assembly = _assembly()
    manifest = create_package_from_assembly(
        assembly, tmp_path / "product.ycpkg", name="product", version="0.1.0"
    )

    entries = export_component_artifacts(
        assembly, manifest, formats=("stl",), strict_step=True, strict_stl=False
    )

    assert len(entries) == 1
    entry = entries[0]
    assert entry["id"] == "component-printed-bracket-stl"
    assert entry["sourceComponent"] == "printed-bracket"
    assert entry["coordinateFrame"] == "component-local"
    assert entry["format"] == "stl"
    assert (manifest.root / entry["path"]).stat().st_size > 84
    assert manifest.get_component("printed-bracket")["exports"] == [entry["id"]]
    assert "exports" not in manifest.get_component("bearing-608")

    ok, messages = validate_package(manifest.root, strict=True)
    assert ok, messages


def test_export_can_include_non_make_dispositions_explicitly(tmp_path: Path):
    assembly = _assembly()
    manifest = create_package_from_assembly(
        assembly, tmp_path / "product.ycpkg", name="product", version="0.1.0"
    )
    entries = export_component_artifacts(
        assembly, manifest, formats=("stl",), dispositions=("make", "buy"),
        strict_stl=False,
    )
    assert {entry["sourceComponent"] for entry in entries} == {
        "printed-bracket", "bearing-608",
    }


def test_export_rejects_unknown_formats_before_writing(tmp_path: Path):
    assembly = _assembly()
    manifest = create_package_from_assembly(
        assembly, tmp_path / "product.ycpkg", name="product", version="0.1.0"
    )
    with pytest.raises(ValueError, match="unsupported component export format"):
        export_component_artifacts(assembly, manifest, formats=("obj",))
    assert manifest.data.get("exports") in (None, [])


def test_export_refuses_to_replace_existing_artifact_without_overwrite(tmp_path: Path):
    assembly = _assembly()
    manifest = create_package_from_assembly(
        assembly, tmp_path / "product.ycpkg", name="product", version="0.1.0"
    )
    export_component_artifacts(
        assembly, manifest, formats=("stl",), strict_stl=False
    )
    with pytest.raises(FileExistsError):
        export_component_artifacts(
            assembly, manifest, formats=("stl",), strict_stl=False
        )
    entries = export_component_artifacts(
        assembly, manifest, formats=("stl",), strict_stl=False, overwrite=True
    )
    assert len(entries) == 1


def test_strict_step_does_not_fall_back_to_faceted_output(tmp_path: Path):
    assembly = _assembly()
    manifest = create_package_from_assembly(
        assembly, tmp_path / "product.ycpkg", name="product", version="0.1.0"
    )
    _clear_brep_data(assembly.geometry["left_bracket"])
    with pytest.raises((RuntimeError, ValueError)):
        export_component_artifacts(
            assembly, manifest, formats=("step",), strict_step=True
        )
    assert not list((manifest.root / "exports").rglob("*.step"))


@pytest.mark.skipif(not occ_available(), reason="pythonocc-core is unavailable")
def test_dsl_packaging_can_emit_strict_component_stl_and_step(tmp_path: Path):
    source = r'''
module printable_package

@meta(component.id="coupon", component.name="Coupon",
      component.disposition="make", component.part_number="TEST-001",
      component.revision="A", manufacturing.process="FDM",
      assembly.datums=[{
          "id": "origin", "kind": "point", "origin_mm": [0.0, 0.0, 0.0]
      }])
command COUPON() -> solid:
    emit box(10.0, 10.0, 4.0)

command BUILD() -> solid:
    let product: assembly = assembly("printable")
    add_part(product, COUPON(), "coupon")
    solve_assembly(product, "coupon")
    emit assembly_compound(product)
'''
    result = package_from_dsl(
        source, "BUILD", {}, tmp_path / "printable.ycpkg",
        name="printable", version="0.1.0",
        component_exports=("stl", "step"), strict_component_step=True,
    )
    assert result.success, result.error_message
    exports = result.manifest.data["exports"]
    assert {entry["format"] for entry in exports} == {"stl", "step"}
    assert next(entry for entry in exports if entry["format"] == "step")[
        "analytic"
    ] is True
    stl_entry = next(entry for entry in exports if entry["format"] == "stl")
    assert stl_entry["brepTessellation"] is True
    import trimesh
    # STL stores triangles independently; merge coincident facet vertices
    # before applying the topological watertightness check.
    mesh = trimesh.load_mesh(result.manifest.root / stl_entry["path"], process=True)
    assert mesh.is_watertight
