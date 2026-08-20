"""Product-definition packaging for mixed fabricated/COTS assemblies."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from yapcad.assembly import Assembly, PartDefinition
from yapcad.dsl import check, parse, tokenize
from yapcad.dsl.packaging import package_from_dsl
from yapcad.geom3d_util import prism
from yapcad.metadata import get_solid_metadata
from yapcad.package import (
    ASSEMBLY_PACKAGE_SCHEMA,
    PackageManifest,
    create_package_from_assembly,
    validate_package,
)


def _translation(x: float):
    result = np.eye(4)
    result[0, 3] = x
    return result


def _definition(component_id: str, *, disposition: str, quantity: float = 1.0):
    part = PartDefinition(component_id)
    part.component_id = component_id
    part.component_name = component_id.replace("-", " ").title()
    part.disposition = disposition
    part.part_number = component_id.upper()
    part.revision = "A"
    part.quantity_per_instance = quantity
    part.unit = "each"
    return part


def _mixed_assembly():
    assembly = Assembly("test_rover")
    wheel = _definition("wheel-hub", disposition="make")
    wheel.material = "petg"
    wheel.manufacturing = {"process": "FDM"}
    bearing = _definition("bearing-608-2rs", disposition="buy", quantity=2)
    bearing.procurement = {"specification": "608-2RS, 8x22x7 mm"}

    assembly.add_part(wheel, "left_wheel", _translation(-20), prism(8, 4, 3))
    assembly.add_part(wheel, "right_wheel", _translation(20), prism(8, 4, 3))
    assembly.add_part(bearing, "left_bearings", _translation(-20), prism(3, 3, 2))
    assembly.add_part(bearing, "right_bearings", _translation(20), prism(3, 3, 2))
    return assembly


def test_create_assembly_package_preserves_components_instances_and_bom(tmp_path: Path):
    manifest = create_package_from_assembly(
        _mixed_assembly(), tmp_path / "rover.ycpkg",
        name="Test rover", version="0.1.0", root_part="left_wheel",
    )

    assert manifest.data["schema"] == ASSEMBLY_PACKAGE_SCHEMA
    assert [c["id"] for c in manifest.data["components"]] == [
        "wheel-hub", "bearing-608-2rs",
    ]
    assert len(manifest.data["instances"]) == 4
    assert {i["component"] for i in manifest.data["instances"]} == {
        "wheel-hub", "bearing-608-2rs",
    }
    assert all(len(i["transform"]) == 4 for i in manifest.data["instances"])

    bom = json.loads((manifest.root / manifest.data["bom"]["path"]).read_text())
    by_component = {item["component"]: item for item in bom["items"]}
    assert by_component["wheel-hub"]["quantity"] == 2
    assert by_component["bearing-608-2rs"]["quantity"] == 4
    assert by_component["wheel-hub"]["disposition"] == "make"
    assert by_component["bearing-608-2rs"]["disposition"] == "buy"

    for component in manifest.data["components"]:
        assert (manifest.root / component["geometry"]["path"]).exists()

    reloaded = PackageManifest.load(manifest.root)
    assert reloaded.get_component("wheel-hub")["disposition"] == "make"
    assert len(reloaded.load_component_geometry("wheel-hub")) == 1
    assert reloaded.load_bom()["items"] == bom["items"]


def test_native_assembly_record_retains_component_references_and_transforms(tmp_path: Path):
    manifest = create_package_from_assembly(
        _mixed_assembly(), tmp_path / "rover.ycpkg",
        name="Test rover", version="0.1.0", root_part="left_wheel",
    )
    snapshot = json.loads(
        (manifest.root / manifest.data["assembly"]["path"]).read_text()
    )
    assert snapshot["name"] == "test_rover"
    assert snapshot["rootPart"] == "left_wheel"
    assert snapshot["parts"]["right_wheel"]["component"] == "wheel-hub"
    assert snapshot["parts"]["right_wheel"]["transform"][0][3] == 20
    assert snapshot["mates"] == []
    assert snapshot["jointCouplings"] == []


def test_strict_validation_checks_product_definition_references_and_bom(tmp_path: Path):
    manifest = create_package_from_assembly(
        _mixed_assembly(), tmp_path / "rover.ycpkg",
        name="Test rover", version="0.1.0", root_part="left_wheel",
    )
    ok, messages = validate_package(manifest.root, strict=True)
    assert ok, messages

    manifest.data["instances"][0]["component"] = "missing-component"
    manifest.save()
    ok, messages = validate_package(manifest.root, strict=True)
    assert not ok
    assert any("missing-component" in message for message in messages)


def test_validation_rejects_edited_bom_quantity(tmp_path: Path):
    manifest = create_package_from_assembly(
        _mixed_assembly(), tmp_path / "rover.ycpkg",
        name="Test rover", version="0.1.0", root_part="left_wheel",
    )
    bom_path = manifest.root / manifest.data["bom"]["path"]
    bom = json.loads(bom_path.read_text())
    bom["items"][0]["quantity"] = 99
    bom_path.write_text(json.dumps(bom))
    manifest.recompute_hashes()
    manifest.save()

    ok, messages = validate_package(manifest.root, strict=True)
    assert not ok
    assert any("BOM quantity" in message for message in messages)


def test_repeated_component_id_rejects_conflicting_local_geometry(tmp_path: Path):
    assembly = Assembly("bad_reuse")
    definition = _definition("same-component", disposition="make")
    assembly.add_part(definition, "one", geometry=prism(2, 2, 2))
    assembly.add_part(definition, "two", geometry=prism(3, 2, 2))
    with pytest.raises(ValueError, match="conflicting local geometry"):
        create_package_from_assembly(
            assembly, tmp_path / "bad.ycpkg", name="bad", version="0.1.0"
        )


def test_buy_component_may_omit_surrogate_geometry(tmp_path: Path):
    assembly = Assembly("mixed_geometry")
    printed = _definition("printed-bracket", disposition="make")
    printed.manufacturing = {"process": "FDM"}
    bearing = _definition("bearing-608-2rs", disposition="buy")
    bearing.procurement = {"specification": "608-2RS, 8x22x7 mm"}
    assembly.add_part(printed, "bracket", geometry=prism(4, 4, 2))
    assembly.add_part(bearing, "bearing")

    manifest = create_package_from_assembly(
        assembly, tmp_path / "mixed.ycpkg", name="mixed", version="0.1.0"
    )

    component = manifest.get_component("bearing-608-2rs")
    assert "geometry" not in component
    assert manifest.load_component_geometry("bearing-608-2rs") == []
    ok, messages = validate_package(manifest.root, strict=True)
    assert ok, messages


DSL_ASSEMBLY = r'''
module packaged_assembly

@meta(
    component.id="block",
    component.name="Printed block",
    component.disposition="make",
    component.part_number="TEST-001",
    component.revision="A",
    manufacturing.process="FDM",
    material="PETG",
    assembly.datums=[{
        "id": "joint", "kind": "axis",
        "origin_mm": [0.0, 0.0, 0.0], "direction": [0.0, 0.0, 1.0]
    }]
)
command BLOCK() -> solid:
    emit box(10.0, 10.0, 10.0)

command BUILD() -> solid:
    let product: assembly = assembly("packaged")
    add_part(product, BLOCK(), "base")
    add_part(product, BLOCK(), "child")
    add_named_mate(product, "fixed", "rigid",
                   "base", "joint", "child", "joint")
    solve_assembly(product, "base")
    emit assembly_compound(product)
'''


def test_component_and_manufacturing_metadata_are_recognised_by_checker():
    result = check(parse(tokenize(DSL_ASSEMBLY), source=DSL_ASSEMBLY))
    warnings = [d for d in result.diagnostics if d.code == "W310"]
    assert warnings == []


def test_dsl_package_uses_retained_assembly_instead_of_flattened_compound(tmp_path: Path):
    result = package_from_dsl(
        DSL_ASSEMBLY, "BUILD", {}, tmp_path / "dsl.ycpkg",
        name="DSL assembly", version="0.1.0",
    )
    assert result.success, result.error_message
    assert result.manifest.data["schema"] == ASSEMBLY_PACKAGE_SCHEMA
    assert len(result.manifest.data["components"]) == 1
    assert len(result.manifest.data["instances"]) == 2
    assert result.manifest.data["components"][0]["manufacturing"] == {
        "process": "FDM"
    }
    assert result.manifest.data["components"][0]["material"] == "PETG"

    primary_path = result.manifest.geometry_primary_path()
    primary = json.loads(primary_path.read_text())
    assert len([e for e in primary["entities"] if e["type"] == "solid"]) == 1
    assembly = json.loads(
        (result.manifest.root / result.manifest.data["assembly"]["path"]).read_text()
    )
    assert set(assembly["parts"]) == {"base", "child"}


def test_metadata_bridge_preserves_component_and_procurement_sections():
    source = DSL_ASSEMBLY.replace(
        'component.revision="A",',
        'component.revision="A", procurement.specification="generic block",',
    )
    from yapcad.dsl import compile_and_run

    result = compile_and_run(source, "BLOCK", {})
    assert result.success
    meta = get_solid_metadata(result.geometry, create=False)
    assert meta["component"]["id"] == "block"
    assert meta["component"]["disposition"] == "make"
    assert meta["procurement"]["specification"] == "generic block"
    assert meta["manufacturing"]["process"] == "FDM"
