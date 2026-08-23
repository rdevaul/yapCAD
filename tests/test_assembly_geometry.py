"""Contract tests for geometry-bearing, positioned assembly instances.

The assembly solver has historically retained only ``PartDefinition`` objects
and transforms.  These tests specify the geometry side of that contract: one
part definition may be instantiated repeatedly, each instance retains its
solid, and callers can retrieve either local or world-positioned geometry or a
compound containing the complete positioned assembly.

Most tests intentionally use yapCAD's tessellated prism representation so this
module runs in the normal, non-OCC test lane.  The final test additionally
checks preservation of analytic BREP data when pythonocc-core is available.
"""

from __future__ import annotations

import copy

import numpy as np
import pytest

from yapcad.assembly.assembly import Assembly, AssemblyError
from yapcad.assembly.datum import PartDefinition
from yapcad.brep import has_brep_data, occ_available
from yapcad.geom3d import bbox, issolid
from yapcad.geom3d_util import prism
from yapcad.io.step import write_step_analytic


def _translation(x: float = 0.0, y: float = 0.0, z: float = 0.0):
    transform = np.eye(4)
    transform[:3, 3] = [x, y, z]
    return transform


def _rotation_z_90_with_translation(x: float, y: float, z: float):
    """Return a rigid transform useful for checking more than translation."""
    return np.array([
        [0.0, -1.0, 0.0, x],
        [1.0, 0.0, 0.0, y],
        [0.0, 0.0, 1.0, z],
        [0.0, 0.0, 0.0, 1.0],
    ])


def _xyz_bounds(solid):
    bounds = bbox(solid)
    assert bounds, "positioned geometry must have a computable bounding box"
    # Tessellated primitive construction may carry sub-femtometre floating
    # noise even before assembly placement; bounds are a geometric, not
    # bit-for-bit, contract.
    return tuple(round(v, 12) for v in bounds[0][:3]), tuple(
        round(v, 12) for v in bounds[1][:3]
    )


@pytest.fixture
def unit_block():
    # Deliberately asymmetric so a 90-degree rotation swaps X/Y extents.
    return prism(2.0, 4.0, 6.0)


def test_part_instance_retains_geometry_and_exposes_local_access(unit_block):
    assembly = Assembly("one_block")
    assembly.add_part(
        PartDefinition("BLOCK"),
        name="block_1",
        geometry=unit_block,
        transform=_translation(10.0, 0.0, 0.0),
    )

    local = assembly.get_part_geometry("block_1", positioned=False)

    assert issolid(local)
    assert _xyz_bounds(local) == ((-1.0, -2.0, -3.0), (1.0, 2.0, 3.0))


def test_positioned_part_geometry_applies_full_rigid_transform(unit_block):
    assembly = Assembly("placed_block")
    assembly.add_part(
        PartDefinition("BLOCK"),
        name="block_1",
        geometry=unit_block,
        transform=_rotation_z_90_with_translation(10.0, 20.0, 30.0),
    )

    placed = assembly.get_part_geometry("block_1", positioned=True)

    # Original 2 x 4 x 6 extents become 4 x 2 x 6 after rotating around Z.
    assert _xyz_bounds(placed) == (
        (pytest.approx(8.0), pytest.approx(19.0), pytest.approx(27.0)),
        (pytest.approx(12.0), pytest.approx(21.0), pytest.approx(33.0)),
    )


def test_positioned_retrieval_does_not_mutate_or_alias_source_geometry(unit_block):
    original = copy.deepcopy(unit_block)
    assembly = Assembly("immutable_source")
    assembly.add_part(
        PartDefinition("BLOCK"),
        name="block_1",
        geometry=unit_block,
        transform=_translation(25.0, 0.0, 0.0),
    )

    first = assembly.get_part_geometry("block_1", positioned=True)
    second = assembly.get_part_geometry("block_1", positioned=True)

    assert unit_block == original
    # Independent BREP copies intentionally receive distinct entity IDs, so
    # compare the geometry payload rather than identity-bearing metadata.
    assert first[1:4] == second[1:4]
    assert _xyz_bounds(first) == _xyz_bounds(second)
    assert first is not unit_block
    assert second is not first
    assert first[1] is not second[1]


def test_repeated_instances_of_same_solid_are_positioned_independently(unit_block):
    assembly = Assembly("repeated_block")
    definition = PartDefinition("BLOCK")
    assembly.add_part(
        definition, name="left", geometry=unit_block,
        transform=_translation(-10.0, 0.0, 0.0),
    )
    assembly.add_part(
        definition, name="right", geometry=unit_block,
        transform=_translation(10.0, 0.0, 0.0),
    )

    positioned = assembly.positioned_parts()

    assert list(positioned) == ["left", "right"]
    assert _xyz_bounds(positioned["left"]) == (
        (-11.0, -2.0, -3.0), (-9.0, 2.0, 3.0),
    )
    assert _xyz_bounds(positioned["right"]) == (
        (9.0, -2.0, -3.0), (11.0, 2.0, 3.0),
    )
    assert positioned["left"] is not positioned["right"]
    assert positioned["left"][1] is not positioned["right"][1]
    assert unit_block == assembly.get_part_geometry("left", positioned=False)
    assert unit_block == assembly.get_part_geometry("right", positioned=False)


def test_positioned_parts_are_returned_in_assembly_insertion_order(unit_block):
    assembly = Assembly("ordered")
    definition = PartDefinition("BLOCK")
    for index, name in enumerate(("rocker", "bogie", "wheel")):
        assembly.add_part(
            definition,
            name=name,
            geometry=unit_block,
            transform=_translation(index * 10.0, 0.0, 0.0),
        )

    assert list(assembly.positioned_parts()) == ["rocker", "bogie", "wheel"]


def test_compound_geometry_contains_every_positioned_instance(unit_block):
    assembly = Assembly("compound")
    definition = PartDefinition("BLOCK")
    assembly.add_part(
        definition, name="left", geometry=unit_block,
        transform=_translation(-10.0, 0.0, 0.0),
    )
    assembly.add_part(
        definition, name="right", geometry=unit_block,
        transform=_translation(10.0, 5.0, 0.0),
    )

    compound = assembly.compound_geometry()

    assert issolid(compound)
    assert len(compound[1]) == 2 * len(unit_block[1])
    assert _xyz_bounds(compound) == (
        (-11.0, -2.0, -3.0), (11.0, 7.0, 3.0),
    )


def test_unknown_part_geometry_access_has_actionable_diagnostic(unit_block):
    assembly = Assembly("diagnostics")
    assembly.add_part(PartDefinition("BLOCK"), name="known", geometry=unit_block)

    with pytest.raises(AssemblyError, match=r"missing.*not found|not found.*missing"):
        assembly.get_part_geometry("missing")


def test_part_without_geometry_has_actionable_diagnostic():
    assembly = Assembly("diagnostics")
    assembly.add_part(PartDefinition("EMPTY"), name="empty")

    with pytest.raises(
        AssemblyError,
        match=r"(?i)(empty.*geometry|geometry.*empty)",
    ):
        assembly.get_part_geometry("empty")

    with pytest.raises(
        AssemblyError,
        match=r"(?i)(empty.*geometry|geometry.*empty)",
    ):
        assembly.compound_geometry(strict=True)


def test_non_strict_compound_skips_missing_geometry_and_keeps_diagnostics(
        unit_block):
    assembly = Assembly("partial")
    assembly.add_part(PartDefinition("BLOCK"), name="present", geometry=unit_block)
    assembly.add_part(PartDefinition("EMPTY"), name="empty")

    compound = assembly.compound_geometry(strict=False)

    assert issolid(compound)
    assert _xyz_bounds(compound) == _xyz_bounds(unit_block)
    assert any(
        "empty" in diagnostic.lower() and "geometry" in diagnostic.lower()
        for diagnostic in assembly.geometry_diagnostics
    )


@pytest.mark.requires_occ
@pytest.mark.skipif(not occ_available(), reason="pythonocc-core is not available")
def test_positioned_compound_preserves_brep_for_strict_analytic_step(
        tmp_path, unit_block):
    assert has_brep_data(unit_block), "OCC prism fixture must carry analytic BREP"
    assembly = Assembly("analytic_compound")
    definition = PartDefinition("BLOCK")
    assembly.add_part(
        definition, name="left", geometry=unit_block,
        transform=_translation(-10.0, 0.0, 0.0),
    )
    assembly.add_part(
        definition, name="right", geometry=unit_block,
        transform=_translation(10.0, 0.0, 0.0),
    )

    compound = assembly.compound_geometry()

    assert has_brep_data(compound)
    output = tmp_path / "positioned_assembly.step"
    assert write_step_analytic(
        compound, str(output), fallback_to_faceted=False,
    ) is True
    assert output.exists()
    assert output.stat().st_size > 0
