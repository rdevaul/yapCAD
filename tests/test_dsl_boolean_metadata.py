"""Boolean, fillet and chamfer results keep the first operand's metadata."""
import textwrap

import pytest

from yapcad.dsl import compile_and_run

SOURCE = textwrap.dedent('''
module meta

@meta(material={"name": "NdFeB N42SH", "color": [0.85, 0.22, 0.2]})
command MAGNET() -> solid:
    emit box(10.0, 5.0, 3.0)

command CUT() -> solid:
    emit difference(MAGNET(), cylinder(1.0, 10.0))

command TRIM() -> solid:
    emit intersection(MAGNET(), box(8.0, 8.0, 8.0))

command JOIN() -> solid:
    emit union(MAGNET(), translate(box(2.0, 2.0, 2.0), 6.0, 0.0, 0.0))

command JOIN_ALL() -> solid:
    emit union_all([MAGNET(), translate(box(2.0, 2.0, 2.0), 6.0, 0.0, 0.0)])

command CUT_ALL() -> solid:
    emit difference_all(MAGNET(), [cylinder(1.0, 10.0)])

command PLAIN() -> solid:
    emit difference(box(4.0, 4.0, 4.0), cylinder(1.0, 10.0))
''')


def metadata(command):
    result = compile_and_run(SOURCE, command, {})
    assert result.success, result.error_message
    solid = result.geometry
    return solid[4] if len(solid) > 4 else {}


@pytest.mark.parametrize("command", ["CUT", "TRIM", "JOIN", "JOIN_ALL", "CUT_ALL"])
def test_first_operand_material_survives(command):
    material = metadata(command)["material"]
    assert material["name"] == "NdFeB N42SH"
    assert material["color"] == [0.85, 0.22, 0.2]


def test_identity_is_not_copied():
    source = metadata("MAGNET")
    derived = metadata("CUT")
    assert "id" not in derived or derived.get("id") != source.get("id")


def test_untagged_operands_stay_untagged():
    assert "material" not in metadata("PLAIN")
