"""Building typed lists in loops, keeping pre-placed parts where they are in
an assembly, and keeping each child's identity in a compound (#80)."""

import textwrap

import numpy as np
import pytest

from yapcad.assembly.assembly import Assembly, PartDefinition
from yapcad.dsl import compile_and_run
from yapcad.geom3d import solidbbox
from yapcad.metadata import get_solid_metadata, get_surface_metadata


def run(source, command="C"):
    return compile_and_run(textwrap.dedent(source), command, {})


# ---------------------------------------------------------------------------
# typed lists
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("element, grow, expected", [
    ("solid", "append(xs, box(1.0, 1.0, 1.0))", "list<solid>"),
    ("solid", "concat(xs, [box(1.0, 1.0, 1.0)])", "list<solid>"),
    ("float", "append(xs, 1.5)", "list<float>"),
    # An int joins a list of float without making it a list of int.
    ("float", "append(xs, 2)", "list<float>"),
    ("string", "append(xs, \"a\")", "list<string>"),
])
def test_a_seeded_list_is_tagged_with_what_it_holds(element, grow, expected):
    """An empty literal has no element to take a type from, and used to be
    a list of int at run time whatever it was declared as -- so a
    ``list<solid>`` built in a loop reached its consumers tagged
    ``list<int>``, and a loop over it handed out solids labelled int."""
    result = run(f"""
        module t
        command C() -> list<{element}>:
            xs: list<{element}> = []
            for i in range(2):
                xs = {grow}
            emit xs
    """)
    assert result.success, result.error_message
    assert str(result.emit_result.value.type) == expected
    assert len(result.emit_result.value.data) == 2


def test_a_list_emptied_by_reassignment_keeps_its_type():
    result = run("""
        module t
        command C() -> list<float>:
            xs: list<float> = [1.0, 2.0]
            xs = []
            xs = append(xs, 3.5)
            emit xs
    """)
    assert result.success, result.error_message
    assert str(result.emit_result.value.type) == "list<float>"
    assert result.emit_result.value.data == [3.5]


def test_an_untyped_empty_list_is_still_a_list_of_int():
    from yapcad.dsl.runtime.builtins import call_builtin
    from yapcad.dsl.runtime.values import Value, float_val, int_val
    from yapcad.dsl.types import INT, ListType
    empty = Value([], ListType(INT))
    assert str(call_builtin("append", [empty, int_val(1)]).type) == "list<int>"
    # With nothing to go on but the item, the item decides.
    assert str(call_builtin("append", [empty, float_val(1.5)]).type) \
        == "list<float>"

def test_an_empty_typed_list_is_filled_with_append():
    result = run("""
        module t
        command C() -> list<solid>:
            parts: list<solid> = []
            for i in range(3):
                parts = append(parts, translate(box(1.0, 1.0, 1.0), 2.0 * i, 0.0, 0.0))
            emit parts
    """)
    assert result.success, result.error_message
    assert len(result.emit_result.value.data) == 3


def test_concat_onto_an_empty_seed_keeps_the_element_type():
    result = run("""
        module t
        command C() -> list<solid>:
            parts: list<solid> = []
            for i in range(2):
                parts = concat(parts, [box(1.0, 1.0, 1.0)])
            emit parts
    """)
    assert result.success, result.error_message
    assert len(result.emit_result.value.data) == 2


def test_append_of_the_wrong_type_is_refused():
    result = run("""
        module t
        command C() -> list<float>:
            xs: list<float> = []
            xs = append(xs, box(1.0, 1.0, 1.0))
            emit xs
    """)
    assert not result.success
    assert "Cannot append 'solid' to 'list<float>'" in result.error_message


def test_append_does_not_change_its_argument():
    result = run("""
        module t
        command C() -> int:
            xs: list<int> = [1, 2]
            ys: list<int> = append(xs, 3)
            emit len(xs) * 10 + len(ys)
    """)
    assert result.success, result.error_message
    assert result.emit_result.value.data == 23


def test_a_non_empty_literal_still_has_to_match():
    result = run("""
        module t
        command C() -> list<solid>:
            parts: list<solid> = [1.0]
            emit parts
    """)
    assert not result.success


# ---------------------------------------------------------------------------
# fixed placement
# ---------------------------------------------------------------------------

ASSEMBLY = """
    module t
    command C() -> assembly:
        a: assembly = assembly("t")
        a = add_part(a, box(1.0, 1.0, 1.0), "base")
        a = add_part(a, translate(box(1.0, 1.0, 1.0), 5.0, 0.0, 0.0), "loose")
        {fix}
        a = solve_assembly(a, "base")
        emit a
"""


def test_an_unmated_part_is_still_refused():
    from yapcad.assembly.assembly import AssemblyError
    with pytest.raises(AssemblyError, match="fix_part"):
        run(ASSEMBLY.format(fix=""))


def test_a_fixed_part_keeps_the_pose_it_was_built_in():
    result = run(ASSEMBLY.format(fix='a = fix_part(a, "loose")'))
    assert result.success, result.error_message
    asm = result.emit_result.value.data
    assert np.allclose(asm.transforms["loose"], np.eye(4))
    lo, hi = solidbbox(asm.get_part_geometry("loose"))
    assert lo[0] == pytest.approx(4.5) and hi[0] == pytest.approx(5.5)


def test_parts_mated_to_a_fixed_part_are_placed_from_it():
    asm = Assembly("t")
    asm.add_part(PartDefinition("base"), name="base")
    shifted = np.eye(4)
    shifted[0, 3] = 10.0
    asm.add_part(PartDefinition("island"), name="island", transform=shifted)
    asm.fix_part("island")
    result = asm.solve("base")
    assert result.success, result.errors
    assert np.allclose(result.transforms["island"], shifted)


def test_fixing_an_unknown_part_fails():
    from yapcad.assembly.assembly import AssemblyError
    with pytest.raises(AssemblyError, match="not found"):
        Assembly("t").fix_part("ghost")


# ---------------------------------------------------------------------------
# compound identity
# ---------------------------------------------------------------------------

def test_a_compound_keeps_each_childs_material():
    result = run("""
        module t
        @meta(material={"name": "red", "color": [1.0, 0.0, 0.0]})
        command RED() -> solid:
            emit box(1.0, 1.0, 1.0)

        @meta(material={"name": "blue", "color": [0.0, 0.0, 1.0]},
              magnet={"M": [0.0, 1.0, 0.0], "Br_T": 1.28})
        command BLUE() -> solid:
            emit translate(box(1.0, 1.0, 1.0), 3.0, 0.0, 0.0)

        command C() -> solid:
            emit compound(RED(), BLUE())
    """)
    assert result.success, result.error_message
    solid_ = result.geometry
    children = get_solid_metadata(solid_)["compound"]["children"]
    assert [c["material"]["name"] for c in children] == ["red", "blue"]
    assert children[1]["magnet"] == {"M": [0.0, 1.0, 0.0], "Br_T": 1.28}
    assert "brep" not in children[0]
    for child in children:
        start, stop = child["surfaces"]
        names = {get_surface_metadata(s).get("material", {}).get("name")
                 for s in solid_[1][start:stop]}
        assert names == {child["material"]["name"]}
