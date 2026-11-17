"""The ivy vine example: a recursive, articulated design built as fields.

``examples/ivy_vine/ivy_vine.dsl`` leans on most of what the DSL offers at
once -- a command that calls itself twice to branch, ``str`` to name what it
adds, the field-only builtins, and revolute mates on declared datums -- so
building it is a good check that those still work together.  Built deferred,
it costs a fraction of a second.
"""

import math
import pathlib

import pytest

from yapcad import sdf
from yapcad.dsl import compile_and_run
from yapcad.rom import RomModel, clearance_sdf

DESIGN = (pathlib.Path(__file__).resolve().parent.parent
          / "examples" / "ivy_vine" / "ivy_vine.dsl")


def build(**params):
    result = compile_and_run(DESIGN.read_text(), "IVY_VINE", params,
                             representation="sdf", sdf_mesh=False,
                             source_path=str(DESIGN))
    assert result.success, result.error_message
    return RomModel.from_result(result)


@pytest.fixture(scope="module")
def vine():
    return build()


def test_the_recursion_grows_every_branch(vine):
    """Four stem links; each link with stem left ahead of it buds a shoot
    one link shorter, down two levels."""
    names = set(vine.parts)
    assert len(names) == 15
    assert {"ivy", "v", "vs", "vss", "vsss"} <= names           # the stem
    assert {"vb", "vbs", "vbss", "vsb", "vsbs", "vssb"} <= names  # level 1
    assert {"vbb", "vbbs", "vbsb", "vsbb"} <= names             # level 2
    assert len(vine.movable_joints()) == 14
    assert not vine.parts_without_fields()


@pytest.mark.parametrize("links,max_level,parts", [
    (1, 2, 2), (3, 0, 4), (3, 1, 7), (4, 1, 11), (4, 2, 15),
])
def test_the_vine_follows_its_parameters(links, max_level, parts):
    assert len(build(links=links, max_level=max_level).parts) == parts


def test_joint_limits_are_the_declared_ranges(vine):
    for joint in vine.movable_joints():
        lo, hi = math.degrees(joint.lower), math.degrees(joint.upper)
        if joint.name.endswith("b"):            # a side shoot
            assert sorted((abs(lo), abs(hi))) == pytest.approx([55.0, 135.0])
            assert lo * hi > 0.0                # it stays on its own side
        else:
            assert (lo, hi) == pytest.approx((-48.0, 48.0))


def test_every_hinge_is_clear_as_printed(vine):
    """A print-in-place part must not fuse: each jointed pair is apart in the
    printed pose, by about the 0.35 mm the design offsets."""
    poses = vine.pose(vine.joint_values())
    for parent, child in vine.mate_pairs():
        gap = clearance_sdf(vine.field(parent), vine.field(child),
                            poses[parent], poses[child], pitch_mm=0.25)
        assert 0.0 < gap.min_distance < 0.5, (parent, child)


def test_a_link_stands_on_the_bed():
    result = compile_and_run(DESIGN.read_text(), "VINE_LINK", {"level": 1},
                             representation="sdf", sdf_mesh=False,
                             source_path=str(DESIGN))
    assert result.success, result.error_message
    assert sdf.is_deferred(result.geometry)
    from yapcad.construction import get_construction
    lo, hi = sdf.from_construction(get_construction(result.geometry)).bounds
    assert lo[2] == pytest.approx(0.0, abs=1e-6)
    assert hi[2] == pytest.approx(6.4, abs=0.5)
