"""Reassigning an outer variable inside for/if blocks updates it."""
import textwrap

import pytest

from yapcad.dsl import compile_and_run


def run(source, command, **params):
    result = compile_and_run(textwrap.dedent(source), command, params)
    assert result.success, result.error_message
    return result.geometry


SOURCE = """
module blocks

command LOOP_SUM(n: int = 4) -> float:
    total: float = 0.0
    for i in range(n):
        total = total + 1.0
    emit total

command NESTED(n: int = 4) -> float:
    total: float = 0.0
    for i in range(n):
        if i > 1:
            total = total + 10.0
        else:
            total = total + 1.0
    emit total

command IF_SET(flag: bool = true) -> float:
    x: float = 0.0
    if flag:
        x = 5.0
    emit x

command TYPED_SHADOWS() -> float:
    x: float = 1.0
    for i in range(3):
        x: float = 100.0
    emit x

command BLOCK_LOCAL() -> float:
    total: float = 0.0
    for i in range(3):
        step = 2.0
        total = total + step
    emit total

command STACK(n: int = 3) -> solid:
    s: solid = box(1.0, 1.0, 1.0)
    for i in range(n):
        s = union(s, translate(box(1.0, 1.0, 1.0), 0.0, 0.0, 1.0 * (i + 1)))
    emit s
"""


@pytest.mark.parametrize(
    "command,expected",
    [("LOOP_SUM", 4.0), ("NESTED", 22.0), ("IF_SET", 5.0), ("TYPED_SHADOWS", 1.0), ("BLOCK_LOCAL", 6.0)],
)
def test_block_assignment(command, expected):
    assert run(SOURCE, command) == expected


def test_loop_accumulates_solids():
    from yapcad.geom3d import solidbbox

    box = solidbbox(run(SOURCE, "STACK"))
    assert box[1][2] == pytest.approx(3.5)
