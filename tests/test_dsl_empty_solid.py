"""Focused DSL runtime regression tests."""

import pytest
from yapcad.dsl import compile_and_run


def run(body, result_type="float", helpers=""):
    result = compile_and_run(
        f"module test\n{helpers}\ncommand TEST() -> {result_type}:\n" + body,
        "TEST",
        {},
    )
    assert result.success, result.error_message
    return result.geometry


def test_empty_solid_can_seed_loop_accumulator():
    pytest.importorskip("manifold3d")
    from yapcad.boolean.manifold_engine import _to_manifold

    body = run(
        "    s: solid = empty_solid()\n    for i in range(3):\n"
        "        s = union(s, translate(box(1.0, 1.0, 1.0), "
        "2.0 * i, 0.0, 0.0))\n    emit s\n",
        "solid",
    )
    assert _to_manifold(body, "test").volume() == pytest.approx(3.0)
