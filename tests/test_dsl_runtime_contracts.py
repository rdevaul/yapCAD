"""Checker-accepted programs obey runtime Value and control-flow contracts."""

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


def test_repeated_assignment_updates_existing_binding():
    assert (
        run(
            "    x: float = 1.0\n    x = x + 1.0\n"
            "    for i in range(3):\n        x = x + 2.0\n"
            "        x = x + 3.0\n    emit x\n"
        )
        == 17.0
    )


def test_assignment_still_checks_type_and_typed_declaration_shadows():
    result = compile_and_run(
        "module t\ncommand TEST() -> float:\n"
        '    x: float = 1.0\n    x = "bad"\n    emit x\n',
        "TEST",
        {},
    )
    assert not result.success
    assert (
        run(
            "    x: float = 1.0\n    for i in range(2):\n"
            "        x: float = 99.0\n    emit x\n"
        )
        == 1.0
    )


def test_emit_exits_nested_blocks_and_loop():
    assert (
        run(
            "    for i in range(3):\n        if i == 1:\n"
            "            emit 12.0\n    require false\n    emit 999.0\n"
        )
        == 12.0
    )


def test_boolean_helper_is_a_runtime_value():
    helper = "command yes() -> bool:\n    emit true\n"
    assert (
        run(
            "    if true and not not yes():\n        emit 2.0\n    emit 1.0\n",
            helpers=helper,
        )
        == 2.0
    )


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


def test_text_solid_accepts_three_arguments():
    body = run('    emit text_solid("ABC", 5.0, 1.0)\n', "solid")
    assert body[1]
