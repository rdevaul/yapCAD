"""Focused DSL runtime regression tests."""

from yapcad.dsl import compile_and_run


def run(body, result_type="float", helpers=""):
    result = compile_and_run(
        f"module test\n{helpers}\ncommand TEST() -> {result_type}:\n" + body,
        "TEST",
        {},
    )
    assert result.success, result.error_message
    return result.geometry


def test_emit_exits_nested_blocks_and_loop():
    assert (
        run(
            "    for i in range(3):\n        if i == 1:\n"
            "            emit 12.0\n    require false\n    emit 999.0\n"
        )
        == 12.0
    )
