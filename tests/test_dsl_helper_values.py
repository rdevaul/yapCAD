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


def test_boolean_helper_is_a_runtime_value():
    helper = "command yes() -> bool:\n    emit true\n"
    assert (
        run(
            "    emit true and not not yes()\n",
            result_type="bool",
            helpers=helper,
        )
        is True
    )
