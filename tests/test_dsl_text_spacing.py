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


def test_text_solid_accepts_three_arguments():
    body = run('    emit text_solid("ABC", 5.0, 1.0)\n', "solid")
    assert body[1]
