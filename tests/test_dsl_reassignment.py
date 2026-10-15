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
