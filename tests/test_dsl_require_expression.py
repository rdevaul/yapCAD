"""A failing require records which condition failed and where."""
import textwrap

from yapcad.dsl import compile_and_run

SOURCE = textwrap.dedent("""
    module plate
    command PLATE(width: float = 60.0, bore: float = 10.0) -> solid:
        require width > 0.0
        require width > bore + 4.0 * 4.5
        emit box(width, 30.0, 6.0)
""")


def test_failing_require_names_its_condition():
    result = compile_and_run(SOURCE, "PLATE", {"width": 10.0})
    assert len(result.require_failures) == 1
    failure = result.require_failures[0]
    assert failure.expression_text == "line 5: width > bore + 4.0 * 4.5"
    assert "width > bore + 4.0 * 4.5" in str(failure)


def test_passing_requires_record_nothing():
    assert compile_and_run(SOURCE, "PLATE", {}).require_failures == []
