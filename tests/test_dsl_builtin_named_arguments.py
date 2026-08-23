"""Contract tests for named arguments on DSL built-in functions."""

import textwrap

import pytest

from yapcad.dsl import check, compile_and_run, parse, tokenize
from yapcad.geom3d import solidbbox


def _check(source):
    source = textwrap.dedent(source)
    return check(parse(tokenize(source), source=source))


def _run_box(expression):
    source = f"""
        module named_builtin
        command BUILD() -> solid:
            emit {expression}
    """
    return compile_and_run(textwrap.dedent(source), "BUILD", {})


def test_builtin_accepts_fully_named_arguments():
    result = _run_box("box(width=10.0, depth=20.0, height=30.0)")
    assert result.success, result.error_message
    bbox = solidbbox(result.geometry)
    assert bbox[1][0] - bbox[0][0] == pytest.approx(10.0)
    assert bbox[1][1] - bbox[0][1] == pytest.approx(20.0)
    assert bbox[1][2] - bbox[0][2] == pytest.approx(30.0)


def test_builtin_accepts_mixed_positional_and_named_arguments():
    result = _run_box("box(10.0, height=30.0, depth=20.0)")
    assert result.success, result.error_message
    bbox = solidbbox(result.geometry)
    assert bbox[1][0] - bbox[0][0] == pytest.approx(10.0)
    assert bbox[1][1] - bbox[0][1] == pytest.approx(20.0)
    assert bbox[1][2] - bbox[0][2] == pytest.approx(30.0)


@pytest.mark.parametrize(
    ("expression", "message"),
    [
        ("box(width=10.0, depth=20.0)", "missing required argument 'height'"),
        ("box(10.0, width=20.0, depth=30.0)", "multiple values for argument 'width'"),
        ("box(width=10.0, depth=20.0, height=30.0, color=\"red\")", "unknown parameter 'color'"),
        ("box(width=\"wide\", depth=20.0, height=30.0)", "argument 'width' expects 'float'"),
    ],
)
def test_builtin_named_argument_errors_are_static(expression, message):
    result = _check(f"""
        module bad_named_builtin
        command BUILD() -> solid:
            emit {expression}
    """)
    assert result.has_errors
    assert any(message in diagnostic.message.lower() for diagnostic in result.diagnostics)
