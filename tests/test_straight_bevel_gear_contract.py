"""Public contract for straight bevel gear generation modes."""

import textwrap

import pytest

from yapcad.dsl import check, compile_and_run, parse, tokenize


def _source(generation_expression):
    return textwrap.dedent(f"""
        module bevel_contract
        command BUILD() -> solid:
            emit straight_bevel_gear(
                teeth=24,
                mate_teeth=24,
                outer_module_mm=1.5,
                face_width_mm=8.0,
                shaft_angle_deg=90.0,
                pressure_angle_deg=20.0,
                backlash_mm=0.25,
                bore_diameter_mm=8.3,
                generation_type={generation_expression},
            )
    """)


def test_unknown_literal_generation_type_is_a_static_error():
    source = _source('"cycloidal"')
    result = check(parse(tokenize(source), source=source))
    assert result.has_errors
    assert any(
        "generation_type" in diagnostic.message
        and "cycloidal" in diagnostic.message
        for diagnostic in result.diagnostics
    )


def test_unknown_positional_generation_type_is_a_static_error():
    source = textwrap.dedent("""
        module positional_bevel_contract
        command BUILD() -> solid:
            emit straight_bevel_gear(
                24, 24, 1.5, 8.0, 90.0, 20.0, 0.25, 8.3, "cycloidal"
            )
    """)
    result = check(parse(tokenize(source), source=source))
    assert result.has_errors
    assert any("cycloidal" in diagnostic.message for diagnostic in result.diagnostics)


def test_dynamic_unknown_generation_type_is_a_runtime_error():
    source = textwrap.dedent("""
        module dynamic_bevel_contract
        command BUILD(generation: string) -> solid:
            emit straight_bevel_gear(
                24, 24, 1.5, 8.0, 90.0, 20.0, 0.25, 8.3, generation
            )
    """)
    checked = check(parse(tokenize(source), source=source))
    assert not checked.has_errors
    result = compile_and_run(source, "BUILD", {"generation": "cycloidal"})
    assert not result.success
    assert "cycloidal" in result.error_message


@pytest.mark.parametrize("generation_type", ["octoid", "gleason"])
def test_reserved_generation_types_are_recognized_then_rejected_at_runtime(
    generation_type,
):
    source = _source(f'"{generation_type}"')
    checked = check(parse(tokenize(source), source=source))
    assert not checked.has_errors, [d.message for d in checked.diagnostics]

    executed = compile_and_run(source, "BUILD", {})
    assert not executed.success
    assert generation_type in executed.error_message
    assert "reserved" in executed.error_message


def test_miter_gear_convenience_has_optional_generation_type():
    source = textwrap.dedent("""
        module miter_contract
        command BUILD() -> solid:
            emit miter_gear(
                teeth=24,
                outer_module_mm=1.5,
                face_width_mm=8.0,
                bore_diameter_mm=8.3,
                generation_type="octoid",
            )
    """)
    checked = check(parse(tokenize(source), source=source))
    assert not checked.has_errors, [d.message for d in checked.diagnostics]
    result = compile_and_run(source, "BUILD", {})
    assert not result.success
    assert "octoid" in result.error_message
