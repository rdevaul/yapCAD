"""A command called twice with the same arguments is built once.

A command's result depends only on its arguments, so the interpreter answers
a repeated call from the first one.  That is what lets a recursive design --
fourteen links of a vine, eight of them distinct -- mesh and simplify each
distinct part once: the repeats are the *same* solid, not equal copies.
"""

import pytest

from yapcad import sdf
from yapcad.dsl.lexer import tokenize
from yapcad.dsl.parser import parse
from yapcad.dsl.runtime.interpreter import Interpreter
from yapcad.dsl import compile_and_run
from yapcad.metadata import get_solid_metadata


def interpret(source, command="TOP", monkeypatch=None):
    interpreter = Interpreter()
    result = interpreter.execute(parse(tokenize(source)), command, {}, source)
    assert result.success, result.error_message
    return interpreter, result


PARTS = '''
module memo

@meta(material="PLA")
@meta(assembly.datums=[
    {"id": "o", "kind": "axis",
     "origin_mm": [0.0, 0.0, 0.0], "direction": [0.0, 0.0, 1.0]}
])
command LINK(level: int, side: float) -> solid:
    emit box(10.0 - level, 4.0, 2.0 + side)

command grow(a: assembly, parent: string, name: string, n: int) -> assembly:
    if n <= 0:
        emit a
    add_part(a, LINK(n % 2, 1.0), name)
    add_named_mate(a, str("j_", name), "rigid", parent, "o", name, "o")
    emit grow(a, name, str(name, "s"), n - 1)

command TOP() -> solid:
    let a: assembly = assembly("t")
    add_part(a, LINK(0, 1.0), "root")
    grow(a, "root", "v", 5)
    emit compound(LINK(0, 1.0), translate(LINK(1, 1.0), 20.0, 0.0, 0.0))
'''


def test_repeated_calls_are_answered_from_the_first():
    interpreter, _ = interpret(PARTS)
    # Eight calls to LINK, two distinct argument lists.
    assert interpreter.memo_hits == 6


def test_repeats_are_the_same_solid():
    source = PARTS.replace(
        'emit compound(LINK(0, 1.0), translate(LINK(1, 1.0), 20.0, 0.0, 0.0))',
        'emit assembly_compound(solve_assembly(a, "root"))')
    _, result = interpret(source)
    assembly = result.emit_result.value.annotations["assembly"]
    solids = assembly.geometry
    assert solids["root"] is solids["vs"] is solids["vsss"]      # LINK(0, 1.0)
    assert solids["v"] is solids["vss"] is solids["vssss"]       # LINK(1, 1.0)
    assert solids["root"] is not solids["v"]


def test_the_environment_turns_it_off(monkeypatch):
    monkeypatch.setenv("YAPCAD_DSL_MEMOIZE", "0")
    interpreter, _ = interpret(PARTS)
    assert interpreter.memo_hits == 0


def test_int_float_and_bool_arguments_are_different_calls():
    source = '''
module memo
command kind(x: float) -> float:
    emit x + 1.0
command flag(b: bool) -> float:
    emit 10.0 if b else 20.0
command TOP() -> float:
    emit kind(1) + kind(1.0) + kind(1.0) + flag(true) + flag(false)
'''
    interpreter, result = interpret(source)
    assert result.emit_result.data == pytest.approx(36.0)
    assert interpreter.memo_hits == 1          # only the second kind(1.0)


def test_a_command_given_an_assembly_always_runs():
    source = '''
module memo
command add(a: assembly, name: string) -> assembly:
    add_part(a, box(1.0, 1.0, 1.0), name)
    emit a
command TOP() -> string:
    let a: assembly = assembly("t")
    add(a, "p")
    let b: assembly = assembly("u")
    add(b, "p")
    emit str(assembly_report(a), assembly_report(b))
'''
    interpreter, result = interpret(source)
    assert interpreter.memo_hits == 0
    assert result.emit_result.data.count("Parts: 1") == 2


def test_a_command_that_makes_an_assembly_makes_a_new_one_each_time():
    source = '''
module memo
command fresh(name: string) -> assembly:
    emit assembly(name)
command TOP() -> string:
    let a: assembly = fresh("t")
    let b: assembly = fresh("t")
    add_part(a, box(1.0, 1.0, 1.0), "only_in_a")
    emit assembly_report(b)
'''
    _, result = interpret(source)
    assert "only_in_a" not in result.emit_result.data


def test_a_command_that_prints_prints_every_time(capsys):
    source = '''
module memo
command noisy(x: float) -> float:
    print("building", x)
    emit x
command TOP() -> float:
    emit noisy(2.0) + noisy(2.0)
'''
    interpreter, _ = interpret(source)
    assert capsys.readouterr().out.count("building") == 2
    assert interpreter.memo_hits == 0


def test_a_native_function_is_not_assumed_pure():
    source = '''
module memo
native python {
calls = [0]

def tick(x):
    calls[0] += 1
    return float(calls[0])
} exports {
    fn tick(x: float) -> float;
}

command wrapped(x: float) -> float:
    emit tick(x)
command TOP() -> float:
    emit wrapped(1.0) + wrapped(1.0)
'''
    result = compile_and_run(source, "TOP", {})
    assert result.success, result.error_message
    assert result.emit_result.data == pytest.approx(3.0)       # 1 + 2


def test_a_hint_on_the_caller_does_not_reach_the_cached_part():
    """WRAPPED emits LINK's cached solid unchanged and adds a hint of its
    own; a later LINK() must not come back carrying it."""
    source = '''
module memo
@meta(material="PLA")
command LINK(n: int) -> solid:
    emit box(4.0, 4.0, 4.0)

@meta(material="steel", tags=["wrapped"])
command WRAPPED(n: int) -> solid:
    emit LINK(n)

command TOP() -> list<solid>:
    let first: solid = LINK(1)
    let wrapped: solid = WRAPPED(1)
    emit [first, wrapped, LINK(1)]
'''
    _, result = interpret(source)
    first, wrapped, again = [
        v.data if hasattr(v, "data") else v for v in result.emit_result.data]
    assert again is first
    assert wrapped is not first
    assert get_solid_metadata(wrapped)["material"] == "steel"
    assert get_solid_metadata(first)["material"] == "PLA"
    assert "wrapped" not in (get_solid_metadata(first).get("tags") or [])


def test_a_shared_part_is_meshed_once_when_deferred():
    source = PARTS.replace(
        'emit compound(LINK(0, 1.0), translate(LINK(1, 1.0), 20.0, 0.0, 0.0))',
        'emit assembly_compound(solve_assembly(a, "root"))')
    result = compile_and_run(source, "TOP", {}, representation="sdf",
                             sdf_mesh=False)
    assert result.success, result.error_message
    from yapcad.rom import assembly_from_result
    solids = assembly_from_result(result).geometry
    assert sdf.is_deferred(solids["root"]) and sdf.is_deferred(solids["vs"])
    sdf.realize(solids["root"])
    # One mesh, built once, seen by every instance of that link.
    assert not sdf.is_deferred(solids["vs"])
    assert solids["vs"][1] is solids["root"][1]
    assert sdf.is_deferred(solids["v"])
