"""Cross-file ``use`` imports in the DSL.

These cover the behaviours that are easy to get subtly wrong and that fail
silently when they are: an aliased library whose commands call each other, a
command name arriving from two different modules, the same module arriving
twice by different routes, and a module edited between two builds in one
process.
"""

import pytest

from yapcad.dsl import compile_and_run
from yapcad.dsl.checker import check
from yapcad.dsl.lexer import tokenize
from yapcad.dsl.modules import clear_module_cache, resolve_uses
from yapcad.dsl.parser import parse
from yapcad.geom3d import volumeof


@pytest.fixture(autouse=True)
def _fresh_cache():
    # The module cache is process-global; keep tests independent of order.
    clear_module_cache()
    yield
    clear_module_cache()


def tiny(name, command, size, uses=""):
    """A one-command module, for fixtures where the body does not matter."""
    header = f"module {name}\n"
    if uses:
        header += f"use {uses}\n"
    return (f"{header}\ncommand {command}() -> solid:\n"
            f"    emit box({size}, {size}, {size})\n")


def write(directory, name, text):
    path = directory / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def run_file(path, command, parameters=None):
    return compile_and_run(
        path.read_text(), command, parameters or {}, source_path=str(path)
    )


def check_file(path):
    module = parse(tokenize(path.read_text()), source=path.read_text())
    return check(module, source_dir=str(path.parent))


LIB = """module lib

command INNER(r: float) -> solid:
    emit cylinder(r, 2.0)

command OUTER(r: float) -> solid:
    emit INNER(r)
"""


# Two modules that both define BOLT. mod_b also calls its own BOLT, which is
# where keeping the first definition silently used to produce mod_a's
# geometry with no diagnostic at all.
MOD_A = """module mod_a

command BOLT() -> solid:
    emit cylinder(1.0, 1.0)
"""

MOD_B = """module mod_b

command BOLT() -> solid:
    emit cylinder(50.0, 1.0)

command B_USES_BOLT() -> solid:
    emit BOLT()
"""


# ---------------------------------------------------------------- resolution

def test_unaliased_use_merges_commands(tmp_path):
    write(tmp_path, "shapes.dsl", """module shapes

command CUBE(s: float) -> solid:
    emit box(s, s, s)
""")
    main = write(tmp_path, "main.dsl", """module main
use shapes

command PART() -> solid:
    emit CUBE(3.0)
""")
    result = run_file(main, "PART")
    assert result.success, result.error_message
    assert volumeof(result.geometry) == pytest.approx(27.0)


def test_dotted_module_path_resolves_a_subdirectory(tmp_path):
    write(tmp_path, "pkg/inner.dsl", """module inner

command INNER_CUBE() -> solid:
    emit box(4.0, 4.0, 4.0)
""")
    main = write(tmp_path, "main.dsl", """module main
use pkg.inner

command PART() -> solid:
    emit INNER_CUBE()
""")
    result = run_file(main, "PART")
    assert result.success, result.error_message
    assert volumeof(result.geometry) == pytest.approx(64.0)


def test_nested_use_chains(tmp_path):
    write(tmp_path, "base.dsl", """module base

command base_unit(s: float) -> solid:
    emit box(s, s, s)
""")
    write(tmp_path, "mid.dsl", """module mid
use base

command MID_PART() -> solid:
    emit base_unit(2.0)
""")
    main = write(tmp_path, "main.dsl", """module main
use mid

command PART() -> solid:
    emit MID_PART()
""")
    result = run_file(main, "PART")
    assert result.success, result.error_message
    assert volumeof(result.geometry) == pytest.approx(8.0)


# ----------------------------------------------------------------- namespace

def test_an_unaliased_import_also_offers_a_namespace(tmp_path):
    """The bare name is a convenience; the namespace is always there."""
    write(tmp_path, "shapes.dsl", """module shapes

command CUBE(s: float) -> solid:
    emit box(s, s, s)
""")
    main = write(tmp_path, "main.dsl", """module main
use shapes

command PART() -> solid:
    emit union(CUBE(3.0), shapes.CUBE(2.0))
""")
    result = run_file(main, "PART")
    assert result.success, result.error_message
    assert volumeof(result.geometry) == pytest.approx(27.0)


def test_a_dotted_path_is_named_by_its_last_component(tmp_path):
    write(tmp_path, "pkg/inner.dsl", """module inner

command INNER_CUBE() -> solid:
    emit box(4.0, 4.0, 4.0)
""")
    main = write(tmp_path, "main.dsl", """module main
use pkg.inner

command PART() -> solid:
    emit inner.INNER_CUBE()
""")
    result = run_file(main, "PART")
    assert result.success, result.error_message
    assert volumeof(result.geometry) == pytest.approx(64.0)


def test_an_imported_command_can_be_the_entry_point(tmp_path):
    """``dsl run main.dsl CUBE`` asks for a command by name, and an
    unaliased import offers that name like any other."""
    write(tmp_path, "shapes.dsl", """module shapes

command CUBE() -> solid:
    emit box(3.0, 3.0, 3.0)
""")
    main = write(tmp_path, "main.dsl", """module main
use shapes

command PART() -> solid:
    emit CUBE()
""")
    result = run_file(main, "CUBE")
    assert result.success, result.error_message
    assert volumeof(result.geometry) == pytest.approx(27.0)


def test_alias_puts_commands_behind_a_dot(tmp_path):
    write(tmp_path, "lib.dsl", LIB)
    main = write(tmp_path, "main.dsl", """module main
use lib as L

command PART() -> solid:
    emit L.INNER(1.0)
""")
    result = run_file(main, "PART")
    assert result.success, result.error_message
    # a faceted cylinder falls slightly short of the ideal 2*pi
    assert volumeof(result.geometry) == pytest.approx(6.2832, rel=1e-2)


def test_an_aliased_library_can_call_its_own_commands(tmp_path):
    """OUTER calls INNER. Renaming the imported commands without rewriting
    the calls between them left OUTER reaching for a name that no longer
    existed, so aliasing broke every library of more than one command."""
    write(tmp_path, "lib.dsl", LIB)
    main = write(tmp_path, "main.dsl", """module main
use lib as L

command PART() -> solid:
    emit L.OUTER(5.0)
""")
    result = run_file(main, "PART")
    assert result.success, result.error_message
    assert volumeof(result.geometry) == pytest.approx(157.08, rel=1e-3)


def test_an_alias_does_not_also_expose_the_bare_name(tmp_path):
    write(tmp_path, "lib.dsl", LIB)
    main = write(tmp_path, "main.dsl", """module main
use lib as L

command PART() -> solid:
    emit OUTER(5.0)
""")
    result = run_file(main, "PART")
    assert not result.success
    assert "OUTER" in result.error_message


def test_a_variable_wins_over_an_import_namespace(tmp_path):
    """An ordinary method call on a value must not be captured by an import
    that happens to be aliased to the same name.

    The library below defines a command called ``translate`` and is aliased
    to ``s``, while ``s`` is also a local solid. ``s.translate(...)`` has to
    be read as a method on the solid. It is asserted through the diagnostic
    because the checker knows a signature for ``solid.translate`` that takes
    a vector: a complaint about that argument proves the method path was
    taken, where capture by the import would have type-checked cleanly
    against ``translate(r: float)`` and quietly produced a cylinder.
    """
    write(tmp_path, "lib.dsl", """module lib

command translate(r: float) -> solid:
    emit cylinder(r, 1.0)
""")
    main = write(tmp_path, "main.dsl", """module main
use lib as s

command PART() -> solid:
    let s: solid = box(2.0, 2.0, 2.0)
    emit s.translate(1.0, 0.0, 0.0)
""")
    result = run_file(main, "PART")
    assert not result.success
    assert "Method argument" in result.error_message


def test_a_namespaced_call_resolves_without_a_shadowing_variable(tmp_path):
    """The other half of the precedence rule: with no variable in the way,
    the same shape of call reaches the imported command."""
    write(tmp_path, "lib.dsl", """module lib

command translate(r: float) -> solid:
    emit cylinder(r, 1.0)
""")
    main = write(tmp_path, "main.dsl", """module main
use lib as s

command PART() -> solid:
    emit s.translate(2.0)
""")
    result = run_file(main, "PART")
    assert result.success, result.error_message
    assert volumeof(result.geometry) == pytest.approx(12.566, rel=1e-2)


# ----------------------------------------------------------------- conflicts

def test_the_same_module_by_two_routes_is_not_a_conflict(tmp_path):
    """A diamond: left and right both use base, top uses both."""
    write(tmp_path, "base.dsl", """module base

command base_unit(s: float) -> solid:
    emit box(s, s, s)
""")
    write(tmp_path, "left.dsl", """module left
use base

command LEFT_P() -> solid:
    emit base_unit(1.0)
""")
    write(tmp_path, "right.dsl", """module right
use base

command RIGHT_P() -> solid:
    emit base_unit(2.0)
""")
    main = write(tmp_path, "top.dsl", """module top
use left
use right

command PART() -> solid:
    emit union(LEFT_P(), RIGHT_P())
""")
    result = run_file(main, "PART")
    assert result.success, result.error_message
    assert volumeof(result.geometry) == pytest.approx(8.0)


def test_a_shared_command_name_is_only_a_problem_where_it_is_called(
    tmp_path,
):
    """Importing two modules that happen to share a command name costs
    nothing on its own. mod_b's own call to BOLT has to reach mod_b's
    BOLT -- radius 50, not mod_a's radius 1."""
    write(tmp_path, "mod_a.dsl", MOD_A)
    write(tmp_path, "mod_b.dsl", MOD_B)
    main = write(tmp_path, "main.dsl", """module main
use mod_a
use mod_b

command PART() -> solid:
    emit B_USES_BOLT()
""")
    result = run_file(main, "PART")
    assert result.success, result.error_message
    assert volumeof(result.geometry) == pytest.approx(7853.98, rel=1e-3)


def test_calling_a_name_two_imports_offer_is_an_error(tmp_path):
    write(tmp_path, "mod_a.dsl", MOD_A)
    write(tmp_path, "mod_b.dsl", MOD_B)
    main = write(tmp_path, "main.dsl", """module main
use mod_a
use mod_b

command PART() -> solid:
    emit BOLT()
""")
    result = run_file(main, "PART")
    assert not result.success
    assert "Ambiguous command 'BOLT'" in result.error_message
    assert "mod_a.BOLT" in result.error_message
    assert "mod_b.BOLT" in result.error_message


def test_an_ambiguous_call_is_reported_by_check(tmp_path):
    """Statically, at the call site -- not as a failure of the import."""
    write(tmp_path, "mod_a.dsl", MOD_A)
    write(tmp_path, "mod_b.dsl", MOD_B)
    main = write(tmp_path, "main.dsl", """module main
use mod_a
use mod_b

command PART() -> solid:
    emit BOLT()
""")
    result = check_file(main)
    assert result.has_errors
    ambiguous = [d for d in result.diagnostics if d.code == "E291"]
    assert len(ambiguous) == 1
    assert ambiguous[0].span.start.line == 6


@pytest.mark.parametrize("namespace,volume", [
    ("mod_a", 3.1416),
    ("mod_b", 7853.98),
])
def test_an_ambiguous_name_resolves_through_its_namespace(
    tmp_path, namespace, volume
):
    write(tmp_path, "mod_a.dsl", MOD_A)
    write(tmp_path, "mod_b.dsl", MOD_B)
    main = write(tmp_path, "main.dsl", f"""module main
use mod_a
use mod_b

command PART() -> solid:
    emit {namespace}.BOLT()
""")
    result = run_file(main, "PART")
    assert result.success, result.error_message
    assert volumeof(result.geometry) == pytest.approx(volume, rel=1e-2)


def test_an_unambiguous_name_is_unaffected_by_an_ambiguous_one(tmp_path):
    """The point of deferring the error: two libraries sharing one command
    name does not make their other commands unreachable."""
    write(tmp_path, "shapes.dsl", """module shapes

command CUBE(s: float) -> solid:
    emit box(s, s, s)
""")
    write(tmp_path, "ice.dsl", """module ice

command CUBE(s: float) -> solid:
    emit box(s, s, s * 2.0)

command PLATE(s: float) -> solid:
    emit box(s, s, 1.0)
""")
    main = write(tmp_path, "main.dsl", """module main
use shapes
use ice

command PART() -> solid:
    emit PLATE(3.0)
""")
    result = run_file(main, "PART")
    assert result.success, result.error_message
    assert volumeof(result.geometry) == pytest.approx(9.0)


def test_a_builtin_wins_over_an_ambiguous_import(tmp_path):
    """Both modules define a command called ``sphere``. The checker resolves
    built-ins before imported commands, so the interpreter has to as well:
    an ambiguity between two imports is no reason for the two layers to
    disagree about what a call means."""
    body = """module {name}

command sphere(r: float) -> solid:
    emit box(r, r, r)
"""
    write(tmp_path, "one.dsl", body.format(name="one"))
    write(tmp_path, "two.dsl", body.format(name="two"))
    main = write(tmp_path, "main.dsl", """module main
use one
use two

command PART() -> solid:
    emit sphere(1.0)
""")
    assert not check_file(main).has_errors
    result = run_file(main, "PART")
    assert result.success, result.error_message
    # the built-in sphere, not either module's box
    assert volumeof(result.geometry) == pytest.approx(4.19, rel=0.2)


def test_a_local_definition_wins_over_an_import_and_says_so(tmp_path):
    """Overriding one command of a library is useful; doing it by accident
    is not, so the local definition wins and the checker warns."""
    write(tmp_path, "shapes.dsl", """module shapes

command CUBE(s: float) -> solid:
    emit box(s, s, s)
""")
    main = write(tmp_path, "main.dsl", """module main
use shapes

command CUBE(s: float) -> solid:
    emit box(s, s, 1.0)

command PART() -> solid:
    emit CUBE(3.0)
""")
    result = run_file(main, "PART")
    assert result.success, result.error_message
    assert volumeof(result.geometry) == pytest.approx(9.0)

    checked = check_file(main)
    assert not checked.has_errors
    shadowed = [d for d in checked.diagnostics if d.code == "W290"]
    assert len(shadowed) == 1
    assert "shapes.CUBE" in shadowed[0].message
    # pointing at the local definition, which is what wins
    assert shadowed[0].span.start.line == 4

    # The library's own command is still reachable.
    main.write_text(main.read_text().replace(
        "emit CUBE(3.0)", "emit shapes.CUBE(3.0)"
    ))
    result = run_file(main, "PART")
    assert result.success, result.error_message
    assert volumeof(result.geometry) == pytest.approx(27.0)


def test_two_module_paths_ending_in_the_same_name_share_a_namespace(
    tmp_path,
):
    """``a.shapes`` and ``b.shapes`` are both ``shapes``, so a command they
    both define cannot be qualified apart. Say that, rather than reporting
    an unknown function."""
    write(tmp_path, "a/shapes.dsl", tiny("shapes", "CUBE", 1.0))
    write(tmp_path, "b/shapes.dsl", tiny("shapes", "CUBE", 4.0))
    main = write(tmp_path, "main.dsl", """module main
use a.shapes
use b.shapes

command PART() -> solid:
    emit shapes.CUBE()
""")
    result = run_file(main, "PART")
    assert not result.success
    assert "shares the namespace 'shapes'" in result.error_message
    assert "alias" in result.error_message.lower()

    main.write_text(main.read_text().replace(
        "use b.shapes", "use b.shapes as b_shapes"
    ).replace("emit shapes.CUBE()", "emit b_shapes.CUBE()"))
    result = run_file(main, "PART")
    assert result.success, result.error_message
    assert volumeof(result.geometry) == pytest.approx(64.0)


def test_aliasing_resolves_a_name_conflict(tmp_path):
    write(tmp_path, "mod_a.dsl", """module mod_a

command BOLT() -> solid:
    emit cylinder(1.0, 1.0)
""")
    write(tmp_path, "mod_b.dsl", """module mod_b

command BOLT() -> solid:
    emit cylinder(50.0, 1.0)

command B_USES_BOLT() -> solid:
    emit BOLT()
""")
    main = write(tmp_path, "main.dsl", """module main
use mod_a
use mod_b as b

command PART() -> solid:
    emit b.B_USES_BOLT()
""")
    result = run_file(main, "PART")
    assert result.success, result.error_message
    # b's own BOLT, radius 50 -- not mod_a's radius 1.
    assert volumeof(result.geometry) == pytest.approx(7853.98, rel=1e-3)


def test_two_imports_cannot_share_an_alias(tmp_path):
    write(tmp_path, "one.dsl", tiny("one", "A", 1.0))
    write(tmp_path, "two.dsl", tiny("two", "B", 2.0))
    main = write(tmp_path, "main.dsl", """module main
use one as x
use two as x

command PART() -> solid:
    emit x.A()
""")
    result = run_file(main, "PART")
    assert not result.success
    assert "alias" in result.error_message.lower()


# ------------------------------------------------------------------- failures

def test_a_circular_import_is_reported(tmp_path):
    write(tmp_path, "a.dsl", tiny("a", "A", 1.0, uses="b"))
    write(tmp_path, "b.dsl", tiny("b", "B", 2.0, uses="a"))
    main = tmp_path / "a.dsl"
    result = run_file(main, "A")
    assert not result.success
    assert "circular" in result.error_message.lower()


def test_a_missing_module_is_reported(tmp_path):
    main = write(tmp_path, "main.dsl", """module main
use nope

command PART() -> solid:
    emit box(1.0, 1.0, 1.0)
""")
    result = run_file(main, "PART")
    assert not result.success
    assert "nope" in result.error_message


@pytest.mark.parametrize("source,needle", [
    (tiny("main", "PART", 1.0, uses="nope"), "nope"),
    (tiny("main", "PART", 1.0, uses="a"), "circular"),
])
def test_check_reports_a_broken_import(tmp_path, source, needle):
    """``dsl check`` used to swallow the import failure and type-check the
    unresolved module, so a file that could not run reported no errors."""
    write(tmp_path, "a.dsl", tiny("a", "A", 1.0, uses="main"))
    main = write(tmp_path, "main.dsl", source)
    result = check_file(main)
    assert result.has_errors
    assert any(needle in d.message.lower() or needle in d.message
               for d in result.diagnostics)


def test_an_import_may_not_escape_the_source_directory(tmp_path):
    write(tmp_path, "secret.dsl", tiny("secret", "S", 7.0))
    main = write(tmp_path / "sub", "main.dsl", """module main
use ...secret

command PART() -> solid:
    emit S()
""")
    result = run_file(main, "PART")
    assert not result.success
    assert "outside the source directory" in result.error_message


def test_export_use_says_it_is_unsupported(tmp_path):
    write(tmp_path, "lib.dsl", LIB)
    main = write(tmp_path, "main.dsl", """module main
export use lib

command PART() -> solid:
    emit OUTER(2.0)
""")
    result = run_file(main, "PART")
    assert not result.success
    assert "export use" in result.error_message
    assert "not supported" in result.error_message


def test_imports_need_a_source_path(tmp_path):
    """Without one there is nothing to resolve against, and the message
    should say so rather than report an unknown function."""
    write(tmp_path, "lib.dsl", LIB)
    main = write(tmp_path, "main.dsl", """module main
use lib

command PART() -> solid:
    emit OUTER(1.0)
""")
    result = compile_and_run(main.read_text(), "PART", {})
    assert not result.success
    assert "source directory" in result.error_message


# ---------------------------------------------------------------------- cache

def test_editing_an_imported_module_is_picked_up(tmp_path):
    """The cache was keyed by path alone, so a long-lived process -- the
    viewer, a pytest session, a build tool -- kept rebuilding from the
    version of the module it first saw."""
    def shapes(size):
        write(tmp_path, "shapes.dsl", f"""module shapes

command CUBE() -> solid:
    emit box({size}, {size}, {size})
""")
    main = write(tmp_path, "main.dsl", """module main
use shapes

command PART() -> solid:
    emit CUBE()
""")
    shapes(1.0)
    first = run_file(main, "PART")
    assert first.success, first.error_message
    assert volumeof(first.geometry) == pytest.approx(1.0)

    shapes(4.0)
    second = run_file(main, "PART")
    assert second.success, second.error_message
    assert volumeof(second.geometry) == pytest.approx(64.0)


def test_resolve_uses_reports_failure_instead_of_raising(tmp_path):
    main = write(tmp_path, "main.dsl", """module main
use nope

command PART() -> solid:
    emit box(1.0, 1.0, 1.0)
""")
    module = parse(tokenize(main.read_text()), source=main.read_text())
    merged = resolve_uses(module, str(tmp_path))
    assert not merged.success
    assert merged.merged_module is None
    assert "nope" in merged.error_message


# ------------------------------------------------------------------ packaging

def test_a_design_with_imports_can_be_packaged(tmp_path):
    """package_from_dsl had no source_path, so any design using an import
    failed to package even though it ran -- which is the path
    tools/build_ycpkg.py depends on."""
    from yapcad.dsl.packaging import package_from_dsl

    write(tmp_path, "lib.dsl", LIB)
    main = write(tmp_path, "main.dsl", """module main
use lib as L

command PART() -> solid:
    emit L.OUTER(3.0)
""")
    result = package_from_dsl(
        main.read_text(), "PART", {}, tmp_path / "pkg",
        name="imports", version="0.1.0", source_path=str(main),
    )
    assert result.success, result.error_message
