"""A call that ends a block still runs.

The parser lifts a block's trailing bare expression into the block's
``final_expression``.  The interpreter used to execute only the block's
``statements``, so a call written for its effect as the *last* line of an
``if`` or ``for`` body -- ``add_part(...)`` ending a loop, a recursive
command ending an ``if`` -- was dropped without a word.
"""

from yapcad.dsl import compile_and_run


def run(source, command="TOP"):
    result = compile_and_run(source, command, {})
    assert result.success, result.error_message
    return result


COUNT_PARTS = '''
module trailing

command add(a: assembly, name: string) -> assembly:
    add_part(a, box(1.0, 1.0, 1.0), name)
    emit a

{body}

command TOP() -> string:
    let a: assembly = assembly("t")
    add_part(a, box(1.0, 1.0, 1.0), "root")
    fill(a)
    emit assembly_report(a)
'''


def report(body):
    return run(COUNT_PARTS.format(body=body)).emit_result.data


def test_a_call_ending_a_for_body_runs():
    text = report('''
command fill(a: assembly) -> assembly:
    for name in ["p0", "p1", "p2"]:
        add(a, name)
    emit a
''')
    assert "Parts: 4" in text
    for name in ("p0", "p1", "p2"):
        assert name in text


def test_a_call_ending_an_if_body_runs():
    text = report('''
command fill(a: assembly) -> assembly:
    if true:
        add(a, "then")
    if false:
        add(a, "never")
    elif true:
        add(a, "elif")
    if false:
        add(a, "never")
    else:
        add(a, "else")
    emit a
''')
    assert "Parts: 4" in text
    for name in ("then", "elif", "else"):
        assert name in text
    assert "never" not in text


def test_a_recursive_call_ending_an_if_body_runs():
    text = report('''
command walk(a: assembly, names: list<string>, i: int) -> assembly:
    add(a, names[i])
    if i < 2:
        walk(a, names, i + 1)
    emit a

command fill(a: assembly) -> assembly:
    emit walk(a, ["one", "two", "three"], 0)
''')
    assert "Parts: 4" in text
    for name in ("one", "two", "three"):
        assert name in text
