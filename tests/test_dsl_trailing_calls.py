"""A call that ends a block still runs, and ``str`` names what it adds.

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
    for i in range(3):
        add(a, str("p", i))
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
    """The shape of a branching recursion: one call continues, and a second,
    inside an ``if``, branches."""
    text = report('''
command branch(a: assembly, name: string, n: int) -> assembly:
    if n <= 0:
        emit a
    add(a, name)
    branch(a, str(name, "s"), n - 1)
    if n > 1:
        branch(a, str(name, "b"), n - 1)
    emit a

command fill(a: assembly) -> assembly:
    emit branch(a, "v", 3)
''')
    # v; vs, vb; vss, vsb, vbs, vbb
    assert "Parts: 8" in text
    for name in ("vsb", "vbs", "vbb"):
        assert name in text


def test_str_joins_the_text_of_its_arguments():
    source = '''
module naming
command TOP() -> string:
    emit str("leaf_", 2, "_", 0.5, "_", true)
'''
    assert run(source).emit_result.data == "leaf_2_0.5_true"
