"""`python -m yapcad.dsl list` prints defaults as written in the source."""
import subprocess
import sys


def test_list_prints_source_defaults(tmp_path):
    dsl = tmp_path / "q.dsl"
    dsl.write_text(
        "module q\n"
        'command BRACKET(width: float = 40.0, name: string = "plate", k: float = 2.0 * 3.0) -> solid:\n'
        "    emit box(width, 1.0, 1.0)\n"
    )
    out = subprocess.run(
        [sys.executable, "-W", "ignore", "-m", "yapcad.dsl", "list", str(dsl)],
        capture_output=True, text=True, check=True,
    ).stdout
    assert 'BRACKET(width: float = 40.0, name: string = "plate", k: float = 2.0 * 3.0) -> solid' in out
    assert "Literal(" not in out
