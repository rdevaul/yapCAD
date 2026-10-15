"""`python -m yapcad.dsl run -o out.dxf` on a solid is refused clearly.

DXF carries 2D geometry; the exporter draws curves, polylines and regions.
The CLI used to hand solids to it anyway and surface the drawable layer's
"bad argument to Drawable.draw(): solid".
"""
import subprocess
import sys


def _run(args):
    return subprocess.run(
        [sys.executable, "-W", "ignore", "-m", "yapcad.dsl", *args],
        capture_output=True, text=True,
    )


def _solid_dsl(tmp_path):
    dsl = tmp_path / "p.dsl"
    dsl.write_text(
        "module p\n"
        "command PLATE() -> solid:\n"
        "    emit box(10.0, 5.0, 2.0)\n"
    )
    return dsl


def test_dxf_export_of_a_solid_is_refused_with_a_useful_message(tmp_path):
    dsl = _solid_dsl(tmp_path)
    out = tmp_path / "p.dxf"
    done = _run(["run", str(dsl), "PLATE", "-o", str(out)])

    assert done.returncode == 1
    assert not out.exists()
    # Names the problem and the formats that do take a solid.
    assert "DXF export needs 2D geometry" in done.stderr
    assert ".stl" in done.stderr and ".step" in done.stderr
    # The old drawable-layer failure must not resurface.
    assert "Drawable.draw" not in done.stderr
    assert "Drawable.draw" not in done.stdout


def test_stl_export_of_the_same_solid_still_works(tmp_path):
    dsl = _solid_dsl(tmp_path)
    out = tmp_path / "p.stl"
    done = _run(["run", str(dsl), "PLATE", "-o", str(out)])

    assert done.returncode == 0, done.stderr
    assert out.exists() and out.stat().st_size > 0
