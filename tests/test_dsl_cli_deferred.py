"""``--defer-mesh`` on ``python -m yapcad.dsl run`` and its mesh report.

An SDF design destined for STL is meshed exactly once, at the end, and
the run reports what that single mesh looks like: triangle count,
watertightness, connected components and bounding box.  The default is
deferred when the representation is SDF and the output is STL; the
report goes to stderr and ``--require-watertight`` can fail the run.
"""
import subprocess
import sys

from yapcad.dsl.__main__ import main
from yapcad.mesh_check import check_mesh, MeshCheck

BOX_DSL = (
    "module defercli\n"
    "command BLOCK() -> solid:\n"
    "    emit box(10.0, 5.0, 2.0)\n"
)


def _dsl(tmp_path):
    dsl = tmp_path / "p.dsl"
    dsl.write_text(BOX_DSL)
    return dsl


def _run(*args):
    return subprocess.run(
        [sys.executable, "-W", "ignore", "-m", "yapcad.dsl", *args],
        capture_output=True, text=True,
    )


def test_stl_output_of_an_sdf_design_defers_by_default(tmp_path):
    dsl = _dsl(tmp_path)
    out = tmp_path / "block.stl"
    done = _run("run", str(dsl), "BLOCK", "--representation", "sdf",
                "--sdf-cell", "2.0", "-o", str(out))
    assert done.returncode == 0, done.stderr
    assert out.exists()
    assert out.stat().st_size > 84  # a binary STL carrying triangles
    # The one-mesh report lands on stderr.
    assert "triangles" in done.stderr
    assert "watertight: yes" in done.stderr
    assert "components: 1" in done.stderr
    assert "bbox" in done.stderr


def test_defer_mesh_off_keeps_the_eager_build(tmp_path):
    dsl = _dsl(tmp_path)
    out = tmp_path / "block.stl"
    done = _run("run", str(dsl), "BLOCK", "--representation", "sdf",
                "--sdf-cell", "2.0", "--no-defer-mesh", "-o", str(out))
    assert done.returncode == 0, done.stderr
    assert out.exists()
    # Eager builds answer the volume up front; the deferred report is absent.
    assert "volume" in done.stdout
    assert "triangles" not in done.stderr


def test_deferred_without_output_builds_no_mesh(tmp_path):
    dsl = _dsl(tmp_path)
    done = _run("run", str(dsl), "BLOCK", "--representation", "sdf",
                "--defer-mesh")
    assert done.returncode == 0, done.stderr
    assert "deferred field" in done.stdout
    # meshing is never asked for, so it is never paid for
    assert "triangles" not in done.stderr


def test_defer_mesh_stays_off_for_mesh_representations(tmp_path):
    dsl = _dsl(tmp_path)
    out = tmp_path / "block.stl"
    done = _run("run", str(dsl), "BLOCK", "-o", str(out))
    assert done.returncode == 0, done.stderr
    assert "triangles" not in done.stderr


def test_require_watertight_passes_a_closed_mesh(tmp_path):
    dsl = _dsl(tmp_path)
    out = tmp_path / "block.stl"
    done = _run("run", str(dsl), "BLOCK", "--representation", "sdf",
                "--sdf-cell", "2.0", "--require-watertight",
                "-o", str(out))
    assert done.returncode == 0, done.stderr
    assert out.exists()


def test_require_watertight_fails_an_open_mesh(tmp_path, monkeypatch,
                                                capsys):
    # The failure path is about the verdict rather than about making a
    # genuinely open mesh at cell 2.0, so the checker's report is forced.
    real = check_mesh

    def fake(target):
        report = real(target)
        return MeshCheck(triangles=report.triangles, watertight=False,
                         components=report.components, bounds=report.bounds)

    monkeypatch.setattr("yapcad.mesh_check.check_mesh", fake)
    dsl = _dsl(tmp_path)
    out = tmp_path / "block.stl"
    code = main(["run", str(dsl), "BLOCK", "--representation", "sdf",
                 "--sdf-cell", "2.0", "--require-watertight",
                 "-o", str(out)])
    assert code == 1
    # refused before the file was written
    assert not out.exists()
    assert "not watertight" in capsys.readouterr().err


def test_report_of_the_deferred_single_mesh():
    from yapcad.dsl import compile_and_run
    from yapcad.sdf.deferred import realize
    result = compile_and_run(BOX_DSL, "BLOCK", {}, representation="sdf",
                             sdf_cell_mm=2.0, sdf_mesh=False)
    assert result.success
    geometry = result.geometry
    if isinstance(geometry, list) and geometry and geometry[0] != "solid":
        geometry = geometry[0]
    realize(geometry)
    report = check_mesh(geometry)
    assert report.watertight is True
    assert report.components == 1
    assert report.bounds is not None
    assert "watertight: yes" in report.summary()
    # the summary carries all four numbers the CLI promised
    for word in ("triangles", "watertight", "components", "bbox"):
        assert word in report.summary()


def test_check_mesh_sees_two_bodies_as_two_components():
    import trimesh
    pair = trimesh.util.concatenate([
        trimesh.creation.box((2.0, 2.0, 2.0)),
        trimesh.creation.box((2.0, 2.0, 2.0)).apply_translation((10.0, 0, 0)),
    ])
    report = check_mesh(pair)
    assert report.watertight is True
    assert report.components == 2


def test_check_mesh_sees_an_open_mesh():
    import numpy as np
    import trimesh
    open_mesh = trimesh.Trimesh(
        vertices=np.array([[0.0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]]),
        faces=np.array([[0, 1, 2], [0, 2, 3]]))
    report = check_mesh(open_mesh)
    assert report.watertight is False
