"""Pure-Python contract tests for OpenCASCADE tessellation controls."""

import pytest

import yapcad.brep as brep_module


def _uninitialized_solid():
    solid = object.__new__(brep_module.BrepSolid)
    solid._shape = object()
    return solid


def test_brep_tessellate_forwards_explicit_meshing_controls(monkeypatch):
    captured = {}

    class MesherReached(RuntimeError):
        pass

    def capture(*args):
        captured["args"] = args
        raise MesherReached

    solid = _uninitialized_solid()
    monkeypatch.setattr(brep_module, "require_occ", lambda: None)
    monkeypatch.setattr(brep_module, "BRepMesh_IncrementalMesh", capture)

    with pytest.raises(MesherReached):
        solid.tessellate(
            0.125,
            angular_deflection=0.25,
            relative=True,
            parallel=False,
        )
    assert captured["args"] == (solid._shape, 0.125, True, 0.25, False)


@pytest.mark.parametrize("name,value", [
    ("deflection", 0),
    ("deflection", float("nan")),
    ("angular_deflection", -0.1),
    ("angular_deflection", float("inf")),
])
def test_brep_tessellate_rejects_invalid_meshing_controls(monkeypatch, name, value):
    solid = _uninitialized_solid()
    monkeypatch.setattr(brep_module, "require_occ", lambda: None)
    with pytest.raises(ValueError, match="positive finite"):
        solid.tessellate(**{name: value})
