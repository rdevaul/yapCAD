"""Which backend meshes a field, and how a mesh says which one did.

numpy is the reference.  The compiled backend makes the same triangles with
vertices that agree to rounding, not bit for bit, and yapCAD signs
packages -- so the choice is recorded in the mesh's meshing record, and
``YAPCAD_SDF_BACKEND=numpy`` forces the reference.
"""

import json
import pathlib

import numpy as np
import pytest

from yapcad import sdf
from yapcad.construction import get_construction
from yapcad.geom3d import issolidclosed, volumeof
from yapcad.sdf import tape
from yapcad.sdf.evaluate import mesh_backend
from yapcad.sdf.node import meshing_from_construction

needs_numba = pytest.mark.skipif(not tape.numba_available(),
                                 reason="numba is not installed")

PART = sdf.subtract(
    sdf.smooth_union(sdf.rounded_box((20.0, 14.0, 8.0), 1.5),
                     sdf.translate(sdf.cylinder(4.0, 12.0), (5.0, 0.0, 4.0)),
                     1.0),
    sdf.translate(sdf.cylinder(2.0, 30.0), (-5.0, 0.0, 0.0)))

#: A part the tape cannot express, so it is always meshed with numpy.
PROFILE = sdf.extrude(sdf.polygon([(0, 0), (8, 0), (8, 5), (0, 5)]), 4.0)


def record(solid):
    return meshing_from_construction(get_construction(solid))


def arrays(solid):
    surface = solid[1][0]
    return (np.array([p[:3] for p in surface[1]], dtype=float),
            np.array(surface[3]))


@pytest.fixture
def setting(monkeypatch):
    def use(value):
        if value is None:
            monkeypatch.delenv("YAPCAD_SDF_BACKEND", raising=False)
        else:
            monkeypatch.setenv("YAPCAD_SDF_BACKEND", value)
    use(None)
    return use


def test_numpy_leaves_the_record_as_it_always_was(setting):
    """An absent key means numpy, so documents meshed with the reference
    backend do not change."""
    setting("numpy")
    assert record(sdf.to_solid(PART, resolution=32)) == {
        "method": "dual-contouring", "resolution": 32}
    assert "backend" not in record(
        sdf.to_solid(PART, resolution=32, backend="numpy"))


@needs_numba
def test_the_compiled_backend_is_the_default_and_says_so(setting):
    solid = sdf.to_solid(PART, resolution=32)
    assert record(solid)["backend"] == "numba"
    assert issolidclosed(solid)


@needs_numba
def test_the_setting_forces_the_reference_backend(setting):
    setting("numpy")
    forced = sdf.to_solid(PART, resolution=32)
    explicit = sdf.to_solid(PART, resolution=32, backend="numpy")
    assert "backend" not in record(forced)
    for a, b in zip(arrays(forced), arrays(explicit)):
        assert np.array_equal(a, b)


@needs_numba
def test_the_compiled_mesh_is_the_same_mesh_to_rounding(setting):
    reference = sdf.to_solid(PART, resolution=48, backend="numpy")
    compiled = sdf.to_solid(PART, resolution=48, backend="numba")
    (rv, rt), (cv, ct) = arrays(reference), arrays(compiled)
    assert np.array_equal(rt, ct)
    assert cv == pytest.approx(rv, abs=1e-8)
    assert volumeof(compiled) == pytest.approx(volumeof(reference),
                                               rel=1e-9)


@needs_numba
def test_an_argument_beats_the_setting(setting):
    setting("numpy")
    assert record(sdf.to_solid(PART, resolution=24,
                               backend="numba"))["backend"] == "numba"


def test_a_part_the_tape_cannot_express_is_meshed_with_numpy(setting):
    solid = sdf.to_solid(PROFILE, resolution=24)
    assert "backend" not in record(solid)
    assert mesh_backend(PROFILE) == "numpy"


def test_without_numba_everything_is_numpy(setting, monkeypatch):
    monkeypatch.setattr(tape, "numba_available", lambda: False)
    assert mesh_backend(PART) == "numpy"
    assert "backend" not in record(sdf.to_solid(PART, resolution=24))
    with pytest.raises(sdf.SdfError, match="numba"):
        sdf.to_solid(PART, resolution=24, backend="numba")


def test_the_gpu_backend_is_refused_for_meshing(setting):
    with pytest.raises(sdf.SdfError, match="single precision"):
        sdf.to_solid(PART, resolution=24, backend="wgpu")
    setting("wgpu")
    with pytest.raises(sdf.SdfError, match="single precision"):
        sdf.to_solid(PART, resolution=24)


def test_an_unknown_backend_is_refused(setting):
    setting("fastest")
    with pytest.raises(sdf.SdfError, match="fastest"):
        sdf.to_solid(PART, resolution=24)


@needs_numba
def test_a_deferred_part_records_the_backend_that_realized_it(setting):
    with sdf.deferring():
        part = sdf.to_solid(PART, resolution=32)
    assert "backend" not in record(part)        # a plan, not a mesh
    sdf.realize(part)
    assert record(part)["backend"] == "numba"


@needs_numba
def test_simplifying_records_its_backend_and_keeps_the_part(setting):
    uniform = sdf.to_solid(PART, resolution=48)
    compiled = sdf.simplify_solid(uniform)
    assert record(compiled)["simplify"]["backend"] == "numba"
    setting("numpy")
    reference = sdf.simplify_solid(sdf.to_solid(PART, resolution=48))
    assert "backend" not in record(reference)["simplify"]
    assert issolidclosed(compiled)
    assert len(arrays(compiled)[1]) < 0.5 * len(arrays(uniform)[1])
    # Collapses that sit on the tolerance can fall either way between
    # backends; both results are within it, so the parts agree closely.
    assert volumeof(compiled) == pytest.approx(volumeof(reference),
                                               rel=2e-3)
    assert len(arrays(compiled)[1]) == pytest.approx(
        len(arrays(reference)[1]), rel=0.1)


@needs_numba
def test_the_backend_survives_a_geometry_document(setting):
    from yapcad.io.geometry_json import geometry_from_json, geometry_to_json
    solid = sdf.to_solid(PART, resolution=24)
    doc = geometry_to_json([solid])
    entry = next(e for e in doc["entities"] if e["type"] == "solid")
    assert entry["representations"]["mesh"]["backend"] == "numba"

    jsonschema = pytest.importorskip("jsonschema")
    schema = json.loads((pathlib.Path(__file__).parent.parent / "docs"
                         / "schemas"
                         / "yapcad-geometry-json-v0.3.schema.json"
                         ).read_text())
    jsonschema.validate(json.loads(json.dumps(doc)), schema)

    back = [g for g in geometry_from_json(doc) if g and g[0] == "solid"][0]
    assert record(back)["backend"] == "numba"
