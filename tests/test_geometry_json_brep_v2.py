"""OCC acceptance tests for authoritative geometry JSON v0.2 BREP."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from yapcad.brep import brep_from_solid, has_brep_data, occ_available
from yapcad.dsl import compile_and_run
from yapcad.geom3d import solidbbox, volumeof
from yapcad.io.geometry_json import geometry_from_json, geometry_to_json
from yapcad.io.step import write_step_analytic


pytestmark = pytest.mark.skipif(
    not occ_available(), reason="pythonocc-core is unavailable"
)

BEVEL_GEAR_SOURCE = r'''
module geometry_json_brep

command TEST_GEAR() -> solid:
    emit straight_bevel_gear(
        teeth=12,
        mate_teeth=12,
        outer_module_mm=1.5,
        face_width_mm=5.0,
        shaft_angle_deg=90.0,
        pressure_angle_deg=20.0,
        backlash_mm=0.2,
        bore_diameter_mm=4.0,
        generation_type="spherical_involute",
    )
'''


def _gear():
    result = compile_and_run(BEVEL_GEAR_SOURCE, "TEST_GEAR", {})
    assert result.success, result.error_message
    return result.geometry


def test_bevel_gear_roundtrips_authoritative_brep(tmp_path: Path):
    original = _gear()
    original_bounds = np.asarray(solidbbox(original), dtype=float)
    original_volume = volumeof(original)

    document = json.loads(json.dumps(geometry_to_json([original], units="mm")))
    solid_entry = next(e for e in document["entities"] if e["type"] == "solid")
    assert solid_entry["representations"]["authoritative"] == "brep"
    assert "brep" not in solid_entry["metadata"]

    restored = geometry_from_json(document)[0]
    brep = brep_from_solid(restored)
    assert has_brep_data(restored)
    assert brep is not None
    from OCC.Core.BRepCheck import BRepCheck_Analyzer
    assert BRepCheck_Analyzer(brep.shape).IsValid()
    np.testing.assert_allclose(solidbbox(restored), original_bounds, atol=1e-9)
    assert volumeof(restored) == pytest.approx(original_volume, rel=1e-9)

    output = tmp_path / "roundtripped-miter-gear.step"
    assert write_step_analytic(
        restored, str(output), fallback_to_faceted=False,
    ) is True
    assert output.stat().st_size > 0


def test_legacy_v01_metadata_brep_still_rehydrates():
    document = geometry_to_json([_gear()])
    document["schema"] = "yapcad-geometry-json-v0.1"
    solid_entry = next(e for e in document["entities"] if e["type"] == "solid")
    explicit = solid_entry.pop("representations")["brep"]
    solid_entry["metadata"]["brep"] = {
        "encoding": "brep-ascii-base64",
        "data": explicit["payload"],
    }

    restored = geometry_from_json(document)[0]
    assert brep_from_solid(restored) is not None
