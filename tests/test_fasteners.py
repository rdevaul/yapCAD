import math
import os
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

os.environ.setdefault("YAPCAD_BOOLEAN_ENGINE", "native")

from yapcad.fasteners import (
    HexCapScrewSpec,
    build_hex_cap_screw,
    metric_hex_cap_catalog,
    metric_hex_cap_screw,
    unified_hex_cap_catalog,
    unified_hex_cap_screw,
)
from yapcad.geom3d import issolid, solidbbox
from yapcad.metadata import get_solid_metadata
from yapcad.threadgen import metric_profile


def test_build_hex_cap_screw_basic():
    profile = metric_profile(6.0, 1.0)
    spec = HexCapScrewSpec(
        diameter=6.0,
        thread_length=12.0,
        shank_length=20.0,
        head_height=6.0,
        head_flat_diameter=10.0,
        washer_thickness=1.0,
    )
    screw = build_hex_cap_screw(profile, spec)
    assert issolid(screw)
    box = solidbbox(screw)
    assert math.isclose(box[0][2], 0.0, abs_tol=1e-6)
    assert box[1][2] > spec.shank_length + spec.head_height
    meta = get_solid_metadata(screw, create=False)
    assert "fastener" in meta.get("tags", [])


def test_metric_helper_defaults_thread_length():
    screw = metric_hex_cap_screw("M8", length=30.0)
    box = solidbbox(screw)
    assert box[1][2] > 30.0
    meta = get_solid_metadata(screw, create=False)
    washer_info = meta["hex_cap_screw"]
    assert math.isclose(washer_info["washer_diameter"], washer_info["head_flat"] * 0.95, rel_tol=0.0, abs_tol=1e-6)
    assert math.isclose(
        washer_info["washer_thickness"],
        metric_hex_cap_catalog()["M8"]["washer_thickness"] / 2.0,
        rel_tol=0.0,
        abs_tol=1e-6,
    )


def test_unified_helper_accepts_inches():
    screw = unified_hex_cap_screw("1/4-20", length_in=1.0)
    box = solidbbox(screw)
    assert box[1][2] > 25.4
    meta = get_solid_metadata(screw, create=False)
    washer_info = meta["hex_cap_screw"]
    assert math.isclose(washer_info["washer_diameter"], washer_info["head_flat"] * 0.95, rel_tol=0.0, abs_tol=1e-6)


def test_tables_expose_washer_dimensions():
    metric_table = metric_hex_cap_catalog()
    assert math.isclose(metric_table["M8"]["washer_thickness"], 0.5, rel_tol=0.0, abs_tol=1e-6)
    
    unified_table = unified_hex_cap_catalog()
    assert "1/4-20" in unified_table
    assert unified_table["1/4-20"]["washer_thickness"] > 0
