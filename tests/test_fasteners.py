import math
import os
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

os.environ.setdefault("YAPCAD_BOOLEAN_ENGINE", "native")

from yapcad.fasteners import (
    HexCapScrewSpec,
    HexNutSpec,
    build_hex_cap_screw,
    build_hex_nut,
    metric_hex_cap_catalog,
    metric_hex_cap_screw,
    metric_hex_nut,
    metric_hex_nut_catalog,
    unified_hex_cap_catalog,
    unified_hex_cap_screw,
    unified_hex_nut,
    unified_hex_nut_catalog,
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


def test_build_hex_nut_basic():
    profile = metric_profile(8.0, 1.25, internal=True)
    spec = HexNutSpec(
        diameter=8.0,
        pitch=1.25,
        width_flat=13.0,
        thickness=6.0,
    )
    nut = build_hex_nut(profile, spec)
    assert issolid(nut)
    bbox = solidbbox(nut)
    assert math.isclose(bbox[0][2], 0.0, abs_tol=1e-6)
    assert math.isclose(bbox[1][2], spec.thickness, abs_tol=1e-6)
    meta = get_solid_metadata(nut, create=False)
    assert "hex_nut" in meta


def test_metric_hex_nut_helper():
    nut = metric_hex_nut("M8")
    bbox = solidbbox(nut)
    assert bbox[1][2] > 6.0
    table = metric_hex_nut_catalog()
    assert "M8" in table
    meta = get_solid_metadata(nut, create=False)
    assert math.isclose(meta["hex_nut"]["width_flat"], table["M8"]["width_flat"], rel_tol=0.0, abs_tol=1e-6)


def test_unified_hex_nut_helper():
    nut = unified_hex_nut("1/4-20")
    table = unified_hex_nut_catalog()
    assert "1/4-20" in table
    bbox = solidbbox(nut)
    assert bbox[1][2] > table["1/4-20"]["thickness"] - 0.1
