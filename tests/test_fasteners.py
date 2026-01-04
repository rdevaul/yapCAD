import math
import os
import pathlib
import sys

import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))


@pytest.fixture(autouse=True)
def use_native_boolean_engine():
    """Use native boolean engine for fastener tests, restore original after."""
    original = os.environ.get("YAPCAD_BOOLEAN_ENGINE")
    os.environ["YAPCAD_BOOLEAN_ENGINE"] = "native"
    yield
    if original is None:
        os.environ.pop("YAPCAD_BOOLEAN_ENGINE", None)
    else:
        os.environ["YAPCAD_BOOLEAN_ENGINE"] = original


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


# ============================================================================
# Tests for the new catalog-based fastener system
# ============================================================================

from yapcad.fasteners import (
    YAPCAD_FASTENER_DATA,
    THREAD_SERIES,
    load_catalog,
    list_available_sizes,
    list_tolerance_classes,
    get_thread_data,
    get_bolt_data,
    get_nut_data,
    clear_cache,
    metric_hex_bolt,
)


class TestCatalogLoading:
    """Test catalog loading functionality."""

    def setup_method(self):
        """Clear cache before each test."""
        clear_cache()

    def test_load_metric_coarse_catalog(self):
        """Test loading the bundled metric coarse catalog."""
        catalog = load_catalog("metric_coarse")
        assert "standard" in catalog
        assert "ISO 262" in catalog["standard"]  # May include other standards
        assert "sizes" in catalog
        assert len(catalog["sizes"]) > 0

    def test_load_unified_coarse_catalog(self):
        """Test loading the bundled unified coarse catalog."""
        catalog = load_catalog("unified_coarse")
        assert "standard" in catalog
        assert "ASME B1.1" in catalog["standard"]  # May include other standards
        assert "sizes" in catalog
        assert len(catalog["sizes"]) > 0

    def test_invalid_thread_series(self):
        """Test that invalid thread series raises ValueError."""
        with pytest.raises(ValueError, match="Unknown thread series"):
            load_catalog("invalid_series")

    def test_thread_series_constant(self):
        """Test that THREAD_SERIES contains expected entries."""
        assert "metric_coarse" in THREAD_SERIES
        assert "unified_coarse" in THREAD_SERIES


class TestMetricCatalogData:
    """Test metric catalog thread data retrieval."""

    def setup_method(self):
        """Clear cache before each test."""
        clear_cache()

    def test_list_metric_sizes(self):
        """Test listing available metric sizes."""
        sizes = list_available_sizes("metric_coarse")
        assert "M8" in sizes
        assert "M10" in sizes
        assert "M12" in sizes

    def test_list_metric_tolerance_classes(self):
        """Test listing available metric tolerance classes."""
        classes = list_tolerance_classes("metric_coarse")
        assert "6g" in classes  # Standard external
        assert "6H" in classes  # Standard internal

    def test_get_metric_thread_data(self):
        """Test retrieving metric thread data."""
        data = get_thread_data("metric_coarse", "M8", "6g")
        assert data["nominal_diameter"] == 8.0
        assert data["pitch"] == 1.25
        assert "tolerances" in data

    def test_metric_tolerance_values(self):
        """Test that metric tolerance values are present."""
        data = get_thread_data("metric_coarse", "M10", "6g")
        tol = data["tolerances"]
        # External thread should have es (upper deviation) and Td2 as keys
        # under the specific tolerance class
        assert "6g" in tol or "es" in tol  # Either nested or flat structure
        if "6g" in tol:
            # Nested structure: tolerances["6g"]["es"]
            assert "es" in tol["6g"]
            assert "Td2" in tol["6g"]
        else:
            # Flat structure: tolerances["es"]
            assert "es" in tol
            assert "Td2" in tol

    def test_invalid_metric_size(self):
        """Test that invalid metric size raises KeyError."""
        with pytest.raises(KeyError, match="M99"):
            get_thread_data("metric_coarse", "M99", "6g")


class TestUnifiedCatalogData:
    """Test unified catalog thread data retrieval."""

    def setup_method(self):
        """Clear cache before each test."""
        clear_cache()

    def test_list_unified_sizes(self):
        """Test listing available unified sizes."""
        sizes = list_available_sizes("unified_coarse")
        assert "1/4-20" in sizes
        assert "1/2-13" in sizes

    def test_get_unified_thread_data(self):
        """Test retrieving unified thread data."""
        data = get_thread_data("unified_coarse", "1/4-20", "2A")
        assert abs(data["nominal_diameter"] - 6.35) < 0.01  # 1/4" = 6.35mm
        assert abs(data["pitch"] - 1.27) < 0.01  # 20 TPI = 1.27mm


class TestCatalogBoltNutData:
    """Test bolt and nut data retrieval from catalog."""

    def setup_method(self):
        """Clear cache before each test."""
        clear_cache()

    def test_get_metric_bolt_data(self):
        """Test retrieving metric bolt data."""
        data = get_bolt_data("metric_coarse", "M10", "6g")
        assert "thread" in data
        assert "head" in data
        assert data["thread"]["nominal_diameter"] == 10.0
        assert "width_across_flats" in data["head"]
        assert "head_height" in data["head"]

    def test_get_metric_nut_data(self):
        """Test retrieving metric nut data."""
        data = get_nut_data("metric_coarse", "M8", "6H")
        assert "thread" in data
        assert "body" in data
        assert data["thread"]["nominal_diameter"] == 8.0
        assert "width_across_flats" in data["body"]
        assert "thickness" in data["body"]


class TestCatalogBasedFastenerGeneration:
    """Test catalog-based fastener solid generation."""

    def setup_method(self):
        """Clear cache before each test."""
        clear_cache()

    def test_metric_hex_bolt_from_catalog(self):
        """Test creating a metric hex bolt using the new API."""
        bolt = metric_hex_bolt("M8", 25.0)
        assert isinstance(bolt, list)
        assert issolid(bolt)

    def test_metric_hex_bolt_different_sizes(self):
        """Test creating metric bolts of different sizes via catalog."""
        for size in ["M6", "M8", "M10"]:
            bolt = metric_hex_bolt(size, 30.0)
            assert issolid(bolt)


class TestDSLFastenerBuiltins:
    """Test DSL integration for fastener builtins."""

    def setup_method(self):
        """Clear cache before each test."""
        clear_cache()

    def test_dsl_metric_bolt(self):
        """Test metric_hex_bolt DSL builtin."""
        from yapcad.dsl import tokenize, parse, check, Interpreter

        source = '''
module test_fastener

command MAKE_BOLT(size: string, length: float) -> solid:
    bolt: solid = metric_hex_bolt(size, length)
    emit bolt
'''
        tokens = tokenize(source)
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(
            module,
            'MAKE_BOLT',
            {'size': 'M8', 'length': 25.0},
            source=source
        )
        assert exec_result.success
        assert exec_result.geometry is not None

    def test_dsl_metric_nut(self):
        """Test metric_hex_nut DSL builtin."""
        from yapcad.dsl import tokenize, parse, check, Interpreter

        source = '''
module test_fastener

command MAKE_NUT(size: string) -> solid:
    nut: solid = metric_hex_nut(size)
    emit nut
'''
        tokens = tokenize(source)
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(
            module,
            'MAKE_NUT',
            {'size': 'M10'},
            source=source
        )
        assert exec_result.success
        assert exec_result.geometry is not None

    def test_dsl_unified_bolt(self):
        """Test unified_hex_bolt DSL builtin."""
        from yapcad.dsl import tokenize, parse, check, Interpreter

        source = '''
module test_fastener

command MAKE_UNC_BOLT(size: string, length: float) -> solid:
    bolt: solid = unified_hex_bolt(size, length)
    emit bolt
'''
        tokens = tokenize(source)
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(
            module,
            'MAKE_UNC_BOLT',
            {'size': '1/4-20', 'length': 1.0},
            source=source
        )
        assert exec_result.success
        assert exec_result.geometry is not None

    def test_dsl_unified_nut(self):
        """Test unified_hex_nut DSL builtin."""
        from yapcad.dsl import tokenize, parse, check, Interpreter

        source = '''
module test_fastener

command MAKE_UNC_NUT(size: string) -> solid:
    nut: solid = unified_hex_nut(size)
    emit nut
'''
        tokens = tokenize(source)
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(
            module,
            'MAKE_UNC_NUT',
            {'size': '1/2-13'},
            source=source
        )
        assert exec_result.success
        assert exec_result.geometry is not None
