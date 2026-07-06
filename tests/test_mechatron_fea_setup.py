"""Tests for the mechatron-canonical FEA setup helper.

These tests validate the FEA-setup CODE against a small synthetic
assembly graph fixture (``tests/fixtures/sample_assembly_graph.json``),
NOT against any real/internal design. The fixture uses generic invented
parts (base-plate, riser-column, top-bracket, ...) and hand-picked
bolt-pattern / load-case values so every assertion is deterministic.
"""
import json
import math
from pathlib import Path

import pytest

from yapcad.package.analysis.mechatron_fea_setup import (
    Bolt,
    FeaSetup,
    LoadAttachResolved,
    generate_bolts_for_interface,
    get_interface,
    get_load_case,
    list_interfaces_with_bolt_patterns,
    list_load_cases,
    load_graph,
    prepare,
    resolve_load_attach,
    stiffness_for_bolt_spec,
    vehicle_to_mesh_direction,
    vehicle_to_mesh_position,
)


# ---------------------------------------------------------------------------
# Stiffness library
# ---------------------------------------------------------------------------

class TestStiffness:
    def test_quarter_20_spec(self):
        k_a, k_s = stiffness_for_bolt_spec("1-4-20-x-1in-button-head-shcs")
        # 1/4-20 = 6.35mm shaft, A = π·6.35²/4 = 31.7 mm² = 3.17e-5 m²
        # k_axial = 200e9 × 3.17e-5 / 0.012 = 5.28e8 N/m = 528 kN/mm
        assert 5e8 < k_a < 7e8, f"k_axial out of range: {k_a:.2e}"
        assert k_s < k_a  # shear is always less than axial for steel

    def test_m5_spec(self):
        k_a, k_s = stiffness_for_bolt_spec("m5-x-16-shcs")
        # 5mm shaft, A = 1.96e-5 m² → k_axial ≈ 3.3e8 N/m
        assert 3e8 < k_a < 4e8

    def test_quarter_20_stiffer_than_m5(self):
        """The 1/4-20 spec should be ~60% stiffer than M5 (A ratio = 1.61)."""
        k_q20, _ = stiffness_for_bolt_spec("1-4-20-x-1in-button-head-shcs")
        k_m5, _ = stiffness_for_bolt_spec("m5-x-16-shcs")
        assert k_q20 > k_m5
        assert 1.5 < (k_q20 / k_m5) < 1.7

    def test_placeholder_returns_zero(self):
        k_a, k_s = stiffness_for_bolt_spec("PLACEHOLDER-engine-chamber-nozzle")
        assert k_a == 0.0
        assert k_s == 0.0

    def test_unknown_spec_falls_back_to_m5(self, capsys):
        k_a, k_s = stiffness_for_bolt_spec("unknown-bolt-spec")
        captured = capsys.readouterr()
        assert "unknown bolt_spec" in captured.err
        # Should be M5-ish
        assert 3e8 < k_a < 4e8


# ---------------------------------------------------------------------------
# Coordinate frame transforms
# ---------------------------------------------------------------------------

class TestCoordinateTransforms:
    def test_vehicle_to_mesh_position(self):
        # vehicle (X=axial=1664.6, Y=0, Z=side=145) → mesh (X=0, Y=145, Z=1664.6)
        m = vehicle_to_mesh_position([1664.6, 0, 145])
        assert m == (0.0, 145.0, 1664.6)

    def test_vehicle_to_mesh_direction(self):
        # vehicle [0,0,1] (radial-out) → mesh [0,1,0] (+Y radial)
        m = vehicle_to_mesh_direction([0, 0, 1])
        assert m == (0.0, 1.0, 0.0)


# ---------------------------------------------------------------------------
# graph.json loading — against the synthetic fixture
# ---------------------------------------------------------------------------

FIXTURE_PATH = str(Path(__file__).parent / "fixtures" / "sample_assembly_graph.json")


@pytest.fixture
def sample_graph():
    """The synthetic sample assembly graph (generic widget parts)."""
    return load_graph(FIXTURE_PATH)


class TestGraphLoading:
    def test_load_sample_graph(self, sample_graph):
        assert "interfaces" in sample_graph
        assert "load_cases" in sample_graph
        assert len(sample_graph["load_cases"]) == 3

    def test_load_graph_raises_on_missing_file(self):
        with pytest.raises(FileNotFoundError):
            load_graph("/no/such/graph.json")

    def test_get_load_case(self, sample_graph):
        lc = get_load_case(sample_graph, "LC-001")
        assert lc["name"] == "Side load on riser bracket (clock 90)"
        assert lc["attach"]["part"] == "top-bracket"

    def test_get_load_case_raises_on_missing(self, sample_graph):
        with pytest.raises(KeyError, match="LoadCase 'LC-999'"):
            get_load_case(sample_graph, "LC-999")

    def test_list_load_cases(self, sample_graph):
        lcs = list_load_cases(sample_graph)
        assert [lc["id"] for lc in lcs] == ["LC-001", "LC-002", "LC-003"]

    def test_list_interfaces_with_bolt_patterns(self, sample_graph):
        ifaces = list_interfaces_with_bolt_patterns(sample_graph)
        # 4 interfaces carry a usable bolt_pattern; the PLACEHOLDER one is
        # skipped and the world↔base one has no bolt_pattern.
        assert len(ifaces) == 4
        ids = [i["id"] for i in ifaces]
        assert "iface_world_base" not in ids  # no bolt_pattern
        assert "iface_placeholder_joint" not in ids  # PLACEHOLDER bolt_spec


# ---------------------------------------------------------------------------
# Bolt enumeration
# ---------------------------------------------------------------------------

class TestBoltEnumeration:
    def test_generate_bolts_for_base_riser(self, sample_graph):
        iface = get_interface(sample_graph, "iface_base_riser")
        bolts = generate_bolts_for_interface(iface, parent_world_origin_mm=(0, 0, 100.0))
        assert len(bolts) == 12  # 12 bolts per the spec
        # First bolt sits at clock=0 (+X on the PCD), R = PCD/2 = 60.
        b0 = bolts[0]
        assert abs(b0.parent_world_mm[0] - 60.0) < 0.01
        assert abs(b0.parent_world_mm[1]) < 0.01  # Y=0
        assert abs(b0.parent_world_mm[2] - 100.0) < 0.01  # world_z
        # Clock spacing is 30° (360/12).
        assert abs(bolts[1].clock_deg - 30.0) < 0.1
        # All bolts share the same spec + a real (non-zero) stiffness.
        for b in bolts:
            assert b.bolt_spec == "1-4-20-x-1in-button-head-shcs"
            assert b.k_axial_N_per_m > 4e8

    def test_clock_90_bolt(self, sample_graph):
        """Find the +Y bolt at clock=90° on the base↔riser interface."""
        iface = get_interface(sample_graph, "iface_base_riser")
        bolts = generate_bolts_for_interface(iface, parent_world_origin_mm=(0, 0, 100.0))
        target = [b for b in bolts if abs(b.clock_deg - 90.0) < 1.0]
        assert len(target) == 1, f"expected one bolt at clock=90°, got {len(target)}"
        b = target[0]
        # World position should be roughly (0, +60, 100) — Y=+60 = PCD/2.
        assert abs(b.parent_world_mm[0]) < 1.0  # X ≈ 0
        assert abs(b.parent_world_mm[1] - 60.0) < 1.0  # Y ≈ +60
        assert abs(b.parent_world_mm[2] - 100.0) < 1.0

    def test_uses_world_z_from_bolt_pattern(self, sample_graph):
        """With no explicit origin override, generator reads world_z_mm."""
        iface = get_interface(sample_graph, "iface_bracket_sensor")
        bolts = generate_bolts_for_interface(iface)  # no override
        assert len(bolts) == 6
        assert all(abs(b.parent_world_mm[2] - 300.0) < 0.01 for b in bolts)


# ---------------------------------------------------------------------------
# Load case resolution
# ---------------------------------------------------------------------------

class TestLoadAttachResolution:
    def test_resolve_lc001(self, sample_graph):
        lc = get_load_case(sample_graph, "LC-001")
        attach = resolve_load_attach(lc)
        assert attach.load_case_id == "LC-001"
        assert attach.part == "top-bracket"
        assert attach.interface_id == "iface_riser_bracket"
        assert attach.clock_deg == 90.0
        assert attach.magnitude_n == 6000.0
        assert attach.coordinate_frame == "vehicle"
        assert attach.direction_unit == (0.0, 0.0, 1.0)

    def test_resolve_whole_assembly_load(self, sample_graph):
        """LC-003 has an empty attach → defaults to the whole assembly."""
        lc = get_load_case(sample_graph, "LC-003")
        attach = resolve_load_attach(lc)
        assert attach.part == "_assembly"
        assert attach.interface_id is None
        assert attach.magnitude_n == 500.0


# ---------------------------------------------------------------------------
# Full prepare() workflow
# ---------------------------------------------------------------------------

class TestPrepare:
    def test_prepare_with_load_case(self, sample_graph):
        setup = prepare(FIXTURE_PATH, load_case_id="LC-001")
        # 4 bolt-pattern interfaces: 12 + 12 + 6 + 8 = 38 bolts.
        assert len(setup.bolts) == 38
        assert setup.load_case_id == "LC-001"
        assert setup.load_attach.part == "top-bracket"

    def test_prepare_without_load_case(self, sample_graph):
        setup = prepare(FIXTURE_PATH)
        assert setup.load_case_id is None
        assert setup.load_attach is None
        assert len(setup.bolts) == 38

    def test_prepare_raises_on_missing_load_case(self, sample_graph):
        with pytest.raises(KeyError):
            prepare(FIXTURE_PATH, load_case_id="LC-bogus")
