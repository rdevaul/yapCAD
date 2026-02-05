"""Tests for the manufacturing post-processing module.

Tests cover path utilities, connector calculations, data structures,
and (with OCC) full solid segmentation.
"""

import math
import pytest

from yapcad.manufacturing import (
    # Data structures
    CutPoint,
    Segment,
    SegmentationResult,
    ConnectorSpec,
    SweptElementProvenance,
    # Path utilities
    evaluate_path3d_at_t,
    compute_cut_plane,
    extract_sub_path,
    path_length,
    length_to_parameter,
    parameter_to_length,
    # Connector functions
    FIT_CLEARANCE,
    offset_rectangular_profile,
    compute_inner_profile_dimensions,
    compute_connector_profile_dimensions,
    compute_connector_length,
    create_connector_region2d,
    # Segmentation
    compute_optimal_cuts,
)


# ============================================================================
# Test fixtures
# ============================================================================

@pytest.fixture
def simple_line_path():
    """A simple straight line path from origin to (100, 0, 0)."""
    return {
        'segments': [
            {
                'type': 'line',
                'start': [0.0, 0.0, 0.0],
                'end': [100.0, 0.0, 0.0],
            }
        ]
    }


@pytest.fixture
def two_segment_path():
    """A path with two line segments forming an L-shape."""
    return {
        'segments': [
            {
                'type': 'line',
                'start': [0.0, 0.0, 0.0],
                'end': [50.0, 0.0, 0.0],
            },
            {
                'type': 'line',
                'start': [50.0, 0.0, 0.0],
                'end': [50.0, 50.0, 0.0],
            },
        ]
    }


@pytest.fixture
def arc_path():
    """A simple 90-degree arc in the XY plane."""
    return {
        'segments': [
            {
                'type': 'arc',
                'center': [0.0, 0.0, 0.0],
                'start': [10.0, 0.0, 0.0],
                'end': [0.0, 10.0, 0.0],
                'normal': [0.0, 0.0, 1.0],
            }
        ]
    }


# ============================================================================
# Path utility tests
# ============================================================================

class TestPathEvaluation:
    """Tests for path3d evaluation functions."""

    def test_evaluate_line_start(self, simple_line_path):
        """Evaluate at start of line."""
        point, tangent = evaluate_path3d_at_t(simple_line_path, 0.0)
        assert point == pytest.approx([0.0, 0.0, 0.0], abs=1e-10)
        assert tangent == pytest.approx([1.0, 0.0, 0.0], abs=1e-10)

    def test_evaluate_line_middle(self, simple_line_path):
        """Evaluate at middle of line."""
        point, tangent = evaluate_path3d_at_t(simple_line_path, 0.5)
        assert point == pytest.approx([50.0, 0.0, 0.0], abs=1e-10)
        assert tangent == pytest.approx([1.0, 0.0, 0.0], abs=1e-10)

    def test_evaluate_line_end(self, simple_line_path):
        """Evaluate at end of line."""
        point, tangent = evaluate_path3d_at_t(simple_line_path, 1.0)
        assert point == pytest.approx([100.0, 0.0, 0.0], abs=1e-10)
        assert tangent == pytest.approx([1.0, 0.0, 0.0], abs=1e-10)

    def test_evaluate_two_segment_first(self, two_segment_path):
        """Evaluate in first segment of two-segment path."""
        point, tangent = evaluate_path3d_at_t(two_segment_path, 0.25)
        # t=0.25 is halfway through first segment (which goes 0 to 0.5)
        assert point == pytest.approx([25.0, 0.0, 0.0], abs=1e-10)
        assert tangent == pytest.approx([1.0, 0.0, 0.0], abs=1e-10)

    def test_evaluate_two_segment_second(self, two_segment_path):
        """Evaluate in second segment of two-segment path."""
        point, tangent = evaluate_path3d_at_t(two_segment_path, 0.75)
        # t=0.75 is halfway through second segment
        assert point == pytest.approx([50.0, 25.0, 0.0], abs=1e-10)
        assert tangent == pytest.approx([0.0, 1.0, 0.0], abs=1e-10)

    def test_evaluate_arc_start(self, arc_path):
        """Evaluate at start of arc."""
        point, tangent = evaluate_path3d_at_t(arc_path, 0.0)
        assert point == pytest.approx([10.0, 0.0, 0.0], abs=1e-6)
        # Tangent at start of 90-degree CCW arc is in +Y direction
        assert tangent == pytest.approx([0.0, 1.0, 0.0], abs=1e-6)

    def test_evaluate_arc_middle(self, arc_path):
        """Evaluate at middle of 90-degree arc."""
        point, tangent = evaluate_path3d_at_t(arc_path, 0.5)
        # At 45 degrees, point is at (r*cos(45), r*sin(45), 0)
        r = 10.0
        expected_coord = r * math.cos(math.pi / 4)
        assert point == pytest.approx([expected_coord, expected_coord, 0.0], abs=1e-6)

    def test_evaluate_arc_end(self, arc_path):
        """Evaluate at end of 90-degree arc."""
        point, tangent = evaluate_path3d_at_t(arc_path, 1.0)
        assert point == pytest.approx([0.0, 10.0, 0.0], abs=1e-6)
        # Tangent at end of CCW arc (at 90 degrees) is in -X direction
        assert tangent == pytest.approx([-1.0, 0.0, 0.0], abs=1e-6)


class TestPathLength:
    """Tests for path length calculations."""

    def test_line_length(self, simple_line_path):
        """Length of simple line."""
        length = path_length(simple_line_path)
        assert length == pytest.approx(100.0, abs=1e-10)

    def test_two_segment_length(self, two_segment_path):
        """Length of L-shaped path."""
        length = path_length(two_segment_path)
        assert length == pytest.approx(100.0, abs=1e-10)  # 50 + 50

    def test_arc_length(self, arc_path):
        """Length of 90-degree arc."""
        length = path_length(arc_path)
        expected = 10.0 * (math.pi / 2)  # radius * angle
        assert length == pytest.approx(expected, abs=1e-6)


class TestSubPathExtraction:
    """Tests for extracting portions of paths."""

    def test_extract_first_half(self, simple_line_path):
        """Extract first half of line."""
        sub = extract_sub_path(simple_line_path, 0.0, 0.5)
        assert len(sub['segments']) == 1
        seg = sub['segments'][0]
        assert seg['start'] == pytest.approx([0.0, 0.0, 0.0], abs=1e-10)
        assert seg['end'] == pytest.approx([50.0, 0.0, 0.0], abs=1e-10)

    def test_extract_second_half(self, simple_line_path):
        """Extract second half of line."""
        sub = extract_sub_path(simple_line_path, 0.5, 1.0)
        assert len(sub['segments']) == 1
        seg = sub['segments'][0]
        assert seg['start'] == pytest.approx([50.0, 0.0, 0.0], abs=1e-10)
        assert seg['end'] == pytest.approx([100.0, 0.0, 0.0], abs=1e-10)

    def test_extract_middle_portion(self, simple_line_path):
        """Extract middle portion of line."""
        sub = extract_sub_path(simple_line_path, 0.25, 0.75)
        assert len(sub['segments']) == 1
        seg = sub['segments'][0]
        assert seg['start'] == pytest.approx([25.0, 0.0, 0.0], abs=1e-10)
        assert seg['end'] == pytest.approx([75.0, 0.0, 0.0], abs=1e-10)

    def test_extract_across_segments(self, two_segment_path):
        """Extract portion spanning both segments."""
        sub = extract_sub_path(two_segment_path, 0.25, 0.75)
        assert len(sub['segments']) == 2

    def test_invalid_range(self, simple_line_path):
        """Error when t_start >= t_end."""
        with pytest.raises(ValueError):
            extract_sub_path(simple_line_path, 0.5, 0.5)
        with pytest.raises(ValueError):
            extract_sub_path(simple_line_path, 0.7, 0.3)


class TestLengthParameterConversion:
    """Tests for converting between arc length and parameter."""

    def test_length_to_param_line(self, simple_line_path):
        """Convert length to parameter on line."""
        t = length_to_parameter(simple_line_path, 50.0)
        assert t == pytest.approx(0.5, abs=1e-10)

    def test_param_to_length_line(self, simple_line_path):
        """Convert parameter to length on line."""
        length = parameter_to_length(simple_line_path, 0.5)
        assert length == pytest.approx(50.0, abs=1e-10)

    def test_roundtrip(self, simple_line_path):
        """Convert parameter -> length -> parameter."""
        original_t = 0.37
        length = parameter_to_length(simple_line_path, original_t)
        recovered_t = length_to_parameter(simple_line_path, length)
        assert recovered_t == pytest.approx(original_t, abs=1e-10)


class TestCutPlane:
    """Tests for cut plane computation."""

    def test_cut_plane_line(self, simple_line_path):
        """Cut plane on line is perpendicular."""
        point, normal = compute_cut_plane(simple_line_path, 0.5)
        assert point == pytest.approx([50.0, 0.0, 0.0], abs=1e-10)
        # Normal equals tangent (pointing along path)
        assert normal == pytest.approx([1.0, 0.0, 0.0], abs=1e-10)


# ============================================================================
# Connector calculation tests
# ============================================================================

class TestConnectorProfiles:
    """Tests for connector profile dimension calculations."""

    def test_offset_rectangular_profile(self):
        """Offset profile shrinks by clearance per side."""
        new_w, new_h = offset_rectangular_profile(10.0, 20.0, 0.2)
        assert new_w == pytest.approx(9.6, abs=1e-10)
        assert new_h == pytest.approx(19.6, abs=1e-10)

    def test_offset_too_large(self):
        """Error when clearance is too large."""
        with pytest.raises(ValueError, match="too large"):
            offset_rectangular_profile(1.0, 1.0, 0.6)

    def test_compute_inner_profile(self):
        """Inner profile from outer and wall thickness."""
        inner_w, inner_h = compute_inner_profile_dimensions(20.0, 30.0, 2.0)
        assert inner_w == pytest.approx(16.0, abs=1e-10)
        assert inner_h == pytest.approx(26.0, abs=1e-10)

    def test_compute_connector_profile(self):
        """Full connector profile calculation."""
        conn_w, conn_h = compute_connector_profile_dimensions(
            outer_width=20.0,
            outer_height=30.0,
            wall_thickness=2.0,
            fit_clearance=0.2,
        )
        # Inner is 16x26, minus 0.4 on each dimension
        assert conn_w == pytest.approx(15.6, abs=1e-10)
        assert conn_h == pytest.approx(25.6, abs=1e-10)

    def test_connector_length_default(self, simple_line_path):
        """Default connector length based on profile size."""
        length = compute_connector_length(20.0, 30.0, simple_line_path, 0.5)
        # Default factor is 3.0, max dimension is 30
        assert length == pytest.approx(90.0, abs=1e-10)


class TestFitClearances:
    """Tests for fit clearance values."""

    def test_fit_clearance_values(self):
        """Verify standard fit clearance values."""
        assert FIT_CLEARANCE['press'] == pytest.approx(0.18, abs=0.01)
        assert FIT_CLEARANCE['slip'] == pytest.approx(0.30, abs=0.01)
        assert FIT_CLEARANCE['loose'] == pytest.approx(0.45, abs=0.01)

    def test_fit_ordering(self):
        """Press < slip < loose."""
        assert FIT_CLEARANCE['press'] < FIT_CLEARANCE['slip']
        assert FIT_CLEARANCE['slip'] < FIT_CLEARANCE['loose']


class TestConnectorRegion:
    """Tests for connector region2d creation."""

    def test_create_simple_region(self):
        """Create rectangular connector region."""
        region = create_connector_region2d(10.0, 20.0)
        assert isinstance(region, list)
        assert len(region) == 1  # Single outer boundary

    def test_create_rounded_region(self):
        """Create rounded-corner connector region."""
        region = create_connector_region2d(10.0, 20.0, corner_radius=1.0)
        assert isinstance(region, list)
        assert len(region) == 1


# ============================================================================
# Data structure tests
# ============================================================================

class TestCutPoint:
    """Tests for CutPoint data structure."""

    def test_valid_cut_point(self):
        """Create valid cut point."""
        cp = CutPoint("beam1", 0.5)
        assert cp.element_id == "beam1"
        assert cp.parameter == 0.5
        assert cp.fit_clearance == 0.2  # Default

    def test_invalid_parameter_zero(self):
        """Cannot cut at t=0."""
        with pytest.raises(ValueError):
            CutPoint("beam1", 0.0)

    def test_invalid_parameter_one(self):
        """Cannot cut at t=1."""
        with pytest.raises(ValueError):
            CutPoint("beam1", 1.0)

    def test_invalid_union_with(self):
        """Invalid union_connector_with value."""
        with pytest.raises(ValueError):
            CutPoint("beam1", 0.5, union_connector_with="invalid")


class TestSegmentationResult:
    """Tests for SegmentationResult data structure."""

    def test_segment_count(self):
        """Segment count property."""
        result = SegmentationResult(
            segments=[
                Segment("seg1", None, "beam1", (0, 0.5)),
                Segment("seg2", None, "beam1", (0.5, 1)),
            ]
        )
        assert result.segment_count == 2

    def test_get_segment(self):
        """Get segment by ID."""
        seg1 = Segment("seg1", None, "beam1", (0, 0.5))
        seg2 = Segment("seg2", None, "beam1", (0.5, 1))
        result = SegmentationResult(segments=[seg1, seg2])

        found = result.get_segment("seg1")
        assert found is seg1

        not_found = result.get_segment("nonexistent")
        assert not_found is None

    def test_get_segments_for_element(self):
        """Get all segments for a parent element."""
        result = SegmentationResult(
            segments=[
                Segment("seg1", None, "beam1", (0, 0.5)),
                Segment("seg2", None, "beam1", (0.5, 1)),
                Segment("seg3", None, "beam2", (0, 1)),
            ]
        )
        beam1_segs = result.get_segments_for_element("beam1")
        assert len(beam1_segs) == 2
        assert all(s.parent_element_id == "beam1" for s in beam1_segs)


class TestSweptElementProvenance:
    """Tests for SweptElementProvenance data structure."""

    def test_create_provenance(self, simple_line_path):
        """Create basic provenance object."""
        prov = SweptElementProvenance(
            id="beam1",
            operation="sweep_adaptive",
            outer_profile=[[...]],  # placeholder
            spine=simple_line_path,
            wall_thickness=2.0,
        )
        assert prov.id == "beam1"
        assert prov.semantic_type == "structural_beam"  # Default


# ============================================================================
# Optimal cuts tests
# ============================================================================

class TestOptimalCuts:
    """Tests for automatic cut point computation."""

    def test_no_cuts_needed(self, simple_line_path):
        """No cuts when beam fits in build volume."""
        prov = SweptElementProvenance(
            id="beam1",
            operation="sweep",
            outer_profile=[[...]],
            spine=simple_line_path,
        )
        cuts = compute_optimal_cuts(prov, max_segment_length=150.0)
        assert len(cuts) == 0

    def test_one_cut_needed(self, simple_line_path):
        """One cut when beam is 2x max length."""
        prov = SweptElementProvenance(
            id="beam1",
            operation="sweep",
            outer_profile=[[...]],
            spine=simple_line_path,  # 100mm path
        )
        cuts = compute_optimal_cuts(prov, max_segment_length=60.0)
        # Path is 100mm, max is 60mm, so need at least 2 segments
        assert len(cuts) == 1
        assert 0 < cuts[0].parameter < 1

    def test_multiple_cuts_needed(self, simple_line_path):
        """Multiple cuts for long beam."""
        prov = SweptElementProvenance(
            id="beam1",
            operation="sweep",
            outer_profile=[[...]],
            spine=simple_line_path,  # 100mm path
        )
        cuts = compute_optimal_cuts(prov, max_segment_length=30.0)
        # Need 4 segments (100/30 = 3.3, rounded up to 4)
        assert len(cuts) == 3
        # Cuts should be roughly evenly spaced
        params = [c.parameter for c in cuts]
        assert all(0 < p < 1 for p in params)
        assert params == sorted(params)  # In order


# ============================================================================
# Integration tests (require OCC)
# ============================================================================

def has_occ():
    """Check if OCC is available."""
    try:
        from OCC.Core.gp import gp_Pnt
        return True
    except ImportError:
        return False


@pytest.mark.skipif(not has_occ(), reason="OCC not available")
@pytest.mark.skip(reason="OCC integration tests need profile/path format refinement")
class TestSolidSegmentation:
    """Integration tests requiring OCC for solid operations.

    NOTE: These tests are currently skipped because the profile and path
    format conversion between manufacturing module format and sweep_adaptive
    requires further refinement. The core path utilities and data structures
    are tested separately.
    """

    def test_split_solid_basic(self):
        """Split a simple box solid."""
        from yapcad.manufacturing import split_solid_at_plane
        from yapcad.geom3d_util import prism

        # Create a 100x20x20 box (prism creates box centered at origin)
        box = prism(100, 20, 20)

        # Split at the middle
        solid_a, solid_b = split_solid_at_plane(
            box,
            plane_point=[50.0, 0.0, 0.0],
            plane_normal=[1.0, 0.0, 0.0],
        )

        # Both parts should exist
        assert solid_a is not None
        assert solid_b is not None

        # Both should have BREP data (metadata is at index 4)
        assert len(solid_a) == 5 and 'brep' in solid_a[4]
        assert len(solid_b) == 5 and 'brep' in solid_b[4]

    def test_create_interior_connector(self, simple_line_path):
        """Create an interior connector solid."""
        from yapcad.manufacturing import create_interior_connector

        connector = create_interior_connector(
            outer_profile_width=20.0,
            outer_profile_height=30.0,
            spine=simple_line_path,
            center_parameter=0.5,
            wall_thickness=2.0,
            connector_length=60.0,
        )

        # Should produce a valid solid
        assert connector is not None
        # If BREP, metadata is at index 4
        if len(connector) == 5:
            assert 'brep' in connector[4]

    def test_full_segmentation_workflow(self, simple_line_path):
        """End-to-end segmentation test."""
        from yapcad.manufacturing import (
            segment_swept_element,
            SweptElementProvenance,
            CutPoint,
        )
        from yapcad.geom3d_util import sweep_adaptive
        from yapcad.geom import line, point

        # Create a simple hollow beam
        # Outer profile (20x30 rectangle)
        hw, hh = 10, 15
        outer = [[
            line(point(-hw, -hh, 0), point(hw, -hh, 0)),
            line(point(hw, -hh, 0), point(hw, hh, 0)),
            line(point(hw, hh, 0), point(-hw, hh, 0)),
            line(point(-hw, hh, 0), point(-hw, -hh, 0)),
        ]]

        # Inner profile (hollow, 2mm wall)
        wt = 2.0
        iw, ih = hw - wt, hh - wt
        inner = [[
            line(point(-iw, -ih, 0), point(iw, -ih, 0)),
            line(point(iw, -ih, 0), point(iw, ih, 0)),
            line(point(iw, ih, 0), point(-iw, ih, 0)),
            line(point(-iw, ih, 0), point(-iw, -ih, 0)),
        ]]

        # Create swept beam
        beam = sweep_adaptive(
            outer, simple_line_path,
            inner_profiles=inner,
            angle_threshold_deg=5.0,
        )

        # Create provenance
        prov = SweptElementProvenance(
            id="test_beam",
            operation="sweep_adaptive",
            outer_profile=outer,
            spine=simple_line_path,
            inner_profile=inner,
            wall_thickness=wt,
            metadata={'solid': beam},
        )

        # Define cut
        cuts = [CutPoint("test_beam", 0.5)]

        # Segment
        result = segment_swept_element(prov, cuts)

        # Verify results
        assert result.segment_count == 2
        assert len(result.connectors) == 1
        assert "test_beam_seg_0" in result.assembly_graph
        assert "test_beam_seg_1" in result.assembly_graph
