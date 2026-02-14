"""Tests for BREP edge selection helper functions."""

import math
import pytest

from yapcad.brep import (
    BrepSolid,
    BrepEdge,
    occ_available,
    fillet_edges,
    chamfer_edges,
)

try:
    from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder
except ImportError:
    BRepPrimAPI_MakeBox = None
    BRepPrimAPI_MakeCylinder = None

# Skip all tests if OCC is not available
pytestmark = pytest.mark.skipif(not occ_available(), reason="pythonocc-core not available")


@pytest.fixture
def box_brep():
    """Create a 10x10x10 box BREP solid."""
    if BRepPrimAPI_MakeBox is None:
        pytest.skip("pythonocc-core not available")
    shape = BRepPrimAPI_MakeBox(10.0, 10.0, 10.0).Shape()
    return BrepSolid(shape)


@pytest.fixture
def tall_box_brep():
    """Create a 5x5x20 tall box BREP solid."""
    if BRepPrimAPI_MakeBox is None:
        pytest.skip("pythonocc-core not available")
    shape = BRepPrimAPI_MakeBox(5.0, 5.0, 20.0).Shape()
    return BrepSolid(shape)


@pytest.fixture
def cylinder_brep():
    """Create a cylinder BREP solid (radius 5, height 10)."""
    if BRepPrimAPI_MakeCylinder is None:
        pytest.skip("pythonocc-core not available")
    shape = BRepPrimAPI_MakeCylinder(5.0, 10.0).Shape()
    return BrepSolid(shape)


class TestGetAllEdges:
    """Tests for get_all_edges function."""

    def test_box_has_12_edges(self, box_brep):
        from yapcad.brep_edge_select import get_all_edges
        edges = get_all_edges(box_brep)
        assert len(edges) == 12
        assert all(isinstance(e, BrepEdge) for e in edges)

    def test_cylinder_has_edges(self, cylinder_brep):
        from yapcad.brep_edge_select import get_all_edges
        edges = get_all_edges(cylinder_brep)
        # Cylinder has 2 circular edges (top and bottom) + possibly seam edge
        assert len(edges) >= 2


class TestSelectVerticalEdges:
    """Tests for select_vertical_edges function."""

    def test_box_has_4_vertical_edges(self, box_brep):
        from yapcad.brep_edge_select import select_vertical_edges
        vertical_edges = select_vertical_edges(box_brep)
        assert len(vertical_edges) == 4

    def test_tall_box_vertical_edges_are_longer(self, tall_box_brep):
        from yapcad.brep_edge_select import select_vertical_edges, edge_info
        vertical_edges = select_vertical_edges(tall_box_brep)
        assert len(vertical_edges) == 4
        # All vertical edges should be 20 units long
        for edge in vertical_edges:
            info = edge_info(edge)
            assert math.isclose(info['length'], 20.0, rel_tol=0.01)


class TestSelectHorizontalEdges:
    """Tests for select_horizontal_edges function."""

    def test_box_has_8_horizontal_edges(self, box_brep):
        from yapcad.brep_edge_select import select_horizontal_edges
        horizontal_edges = select_horizontal_edges(box_brep)
        assert len(horizontal_edges) == 8

    def test_horizontal_edges_are_perpendicular_to_z(self, box_brep):
        from yapcad.brep_edge_select import select_horizontal_edges, edge_info
        horizontal_edges = select_horizontal_edges(box_brep)
        for edge in horizontal_edges:
            info = edge_info(edge)
            assert info['is_horizontal'] is True


class TestSelectEdgesByDirection:
    """Tests for select_edges_by_direction function."""

    def test_x_direction_edges(self, box_brep):
        from yapcad.brep_edge_select import select_edges_by_direction
        x_edges = select_edges_by_direction(box_brep, (1, 0, 0))
        # A box has 4 edges parallel to X axis
        assert len(x_edges) == 4

    def test_y_direction_edges(self, box_brep):
        from yapcad.brep_edge_select import select_edges_by_direction
        y_edges = select_edges_by_direction(box_brep, (0, 1, 0))
        # A box has 4 edges parallel to Y axis
        assert len(y_edges) == 4

    def test_z_direction_edges(self, box_brep):
        from yapcad.brep_edge_select import select_edges_by_direction
        z_edges = select_edges_by_direction(box_brep, (0, 0, 1))
        # A box has 4 edges parallel to Z axis (vertical edges)
        assert len(z_edges) == 4

    def test_no_diagonal_edges_in_box(self, box_brep):
        from yapcad.brep_edge_select import select_edges_by_direction
        diagonal_edges = select_edges_by_direction(box_brep, (1, 1, 1))
        assert len(diagonal_edges) == 0


class TestSelectEdgesByLength:
    """Tests for select_edges_by_length function."""

    def test_all_edges_of_cube_same_length(self, box_brep):
        from yapcad.brep_edge_select import select_edges_by_length
        # All edges of a 10x10x10 cube are 10 units long
        edges = select_edges_by_length(box_brep, min_length=9.9, max_length=10.1)
        assert len(edges) == 12

    def test_filter_short_edges(self, tall_box_brep):
        from yapcad.brep_edge_select import select_edges_by_length
        # tall_box is 5x5x20
        # Short edges are 5 units, long edges are 20 units
        short_edges = select_edges_by_length(tall_box_brep, max_length=6.0)
        assert len(short_edges) == 8  # 4 top + 4 bottom

    def test_filter_long_edges(self, tall_box_brep):
        from yapcad.brep_edge_select import select_edges_by_length
        long_edges = select_edges_by_length(tall_box_brep, min_length=15.0)
        assert len(long_edges) == 4  # Vertical edges


class TestSelectEdgesAtZ:
    """Tests for select_edges_at_z function."""

    def test_bottom_edges_at_z0(self, box_brep):
        from yapcad.brep_edge_select import select_edges_at_z
        bottom_edges = select_edges_at_z(box_brep, 0.0)
        assert len(bottom_edges) == 4

    def test_top_edges_at_z10(self, box_brep):
        from yapcad.brep_edge_select import select_edges_at_z
        top_edges = select_edges_at_z(box_brep, 10.0)
        assert len(top_edges) == 4

    def test_no_edges_at_z5(self, box_brep):
        from yapcad.brep_edge_select import select_edges_at_z
        # No horizontal edges at Z=5 in a box (only vertical edges pass through)
        mid_edges = select_edges_at_z(box_brep, 5.0)
        assert len(mid_edges) == 0


class TestSelectTopBottomEdges:
    """Tests for select_top_edges and select_bottom_edges."""

    def test_select_top_edges(self, box_brep):
        from yapcad.brep_edge_select import select_top_edges
        top_edges = select_top_edges(box_brep)
        assert len(top_edges) == 4

    def test_select_bottom_edges(self, box_brep):
        from yapcad.brep_edge_select import select_bottom_edges
        bottom_edges = select_bottom_edges(box_brep)
        assert len(bottom_edges) == 4


class TestSelectEdgesCrossingZ:
    """Tests for select_edges_crossing_z function."""

    def test_vertical_edges_cross_z5(self, box_brep):
        from yapcad.brep_edge_select import select_edges_crossing_z
        crossing = select_edges_crossing_z(box_brep, 5.0)
        assert len(crossing) == 4  # The 4 vertical edges


class TestEdgeInfo:
    """Tests for edge_info function."""

    def test_edge_info_returns_dict(self, box_brep):
        from yapcad.brep_edge_select import get_all_edges, edge_info
        edges = get_all_edges(box_brep)
        info = edge_info(edges[0])

        assert 'length' in info
        assert 'endpoints' in info
        assert 'midpoint' in info
        assert 'direction' in info
        assert 'is_linear' in info
        assert 'is_vertical' in info
        assert 'is_horizontal' in info

    def test_linear_edge_has_direction(self, box_brep):
        from yapcad.brep_edge_select import get_all_edges, edge_info
        edges = get_all_edges(box_brep)
        for edge in edges:
            info = edge_info(edge)
            assert info['is_linear'] is True
            assert info['direction'] is not None


class TestFilterFunctions:
    """Tests for filter_curved_edges and filter_linear_edges."""

    def test_filter_linear_preserves_box_edges(self, box_brep):
        from yapcad.brep_edge_select import get_all_edges, filter_linear_edges
        all_edges = get_all_edges(box_brep)
        linear = filter_linear_edges(all_edges)
        assert len(linear) == 12  # All box edges are linear

    def test_filter_curved_on_box_returns_empty(self, box_brep):
        from yapcad.brep_edge_select import get_all_edges, filter_curved_edges
        all_edges = get_all_edges(box_brep)
        curved = filter_curved_edges(all_edges)
        assert len(curved) == 0  # Box has no curved edges

    def test_cylinder_has_curved_edges(self, cylinder_brep):
        from yapcad.brep_edge_select import get_all_edges, filter_curved_edges
        all_edges = get_all_edges(cylinder_brep)
        curved = filter_curved_edges(all_edges)
        # Cylinder has 2 circular edges (top and bottom)
        assert len(curved) >= 2


class TestSetOperations:
    """Tests for union_edges, intersect_edges, subtract_edges."""

    def test_union_removes_duplicates(self, box_brep):
        from yapcad.brep_edge_select import (
            select_vertical_edges, select_horizontal_edges, union_edges
        )
        vertical = select_vertical_edges(box_brep)
        horizontal = select_horizontal_edges(box_brep)
        combined = union_edges(vertical, horizontal)
        assert len(combined) == 12  # 4 + 8 = 12 unique edges

    def test_intersect_finds_common(self, box_brep):
        from yapcad.brep_edge_select import (
            select_vertical_edges, select_edges_by_length, intersect_edges
        )
        vertical = select_vertical_edges(box_brep)
        # All edges of the cube are 10 units
        length_10 = select_edges_by_length(box_brep, min_length=9.9, max_length=10.1)
        common = intersect_edges(vertical, length_10)
        assert len(common) == 4  # Vertical edges that are 10 units long

    def test_subtract_removes_edges(self, box_brep):
        from yapcad.brep_edge_select import (
            get_all_edges, select_vertical_edges, subtract_edges
        )
        all_edges = get_all_edges(box_brep)
        vertical = select_vertical_edges(box_brep)
        non_vertical = subtract_edges(all_edges, vertical)
        assert len(non_vertical) == 8  # 12 - 4 = 8


class TestFilletSelectiveEdges:
    """Test applying fillet only to selected edges."""

    def test_fillet_only_vertical_edges(self, box_brep):
        """Example: fillet only vertical edges of a box."""
        from yapcad.brep_edge_select import select_vertical_edges

        vertical_edges = select_vertical_edges(box_brep)
        assert len(vertical_edges) == 4

        # Apply fillet only to vertical edges
        filleted = fillet_edges(box_brep, vertical_edges, 0.5)

        assert filleted is not None
        assert isinstance(filleted, BrepSolid)

        # Verify tessellation works
        surface = filleted.tessellate()
        assert surface[0] == 'surface'
        assert len(surface[1]) > 0

    def test_fillet_only_top_edges(self, box_brep):
        """Fillet only the top edges of a box."""
        from yapcad.brep_edge_select import select_top_edges

        top_edges = select_top_edges(box_brep)
        assert len(top_edges) == 4

        filleted = fillet_edges(box_brep, top_edges, 0.5)

        assert filleted is not None
        surface = filleted.tessellate()
        assert surface[0] == 'surface'

    def test_chamfer_only_bottom_edges(self, box_brep):
        """Chamfer only the bottom edges of a box."""
        from yapcad.brep_edge_select import select_bottom_edges

        bottom_edges = select_bottom_edges(box_brep)
        assert len(bottom_edges) == 4

        chamfered = chamfer_edges(box_brep, bottom_edges, 0.3)

        assert chamfered is not None
        surface = chamfered.tessellate()
        assert surface[0] == 'surface'


class TestPocketVerticalEdgeExample:
    """Test case demonstrating filleting vertical edges of a pocket."""

    def test_fillet_pocket_vertical_edges(self):
        """
        Example use case: Create a plate with a pocket, then fillet
        only the vertical edges of the pocket.
        """
        if BRepPrimAPI_MakeBox is None:
            pytest.skip("pythonocc-core not available")

        from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut
        from yapcad.brep_edge_select import (
            select_vertical_edges,
            select_edges_in_z_range,
            intersect_edges,
        )

        # Create base plate: 50x50x10
        plate_shape = BRepPrimAPI_MakeBox(50.0, 50.0, 10.0).Shape()

        # Create pocket tool: 30x30x6, positioned at center
        # (starting at x=10, y=10, z=4 to make it a 6mm deep pocket from top)
        from OCC.Core.gp import gp_Vec
        from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_Transform
        from OCC.Core.gp import gp_Trsf

        pocket_shape = BRepPrimAPI_MakeBox(30.0, 30.0, 6.0).Shape()

        # Translate pocket to center of plate, sitting on top surface
        trsf = gp_Trsf()
        trsf.SetTranslation(gp_Vec(10.0, 10.0, 4.0))
        transform = BRepBuilderAPI_Transform(pocket_shape, trsf, True)
        pocket_translated = transform.Shape()

        # Cut pocket from plate
        cutter = BRepAlgoAPI_Cut(plate_shape, pocket_translated)
        result_shape = cutter.Shape()
        pocket_brep = BrepSolid(result_shape)

        # Get vertical edges
        vertical_edges = select_vertical_edges(pocket_brep)
        assert len(vertical_edges) > 4  # Should have plate corners + pocket corners

        # Filter to only edges in the pocket region (Z between 4 and 10)
        pocket_region_edges = select_edges_in_z_range(pocket_brep, 4.0, 10.0)

        # Intersect to get only vertical edges in the pocket region
        pocket_vertical = intersect_edges(vertical_edges, pocket_region_edges)

        # These should be the 4 vertical edges of the pocket
        # (pocket corners go from Z=4 to Z=10)
        # Note: This may include some plate edges depending on geometry
        assert len(pocket_vertical) >= 4

        # Apply fillet to pocket vertical edges
        filleted = fillet_edges(pocket_brep, pocket_vertical, 1.0)
        assert filleted is not None

        # Verify result
        surface = filleted.tessellate()
        assert surface[0] == 'surface'
        assert len(surface[1]) > 0
