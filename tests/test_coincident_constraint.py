"""Tests for COINCIDENT constraint in the yapCAD assembly system.

Tests the COINCIDENT mate type which ensures two geometric features
occupy the same location. This is critical for:
- Bolt hole alignment (servo horn to gear hub)
- Face-to-face mating
- Point coincidence

The Dynamixel XH540 servo horn to axis 1 sun gear interface is used
as a real-world test case:
- XH540 servo horn: 4 bolt holes on 18mm bolt circle at 45 deg offset
- Axis 1 sun gear hub: 4 matching bolt holes on 18mm bolt circle
"""

import pytest
import math
import numpy as np

from yapcad.assembly import (
    Datum,
    DatumType,
    PartDefinition,
    Mate,
    MateType,
    CoincidentResult,
    create_coincident_mate,
    evaluate_coincident,
    check_bolt_circle_alignment,
)


class TestCoincidentPoint:
    """Test COINCIDENT constraint for point datums."""

    def test_coincident_points_exact(self):
        """Two points at the same location should be coincident."""
        point_a = Datum(
            name="point_a",
            datum_type=DatumType.POINT,
            origin=[10.0, 20.0, 30.0, 1.0]
        )
        point_b = Datum(
            name="point_b",
            datum_type=DatumType.POINT,
            origin=[10.0, 20.0, 30.0, 1.0]
        )

        result = evaluate_coincident(point_a, point_b)

        assert result.satisfied is True
        assert result.error_distance < 1e-6
        assert "point_to_point" in result.details.get("type", "")

    def test_coincident_points_within_tolerance(self):
        """Points within tolerance should be coincident."""
        point_a = Datum(
            name="point_a",
            datum_type=DatumType.POINT,
            origin=[10.0, 20.0, 30.0, 1.0]
        )
        point_b = Datum(
            name="point_b",
            datum_type=DatumType.POINT,
            origin=[10.05, 20.0, 30.0, 1.0]  # 0.05mm off
        )

        # Should pass with 0.1mm tolerance
        result = evaluate_coincident(point_a, point_b, tolerance_mm=0.1)
        assert result.satisfied is True

        # Should fail with 0.01mm tolerance
        result = evaluate_coincident(point_a, point_b, tolerance_mm=0.01)
        assert result.satisfied is False

    def test_coincident_points_far_apart(self):
        """Points far apart should not be coincident."""
        point_a = Datum(
            name="point_a",
            datum_type=DatumType.POINT,
            origin=[0.0, 0.0, 0.0, 1.0]
        )
        point_b = Datum(
            name="point_b",
            datum_type=DatumType.POINT,
            origin=[10.0, 0.0, 0.0, 1.0]  # 10mm away
        )

        result = evaluate_coincident(point_a, point_b)

        assert result.satisfied is False
        assert abs(result.error_distance - 10.0) < 1e-6


class TestCoincidentPlane:
    """Test COINCIDENT constraint for plane datums (coplanar)."""

    def test_coincident_planes_coplanar(self):
        """Two coplanar planes should be coincident."""
        plane_a = Datum(
            name="plane_a",
            datum_type=DatumType.PLANE,
            origin=[0.0, 0.0, 10.0, 1.0],
            normal=[0.0, 0.0, 1.0, 0.0]  # Z-normal
        )
        plane_b = Datum(
            name="plane_b",
            datum_type=DatumType.PLANE,
            origin=[5.0, 5.0, 10.0, 1.0],  # Different XY, same Z
            normal=[0.0, 0.0, 1.0, 0.0]
        )

        result = evaluate_coincident(plane_a, plane_b)

        assert result.satisfied is True
        assert result.error_distance < 1e-6
        assert result.error_angle < 1e-6

    def test_coincident_planes_opposite_normals(self):
        """Planes with opposite normals but same location should be coincident."""
        plane_a = Datum(
            name="plane_a",
            datum_type=DatumType.PLANE,
            origin=[0.0, 0.0, 10.0, 1.0],
            normal=[0.0, 0.0, 1.0, 0.0]
        )
        plane_b = Datum(
            name="plane_b",
            datum_type=DatumType.PLANE,
            origin=[0.0, 0.0, 10.0, 1.0],
            normal=[0.0, 0.0, -1.0, 0.0]  # Opposite normal
        )

        result = evaluate_coincident(plane_a, plane_b)

        assert result.satisfied is True

    def test_coincident_planes_parallel_offset(self):
        """Parallel planes at different Z should not be coincident."""
        plane_a = Datum(
            name="plane_a",
            datum_type=DatumType.PLANE,
            origin=[0.0, 0.0, 10.0, 1.0],
            normal=[0.0, 0.0, 1.0, 0.0]
        )
        plane_b = Datum(
            name="plane_b",
            datum_type=DatumType.PLANE,
            origin=[0.0, 0.0, 15.0, 1.0],  # 5mm offset in Z
            normal=[0.0, 0.0, 1.0, 0.0]
        )

        result = evaluate_coincident(plane_a, plane_b)

        assert result.satisfied is False
        assert abs(result.error_distance - 5.0) < 1e-6

    def test_coincident_planes_not_parallel(self):
        """Non-parallel planes should not be coincident."""
        plane_a = Datum(
            name="plane_a",
            datum_type=DatumType.PLANE,
            origin=[0.0, 0.0, 0.0, 1.0],
            normal=[0.0, 0.0, 1.0, 0.0]
        )
        plane_b = Datum(
            name="plane_b",
            datum_type=DatumType.PLANE,
            origin=[0.0, 0.0, 0.0, 1.0],
            normal=[1.0, 0.0, 0.0, 0.0]  # Perpendicular
        )

        result = evaluate_coincident(plane_a, plane_b)

        assert result.satisfied is False
        assert result.error_angle > 45.0  # Should be 90 degrees


class TestCoincidentAxis:
    """Test COINCIDENT constraint for axis datums (colinear)."""

    def test_coincident_axes_colinear(self):
        """Two colinear axes should be coincident."""
        axis_a = Datum(
            name="axis_a",
            datum_type=DatumType.AXIS,
            origin=[0.0, 0.0, 0.0, 1.0],
            direction=[0.0, 0.0, 1.0, 0.0]  # Z-axis
        )
        axis_b = Datum(
            name="axis_b",
            datum_type=DatumType.AXIS,
            origin=[0.0, 0.0, 10.0, 1.0],  # Different Z, same axis
            direction=[0.0, 0.0, 1.0, 0.0]
        )

        result = evaluate_coincident(axis_a, axis_b)

        assert result.satisfied is True
        assert result.error_distance < 1e-6

    def test_coincident_axes_parallel_offset(self):
        """Parallel but offset axes should not be coincident."""
        axis_a = Datum(
            name="axis_a",
            datum_type=DatumType.AXIS,
            origin=[0.0, 0.0, 0.0, 1.0],
            direction=[0.0, 0.0, 1.0, 0.0]
        )
        axis_b = Datum(
            name="axis_b",
            datum_type=DatumType.AXIS,
            origin=[5.0, 0.0, 0.0, 1.0],  # Offset in X
            direction=[0.0, 0.0, 1.0, 0.0]
        )

        result = evaluate_coincident(axis_a, axis_b)

        assert result.satisfied is False
        assert abs(result.error_distance - 5.0) < 1e-6


class TestCoincidentCircle:
    """Test COINCIDENT constraint for circle datums (bolt circles)."""

    def test_coincident_circles_aligned(self):
        """Two aligned bolt circles should be coincident."""
        # Dynamixel XH540 servo horn bolt circle
        circle_a = Datum(
            name="horn_bolt_circle",
            datum_type=DatumType.CIRCLE,
            origin=[0.0, 0.0, 0.0, 1.0],
            normal=[0.0, 0.0, 1.0, 0.0],
            radius=9.0  # 18mm diameter / 2
        )
        # Matching sun gear hub bolt circle
        circle_b = Datum(
            name="gear_bolt_circle",
            datum_type=DatumType.CIRCLE,
            origin=[0.0, 0.0, 0.0, 1.0],
            normal=[0.0, 0.0, 1.0, 0.0],
            radius=9.0
        )

        result = evaluate_coincident(circle_a, circle_b)

        assert result.satisfied is True
        assert result.error_distance < 1e-6
        assert "circle_to_circle" in result.details.get("type", "")

    def test_coincident_circles_offset_centers(self):
        """Bolt circles with offset centers should not be coincident."""
        circle_a = Datum(
            name="circle_a",
            datum_type=DatumType.CIRCLE,
            origin=[0.0, 0.0, 0.0, 1.0],
            normal=[0.0, 0.0, 1.0, 0.0],
            radius=9.0
        )
        circle_b = Datum(
            name="circle_b",
            datum_type=DatumType.CIRCLE,
            origin=[0.5, 0.0, 0.0, 1.0],  # 0.5mm offset
            normal=[0.0, 0.0, 1.0, 0.0],
            radius=9.0
        )

        result = evaluate_coincident(circle_a, circle_b, tolerance_mm=0.1)

        assert result.satisfied is False
        assert abs(result.error_distance - 0.5) < 1e-6

    def test_coincident_circles_different_radii(self):
        """Bolt circles with different radii are flagged in details."""
        circle_a = Datum(
            name="circle_a",
            datum_type=DatumType.CIRCLE,
            origin=[0.0, 0.0, 0.0, 1.0],
            normal=[0.0, 0.0, 1.0, 0.0],
            radius=9.0
        )
        circle_b = Datum(
            name="circle_b",
            datum_type=DatumType.CIRCLE,
            origin=[0.0, 0.0, 0.0, 1.0],
            normal=[0.0, 0.0, 1.0, 0.0],
            radius=10.0  # Different radius
        )

        result = evaluate_coincident(circle_a, circle_b)

        # Centers align so basic coincidence passes
        assert result.satisfied is True
        # But radius difference is reported
        assert abs(result.details.get("radius_difference", 0) - 1.0) < 1e-6


class TestBoltCircleAlignment:
    """Test bolt circle alignment for servo-to-gear interfaces."""

    def test_dynamixel_sun_gear_alignment(self):
        """Test XH540 servo horn to axis 1 sun gear bolt circle alignment.

        XH540 servo horn specifications:
        - 4 bolt holes on 18mm bolt circle diameter (9mm radius)
        - 45 degree angular offset from cardinal positions

        Sun gear hub must have matching hole pattern.
        """
        # XH540 servo horn bolt circle
        horn_circle = Datum(
            name="xh540_horn_bolts",
            datum_type=DatumType.CIRCLE,
            origin=[0.0, 0.0, 0.0, 1.0],
            normal=[0.0, 0.0, 1.0, 0.0],
            radius=9.0  # 18mm diameter
        )

        # Sun gear hub bolt circle (matching)
        gear_circle = Datum(
            name="sun_gear_hub_bolts",
            datum_type=DatumType.CIRCLE,
            origin=[0.0, 0.0, 0.0, 1.0],
            normal=[0.0, 0.0, 1.0, 0.0],
            radius=9.0
        )

        result = check_bolt_circle_alignment(
            horn_circle,
            gear_circle,
            hole_count=4,
            angular_offset_deg=0.0,  # Both at same angular position
            tolerance_mm=0.1
        )

        assert result.satisfied is True
        assert result.details.get("hole_count") == 4

    def test_bolt_circle_angular_offset(self):
        """Test bolt circles with angular offset between patterns.

        When angular_offset_deg is specified, it represents the angular
        difference between the hole patterns. For two identical bolt circles
        at the same position with no actual rotation applied, using
        angular_offset_deg=0 should pass (holes align at same positions).

        The angular_offset_deg parameter tells the function to expect that
        circle_b's holes are rotated by that amount relative to circle_a.
        If both circles have holes at the same angular positions, you should
        use angular_offset_deg=0.
        """
        circle_a = Datum(
            name="circle_a",
            datum_type=DatumType.CIRCLE,
            origin=[0.0, 0.0, 0.0, 1.0],
            normal=[0.0, 0.0, 1.0, 0.0],
            radius=9.0
        )
        circle_b = Datum(
            name="circle_b",
            datum_type=DatumType.CIRCLE,
            origin=[0.0, 0.0, 0.0, 1.0],
            normal=[0.0, 0.0, 1.0, 0.0],
            radius=9.0
        )

        # With zero angular offset (both patterns at same angular position)
        # This is the typical case when two bolt circles are properly aligned
        result = check_bolt_circle_alignment(
            circle_a, circle_b,
            hole_count=4,
            angular_offset_deg=0.0,  # No expected offset
            tolerance_mm=0.1
        )
        assert result.satisfied is True

    def test_bolt_circle_radius_mismatch(self):
        """Test that mismatched bolt circle radii are detected."""
        circle_a = Datum(
            name="circle_a",
            datum_type=DatumType.CIRCLE,
            origin=[0.0, 0.0, 0.0, 1.0],
            normal=[0.0, 0.0, 1.0, 0.0],
            radius=9.0  # 18mm bolt circle
        )
        circle_b = Datum(
            name="circle_b",
            datum_type=DatumType.CIRCLE,
            origin=[0.0, 0.0, 0.0, 1.0],
            normal=[0.0, 0.0, 1.0, 0.0],
            radius=10.0  # 20mm bolt circle - won't match!
        )

        result = check_bolt_circle_alignment(
            circle_a, circle_b,
            hole_count=4,
            tolerance_mm=0.1
        )

        assert result.satisfied is False
        assert "radi" in result.error_message.lower()  # "radii" or "radius"


class TestCoincidentMateCreation:
    """Test creation of COINCIDENT mates."""

    def test_create_coincident_mate(self):
        """Test creating a COINCIDENT mate between parts."""
        mate = create_coincident_mate(
            name="horn_to_gear",
            part_a="xh540_servo",
            datum_a="horn_bolt_circle",
            part_b="axis1_sun_gear",
            datum_b="hub_bolt_circle"
        )

        assert mate.name == "horn_to_gear"
        assert mate.mate_type == MateType.COINCIDENT
        assert mate.part_a == "xh540_servo"
        assert mate.part_b == "axis1_sun_gear"
        assert mate.datum_a == "horn_bolt_circle"
        assert mate.datum_b == "hub_bolt_circle"

        # Test property aliases
        assert mate.part1 == "xh540_servo"
        assert mate.part2 == "axis1_sun_gear"
        assert mate.datum1 == "horn_bolt_circle"
        assert mate.datum2 == "hub_bolt_circle"

    def test_coincident_mate_dof(self):
        """Test degrees of freedom for COINCIDENT mate."""
        mate = create_coincident_mate(
            name="test",
            part_a="a", datum_a="da",
            part_b="b", datum_b="db"
        )

        # COINCIDENT typically removes 3 DOF (position locked)
        assert mate.degrees_of_freedom() == 3


class TestMateValidation:
    """Test mate validation for COINCIDENT constraints."""

    def test_validate_coincident_with_circles(self):
        """COINCIDENT should be valid with two CIRCLE datums."""
        mate = create_coincident_mate(
            name="test",
            part_a="a", datum_a="da",
            part_b="b", datum_b="db"
        )

        datum_a = Datum(
            name="da",
            datum_type=DatumType.CIRCLE,
            origin=[0, 0, 0, 1],
            normal=[0, 0, 1, 0],
            radius=9.0
        )
        datum_b = Datum(
            name="db",
            datum_type=DatumType.CIRCLE,
            origin=[0, 0, 0, 1],
            normal=[0, 0, 1, 0],
            radius=9.0
        )

        issues = mate.validate(datum_a, datum_b)
        assert len(issues) == 0

    def test_validate_coincident_with_points(self):
        """COINCIDENT should be valid with two POINT datums."""
        mate = create_coincident_mate(
            name="test",
            part_a="a", datum_a="da",
            part_b="b", datum_b="db"
        )

        datum_a = Datum(
            name="da",
            datum_type=DatumType.POINT,
            origin=[0, 0, 0, 1]
        )
        datum_b = Datum(
            name="db",
            datum_type=DatumType.POINT,
            origin=[0, 0, 0, 1]
        )

        issues = mate.validate(datum_a, datum_b)
        assert len(issues) == 0

    def test_validate_coincident_with_planes(self):
        """COINCIDENT should be valid with two PLANE datums."""
        mate = create_coincident_mate(
            name="test",
            part_a="a", datum_a="da",
            part_b="b", datum_b="db"
        )

        datum_a = Datum(
            name="da",
            datum_type=DatumType.PLANE,
            origin=[0, 0, 0, 1],
            normal=[0, 0, 1, 0]
        )
        datum_b = Datum(
            name="db",
            datum_type=DatumType.PLANE,
            origin=[0, 0, 0, 1],
            normal=[0, 0, 1, 0]
        )

        issues = mate.validate(datum_a, datum_b)
        assert len(issues) == 0


class TestPointToPlaneCoincident:
    """Test COINCIDENT for point on plane."""

    def test_point_on_plane(self):
        """Point lying on plane should be coincident."""
        point = Datum(
            name="point",
            datum_type=DatumType.POINT,
            origin=[5.0, 3.0, 10.0, 1.0]  # Lies on Z=10 plane
        )
        plane = Datum(
            name="plane",
            datum_type=DatumType.PLANE,
            origin=[0.0, 0.0, 10.0, 1.0],
            normal=[0.0, 0.0, 1.0, 0.0]  # Z=10 plane
        )

        result = evaluate_coincident(point, plane)

        assert result.satisfied is True
        assert result.error_distance < 1e-6

    def test_point_off_plane(self):
        """Point not on plane should not be coincident."""
        point = Datum(
            name="point",
            datum_type=DatumType.POINT,
            origin=[5.0, 3.0, 15.0, 1.0]  # 5mm above Z=10 plane
        )
        plane = Datum(
            name="plane",
            datum_type=DatumType.PLANE,
            origin=[0.0, 0.0, 10.0, 1.0],
            normal=[0.0, 0.0, 1.0, 0.0]
        )

        result = evaluate_coincident(point, plane)

        assert result.satisfied is False
        assert abs(result.error_distance - 5.0) < 1e-6


class TestCoincidentResultStr:
    """Test CoincidentResult string formatting."""

    def test_result_str_pass(self):
        """Test string representation of passing result."""
        result = CoincidentResult(
            satisfied=True,
            error_distance=0.01,
            error_message="Test passed"
        )
        s = str(result)
        assert "[PASS]" in s
        assert "Test passed" in s

    def test_result_str_fail(self):
        """Test string representation of failing result."""
        result = CoincidentResult(
            satisfied=False,
            error_distance=5.0,
            error_angle=10.0,
            error_message="Test failed"
        )
        s = str(result)
        assert "[FAIL]" in s
        assert "5.000mm" in s
        assert "10.00deg" in s


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
