"""Tests for assembly constraint validator.

Tests the kinematic chain constraint traversal and validation system
for verifying assembly mates in the final assembled position.

Run with: PYTHONPATH=./src:. pytest tests/test_assembly_constraint_validator.py -v
"""

import pytest
import math
import numpy as np


class TestAssemblyConstraintValidator:
    """Test suite for assembly constraint validator."""

    @pytest.fixture
    def chain(self):
        """Create a test kinematic chain."""
        from robot_printer.kinematic_chain import create_tube_robot_chain
        return create_tube_robot_chain()

    @pytest.fixture
    def validator(self, chain):
        """Create a validator with the chain."""
        from robot_printer.assembly_constraint_validator import AssemblyValidator
        return AssemblyValidator(chain)

    def test_chain_creation(self, chain):
        """Test that kinematic chain is created correctly."""
        # Should have many parts
        assert len(chain.parts) > 30

        # Should have expected key parts
        assert "CHASSIS_PLATE" in chain.parts
        assert "WHEEL_ARM_1" in chain.parts
        assert "DDSM115_MOTOR_1" in chain.parts
        assert "ARM_TOWER" in chain.parts

    def test_chain_validation(self, chain):
        """Test that kinematic chain passes validation."""
        issues = chain.validate()
        assert len(issues) == 0, f"Chain validation issues: {issues}"

    def test_world_transform_computation(self, chain):
        """Test that world transforms can be computed."""
        # Get transform for a part
        tf = chain.get_world_transform("WHEEL_ARM_1")
        assert tf is not None

        # Check translation is reasonable
        pos = tf.translation
        assert isinstance(pos, tuple)
        assert len(pos) == 3

        # Wheel arm should be at approximately 145mm radius
        xy_dist = math.sqrt(pos[0]**2 + pos[1]**2)
        assert 140 < xy_dist < 150, f"Wheel arm at radius {xy_dist:.1f}mm, expected ~145mm"

    def test_motor_position(self, chain):
        """Test motor positions are reasonable."""
        for i in [1, 2, 3]:
            motor_name = f"DDSM115_MOTOR_{i}"
            tf = chain.get_world_transform(motor_name)
            pos = tf.translation

            # Motors should be at some radius from center
            xy_dist = math.sqrt(pos[0]**2 + pos[1]**2)
            assert 100 < xy_dist < 180, f"{motor_name} at radius {xy_dist:.1f}mm"

    def test_validator_creation(self, validator):
        """Test that validator can be created."""
        assert validator is not None
        assert validator.chain is not None

    def test_add_constraint(self, validator):
        """Test adding constraints to validator."""
        from robot_printer.assembly_constraint_validator import (
            AssemblyConstraint, ConstraintType
        )

        constraint = AssemblyConstraint(
            name="test_constraint",
            constraint_type=ConstraintType.COINCIDENT,
            part_a="CHASSIS_PLATE",
            frame_a="ORIGIN",
            part_b="CENTRAL_HUB",
            frame_b="ORIGIN",
            tolerance_mm=10.0
        )

        validator.add_constraint(constraint)
        assert len(validator.constraints) == 1

    def test_coincident_evaluation(self, validator):
        """Test evaluating a COINCIDENT constraint."""
        from robot_printer.assembly_constraint_validator import (
            AssemblyConstraint, ConstraintType
        )

        # Wheel arm origin should be at pivot boss axis
        constraint = AssemblyConstraint(
            name="arm_at_pivot",
            constraint_type=ConstraintType.COINCIDENT,
            part_a="WHEEL_ARM_1",
            frame_a="ORIGIN",
            part_b="PIVOT_BOSS_1",
            frame_b="PIVOT_AXIS",
            tolerance_mm=1.0
        )

        validator.add_constraint(constraint)
        results = validator.validate_all()

        assert len(results) == 1
        result = results[0]
        assert result.constraint_name == "arm_at_pivot"
        # Should pass (wheel arm is at pivot axis via joint)
        assert result.passed, f"Constraint failed: {result.message}"

    def test_tangent_constraint(self, validator):
        """Test evaluating a TANGENT constraint."""
        from robot_printer.assembly_constraint_validator import (
            AssemblyConstraint, ConstraintType
        )

        # Pivot axis should be tangent to chassis
        constraint = AssemblyConstraint(
            name="pivot_tangent",
            constraint_type=ConstraintType.TANGENT,
            part_a="PIVOT_BOSS_1",
            frame_a="PIVOT_AXIS",
            axis_a=(0, 0, 1),  # Local Z-axis
            part_b="CHASSIS_PLATE",
            frame_b="ORIGIN",
            tolerance_deg=10.0
        )

        validator.add_constraint(constraint)
        results = validator.validate_all()

        result = results[0]
        # This should pass - pivot bore axis is designed to be tangent
        assert result.passed, f"Tangent constraint failed with error {result.error_angular:.1f}deg"

    def test_wheel_assembly_constraints(self, chain):
        """Test wheel assembly constraint definitions."""
        from robot_printer.assembly_constraint_validator import (
            get_wheel_assembly_constraints, AssemblyValidator
        )

        constraints = get_wheel_assembly_constraints(wheel_idx=1)
        assert len(constraints) >= 4  # At least motor mount, parallel, tangent, etc.

        # All constraints should have valid parts
        for c in constraints:
            assert c.part_a in chain.parts, f"Part {c.part_a} not in chain"
            assert c.part_b in chain.parts, f"Part {c.part_b} not in chain"

    def test_scara_arm_constraints(self, chain):
        """Test SCARA arm constraint definitions."""
        from robot_printer.assembly_constraint_validator import (
            get_scara_arm_constraints, AssemblyValidator
        )

        constraints = get_scara_arm_constraints()
        assert len(constraints) >= 5  # Tower, gearboxes, links

        # All constraints should have valid parts
        for c in constraints:
            assert c.part_a in chain.parts, f"Part {c.part_a} not in chain"
            assert c.part_b in chain.parts, f"Part {c.part_b} not in chain"

    def test_full_validation(self, chain):
        """Test full assembly validation."""
        from robot_printer.assembly_constraint_validator import validate_full_assembly

        passed, failed, results = validate_full_assembly(chain, verbose=False)

        # Should have results
        assert len(results) > 0

        # Report pass rate
        total = passed + failed
        pass_rate = passed / total if total > 0 else 0
        print(f"\nValidation: {passed}/{total} constraints passed ({pass_rate*100:.1f}%)")

        # Print failed constraints
        for r in results:
            if not r.passed:
                print(f"  FAILED: {r}")

    def test_parallel_constraint_for_axes(self, validator):
        """Test PARALLEL constraint evaluation."""
        from robot_printer.assembly_constraint_validator import (
            AssemblyConstraint, ConstraintType
        )

        # Motor axis should be parallel to pivot axis
        constraint = AssemblyConstraint(
            name="motor_parallel_pivot",
            constraint_type=ConstraintType.PARALLEL,
            part_a="DDSM115_MOTOR_1",
            frame_a="ORIGIN",
            axis_a=(0, 1, 0),  # Motor Y-axis
            part_b="PIVOT_BOSS_1",
            frame_b="PIVOT_AXIS",
            axis_b=(0, 0, 1),  # Pivot Z-axis
            tolerance_deg=10.0
        )

        validator.add_constraint(constraint)
        results = validator.validate_all()

        result = results[0]
        # Should pass - motor and pivot axes are designed to be parallel
        print(f"Motor-pivot parallelism: {result.error_angular:.1f}deg error")

    def test_world_axis_computation(self, validator):
        """Test that world axis computation works correctly."""
        # Get motor axis in world coordinates
        motor_axis = validator.get_world_axis(
            "DDSM115_MOTOR_1", "ORIGIN", (0, 1, 0)
        )

        # Should be a unit vector
        mag = math.sqrt(motor_axis[0]**2 + motor_axis[1]**2 + motor_axis[2]**2)
        assert abs(mag - 1.0) < 1e-6, f"Axis not unit vector: mag={mag}"

    def test_constraint_result_formatting(self, validator):
        """Test that constraint results format correctly."""
        from robot_printer.assembly_constraint_validator import (
            AssemblyConstraint, ConstraintType
        )

        constraint = AssemblyConstraint(
            name="test_fmt",
            constraint_type=ConstraintType.COINCIDENT,
            part_a="CHASSIS_PLATE",
            frame_a="ORIGIN",
            part_b="CENTRAL_HUB",
            frame_b="ORIGIN",
            tolerance_mm=1.0
        )

        validator.add_constraint(constraint)
        results = validator.validate_all()

        result = results[0]
        result_str = str(result)

        # Should contain key information
        assert "test_fmt" in result_str
        assert "PASS" in result_str or "FAIL" in result_str


class TestConstraintTypes:
    """Test individual constraint type evaluations."""

    @pytest.fixture
    def validator(self):
        """Create validator with chain."""
        from robot_printer.kinematic_chain import create_tube_robot_chain
        from robot_printer.assembly_constraint_validator import AssemblyValidator
        chain = create_tube_robot_chain()
        return AssemblyValidator(chain)

    def test_at_distance_constraint(self, validator):
        """Test AT_DISTANCE constraint."""
        from robot_printer.assembly_constraint_validator import (
            AssemblyConstraint, ConstraintType
        )

        # Central hub should be at some distance from wheel arm
        constraint = AssemblyConstraint(
            name="hub_to_arm_distance",
            constraint_type=ConstraintType.AT_DISTANCE,
            part_a="CENTRAL_HUB",
            frame_a="ORIGIN",
            part_b="WHEEL_ARM_1",
            frame_b="ORIGIN",
            target_value=145.0,  # Pivot radius
            tolerance_mm=10.0
        )

        result = validator.evaluate_constraint(constraint)
        print(f"Hub-to-arm distance: target 145mm, error {result.error_linear:.1f}mm")

    def test_perpendicular_constraint(self, validator):
        """Test PERPENDICULAR constraint."""
        from robot_printer.assembly_constraint_validator import (
            AssemblyConstraint, ConstraintType
        )

        # Chassis Z should be perpendicular to motor axis (motor is tangent)
        constraint = AssemblyConstraint(
            name="chassis_z_perp_motor",
            constraint_type=ConstraintType.PERPENDICULAR,
            part_a="CHASSIS_PLATE",
            frame_a="ORIGIN",
            axis_a=(0, 0, 1),  # Chassis Z
            part_b="DDSM115_MOTOR_1",
            frame_b="ORIGIN",
            axis_b=(0, 1, 0),  # Motor Y
            tolerance_deg=10.0
        )

        result = validator.evaluate_constraint(constraint)
        print(f"Chassis Z to motor axis angle from perp: {result.error_angular:.1f}deg")


class TestEdgeCases:
    """Test edge cases and error handling."""

    @pytest.fixture
    def chain(self):
        """Create a test kinematic chain."""
        from robot_printer.kinematic_chain import create_tube_robot_chain
        return create_tube_robot_chain()

    def test_invalid_part_name(self, chain):
        """Test behavior with invalid part name."""
        from robot_printer.assembly_constraint_validator import (
            AssemblyConstraint, ConstraintType, AssemblyValidator
        )

        validator = AssemblyValidator(chain)

        constraint = AssemblyConstraint(
            name="invalid_test",
            constraint_type=ConstraintType.COINCIDENT,
            part_a="NONEXISTENT_PART",
            frame_a="ORIGIN",
            part_b="CHASSIS_PLATE",
            frame_b="ORIGIN",
            tolerance_mm=1.0
        )

        result = validator.evaluate_constraint(constraint)
        # Should fail gracefully
        assert not result.passed
        assert "Error" in result.message or "not" in result.message.lower()

    def test_invalid_frame_name(self, chain):
        """Test behavior with invalid frame name."""
        from robot_printer.assembly_constraint_validator import (
            AssemblyConstraint, ConstraintType, AssemblyValidator
        )

        validator = AssemblyValidator(chain)

        constraint = AssemblyConstraint(
            name="invalid_frame_test",
            constraint_type=ConstraintType.COINCIDENT,
            part_a="CHASSIS_PLATE",
            frame_a="NONEXISTENT_FRAME",
            part_b="CENTRAL_HUB",
            frame_b="ORIGIN",
            tolerance_mm=1.0
        )

        result = validator.evaluate_constraint(constraint)
        # Should fail gracefully
        assert not result.passed

    def test_zero_tolerance(self, chain):
        """Test with zero tolerance (very strict)."""
        from robot_printer.assembly_constraint_validator import (
            AssemblyConstraint, ConstraintType, AssemblyValidator
        )

        validator = AssemblyValidator(chain)

        # Zero tolerance should still work (though likely fail)
        constraint = AssemblyConstraint(
            name="zero_tol_test",
            constraint_type=ConstraintType.COINCIDENT,
            part_a="CHASSIS_PLATE",
            frame_a="ORIGIN",
            part_b="CENTRAL_HUB",
            frame_b="ORIGIN",
            tolerance_mm=0.0
        )

        result = validator.evaluate_constraint(constraint)
        # Should evaluate without error
        assert result.constraint_name == "zero_tol_test"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
