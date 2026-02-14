#!/usr/bin/env python3
"""Demonstration of yapCAD assembly constraint system.

This example shows how to use design constraints to validate that an
assembly meets high-level design intent. Constraints check orientation
and positioning requirements that are difficult to express with mates alone.

Example use case: A tube-climbing robot with motors on a wheel path.
The motors must:
1. Be positioned at a specific radius from the chassis center (AT_RADIUS)
2. Have their rotation axis tangent to the wheel path (TANGENT_TO_CIRCLE)
3. Have their tire treads facing outward toward the tube wall (RADIAL_FROM_CENTER)
4. Have their stator mounting faces pointing inward (RADIAL_FROM_CENTER inward)
"""

from yapcad.assembly import (
    Datum, DatumType,
    Constraint, ConstraintType, ConstraintResult
)

# Define the robot chassis specification
CHASSIS_CENTER = (0.0, 0.0, 0.0)
WHEEL_PATH_RADIUS = 124.5  # mm - radius at which motors are positioned
NUM_MOTORS = 3  # Three motors equally spaced around chassis

print("=" * 80)
print("Assembly Constraint Validation Demo: Tube-Climbing Robot")
print("=" * 80)


def validate_motor_placement(motor_position, motor_axis_direction,
                             tire_normal, stator_normal, motor_id):
    """Validate that a motor is correctly positioned and oriented.

    Args:
        motor_position: (x, y, z) position of motor center
        motor_axis_direction: (x, y, z) direction of motor rotation axis
        tire_normal: (x, y, z) normal vector of tire tread surface
        stator_normal: (x, y, z) normal vector of stator mounting face
        motor_id: Motor identifier (1, 2, 3)

    Returns:
        List of constraint results
    """
    print(f"\n{'='*80}")
    print(f"Motor {motor_id} Validation")
    print(f"{'='*80}")
    print(f"Position: ({motor_position[0]:.1f}, {motor_position[1]:.1f}, "
          f"{motor_position[2]:.1f})")
    print(f"Axis direction: ({motor_axis_direction[0]:.1f}, "
          f"{motor_axis_direction[1]:.1f}, {motor_axis_direction[2]:.1f})")

    results = []

    # Constraint 1: Motor must be at correct radius from center
    print(f"\n1. Checking position (AT_RADIUS)...")
    motor_position_datum = Datum(
        name=f"motor_{motor_id}_center",
        datum_type=DatumType.POINT,
        origin=[motor_position[0], motor_position[1], motor_position[2], 1],
        description=f"Motor {motor_id} center point"
    )

    at_radius_constraint = Constraint(
        name=f"motor_{motor_id}_at_radius",
        constraint_type=ConstraintType.AT_RADIUS,
        part=f"DDSM115_MOTOR_{motor_id}",
        datum=f"motor_{motor_id}_center",
        center=CHASSIS_CENTER,
        radius=WHEEL_PATH_RADIUS,
        tolerance_mm=0.5,
        description=f"Motor {motor_id} positioned on wheel path circle"
    )

    result = at_radius_constraint.evaluate(motor_position_datum)
    print(f"   {result}")
    results.append(result)

    # Constraint 2: Motor axis must be tangent (perpendicular to radial)
    print(f"\n2. Checking axis tangency (TANGENT_TO_CIRCLE)...")
    motor_axis_datum = Datum(
        name=f"motor_{motor_id}_axis",
        datum_type=DatumType.AXIS,
        origin=[motor_position[0], motor_position[1], motor_position[2], 1],
        direction=[motor_axis_direction[0], motor_axis_direction[1],
                  motor_axis_direction[2], 0],
        description=f"Motor {motor_id} rotation axis"
    )

    tangent_constraint = Constraint(
        name=f"motor_{motor_id}_tangent",
        constraint_type=ConstraintType.TANGENT_TO_CIRCLE,
        part=f"DDSM115_MOTOR_{motor_id}",
        datum=f"motor_{motor_id}_axis",
        center=CHASSIS_CENTER,
        radius=WHEEL_PATH_RADIUS,
        tolerance_deg=1.0,
        description=f"Motor {motor_id} axis tangent for rolling motion"
    )

    result = tangent_constraint.evaluate(motor_axis_datum)
    print(f"   {result}")
    results.append(result)

    # Constraint 3: Tire tread must face outward (toward tube wall)
    print(f"\n3. Checking tire orientation (RADIAL_FROM_CENTER outward)...")
    tire_face_datum = Datum(
        name=f"motor_{motor_id}_tire_face",
        datum_type=DatumType.PLANE,
        origin=[motor_position[0], motor_position[1], motor_position[2], 1],
        normal=[tire_normal[0], tire_normal[1], tire_normal[2], 0],
        description=f"Motor {motor_id} tire tread surface"
    )

    tire_radial_constraint = Constraint(
        name=f"motor_{motor_id}_tire_outward",
        constraint_type=ConstraintType.RADIAL_FROM_CENTER,
        part=f"DDSM115_MOTOR_{motor_id}",
        datum=f"motor_{motor_id}_tire_face",
        center=CHASSIS_CENTER,
        direction="outward",
        tolerance_deg=5.0,
        description=f"Motor {motor_id} tire faces tube wall for traction"
    )

    result = tire_radial_constraint.evaluate(tire_face_datum)
    print(f"   {result}")
    results.append(result)

    # Constraint 4: Stator must face inward (toward mounting bracket)
    print(f"\n4. Checking stator orientation (RADIAL_FROM_CENTER inward)...")
    stator_face_datum = Datum(
        name=f"motor_{motor_id}_stator_face",
        datum_type=DatumType.PLANE,
        origin=[motor_position[0], motor_position[1], motor_position[2], 1],
        normal=[stator_normal[0], stator_normal[1], stator_normal[2], 0],
        description=f"Motor {motor_id} stator mounting surface"
    )

    stator_radial_constraint = Constraint(
        name=f"motor_{motor_id}_stator_inward",
        constraint_type=ConstraintType.RADIAL_FROM_CENTER,
        part=f"DDSM115_MOTOR_{motor_id}",
        datum=f"motor_{motor_id}_stator_face",
        center=CHASSIS_CENTER,
        direction="inward",
        tolerance_deg=5.0,
        description=f"Motor {motor_id} stator faces bracket for mounting"
    )

    result = stator_radial_constraint.evaluate(stator_face_datum)
    print(f"   {result}")
    results.append(result)

    return results


# Example 1: Correctly oriented motor at 0° (along +X axis)
print("\n" + "#" * 80)
print("# EXAMPLE 1: Correctly Oriented Motor")
print("#" * 80)

results_correct = validate_motor_placement(
    motor_position=(124.5, 0.0, 50.0),      # At correct radius
    motor_axis_direction=(0.0, 1.0, 0.0),   # Tangent (perpendicular to radial)
    tire_normal=(1.0, 0.0, 0.0),            # Outward
    stator_normal=(-1.0, 0.0, 0.0),         # Inward
    motor_id=1
)

all_passed = all(r.passed for r in results_correct)
print(f"\n{'='*80}")
print(f"Motor 1 Overall: {'✅ PASS' if all_passed else '❌ FAIL'}")
print(f"{'='*80}")


# Example 2: Incorrectly oriented motor (common mistake: radial instead of tangent)
print("\n\n" + "#" * 80)
print("# EXAMPLE 2: Incorrectly Oriented Motor (COMMON ERROR)")
print("#" * 80)

results_incorrect = validate_motor_placement(
    motor_position=(124.5, 0.0, 50.0),      # Position is correct
    motor_axis_direction=(1.0, 0.0, 0.0),   # ❌ RADIAL (should be tangent!)
    tire_normal=(1.0, 0.0, 0.0),            # Outward - correct
    stator_normal=(-1.0, 0.0, 0.0),         # Inward - correct
    motor_id=2
)

all_passed = all(r.passed for r in results_incorrect)
print(f"\n{'='*80}")
print(f"Motor 2 Overall: {'✅ PASS' if all_passed else '❌ FAIL - Axis must be TANGENT, not RADIAL!'}")
print(f"{'='*80}")


# Example 3: Motor at wrong radius
print("\n\n" + "#" * 80)
print("# EXAMPLE 3: Motor at Wrong Radius")
print("#" * 80)

results_wrong_radius = validate_motor_placement(
    motor_position=(120.0, 0.0, 50.0),      # ❌ Wrong radius!
    motor_axis_direction=(0.0, 1.0, 0.0),   # Tangent - correct
    tire_normal=(1.0, 0.0, 0.0),            # Outward - correct
    stator_normal=(-1.0, 0.0, 0.0),         # Inward - correct
    motor_id=3
)

all_passed = all(r.passed for r in results_wrong_radius)
print(f"\n{'='*80}")
print(f"Motor 3 Overall: {'✅ PASS' if all_passed else '❌ FAIL - Wrong radius!'}")
print(f"{'='*80}")


# Summary
print("\n\n" + "=" * 80)
print("SUMMARY: Constraint System Benefits")
print("=" * 80)
print("""
The constraint system provides:

1. ✅ EXPLICIT DESIGN INTENT
   Constraints document what orientations are required (tangent vs radial)
   and why ("for rolling motion", "for traction")

2. ✅ AUTOMATIC VALIDATION
   Invalid assemblies are detected immediately at definition time,
   not after hours of debugging incorrect renders

3. ✅ CLEAR ERROR MESSAGES
   When constraints fail, you get precise geometric errors:
   "Tangent angle error: 90.00°" clearly indicates radial vs tangent issue

4. ✅ PARAMETRIC CORRECTNESS
   As design parameters change (wheel radius, motor positions),
   constraints ensure all parts remain correctly oriented

5. ✅ COMPOSABILITY
   Complex assemblies validate sub-assemblies, preventing integration issues

Implemented Constraint Types:
""")

for ct in ConstraintType:
    print(f"  - {ct.value}")

print("\n" + "=" * 80)
