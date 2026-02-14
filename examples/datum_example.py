#!/usr/bin/env python3
"""Example demonstrating datum features for the robot printer assembly.

This example shows how to define datum features on parts that can later
be used to create assembly mates and constraints. We define datums for
the motor and wheel bracket that will eventually be used to ensure:
  - Motor axis is TANGENT to chassis (not radial)
  - Stator face is FLUSH with bracket mounting surface
  - Motor axis is CONCENTRIC with bracket bore
"""

import sys
from pathlib import Path

# Add src to path for local development
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from yapcad.assembly.datum import Datum, DatumType, PartDefinition
from yapcad.geom import point
from yapcad.xform import Translation, Rotation


def define_motor_part():
    """Define the DDSM115 motor with its datum features."""
    print("=" * 60)
    print("Defining DDSM115 Motor Part")
    print("=" * 60)

    motor = PartDefinition(
        name="DDSM115_MOTOR",
        geometry_source="robot_printer/cots/ddsm115_motor.dsl",
        is_printable=False,
        material="Aluminum",
        description="DDSM115 servo motor with integrated encoder"
    )

    # Motor rotation axis (along Y in motor's local coordinates)
    motor.add_datum(Datum(
        name="motor_axis",
        datum_type=DatumType.AXIS,
        origin=point(0, 0, 0),
        direction=[0, 1, 0, 0],  # Y-axis
        description="Motor rotation axis - must be tangent to chassis",
        tags=["kinematic", "critical_orientation"]
    ))

    # Stator mounting face (facing +Y)
    motor.add_datum(Datum(
        name="stator_face",
        datum_type=DatumType.PLANE,
        origin=point(0, 10, 0),
        normal=[0, 1, 0, 0],  # Normal pointing +Y
        description="Stator mounting surface - must face bracket",
        tags=["mounting", "critical_orientation"]
    ))

    # Mounting hole bolt circle
    motor.add_datum(Datum(
        name="mounting_holes",
        datum_type=DatumType.CIRCLE,
        origin=point(0, 10, 0),
        normal=[0, 1, 0, 0],
        radius=15.2,
        description="M2.5 mounting hole pattern at 15.2mm radius",
        tags=["mounting", "fastener"]
    ))

    # Output shaft center point
    motor.add_datum(Datum(
        name="shaft_center",
        datum_type=DatumType.POINT,
        origin=point(0, -5, 0),
        description="Output shaft center for wheel mounting",
        tags=["kinematic"]
    ))

    print(f"\n✓ Created motor part with {len(motor.datums)} datums:")
    for name, datum in motor.datums.items():
        print(f"  - {name}: {datum.datum_type.value}")
        if datum.description:
            print(f"    {datum.description}")

    # Validate datums
    issues = motor.validate_datums()
    if issues:
        print("\n⚠ Validation warnings:")
        for issue in issues:
            print(f"  {issue}")
    else:
        print("\n✓ All datums valid")

    return motor


def define_bracket_part():
    """Define the wheel bracket with its datum features."""
    print("\n" + "=" * 60)
    print("Defining Wheel Bracket Part")
    print("=" * 60)

    bracket = PartDefinition(
        name="WHEEL_BRACKET",
        geometry_source="robot_printer/chassis/wheel_bracket.dsl",
        is_printable=True,
        material="PETG",
        description="Mounting bracket for drive wheel motor"
    )

    # Motor mounting face (facing -Z)
    bracket.add_datum(Datum(
        name="motor_interface",
        datum_type=DatumType.PLANE,
        origin=point(0, 0, 0),
        normal=[0, 0, -1, 0],  # Normal pointing -Z
        description="Motor mounting surface - mates with stator face",
        tags=["mounting"]
    ))

    # Motor bore axis (along Z)
    bracket.add_datum(Datum(
        name="bore_axis",
        datum_type=DatumType.AXIS,
        origin=point(0, 0, 0),
        direction=[0, 0, 1, 0],  # Z-axis
        description="Motor shaft bore - must be concentric with motor axis",
        tags=["kinematic"]
    ))

    # Mounting hole bolt circle (matches motor)
    bracket.add_datum(Datum(
        name="motor_bolt_holes",
        datum_type=DatumType.CIRCLE,
        origin=point(0, 0, 0),
        normal=[0, 0, -1, 0],
        radius=15.2,
        description="M2.5 bolt holes for motor mounting",
        tags=["mounting", "fastener"]
    ))

    # Chassis mounting plane
    bracket.add_datum(Datum(
        name="chassis_mount",
        datum_type=DatumType.PLANE,
        origin=point(20, 0, 0),
        normal=[1, 0, 0, 0],  # Normal pointing +X
        description="Surface that mounts to chassis pivot boss",
        tags=["mounting"]
    ))

    print(f"\n✓ Created bracket part with {len(bracket.datums)} datums:")
    for name, datum in bracket.datums.items():
        print(f"  - {name}: {datum.datum_type.value}")
        if datum.description:
            print(f"    {datum.description}")

    # Validate datums
    issues = bracket.validate_datums()
    if issues:
        print("\n⚠ Validation warnings:")
        for issue in issues:
            print(f"  {issue}")
    else:
        print("\n✓ All datums valid")

    return bracket


def demonstrate_transformation():
    """Demonstrate how datums transform with parts."""
    print("\n" + "=" * 60)
    print("Demonstrating Datum Transformation")
    print("=" * 60)

    # Create a simple axis datum
    original = Datum(
        name="test_axis",
        datum_type=DatumType.AXIS,
        origin=point(0, 0, 0),
        direction=[1, 0, 0, 0],  # X-axis
        description="Original axis along X"
    )

    print(f"\nOriginal datum:")
    print(f"  Origin: {[round(x, 2) for x in original.origin[:3]]}")
    print(f"  Direction: {[round(x, 2) for x in original.direction[:3]]}")

    # Apply translation
    T = Translation(point(10, 20, 30))
    translated = original.transform(T)

    print(f"\nAfter translation by [10, 20, 30]:")
    print(f"  Origin: {[round(x, 2) for x in translated.origin[:3]]}")
    print(f"  Direction: {[round(x, 2) for x in translated.direction[:3]]} (unchanged)")

    # Apply rotation (90° around Z-axis)
    R = Rotation([0, 0, 1, 0], 90)
    rotated = original.transform(R)

    print(f"\nAfter 90° rotation around Z-axis:")
    print(f"  Origin: {[round(x, 2) for x in rotated.origin[:3]]}")
    print(f"  Direction: {[round(x, 4) for x in rotated.direction[:3]]} (X->Y)")

    # Combined transformation
    combined = T.mul(R)
    both = original.transform(combined)

    print(f"\nAfter rotate then translate:")
    print(f"  Origin: {[round(x, 2) for x in both.origin[:3]]}")
    print(f"  Direction: {[round(x, 4) for x in both.direction[:3]]}")

    print("\n✓ Transformations preserve datum relationships")


def demonstrate_future_usage():
    """Show how datums will be used in future mate definitions."""
    print("\n" + "=" * 60)
    print("Future Usage: Assembly Mates (Not Yet Implemented)")
    print("=" * 60)

    print("""
Once the mate system is implemented, datums will be used like this:

    from yapcad.assembly import Assembly, Mate, MateType

    # Create assembly
    assembly = Assembly("wheel_drive_assembly")
    assembly.add_part(motor)
    assembly.add_part(bracket)

    # Define mates using datum names
    assembly.add_mate(Mate(
        name="stator_flush",
        mate_type=MateType.FLUSH,
        part1="DDSM115_MOTOR", datum1="stator_face",
        part2="WHEEL_BRACKET", datum2="motor_interface",
        description="Motor stator flush with bracket face"
    ))

    assembly.add_mate(Mate(
        name="axis_concentric",
        mate_type=MateType.CONCENTRIC,
        part1="DDSM115_MOTOR", datum1="motor_axis",
        part2="WHEEL_BRACKET", datum2="bore_axis",
        description="Motor shaft concentric with bore"
    ))

    # Add design constraint
    assembly.add_constraint(Constraint(
        constraint_type=ConstraintType.TANGENT_TO_CIRCLE,
        part1="DDSM115_MOTOR", datum1="motor_axis",
        part2="CHASSIS", datum2="wheel_orbit",
        description="Motor axis must be tangent for proper wheel rolling"
    ))

    # Solve and validate
    result = assembly.solve()
    if result.is_valid:
        print("✓ Assembly valid")
        transforms = result.transforms
    else:
        print("✗ Assembly invalid:")
        for error in result.errors:
            print(f"  - {error}")
    """)


def main():
    """Run the complete datum feature demonstration."""
    print("\n" + "=" * 60)
    print("yapCAD Datum Feature Example")
    print("Robot Printer Assembly System")
    print("=" * 60)

    # Define parts with datums
    motor = define_motor_part()
    bracket = define_bracket_part()

    # Demonstrate transformation
    demonstrate_transformation()

    # Show future usage
    demonstrate_future_usage()

    print("\n" + "=" * 60)
    print("Example Complete")
    print("=" * 60)
    print(f"\nDefined {len(motor.datums) + len(bracket.datums)} total datums")
    print("These datums will be used for mate and constraint definitions")
    print("once the complete assembly system is implemented.")
    print()


if __name__ == "__main__":
    main()
