"""
Example usage of the yapCAD assembly mate system.

This file demonstrates how to define kinematic constraints for mechanical
assemblies using the mate module.
"""

import math
from yapcad.assembly.mate import (
    Mate, MateType, MateLimits, MateDynamics,
    create_revolute_mate, create_prismatic_mate, create_gear_mate
)


def example_robot_arm():
    """
    Define a simple 2-DOF robot arm with revolute joints.

    The arm consists of:
    - Base (fixed)
    - Shoulder joint (revolute, ±90°)
    - Upper arm
    - Elbow joint (revolute, 0-180°)
    - Forearm
    """
    print("Example 1: Robot Arm Assembly")
    print("=" * 60)

    # Shoulder joint - allows ±90 degree pitch motion
    shoulder = create_revolute_mate(
        name="shoulder_pitch",
        part_a="robot_base",
        datum_a="shoulder_mount_axis",
        part_b="upper_arm",
        datum_b="arm_root_axis",
        axis=[0, 1, 0, 0],  # Y-axis rotation
        min_angle=-math.pi/2,
        max_angle=math.pi/2,
        max_velocity=1.5,
        max_torque=100.0,
        friction=0.05,
        damping=0.2
    )

    print(f"Shoulder joint: {shoulder.name}")
    print(f"  Type: {shoulder.mate_type.value}")
    print(f"  DOF: {shoulder.degrees_of_freedom()}")
    print(f"  Range: {math.degrees(shoulder.limits.min_value):.1f}° to "
          f"{math.degrees(shoulder.limits.max_value):.1f}°")
    print(f"  Max velocity: {shoulder.limits.max_velocity} rad/s")
    print(f"  Friction: {shoulder.dynamics.friction_static}")
    print()

    # Elbow joint - allows 0-180 degree flexion
    elbow = create_revolute_mate(
        name="elbow_flex",
        part_a="upper_arm",
        datum_a="elbow_axis",
        part_b="forearm",
        datum_b="forearm_root_axis",
        axis=[0, 1, 0, 0],  # Y-axis rotation
        min_angle=0.0,
        max_angle=math.pi,
        max_velocity=2.0,
        max_torque=50.0,
        friction=0.03,
        damping=0.15
    )

    print(f"Elbow joint: {elbow.name}")
    print(f"  Type: {elbow.mate_type.value}")
    print(f"  DOF: {elbow.degrees_of_freedom()}")
    print(f"  Range: {math.degrees(elbow.limits.min_value):.1f}° to "
          f"{math.degrees(elbow.limits.max_value):.1f}°")
    print()

    # Test constraint evaluation
    test_position = math.pi / 4  # 45 degrees
    result = shoulder.evaluate(current_position=test_position)
    print(f"Testing shoulder at {math.degrees(test_position):.1f}°:")
    print(f"  Satisfied: {result['satisfied']}")
    print(f"  DOF remaining: {result['dof_remaining']}")
    print()

    # Test limit violation
    test_position = math.pi  # 180 degrees - exceeds limit
    result = shoulder.evaluate(current_position=test_position)
    print(f"Testing shoulder at {math.degrees(test_position):.1f}° (beyond limit):")
    print(f"  Satisfied: {result['satisfied']}")
    print(f"  Error: {result['error']:.4f} radians")
    print(f"  Violations: {result['limit_violations']}")
    print()


def example_gearbox():
    """
    Define a 3-stage planetary gearbox with gear mates.

    Each stage has a 3:1 reduction, resulting in 27:1 overall.
    """
    print("Example 2: Planetary Gearbox")
    print("=" * 60)

    # Stage 1: Motor to first sun gear (direct connection, 1:1)
    motor_to_sun1 = Mate(
        name="motor_to_sun1",
        mate_type=MateType.RIGID,
        part_a="motor_shaft",
        datum_a="output_axis",
        part_b="sun_gear_1",
        datum_b="gear_axis"
    )
    print(f"Motor connection: {motor_to_sun1.name}")
    print(f"  Type: {motor_to_sun1.mate_type.value}")
    print(f"  DOF: {motor_to_sun1.degrees_of_freedom()}")
    print()

    # Stage 1: Sun to carrier (3:1 reduction)
    stage1 = create_gear_mate(
        name="sun1_to_carrier1",
        part_a="sun_gear_1",
        datum_a="gear_axis",
        part_b="carrier_1",
        datum_b="carrier_axis",
        ratio=1.0/3.0  # Carrier rotates 1/3 speed of sun
    )
    print(f"Stage 1 reduction: {stage1.name}")
    print(f"  Ratio: {stage1.coupling_ratio} (3:1 reduction)")
    print(f"  Is coupled: {stage1.is_coupled()}")

    # Test coupling computation
    input_angle = math.pi  # 180 degrees input
    output_angle = stage1.compute_coupled_motion(input_angle)
    print(f"  Input: {math.degrees(input_angle):.1f}° -> "
          f"Output: {math.degrees(output_angle):.1f}°")
    print()

    # Stage 2: Carrier 1 to sun 2 (rigid connection)
    carrier1_to_sun2 = create_gear_mate(
        name="carrier1_to_sun2",
        part_a="carrier_1",
        datum_a="output_axis",
        part_b="sun_gear_2",
        datum_b="gear_axis",
        ratio=1.0
    )

    # Stage 2: Sun 2 to carrier 2 (3:1 reduction)
    stage2 = create_gear_mate(
        name="sun2_to_carrier2",
        part_a="sun_gear_2",
        datum_a="gear_axis",
        part_b="carrier_2",
        datum_b="carrier_axis",
        ratio=1.0/3.0
    )
    print(f"Stage 2 reduction: {stage2.name}")
    print(f"  Ratio: {stage2.coupling_ratio} (3:1 reduction)")
    print()

    # Stage 3: Similar to stages 1 and 2
    stage3 = create_gear_mate(
        name="sun3_to_carrier3",
        part_a="sun_gear_3",
        datum_a="gear_axis",
        part_b="carrier_3",
        datum_b="carrier_axis",
        ratio=1.0/3.0
    )
    print(f"Stage 3 reduction: {stage3.name}")
    print(f"  Ratio: {stage3.coupling_ratio} (3:1 reduction)")

    # Calculate overall reduction
    overall_ratio = (1.0/3.0) ** 3
    print(f"  Overall reduction: {1.0/overall_ratio:.1f}:1")
    print()


def example_linear_actuator():
    """
    Define a linear actuator with prismatic joint and screw coupling.
    """
    print("Example 3: Linear Actuator")
    print("=" * 60)

    # Motor to leadscrew (rigid connection)
    motor_coupling = Mate(
        name="motor_to_leadscrew",
        mate_type=MateType.RIGID,
        part_a="motor",
        datum_a="output_shaft",
        part_b="leadscrew",
        datum_b="screw_axis"
    )
    print(f"Motor coupling: {motor_coupling.name}")
    print(f"  Type: {motor_coupling.mate_type.value}")
    print()

    # Leadscrew to carriage (screw coupling, 2mm pitch)
    screw_coupling = Mate(
        name="leadscrew_to_carriage",
        mate_type=MateType.SCREW,
        part_a="leadscrew",
        datum_a="screw_axis",
        part_b="carriage",
        datum_b="nut_axis",
        coupling_pitch=2.0  # 2mm per revolution
    )
    print(f"Screw coupling: {screw_coupling.name}")
    print(f"  Type: {screw_coupling.mate_type.value}")
    print(f"  Pitch: {screw_coupling.coupling_pitch} mm/rev")
    print(f"  Is coupled: {screw_coupling.is_coupled()}")

    # Test motion calculation
    rotations = 10  # 10 revolutions
    input_angle = rotations * 2.0 * math.pi
    linear_motion = screw_coupling.compute_coupled_motion(input_angle)
    print(f"  {rotations} revolutions -> {linear_motion:.1f} mm linear travel")
    print()

    # Carriage to rail (prismatic joint with limits)
    rail_guide = create_prismatic_mate(
        name="carriage_on_rail",
        part_a="base",
        datum_a="rail_axis",
        part_b="carriage",
        datum_b="slider_axis",
        axis=[0, 0, 1, 0],  # Z-axis travel
        min_position=0.0,
        max_position=200.0,  # 200mm stroke
        max_velocity=100.0,  # 100 mm/s
        max_force=500.0,     # 500N force limit
        friction=0.02,
        damping=10.0
    )
    print(f"Rail guide: {rail_guide.name}")
    print(f"  Type: {rail_guide.mate_type.value}")
    print(f"  DOF: {rail_guide.degrees_of_freedom()}")
    print(f"  Stroke: {rail_guide.limits.max_value} mm")
    print(f"  Max speed: {rail_guide.limits.max_velocity} mm/s")
    print()


def example_spherical_joint():
    """
    Define a spherical (ball-and-socket) joint with multi-axis limits.
    """
    print("Example 4: Ball-and-Socket Joint")
    print("=" * 60)

    ball_joint = Mate(
        name="shoulder_ball",
        mate_type=MateType.SPHERICAL,
        part_a="torso",
        datum_a="shoulder_socket",
        part_b="upper_arm",
        datum_b="ball_center",
        limits=MateLimits(
            min_value=-math.pi/3,      # Primary axis: ±60°
            max_value=math.pi/3,
            min_secondary=-math.pi/4,  # Secondary axis: ±45°
            max_secondary=math.pi/4,
            max_velocity=3.0
        ),
        dynamics=MateDynamics(
            friction_static=0.1,
            damping=0.5
        )
    )

    print(f"Spherical joint: {ball_joint.name}")
    print(f"  Type: {ball_joint.mate_type.value}")
    print(f"  DOF: {ball_joint.degrees_of_freedom()}")
    print(f"  Primary range: ±{math.degrees(ball_joint.limits.max_value):.1f}°")
    print(f"  Secondary range: ±{math.degrees(ball_joint.limits.max_secondary):.1f}°")
    print()


def example_serialization():
    """
    Demonstrate mate serialization to/from dictionary for JSON export.
    """
    print("Example 5: Mate Serialization")
    print("=" * 60)

    # Create a mate with full configuration
    original = create_revolute_mate(
        name="test_joint",
        part_a="part_a",
        datum_a="axis_a",
        part_b="part_b",
        datum_b="axis_b",
        min_angle=-math.pi/2,
        max_angle=math.pi/2,
        max_velocity=2.0,
        max_torque=50.0,
        friction=0.05,
        damping=0.2
    )

    print("Original mate:")
    print(f"  Name: {original.name}")
    print(f"  Type: {original.mate_type.value}")
    print()

    # Serialize to dictionary
    mate_dict = original.to_dict()
    print("Serialized dictionary keys:")
    for key in mate_dict.keys():
        print(f"  - {key}")
    print()

    # Deserialize from dictionary
    restored = Mate.from_dict(mate_dict)
    print("Restored mate:")
    print(f"  Name: {restored.name}")
    print(f"  Type: {restored.mate_type.value}")
    matches = (restored.name == original.name and
               restored.mate_type == original.mate_type)
    print(f"  Matches original: {matches}")
    print()


def main():
    """Run all examples."""
    example_robot_arm()
    print("\n")
    example_gearbox()
    print("\n")
    example_linear_actuator()
    print("\n")
    example_spherical_joint()
    print("\n")
    example_serialization()


if __name__ == "__main__":
    main()
