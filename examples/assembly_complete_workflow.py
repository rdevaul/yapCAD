#!/usr/bin/env python3
"""Complete assembly workflow demonstrating datum, mate, and constraint integration.

This example shows how the three assembly modules work together:
1. Datum - Define geometric features on parts
2. Mate - Position parts relative to each other
3. Constraint - Validate that assembly meets design intent

The example creates a simplified wheel assembly with:
- Motor positioned at wheel path radius
- Motor axis tangent to wheel path for rolling
- Tire facing outward for traction
"""

from yapcad.assembly import (
    Datum, DatumType, PartDefinition,
    Mate, MateType,
    Constraint, ConstraintType
)

print("=" * 80)
print("Complete Assembly Workflow: Datum + Mate + Constraint")
print("=" * 80)

# Step 1: Define parts with datum features
print("\n1. Defining parts with datum features...")
print("-" * 80)

# Motor part definition
motor = PartDefinition(
    name="DDSM115_MOTOR",
    geometry_source="robot_printer/cots/ddsm115_motor.dsl",
    is_printable=False,
    material="Aluminum"
)

# Add motor datums
motor.add_datum(Datum(
    name="motor_axis",
    datum_type=DatumType.AXIS,
    origin=[0.0, 0.0, 0.0, 1],
    direction=[0.0, 1.0, 0.0, 0],
    description="Motor rotation axis"
))

motor.add_datum(Datum(
    name="stator_face",
    datum_type=DatumType.PLANE,
    origin=[0.0, -5.0, 0.0, 1],
    normal=[0.0, -1.0, 0.0, 0],
    description="Stator mounting surface"
))

motor.add_datum(Datum(
    name="tire_face",
    datum_type=DatumType.PLANE,
    origin=[0.0, 10.0, 0.0, 1],
    normal=[0.0, 1.0, 0.0, 0],
    description="Tire tread surface"
))

print(f"✓ Motor defined with {len(motor.datums)} datums:")
for datum_name in motor.datums:
    datum = motor.get_datum(datum_name)
    print(f"    - {datum_name}: {datum.datum_type.value}")

# Bracket part definition
bracket = PartDefinition(
    name="WHEEL_BRACKET",
    geometry_source="robot_printer/chassis/wheel_arm.dsl",
    is_printable=True,
    material="PLA"
)

bracket.add_datum(Datum(
    name="motor_mount_face",
    datum_type=DatumType.PLANE,
    origin=[124.5, 0.0, 50.0, 1],  # At wheel path radius
    normal=[1.0, 0.0, 0.0, 0],     # Facing outward
    description="Motor mounting surface"
))

bracket.add_datum(Datum(
    name="motor_bore_axis",
    datum_type=DatumType.AXIS,
    origin=[124.5, 0.0, 50.0, 1],
    direction=[0.0, 1.0, 0.0, 0],  # Tangent direction
    description="Motor shaft bore"
))

print(f"✓ Bracket defined with {len(bracket.datums)} datums:")
for datum_name in bracket.datums:
    datum = bracket.get_datum(datum_name)
    print(f"    - {datum_name}: {datum.datum_type.value}")


# Step 2: Define mates to position parts
print("\n2. Defining mates to position motor on bracket...")
print("-" * 80)

# Mate 1: Stator face coincident with bracket mount face
flush_mate = Mate(
    name="motor_flush_mount",
    mate_type=MateType.COINCIDENT,
    part_a="DDSM115_MOTOR",
    datum_a="stator_face",
    part_b="WHEEL_BRACKET",
    datum_b="motor_mount_face"
)

print(f"✓ {flush_mate.name}: {flush_mate.mate_type.value}")
print(f"    {flush_mate.part_a}.{flush_mate.datum_a} ↔ {flush_mate.part_b}.{flush_mate.datum_b}")

# Mate 2: Motor axis concentric with bracket bore
concentric_mate = Mate(
    name="motor_axis_align",
    mate_type=MateType.CONCENTRIC,
    part_a="DDSM115_MOTOR",
    datum_a="motor_axis",
    part_b="WHEEL_BRACKET",
    datum_b="motor_bore_axis"
)

print(f"✓ {concentric_mate.name}: {concentric_mate.mate_type.value}")
print(f"    {concentric_mate.part_a}.{concentric_mate.datum_a} ↔ "
      f"{concentric_mate.part_b}.{concentric_mate.datum_b}")


# Step 3: Define constraints to validate design intent
print("\n3. Defining constraints to validate assembly...")
print("-" * 80)

constraints = []

# Constraint 1: Motor at correct radius
c1 = Constraint(
    name="motor_at_wheel_radius",
    constraint_type=ConstraintType.AT_RADIUS,
    part="DDSM115_MOTOR",
    datum="motor_axis",
    center=(0.0, 0.0, 0.0),
    radius=124.5,
    tolerance_mm=0.5,
    description="Motor positioned on wheel path circle"
)
constraints.append(c1)
print(f"✓ {c1.name}: {c1.constraint_type.value}")
print(f"    Validates: Motor at radius {c1.radius}mm from center")

# Constraint 2: Motor axis tangent to wheel path
c2 = Constraint(
    name="motor_axis_tangent",
    constraint_type=ConstraintType.TANGENT_TO_CIRCLE,
    part="DDSM115_MOTOR",
    datum="motor_axis",
    center=(0.0, 0.0, 0.0),
    radius=124.5,
    tolerance_deg=1.0,
    description="Motor axis must be tangent for rolling motion"
)
constraints.append(c2)
print(f"✓ {c2.name}: {c2.constraint_type.value}")
print(f"    Validates: Axis tangent (not radial) to wheel path")

# Constraint 3: Tire faces outward toward tube wall
c3 = Constraint(
    name="tire_faces_outward",
    constraint_type=ConstraintType.RADIAL_FROM_CENTER,
    part="DDSM115_MOTOR",
    datum="tire_face",
    center=(0.0, 0.0, 0.0),
    direction="outward",
    tolerance_deg=5.0,
    description="Tire tread must face tube wall for traction"
)
constraints.append(c3)
print(f"✓ {c3.name}: {c3.constraint_type.value}")
print(f"    Validates: Tire faces {c3.direction} toward tube wall")


# Step 4: Simulate assembly and validation
print("\n4. Simulating assembly with constraint validation...")
print("-" * 80)

# In a real assembly, mates would compute transforms via solver
# For this demo, we manually create a "solved" motor position

# Motor positioned at wheel path with correct orientation
motor_world = Datum(
    name="motor_axis",
    datum_type=DatumType.AXIS,
    origin=[124.5, 0.0, 50.0, 1],  # At wheel path radius
    direction=[0.0, 1.0, 0.0, 0],   # Tangent direction
    description="Motor axis in world coordinates"
)

tire_world = Datum(
    name="tire_face",
    datum_type=DatumType.PLANE,
    origin=[124.5, 10.0, 50.0, 1],
    normal=[1.0, 0.0, 0.0, 0],      # Outward
    description="Tire face in world coordinates"
)

print("\nEvaluating constraints...")
all_passed = True

# Evaluate AT_RADIUS constraint
result = c1.evaluate(motor_world)
print(f"\n  {c1.name}:")
print(f"    {result}")
if not result.passed:
    all_passed = False

# Evaluate TANGENT_TO_CIRCLE constraint
result = c2.evaluate(motor_world)
print(f"\n  {c2.name}:")
print(f"    {result}")
if not result.passed:
    all_passed = False

# Evaluate RADIAL_FROM_CENTER constraint
result = c3.evaluate(tire_world)
print(f"\n  {c3.name}:")
print(f"    {result}")
if not result.passed:
    all_passed = False

# Summary
print("\n" + "=" * 80)
if all_passed:
    print("✅ Assembly VALID - All constraints satisfied!")
else:
    print("❌ Assembly INVALID - Some constraints failed!")
print("=" * 80)


# Step 5: Show the value proposition
print("\n5. Value Proposition: Why This System Matters")
print("-" * 80)
print("""
Before constraint system:
  ❌ Orientation errors caught only after render/export
  ❌ Design intent documented only in comments/markdown
  ❌ No automatic validation of assembly correctness
  ❌ Easy to make radial vs. tangent mistakes

After constraint system:
  ✅ Orientation errors caught at definition time
  ✅ Design intent is executable code
  ✅ Automatic validation with clear error messages
  ✅ Constraints enforce correct orientations

Example error detection:
  "Tangent angle error: 90.00° (tolerance: 1.0°)"
  → Immediately tells you axis is RADIAL instead of TANGENT

This prevents hours of debugging incorrect assemblies!
""")

print("=" * 80)
print("Workflow Complete: Datum → Mate → Constraint")
print("=" * 80)
