# yapCAD Assembly System

This package provides a constraint-based assembly system for yapCAD, enabling declarative part positioning through datum features, mates, and design constraints.

## Overview

The assembly system is inspired by modern CAD tools (SolidWorks mates, Fusion 360 joints) but tailored for programmatic 3D generative design workflows. It prevents common orientation errors by making geometric constraints first-class objects rather than documentation.

## Key Concepts

### Datum Features

**Datums** are named geometric references defined on parts in their local coordinate systems. They provide explicit geometric features that can be used to establish mates and constraints.

**Datum Types:**
- `POINT`: A single point in 3D space
- `AXIS`: An infinite line (origin + direction)
- `PLANE`: An infinite plane (origin + normal)
- `FRAME`: A full coordinate frame (origin + x_axis + y_axis)
- `CIRCLE`: A circle (origin + normal + radius)

### Part Definitions

A `PartDefinition` is a part with named datum features and metadata (geometry source, material, printability).

### Assembly System (To Be Implemented)

The complete system will include:
- **Mates**: Position and orient parts relative to each other (FLUSH, CONCENTRIC, etc.)
- **Constraints**: High-level design rules that validate intent
- **Solver**: Computes transforms and validates all constraints

## Module Status

### Implemented (Current)

- ✅ `datum.py`: Datum features and part definitions
  - `Datum` dataclass with full validation
  - `DatumType` enum
  - `PartDefinition` with datum management
  - Transform support using yapCAD's Matrix class

### To Be Implemented

- ⏳ `mate.py`: Assembly mates (FLUSH, CONCENTRIC, PARALLEL, etc.)
- ⏳ `constraint.py`: Design constraints (TANGENT_TO_CIRCLE, POINTS_TOWARD, etc.)
- ⏳ `assembly.py`: Assembly container and result classes
- ⏳ `solver.py`: Constraint solver and validation engine

## Usage Examples

### Basic Datum Creation

```python
from yapcad.assembly.datum import Datum, DatumType, PartDefinition
from yapcad.geom import point

# Create a point datum
center = Datum(
    name="center",
    datum_type=DatumType.POINT,
    origin=point(0, 0, 0),
    description="Part center point"
)

# Create an axis datum (for rotation or bore alignment)
axis = Datum(
    name="motor_axis",
    datum_type=DatumType.AXIS,
    origin=point(0, 0, 0),
    direction=[0, 1, 0, 0],  # Y-axis direction (w=0)
    description="Motor rotation axis"
)

# Create a plane datum (for mounting surfaces)
mount_face = Datum(
    name="stator_face",
    datum_type=DatumType.PLANE,
    origin=point(0, 10, 0),
    normal=[0, 1, 0, 0],  # Normal pointing +Y (w=0)
    description="Stator mounting surface"
)

# Create a circle datum (for bolt hole patterns)
bolt_pattern = Datum(
    name="mounting_holes",
    datum_type=DatumType.CIRCLE,
    origin=point(0, 10, 0),
    normal=[0, 1, 0, 0],
    radius=15.2,
    description="M2.5 mounting holes at 15.2mm radius"
)
```

### Part Definition with Datums

```python
# Create a motor part definition
motor = PartDefinition(
    name="DDSM115_MOTOR",
    geometry_source="cots/motor.dsl",
    is_printable=False,
    material="Aluminum",
    description="DDSM115 servo motor"
)

# Add datum features
motor.add_datum(axis)
motor.add_datum(mount_face)
motor.add_datum(bolt_pattern)

# Retrieve a datum by name
shaft_axis = motor.get_datum("motor_axis")

# Validate all datums
issues = motor.validate_datums()
if issues:
    print("Validation warnings:")
    for issue in issues:
        print(f"  - {issue}")
```

### Datum Transformation

Datums transform with parts when assembled:

```python
from yapcad.xform import Translation, Rotation

# Original datum at origin
datum = Datum(
    "axis",
    DatumType.AXIS,
    origin=point(0, 0, 0),
    direction=[1, 0, 0, 0]
)

# Translate by [10, 20, 30]
T = Translation(point(10, 20, 30))
translated_datum = datum.transform(T)
# translated_datum.origin == [10, 20, 30, 1]
# translated_datum.direction == [1, 0, 0, 0] (unchanged)

# Rotate 90° around Z-axis
R = Rotation([0, 0, 1, 0], 90)
rotated_datum = datum.transform(R)
# rotated_datum.direction rotates from X-axis to Y-axis
```

## Design Principles

1. **Constraints are First-Class Objects**: Not comments, but executable specifications
2. **Fail-Fast Validation**: Invalid assemblies detected at definition time
3. **Separation of Concerns**: Mates position parts, constraints verify intent
4. **Composability**: Complex assemblies built from validated subassemblies
5. **Integration**: Augments existing kinematic chain transform system

## Coordinate System

All geometric data uses yapCAD's homogeneous coordinate convention:
- **Points**: `[x, y, z, 1]` where w=1 indicates a position
- **Direction vectors**: `[x, y, z, 0]` where w=0 indicates a direction (no translation)

Direction vectors are automatically normalized during transformation to maintain unit length.

## Validation

The module provides extensive validation:

### Datum Validation (at creation)
- Required attributes present for each datum type
- Vectors have proper 4-component format
- Direction/normal vectors are non-zero
- Radius is positive (for circles)

### Part Validation
```python
issues = part.validate_datums()
```
Checks for:
- Duplicate datum definitions (same origin and type)
- Overlapping features that might indicate copy-paste errors
- Parallel/coplanar datums at same location

## Error Handling

The module raises appropriate exceptions:
- `ValueError`: Invalid datum parameters or duplicate names
- `KeyError`: Datum not found in part definition

Example error messages include context:
```
ValueError: Axis datum 'motor_shaft' requires direction vector
ValueError: Circle datum 'bolt_pattern' requires normal and radius
KeyError: Datum 'nonexistent' not found on part 'MOTOR'. Available datums: axis, mount_face
```

## Testing

Test scripts are provided in the repository root:
- `test_datum_transform.py`: Tests transformation functionality
- `test_datum_validation.py`: Tests validation and error handling

Run tests:
```bash
python3 test_datum_transform.py
python3 test_datum_validation.py
```

## See Also

- **Usage Patterns**: See `geom.py` and `xform.py` for coordinate conventions
- **Related Modules**: `yapcad.geom`, `yapcad.xform`, `yapcad.geom3d`

## Future Development

Once the complete assembly system is implemented, typical workflow will be:

```python
# Define parts with datums
motor = PartDefinition("motor")
motor.add_datum(Datum(...))

bracket = PartDefinition("bracket")
bracket.add_datum(Datum(...))

# Create assembly
assembly = Assembly("wheel_drive")
assembly.add_part(motor)
assembly.add_part(bracket)

# Add mates to position parts
assembly.add_mate(Mate(
    MateType.FLUSH,
    "motor", "stator_face",
    "bracket", "motor_mount"
))

# Add constraints to validate design intent
assembly.add_constraint(Constraint(
    ConstraintType.TANGENT_TO_CIRCLE,
    "motor", "axis",
    "chassis", "wheel_orbit",
    description="Motor axis tangent for proper wheel rolling"
))

# Solve and validate
result = assembly.solve()
if result.is_valid:
    # Use transforms for rendering
    transforms = result.transforms
else:
    # Show errors
    print(result.errors)
```
