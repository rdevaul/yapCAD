# yapCAD DSL Reference

A domain-specific language for parametric CAD design with full type safety and provenance tracking.

## Quick Start

```python
module my_design

command MAKE_BOX(width: float, height: float, depth: float) -> solid:
    result: solid = box(width, height, depth)
    emit result
```

## CLI Usage

```bash
# Check DSL syntax and types
python -m yapcad.dsl check myfile.dsl

# List available commands
python -m yapcad.dsl list myfile.dsl

# Run a command
python -m yapcad.dsl run myfile.dsl COMMAND_NAME

# Run with parameters
python -m yapcad.dsl run myfile.dsl COMMAND_NAME --param width=10.0 --param height=5.0

# Export to STEP file
python -m yapcad.dsl run myfile.dsl COMMAND_NAME --output result.step

# Export to package
python -m yapcad.dsl run myfile.dsl COMMAND_NAME --package output.ycpkg
```

## Syntax Overview

yapCAD DSL uses Python-like syntax:
- Colons (`:`) after command signatures and control flow
- Indentation defines blocks
- `#` for comments
- Strong static typing with explicit type annotations
- Both `let name: type = value` and `name: type = value` for variable declaration

## Module Structure

```python
module module_name

# Optional imports (future feature)
use other_module
use package.submodule as alias

# Helper command (lowercase name, not exported)
command make_helper(x: float) -> solid:
    # ...
    emit result

# Exported command (UPPERCASE name, appears in CLI)
command MAKE_PART(param: float, param2: float = 10.0) -> solid:
    # body
    emit result
```

**Note:** Commands with UPPERCASE names are exported and visible to `dsl list`. Commands with lowercase names are helpers usable within the module but not directly callable from CLI.

## Types

### Primitive Types

| Type | Description | Examples |
|------|-------------|----------|
| `int` | Integer number | `42`, `-10`, `0` |
| `float` | Floating-point number | `3.14`, `-1.5`, `0.0` |
| `bool` | Boolean | `true`, `false` |
| `string` | Text string | `"hello"`, `"gear_1"` |

### Geometric Types

| Type | Description | Constructor |
|------|-------------|-------------|
| `point` | 2D or 3D point | `point(x, y)` or `point(x, y, z)` |
| `point2d` | 2D point | `point2d(x, y)` |
| `point3d` | 3D point | `point(x, y, z)` |
| `vector` | 2D or 3D direction | `vector(dx, dy)` or `vector(dx, dy, dz)` |
| `vector2d` | 2D vector | `vector2d(dx, dy)` |
| `vector3d` | 3D vector | `vector(dx, dy, dz)` |
| `transform` | Transformation matrix | `translate_xform()`, `rotate_xform()`, etc. |

### Curve Types

| Type | Description | Constructor |
|------|-------------|-------------|
| `line_segment` | Straight line | `line(start, end)` |
| `arc` | Circular arc | `arc(center, radius, start_angle, end_angle)` |
| `circle` | Full circle | `circle(center, radius)` |

### Compound Types

| Type | Description | Constructor |
|------|-------------|-------------|
| `path3d` | 3D path of segments | `make_path3d(segments)`, `path3d_line()`, `path3d_arc()` |
| `region2d` | Closed 2D region | `rectangle()`, `regular_polygon()` |
| `solid` | 3D solid volume | `box()`, `cylinder()`, `sphere()`, etc. |

### Generic Types

| Type | Description | Example |
|------|-------------|---------|
| `list<T>` | List of elements | `list<float>`, `list<point>`, `list<solid>` |

## Built-in Functions

### Math Functions

```python
# Trigonometry (argument in radians)
sin(x: float) -> float
cos(x: float) -> float
tan(x: float) -> float
asin(x: float) -> float
acos(x: float) -> float
atan(x: float) -> float
atan2(y: float, x: float) -> float

# General math
sqrt(x: float) -> float
abs(x: float) -> float
pow(base: float, exp: float) -> float
floor(x: float) -> int
ceil(x: float) -> int
round(x: float) -> int
min(a: float, b: float, ...) -> float  # variadic
max(a: float, b: float, ...) -> float  # variadic

# Angle conversion
radians(degrees: float) -> float
degrees(radians: float) -> float

# Constants
pi() -> float   # 3.14159...
tau() -> float  # 2 * pi
```

### Point and Vector Constructors

```python
# Points
point(x: float, y: float) -> point2d
point(x: float, y: float, z: float) -> point3d
point2d(x: float, y: float) -> point2d

# Vectors
vector(dx: float, dy: float) -> vector2d
vector(dx: float, dy: float, dz: float) -> vector3d
vector2d(dx: float, dy: float) -> vector2d
```

### 2D Shape Constructors

```python
# Rectangle centered at a point (or origin if center omitted)
rectangle(width: float, height: float) -> region2d
rectangle(width: float, height: float, center: point2d) -> region2d

# Regular polygon
regular_polygon(sides: int, radius: float) -> region2d
regular_polygon(sides: int, radius: float, center: point2d) -> region2d
```

### Curve Constructors

```python
# Line segment between two points
line(start: point, end: point) -> line_segment

# Arc from center point
arc(center: point, radius: float, start_angle: float, end_angle: float) -> arc

# Full circle
circle(center: point, radius: float) -> circle
```

### Path3D Constructors (for sweep operations)

```python
# Create a path3d from segment list
make_path3d(segments: list<path3d>) -> path3d

# Line segment for path3d
path3d_line(start: point3d, end: point3d) -> path3d

# Arc segment for path3d
path3d_arc(center: point3d, start: point3d, end: point3d, normal: vector3d) -> path3d
```

### Solid Constructors

```python
# Box: width (X), depth (Y), height (Z) - centered at origin
box(width: float, depth: float, height: float) -> solid

# Cylinder: centered at origin, extends from Z=0 to Z=height
cylinder(radius: float, height: float) -> solid

# Sphere: centered at origin
sphere(radius: float) -> solid

# Cone/frustum: radius1 at bottom, radius2 at top
cone(radius1: float, radius2: float, height: float) -> solid

# Involute spur gear
involute_gear(teeth: int, module_mm: float, pressure_angle: float, face_width: float) -> solid
```

### Solid from 2D Operations

```python
# Extrude a 2D region along Z axis
extrude(profile: region2d, height: float) -> solid

# Revolve a 2D region around an axis
revolve(profile: region2d, axis: vector3d, angle: float) -> solid

# Sweep a 2D profile along a 3D path
sweep(profile: region2d, spine: path3d) -> solid

# Sweep with inner void (hollow tube)
sweep_hollow(outer_profile: region2d, inner_profile: region2d, spine: path3d) -> solid

# Adaptive sweep - profile rotates to track path tangent
# Creates ruled surfaces for straight segments
sweep_adaptive(profile: region2d, spine: path3d, threshold_deg: float) -> solid

# Adaptive sweep with hollow profile
sweep_adaptive_hollow(
    outer_profile: region2d,
    inner_profile: region2d,
    spine: path3d,
    threshold_deg: float
) -> solid
```

### Boolean Operations

```python
# Union (combine) solids
union(a: solid, b: solid) -> solid
union(solids: list<solid>) -> solid  # variadic

# Difference (subtract b from a)
difference(a: solid, b: solid) -> solid
difference(a: solid, tools: list<solid>) -> solid  # variadic

# Intersection (keep overlapping volume)
intersection(a: solid, b: solid) -> solid
intersection(solids: list<solid>) -> solid  # variadic
```

### Transformation Functions

```python
# Transform solids directly (angles in degrees)
translate(s: solid, x: float, y: float, z: float) -> solid
rotate(s: solid, rx: float, ry: float, rz: float) -> solid  # Euler angles
scale(s: solid, sx: float, sy: float, sz: float) -> solid

# Create transform matrices (for advanced use)
translate_xform(v: vector) -> transform
rotate_xform(axis: vector3d, angle: float) -> transform
scale_xform(factors: vector) -> transform
identity_transform() -> transform

# Apply transform to geometry
apply(t: transform, shape: solid) -> solid
apply_surface(t: transform, surf: surface) -> surface
apply_point(t: transform, p: point) -> point
apply_vector(t: transform, v: vector3d) -> vector3d
```

### Query Functions

```python
# Solid queries
volume(s: solid) -> float
surface_area(s: solid) -> float

# Region2D queries
area(r: region2d) -> float
perimeter(r: region2d) -> float
```

### List Functions

```python
# List operations
len(lst: list<T>) -> int
range(end: int) -> list<int>                       # [0, 1, ..., end-1]
range(start: int, end: int) -> list<int>           # [start, ..., end-1]
range(start: int, end: int, step: int) -> list<int>
concat(list1: list<T>, list2: list<T>) -> list<T>
reverse(lst: list<T>) -> list<T>
flatten(nested: list<list<T>>) -> list<T>
```

### Utility Functions

```python
# Debug output
print(value, ...) -> bool  # variadic, returns true
```

## Statements

### Variable Declaration

```python
# Explicit type annotation (preferred)
width: float = 10.0
gear: solid = involute_gear(24, 2.0, 20.0, 10.0)

# With 'let' keyword (also supported)
let height: float = 5.0
```

### Require (Assertions)

```python
# Validate parameters - raises error if false
require width > 0.0
require height > 0.0 and depth > 0.0
```

### Emit (Output)

```python
# Return the result from a command
emit result
```

## Common Patterns

### Box with Hole

```python
module bracket

command MAKE_BRACKET(
    width: float,
    height: float,
    thickness: float,
    hole_radius: float
) -> solid:
    # Create main plate
    plate: solid = box(width, height, thickness)

    # Create hole cylinder (slightly longer for clean cut)
    hole: solid = cylinder(hole_radius, thickness + 1.0)

    # Position hole at center of plate (adjust Z to cut through)
    hole_pos: solid = translate(hole, 0.0, 0.0, -0.5)

    # Subtract hole from plate
    result: solid = difference(plate, hole_pos)
    emit result
```

### Hollow Box (Shell)

```python
module enclosure

command MAKE_ENCLOSURE(
    width: float,
    depth: float,
    height: float,
    wall_thickness: float
) -> solid:
    require wall_thickness < width / 2.0
    require wall_thickness < depth / 2.0

    outer: solid = box(width, depth, height)

    inner_w: float = width - 2.0 * wall_thickness
    inner_d: float = depth - 2.0 * wall_thickness
    inner_h: float = height - wall_thickness

    inner: solid = box(inner_w, inner_d, inner_h)
    inner_positioned: solid = translate(inner, 0.0, 0.0, wall_thickness)

    shell: solid = difference(outer, inner_positioned)
    emit shell
```

### Swept Path with Bends

```python
module pipe_jig

# Helper: create a path with two bends
command make_bent_path(
    seg1_length: float,
    seg2_length: float,
    seg3_length: float,
    bend_angle: float
) -> path3d:
    angle_rad: float = radians(bend_angle)
    double_rad: float = radians(2.0 * bend_angle)

    # First segment along Y
    p0: point3d = point(0.0, 0.0, 0.0)
    p1: point3d = point(0.0, seg1_length, 0.0)

    # Direction after first bend
    dir1_x: float = sin(angle_rad)
    dir1_y: float = cos(angle_rad)

    # Second bend point
    p2_x: float = seg2_length * dir1_x
    p2_y: float = seg1_length + seg2_length * dir1_y
    p2: point3d = point(p2_x, p2_y, 0.0)

    # Direction after second bend
    dir2_x: float = sin(double_rad)
    dir2_y: float = cos(double_rad)

    # End point
    p3_x: float = p2_x + seg3_length * dir2_x
    p3_y: float = p2_y + seg3_length * dir2_y
    p3: point3d = point(p3_x, p3_y, 0.0)

    # Build path segments
    seg1: path3d = path3d_line(p0, p1)
    seg2: path3d = path3d_line(p1, p2)
    seg3: path3d = path3d_line(p2, p3)

    segments: list<path3d> = [seg1, seg2, seg3]
    spine: path3d = make_path3d(segments)
    emit spine

# Main command: sweep a profile along the bent path
command MAKE_BENT_TUBE(
    profile_width: float = 10.0,
    profile_height: float = 10.0,
    bend_angle: float = 15.0
) -> solid:
    # Create profile
    profile: region2d = rectangle(profile_width, profile_height)

    # Create path
    spine: path3d = make_bent_path(50.0, 200.0, 50.0, bend_angle)

    # Sweep with adaptive tangent tracking (5 degree threshold)
    result: solid = sweep_adaptive(profile, spine, 5.0)
    emit result
```

### Gear Creation

```python
module gears

command MAKE_SPUR_GEAR(
    teeth: int,
    module_mm: float,
    face_width: float
) -> solid:
    require teeth >= 6
    require module_mm > 0.0

    # Standard 20-degree pressure angle
    gear: solid = involute_gear(teeth, module_mm, 20.0, face_width)
    emit gear
```

### Combining Multiple Parts

```python
module assembly

command MAKE_ASSEMBLY(count: int) -> solid:
    # Create base
    base: solid = box(100.0, 100.0, 10.0)

    # Create a cylinder to add multiple times
    peg: solid = cylinder(5.0, 20.0)

    # Position and combine
    result: solid = base
    spacing: float = 20.0

    # Note: for loops with range work
    i: int = 0
    # Manual loop unrolling for now
    peg1: solid = translate(peg, -30.0, -30.0, 10.0)
    result = union(result, peg1)

    peg2: solid = translate(peg, 30.0, -30.0, 10.0)
    result = union(result, peg2)

    peg3: solid = translate(peg, -30.0, 30.0, 10.0)
    result = union(result, peg3)

    peg4: solid = translate(peg, 30.0, 30.0, 10.0)
    result = union(result, peg4)

    emit result
```

## Type System Notes

1. **Int/Float compatibility**: `int` is assignable to `float` parameters
2. **Point polymorphism**: `point2d` and `point3d` are subtypes of `point`
3. **Vector polymorphism**: `vector2d` and `vector3d` are subtypes of `vector`
4. **Explicit typing required**: All variables must have explicit type annotations

## Error Messages

The DSL provides detailed error messages with source locations:

```
error[E201]: Type mismatch: expected float, got string
  --> design.dsl:5:12
   |
 5 |     x: float = "hello"
   |            ^^^^^
```

Common error codes:
- `E101`: Parser errors (syntax)
- `E201`: Type errors
- `E301`: Undefined variables/functions
- `E302`: Duplicate definitions

## Package Integration

When using `--package`, the DSL automatically:
- Creates a `.ycpkg` directory structure
- Exports geometry to JSON format
- Generates STEP export
- Records provenance metadata (DSL command, parameters, version)

```bash
# Create a package with full provenance
python -m yapcad.dsl run design.dsl MAKE_PART --package output.ycpkg

# View the package
python tools/ycpkg_viewer.py output.ycpkg
```

## Programmatic API

For scripts and automation:

```python
from yapcad.dsl import compile_and_run

# Run a DSL command programmatically
source = open("design.dsl").read()
result = compile_and_run(source, "MAKE_PART", {"width": 10.0, "height": 5.0})

if result.success:
    solid = result.geometry
    # Use the solid...
else:
    print(f"Error: {result.error_message}")
```

## See Also

- `docs/yapCADone.rst` - yapCAD 1.0 roadmap and feature status
- `docs/ycpkg_spec.rst` - Package specification
- `examples/` - Example DSL files
