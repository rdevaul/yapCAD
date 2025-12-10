# yapCAD DSL Reference

A domain-specific language for parametric CAD design with full type safety and provenance tracking.

## Quick Start

```python
module my_design

command MAKE_BOX(width: float, height: float, depth: float) -> solid:
    result: solid = box(width, height, depth)
    emit result
```

## Syntax Overview

yapCAD DSL uses Python-like indentation-based syntax:
- No semicolons required (optional for legacy compatibility)
- Indentation defines blocks
- `#` for comments
- Strong static typing with type inference

## Module Structure

```python
module module_name

# Optional imports
use other_module
use package.submodule as alias

# Command definitions
command COMMAND_NAME(param: type, param2: type = default) -> return_type:
    # body
    emit result
```

## Types

### Primitive Types (Tier 0)

| Type | Description | Examples |
|------|-------------|----------|
| `int` | Integer number | `42`, `-10`, `0` |
| `float` | Floating-point number | `3.14`, `-1.5`, `0.0` |
| `bool` | Boolean | `true`, `false` |
| `string` | Text string | `"hello"`, `"gear_1"` |

### Geometric Primitives (Tier 1)

| Type | Description | Constructor |
|------|-------------|-------------|
| `point` | 2D or 3D point | `point(x, y)` or `point(x, y, z)` |
| `point2d` | 2D point | `point(x, y)` |
| `point3d` | 3D point | `point(x, y, z)` |
| `vector` | 2D or 3D direction | `vector(dx, dy)` or `vector(dx, dy, dz)` |
| `vector2d` | 2D vector | `vector(dx, dy)` |
| `vector3d` | 3D vector | `vector(dx, dy, dz)` |
| `transform` | Transformation matrix | `translate_xform()`, `rotate_xform()`, etc. |

### Curves (Tier 2)

| Type | Description | Constructor |
|------|-------------|-------------|
| `line_segment` | Straight line | `line(start, end)` |
| `arc` | Circular arc | `arc(center, radius, start_angle, end_angle)` |
| `circle` | Full circle | `circle(center, radius)` |
| `ellipse` | Ellipse | `ellipse(center, major, minor, rotation?)` |
| `bezier` | Bezier curve | `bezier(control_points)` |
| `catmullrom` | Catmull-Rom spline | `catmullrom(control_points, tension?)` |
| `nurbs` | NURBS curve | `nurbs(control_points, weights, knots, degree)` |

**Curve Methods:**
- `.at(t: float) -> point` - Point at parameter t (0-1)
- `.tangent_at(t: float) -> vector` - Tangent vector at t
- `.normal_at(t: float) -> vector` - Normal vector at t
- `.curvature_at(t: float) -> float` - Curvature at t
- `.length() -> float` - Total arc length

### Compound Curves (Tier 3)

| Type | Description | Constructor |
|------|-------------|-------------|
| `path2d` | 2D path of curves | `path(segments)` |
| `path3d` | 3D path of curves | `path(segments)` |
| `profile2d` | Open 2D profile | `path(segments)` |
| `region2d` | Closed 2D region | `rectangle()`, `regular_polygon()`, `close()` |
| `loop3d` | Closed 3D loop | `close(path3d)` |

**Region2D Methods:**
- `.union(other: region2d) -> region2d` - Boolean union
- `.difference(other: region2d) -> region2d` - Boolean subtraction
- `.intersection(other: region2d) -> region2d` - Boolean intersection

### Surfaces (Tier 4)

| Type | Description | Constructor |
|------|-------------|-------------|
| `surface` | 3D surface | `planar_surface()`, `cylindrical_surface()`, `loft_surface()` |
| `shell` | Collection of surfaces | `shell(surfaces)` |

### Solids (Tier 5)

| Type | Description | Constructor |
|------|-------------|-------------|
| `solid` | 3D solid volume | `box()`, `cylinder()`, `sphere()`, `extrude()`, etc. |

**Solid Methods:**
- `.union(other: solid) -> solid` - Boolean union
- `.difference(other: solid) -> solid` - Boolean subtraction
- `.intersection(other: solid) -> solid` - Boolean intersection
- `.translate(v: vector) -> solid` - Move by vector
- `.rotate(axis: vector3d, angle: float) -> solid` - Rotate around axis
- `.scale(factors: vector) -> solid` - Scale by factors
- `.apply(t: transform) -> solid` - Apply transformation

### Generic Types

| Type | Description | Example |
|------|-------------|---------|
| `list<T>` | List of elements | `list<float>`, `list<point>` |
| `T?` | Optional type | `float?` (nullable float) |

## Built-in Functions

### Solid Constructors

```python
# Box: width (x), depth (y), height (z)
box(width: float, depth: float, height: float) -> solid

# Cylinder centered at origin
cylinder(radius: float, height: float) -> solid

# Sphere centered at origin
sphere(radius: float) -> solid

# Cone/frustum
cone(radius1: float, radius2: float, height: float) -> solid

# Involute spur gear
involute_gear(teeth: int, module_mm: float, pressure_angle: float, face_width: float) -> solid
```

### Solid from 2D

```python
# Extrude a 2D region
extrude(profile: region2d, height: float, direction: vector3d) -> solid

# Revolve around an axis
revolve(profile: region2d, axis: vector3d, angle: float) -> solid

# Sweep along a path
sweep(profile: region2d, spine: path3d) -> solid

# Loft between profiles
loft(profiles: list<region2d>) -> solid
```

### Boolean Operations

```python
# Combine solids
union(a: solid, b: solid) -> solid

# Subtract b from a
difference(a: solid, b: solid) -> solid

# Keep overlapping volume
intersection(a: solid, b: solid) -> solid
```

### Transformations

```python
# Transform solids directly
translate(s: solid, x: float, y: float, z: float) -> solid
rotate(s: solid, x: float, y: float, z: float) -> solid  # angles in degrees
scale(s: solid, x: float, y: float, z: float) -> solid

# Create transform matrices
translate_xform(v: vector) -> transform
rotate_xform(axis: vector3d, angle: float) -> transform
scale_xform(factors: vector) -> transform
scale_uniform(factor: float) -> transform
mirror(plane_normal: vector3d) -> transform
identity_transform() -> transform

# Apply transforms
apply(t: transform, shape: solid) -> solid
```

### 2D Shapes

```python
# Rectangle centered at a point
rectangle(width: float, height: float, center: point2d) -> region2d

# Regular polygon
regular_polygon(n: int, radius: float, center: point2d) -> region2d

# Circle
circle(center: point, radius: float) -> circle
```

### Pattern Operations

```python
# Radial/circular pattern
radial_pattern(shape, count: int, axis: vector3d, center: point)

# Linear pattern
linear_pattern(shape, count: int, spacing: vector)
```

### Query Operations

```python
volume(s: solid) -> float
surface_area(s: solid) -> float
area(r: region2d) -> float
perimeter(r: region2d) -> float
centroid(s: solid) -> point3d
distance(a: point, b: point, tolerance: float) -> float
is_empty(s: solid) -> bool
```

### Math Functions

```python
# Trigonometry (radians)
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
min(a: float, b: float) -> float
max(a: float, b: float) -> float

# Angle conversion
radians(degrees: float) -> float
degrees(radians: float) -> float
```

### List Operations

```python
len(list: list<T>) -> int
range(end: int) -> list<int>              # [0, 1, ..., end-1]
range(start: int, end: int) -> list<int>  # [start, ..., end-1]
```

### Utility

```python
empty_solid() -> solid      # Empty solid for accumulating unions
empty_region() -> region2d  # Empty region
print(value) -> bool        # Debug printing
```

## Statements

### Variable Declaration

```python
# With explicit type
width: float = 10.0
gear: solid = involute_gear(24, 2.0, 20.0, 10.0)

# Legacy 'let' syntax (still supported)
let height: float = 5.0
```

### Require (Assertions)

```python
# Validate parameters
require width > 0.0, "width must be positive"
require height > 0.0 and depth > 0.0, "dimensions must be positive"
```

### Emit (Output)

```python
# Emit the result of the command
emit result

# Emit with metadata
emit result, name="gear_24t", material="steel"
```

### Control Flow

```python
# For loop with list comprehension
squares: list<int> = [x * x for x in range(10)]

# Note: Full if/while statements work in pure indentation mode
```

## Common Patterns

### Box with Hole

```python
command BRACKET(width: float, height: float, thickness: float, hole_radius: float) -> solid:
    # Create main plate
    plate: solid = box(width, height, thickness)

    # Create hole
    hole: solid = cylinder(hole_radius, thickness + 1.0)
    hole_pos: solid = translate(hole, width/2.0, height/2.0, -0.5)

    # Subtract hole from plate
    result: solid = difference(plate, hole_pos)
    emit result
```

### Gear Creation

```python
command MAKE_GEAR(teeth: int, module_mm: float) -> solid:
    # Standard 20-degree pressure angle, 10mm face width
    gear: solid = involute_gear(teeth, module_mm, 20.0, 10.0)
    emit gear
```

### Parametric Design

```python
command ENCLOSURE(
    width: float,
    depth: float,
    height: float,
    wall_thickness: float
) -> solid:
    require wall_thickness < width / 2.0, "walls too thick"
    require wall_thickness < depth / 2.0, "walls too thick"

    outer: solid = box(width, depth, height)
    inner: solid = box(
        width - 2.0 * wall_thickness,
        depth - 2.0 * wall_thickness,
        height - wall_thickness
    )
    inner_positioned: solid = translate(inner, wall_thickness, wall_thickness, wall_thickness)

    shell: solid = difference(outer, inner_positioned)
    emit shell
```

### Transform Chain

```python
command POSITIONED_PART(angle: float, offset_x: float) -> solid:
    base: solid = box(10.0, 10.0, 10.0)
    rotated: solid = rotate(base, 0.0, 0.0, angle)
    moved: solid = translate(rotated, offset_x, 0.0, 0.0)
    emit moved
```

## Programmatic API

For agentic tools and programmatic access:

```python
from yapcad.dsl import (
    get_api_reference,      # Get complete API as dict
    get_function_info,      # Get info about a function
    list_functions,         # List all functions
    list_types,             # List all types
    describe_function,      # Human-readable description
    get_common_pattern,     # Get example code patterns
)

# Get machine-readable API reference
api = get_api_reference()
print(api["functions"]["box"]["signature"])
# -> "box(width: float, depth: float, height: float) -> solid"

# Get function details
info = get_function_info("cylinder")
print(info["description"])
# -> "Create a cylinder solid from radius and height"

# List solid-creating functions
solid_funcs = list_functions(category="solid")
# -> ['apply', 'box', 'cone', 'cylinder', ...]
```

## Type System Notes

1. **Int/Float compatibility**: `int` is assignable to `float` parameters
2. **Point polymorphism**: `point2d` and `point3d` are subtypes of `point`
3. **Vector polymorphism**: `vector2d` and `vector3d` are subtypes of `vector`
4. **List covariance**: `list<int>` is assignable to `list<float>` (elements coerced)

## Module System and Namespaces

### Import Syntax

```python
module my_design

# Import entire module
use other_module

# Import with alias
use package.submodule as alias

# Import from standard library
use yapcad.stdlib.gears
use yapcad.stdlib.fasteners as fasteners
```

### Namespace Conventions

yapCAD DSL organizes functions into three tiers:

1. **Core Builtins** (no import required)
   - Primitives: `point()`, `vector()`, `box()`, `cylinder()`, `sphere()`
   - Operations: `union()`, `difference()`, `intersection()`
   - Transforms: `translate()`, `rotate()`, `scale()`
   - Math: `sin()`, `cos()`, `sqrt()`, etc.
   - Queries: `volume()`, `area()`, `perimeter()`

2. **Standard Library** (`yapcad.stdlib.*`)
   - `yapcad.stdlib.gears` - Parametric gear generation
     - `involute_spur(teeth, module_mm, pressure_angle, face_width)`
     - `involute_helical(teeth, module_mm, helix_angle, face_width)`
     - `bevel_gear(teeth, module_mm, pitch_angle)`
   - `yapcad.stdlib.fasteners` - Standard fasteners and holes
     - `hex_bolt(standard, size, length)` - ISO/ANSI hex bolts
     - `hex_nut(standard, size)` - Hex nuts
     - `washer(standard, size)` - Flat washers
     - `threaded_hole(standard, size, depth, tapped?)` - Threaded holes
   - `yapcad.stdlib.profiles` - Common 2D profiles
     - `channel(width, height, flange, web)`
     - `angle(leg1, leg2, thickness)`
     - `i_beam(width, height, flange, web)`

3. **User Packages** (project-specific)
   - Custom DSL modules in `.ycpkg` packages
   - Imported via `use package_name.module`

### Example: Using Standard Library

```python
module assembly

use yapcad.stdlib.gears
use yapcad.stdlib.fasteners as f

command GEAR_ASSEMBLY(teeth: int, bolt_count: int) -> solid:
    # Create gear using stdlib
    gear: solid = gears.involute_spur(teeth, 2.0, 20.0, 10.0)

    # Add mounting holes using fasteners stdlib
    bolt_circle_radius: float = (teeth * 2.0) / 4.0

    # Create bolt holes
    holes: list<solid> = []
    for i in range(bolt_count):
        angle: float = 360.0 * i / bolt_count
        x: float = bolt_circle_radius * cos(radians(angle))
        y: float = bolt_circle_radius * sin(radians(angle))
        hole: solid = f.threaded_hole("ISO", "M6", 15.0, true)
        positioned: solid = translate(hole, x, y, 0.0)
        holes = concat(holes, [positioned])

    # Subtract all holes
    result: solid = difference(gear, union(holes))
    emit result
```

### Future: Package Registry

A package registry will allow sharing and discovering DSL modules:

```python
# Future syntax for external packages
use external:gears-advanced.planetary
use external:aerospace-fasteners.an_series
```

## Error Handling

The DSL provides detailed error messages with source locations:

```
error[E201]: Type mismatch: expected float, got string
  --> design.dsl:5:12
   |
 5 |     let x: float = "hello";
   |            ^^^^^
```

Common error codes:
- `E101`: Parser errors (syntax)
- `E201`: Type errors
- `E301`: Undefined variables/functions
- `E302`: Duplicate definitions
