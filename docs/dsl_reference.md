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
| `ellipse` | Ellipse or elliptical arc | `ellipse(center, semi_major, semi_minor, ...)` |
| `bezier` | Bezier curve | `bezier(control_points)` |
| `catmullrom` | Catmull-Rom spline | `catmullrom(points, closed?, alpha?)` |
| `nurbs` | NURBS curve | `nurbs(points, weights?, degree?)` |

### Compound Types

| Type | Description | Constructor |
|------|-------------|-------------|
| `path2d` | 2D path of segments | `make_path2d(curves)` |
| `path3d` | 3D path of segments | `make_path3d(segments)`, `path3d_line()`, `path3d_arc()` |
| `region2d` | Closed 2D region | `rectangle()`, `regular_polygon()`, `polygon()`, `disk()` |
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
exp(x: float) -> float          # e^x
log(x: float) -> float          # natural logarithm
log10(x: float) -> float        # base-10 logarithm
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

# Polygon from list of points
polygon(points: list<point>) -> region2d

# Disk (filled circle as polygon approximation)
disk(center: point, radius: float) -> region2d
disk(center: point, radius: float, segments: int) -> region2d  # default 64 segments
```

### Curve Constructors

```python
# Line segment between two points
line(start: point, end: point) -> line_segment

# Arc from center point (angles in degrees)
arc(center: point, radius: float, start_angle: float, end_angle: float) -> arc

# Full circle
circle(center: point, radius: float) -> circle

# Ellipse (angles in degrees, rotation in degrees)
ellipse(center: point, semi_major: float, semi_minor: float) -> ellipse
ellipse(center: point, semi_major: float, semi_minor: float,
        rotation: float, start: float, end: float) -> ellipse

# Bezier curve from control points
bezier(control_points: list<point>) -> bezier

# Catmull-Rom spline
catmullrom(points: list<point>) -> catmullrom
catmullrom(points: list<point>, closed: bool, alpha: float) -> catmullrom
# alpha: 0.0=uniform, 0.5=centripetal (default), 1.0=chordal

# NURBS curve
nurbs(points: list<point>) -> nurbs
nurbs(points: list<point>, weights: list<float>, degree: int) -> nurbs
# degree default: 3
```

### Curve Sampling Functions

```python
# Sample a point on a curve at parameter t in [0, 1]
sample_curve(curve, t: float) -> point

# Sample n points along a curve
sample_curve_n(curve, n: int) -> list<point>

# Get the length of a curve
curve_length(curve) -> float
```

### Path Constructors

#### 2D Paths

```python
# Create a path from a list of curves
make_path2d(curves: list<curve>) -> path2d

# Close an open path to create a region
close_path(path: path2d) -> region2d

# Convert a spline to a region (polygon approximation)
region_from_spline(spline, segments: int = 64) -> region2d
```

#### 3D Paths (for sweep operations)

```python
# Create a path3d from segments
make_path3d(segments...) -> path3d  # variadic

# Line segment for path3d
path3d_line(start: point3d, end: point3d) -> path3d

# Arc segment for path3d (explicit normal)
path3d_arc(center: point3d, start: point3d, end: point3d, normal: vector3d) -> path3d

# Arc segment with auto-computed normal
# flip=false: shorter arc, flip=true: longer arc (opposite direction)
path3d_arc_auto(center: point3d, start: point3d, end: point3d, flip: bool) -> path3d
```

### 2D Boolean Operations

```python
# Union of two 2D regions
union2d(a: region2d, b: region2d) -> region2d

# Difference (subtract b from a)
difference2d(a: region2d, b: region2d) -> region2d

# Intersection (keep overlapping area)
intersection2d(a: region2d, b: region2d) -> region2d

# Aggregation operations (for lists)
union2d_all(regions: list<region2d>) -> region2d        # Union all regions
difference2d_all(base: region2d, tools: list<region2d>) -> region2d  # Subtract all tools from base
intersection2d_all(regions: list<region2d>) -> region2d  # Intersect all regions
```

### Solid Constructors

```python
# Box: width (X), depth (Y), height (Z) - centered at origin
box(width: float, depth: float, height: float) -> solid

# Cylinder: centered at origin, extends from Z=0 to Z=height
cylinder(radius: float, height: float) -> solid

# Sphere: centered at origin
sphere(radius: float) -> solid

# Oblate spheroid (flattened sphere like Earth/Mars)
oblate_spheroid(equatorial_diameter: float, oblateness: float) -> solid
# oblateness: 0=sphere, typical values: Earth~0.00335, Mars~0.00648

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
# Uses minimal-twist frame (default) to avoid unwanted rotation
sweep_adaptive(profile: region2d, spine: path3d, threshold_deg: float) -> solid

# Adaptive sweep with hollow profile
sweep_adaptive_hollow(
    outer_profile: region2d,
    inner_profile: region2d,
    spine: path3d,
    threshold_deg: float
) -> solid

# Frenet frame variants - profile follows natural curve curvature
# Appropriate for paths like helices where you want natural twisting
sweep_adaptive_frenet(profile: region2d, spine: path3d, threshold_deg: float) -> solid

sweep_adaptive_hollow_frenet(
    outer_profile: region2d,
    inner_profile: region2d,
    spine: path3d,
    threshold_deg: float
) -> solid

# Loft between multiple profiles
loft(profiles: list<region2d>) -> solid
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

# Compound - combine without merging (for multi-body assemblies)
compound(a: solid, b: solid) -> solid
compound(solids: list<solid>) -> solid  # variadic

# Aggregation operations (for lists) - cleaner syntax than variadic
union_all(solids: list<solid>) -> solid        # Union all solids in list
difference_all(base: solid, tools: list<solid>) -> solid  # Subtract all tools from base
intersection_all(solids: list<solid>) -> solid  # Intersect all solids in list
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
rotate_2d(angle: float) -> transform
scale_xform(factors: vector) -> transform
scale_uniform(factor: float) -> transform
mirror(plane_normal: vector3d) -> transform
mirror_2d(axis: vector2d) -> transform
mirror_y() -> transform  # Convenience: mirror across Y axis
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
centroid(s: solid) -> point3d
is_empty(s: solid) -> bool

# Region2D queries
area(r: region2d) -> float
perimeter(r: region2d) -> float

# Distance between points
distance(a: point, b: point, tolerance: float) -> float
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

### Aggregation Functions

```python
# Numeric aggregation
sum(values: list<float>) -> float       # Sum all values
product(values: list<float>) -> float   # Multiply all values
min_of(values: list<float>) -> float    # Find minimum value
max_of(values: list<float>) -> float    # Find maximum value

# Boolean aggregation
any_true(values: list<bool>) -> bool    # True if any element is true
all_true(values: list<bool>) -> bool    # True if all elements are true
```

### Utility Functions

```python
# Debug output
print(value, ...) -> bool  # variadic, returns true

# Empty geometry constructors
empty_solid() -> solid
empty_region() -> region2d
```

## Method Syntax

Some types support method-style calls as an alternative to function calls:

### Solid Methods

```python
# These are equivalent:
result = translate(my_solid, 10.0, 0.0, 0.0)
result = my_solid.translate(vector(10.0, 0.0, 0.0))

# Available methods on solids:
solid.union(other: solid) -> solid
solid.difference(other: solid) -> solid
solid.intersection(other: solid) -> solid
solid.translate(v: vector) -> solid
solid.rotate(axis: vector3d, angle: float) -> solid
solid.scale(factors: vector) -> solid
solid.apply(t: transform) -> solid
```

### Region2D Methods

```python
region.union(other: region2d) -> region2d
region.difference(other: region2d) -> region2d
region.intersection(other: region2d) -> region2d
```

### Transform Methods

```python
transform.compose(other: transform) -> transform
transform.inverse() -> transform
transform.translation() -> vector3d
transform.is_rigid() -> bool
```

### Curve Methods

```python
curve.at(t: float) -> point
curve.tangent_at(t: float) -> vector
curve.normal_at(t: float) -> vector
curve.curvature_at(t: float) -> float
curve.length() -> float
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

### For Loops

```python
# Iterate over a range
for i in range(10):
    # body

# Iterate over a list
for item in my_list:
    # body
```

### Conditionals

```python
if condition:
    # body
elif other_condition:
    # body
else:
    # body
```

### Conditional Expressions (Ternary)

The DSL supports Python-style conditional expressions for inline value selection:

```python
# Basic syntax: value_if_true if condition else value_if_false
result: float = 10.0 if use_metric else 25.4 * value

# With comparison
status: string = "hot" if temperature > 100.0 else "cold"

# Nested (chained) conditionals
grade: string = "A" if score >= 90.0 else ("B" if score >= 80.0 else "C")

# Selecting geometry
shape: solid = box(10.0, 10.0, 10.0) if use_cube else cylinder(5.0, 10.0)
```

Conditional expressions are useful for:
- Selecting between two values based on a boolean parameter
- Unit conversion (metric vs imperial)
- Choosing geometry based on configuration
- Inline computation without separate if/else blocks

**Note:** Both branches must have compatible types, and the condition must be a boolean.

### List Comprehensions

Create lists using comprehension syntax:

```python
# Map: transform each element
squares: list<float> = [x * x for x in values]

# Generate from range
angles: list<float> = [i * 30.0 for i in range(12)]

# Filter: select elements matching condition
positives: list<float> = [x for x in values if x > 0.0]

# Combined map and filter
big_squares: list<float> = [x * x for x in values if x > 10.0]
```

#### Nested Comprehensions

Multiple `for` clauses create nested iterations (cartesian products):

```python
# Nested comprehension - generates all (x, y) combinations
sums: list<int> = [x + y for x in xs for y in ys]
# Equivalent to: for x in xs: for y in ys: append(x + y)

# With conditions on outer loop
filtered_outer: list<int> = [x + y for x in xs if x > 0 for y in ys]

# With conditions on inner loop
filtered_inner: list<int> = [x + y for x in xs for y in ys if y < 10]

# With conditions on both loops
filtered_both: list<int> = [x + y for x in xs if x > 0 for y in ys if y < 10]

# Triple nesting
products: list<int> = [x + y + z for x in xs for y in ys for z in zs]
```

**Example: Creating patterns with symmetric geometry**

```python
# Generate holes at 3 sectors (0°, 120°, 240°) with multiple angles per sector
sector_offsets: list<float> = [0.0, 120.0, 240.0]
base_angles: list<float> = [50.0, 60.0, 70.0]

all_holes: list<solid> = [
    make_hole(radius, thickness, base_angle + offset)
    for offset in sector_offsets
    for base_angle in base_angles
]
# Creates 9 holes: 3 sectors × 3 angles per sector
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

    spine: path3d = make_path3d(seg1, seg2, seg3)
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

### Pattern with For Loop

```python
module assembly

command MAKE_PEGBOARD(count: int) -> solid:
    require count > 0 and count <= 10

    # Create base
    base: solid = box(100.0, 100.0, 10.0)
    peg: solid = cylinder(5.0, 20.0)

    # Create pegs in a grid
    result: solid = base
    spacing: float = 80.0 / (count + 1)

    for i in range(count):
        for j in range(count):
            x: float = -40.0 + (i + 1) * spacing
            y: float = -40.0 + (j + 1) * spacing
            peg_pos: solid = translate(peg, x, y, 10.0)
            result = union(result, peg_pos)

    emit result
```

### Pattern with Aggregation Functions

Using `union_all` and `difference_all` for cleaner multi-body operations:

```python
module thrust_plate

command MAKE_PLATE_WITH_HOLES(
    diameter: float,
    thickness: float,
    hole_radius: float
) -> solid:
    # Create base plate
    plate: solid = cylinder(diameter / 2.0, thickness)

    # Create holes at 120° intervals using list comprehension
    hole_angles: list<float> = [i * 120.0 for i in range(3)]
    holes: list<solid> = [
        translate(
            cylinder(hole_radius, thickness + 2.0),
            (diameter / 3.0) * cos(radians(angle)),
            (diameter / 3.0) * sin(radians(angle)),
            -1.0
        )
        for angle in hole_angles
    ]

    # Subtract all holes at once using difference_all
    result: solid = difference_all(plate, holes)
    emit result
```

Using `union_all` to combine multiple parts:

```python
module assembly

command MAKE_FRAME() -> solid:
    # Create individual parts
    bottom_rail: solid = box(100.0, 10.0, 10.0)
    top_rail: solid = translate(box(100.0, 10.0, 10.0), 0.0, 0.0, 50.0)
    left_post: solid = translate(box(10.0, 10.0, 50.0), -45.0, 0.0, 25.0)
    right_post: solid = translate(box(10.0, 10.0, 50.0), 45.0, 0.0, 25.0)

    # Union all parts at once
    parts: list<solid> = [bottom_rail, top_rail, left_post, right_post]
    frame: solid = union_all(parts)
    emit frame
```

### Spline-Based Profile

```python
module organic

command MAKE_BLOB(scale: float = 10.0) -> solid:
    # Create organic shape using Catmull-Rom spline
    pts: list<point> = [
        point(1.0 * scale, 0.0),
        point(0.8 * scale, 0.6 * scale),
        point(0.0, 1.0 * scale),
        point(-0.8 * scale, 0.6 * scale),
        point(-1.0 * scale, 0.0),
        point(-0.8 * scale, -0.6 * scale),
        point(0.0, -1.0 * scale),
        point(0.8 * scale, -0.6 * scale)
    ]

    # Create closed Catmull-Rom spline
    spline: catmullrom = catmullrom(pts, true, 0.5)

    # Convert to region and extrude
    profile: region2d = region_from_spline(spline, 64)
    result: solid = extrude(profile, scale * 0.5)
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
    solid = result.emit_result.data  # The yapCAD solid
    # Use the solid...
else:
    print(f"Error: {result.error_message}")
```

## See Also

- `docs/yapCADone.rst` - yapCAD 1.0 roadmap and feature status
- `docs/ycpkg_spec.rst` - Package specification
- `examples/` - Example DSL files
