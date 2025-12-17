# yapCAD DSL Tutorial

This tutorial walks through creating a parametric part using the yapCAD DSL, from design to export.

## Prerequisites

- Python 3.10+ with yapCAD installed
- For STEP export: pythonocc-core (via conda)

```bash
# Install yapCAD (if not already)
pip install -e .

# For BREP/STEP support (optional but recommended)
conda install -c conda-forge pythonocc-core
```

## Tutorial: Creating a Parametric Pipe Fitting

We'll create a pipe fitting with:
- A main body (hollow cylinder)
- Mounting flanges with bolt holes
- Parametric dimensions

### Step 1: Create the DSL File

Create a file `pipe_fitting.dsl`:

```python
# Pipe Fitting - Parametric Design
# A hollow cylindrical body with mounting flanges

module pipe_fitting

# Helper: create a mounting flange with bolt holes
command make_flange(
    outer_diameter: float,
    inner_diameter: float,
    thickness: float,
    bolt_circle_diameter: float,
    bolt_hole_diameter: float,
    bolt_count: int
) -> solid:
    # Create the base flange disc
    outer_radius: float = outer_diameter / 2.0
    flange_cyl: solid = cylinder(outer_radius, thickness)

    # Create the center hole
    inner_radius: float = inner_diameter / 2.0
    center_hole: solid = cylinder(inner_radius, thickness + 2.0)
    center_hole_pos: solid = translate(center_hole, 0.0, 0.0, -1.0)

    # Start with flange minus center hole
    flange: solid = difference(flange_cyl, center_hole_pos)

    # Add bolt holes around the bolt circle
    bolt_radius: float = bolt_hole_diameter / 2.0
    bolt_circle_radius: float = bolt_circle_diameter / 2.0

    # Create bolt holes (manual placement for 4 bolts)
    # For 4 bolts at 0, 90, 180, 270 degrees
    bolt_hole: solid = cylinder(bolt_radius, thickness + 2.0)
    bolt_hole_z: solid = translate(bolt_hole, 0.0, 0.0, -1.0)

    # Position at 0 degrees (positive X)
    hole1: solid = translate(bolt_hole_z, bolt_circle_radius, 0.0, 0.0)
    flange = difference(flange, hole1)

    # Position at 90 degrees (positive Y)
    hole2: solid = translate(bolt_hole_z, 0.0, bolt_circle_radius, 0.0)
    flange = difference(flange, hole2)

    # Position at 180 degrees (negative X)
    neg_bcr: float = 0.0 - bolt_circle_radius
    hole3: solid = translate(bolt_hole_z, neg_bcr, 0.0, 0.0)
    flange = difference(flange, hole3)

    # Position at 270 degrees (negative Y)
    hole4: solid = translate(bolt_hole_z, 0.0, neg_bcr, 0.0)
    flange = difference(flange, hole4)

    emit flange

# Main command: create the complete pipe fitting
command MAKE_PIPE_FITTING(
    body_length: float = 100.0,
    body_outer_diameter: float = 50.0,
    body_inner_diameter: float = 40.0,
    flange_outer_diameter: float = 80.0,
    flange_thickness: float = 10.0,
    bolt_circle_diameter: float = 65.0,
    bolt_hole_diameter: float = 8.0
) -> solid:
    # Validate parameters
    require body_outer_diameter > body_inner_diameter
    require flange_outer_diameter > body_outer_diameter
    require bolt_circle_diameter < flange_outer_diameter
    require bolt_circle_diameter > body_outer_diameter

    # Create the main body (hollow cylinder)
    body_outer_radius: float = body_outer_diameter / 2.0
    body_inner_radius: float = body_inner_diameter / 2.0

    outer_cyl: solid = cylinder(body_outer_radius, body_length)
    inner_cyl: solid = cylinder(body_inner_radius, body_length + 2.0)
    inner_cyl_pos: solid = translate(inner_cyl, 0.0, 0.0, -1.0)

    body: solid = difference(outer_cyl, inner_cyl_pos)

    # Create flanges
    flange: solid = make_flange(
        flange_outer_diameter,
        body_inner_diameter,
        flange_thickness,
        bolt_circle_diameter,
        bolt_hole_diameter,
        4
    )

    # Position bottom flange (already at Z=0)
    bottom_flange: solid = flange

    # Position top flange at top of body
    top_flange: solid = translate(flange, 0.0, 0.0, body_length - flange_thickness)

    # Combine all parts
    fitting: solid = union(body, bottom_flange)
    fitting = union(fitting, top_flange)

    emit fitting

# Variant: just the body without flanges
command MAKE_PIPE_BODY(
    length: float = 100.0,
    outer_diameter: float = 50.0,
    wall_thickness: float = 5.0
) -> solid:
    require wall_thickness > 0.0
    require wall_thickness < outer_diameter / 2.0

    outer_radius: float = outer_diameter / 2.0
    inner_radius: float = outer_radius - wall_thickness

    outer_cyl: solid = cylinder(outer_radius, length)
    inner_cyl: solid = cylinder(inner_radius, length + 2.0)
    inner_cyl_pos: solid = translate(inner_cyl, 0.0, 0.0, -1.0)

    body: solid = difference(outer_cyl, inner_cyl_pos)
    emit body
```

### Step 2: Verify the DSL Syntax

Check the file for syntax and type errors:

```bash
python -m yapcad.dsl check pipe_fitting.dsl
```

Expected output:
```
OK: pipe_fitting.dsl passed type checking
```

### Step 3: List Available Commands

See what commands are available:

```bash
python -m yapcad.dsl list pipe_fitting.dsl
```

Expected output:
```
Commands in pipe_fitting.dsl:
  MAKE_PIPE_FITTING(body_length: float = 100.0, body_outer_diameter: float = 50.0, ...) -> solid
  MAKE_PIPE_BODY(length: float = 100.0, outer_diameter: float = 50.0, wall_thickness: float = 5.0) -> solid
```

### Step 4: Generate Geometry

Run the command with default parameters:

```bash
python -m yapcad.dsl run pipe_fitting.dsl MAKE_PIPE_FITTING --output pipe_fitting.step
```

Or with custom parameters:

```bash
python -m yapcad.dsl run pipe_fitting.dsl MAKE_PIPE_FITTING \
    --param body_length=150.0 \
    --param body_outer_diameter=60.0 \
    --output pipe_fitting_large.step
```

### Step 5: Create a Package

Create a full yapCAD package with provenance tracking:

```bash
python -m yapcad.dsl run pipe_fitting.dsl MAKE_PIPE_FITTING \
    --package pipe_fitting.ycpkg
```

This creates:
```
pipe_fitting.ycpkg/
├── manifest.yaml          # Package metadata
├── geometry/
│   └── primary.json       # yapCAD geometry format
├── exports/
│   └── model.step         # STEP export
└── src/
    └── pipe_fitting.dsl   # Source DSL (if packaged)
```

### Step 6: View the Package

View the generated geometry:

```bash
python tools/ycpkg_viewer.py pipe_fitting.ycpkg
```

Viewer controls:
- Mouse drag: Rotate view
- Scroll: Zoom
- Arrow keys: Pan
- `1-9`: Toggle layers
- `H` or `F1`: Help overlay

### Step 7: Validate the Package

Check package integrity:

```bash
python tools/ycpkg_validate.py pipe_fitting.ycpkg
```

## Advanced: Swept Geometry

For parts with non-straight paths, use sweep operations:

```python
module bent_tube

# Create a path with a 90-degree bend
command make_elbow_path(
    straight_length: float,
    bend_radius: float
) -> path3d:
    # First straight segment along Y
    p0: point3d = point(0.0, 0.0, 0.0)
    p1: point3d = point(0.0, straight_length, 0.0)

    # After bend, go along X
    p2: point3d = point(straight_length, straight_length + bend_radius, 0.0)

    seg1: path3d = path3d_line(p0, p1)
    seg2: path3d = path3d_line(p1, p2)

    segments: list<path3d> = [seg1, seg2]
    emit make_path3d(segments)

command MAKE_ELBOW(
    outer_diameter: float = 50.0,
    wall_thickness: float = 3.0,
    straight_length: float = 50.0
) -> solid:
    # Create profiles
    outer_radius: float = outer_diameter / 2.0
    inner_radius: float = outer_radius - wall_thickness

    outer_profile: region2d = regular_polygon(32, outer_radius)
    inner_profile: region2d = regular_polygon(32, inner_radius)

    # Create path
    spine: path3d = make_elbow_path(straight_length, outer_diameter)

    # Sweep with adaptive tangent tracking
    result: solid = sweep_adaptive_hollow(outer_profile, inner_profile, spine, 5.0)
    emit result
```

## Workflow Summary

1. **Design**: Write DSL in `.dsl` file
2. **Check**: `python -m yapcad.dsl check file.dsl`
3. **List**: `python -m yapcad.dsl list file.dsl`
4. **Run**: `python -m yapcad.dsl run file.dsl COMMAND --output file.step`
5. **Package**: `python -m yapcad.dsl run file.dsl COMMAND --package output.ycpkg`
6. **View**: `python tools/ycpkg_viewer.py output.ycpkg`
7. **Validate**: `python tools/ycpkg_validate.py output.ycpkg`

## Tips

### Debugging

Use `print()` for debugging intermediate values:

```python
command DEBUG_EXAMPLE(x: float) -> solid:
    intermediate: float = x * 2.0
    print("intermediate value:", intermediate)
    # ...
```

### Parameter Validation

Always validate parameters with `require`:

```python
command SAFE_CYLINDER(radius: float, height: float) -> solid:
    require radius > 0.0
    require height > 0.0
    # ...
```

### Naming Conventions

- `UPPERCASE_NAMES`: Exported commands (visible in CLI)
- `lowercase_names`: Helper commands (internal use)

### Boolean Operations

When combining many parts, chain unions efficiently:

```python
# Instead of:
# result = union(union(union(a, b), c), d)

# Do:
result: solid = a
result = union(result, b)
result = union(result, c)
result = union(result, d)
```

## Next Steps

- See `docs/dsl_reference.md` for complete function reference
- See `examples/` for more example DSL files
- See `docs/ycpkg_spec.rst` for package format details
