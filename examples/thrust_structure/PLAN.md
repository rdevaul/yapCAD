# Thrust Structure FEA Demo - Implementation Plan

## Overview

This demo showcases the complete yapCAD workflow for designing, analyzing, and optimizing a liquid fuel rocket thrust structure. The workflow demonstrates:

1. **Parametric DSL Design** - Geometry defined in yapCAD DSL with optimization parameters
2. **Requirements-Driven Analysis** - FEA validation against acceptance criteria
3. **Package-Based Workflow** - `.ycpkg` format bundling geometry, requirements, and results
4. **Design Optimization** - Iterative search for minimum-weight solution

## Design Specifications

### Geometry
- **Outer diameter**: 12 inches (304.8 mm)
- **Thickness**: 8 mm
- **Material**: 6061-T6 Aluminum
- **Motor mount**: 4 inch (101.6 mm) diameter central hole
- **Stringer notches**: 3x rectangular (2" × 1" / 50.8 × 25.4 mm) at 120° intervals
- **Lightening holes**: Parameterized pattern between motor mount and edge
- **Symmetry**: Trilateral (120° rotational)

### Requirements
- **Thrust load**: 1500 lbf (6672 N) applied at motor mount
- **Max deflection**: 10 mm at any point
- **Constraint**: Fixed at stringer notch locations (simulating stringer attachment)

### Optimization Goal
Minimize mass while meeting deflection requirement.

## Architecture

```
thrust_structure/
├── PLAN.md                          # This file
├── thrust_structure.dsl             # Parametric geometry definition
├── thrust_structure.ycpkg/          # yapCAD package
│   ├── manifest.yaml                # Package manifest
│   ├── geometry/
│   │   └── primary.json             # Generated geometry
│   ├── validation/
│   │   ├── plans/
│   │   │   └── thrust-fea.yaml      # FEA analysis plan
│   │   └── results/
│   │       └── thrust-fea/
│   │           ├── summary.json     # Analysis results
│   │           ├── mesh.msh         # Gmsh mesh
│   │           └── displacement.vtu # ParaView visualization
│   └── attachments/
│       └── source.dsl               # Original DSL source
├── optimize.py                      # Optimization driver script
└── visualize.py                     # Viewer/export script
```

## Implementation Steps

### Phase 1: DSL Geometry Module

Create `thrust_structure.dsl` with:

```
module thrust_structure

# Parameters for optimization
command MAKE_THRUST_PLATE(
    # Fixed geometry
    outer_diameter_mm: float = 304.8,    # 12 inches
    thickness_mm: float = 8.0,
    motor_mount_diameter_mm: float = 101.6,  # 4 inches
    stringer_width_mm: float = 50.8,     # 2 inches
    stringer_depth_mm: float = 25.4,     # 1 inch

    # Optimization parameters
    num_holes_per_sector: int = 3,       # Holes in each 120° sector
    hole_radius_mm: float = 20.0,        # Lightening hole radius
    hole_radial_position: float = 0.6,   # 0-1 fraction from center to edge
    hole_angular_spread: float = 0.5     # 0-1 fraction of sector used
) -> solid
```

Key features:
- Trilateral symmetry via 120° pattern
- Parameterized lightening hole pattern
- Material assignment for density calculation

### Phase 2: Analysis Plan

Create `thrust-fea.yaml`:

```yaml
id: thrust-fea
kind: structural
backend: fenics
name: Thrust Structure Static Analysis

geometry:
  source: geometry/primary.json

materials:
  default:
    name: 6061-T6
    youngs_modulus_pa: 68.9e9
    poisson_ratio: 0.33
    density_kgm3: 2700

loads:
  - id: thrust
    type: pressure
    surfaces: motor_mount_face
    magnitude_pa: calculated  # 6672N / area

boundary_conditions:
  - id: stringer_fixed
    type: fixed
    surfaces: stringer_notch_faces

acceptance:
  displacement.max_mm:
    limit: 10.0
    comparison: "<="

backendOptions:
  mesh:
    element_size: 5.0
    min_element_size: 2.0
    element_order: 1
```

### Phase 3: Package Assembly

Script to:
1. Compile DSL → geometry
2. Create package structure
3. Embed analysis plan
4. Compute mass from volume × density

### Phase 4: Optimization Loop

```python
def optimize_thrust_structure():
    """Find minimum-mass design meeting deflection requirement."""

    # Parameter space
    params = {
        'num_holes_per_sector': [2, 3, 4, 5],
        'hole_radius_mm': [15, 20, 25, 30],
        'hole_radial_position': [0.5, 0.6, 0.7],
        'hole_angular_spread': [0.4, 0.5, 0.6]
    }

    best_mass = float('inf')
    best_config = None

    for config in parameter_combinations(params):
        # Generate geometry
        # Run FEA
        # Check acceptance
        # Track best

    return best_config
```

### Phase 5: Visualization & Export

- **Viewer**: 4-view display with material colors
- **STEP export**: For CAM/manufacturing
- **STL export**: For 3D printing prototypes
- **ParaView**: Load `.vtu` for stress/displacement visualization

## DSL Design Details

### Geometry Construction Approach

```
# Single sector (120°) construction:
1. Start with full disk
2. Subtract motor mount hole
3. Subtract stringer notch
4. Subtract lightening holes (if any)

# Full plate:
1. Create single sector
2. Rotate copy at 120°
3. Rotate copy at 240°
4. Union all three sectors
```

### Face Naming for BCs

```
face_names:
  motor_mount_face: inner cylindrical surface of motor hole
  stringer_notch_faces: faces of rectangular notches
  top_face: +Z face of plate
  bottom_face: -Z face of plate
```

### Material Properties

| Property | Value | Units |
|----------|-------|-------|
| Young's Modulus | 68.9 | GPa |
| Poisson's Ratio | 0.33 | - |
| Density | 2700 | kg/m³ |
| Yield Strength | 276 | MPa |

## Validation Approach

### Load Application
- Thrust force: 6672 N (1500 lbf)
- Applied as pressure on motor mount inner face
- Pressure = Force / Area = 6672 / (π × 50.8 × 8) = ~5.2 MPa

### Boundary Conditions
- Stringer notches fixed (all DOF = 0)
- Simulates welded/bolted stringer attachment

### Acceptance Criteria
- Max displacement anywhere < 10 mm
- (Optional) Max von Mises stress < yield strength with safety factor

## CLI Usage

```bash
# Check DSL syntax
python -m yapcad.dsl check thrust_structure.dsl

# Generate geometry with default parameters
python -m yapcad.dsl run thrust_structure.dsl MAKE_THRUST_PLATE \
    --output thrust_plate.step

# Generate with custom parameters
python -m yapcad.dsl run thrust_structure.dsl MAKE_THRUST_PLATE \
    --param num_holes_per_sector=4 \
    --param hole_radius_mm=25.0 \
    --output thrust_plate_optimized.step

# Create package
python -m yapcad.dsl run thrust_structure.dsl MAKE_THRUST_PLATE \
    --package thrust_structure.ycpkg \
    --name "Thrust Structure" \
    --version "1.0.0"

# Run FEA analysis
python -m yapcad.package.analysis.cli thrust_structure.ycpkg \
    --plan validation/plans/thrust-fea.yaml

# View results
python visualize.py thrust_structure.ycpkg
```

## Optimization Strategy

### Approach 1: Grid Search
- Enumerate parameter combinations
- Run FEA on each
- Select minimum mass that passes

### Approach 2: Bayesian Optimization (Future)
- Use Gaussian process surrogate
- Balance exploration/exploitation
- More efficient for larger parameter spaces

### Approach 3: Gradient-Based (Future)
- Sensitivity analysis via FEA
- Shape optimization
- Requires adjoint solver support

## Expected Results

### Baseline (No Lightening Holes)
- Mass: ~1.6 kg
- Max deflection: ~2-3 mm (well under limit)
- Over-engineered, room for optimization

### Optimized
- Target: <1.0 kg
- Max deflection: ~8-9 mm (near limit)
- Balance of structural efficiency

## Files to Create

1. `thrust_structure.dsl` - Parametric geometry
2. `thrust-fea.yaml` - Analysis plan
3. `optimize.py` - Optimization driver
4. `visualize.py` - Viewer/export helper
5. `run_demo.sh` - End-to-end demo script

## Dependencies

- yapCAD with OCC BREP support
- Gmsh (for meshing)
- FEniCSx/DOLFINx (for FEA) - optional, can use CalculiX fallback
- ParaView (for result visualization) - optional

## Notes

### Unit Consistency
All internal calculations in SI (meters, Pascals, Newtons).
User-facing parameters can use mm for convenience.

### Symmetry Exploitation
Could mesh only 1/3 of structure with symmetry BCs for faster analysis.
Full model used here for clarity.

### Mesh Convergence
Should verify results with finer mesh to ensure convergence.
Element size ~5mm is starting point.
