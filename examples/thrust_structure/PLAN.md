# Thrust Structure FEA Demo

## Overview

This demo showcases the complete yapCAD workflow for designing, analyzing, and optimizing a liquid fuel rocket thrust structure. The workflow demonstrates:

1. **Parametric DSL Design** - Geometry defined in yapCAD DSL with optimization parameters
2. **Requirements-Driven Analysis** - FEA validation against acceptance criteria
3. **Mass Budget Check** - Native yapCAD geometric validation (no external deps)
4. **Package-Based Workflow** - `.ycpkg` format bundling geometry, requirements, and results

## Design Specifications

### Geometry
- **Outer diameter**: 12 inches (304.8 mm)
- **Thickness**: 8 mm
- **Material**: 6061-T6 Aluminum
- **Motor mount**: 4 inch (101.6 mm) diameter central hole
- **Stringer notches**: 3x rectangular (2" x 1" / 50.8 x 25.4 mm) at 120 degree intervals
- **Lightening holes**: Parameterized pattern between motor mount and edge
- **Symmetry**: Trilateral (120 degree rotational)

### Requirements
| Requirement | Limit | Actual (Optimized) | Margin |
|-------------|-------|-------------------|--------|
| **Mass Budget** | < 1.0 kg | 0.929 kg | 7.1% |
| **Max Deflection** | < 1.0 mm | 0.512 mm | 48.8% |
| **Thrust Load** | 1500 lbf (6672 N) | - | - |

### Constraint
Fixed at stringer notch locations (simulating welded stringer attachment)

### Optimization Goal
Minimize mass while meeting deflection and mass budget requirements.

## Architecture

```
thrust_structure/
├── PLAN.md                          # This file
├── thrust_structure.dsl             # Parametric geometry definition
├── run_validation.py                # Primary validation script
├── run_demo.sh                      # End-to-end demo script
├── optimize.py                      # Optimization driver (for exploration)
├── visualize.py                     # Viewer/export helper
├── validation/
│   └── plans/
│       ├── thrust-fea.yaml          # FEA analysis plan
│       └── mass-check.yaml          # Mass budget check plan
└── thrust_structure.ycpkg/          # Output yapCAD package
    ├── manifest.yaml                # Package manifest
    ├── geometry/
    │   └── primary.json             # Generated geometry
    └── validation/
        ├── plans/
        │   ├── thrust-fea.yaml
        │   └── mass-check.yaml
        └── results/
            ├── mass-check/
            │   └── summary.json     # Mass check results
            └── thrust-fea/
                ├── summary.json     # FEA results
                ├── mesh.msh         # Gmsh mesh
                ├── mesh.xdmf        # XDMF mesh for DOLFINx
                ├── displacement.xdmf # Displacement field
                └── stress.xdmf      # Von Mises stress field
```

## Usage

### Quick Start (Full Validation)
```bash
# Run both mass check and FEA analysis
./run_demo.sh

# Or run validation directly
python run_validation.py --design optimized
```

### Mass Check Only (No External Dependencies)
```bash
python run_validation.py --mass-only
```

### Generate STEP Files for CAD
```bash
python -m yapcad.dsl run thrust_structure.dsl MAKE_OPTIMIZED_PLATE \
    --output output/optimized_plate.step
```

### Available Design Commands
| Command | Description | Holes |
|---------|-------------|-------|
| `MAKE_BASELINE_PLATE` | Solid plate, no lightening | 0 |
| `MAKE_LIGHTENED_N(n)` | Parameterized lightening | n per sector |
| `MAKE_LIGHTENED_PLATE_3` | Convenience wrapper | 9 total |
| `MAKE_LIGHTENED_PLATE_4` | Convenience wrapper | 12 total |
| `MAKE_OPTIMIZED_PLATE` | Best configuration found | 12 total |

## Validation Plans

### Mass Budget Check (`mass-check.yaml`)
- **Backend**: yapcad (native, no external dependencies)
- **Check**: Volume * density < 1.0 kg
- **Material**: 6061-T6 Aluminum (2700 kg/m^3)

### FEA Analysis (`thrust-fea.yaml`)
- **Backend**: fenics (DOLFINx)
- **Analysis**: Linear elastic static
- **Load**: 2.61 MPa pressure on motor mount inner surface
- **Boundary**: Fixed at stringer notches
- **Acceptance**: Max displacement < 1.0 mm

## Material Properties (6061-T6 Aluminum)

| Property | Value | Units |
|----------|-------|-------|
| Young's Modulus | 68.9 | GPa |
| Poisson's Ratio | 0.33 | - |
| Density | 2700 | kg/m^3 |
| Yield Strength | 276 | MPa |
| Ultimate Strength | 310 | MPa |

## Results Summary

### Optimized Design Performance
- **Mass**: 0.929 kg (vs 1.0 kg budget)
- **Max Displacement**: 0.512 mm (vs 1.0 mm limit)
- **Max Von Mises Stress**: 103.1 MPa (vs 276 MPa yield)
- **Safety Factor**: ~2.7x on stress

### Mesh Statistics
- Nodes: ~6,400
- Elements: ~19,000 tetrahedra
- Element size: 5.0 mm target

## Dependencies

### Required (for mass check only)
- yapCAD with OCC BREP support

### Optional (for FEA analysis)
```bash
conda install -c conda-forge fenics-dolfinx gmsh meshio pyvista
```

### Optional (for visualization)
- ParaView: `paraview displacement.xdmf`

## Notes

### Unit Consistency
All internal calculations in SI (meters, Pascals, Newtons).
User-facing parameters use mm for convenience.

### Gmsh Healing
The Gmsh mesher is configured to skip automatic healing during BREP import,
as this can break valid geometry topology. OCC healing is applied at the
BREP level before export.

### Design Margin
The optimized design has significant margin on both the deflection limit
(64% margin) and stress (3.7x safety factor). Further mass reduction is
possible if requirements allow.
