"""Analytic surfaces demo: exports one of each surface type to STEP.

This example creates one solid for each class of analytic surface supported
by the native BREP system and exports them to STEP files. The user can choose
between:
1. Analytic form (native BREP with exact geometric definitions)
2. Faceted form (triangle mesh approximation)

Supported surface types:
- Plane (box has 6 planar faces)
- Cylinder (cylindrical surface with circular end caps)
- Sphere (spherical surface)
- Cone (conical surface with circular base)
- Torus (toroidal surface - donut shape)

Usage:
    # Export all surfaces using analytic BREP (default)
    python analytic_surfaces_demo.py --mode analytic --output build/analytic

    # Export all surfaces using faceted mesh
    python analytic_surfaces_demo.py --mode faceted --output build/faceted

    # Export both for comparison
    python analytic_surfaces_demo.py --mode both --output build/surfaces

    # View in 3D viewer
    python analytic_surfaces_demo.py --mode view
"""

import argparse
import sys
from pathlib import Path

# Check for OCC availability early
try:
    from OCC.Core.BRepPrimAPI import (
        BRepPrimAPI_MakeBox,
        BRepPrimAPI_MakeCylinder,
        BRepPrimAPI_MakeSphere,
        BRepPrimAPI_MakeCone,
        BRepPrimAPI_MakeTorus,
    )
    from OCC.Core.gp import gp_Pnt, gp_Ax2, gp_Dir
    OCC_AVAILABLE = True
except ImportError:
    OCC_AVAILABLE = False
    print("Warning: pythonocc-core not available. Install it for full functionality.")
    print("         conda install -c conda-forge pythonocc-core")


def create_box(size=10.0, offset=(0, 0, 0)):
    """Create a box with 6 planar faces.

    Surface types: 6x plane_surface
    """
    if not OCC_AVAILABLE:
        return None

    corner = gp_Pnt(offset[0], offset[1], offset[2])
    shape = BRepPrimAPI_MakeBox(corner, size, size, size).Shape()
    return shape, "box", f"Box {size}x{size}x{size} (6 planar faces)"


def create_cylinder(radius=5.0, height=15.0, offset=(0, 0, 0)):
    """Create a cylinder with cylindrical surface and planar end caps.

    Surface types: 1x cylinder_surface, 2x plane_surface (end caps)
    """
    if not OCC_AVAILABLE:
        return None

    axis = gp_Ax2(gp_Pnt(offset[0], offset[1], offset[2]), gp_Dir(0, 0, 1))
    shape = BRepPrimAPI_MakeCylinder(axis, radius, height).Shape()
    return shape, "cylinder", f"Cylinder r={radius}, h={height} (1 cylindrical, 2 planar)"


def create_sphere(radius=8.0, offset=(0, 0, 0)):
    """Create a sphere with a single spherical surface.

    Surface types: 1x sphere_surface
    """
    if not OCC_AVAILABLE:
        return None

    center = gp_Pnt(offset[0], offset[1], offset[2])
    shape = BRepPrimAPI_MakeSphere(center, radius).Shape()
    return shape, "sphere", f"Sphere r={radius} (1 spherical surface)"


def create_cone(base_radius=8.0, top_radius=2.0, height=12.0, offset=(0, 0, 0)):
    """Create a truncated cone with conical surface and planar end caps.

    Surface types: 1x cone_surface, 2x plane_surface (end caps)
    """
    if not OCC_AVAILABLE:
        return None

    axis = gp_Ax2(gp_Pnt(offset[0], offset[1], offset[2]), gp_Dir(0, 0, 1))
    shape = BRepPrimAPI_MakeCone(axis, base_radius, top_radius, height).Shape()
    return shape, "cone", f"Cone r1={base_radius}, r2={top_radius}, h={height} (1 conical, 2 planar)"


def create_torus(major_radius=10.0, minor_radius=3.0, offset=(0, 0, 0)):
    """Create a torus (donut shape) with toroidal surface.

    Surface types: 1x torus_surface
    """
    if not OCC_AVAILABLE:
        return None

    axis = gp_Ax2(gp_Pnt(offset[0], offset[1], offset[2]), gp_Dir(0, 0, 1))
    shape = BRepPrimAPI_MakeTorus(axis, major_radius, minor_radius).Shape()
    return shape, "torus", f"Torus R={major_radius}, r={minor_radius} (1 toroidal surface)"


def export_analytic(shape, output_path, name="shape"):
    """Export shape using analytic BREP (exact geometry).

    The STEP file will contain true geometric definitions:
    PLANE(), CYLINDRICAL_SURFACE(), SPHERICAL_SURFACE(), etc.
    """
    from OCC.Core.STEPControl import STEPControl_Writer, STEPControl_AsIs
    from OCC.Core.IFSelect import IFSelect_RetDone
    from OCC.Core.Interface import Interface_Static

    writer = STEPControl_Writer()
    Interface_Static.SetCVal("write.step.product.name", name)

    status = writer.Transfer(shape, STEPControl_AsIs)
    if status != IFSelect_RetDone:
        raise RuntimeError(f"Failed to transfer shape: {status}")

    status = writer.Write(str(output_path))
    if status != IFSelect_RetDone:
        raise RuntimeError(f"Failed to write STEP file: {status}")

    return True


def export_faceted(shape, output_path, name="shape", linear_deflection=0.1):
    """Export shape using faceted BREP (triangle mesh).

    The STEP file will contain triangulated geometry.
    """
    from yapcad.brep import BrepSolid
    from yapcad.io.step import write_step
    from yapcad.geom3d import solid

    # Wrap in BrepSolid and tessellate
    brep_solid = BrepSolid(shape)
    surf = brep_solid.tessellate(deflection=linear_deflection)

    # Create a solid from the surface for STEP export
    # The surface from tessellate() is already in the right format
    write_step(surf, output_path, name=name)
    return True


def count_surface_types(step_file):
    """Parse a STEP file and count surface types present."""
    with open(step_file, 'r') as f:
        content = f.read()

    surface_types = {
        'PLANE': content.count('PLANE('),
        'CYLINDRICAL_SURFACE': content.count('CYLINDRICAL_SURFACE('),
        'SPHERICAL_SURFACE': content.count('SPHERICAL_SURFACE('),
        'CONICAL_SURFACE': content.count('CONICAL_SURFACE('),
        'TOROIDAL_SURFACE': content.count('TOROIDAL_SURFACE('),
        'B_SPLINE_SURFACE': content.count('B_SPLINE_SURFACE'),
    }

    # Also check for faceted indicators
    surface_types['TRIANGLES'] = content.count('ADVANCED_FACE')

    return {k: v for k, v in surface_types.items() if v > 0}


def view_shapes(shapes_dict):
    """View shapes using yapCAD's pyglet viewer."""
    from yapcad.brep import BrepSolid
    from yapcad.geometry import Geometry
    from yapcad.pyglet_drawable import pygletDraw

    viewer = pygletDraw()

    colors = ['cyan', 'yellow', 'magenta', 'green', 'orange']

    for idx, (name, (shape, _, desc)) in enumerate(shapes_dict.items()):
        if shape is None:
            continue

        color = colors[idx % len(colors)]
        brep = BrepSolid(shape)
        geom = Geometry(brep)

        viewer.make_object(name, lighting=True, linecolor=color, material="pearl")

        # Get faceted surface for display
        surf = geom.surface()
        viewer.draw_surface(surf, name=name)

    viewer.display()


def main():
    parser = argparse.ArgumentParser(
        description="Export analytic surface types to STEP for comparison testing"
    )
    parser.add_argument(
        "--mode",
        choices=["analytic", "faceted", "both", "view"],
        default="analytic",
        help="Export mode: analytic (exact), faceted (mesh), both, or view"
    )
    parser.add_argument(
        "--output",
        default="build/analytic_surfaces",
        help="Output directory for STEP files"
    )
    parser.add_argument(
        "--deflection",
        type=float,
        default=0.1,
        help="Linear deflection for faceted meshing (smaller = finer)"
    )
    args = parser.parse_args()

    if not OCC_AVAILABLE:
        print("Error: pythonocc-core is required for this demo.")
        print("Install with: conda install -c conda-forge pythonocc-core")
        sys.exit(1)

    # Create all primitive shapes with offsets so they don't overlap in viewer
    shapes = {
        'box': create_box(10.0, offset=(0, 0, 0)),
        'cylinder': create_cylinder(5.0, 15.0, offset=(25, 0, 0)),
        'sphere': create_sphere(8.0, offset=(50, 0, 8)),
        'cone': create_cone(8.0, 2.0, 12.0, offset=(75, 0, 0)),
        'torus': create_torus(10.0, 3.0, offset=(100, 0, 3)),
    }

    print("Analytic Surface Types Demo")
    print("=" * 60)
    print()

    for name, (shape, short_name, desc) in shapes.items():
        print(f"  {name}: {desc}")
    print()

    if args.mode == "view":
        print("Opening 3D viewer...")
        view_shapes(shapes)
        return

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.mode in ("analytic", "both"):
        analytic_dir = output_dir / "analytic" if args.mode == "both" else output_dir
        analytic_dir.mkdir(parents=True, exist_ok=True)

        print("Exporting ANALYTIC (exact geometry) STEP files:")
        print("-" * 50)

        for name, (shape, short_name, desc) in shapes.items():
            if shape is None:
                continue

            step_path = analytic_dir / f"{short_name}_analytic.step"
            try:
                export_analytic(shape, step_path, name=short_name)
                surface_types = count_surface_types(step_path)
                print(f"  {step_path}")
                print(f"    Surfaces: {surface_types}")
            except Exception as e:
                print(f"  {step_path} - FAILED: {e}")

        print()

    if args.mode in ("faceted", "both"):
        faceted_dir = output_dir / "faceted" if args.mode == "both" else output_dir
        faceted_dir.mkdir(parents=True, exist_ok=True)

        print("Exporting FACETED (triangle mesh) STEP files:")
        print("-" * 50)

        for name, (shape, short_name, desc) in shapes.items():
            if shape is None:
                continue

            step_path = faceted_dir / f"{short_name}_faceted.step"
            try:
                export_faceted(shape, step_path, name=short_name,
                             linear_deflection=args.deflection)
                surface_types = count_surface_types(step_path)
                print(f"  {step_path}")
                print(f"    Surfaces: {surface_types}")
            except Exception as e:
                print(f"  {step_path} - FAILED: {e}")

        print()

    print("=" * 60)
    print("Export complete!")
    print()
    print("To test in FreeCAD:")
    print("  1. Open FreeCAD")
    print("  2. File > Import... > select a STEP file")
    print("  3. In Model panel, expand the shape tree")
    print("  4. Check face types - analytic should show 'Plane', 'Cylinder', etc.")
    print("     while faceted will show triangulated faces")


if __name__ == "__main__":
    main()
