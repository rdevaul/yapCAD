"""Conic curves demo: parabola and hyperbola primitives.

This example demonstrates:
1. Creating parabola and hyperbola curves using yapCAD's conic primitives
2. Exporting 2D curves to DXF format
3. Revolving the curves to create 3D solids (paraboloid, hyperboloid)
4. Exporting the 3D solids to STEP format

Usage:
    python conic_demo.py                       # View 3D solids
    python conic_demo.py --mode dxf            # Export 2D curves to DXF
    python conic_demo.py --mode step           # Export 3D solids to STEP
    python conic_demo.py --mode all --output build/conic  # All outputs
"""

import argparse
import math
from pathlib import Path

import ezdxf

from yapcad.geom import (
    point, sample, parabola, hyperbola, isparabola, ishyperbola
)
from yapcad.geom3d_util import makeRevolutionSolid
from yapcad.io import write_step


def create_parabola_curve(focal_length=1.0, t_range=(-3.0, 3.0)):
    """Create a parabola curve with specified focal length.

    The parabola opens along positive X with vertex at origin.
    Parametric form: x = t²/(4f), y = t
    """
    return parabola(
        point(0, 0, 0),
        focal_length,
        start=t_range[0],
        end=t_range[1]
    )


def create_hyperbola_curve(semi_major=2.0, semi_minor=1.5, t_range=(-1.5, 1.5)):
    """Create a hyperbola curve (right branch).

    Standard form: x²/a² - y²/b² = 1
    Parametric form: x = a*cosh(t), y = b*sinh(t)
    """
    return hyperbola(
        point(0, 0, 0),
        semi_major,
        semi_minor,
        start=t_range[0],
        end=t_range[1],
        branch=1  # Right branch (positive x)
    )


def sample_curve_to_points(curve, num_samples=64):
    """Sample a curve to a list of 2D points for DXF export."""
    pts = []
    for i in range(num_samples + 1):
        u = i / num_samples
        pt = sample(curve, u)
        pts.append((pt[0], pt[1]))
    return pts


def export_curves_to_dxf(output_path, curves_dict):
    """Export curves to DXF file with labeled layers.

    Args:
        output_path: Path to the output DXF file
        curves_dict: Dict mapping layer names to curves
    """
    doc = ezdxf.new(dxfversion="R2010", setup=True)
    msp = doc.modelspace()

    colors = {"parabola": 1, "hyperbola": 3, "axis": 7}  # DXF color indices

    for name, curve in curves_dict.items():
        layer_name = name.upper()
        if layer_name not in doc.layers:
            doc.layers.add(name=layer_name, color=colors.get(name, 7))

        pts = sample_curve_to_points(curve, num_samples=100)
        msp.add_lwpolyline(
            pts,
            format="xy",
            dxfattribs={"layer": layer_name}
        )

    # Add coordinate axes for reference
    axis_layer = "AXES"
    if axis_layer not in doc.layers:
        doc.layers.add(name=axis_layer, color=7)

    # X and Y axes
    msp.add_line((-5, 0), (5, 0), dxfattribs={"layer": axis_layer})
    msp.add_line((0, -5), (0, 5), dxfattribs={"layer": axis_layer})

    # Axis labels
    msp.add_text("X", height=0.3, dxfattribs={"layer": axis_layer}).set_placement((5.2, 0))
    msp.add_text("Y", height=0.3, dxfattribs={"layer": axis_layer}).set_placement((0.2, 5.2))

    doc.saveas(output_path)
    print(f"Wrote DXF: {output_path}")


def create_paraboloid(focal_length=1.0, z_max=4.0, steps=48, arc_samples=64):
    """Create a paraboloid solid of revolution.

    The paraboloid is created by rotating a parabola around the Z axis.
    For a parabola y² = 4fx (opening along X), rotating around Y gives
    a paraboloid. We use z as the height and r as the radius.

    Equation: r² = 4*f*z, so r(z) = sqrt(4*f*z) for z >= 0
    """
    def contour(z):
        if z < 0:
            return 0.0
        return math.sqrt(4.0 * focal_length * z)

    return makeRevolutionSolid(contour, 0.0, z_max, steps, arcSamples=arc_samples)


def create_hyperboloid_one_sheet(a=2.0, b=1.5, z_range=(-3.0, 3.0),
                                  steps=48, arc_samples=64):
    """Create a hyperboloid of one sheet solid of revolution.

    The hyperboloid of one sheet is created by rotating a hyperbola
    around its conjugate axis. For a hyperbola x²/a² - y²/b² = 1,
    rotating around the Y axis gives:

    Equation: x²/a² + z²/a² - y²/b² = 1
    Which means: r(z) = a*sqrt(1 + z²/b²)

    This creates the classic "cooling tower" shape.
    """
    def contour(z):
        return a * math.sqrt(1.0 + (z * z) / (b * b))

    z_start, z_end = z_range
    return makeRevolutionSolid(contour, z_start, z_end, steps, arcSamples=arc_samples)


def create_hyperboloid_two_sheets(a=2.0, b=1.5, z_range=(2.0, 5.0),
                                   steps=48, arc_samples=64):
    """Create one sheet of a hyperboloid of two sheets.

    For a hyperbola x²/a² - y²/b² = 1 rotated around the transverse axis (X),
    we get a hyperboloid of two sheets:

    Equation: x²/a² - y²/b² - z²/b² = 1
    Which means: r(z) = b*sqrt((x²/a²) - 1) where x > a

    We parameterize by x (height) and compute radius.
    """
    def contour(x):
        if x < a:
            return 0.001  # Minimum radius at vertex
        return b * math.sqrt((x * x) / (a * a) - 1.0)

    x_start, x_end = z_range
    return makeRevolutionSolid(contour, x_start, x_end, steps, arcSamples=arc_samples)


def export_solids_to_step(output_dir, solids_dict):
    """Export solids to STEP files.

    Args:
        output_dir: Directory for output files
        solids_dict: Dict mapping names to solid objects
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for name, solid in solids_dict.items():
        step_path = output_dir / f"{name}.step"
        write_step(solid, step_path, name=name)
        print(f"Wrote STEP: {step_path}")


def view_solids(solids_dict):
    """View solids using pyglet viewer."""
    from yapcad.pyglet_drawable import pygletDraw

    viewer = pygletDraw()

    # Colors for different solids
    colors = {
        "paraboloid": "cyan",
        "hyperboloid_one_sheet": "yellow",
        "hyperboloid_two_sheets": "magenta"
    }

    for name, solid in solids_dict.items():
        viewer.make_object(name, lighting=True,
                          linecolor=colors.get(name, "white"),
                          material="pearl")
        viewer.draw_solid(solid, name=name)

    viewer.display()


def main():
    parser = argparse.ArgumentParser(
        description="Conic curves demo: parabola and hyperbola primitives"
    )
    parser.add_argument(
        "--mode",
        choices=["view", "dxf", "step", "all"],
        default="view",
        help="Output mode: view 3D, export DXF, export STEP, or all"
    )
    parser.add_argument(
        "--output",
        default="build/conic",
        help="Output directory for exports"
    )
    parser.add_argument(
        "--focal-length", "-f",
        type=float,
        default=1.0,
        help="Parabola focal length"
    )
    parser.add_argument(
        "--semi-major", "-a",
        type=float,
        default=2.0,
        help="Hyperbola semi-major axis"
    )
    parser.add_argument(
        "--semi-minor", "-b",
        type=float,
        default=1.5,
        help="Hyperbola semi-minor axis"
    )
    args = parser.parse_args()

    output_dir = Path(args.output)

    # Create 2D curves
    print("Creating 2D conic curves...")
    parabola_curve = create_parabola_curve(
        focal_length=args.focal_length,
        t_range=(-3.0, 3.0)
    )
    hyperbola_curve = create_hyperbola_curve(
        semi_major=args.semi_major,
        semi_minor=args.semi_minor,
        t_range=(-1.5, 1.5)
    )

    curves = {
        "parabola": parabola_curve,
        "hyperbola": hyperbola_curve
    }

    # Create 3D revolution solids
    print("Creating 3D revolution solids...")
    paraboloid = create_paraboloid(
        focal_length=args.focal_length,
        z_max=4.0
    )
    hyperboloid_1 = create_hyperboloid_one_sheet(
        a=args.semi_major,
        b=args.semi_minor,
        z_range=(-3.0, 3.0)
    )

    solids = {
        "paraboloid": paraboloid,
        "hyperboloid_one_sheet": hyperboloid_1,
    }

    # Execute based on mode
    if args.mode in ("dxf", "all"):
        output_dir.mkdir(parents=True, exist_ok=True)
        dxf_path = output_dir / "conic_curves.dxf"
        export_curves_to_dxf(dxf_path, curves)

    if args.mode in ("step", "all"):
        export_solids_to_step(output_dir, solids)

    if args.mode in ("view", "all"):
        view_solids(solids)

    print("\nDemo complete!")
    print(f"  Parabola: vertex at origin, focal length = {args.focal_length}")
    print(f"  Hyperbola: center at origin, a = {args.semi_major}, b = {args.semi_minor}")
    print(f"  Paraboloid: r(z) = sqrt(4 * {args.focal_length} * z)")
    print(f"  Hyperboloid of one sheet: r(z) = {args.semi_major} * sqrt(1 + z^2/{args.semi_minor}^2)")


if __name__ == "__main__":
    main()
