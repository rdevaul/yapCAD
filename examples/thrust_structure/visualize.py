#!/usr/bin/env python3
"""
Visualization and export utility for thrust structure designs.

Usage:
    python visualize.py [OPTIONS]

Commands:
    python visualize.py                      # View baseline design
    python visualize.py --light              # View 9-hole design
    python visualize.py --optimized          # View 12-hole design
    python visualize.py --export STEP|STL    # Export to file
    python visualize.py --info               # Show design info only

Options:
    --baseline      View baseline design (no holes)
    --light         View light design (9 holes)
    --optimized     View optimized design (12 holes)
    --export FMT    Export to format (STEP, STL)
    --output PATH   Output file path
    --info          Print design info without viewing
    --no-viewer     Skip opening viewer
"""

import argparse
import sys
from pathlib import Path

# Add src to path for development
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))


DESIGNS = {
    "baseline": {
        "name": "Baseline (no holes)",
        "command": "MAKE_BASELINE_PLATE",
        "num_holes": 0,
    },
    "light": {
        "name": "Light (9 holes)",
        "command": "MAKE_LIGHTENED_PLATE_3",
        "num_holes": 9,
    },
    "optimized": {
        "name": "Optimized (12 holes)",
        "command": "MAKE_OPTIMIZED_PLATE",
        "num_holes": 12,
    },
}


def generate_thrust_plate(design_key: str):
    """Generate thrust plate geometry for the specified design."""
    from yapcad.dsl import compile_and_run

    design = DESIGNS[design_key]

    # Read DSL source
    dsl_path = Path(__file__).parent / "thrust_structure.dsl"
    source = dsl_path.read_text()

    result = compile_and_run(source, design["command"], {})

    if not result.success:
        raise RuntimeError(f"Failed to generate geometry: {result.error_message}")

    return result.geometry, design


def calculate_mass(solid):
    """Calculate mass of solid (assuming 6061-T6 aluminum)."""
    from yapcad.geom3d import volumeof

    volume_mm3 = volumeof(solid)
    volume_m3 = volume_mm3 / 1e9
    mass_kg = volume_m3 * 2700  # 6061-T6 density
    return mass_kg, volume_mm3


def print_design_info(design: dict, mass_kg: float, volume_mm3: float):
    """Print design information."""
    print("\n" + "=" * 50)
    print("THRUST STRUCTURE DESIGN")
    print("=" * 50)
    print(f"\nDesign: {design['name']}")
    print(f"\nFixed Parameters:")
    print(f"  Outer diameter:     304.8 mm (12 in)")
    print(f"  Thickness:          8.0 mm")
    print(f"  Motor mount:        101.6 mm (4 in)")
    print(f"  Stringer notches:   3 @ 120° intervals")
    print(f"  Notch size:         50.8 x 25.4 mm (2x1 in)")

    print(f"\nLightening Holes:")
    print(f"  Total holes:        {design['num_holes']}")

    print(f"\nResults:")
    print(f"  Volume:             {volume_mm3:.0f} mm³")
    print(f"  Mass:               {mass_kg:.3f} kg ({mass_kg * 2.205:.3f} lbs)")
    print(f"  Material:           6061-T6 Aluminum")
    print("=" * 50)


def export_geometry(solid, output_path: Path, format: str):
    """Export geometry to specified format."""
    format = format.upper()

    if format == "STEP" or format == "STP":
        from yapcad.io import write_step
        write_step(solid, str(output_path))
    elif format == "STL":
        from yapcad.io import write_stl
        write_stl(solid, str(output_path))
    else:
        raise ValueError(f"Unknown format: {format}")

    print(f"Exported to: {output_path}")


def view_geometry(solid, design_name: str = "Thrust Structure"):
    """Open the yapCAD viewer to display the geometry."""
    import tempfile

    try:
        from yapcad.package import view_package, create_package_from_entities

        # Create a temporary package for viewing
        with tempfile.TemporaryDirectory() as tmpdir:
            pkg_path = Path(tmpdir) / "temp_view.ycpkg"

            # Create package with the solid
            create_package_from_entities(
                [solid],  # entities as a list
                pkg_path,
                name=design_name,
                version="1.0.0",
                description="Thrust structure visualization",
                units="mm",
            )

            print("\nOpening viewer...")
            print("  Controls:")
            print("    - Click and drag to rotate")
            print("    - Scroll to zoom")
            print("    - Press 'H' for help")
            print("    - Press 'Q' or close window to exit")

            view_package(pkg_path)

    except ImportError as e:
        print(f"\nViewer not available: {e}")
        print("Export to STEP/STL and view in external application.")
        print("  Example: python visualize.py --optimized --export STEP --output plate.step")
    except Exception as e:
        import traceback
        print(f"\nViewer error: {e}")
        traceback.print_exc()
        print("This may occur in headless environments or if pyglet is not installed.")


def main():
    parser = argparse.ArgumentParser(
        description="Visualize and export thrust structure designs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python visualize.py                          # View baseline
    python visualize.py --light                  # View 9-hole design
    python visualize.py --optimized              # View 12-hole design
    python visualize.py --export STEP --output plate.step
        """
    )

    # Design selection
    design_group = parser.add_mutually_exclusive_group()
    design_group.add_argument(
        "--baseline", action="store_true",
        help="Use baseline configuration (no holes)"
    )
    design_group.add_argument(
        "--light", action="store_true",
        help="Use light configuration (9 holes)"
    )
    design_group.add_argument(
        "--optimized", action="store_true",
        help="Use optimized configuration (12 holes)"
    )

    # Export options
    parser.add_argument(
        "--export", type=str, choices=["STEP", "STL", "step", "stl"],
        help="Export format"
    )
    parser.add_argument(
        "--output", type=Path,
        help="Output file path (default: thrust_plate.{format})"
    )

    # Display options
    parser.add_argument(
        "--info", action="store_true",
        help="Print design info only (no viewer)"
    )
    parser.add_argument(
        "--no-viewer", action="store_true",
        help="Skip opening viewer"
    )

    args = parser.parse_args()

    # Determine which design to use
    if args.light:
        design_key = "light"
    elif args.optimized:
        design_key = "optimized"
    else:
        design_key = "baseline"

    print(f"Using {DESIGNS[design_key]['name']} design")

    # Generate geometry
    print("\nGenerating geometry...")
    solid, design = generate_thrust_plate(design_key)

    # Calculate mass
    mass_kg, volume_mm3 = calculate_mass(solid)

    # Print info
    print_design_info(design, mass_kg, volume_mm3)

    # Export if requested
    if args.export:
        output_path = args.output
        if output_path is None:
            ext = args.export.lower()
            output_path = Path(f"thrust_plate_{design_key}.{ext}")
        export_geometry(solid, output_path, args.export)

    # View if not suppressed
    if not args.info and not args.no_viewer:
        view_geometry(solid, design["name"])

    return 0


if __name__ == "__main__":
    sys.exit(main())
