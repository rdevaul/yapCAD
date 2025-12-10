#!/usr/bin/env python3
"""
CLI for yapCAD DSL compiler and runner.

Usage:
    python -m yapcad.dsl run FILE.dsl COMMAND [--param NAME=VALUE ...] [--output FILE]
    python -m yapcad.dsl check FILE.dsl
    python -m yapcad.dsl list FILE.dsl

TODO: Support 2D region export
  Currently the CLI assumes solid (3D) output and exports to STEP/STL/DXF
  as 3D geometry. Commands that return region2d should be exportable as:
  - DXF: Native 2D entities (LINE, ARC, CIRCLE, ELLIPSE, SPLINE)
  - SVG: For web/documentation use
  Detection: Check if result type is region2d and route to 2D exporter.
  See: yapcad.ezdxf_exporter for existing 2D DXF support

Examples:
    # Check syntax and types
    python -m yapcad.dsl check examples/spur_gears.dsl

    # List available commands
    python -m yapcad.dsl list examples/spur_gears.dsl

    # Run a command with parameters
    python -m yapcad.dsl run examples/spur_gears.dsl MAKE_GEAR_PAIR \
        --param teeth1=24 --param teeth2=36 --param module_mm=2.0 \
        --param pressure_angle=20.0 --param face_width=10.0 \
        --param hub_diameter=15.0 --param hub_height=5.0

    # Run and export to STEP file
    python -m yapcad.dsl run examples/spur_gears.dsl MAKE_GEAR_PAIR \
        --param teeth1=24 --param teeth2=36 --param module_mm=2.0 \
        --param pressure_angle=20.0 --param face_width=10.0 \
        --param hub_diameter=15.0 --param hub_height=5.0 \
        --output gears.step

    # Run and create a package
    python -m yapcad.dsl run examples/spur_gears.dsl MAKE_GEAR_PAIR \
        --param teeth1=24 --param teeth2=36 --param module_mm=2.0 \
        --param pressure_angle=20.0 --param face_width=10.0 \
        --param hub_diameter=15.0 --param hub_height=5.0 \
        --package output_pkg --name "gear_pair" --version "1.0.0"
"""

import argparse
import os
import sys
from pathlib import Path
from typing import Any, Dict


def format_type(type_node) -> str:
    """Format a TypeNode for display."""
    from .ast import SimpleType, GenericType, OptionalType
    if isinstance(type_node, SimpleType):
        return type_node.name
    elif isinstance(type_node, GenericType):
        args = ", ".join(format_type(a) for a in type_node.type_args)
        return f"{type_node.name}<{args}>"
    elif isinstance(type_node, OptionalType):
        return f"{format_type(type_node.inner)}?"
    else:
        return str(type_node)


def parse_param(param_str: str) -> tuple:
    """Parse a parameter string like 'name=value' into (name, typed_value)."""
    if '=' not in param_str:
        raise ValueError(f"Invalid parameter format: {param_str} (expected name=value)")

    name, value_str = param_str.split('=', 1)
    name = name.strip()
    value_str = value_str.strip()

    # Try to parse as int, float, bool, or string
    if value_str.lower() == 'true':
        return (name, True)
    elif value_str.lower() == 'false':
        return (name, False)

    try:
        return (name, int(value_str))
    except ValueError:
        pass

    try:
        return (name, float(value_str))
    except ValueError:
        pass

    # Strip quotes if present
    if (value_str.startswith('"') and value_str.endswith('"')) or \
       (value_str.startswith("'") and value_str.endswith("'")):
        value_str = value_str[1:-1]

    return (name, value_str)


def cmd_check(args):
    """Check a DSL file for syntax and type errors."""
    from . import tokenize, parse, check

    source_path = Path(args.file)
    if not source_path.exists():
        print(f"Error: File not found: {source_path}", file=sys.stderr)
        return 1

    source = source_path.read_text()

    try:
        tokens = tokenize(source)
        module = parse(tokens, source=source)  # Pass source for native block extraction
        result = check(module)

        if result.has_errors:
            print(f"Type checking failed with {len(result.diagnostics)} error(s):")
            for diag in result.diagnostics:
                print(f"  Line {diag.span.start.line}: {diag.message}")
            return 1

        print(f"OK: {source_path.name} - {len(module.commands)} command(s), no errors")
        if result.has_warnings:
            print(f"  {sum(1 for d in result.diagnostics if d.severity.value == 'warning')} warning(s)")

        return 0

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


def cmd_list(args):
    """List commands in a DSL file."""
    from . import tokenize, parse, check

    source_path = Path(args.file)
    if not source_path.exists():
        print(f"Error: File not found: {source_path}", file=sys.stderr)
        return 1

    source = source_path.read_text()

    try:
        tokens = tokenize(source)
        module = parse(tokens, source=source)  # Pass source for native block extraction

        print(f"Module: {module.name}")
        print(f"Commands ({len(module.commands)}):")

        for cmd in module.commands:
            params = ", ".join(
                f"{p.name}: {format_type(p.type_annotation)}" + (f" = {p.default_value}" if p.default_value else "")
                for p in cmd.parameters
            )
            print(f"  {cmd.name}({params}) -> {format_type(cmd.return_type)}")

        return 0

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


def cmd_run(args):
    """Run a command from a DSL file."""
    from . import compile_and_run
    from .packaging import package_from_dsl

    source_path = Path(args.file)
    if not source_path.exists():
        print(f"Error: File not found: {source_path}", file=sys.stderr)
        return 1

    source = source_path.read_text()

    # Parse parameters
    parameters: Dict[str, Any] = {}
    for param_str in args.param or []:
        try:
            name, value = parse_param(param_str)
            parameters[name] = value
        except ValueError as e:
            print(f"Error: {e}", file=sys.stderr)
            return 1

    print(f"Executing {args.command} with parameters: {parameters}")

    # Check if creating a package
    if args.package:
        result = package_from_dsl(
            source,
            args.command,
            parameters,
            args.package,
            name=args.name or args.command.lower(),
            version=args.version or "1.0.0",
            description=args.description,
            overwrite=args.force,
        )

        if not result.success:
            print(f"Error: {result.error_message}", file=sys.stderr)
            return 1

        print(f"Package created at: {result.manifest.root}")
        return 0

    # Regular execution
    result = compile_and_run(source, args.command, parameters)

    if not result.success:
        print(f"Error: {result.error_message}", file=sys.stderr)
        return 1

    print("Execution successful!")

    # Check geometry type
    from yapcad.geom3d import issolid, volumeof

    geometry = result.geometry
    if issolid(geometry):
        try:
            vol = volumeof(geometry)
            print(f"Result: solid with volume {vol:.2f}")
        except ValueError:
            print("Result: solid (volume calculation not available)")
    elif isinstance(geometry, list) and geometry and issolid(geometry[0]):
        try:
            total_vol = sum(volumeof(g) for g in geometry if issolid(g))
            print(f"Result: {len(geometry)} solid(s) with total volume {total_vol:.2f}")
        except ValueError:
            print(f"Result: {len(geometry)} solid(s)")
    else:
        print(f"Result: {type(geometry).__name__}")

    # Export if output specified
    if args.output:
        output_path = Path(args.output)
        suffix = output_path.suffix.lower()

        # Ensure we have a solid to export
        if issolid(geometry):
            solid = geometry
        elif isinstance(geometry, list) and geometry and issolid(geometry[0]):
            # For multiple solids, union them together
            if len(geometry) == 1:
                solid = geometry[0]
            else:
                from yapcad.geom3d import solid_boolean
                solid = geometry[0]
                for i in range(1, len(geometry)):
                    solid = solid_boolean(solid, geometry[i], 'union')
        else:
            print(f"Warning: Cannot export non-solid geometry", file=sys.stderr)
            return 0

        if suffix == '.step' or suffix == '.stp':
            step_format = os.environ.get('YAPCAD_STEP_FORMAT', 'faceted').lower()
            if step_format == 'analytic':
                from yapcad.io.step import write_step_analytic
                analytic_ok = write_step_analytic(solid, str(output_path))
                if analytic_ok:
                    print(f"Exported to: {output_path} (analytic BREP)")
                else:
                    print(f"Exported to: {output_path} (faceted fallback)")
            else:
                from yapcad.io import write_step
                write_step(solid, str(output_path))
                print(f"Exported to: {output_path}")
        elif suffix == '.stl':
            from yapcad.io import write_stl
            write_stl(solid, str(output_path))
            print(f"Exported to: {output_path}")
        elif suffix == '.dxf':
            from yapcad.ezdxf_exporter import write_dxf
            write_dxf(solid, str(output_path))
            print(f"Exported to: {output_path}")
        else:
            print(f"Warning: Unknown output format: {suffix}", file=sys.stderr)

    return 0


def main():
    parser = argparse.ArgumentParser(
        prog='python -m yapcad.dsl',
        description='yapCAD DSL compiler and runner',
    )

    subparsers = parser.add_subparsers(dest='action', required=True)

    # check command
    check_parser = subparsers.add_parser('check', help='Check DSL file for errors')
    check_parser.add_argument('file', help='DSL source file')

    # list command
    list_parser = subparsers.add_parser('list', help='List commands in DSL file')
    list_parser.add_argument('file', help='DSL source file')

    # run command
    run_parser = subparsers.add_parser('run', help='Run a command from DSL file')
    run_parser.add_argument('file', help='DSL source file')
    run_parser.add_argument('command', help='Command name to execute')
    run_parser.add_argument('-p', '--param', action='append', metavar='NAME=VALUE',
                          help='Parameter value (can be repeated)')
    run_parser.add_argument('-o', '--output', metavar='FILE',
                          help='Output file (STEP, STL, or DXF)')
    run_parser.add_argument('--package', metavar='DIR',
                          help='Create a yapCAD package instead of raw export')
    run_parser.add_argument('--name', help='Package name (for --package)')
    run_parser.add_argument('--version', default='1.0.0', help='Package version')
    run_parser.add_argument('--description', help='Package description')
    run_parser.add_argument('-f', '--force', action='store_true',
                          help='Overwrite existing output')

    args = parser.parse_args()

    if args.action == 'check':
        return cmd_check(args)
    elif args.action == 'list':
        return cmd_list(args)
    elif args.action == 'run':
        return cmd_run(args)
    else:
        parser.print_help()
        return 1


if __name__ == '__main__':
    sys.exit(main())
