#!/usr/bin/env python3
"""
Render a yapCAD assembly from DSL as multi-view PNG.

A thin wrapper over :mod:`yapcad.render`: each layer command is built as
an SDF field (meshing deferred -- the renderer sphere-traces fields and
never looks at triangles), every part keeps its own colour, and the sheet
comes from ``yapcad.render.four_up`` or ``yapcad.render.image``.

Usage:
    python tools/render_assembly.py <dsl_file> [--commands CMD1,CMD2,...]
                                    [--output FILE] [--views 4|6]
                                    [--view iso|front|top|...] [--title T]

The old matplotlib styles (solid|wireframe|xray) and --dpi are gone; the
renderer shades fields and does not rasterise triangle soup.

The script:
1. Runs each layer command as a deferred SDF build
2. Traces every part against its own field
3. Saves one PNG (4-up sheet, or one image per named view with --views 6)
"""

import argparse
import sys
from pathlib import Path
from typing import List, Tuple

sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from yapcad.render.api import (  # noqa: E402
    Part, as_parts, four_up, image, scene_bounds)
from yapcad.render.camera import VIEWS  # noqa: E402

SIX_VIEWS = ("iso", "front", "back", "left", "right", "top")


def parse_commands(spec: str) -> List[Tuple[str, str]]:
    """Parse a COMMAND:layer,... list into (command, layer) pairs."""
    commands = []
    for pair in spec.split(','):
        parts = pair.strip().split(':')
        if not parts or not parts[0]:
            continue
        cmd = parts[0]
        layer = parts[1] if len(parts) > 1 else cmd.lower()
        commands.append((cmd, layer))
    return commands


def detect_commands(source: str) -> List[Tuple[str, str]]:
    """Auto-detect uppercase commands in the DSL source."""
    import re
    all_cmds = re.findall(r'^command\s+([A-Z][A-Z_0-9]*)\s*\(',
                          source, re.MULTILINE)
    # ASSEMBLY is the combined model, not a layer
    commands = [(c, c.lower()) for c in all_cmds if c != 'ASSEMBLY']
    if not commands and all_cmds:
        commands = [(all_cmds[0], all_cmds[0].lower())]
    return commands


def build_parts(dsl_path: str,
                commands: List[Tuple[str, str]]) -> List[dict]:
    """Run each layer command as a deferred SDF build.

    Returns one record per command: the field wrapped as render
    :class:`~yapcad.render.api.Part` plus its layer name, ready to trace.
    """
    from yapcad.dsl import compile_and_run

    source = Path(dsl_path).read_text()
    parts = []
    for idx, (cmd, layer) in enumerate(commands):
        print(f"  [{idx+1}/{len(commands)}] Building {cmd}...",
              end=' ', flush=True)
        result = compile_and_run(source, cmd, {}, representation="sdf",
                                 sdf_mesh=False, source_path=str(dsl_path))
        if not result.success:
            print(f"FAILED: {result.error_message}")
            continue
        field = as_parts(result.geometry)
        if not field:
            print("SKIP (no SDF field)")
            continue
        part = Part(field[0].node, name=cmd, group=layer)
        parts.append({'part': part, 'layer': layer, 'name': cmd})
        lo, hi = part.node.bounds
        dims = ', '.join(f"{hi[i] - lo[i]:.1f}" for i in range(3))
        print(f"OK ({dims} mm)")
    return parts


def main():
    parser = argparse.ArgumentParser(
        description='Render a yapCAD assembly multi-view (SDF fields)')
    parser.add_argument('dsl_file', help='Path to DSL source file')
    parser.add_argument('--commands', '-c', default=None,
                        help='Comma-separated COMMAND:layer pairs '
                             '(e.g., BODY:body,ARMS:arm)')
    parser.add_argument('--output', '-o', default=None,
                        help='Output PNG path (default: <dsl>_render.png)')
    parser.add_argument('--views', '-v', type=int, choices=[4, 6],
                        default=4, help='4 = one 4-up sheet, '
                                        '6 = one image per named view')
    parser.add_argument('--view', default=None,
                        help='Named view for --views 6; the default draws '
                             'every view (one of '
                             + ', '.join(sorted(VIEWS)) + ')')
    parser.add_argument('--theme', default='light',
                        choices=('light', 'dark'))
    parser.add_argument('--title', '-t', default=None, help='Figure title')

    args = parser.parse_args()

    dsl_path = Path(args.dsl_file)
    if not dsl_path.exists():
        print(f"ERROR: DSL file not found: {dsl_path}", file=sys.stderr)
        return 1

    if args.commands:
        commands = parse_commands(args.commands)
    else:
        commands = detect_commands(dsl_path.read_text())
        if not commands:
            print("ERROR: No uppercase commands found. "
                  "Use --commands to specify.", file=sys.stderr)
            return 1
        print(f"Auto-detected commands: {[c[0] for c in commands]}")

    title = args.title or f"{dsl_path.stem} Assembly"

    print(f"\n[1/2] Building {len(commands)} layer(s) as deferred fields...")
    records = build_parts(str(dsl_path), commands)
    if not records:
        print("\nERROR: No renderable geometry produced!", file=sys.stderr)
        return 1

    parts = [r['part'] for r in records]
    print(f"\n  Total: {len(parts)} part(s) traced from fields, "
          f"bounds {scene_bounds(parts)}")

    print(f"\n[2/2] Rendering "
          f"{'4-up sheet' if args.views == 4 else '6 views'}...")
    if args.views == 4:
        out = args.output or f"/tmp/{dsl_path.stem}_render.png"
        four_up(parts, theme=args.theme, colour_by="part", path=out)
        print(f"\nRENDER COMPLETE ({title}): {out}")
    else:
        views = [args.view] if args.view else SIX_VIEWS
        stem = args.output or f"/tmp/{dsl_path.stem}_render.png"
        base = Path(stem)
        for view in views:
            out = base.with_name(f"{base.stem}_{view}{base.suffix}")
            image(parts, view=view, title=title, theme=args.theme,
                  colour_by="part", path=str(out))
            print(f"RENDER COMPLETE: {out}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
