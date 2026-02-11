#!/usr/bin/env python3
"""
Render a yapCAD assembly from DSL as multi-view PNG.

Adapted from Jeremy's assembly visualization scripts for the Roswell project.
Generalized to work with any yapCAD DSL module that has per-layer commands.

Usage:
    python render_assembly.py <dsl_file> [--commands CMD1,CMD2,...] [--output FILE]
                              [--style solid|wireframe|xray] [--views 4|6]

The script:
1. Runs each layer command to get geometry
2. Converts solids to trimesh for rendering
3. Generates multi-view matplotlib render
4. Saves PNG for inspection

Requires: conda yapcad-brep environment with trimesh + matplotlib
"""

import sys
import json
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Headless rendering
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from pathlib import Path
from typing import List, Dict, Tuple, Optional

sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from yapcad.dsl.runtime import compile_and_run
from yapcad.geom3d import issolid

# Try importing trimesh-based conversion
try:
    import trimesh
    HAS_TRIMESH = True
except ImportError:
    HAS_TRIMESH = False


# Default layer colors (RGBA)
LAYER_COLORS = {
    'body':      ((0.4, 0.5, 0.7, 0.7),  '#6688aa'),
    'connector': ((0.8, 0.5, 0.2, 0.8),  '#cc8844'),
    'arm':       ((0.6, 0.6, 0.6, 0.7),  '#999999'),
    'effector':  ((0.3, 0.7, 0.4, 0.8),  '#44aa66'),
    'actuator':  ((0.8, 0.2, 0.2, 0.8),  '#cc3333'),
    'power':     ((0.7, 0.6, 0.2, 0.7),  '#bb9933'),
    'default':   ((0.5, 0.5, 0.5, 0.6),  '#888888'),
}

# Color cycle for unknown layers
COLOR_CYCLE = [
    (0.2, 0.4, 0.8, 0.7),
    (0.8, 0.4, 0.2, 0.7),
    (0.3, 0.7, 0.3, 0.7),
    (0.7, 0.3, 0.7, 0.7),
    (0.2, 0.7, 0.7, 0.7),
    (0.7, 0.7, 0.2, 0.7),
]


def get_layer_color(layer_name: str, idx: int = 0) -> Tuple[float, float, float, float]:
    """Get RGBA color for a layer."""
    if layer_name in LAYER_COLORS:
        return LAYER_COLORS[layer_name][0]
    return COLOR_CYCLE[idx % len(COLOR_CYCLE)]


def solid_to_trimesh(solid) -> Optional[trimesh.Trimesh]:
    """Convert a yapCAD solid to a trimesh object.
    
    Uses the trimesh engine's internal converter which extracts
    triangulated surfaces from the solid's native representation.
    Falls back to BRep tessellation if available.
    """
    # Primary: use trimesh engine's solid-to-mesh converter
    try:
        from yapcad.boolean.trimesh_engine import _solid_to_mesh
        mesh = _solid_to_mesh(solid)
        if mesh is not None and len(mesh.vertices) > 0:
            return mesh
    except Exception as e:
        print(f"    [trimesh convert failed: {e}]")
    
    # Fallback: BRep tessellation via OCC
    try:
        from yapcad.brep import BrepSolid
        brep = BrepSolid.from_yapcad(solid)
        verts, faces = brep.tessellate(deflection=0.5)
        if len(verts) > 0 and len(faces) > 0:
            return trimesh.Trimesh(vertices=np.array(verts), faces=np.array(faces))
    except Exception:
        pass
    
    return None


def render_meshes_multiview(
    meshes: List[Dict],
    output_path: str,
    title: str = "Assembly Render",
    style: str = "xray",
    num_views: int = 4,
    dpi: int = 150,
):
    """Render meshes as multi-view PNG.
    
    Args:
        meshes: List of {'name': str, 'mesh': trimesh.Trimesh, 'color': RGBA tuple, 'layer': str}
        output_path: Path to save PNG
        title: Figure title
        style: 'solid', 'wireframe', or 'xray'
        num_views: 4 or 6
        dpi: Output resolution
    """
    # Calculate bounding box
    all_vertices = np.vstack([m['mesh'].vertices for m in meshes])
    bounds_min = all_vertices.min(axis=0)
    bounds_max = all_vertices.max(axis=0)
    bounds_center = (bounds_min + bounds_max) / 2
    bounds_size = (bounds_max - bounds_min).max()
    
    if bounds_size == 0:
        bounds_size = 1.0

    # Define views
    if num_views == 6:
        fig = plt.figure(figsize=(24, 16), facecolor='white')
        views = [
            (0, -90, 'Front (Y-)'),
            (0, 90, 'Back (Y+)'),
            (0, 0, 'Left (X+)'),
            (0, 180, 'Right (X-)'),
            (90, -90, 'Top (Z+)'),
            (30, -60, 'Isometric'),
        ]
        rows, cols = 2, 3
    else:
        fig = plt.figure(figsize=(16, 16), facecolor='white')
        views = [
            (30, -60, 'Perspective'),
            (0, -90, 'Front (Y-)'),
            (90, -90, 'Top (Z+)'),
            (0, 0, 'Right (X+)'),
        ]
        rows, cols = 2, 2

    fig.suptitle(title, fontsize=16, fontweight='bold', y=0.98)

    for idx, (elev, azim, vtitle) in enumerate(views, 1):
        ax = fig.add_subplot(rows, cols, idx, projection='3d')

        for mesh_data in meshes:
            mesh = mesh_data['mesh']
            color = mesh_data['color']
            vertices = mesh.vertices
            faces = mesh.faces

            if style == 'wireframe':
                edges = mesh.edges_unique
                segments = [[vertices[e[0]], vertices[e[1]]] for e in edges]
                line_color = color[:3]  # RGB only
                lc = Line3DCollection(segments, colors=[line_color], linewidths=0.3, alpha=0.7)
                ax.add_collection3d(lc)

            elif style == 'xray':
                polys = Poly3DCollection(
                    vertices[faces],
                    alpha=color[3] * 0.4,
                    edgecolor='black',
                    linewidths=0.1,
                )
                polys.set_facecolor(color)
                ax.add_collection3d(polys)

            else:  # solid
                polys = Poly3DCollection(
                    vertices[faces],
                    alpha=color[3],
                    edgecolor=(0, 0, 0, 0.05),
                    linewidths=0.05,
                )
                polys.set_facecolor(color)
                ax.add_collection3d(polys)

        # Configure view
        ax.view_init(elev=elev, azim=azim)
        ax.set_title(vtitle, fontsize=11, fontweight='bold', pad=10)
        ax.set_box_aspect([1, 1, 1])
        ax.set_xlabel('X', fontsize=8)
        ax.set_ylabel('Y', fontsize=8)
        ax.set_zlabel('Z', fontsize=8)
        ax.grid(True, alpha=0.2, linestyle='--')
        ax.tick_params(labelsize=7)

        if style == 'xray':
            ax.set_facecolor('#f8f8f8')

        # Set limits
        margin = bounds_size * 0.15
        for setter, center_val in [
            (ax.set_xlim, bounds_center[0]),
            (ax.set_ylim, bounds_center[1]),
            (ax.set_zlim, bounds_center[2]),
        ]:
            setter([center_val - bounds_size/2 - margin,
                    center_val + bounds_size/2 + margin])

        # Coordinate axes
        scale = bounds_size * 0.1
        o = bounds_min - margin * 0.5
        ax.plot([o[0], o[0]+scale], [o[1], o[1]], [o[2], o[2]], 'r-', lw=2)
        ax.plot([o[0], o[0]], [o[1], o[1]+scale], [o[2], o[2]], 'g-', lw=2)
        ax.plot([o[0], o[0]], [o[1], o[1]], [o[2], o[2]+scale], 'b-', lw=2)

    # Layer legend
    from matplotlib.patches import Patch
    seen_layers = {}
    for m in meshes:
        layer = m.get('layer', 'default')
        if layer not in seen_layers:
            seen_layers[layer] = m['color']
    
    legend_elements = [
        Patch(facecolor=color[:3], alpha=color[3], label=layer)
        for layer, color in seen_layers.items()
    ]
    if legend_elements:
        fig.legend(handles=legend_elements, loc='lower center',
                   bbox_to_anchor=(0.5, 0.0), ncol=min(8, len(legend_elements)),
                   fontsize=10, frameon=True)

    plt.tight_layout(rect=[0, 0.03, 1, 0.96])
    
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close()
    
    return output_path


def run_dsl_layers(dsl_path: str, commands: List[Tuple[str, str]]) -> List[Dict]:
    """Run DSL commands and convert results to renderable meshes.
    
    Args:
        dsl_path: Path to DSL file
        commands: List of (command_name, layer_name) tuples
    
    Returns:
        List of mesh dicts for rendering
    """
    source = Path(dsl_path).read_text()
    meshes = []
    
    for idx, (cmd, layer) in enumerate(commands):
        print(f"  [{idx+1}/{len(commands)}] Running {cmd}...", end=' ', flush=True)
        
        result = compile_and_run(source, cmd, {})
        if not result.success:
            print(f"FAILED: {result.error_message}")
            continue
        
        geom = result.geometry
        meta = result.metadata
        layer_name = meta.get('layer', layer)
        
        # Convert to trimesh
        mesh = solid_to_trimesh(geom)
        
        if mesh is None:
            print(f"SKIP (no mesh)")
            continue
        
        color = get_layer_color(layer_name, idx)
        
        meshes.append({
            'name': meta.get('name', cmd),
            'mesh': mesh,
            'color': color,
            'layer': layer_name,
        })
        
        print(f"OK ({len(mesh.vertices)} verts, {len(mesh.faces)} faces)")
    
    return meshes


def main():
    parser = argparse.ArgumentParser(description='Render yapCAD assembly multi-view')
    parser.add_argument('dsl_file', help='Path to DSL source file')
    parser.add_argument('--commands', '-c', default=None,
                        help='Comma-separated COMMAND:layer pairs (e.g., BODY:body,ARMS:arm)')
    parser.add_argument('--output', '-o', default=None,
                        help='Output PNG path (default: <dsl_name>_render.png)')
    parser.add_argument('--style', '-s', choices=['solid', 'wireframe', 'xray'],
                        default='xray', help='Render style')
    parser.add_argument('--views', '-v', type=int, choices=[4, 6],
                        default=4, help='Number of views')
    parser.add_argument('--dpi', type=int, default=150, help='Output DPI')
    parser.add_argument('--title', '-t', default=None, help='Figure title')
    
    args = parser.parse_args()
    
    dsl_path = Path(args.dsl_file)
    if not dsl_path.exists():
        print(f"ERROR: DSL file not found: {dsl_path}")
        return 1
    
    # Parse commands
    if args.commands:
        commands = []
        for pair in args.commands.split(','):
            parts = pair.strip().split(':')
            cmd = parts[0]
            layer = parts[1] if len(parts) > 1 else cmd.lower()
            commands.append((cmd, layer))
    else:
        # Auto-detect: look for uppercase commands in the DSL source
        source = dsl_path.read_text()
        import re
        all_cmds = re.findall(r'^command\s+([A-Z][A-Z_0-9]*)\s*\(', source, re.MULTILINE)
        if not all_cmds:
            print("ERROR: No uppercase commands found. Use --commands to specify.")
            return 1
        # Filter out ASSEMBLY (it's the combined one)
        commands = [(cmd, cmd.lower()) for cmd in all_cmds if cmd != 'ASSEMBLY']
        if not commands:
            commands = [(all_cmds[0], all_cmds[0].lower())]
        print(f"Auto-detected commands: {[c[0] for c in commands]}")
    
    # Output path
    output_path = args.output or f"/tmp/{dsl_path.stem}_render.png"
    title = args.title or f"{dsl_path.stem} Assembly"
    
    print(f"\n{'='*60}")
    print(f"yapCAD Assembly Renderer")
    print(f"{'='*60}")
    print(f"DSL:      {dsl_path}")
    print(f"Commands: {[c[0] for c in commands]}")
    print(f"Style:    {args.style}")
    print(f"Views:    {args.views}")
    print(f"Output:   {output_path}")
    print(f"{'='*60}\n")
    
    # Run DSL and get meshes
    print("[1/2] Building geometry from DSL...")
    meshes = run_dsl_layers(str(dsl_path), commands)
    
    if not meshes:
        print("\nERROR: No renderable geometry produced!")
        return 1
    
    print(f"\n  Total: {len(meshes)} layers, "
          f"{sum(len(m['mesh'].vertices) for m in meshes)} vertices, "
          f"{sum(len(m['mesh'].faces) for m in meshes)} faces")
    
    # Render
    print(f"\n[2/2] Rendering {args.views}-view {args.style}...")
    render_meshes_multiview(
        meshes, output_path, title=title,
        style=args.style, num_views=args.views, dpi=args.dpi,
    )
    
    print(f"\n{'='*60}")
    print(f"RENDER COMPLETE: {output_path}")
    print(f"{'='*60}")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
