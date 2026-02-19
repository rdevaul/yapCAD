#!/usr/bin/env python3
"""Fast matplotlib multiview render using wireframes instead of full polygon rendering.

Generates a 2x3 grid showing:
- Front (looking along -Y)
- Back (looking along +Y)
- Left (looking along +X)
- Right (looking along -X)
- Top (looking along -Z)
- Isometric (3/4 view)

Uses kinematic chain positions from positions_home.json and loads all STL files.
Renders as wireframe/edge representation for faster rendering.
"""

import json
import os
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import trimesh

# Color scheme by subsystem (simplified for wireframe)
COLORS = {
    'chassis': '#3333CC',
    'wheel_arm': '#888888',
    'pivot': '#5555AA',
    'motor': '#1a1a1a',
    'gearbox_axis1': '#CC6633',
    'gearbox_axis2': '#DD7744',
    'gearbox_axis3': '#EE8855',
    'scara_link': '#33AA33',
    'scara_tower': '#4d4d7f',
    'battery': '#4B0082',
    'servo': '#1a1a1a',
    'default': '#666666'
}

def load_positions(json_path):
    """Load assembly positions from JSON file."""
    with open(json_path, 'r') as f:
        return json.load(f)

def find_stl_for_part(part_name, stl_base_dir, positions=None):
    """Find the STL file for a given part name.

    First checks positions JSON for stl_path (source of truth),
    then falls back to hardcoded mappings for legacy compatibility.
    """
    stl_base = Path(stl_base_dir)

    # First try: Use stl_path from positions JSON (source of truth)
    if positions and part_name in positions:
        stl_rel = positions[part_name].get('stl_path')
        if stl_rel and stl_rel != 'None':
            stl_path = stl_base / stl_rel
            if stl_path.exists():
                return stl_path

    # Fallback: hardcoded mappings for parts without stl_path in JSON
    mappings = {
        # Chassis (CENTRAL_HUB is unioned into chassis_with_battery_mounts.stl)
        'CHASSIS_PLATE': 'chassis/chassis_with_battery_mounts.stl',
        'ARM_TOWER_MOUNT': 'chassis/arm_tower_mount.stl',
        'PIVOT_BOSS_1': 'chassis/pivot_boss.stl',
        'PIVOT_BOSS_2': 'chassis/pivot_boss.stl',
        'PIVOT_BOSS_3': 'chassis/pivot_boss.stl',
    }

    if part_name in mappings:
        stl_path = stl_base / mappings[part_name]
        if stl_path.exists():
            return stl_path

    return None

def get_part_color(part_name):
    """Get color for a part based on subsystem."""
    if 'CHASSIS' in part_name:
        return COLORS['chassis']
    elif 'WHEEL_ARM' in part_name:
        return COLORS['wheel_arm']
    elif 'PIVOT_BOSS' in part_name:
        return COLORS['pivot']
    elif 'MOTOR' in part_name or 'DDSM115' in part_name:
        return COLORS['motor']
    elif 'BATTERY' in part_name:
        return COLORS['battery']
    elif 'AXIS1' in part_name:
        return COLORS['gearbox_axis1']
    elif 'AXIS2' in part_name:
        return COLORS['gearbox_axis2']
    elif 'AXIS3' in part_name or 'AXIS4' in part_name:
        return COLORS['gearbox_axis3']
    elif 'SCARA_LINK' in part_name or 'LINK' in part_name or 'EXTRUDER' in part_name:
        return COLORS['scara_link']
    elif 'ARM_TOWER' in part_name:
        return COLORS['scara_tower']
    elif 'SERVO' in part_name:
        return COLORS['servo']
    else:
        return COLORS['default']

def load_and_transform_mesh(stl_path, transform_matrix, simplify=True):
    """Load STL and apply transform. Optionally simplify for faster rendering."""
    mesh = trimesh.load(str(stl_path))

    # Simplify mesh to reduce face count for faster rendering
    if simplify and len(mesh.faces) > 1000:
        # Reduce to ~500 faces for complex meshes
        target_faces = min(500, len(mesh.faces) // 2)
        try:
            mesh = mesh.simplify_quadric_decimation(target_faces)
        except:
            pass  # Keep original if simplification fails

    # Apply transform
    mesh.apply_transform(transform_matrix)

    return mesh

def mesh_to_wireframe(mesh, color):
    """Convert trimesh to wireframe Line3DCollection (edges only)."""
    # Get unique edges
    edges_sparse = mesh.edges_unique
    vertices = mesh.vertices

    # Create line segments
    segments = [[vertices[edge[0]], vertices[edge[1]]] for edge in edges_sparse]

    return Line3DCollection(segments, colors=color, linewidths=0.3, alpha=0.7)

def setup_view(ax, elev, azim, title):
    """Set up a specific view of the assembly."""
    ax.view_init(elev=elev, azim=azim)
    ax.set_title(title, fontsize=12, fontweight='bold')

    # Set equal aspect ratio
    ax.set_box_aspect([1, 1, 1])

    # Labels
    ax.set_xlabel('X (mm)', fontsize=8)
    ax.set_ylabel('Y (mm)', fontsize=8)
    ax.set_zlabel('Z (mm)', fontsize=8)

    # Grid
    ax.grid(True, alpha=0.3)
    ax.tick_params(labelsize=7)

def main():
    # Paths
    script_dir = Path(__file__).parent
    positions_path = script_dir / 'output' / 'chain_assembly' / 'positions_home.json'
    stl_base_dir = script_dir / 'output' / 'stl'
    output_path = script_dir / 'output' / 'renders' / 'assembly_matplotlib_multiview.png'

    print("Loading assembly positions...")
    positions = load_positions(positions_path)

    print(f"Found {len(positions)} parts in kinematic chain")

    # Load all meshes
    print("\nLoading STL files (simplified for fast rendering)...")
    meshes = []
    part_count = 0
    missing_count = 0

    for part_name, part_data in positions.items():
        stl_path = find_stl_for_part(part_name, stl_base_dir, positions)

        if stl_path is None:
            missing_count += 1
            continue

        try:
            transform_matrix = np.array(part_data['transform_matrix'])
            mesh = load_and_transform_mesh(stl_path, transform_matrix, simplify=True)
            color = get_part_color(part_name)

            meshes.append({
                'name': part_name,
                'mesh': mesh,
                'color': color
            })

            part_count += 1
            if part_count % 5 == 0:
                print(f"  Loaded {part_count} parts...")

        except Exception as e:
            print(f"  [ERROR] {part_name}: {e}")
            missing_count += 1

    print(f"\nLoaded {part_count} parts, {missing_count} missing")

    if part_count == 0:
        print("ERROR: No parts loaded!")
        return 1

    # Calculate bounding box
    print("\nCalculating assembly bounds...")
    all_vertices = np.vstack([m['mesh'].vertices for m in meshes])
    bounds_min = all_vertices.min(axis=0)
    bounds_max = all_vertices.max(axis=0)
    bounds_center = (bounds_min + bounds_max) / 2
    bounds_size = (bounds_max - bounds_min).max()

    print(f"Assembly bounds:")
    print(f"  X: [{bounds_min[0]:.1f}, {bounds_max[0]:.1f}] mm")
    print(f"  Y: [{bounds_min[1]:.1f}, {bounds_max[1]:.1f}] mm")
    print(f"  Z: [{bounds_min[2]:.1f}, {bounds_max[2]:.1f}] mm")
    print(f"  Max dimension: {bounds_size:.1f} mm")

    # Create figure with 2x3 grid
    print("\nGenerating multiview wireframe render...")
    fig = plt.figure(figsize=(24, 16))
    fig.suptitle('Robot Assembly - Multiview Wireframe Render', fontsize=18, fontweight='bold')

    # Define views: (elev, azim, title)
    views = [
        (0, -90, 'Front View (looking along -Y)'),
        (0, 90, 'Back View (looking along +Y)'),
        (0, 0, 'Left View (looking along +X)'),
        (0, 180, 'Right View (looking along -X)'),
        (90, -90, 'Top View (looking along -Z)'),
        (30, -60, 'Isometric View (3/4)')
    ]

    # Create subplots
    for idx, (elev, azim, title) in enumerate(views, 1):
        print(f"  Rendering view {idx}/6: {title}")
        ax = fig.add_subplot(2, 3, idx, projection='3d')

        # Add all meshes as wireframes
        for mesh_data in meshes:
            wireframe = mesh_to_wireframe(mesh_data['mesh'], mesh_data['color'])
            ax.add_collection3d(wireframe)

        # Set view
        setup_view(ax, elev, azim, title)

        # Set limits
        margin = bounds_size * 0.1
        ax.set_xlim([bounds_center[0] - bounds_size/2 - margin,
                     bounds_center[0] + bounds_size/2 + margin])
        ax.set_ylim([bounds_center[1] - bounds_size/2 - margin,
                     bounds_center[1] + bounds_size/2 + margin])
        ax.set_zlim([bounds_center[2] - bounds_size/2 - margin,
                     bounds_center[2] + bounds_size/2 + margin])

    # Adjust layout
    plt.tight_layout(rect=[0, 0, 1, 0.97])

    # Save
    print(f"\nSaving to {output_path}...")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved: {output_path}")

    # Summary
    print("\n" + "="*70)
    print("ASSEMBLY MULTIVIEW RENDER COMPLETE")
    print("="*70)
    print(f"Parts included: {part_count}")
    print(f"Parts missing: {missing_count}")
    print(f"Output: {output_path}")
    print("\nKey subsystems (wireframe rendering):")
    print("  - Chassis (blue): Plate, hub, tower mount, pivot bosses")
    print("  - Wheel assemblies (gray/black): 3x wheel arms + DDSM115 motors")
    print("  - Gearboxes (orange): Axis 1/2/3 planetary gearboxes")
    print("  - SCARA arm (green): Links 1 & 2, extruder mount")
    print("  - Battery (purple): Battery cage")
    print("  - Servos (black): XH540/XH430 servos (COTS)")
    print("="*70)

    return 0

if __name__ == '__main__':
    sys.exit(main())
