#!/usr/bin/env python3
"""Matplotlib multiview render of the entire robot assembly.

Generates a 2x3 grid showing:
- Front (looking along -Y)
- Back (looking along +Y)
- Left (looking along +X)
- Right (looking along -X)
- Top (looking along -Z)
- Isometric (3/4 view)

Uses kinematic chain positions from positions_home.json and loads all STL files.
"""

import json
import os
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import trimesh

# Color scheme by subsystem
COLORS = {
    'chassis': '#3333CC',      # Blue
    'wheel_arm': '#888888',    # Gray
    'pivot': '#5555AA',        # Medium blue
    'motor': '#1a1a1a',        # Black
    'wheel': '#4d4d4d',        # Dark gray
    'gearbox_axis1': '#CC6633',  # Orange (Axis 1)
    'gearbox_axis2': '#DD7744',  # Lighter orange (Axis 2)
    'gearbox_axis3': '#EE8855',  # Even lighter orange (Axis 3)
    'scara_link': '#33AA33',   # Green
    'scara_tower': '#4d4d7f',  # Blue-purple
    'battery': '#4B0082',      # Indigo (purple)
    'servo': '#1a1a1a',        # Black (COTS)
    'default': '#666666'       # Medium gray
}

def load_positions(json_path):
    """Load assembly positions from JSON file."""
    with open(json_path, 'r') as f:
        return json.load(f)

def find_stl_for_part(part_name, stl_base_dir):
    """Find the STL file for a given part name.

    Args:
        part_name: Name from kinematic chain (e.g., "CHASSIS_PLATE")
        stl_base_dir: Base directory for STL files

    Returns:
        Path to STL file or None if not found
    """
    stl_base = Path(stl_base_dir)

    # Map part names to STL paths
    mappings = {
        # Chassis (CENTRAL_HUB is unioned into chassis_plate.stl)
        'CHASSIS_PLATE': 'chassis/chassis_plate.stl',
        'ARM_TOWER_MOUNT': 'chassis/arm_tower_mount.stl',
        'PIVOT_BOSS_1': 'chassis/pivot_boss.stl',
        'PIVOT_BOSS_2': 'chassis/pivot_boss.stl',
        'PIVOT_BOSS_3': 'chassis/pivot_boss.stl',

        # Wheel assemblies
        'WHEEL_ARM_1': 'drive/wheel_arm.stl',
        'WHEEL_ARM_2': 'drive/wheel_arm.stl',
        'WHEEL_ARM_3': 'drive/wheel_arm.stl',
        'DDSM115_MOTOR_1': 'cots/ddsm115_motor_official.stl',
        'DDSM115_MOTOR_2': 'cots/ddsm115_motor_official.stl',
        'DDSM115_MOTOR_3': 'cots/ddsm115_motor_official.stl',

        # SCARA arm tower
        'ARM_TOWER': 'scara/arm_tower.stl',

        # Axis 1 gearbox
        'AXIS1_RING_HOUSING': 'gearbox/axis1_ring_housing.stl',
        'AXIS1_SUN_GEAR': 'gearbox/axis1_sun_gear.stl',
        'AXIS1_PLANET_GEAR_1': 'gearbox/axis1_planet_gear.stl',
        'AXIS1_PLANET_GEAR_2': 'gearbox/axis1_planet_gear.stl',
        'AXIS1_PLANET_GEAR_3': 'gearbox/axis1_planet_gear.stl',
        'AXIS1_CARRIER_TOP': 'gearbox/axis1_carrier_top.stl',

        # Axis 2 gearbox
        'AXIS2_RING_HOUSING': 'gearbox/axis2_ring_housing.stl',
        'AXIS2_SUN_GEAR': 'gearbox/axis2_sun_gear.stl',
        'AXIS2_PLANET_GEAR_1': 'gearbox/axis2_planet_gear.stl',
        'AXIS2_PLANET_GEAR_2': 'gearbox/axis2_planet_gear.stl',
        'AXIS2_PLANET_GEAR_3': 'gearbox/axis2_planet_gear.stl',
        'AXIS2_CARRIER_TOP': 'gearbox/axis2_carrier_top.stl',

        # Axis 3 gearbox
        'AXIS3_RING_HOUSING': 'gearbox/axis3_ring_housing.stl',
        'AXIS3_SUN_GEAR': 'gearbox/axis3_sun_gear.stl',
        'AXIS3_PLANET_GEAR_1': 'gearbox/axis3_planet_gear.stl',
        'AXIS3_PLANET_GEAR_2': 'gearbox/axis3_planet_gear.stl',
        'AXIS3_PLANET_GEAR_3': 'gearbox/axis3_planet_gear.stl',
        'AXIS3_CARRIER_TOP': 'gearbox/axis3_carrier_top.stl',

        # SCARA links
        'LINK_2_3': 'scara/link1_attractive.stl',  # Use attractive version
        'LINK_3_4': 'scara/link2_attractive.stl',

        # Axis 4 (wrist)
        'AXIS4_GEARBOX': 'scara/axis4_housing.stl',
        'EXTRUDER_MOUNT': 'scara/axis4_output.stl',

        # Battery
        'BATTERY_CAGE': 'power/battery_cage.stl',

        # COTS servos (use xh430/xh540 as available)
        'AXIS1_SERVO_XH540': 'cots/xh540_cots.stl',
        'AXIS2_SERVO_XH430': 'cots/xh430_cots.stl',
        'AXIS3_SERVO_XH430': 'cots/xh430_cots.stl',
        'AXIS4_SERVO_XL330': 'cots/xh430_cots.stl',  # Use xh430 as placeholder
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

def load_and_transform_mesh(stl_path, transform_matrix):
    """Load STL and apply transform."""
    mesh = trimesh.load(str(stl_path))

    # Apply transform
    mesh.apply_transform(transform_matrix)

    return mesh

def mesh_to_poly3d(mesh, color, alpha=0.8):
    """Convert trimesh to matplotlib Poly3DCollection."""
    vertices = mesh.vertices
    faces = mesh.faces

    # Create polygon collection
    poly = [[vertices[face[i]] for i in range(3)] for face in faces]

    return Poly3DCollection(poly, alpha=alpha, facecolor=color, edgecolor='k', linewidths=0.05)

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
    project_root = script_dir.parent
    positions_path = script_dir / 'output' / 'chain_assembly' / 'positions_home.json'
    stl_base_dir = script_dir / 'output' / 'stl'
    output_path = script_dir / 'output' / 'renders' / 'assembly_matplotlib_multiview.png'

    print("Loading assembly positions...")
    positions = load_positions(positions_path)

    print(f"Found {len(positions)} parts in kinematic chain")

    # Load all meshes
    print("\nLoading STL files...")
    meshes = []
    part_count = 0
    missing_count = 0

    for part_name, part_data in positions.items():
        stl_path = find_stl_for_part(part_name, stl_base_dir)

        if stl_path is None:
            print(f"  [SKIP] {part_name}: STL not found")
            missing_count += 1
            continue

        try:
            transform_matrix = np.array(part_data['transform_matrix'])
            mesh = load_and_transform_mesh(stl_path, transform_matrix)
            color = get_part_color(part_name)

            meshes.append({
                'name': part_name,
                'mesh': mesh,
                'color': color
            })

            part_count += 1
            print(f"  [OK] {part_name}")

        except Exception as e:
            print(f"  [ERROR] {part_name}: {e}")
            missing_count += 1

    print(f"\nLoaded {part_count} parts, {missing_count} missing")

    if part_count == 0:
        print("ERROR: No parts loaded!")
        return 1

    # Calculate bounding box for all meshes
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
    print("\nGenerating multiview render...")
    fig = plt.figure(figsize=(24, 16))
    fig.suptitle('Robot Assembly - Multiview Render', fontsize=18, fontweight='bold')

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

        # Add all meshes to this view
        for mesh_data in meshes:
            poly = mesh_to_poly3d(mesh_data['mesh'], mesh_data['color'], alpha=0.8)
            ax.add_collection3d(poly)

        # Set view
        setup_view(ax, elev, azim, title)

        # Set limits (same for all views for consistency)
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
    print("\nKey subsystems:")
    print("  - Chassis (blue): Main plate, central hub, arm tower mount")
    print("  - Wheel assemblies (gray/black): 3x pivot bosses, wheel arms, DDSM115 motors")
    print("  - Gearboxes (orange): Axis 1/2/3 planetary gearboxes")
    print("  - SCARA arm (green): Links 1 & 2, extruder mount")
    print("  - Battery (purple): Battery cage")
    print("="*70)

    return 0

if __name__ == '__main__':
    sys.exit(main())
