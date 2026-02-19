#!/usr/bin/env python3
"""
6-View Matplotlib X-Ray Visualization of Robot Assembly

Generates a 2x3 grid showing all 6 orthographic views:
- Front (looking at +Y face, viewing along -Y axis)
- Back (looking at -Y face, viewing along +Y axis)
- Left (looking at +X face, viewing along -X axis)
- Right (looking at -X face, viewing along +X axis)
- Top (looking at +Z face, viewing along -Z axis)
- Bottom (looking at -Z face, viewing along +Z axis)

Uses transparent/wireframe rendering for x-ray effect with part labels.
Loads all STL files from kinematic chain and applies world transforms.
"""

import json
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import trimesh
from pathlib import Path
from typing import Dict, Tuple, Optional
import sys

# Import viz utilities
sys.path.insert(0, str(Path(__file__).parent))
from viz_utils import set_equal_3d_aspect, load_stl_safe


# X-ray color scheme - transparent with visible edges
SUBSYSTEM_COLORS = {
    'chassis': (0.5, 0.5, 0.5, 0.2),      # Gray, very transparent
    'wheel': (1.0, 0.6, 0.2, 0.25),       # Orange
    'scara_tower': (0.2, 0.4, 0.8, 0.3),  # Blue
    'scara_link': (0.3, 0.8, 0.3, 0.3),   # Green
    'gearbox_axis1': (0.9, 0.4, 0.2, 0.3), # Red-orange (Axis 1)
    'gearbox_axis2': (0.9, 0.5, 0.3, 0.3), # Orange (Axis 2)
    'gearbox_axis3': (0.95, 0.6, 0.4, 0.3), # Light orange (Axis 3)
    'gearbox_axis4': (1.0, 0.7, 0.5, 0.3), # Very light orange (Axis 4)
    'servo': (0.0, 0.7, 0.0, 0.4),        # Green - highlight servos
    'power': (0.6, 0.2, 0.8, 0.25),       # Purple
    'default': (0.6, 0.6, 0.6, 0.2),      # Medium gray
}


def get_part_color(part_name: str) -> Tuple[float, float, float, float]:
    """Determine color based on part name/subsystem."""
    part_lower = part_name.lower()

    if 'chassis' in part_lower or 'pivot' in part_lower or 'hub' in part_lower or 'tower_mount' in part_lower:
        return SUBSYSTEM_COLORS['chassis']
    elif 'wheel' in part_lower:
        return SUBSYSTEM_COLORS['wheel']
    elif 'arm_tower' in part_lower:
        return SUBSYSTEM_COLORS['scara_tower']
    elif 'link' in part_lower or 'extruder_mount' in part_lower:
        return SUBSYSTEM_COLORS['scara_link']
    elif 'axis1' in part_lower and 'gear' in part_lower or 'axis1' in part_lower and ('ring' in part_lower or 'carrier' in part_lower):
        return SUBSYSTEM_COLORS['gearbox_axis1']
    elif 'axis2' in part_lower and 'gear' in part_lower or 'axis2' in part_lower and ('ring' in part_lower or 'carrier' in part_lower or 'mount' in part_lower):
        return SUBSYSTEM_COLORS['gearbox_axis2']
    elif 'axis3' in part_lower and 'gear' in part_lower or 'axis3' in part_lower and ('ring' in part_lower or 'carrier' in part_lower):
        return SUBSYSTEM_COLORS['gearbox_axis3']
    elif 'axis4' in part_lower:
        return SUBSYSTEM_COLORS['gearbox_axis4']
    elif 'servo' in part_lower or 'motor' in part_lower:
        return SUBSYSTEM_COLORS['servo']
    elif 'battery' in part_lower or 'cage' in part_lower:
        return SUBSYSTEM_COLORS['power']
    else:
        return SUBSYSTEM_COLORS['default']


def find_stl_file(part_name: str, stl_base_dir: Path) -> Optional[Path]:
    """Find the STL file for a given part name."""
    # Map part names to STL files
    stl_mappings = {
        # Chassis (CENTRAL_HUB is unioned into chassis_with_battery_mounts.stl)
        'CHASSIS_PLATE': 'chassis/chassis_with_battery_mounts.stl',
        'PIVOT_BOSS_1': 'chassis/pivot_boss.stl',
        'PIVOT_BOSS_2': 'chassis/pivot_boss.stl',
        'PIVOT_BOSS_3': 'chassis/pivot_boss.stl',
        'ARM_TOWER_MOUNT': 'chassis/arm_tower_mount.stl',
        'ARM_TOWER': 'scara/arm_tower.stl',

        # Wheel assemblies
        'WHEEL_ARM_1': 'drive/wheel_arm.stl',
        'WHEEL_ARM_2': 'drive/wheel_arm.stl',
        'WHEEL_ARM_3': 'drive/wheel_arm.stl',

        # COTS motors
        'DDSM115_MOTOR_1': 'cots/ddsm115_motor_official.stl',
        'DDSM115_MOTOR_2': 'cots/ddsm115_motor_official.stl',
        'DDSM115_MOTOR_3': 'cots/ddsm115_motor_official.stl',
        'AXIS1_SERVO_XH540': 'cots/xh540_cots.stl',
        'AXIS2_SERVO_XH430': 'cots/xh430_cots.stl',
        'AXIS3_SERVO_XH430': 'cots/xh430_cots.stl',
        'AXIS4_SERVO_XL330': 'cots/xh430_cots.stl',  # Placeholder

        # Gearbox axis 1
        'AXIS1_RING_HOUSING': 'gearbox/axis1_ring_housing.stl',
        'AXIS1_SUN_GEAR': 'gearbox/axis1_sun_gear.stl',
        'AXIS1_PLANET_GEAR_1': 'gearbox/axis1_planet_gear.stl',
        'AXIS1_PLANET_GEAR_2': 'gearbox/axis1_planet_gear.stl',
        'AXIS1_PLANET_GEAR_3': 'gearbox/axis1_planet_gear.stl',
        'AXIS1_CARRIER_TOP': 'gearbox/axis1_carrier_top.stl',
        'AXIS1_CARRIER_BOTTOM': 'gearbox/axis1_carrier_bottom.stl',

        # Gearbox axis 2
        'AXIS2_MOUNT': 'scara/axis2_mount.stl',
        'AXIS2_RING_HOUSING': 'gearbox/axis2_ring_housing.stl',
        'AXIS2_SUN_GEAR': 'gearbox/axis2_sun_gear.stl',
        'AXIS2_PLANET_GEAR_1': 'gearbox/axis2_planet_gear.stl',
        'AXIS2_PLANET_GEAR_2': 'gearbox/axis2_planet_gear.stl',
        'AXIS2_PLANET_GEAR_3': 'gearbox/axis2_planet_gear.stl',
        'AXIS2_CARRIER_TOP': 'gearbox/axis2_carrier_top.stl',
        'AXIS2_CARRIER_BOTTOM': 'gearbox/axis2_carrier_bottom.stl',

        # Gearbox axis 3
        'AXIS3_RING_HOUSING': 'gearbox/axis3_ring_housing.stl',
        'AXIS3_SUN_GEAR': 'gearbox/axis3_sun_gear.stl',
        'AXIS3_PLANET_GEAR_1': 'gearbox/axis3_planet_gear.stl',
        'AXIS3_PLANET_GEAR_2': 'gearbox/axis3_planet_gear.stl',
        'AXIS3_PLANET_GEAR_3': 'gearbox/axis3_planet_gear.stl',
        'AXIS3_CARRIER_TOP': 'gearbox/axis3_carrier_top.stl',
        'AXIS3_CARRIER_BOTTOM': 'gearbox/axis3_carrier_bottom.stl',

        # SCARA links
        'LINK_2_3': 'scara/link1_attractive.stl',
        'LINK_3_4': 'scara/link2_attractive.stl',

        # Axis 4 (wrist)
        'AXIS4_GEARBOX': 'scara/axis4_housing.stl',
        'EXTRUDER_MOUNT': 'scara/axis4_output.stl',
        'BIQU_H2O_EXTRUDER': 'scara/h2o_extruder.stl',

        # Power
        'BATTERY_CAGE': 'power/battery_cage.stl',
    }

    if part_name in stl_mappings:
        stl_path = stl_base_dir / stl_mappings[part_name]
        if stl_path.exists():
            return stl_path

    return None


def render_mesh_xray(ax, mesh, transform_matrix, color, edgecolor='black', linewidth=0.2):
    """Render mesh with x-ray style (transparent with visible edges)."""
    # Apply transform
    vertices = mesh.vertices.copy()
    ones = np.ones((vertices.shape[0], 1))
    verts_h = np.hstack([vertices, ones])
    verts_transformed = (np.array(transform_matrix) @ verts_h.T).T
    vertices = verts_transformed[:, :3]

    faces = mesh.faces

    # Create transparent polygon collection with visible edges
    poly = Poly3DCollection(
        vertices[faces],
        alpha=color[3],
        edgecolor=edgecolor,
        linewidths=linewidth
    )
    poly.set_facecolor(color)
    ax.add_collection3d(poly)

    return vertices


def add_coordinate_axes(ax, origin, scale=50):
    """Add RGB XYZ coordinate axes at origin."""
    # X axis - Red
    ax.plot([origin[0], origin[0] + scale], [origin[1], origin[1]], [origin[2], origin[2]],
            'r-', linewidth=2, label='X')
    # Y axis - Green
    ax.plot([origin[0], origin[0]], [origin[1], origin[1] + scale], [origin[2], origin[2]],
            'g-', linewidth=2, label='Y')
    # Z axis - Blue
    ax.plot([origin[0], origin[0]], [origin[1], origin[1]], [origin[2], origin[2] + scale],
            'b-', linewidth=2, label='Z')


def setup_view(ax, elev, azim, title):
    """Configure a view with consistent styling."""
    ax.view_init(elev=elev, azim=azim)
    ax.set_title(title, fontsize=11, fontweight='bold', pad=10)
    ax.set_xlabel('X (mm)', fontsize=8)
    ax.set_ylabel('Y (mm)', fontsize=8)
    ax.set_zlabel('Z (mm)', fontsize=8)
    ax.grid(True, alpha=0.2, linestyle='--')
    ax.tick_params(labelsize=7)

    # Light background for x-ray effect
    ax.set_facecolor('#f8f8f8')
    ax.xaxis.pane.fill = True
    ax.yaxis.pane.fill = True
    ax.zaxis.pane.fill = True
    ax.xaxis.pane.set_alpha(0.1)
    ax.yaxis.pane.set_alpha(0.1)
    ax.zaxis.pane.set_alpha(0.1)


def annotate_key_parts(ax, positions, elev, azim):
    """Add labels for key components (view-dependent positioning)."""
    # Define key parts to label
    key_parts = {
        'ARM_TOWER': 'ARM_TOWER',
        'LINK_2_3': 'LINK_2_3',
        'LINK_3_4': 'LINK_3_4',
        'AXIS1_SERVO_XH540': 'SERVO_1',
        'AXIS2_SERVO_XH430': 'SERVO_2',
        'AXIS3_SERVO_XH430': 'SERVO_3',
        'AXIS1_RING_HOUSING': 'GEARBOX_1',
        'AXIS2_RING_HOUSING': 'GEARBOX_2',
        'AXIS3_RING_HOUSING': 'GEARBOX_3',
    }

    for part_name, label in key_parts.items():
        if part_name in positions:
            # Get world position from transform matrix
            matrix = np.array(positions[part_name]['transform_matrix'])
            pos = matrix[:3, 3]

            # Add small text annotation
            ax.text(pos[0], pos[1], pos[2], label,
                   fontsize=6, color='darkblue',
                   bbox=dict(boxstyle='round,pad=0.2', facecolor='yellow', alpha=0.7))


def main():
    """Main entry point."""
    # Paths
    base_dir = Path(__file__).parent
    positions_path = base_dir / 'output' / 'chain_assembly' / 'positions_home.json'
    stl_base_dir = base_dir / 'output' / 'stl'
    output_path = base_dir / 'output' / 'renders' / 'assembly_xray_6view.png'

    print("="*70)
    print("6-VIEW X-RAY VISUALIZATION - Robot Assembly")
    print("="*70)

    # Load positions
    print("\n[1/4] Loading assembly positions...")
    with open(positions_path) as f:
        positions = json.load(f)
    print(f"      Found {len(positions)} parts in kinematic chain")

    # Load all meshes
    print("\n[2/4] Loading STL files and applying transforms...")
    meshes = []
    all_vertices = []
    part_count = 0
    missing_count = 0

    for part_name, part_data in positions.items():
        stl_path = find_stl_file(part_name, stl_base_dir)

        if stl_path is None:
            print(f"      [SKIP] {part_name}: STL not found")
            missing_count += 1
            continue

        try:
            mesh = load_stl_safe(stl_path)
            if mesh is None:
                missing_count += 1
                continue

            transform_matrix = part_data['transform_matrix']
            color = get_part_color(part_name)

            meshes.append({
                'name': part_name,
                'mesh': mesh,
                'transform': transform_matrix,
                'color': color
            })

            # Track vertices for bounds
            vertices = mesh.vertices.copy()
            ones = np.ones((vertices.shape[0], 1))
            verts_h = np.hstack([vertices, ones])
            verts_transformed = (np.array(transform_matrix) @ verts_h.T).T
            all_vertices.extend(verts_transformed[:, :3])

            part_count += 1
            print(f"      [OK] {part_name}")

        except Exception as e:
            print(f"      [ERROR] {part_name}: {e}")
            missing_count += 1

    print(f"\n      Loaded: {part_count} parts")
    print(f"      Missing: {missing_count} parts")

    if part_count == 0:
        print("\nERROR: No parts loaded!")
        return 1

    # Calculate bounds
    print("\n[3/4] Calculating assembly bounds...")
    all_vertices = np.array(all_vertices)
    bounds_min = all_vertices.min(axis=0)
    bounds_max = all_vertices.max(axis=0)
    bounds_center = (bounds_min + bounds_max) / 2
    bounds_size = (bounds_max - bounds_min).max()

    print(f"      X: [{bounds_min[0]:.1f}, {bounds_max[0]:.1f}] mm")
    print(f"      Y: [{bounds_min[1]:.1f}, {bounds_max[1]:.1f}] mm")
    print(f"      Z: [{bounds_min[2]:.1f}, {bounds_max[2]:.1f}] mm")
    print(f"      Max dimension: {bounds_size:.1f} mm")

    # Create figure with 2x3 grid for 6 orthographic views
    print("\n[4/4] Rendering 6-view x-ray visualization...")
    fig = plt.figure(figsize=(24, 16), facecolor='white')
    fig.suptitle('Robot Assembly - 6-View X-Ray Visualization',
                 fontsize=18, fontweight='bold', y=0.98)

    # Define views: (elev, azim, title)
    # Note: matplotlib's view_init uses elev (elevation from XY plane) and azim (azimuth rotation)
    views = [
        (0, 0, 'Front View\n(looking along +X)'),
        (0, 180, 'Back View\n(looking along -X)'),
        (0, 90, 'Left View\n(looking along +Y)'),
        (0, -90, 'Right View\n(looking along -Y)'),
        (90, 0, 'Top View\n(looking along +Z)'),
        (-90, 0, 'Bottom View\n(looking along -Z)'),
    ]

    # Render each view
    for idx, (elev, azim, title) in enumerate(views, 1):
        print(f"      Rendering view {idx}/6: {title.split('(')[0].strip()}")
        ax = fig.add_subplot(2, 3, idx, projection='3d')

        # Render all meshes
        for mesh_data in meshes:
            render_mesh_xray(
                ax,
                mesh_data['mesh'],
                mesh_data['transform'],
                mesh_data['color'],
                edgecolor='black',
                linewidth=0.15
            )

        # Add coordinate axes at origin
        add_coordinate_axes(ax, origin=[0, 0, 0], scale=30)

        # Set up view
        setup_view(ax, elev, azim, title)

        # Set equal aspect ratio with padding
        margin = bounds_size * 0.15
        ax.set_xlim([bounds_center[0] - bounds_size/2 - margin,
                     bounds_center[0] + bounds_size/2 + margin])
        ax.set_ylim([bounds_center[1] - bounds_size/2 - margin,
                     bounds_center[1] + bounds_size/2 + margin])
        ax.set_zlim([bounds_center[2] - bounds_size/2 - margin,
                     bounds_center[2] + bounds_size/2 + margin])

        # Force equal box aspect
        ax.set_box_aspect([1, 1, 1])

        # Add labels for key parts (only on some views for clarity)
        if idx in [1, 5]:  # Front and top views
            annotate_key_parts(ax, positions, elev, azim)

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=SUBSYSTEM_COLORS['chassis'], label='Chassis'),
        Patch(facecolor=SUBSYSTEM_COLORS['scara_tower'], label='SCARA Tower'),
        Patch(facecolor=SUBSYSTEM_COLORS['scara_link'], label='SCARA Links'),
        Patch(facecolor=SUBSYSTEM_COLORS['gearbox_axis1'], label='Gearbox Axis 1'),
        Patch(facecolor=SUBSYSTEM_COLORS['gearbox_axis2'], label='Gearbox Axis 2'),
        Patch(facecolor=SUBSYSTEM_COLORS['gearbox_axis3'], label='Gearbox Axis 3'),
        Patch(facecolor=SUBSYSTEM_COLORS['servo'], label='Servos/Motors'),
        Patch(facecolor=SUBSYSTEM_COLORS['power'], label='Battery'),
    ]
    fig.legend(handles=legend_elements, loc='lower center',
               bbox_to_anchor=(0.5, -0.01), ncol=8, fontsize=10, frameon=True)

    # Adjust layout
    plt.tight_layout(rect=[0, 0.02, 1, 0.97])

    # Save
    print(f"\n[SAVE] Writing to {output_path}...")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')

    # Summary
    print("\n" + "="*70)
    print("6-VIEW X-RAY VISUALIZATION COMPLETE")
    print("="*70)
    print(f"Parts rendered: {part_count}")
    print(f"Parts missing: {missing_count}")
    print(f"Output: {output_path}")
    print("\nViews generated:")
    print("  1. Front (X+)  - Looking at YZ plane")
    print("  2. Back (X-)   - Opposite view")
    print("  3. Left (Y+)   - Looking at XZ plane")
    print("  4. Right (Y-)  - Opposite view")
    print("  5. Top (Z+)    - Looking at XY plane from above")
    print("  6. Bottom (Z-) - Looking at XY plane from below")
    print("\nKey features:")
    print("  - Transparent x-ray rendering with visible edges")
    print("  - Color-coded by subsystem")
    print("  - Coordinate axes shown (RGB = XYZ)")
    print("  - Key parts labeled on select views")
    print("="*70)

    return 0


if __name__ == '__main__':
    sys.exit(main())
