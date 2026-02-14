#!/usr/bin/env python3
"""
Simple Assembly Example - Demonstrating the yapCAD Assembly System

This example shows how to use yapCAD's datum-driven assembly system to
position parts using geometric constraints rather than hardcoded transforms.

Key concepts demonstrated:
1. Define parts with datum features (geometric reference points/planes/axes)
2. Create mate constraints between datums
3. Use the MateConstraintSolver to compute transforms
4. Build a KinematicChain from the solved positions
5. Export positions to JSON

The example creates a simple assembly:
- A base plate with a mounting face
- A bracket that mounts to the base plate
- Datums on each part define the mounting interface
- A COINCIDENT mate constraint aligns the parts face-to-face

No external files needed - all geometry and datums are defined programmatically.

Copyright (c) 2026 yapCAD contributors
License: MIT
"""

import sys
import json
from pathlib import Path

# Add src to path so we can import yapcad
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import numpy as np
from yapcad.geom import point, vect
from yapcad.geom3d_util import prism
from yapcad.assembly.datum import Datum, DatumType, PartDefinition
from yapcad.assembly.datum_registry import DatumRegistry
from yapcad.assembly.mate import Mate, MateType
from yapcad.assembly.solver import MateConstraintSolver
from yapcad.kinematics import KinematicChain, KinematicPart, Joint, JointType, Transform


def create_base_plate():
    """Create a base plate part with mounting datums.

    The base plate has:
    - A top mounting face (where the bracket will attach)
    - A bolt hole circle pattern (4 holes at 45mm radius)
    - A central mounting axis

    Returns:
        PartDefinition with datums
    """
    print("\n[1] Creating BASE_PLATE part...")

    # Create the part definition
    base = PartDefinition(
        name="BASE_PLATE",
        geometry_source="prism(100, 100, 10)",
        is_printable=True,
        material="PETG",
        description="Base mounting plate with 4 bolt holes"
    )

    # Add top mounting face datum (Z=5mm, facing +Z)
    # This is where the bracket will mate to the base
    base.add_datum(Datum(
        name="top_face",
        datum_type=DatumType.PLANE,
        origin=point(0, 0, 5),      # Top surface of 10mm thick plate
        normal=vect(0, 0, 1, 0),     # Normal pointing up (+Z)
        description="Top mounting surface"
    ))

    # Add bolt hole circle datum (4 x M3 holes at 45mm radius)
    # This ensures bolt holes align with the bracket
    base.add_datum(Datum(
        name="bolt_circle",
        datum_type=DatumType.CIRCLE,
        origin=point(0, 0, 5),       # Center at top face
        normal=vect(0, 0, 1, 0),     # Circle in XY plane
        radius=45.0,                 # 45mm bolt circle diameter
        description="4 x M3 mounting holes at 45mm radius"
    ))

    # Add central axis (for alignment reference)
    base.add_datum(Datum(
        name="central_axis",
        datum_type=DatumType.AXIS,
        origin=point(0, 0, 0),
        direction=vect(0, 0, 1, 0),  # Vertical axis
        description="Central vertical axis"
    ))

    print(f"  - Added {len(base.datums)} datums: {list(base.datums.keys())}")
    return base


def create_mounting_bracket():
    """Create a mounting bracket part with interface datums.

    The bracket has:
    - A bottom mounting face (mates to base plate top face)
    - A bolt hole circle (must align with base plate holes)
    - A vertical mounting post

    Returns:
        PartDefinition with datums
    """
    print("\n[2] Creating MOUNTING_BRACKET part...")

    # Create the part definition
    bracket = PartDefinition(
        name="MOUNTING_BRACKET",
        geometry_source="prism(80, 80, 15)",
        is_printable=True,
        material="PETG",
        description="L-shaped mounting bracket"
    )

    # Add bottom mounting face datum (Z=0, facing -Z)
    # This mates with the base plate's top face
    # The normals will point in opposite directions (face-to-face contact)
    bracket.add_datum(Datum(
        name="mount_face",
        datum_type=DatumType.PLANE,
        origin=point(0, 0, 0),       # Bottom of bracket
        normal=vect(0, 0, -1, 0),    # Normal pointing down (-Z)
        description="Bottom mounting surface (mates to base)"
    ))

    # Add bolt hole circle (must match base plate pattern)
    bracket.add_datum(Datum(
        name="bolt_circle",
        datum_type=DatumType.CIRCLE,
        origin=point(0, 0, 0),       # Center at bottom face
        normal=vect(0, 0, -1, 0),    # Circle in XY plane
        radius=45.0,                 # Same 45mm radius as base
        description="4 x M3 mounting holes matching base pattern"
    ))

    # Add tool mounting point at top of bracket
    bracket.add_datum(Datum(
        name="tool_mount",
        datum_type=DatumType.POINT,
        origin=point(0, 0, 15),
        description="Tool attachment point"
    ))

    print(f"  - Added {len(bracket.datums)} datums: {list(bracket.datums.keys())}")
    return bracket


def register_parts_in_datum_registry(base, bracket):
    """Register parts in the DatumRegistry for solver lookup.

    The MateConstraintSolver uses the DatumRegistry to find datum features
    when solving constraints.

    Args:
        base: BASE_PLATE PartDefinition
        bracket: MOUNTING_BRACKET PartDefinition
    """
    print("\n[3] Registering parts in DatumRegistry...")

    # Register base plate datums
    DatumRegistry.register_source(
        source_id="BASE_PLATE",
        datums=base.datums,
        source_type="programmatic",
        metadata={"description": base.description}
    )

    # Register bracket datums
    DatumRegistry.register_source(
        source_id="MOUNTING_BRACKET",
        datums=bracket.datums,
        source_type="programmatic",
        metadata={"description": bracket.description}
    )

    print(f"  - Registered sources: {DatumRegistry.list_sources()}")


def create_mate_constraint():
    """Create a COINCIDENT mate between base and bracket.

    This constraint specifies that:
    - The bracket's mount_face should coincide with the base's top_face
    - Coincident means: faces touch, normals point in opposite directions
    - This is a "face-to-face" constraint, common in mechanical assemblies

    Returns:
        Mate object defining the constraint
    """
    print("\n[4] Creating COINCIDENT mate constraint...")

    mate = Mate(
        name="base_to_bracket",
        mate_type=MateType.COINCIDENT,
        part_a="BASE_PLATE",         # Parent part (fixed)
        datum_a="top_face",           # Base's top mounting face
        part_b="MOUNTING_BRACKET",    # Child part (will be positioned)
        datum_b="mount_face",         # Bracket's bottom mounting face
        offset=0.0                    # No gap between faces
    )

    print(f"  - Mate: {mate.part_a}.{mate.datum_a} <-> {mate.part_b}.{mate.datum_b}")
    print(f"  - Type: {mate.mate_type.value}")
    return mate


def solve_mate(mate):
    """Solve the mate constraint to compute the bracket's transform.

    The MateConstraintSolver uses the datum features to compute the 4x4
    transformation matrix that positions the child part (bracket) relative
    to the parent part (base plate) such that the constraint is satisfied.

    Args:
        mate: Mate constraint to solve

    Returns:
        4x4 numpy array transform matrix
    """
    print("\n[5] Solving mate constraint...")

    # Create solver instance
    solver = MateConstraintSolver(
        tolerance=0.001,        # Position tolerance (mm)
        angle_tolerance=0.1     # Angular tolerance (degrees)
    )

    # Solve the constraint
    result = solver.solve_mate(mate)

    if not result.success:
        print(f"  ERROR: Mate solve failed: {result.error_message}")
        sys.exit(1)

    print(f"  - Solve succeeded!")
    print(f"  - Residual error: {result.residual:.6f}")

    # Extract the transform
    transform = result.transform
    print(f"  - Computed transform shape: {transform.shape}")

    # Display the transform
    print("\n  Transform matrix:")
    for row in transform:
        print(f"    [{row[0]:7.3f} {row[1]:7.3f} {row[2]:7.3f} {row[3]:7.3f}]")

    return transform


def build_kinematic_chain(bracket_transform):
    """Build a kinematic chain from the solved assembly.

    The KinematicChain represents the assembly as a tree of parts with
    parent-child relationships and joints. This enables:
    - Computing world transforms for all parts
    - Animating assemblies
    - Exporting to URDF, SDF, or other formats

    Args:
        bracket_transform: 4x4 transform matrix for bracket (from solver)

    Returns:
        KinematicChain with base and bracket parts
    """
    print("\n[6] Building kinematic chain...")

    # Create the chain
    chain = KinematicChain("simple_assembly")

    # Add base plate as root (attached to world origin)
    # The base plate has a fixed joint and stays at the world origin
    base_part = KinematicPart(
        name="BASE_PLATE",
        parent=None,                          # Root part (no parent)
        joint=Joint("base_fixed", JointType.FIXED)
    )
    chain.add_part(base_part)
    print(f"  - Added BASE_PLATE (root)")

    # Add bracket as child of base
    # Create a frame on the base plate for the bracket attachment point
    # This frame is positioned by the mate solver result
    base_part.add_frame(
        "BRACKET_MOUNT",
        Transform(bracket_transform),
        "Bracket mounting location from mate solver"
    )

    # Add the bracket, attached to the BRACKET_MOUNT frame
    bracket_part = KinematicPart(
        name="MOUNTING_BRACKET",
        parent="BASE_PLATE",
        parent_frame="BRACKET_MOUNT",           # Attach to the mount frame
        joint=Joint("bracket_fixed", JointType.FIXED)
    )
    chain.add_part(bracket_part)
    print(f"  - Added MOUNTING_BRACKET (child of BASE_PLATE at BRACKET_MOUNT frame)")

    return chain


def print_world_transforms(chain):
    """Compute and print world transforms for all parts.

    World transforms show the absolute position and orientation of each
    part in the global coordinate system.

    Args:
        chain: KinematicChain to query
    """
    print("\n[7] Computing world transforms...")

    for part_name in ["BASE_PLATE", "MOUNTING_BRACKET"]:
        world_tf = chain.get_world_transform(part_name)
        print(f"\n  {part_name} world transform:")
        for row in world_tf.matrix:
            print(f"    [{row[0]:7.3f} {row[1]:7.3f} {row[2]:7.3f} {row[3]:7.3f}]")

        # Extract position
        position = world_tf.matrix[:3, 3]
        print(f"    Position: [{position[0]:.3f}, {position[1]:.3f}, {position[2]:.3f}]")


def export_to_json(chain, output_path):
    """Export the assembly positions to JSON.

    This creates a JSON file with all part transforms, which can be used
    for rendering, visualization, or import into other tools.

    Args:
        chain: KinematicChain to export
        output_path: Path to output JSON file
    """
    print(f"\n[8] Exporting positions to JSON...")

    chain.export_json(output_path)
    print(f"  - Exported to: {output_path}")

    # Also print the contents
    with open(output_path, 'r') as f:
        data = json.load(f)

    print("\n  JSON contents:")
    print(json.dumps(data, indent=2))


def main():
    """Main assembly demonstration."""
    print("=" * 70)
    print("yapCAD Simple Assembly Example")
    print("=" * 70)
    print("\nDemonstrating datum-driven assembly with mate constraints")

    # Step 1: Create parts with datums
    base = create_base_plate()
    bracket = create_mounting_bracket()

    # Validate datums
    base_issues = base.validate_datums()
    bracket_issues = bracket.validate_datums()
    if base_issues or bracket_issues:
        print("\nWARNING: Datum validation issues:")
        for issue in base_issues + bracket_issues:
            print(f"  - {issue}")

    # Step 2: Register parts in datum registry
    register_parts_in_datum_registry(base, bracket)

    # Step 3: Create mate constraint
    mate = create_mate_constraint()

    # Step 4: Solve mate to get bracket transform
    bracket_transform = solve_mate(mate)

    # Step 5: Build kinematic chain
    chain = build_kinematic_chain(bracket_transform)

    # Step 6: Print world transforms
    print_world_transforms(chain)

    # Step 7: Export to JSON
    output_dir = Path(__file__).parent / "output"
    output_dir.mkdir(exist_ok=True)
    output_path = output_dir / "simple_assembly_positions.json"
    export_to_json(chain, str(output_path))

    print("\n" + "=" * 70)
    print("Assembly complete!")
    print("=" * 70)
    print("\nKey takeaways:")
    print("  1. Datums define explicit geometric references on parts")
    print("  2. Mates specify constraints between datums (not hardcoded offsets)")
    print("  3. The solver computes transforms from constraints")
    print("  4. KinematicChain manages the assembly tree")
    print("  5. Positions can be exported to JSON for rendering/visualization")
    print("\nThis datum-driven approach eliminates hardcoded transforms and")
    print("makes assemblies easier to maintain and modify.")
    print("=" * 70)


if __name__ == "__main__":
    main()
