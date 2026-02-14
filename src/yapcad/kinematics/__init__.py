"""Kinematic Chain System for Multi-Body Assemblies.

====================
OVERVIEW
====================

This package provides a complete kinematic tree system for representing
articulated assemblies with joints, frames, and transforms. It forms the
foundation for computing world positions of parts in multi-body mechatronic
assemblies.

Key Concepts:
    - **Transform**: 4x4 homogeneous matrix for 6DOF positioning
    - **Joint**: Connection between parts (fixed, revolute, prismatic, etc.)
    - **CoordinateFrame**: Named reference frame on a part
    - **KinematicPart**: Tree node with parent, frames, and joint
    - **KinematicChain**: Tree container with world transform computation

Quick Start::

    from yapcad.kinematics import (
        KinematicChain, KinematicPart, Joint, JointType, Transform
    )

    # Create a kinematic chain
    chain = KinematicChain("my_robot")

    # Add base part (attached to world origin)
    base = KinematicPart(
        name="BASE",
        parent=None,
        joint=Joint("base_joint", JointType.FIXED),
    )
    chain.add_part(base)

    # Add link with revolute joint
    link = KinematicPart(
        name="LINK1",
        parent="BASE",
        parent_frame="MOUNT",
        joint=Joint("joint1", JointType.REVOLUTE, axis=(0, 0, 1)),
    )
    chain.add_part(link)

    # Set joint angle and get world transform
    chain.set_joint_value("LINK1", 45.0)  # degrees
    world_tf = chain.get_world_transform("LINK1")

    # Export to JSON
    chain.export_json("positions.json")

Kinematic chain system contributed by Jeremy Mika.

Copyright (c) 2026 yapCAD contributors
License: MIT
"""

from __future__ import annotations

from .transform import Transform
from .joint import Joint, JointType
from .frame import CoordinateFrame
from .part import KinematicPart
from .chain import KinematicChain

__all__ = [
    # Core transform
    "Transform",
    # Joint system
    "Joint",
    "JointType",
    # Frame system
    "CoordinateFrame",
    # Part and chain
    "KinematicPart",
    "KinematicChain",
]

__version__ = "0.1.0"
