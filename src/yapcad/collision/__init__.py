"""Collision detection system for yapCAD assemblies.

This package provides collision detection capabilities for validating assembly
configurations. It supports multiple detection methods with automatic fallback:

1. **BREP-based detection** (requires pythonocc-core):
   Uses exact boolean intersection via BRepAlgoAPI_Common for precise
   collision volumes. This is the most accurate method.

2. **Mesh-based detection** (requires trimesh):
   Uses point sampling and containment checks for collision detection.
   Works with STL/mesh geometry when BREP is unavailable.

3. **AABB-based detection**:
   Fast bounding box checks for quick rejection of non-colliding pairs.
   Always available as a first-pass filter.

Key Features:
    - Pluggable geometry source via GeometryProvider protocol
    - Interface volume system for allowed overlaps (gear mesh, threads, etc.)
    - Detailed collision results with volumes, penetration depths, contact points
    - Integration with assembly constraint system

Design Principles:
    - No hardcoded part-to-file mappings (use GeometryProvider)
    - No project-specific exclusion lists (use InterfaceRegistry)
    - Works with any assembly via pluggable interfaces
    - Graceful degradation when optional dependencies unavailable

Quick Start:
    >>> from yapcad.collision import (
    ...     CollisionDetector,
    ...     CollisionResult,
    ...     InterfaceRegistry,
    ...     GearMeshInterface,
    ... )
    >>>
    >>> # Create detector with geometry provider
    >>> class MyGeometryProvider:
    ...     def get_geometry(self, part_name: str) -> Optional[Any]:
    ...         # Return STEP/STL path or mesh object
    ...         return f"parts/{part_name}.step"
    ...
    >>> detector = CollisionDetector(MyGeometryProvider())
    >>>
    >>> # Register allowed overlaps for gear meshing
    >>> registry = InterfaceRegistry()
    >>> registry.register(GearMeshInterface(
    ...     name="sun_teeth", part_name="SUN_GEAR",
    ...     module=0.75, teeth=18, pressure_angle=20.0, face_width=10.0
    ... ))
    >>> detector.set_interface_registry(registry)
    >>>
    >>> # Check assembly for collisions
    >>> world_transforms = {"SUN_GEAR": tf1, "PLANET_GEAR": tf2, ...}
    >>> results = detector.check_assembly(world_transforms)
    >>> for result in results:
    ...     if result.collides and not result.compatible_interface:
    ...         print(f"COLLISION: {result.part_a} <-> {result.part_b}")

Available Classes:
    Detection:
        CollisionDetector - Main collision detection class
        CollisionResult - Detailed collision result data
        GeometryProvider - Protocol for geometry source

    Interface Volumes (allowed overlaps):
        InterfaceVolume - Base class for interface volumes
        InterfaceType - Enum of interface types
        GearMeshInterface - Gear teeth meshing interface
        ThreadInterface - Screw thread engagement interface
        InterfaceRegistry - Registry for managing interfaces

See Also:
    - yapcad.assembly: Constraint-based assembly system
    - yapcad.octtree: Spatial indexing for acceleration

Collision detection system contributed by Jeremy Mika.

Copyright (c) 2026 yapCAD contributors
License: MIT
"""

from __future__ import annotations

# Result classes (no dependencies)
from .result import (
    CollisionResult,
    CollisionMethod,
)

# Interface volume system
from .interface import (
    InterfaceVolume,
    InterfaceType,
    CompatibilityResult,
)

from .gear_interface import (
    GearMeshInterface,
    check_gear_mesh_collision,
    create_planetary_gearbox_interfaces,
)

from .thread_interface import (
    ThreadInterface,
    ThreadType,
    ThreadClass,
)

from .registry import (
    InterfaceRegistry,
)

# Main detector
from .detector import (
    CollisionDetector,
    GeometryProvider,
)

__all__ = [
    # Detection
    "CollisionDetector",
    "CollisionResult",
    "CollisionMethod",
    "GeometryProvider",
    # Interface volumes
    "InterfaceVolume",
    "InterfaceType",
    "CompatibilityResult",
    "GearMeshInterface",
    "ThreadInterface",
    "ThreadType",
    "ThreadClass",
    "InterfaceRegistry",
    # Utility functions
    "check_gear_mesh_collision",
    "create_planetary_gearbox_interfaces",
]

__version__ = "0.1.0"
