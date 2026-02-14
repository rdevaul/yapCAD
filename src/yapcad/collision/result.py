"""Collision detection result data structures.

This module defines the CollisionResult dataclass and CollisionMethod enum
that represent the output of collision detection operations.

Copyright (c) 2026 yapCAD contributors
License: MIT
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum, auto
from typing import List, Tuple, Optional, Any, Dict


class CollisionMethod(Enum):
    """Method used for collision detection.

    Attributes:
        BREP: Exact BREP boolean intersection (pythonocc)
        MESH: Mesh-based sampling and containment (trimesh)
        AABB: Axis-aligned bounding box check only
        FCL: Flexible Collision Library (via trimesh)
        UNKNOWN: Method not specified or hybrid approach
    """
    BREP = auto()
    MESH = auto()
    AABB = auto()
    FCL = auto()
    UNKNOWN = auto()


@dataclass
class CollisionResult:
    """Result of collision detection between two parts.

    This dataclass captures all information about a collision check between
    two assembly parts, including whether they collide, the detection method
    used, collision metrics, and interface compatibility status.

    Attributes:
        part_a: Name/identifier of the first part
        part_b: Name/identifier of the second part
        collides: True if the parts geometrically intersect
        method: Detection method used (BREP, MESH, AABB, FCL)
        intersection_volume: Volume of intersection region in mm^3 (if computed)
        penetration_depth: Maximum penetration depth in mm (if computed)
        contact_points: List of contact/intersection points as (x, y, z) tuples
        compatible_interface: True if overlap is due to compatible interfaces
            (e.g., meshing gears with matching module/pressure angle)
        interface_names: Names of compatible interfaces that explain the overlap
        error_message: Error message if detection failed
        metadata: Additional method-specific data

    Example:
        >>> result = CollisionResult(
        ...     part_a="SUN_GEAR",
        ...     part_b="PLANET_GEAR_1",
        ...     collides=True,
        ...     method=CollisionMethod.BREP,
        ...     intersection_volume=15.3,
        ...     compatible_interface=True,
        ...     interface_names=["sun_teeth", "planet_1_teeth"]
        ... )
        >>> if result.collides and not result.compatible_interface:
        ...     print(f"ERROR: Unintended collision between {result.part_a} and {result.part_b}")
        ... elif result.collides:
        ...     print(f"OK: Expected overlap (gear mesh) between {result.part_a} and {result.part_b}")

    Notes:
        - `collides=True` with `compatible_interface=True` indicates expected overlap
          (e.g., gear teeth meshing, threads engaging)
        - `collides=True` with `compatible_interface=False` indicates an actual
          collision that needs to be resolved
        - `intersection_volume` and `penetration_depth` may be None/0 if not
          computed by the detection method used
    """
    part_a: str
    part_b: str
    collides: bool
    method: CollisionMethod = CollisionMethod.UNKNOWN
    intersection_volume: Optional[float] = None
    penetration_depth: float = 0.0
    contact_points: List[Tuple[float, float, float]] = field(default_factory=list)
    compatible_interface: bool = False
    interface_names: List[str] = field(default_factory=list)
    error_message: str = ""
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __str__(self) -> str:
        """Human-readable string representation."""
        if self.error_message:
            return f"[ERROR] {self.part_a} <-> {self.part_b}: {self.error_message}"

        if not self.collides:
            return f"[OK] {self.part_a} <-> {self.part_b}: No collision"

        if self.compatible_interface:
            ifaces = ", ".join(self.interface_names) if self.interface_names else "compatible"
            return f"[INTERFACE] {self.part_a} <-> {self.part_b}: Expected overlap ({ifaces})"

        details = []
        if self.intersection_volume is not None and self.intersection_volume > 0:
            details.append(f"volume={self.intersection_volume:.3f}mm^3")
        if self.penetration_depth > 0:
            details.append(f"depth={self.penetration_depth:.3f}mm")
        if self.contact_points:
            details.append(f"contacts={len(self.contact_points)}")

        detail_str = ", ".join(details) if details else "detected"
        return f"[COLLISION] {self.part_a} <-> {self.part_b}: {detail_str}"

    def __repr__(self) -> str:
        """Detailed repr for debugging."""
        return (
            f"CollisionResult(part_a={self.part_a!r}, part_b={self.part_b!r}, "
            f"collides={self.collides}, method={self.method.name}, "
            f"intersection_volume={self.intersection_volume}, "
            f"penetration_depth={self.penetration_depth}, "
            f"compatible_interface={self.compatible_interface})"
        )

    @property
    def is_error(self) -> bool:
        """True if collision represents an actual error (not expected interface)."""
        return self.collides and not self.compatible_interface

    @property
    def is_interface_overlap(self) -> bool:
        """True if collision is due to compatible interface (expected overlap)."""
        return self.collides and self.compatible_interface

    @property
    def pair_key(self) -> Tuple[str, str]:
        """Canonical pair key (sorted alphabetically) for deduplication."""
        return tuple(sorted([self.part_a, self.part_b]))

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization.

        Returns:
            Dictionary representation suitable for JSON export
        """
        return {
            "part_a": self.part_a,
            "part_b": self.part_b,
            "collides": self.collides,
            "method": self.method.name,
            "intersection_volume": self.intersection_volume,
            "penetration_depth": self.penetration_depth,
            "contact_points": self.contact_points,
            "compatible_interface": self.compatible_interface,
            "interface_names": self.interface_names,
            "error_message": self.error_message,
            "metadata": self.metadata,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'CollisionResult':
        """Create CollisionResult from dictionary.

        Args:
            data: Dictionary with collision result data

        Returns:
            CollisionResult instance
        """
        method_name = data.get("method", "UNKNOWN")
        try:
            method = CollisionMethod[method_name]
        except KeyError:
            method = CollisionMethod.UNKNOWN

        return cls(
            part_a=data["part_a"],
            part_b=data["part_b"],
            collides=data["collides"],
            method=method,
            intersection_volume=data.get("intersection_volume"),
            penetration_depth=data.get("penetration_depth", 0.0),
            contact_points=[tuple(p) for p in data.get("contact_points", [])],
            compatible_interface=data.get("compatible_interface", False),
            interface_names=data.get("interface_names", []),
            error_message=data.get("error_message", ""),
            metadata=data.get("metadata", {}),
        )

    @classmethod
    def no_collision(cls, part_a: str, part_b: str,
                     method: CollisionMethod = CollisionMethod.UNKNOWN) -> 'CollisionResult':
        """Factory for creating a no-collision result.

        Args:
            part_a: First part name
            part_b: Second part name
            method: Detection method used

        Returns:
            CollisionResult with collides=False
        """
        return cls(
            part_a=part_a,
            part_b=part_b,
            collides=False,
            method=method,
        )

    @classmethod
    def error(cls, part_a: str, part_b: str, message: str) -> 'CollisionResult':
        """Factory for creating an error result.

        Args:
            part_a: First part name
            part_b: Second part name
            message: Error message describing what went wrong

        Returns:
            CollisionResult with error_message set
        """
        return cls(
            part_a=part_a,
            part_b=part_b,
            collides=False,
            error_message=message,
        )
