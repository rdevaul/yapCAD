"""Base interface volume classes for controlled overlap regions.

This module defines the base InterfaceVolume class and InterfaceType enum
for describing regions where parts are designed to overlap (gear teeth,
threads, press fits, etc.).

Copyright (c) 2026 yapCAD contributors
License: MIT
"""

from __future__ import annotations

import math
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Dict, List, Optional, Tuple, Any

# Try to import numpy for geometric calculations
try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None


class InterfaceType(Enum):
    """Types of interface volumes for controlled overlap.

    Each type represents a specific mechanical interface where
    geometric overlap is expected and designed.

    Attributes:
        GEAR_MESH: Involute gear teeth meshing
        THREAD: Screw thread engagement (bolt/nut, tap/screw)
        PRESS_FIT: Interference fit (shaft/bore)
        BEARING: Bearing race contact regions
        SPLINE: Spline coupling interface
        KEY: Key/keyway interface
        SNAP_FIT: Snap-fit detent mechanism
        CUSTOM: User-defined interface type
    """
    GEAR_MESH = auto()
    THREAD = auto()
    PRESS_FIT = auto()
    BEARING = auto()
    SPLINE = auto()
    KEY = auto()
    SNAP_FIT = auto()
    CUSTOM = auto()


@dataclass
class CompatibilityResult:
    """Result of checking interface compatibility.

    When two interface volumes overlap, this result indicates whether
    they are designed to mate properly (compatible) or represent an
    unintended collision (incompatible).

    Attributes:
        is_compatible: True if interfaces can mate correctly
        reason: Explanation of compatibility status
        warnings: Non-fatal compatibility concerns (e.g., fit tolerance warnings)
        overlap_volume: Estimated overlap volume in mm^3 (if compatible)
        required_clearance: Minimum clearance needed for assembly in mm

    Example:
        >>> result = sun_gear.check_compatibility(planet_gear)
        >>> if result.is_compatible:
        ...     print(f"Gears can mesh: {result.reason}")
        ...     if result.warnings:
        ...         print(f"Warnings: {result.warnings}")
        ... else:
        ...     print(f"Incompatible: {result.reason}")
    """
    is_compatible: bool
    reason: str
    warnings: List[str] = field(default_factory=list)
    overlap_volume: float = 0.0
    required_clearance: float = 0.0

    def __str__(self) -> str:
        status = "COMPATIBLE" if self.is_compatible else "INCOMPATIBLE"
        return f"{status}: {self.reason}"

    def __repr__(self) -> str:
        return (
            f"CompatibilityResult(is_compatible={self.is_compatible}, "
            f"reason={self.reason!r}, warnings={self.warnings}, "
            f"overlap_volume={self.overlap_volume:.3f})"
        )


@dataclass
class InterfaceVolume(ABC):
    """Base class for interface volumes that allow controlled overlap.

    Interface volumes define regions where parts are DESIGNED to overlap,
    such as gear teeth meshing, threads engaging, or press fits. The
    collision detector uses these to distinguish expected overlap from
    actual collisions.

    Subclasses must implement:
        - get_bounding_cylinder(): Return (radius, height) of enclosing cylinder
        - check_compatibility(): Check if another interface can mate properly
        - get_engagement_depth(): Return the engagement/overlap depth

    Attributes:
        name: Unique identifier for this interface
        part_name: Name of the part this interface belongs to
        interface_type: Type of interface (GEAR_MESH, THREAD, etc.)
        center: Center point of interface volume as (x, y, z) tuple
        axis: Primary axis direction as (x, y, z) tuple (rotation/engagement axis)
        description: Human-readable description
        metadata: Additional type-specific parameters

    Example:
        >>> class MyInterface(InterfaceVolume):
        ...     def get_bounding_cylinder(self):
        ...         return (10.0, 5.0)  # radius=10mm, height=5mm
        ...     def check_compatibility(self, other):
        ...         if isinstance(other, MyInterface):
        ...             return CompatibilityResult(True, "Compatible")
        ...         return CompatibilityResult(False, "Type mismatch")
        ...     def get_engagement_depth(self):
        ...         return 5.0
    """
    name: str
    part_name: str
    interface_type: InterfaceType
    center: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    axis: Tuple[float, float, float] = (0.0, 0.0, 1.0)
    description: str = ""
    metadata: Dict[str, Any] = field(default_factory=dict)

    @abstractmethod
    def get_bounding_cylinder(self) -> Tuple[float, float]:
        """Get cylindrical bounding volume for the interface.

        Returns:
            Tuple of (radius, height) in mm defining the minimum enclosing
            cylinder aligned with the interface axis.
        """
        pass

    @abstractmethod
    def check_compatibility(self, other: 'InterfaceVolume') -> CompatibilityResult:
        """Check if this interface is compatible with another.

        Two interfaces are compatible if they are designed to mate
        properly when overlapping (e.g., matching gear modules,
        matching thread pitches).

        Args:
            other: Another interface volume to check against

        Returns:
            CompatibilityResult indicating if interfaces can properly mate
        """
        pass

    @abstractmethod
    def get_engagement_depth(self) -> float:
        """Get the depth/length of interface engagement in mm.

        For different interface types:
            - Gears: face width
            - Threads: engagement length
            - Press fits: engagement length
            - Bearings: bearing width

        Returns:
            Engagement depth in mm
        """
        pass

    def overlaps_with(self, other: 'InterfaceVolume', tolerance: float = 0.1) -> bool:
        """Check if two interface volumes spatially overlap.

        This is a geometric overlap check, not a compatibility check.
        Two interfaces can overlap geometrically but be incompatible
        (e.g., mismatched gear modules).

        Args:
            other: Another interface volume to check
            tolerance: Distance tolerance for overlap detection in mm

        Returns:
            True if bounding volumes overlap within tolerance
        """
        if not HAS_NUMPY:
            # Fallback: simple center distance check
            dx = self.center[0] - other.center[0]
            dy = self.center[1] - other.center[1]
            dz = self.center[2] - other.center[2]
            dist = math.sqrt(dx*dx + dy*dy + dz*dz)

            r1, h1 = self.get_bounding_cylinder()
            r2, h2 = other.get_bounding_cylinder()

            # Simple check: combined radius plus tolerance
            return dist < (r1 + r2 + tolerance)

        # Get bounding cylinders
        r1, h1 = self.get_bounding_cylinder()
        r2, h2 = other.get_bounding_cylinder()

        # Convert to numpy arrays
        c1 = np.array(self.center)
        c2 = np.array(other.center)
        a1 = np.array(self.axis)
        a2 = np.array(other.axis)

        # Normalize axes
        a1_norm = np.linalg.norm(a1)
        a2_norm = np.linalg.norm(a2)
        if a1_norm > 1e-10:
            a1 = a1 / a1_norm
        if a2_norm > 1e-10:
            a2 = a2 / a2_norm

        # Check axis alignment
        axis_dot = abs(np.dot(a1, a2))
        if axis_dot < 0.99:
            # Axes not parallel - use center distance
            center_dist = np.linalg.norm(c2 - c1)
            return center_dist < (r1 + r2 + tolerance)

        # For parallel axes, check cylinder overlap
        center_vec = c2 - c1
        axial_dist = abs(np.dot(center_vec, a1))

        # Radial distance (perpendicular to axis)
        radial_vec = center_vec - np.dot(center_vec, a1) * a1
        radial_dist = np.linalg.norm(radial_vec)

        # Check axial overlap
        axial_overlap = axial_dist < (h1/2 + h2/2 + tolerance)

        # Check radial overlap
        radial_overlap = radial_dist < (r1 + r2 + tolerance)

        return axial_overlap and radial_overlap

    def transform_center(self, transform_matrix: Any) -> Tuple[float, float, float]:
        """Transform the interface center by a 4x4 transformation matrix.

        Args:
            transform_matrix: 4x4 transformation matrix (numpy array or list)

        Returns:
            Transformed center point as (x, y, z) tuple
        """
        if not HAS_NUMPY:
            # Can't transform without numpy
            return self.center

        tf = np.array(transform_matrix)
        if tf.shape != (4, 4):
            return self.center

        # Apply transform to center point
        pt = np.array([self.center[0], self.center[1], self.center[2], 1.0])
        pt_transformed = tf @ pt
        return (pt_transformed[0], pt_transformed[1], pt_transformed[2])

    def transform_axis(self, transform_matrix: Any) -> Tuple[float, float, float]:
        """Transform the interface axis by a 4x4 transformation matrix.

        Args:
            transform_matrix: 4x4 transformation matrix (numpy array or list)

        Returns:
            Transformed axis direction as (x, y, z) tuple (normalized)
        """
        if not HAS_NUMPY:
            return self.axis

        tf = np.array(transform_matrix)
        if tf.shape != (4, 4):
            return self.axis

        # Extract rotation part (upper-left 3x3)
        rot = tf[:3, :3]

        # Transform axis direction
        axis = np.array(self.axis)
        axis_transformed = rot @ axis

        # Normalize
        norm = np.linalg.norm(axis_transformed)
        if norm > 1e-10:
            axis_transformed = axis_transformed / norm

        return tuple(axis_transformed)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization.

        Returns:
            Dictionary representation suitable for JSON export
        """
        return {
            "name": self.name,
            "part_name": self.part_name,
            "interface_type": self.interface_type.name,
            "center": list(self.center),
            "axis": list(self.axis),
            "description": self.description,
            "metadata": self.metadata,
        }
