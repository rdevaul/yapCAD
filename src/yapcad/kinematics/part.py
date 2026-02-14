"""Kinematic part (tree node) for kinematic chains.

Copyright (c) 2026 yapCAD contributors
License: MIT
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional, Tuple, Any, List

from .transform import Transform
from .joint import Joint, JointType
from .frame import CoordinateFrame


@dataclass
class KinematicPart:
    """Node in a kinematic tree.

    A KinematicPart represents a rigid body in an articulated assembly.
    It connects to a parent part via a joint and can have multiple
    named coordinate frames for attaching children or defining features.

    Attributes:
        name: Unique identifier for this part
        parent: Name of parent part (None for root)
        parent_frame: Frame on parent to attach to (default "ORIGIN")
        joint: Joint connecting to parent
        frames: Dict of named coordinate frames on this part
        stl_path: Optional path to STL geometry file
        is_printable: Whether this is a 3D printable part
        material: Material specification (e.g., "PETG", "aluminum")
        color: RGB color tuple for visualization
        description: Human-readable description

    Example::

        # Create a part with frames
        link = KinematicPart(
            name="LINK1",
            parent="BASE",
            parent_frame="MOUNT",
            joint=Joint("joint1", JointType.REVOLUTE),
        )

        # Add frames for child attachments
        link.add_frame(
            "END_EFFECTOR",
            Transform.from_translation(100, 0, 0),
            "Tool attachment point"
        )
    """

    name: str
    parent: Optional[str] = None
    parent_frame: str = "ORIGIN"
    joint: Joint = field(default_factory=lambda: Joint("default", JointType.FIXED))
    frames: Dict[str, CoordinateFrame] = field(default_factory=dict)

    # Geometry and visualization
    stl_path: Optional[str] = None
    is_printable: bool = True
    material: str = "PETG"
    color: Tuple[float, float, float] = (0.5, 0.5, 0.5)
    description: str = ""

    # Internal caching
    _world_transform: Optional[Transform] = field(
        default=None, repr=False, compare=False
    )
    _transform_dirty: bool = field(default=True, repr=False, compare=False)

    def __post_init__(self):
        """Ensure ORIGIN frame exists."""
        if "ORIGIN" not in self.frames:
            self.frames["ORIGIN"] = CoordinateFrame(
                name="ORIGIN",
                transform=Transform.identity(),
                description="Part origin"
            )

    def add_frame(
        self,
        name: str,
        transform: Transform,
        description: str = "",
    ) -> CoordinateFrame:
        """Add a coordinate frame to this part.

        :param name: Unique frame name
        :param transform: Transform from part origin to frame
        :param description: Human-readable description
        :returns: The created frame
        """
        frame = CoordinateFrame(
            name=name,
            transform=transform,
            description=description,
        )
        self.frames[name] = frame
        return frame

    def get_frame(self, name: str) -> Optional[CoordinateFrame]:
        """Get a coordinate frame by name.

        :param name: Frame name
        :returns: CoordinateFrame or None if not found
        """
        return self.frames.get(name)

    def invalidate_cache(self) -> None:
        """Mark world transform cache as dirty."""
        object.__setattr__(self, '_transform_dirty', True)
        object.__setattr__(self, '_world_transform', None)

    def set_cached_world_transform(self, transform: Transform) -> None:
        """Set cached world transform (called by KinematicChain)."""
        object.__setattr__(self, '_world_transform', transform)
        object.__setattr__(self, '_transform_dirty', False)

    def get_cached_world_transform(self) -> Optional[Transform]:
        """Get cached world transform if valid."""
        if self._transform_dirty:
            return None
        return self._world_transform

    def is_cache_valid(self) -> bool:
        """Check if world transform cache is valid."""
        return not self._transform_dirty

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "name": self.name,
            "parent": self.parent,
            "parent_frame": self.parent_frame,
            "joint": self.joint.to_dict(),
            "frames": {
                name: frame.to_dict()
                for name, frame in self.frames.items()
            },
            "stl_path": self.stl_path,
            "is_printable": self.is_printable,
            "material": self.material,
            "color": list(self.color),
            "description": self.description,
        }

    @classmethod
    def from_dict(cls, data: dict) -> KinematicPart:
        """Create part from dictionary."""
        frames = {}
        for name, frame_data in data.get("frames", {}).items():
            frames[name] = CoordinateFrame.from_dict(frame_data)

        return cls(
            name=data["name"],
            parent=data.get("parent"),
            parent_frame=data.get("parent_frame", "ORIGIN"),
            joint=Joint.from_dict(data["joint"]) if "joint" in data else Joint("default", JointType.FIXED),
            frames=frames,
            stl_path=data.get("stl_path"),
            is_printable=data.get("is_printable", True),
            material=data.get("material", "PETG"),
            color=tuple(data.get("color", [0.5, 0.5, 0.5])),
            description=data.get("description", ""),
        )


__all__ = ["KinematicPart"]
