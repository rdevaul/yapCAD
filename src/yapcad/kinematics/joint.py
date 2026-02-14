"""Joint types and Joint class for kinematic connections.

Copyright (c) 2026 yapCAD contributors
License: MIT
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum
from typing import Tuple, Optional

from .transform import Transform


class JointType(Enum):
    """Types of kinematic joints.

    Attributes:
        FIXED: No relative motion (rigid attachment)
        REVOLUTE: Rotation about single axis
        PRISMATIC: Translation along single axis
        CYLINDRICAL: Rotation + translation about same axis
        SPHERICAL: Ball joint (3 rotational DOF)
    """

    FIXED = "fixed"
    REVOLUTE = "revolute"
    PRISMATIC = "prismatic"
    CYLINDRICAL = "cylindrical"
    SPHERICAL = "spherical"


@dataclass
class Joint:
    """Connection between kinematic parts.

    A joint defines how a child part connects to its parent part,
    including the base transform (fixed offset) and any articulated
    motion (rotation or translation along an axis).

    Attributes:
        name: Unique identifier for this joint
        joint_type: Type of joint (FIXED, REVOLUTE, etc.)
        base_transform: Fixed transform from parent frame to joint origin
        axis: Motion axis for revolute/prismatic joints
        value: Current joint position (angle in deg or distance in mm)
        min_limit: Minimum joint value
        max_limit: Maximum joint value

    Example::

        # Fixed joint (no motion)
        fixed = Joint("mount", JointType.FIXED)

        # Revolute joint about Z axis
        revolute = Joint(
            name="shoulder",
            joint_type=JointType.REVOLUTE,
            axis=(0, 0, 1),
            min_limit=-90,
            max_limit=90,
        )
        revolute.value = 45.0  # Set current angle

        # Get total transform
        tf = revolute.get_transform()  # base_transform @ rotation(45°)
    """

    name: str
    joint_type: JointType = JointType.FIXED
    base_transform: Transform = field(default_factory=Transform.identity)
    axis: Tuple[float, float, float] = (0.0, 0.0, 1.0)
    value: float = 0.0
    min_limit: float = -180.0
    max_limit: float = 180.0

    def get_transform(self) -> Transform:
        """Compute total joint transform (base + motion).

        :returns: Transform from parent frame to child origin

        For FIXED joints, returns just the base_transform.
        For REVOLUTE joints, adds rotation about axis.
        For PRISMATIC joints, adds translation along axis.
        """
        if self.joint_type == JointType.FIXED:
            return self.base_transform

        elif self.joint_type == JointType.REVOLUTE:
            # Rotation about axis by value degrees
            motion = Transform.from_axis_angle(self.axis, self.value)
            return self.base_transform @ motion

        elif self.joint_type == JointType.PRISMATIC:
            # Translation along axis by value mm
            ax = self.axis
            norm = math.sqrt(ax[0]**2 + ax[1]**2 + ax[2]**2)
            if norm > 1e-10:
                ax = (ax[0]/norm, ax[1]/norm, ax[2]/norm)
            motion = Transform.from_translation(
                ax[0] * self.value,
                ax[1] * self.value,
                ax[2] * self.value,
            )
            return self.base_transform @ motion

        elif self.joint_type == JointType.CYLINDRICAL:
            # Both rotation and translation about same axis
            # value[0] = rotation (deg), value[1] = translation (mm)
            # For simplicity, use value as rotation, assume no translation
            motion = Transform.from_axis_angle(self.axis, self.value)
            return self.base_transform @ motion

        elif self.joint_type == JointType.SPHERICAL:
            # Ball joint - value is ignored (orientation set externally)
            return self.base_transform

        return self.base_transform

    def clamp_value(self, value: float) -> float:
        """Clamp value to joint limits.

        :param value: Unclamped joint value
        :returns: Value clamped to [min_limit, max_limit]
        """
        return max(self.min_limit, min(self.max_limit, value))

    def set_value(self, value: float, clamp: bool = True) -> None:
        """Set joint value with optional clamping.

        :param value: New joint value
        :param clamp: If True, clamp to limits
        """
        if clamp:
            self.value = self.clamp_value(value)
        else:
            self.value = value

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "name": self.name,
            "type": self.joint_type.value,
            "axis": list(self.axis),
            "value": self.value,
            "limits": [self.min_limit, self.max_limit],
            "base_transform": self.base_transform.to_dict(),
        }

    @classmethod
    def from_dict(cls, data: dict) -> Joint:
        """Create joint from dictionary."""
        return cls(
            name=data["name"],
            joint_type=JointType(data.get("type", "fixed")),
            axis=tuple(data.get("axis", [0, 0, 1])),
            value=data.get("value", 0.0),
            min_limit=data.get("limits", [-180, 180])[0],
            max_limit=data.get("limits", [-180, 180])[1],
            base_transform=Transform.from_dict(data["base_transform"])
            if "base_transform" in data else Transform.identity(),
        )


__all__ = ["Joint", "JointType"]
