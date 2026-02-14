"""Kinematic chain (tree of parts) with world transform computation.

Copyright (c) 2026 yapCAD contributors
License: MIT
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Any

from .transform import Transform
from .joint import Joint, JointType
from .frame import CoordinateFrame
from .part import KinematicPart

logger = logging.getLogger(__name__)


class KinematicChain:
    """Tree of kinematic parts with world transform computation.

    A KinematicChain manages a hierarchy of parts connected by joints.
    It provides efficient world transform computation with caching and
    supports JSON import/export.

    The chain has a virtual root node called GLOBAL_DATUM at the world origin.
    All root parts (parent=None) are implicitly attached to GLOBAL_DATUM.

    Attributes:
        name: Name of this kinematic chain/assembly
        parts: Dict mapping part names to KinematicPart objects
        DATUM_NAME: Name of the virtual root node ("GLOBAL_DATUM")

    Example::

        chain = KinematicChain("robot_arm")

        # Add parts
        chain.add_part(KinematicPart(name="BASE"))
        chain.add_part(KinematicPart(
            name="LINK1",
            parent="BASE",
            joint=Joint("j1", JointType.REVOLUTE),
        ))

        # Set joint values
        chain.set_joint_value("LINK1", 45.0)

        # Get world transforms
        world_tf = chain.get_world_transform("LINK1")

        # Export to JSON
        chain.export_json("positions.json")
    """

    DATUM_NAME = "GLOBAL_DATUM"

    def __init__(self, name: str = "assembly"):
        """Initialize kinematic chain.

        :param name: Name for this chain/assembly
        """
        self.name = name
        self.parts: Dict[str, KinematicPart] = {}
        self._children: Dict[str, List[str]] = {}

        # Create global datum as virtual root
        self._create_global_datum()

    def _create_global_datum(self) -> None:
        """Create the virtual global datum root node."""
        datum = KinematicPart(
            name=self.DATUM_NAME,
            parent=None,
            joint=Joint("global_joint", JointType.FIXED),
            description="Global world origin",
        )
        datum.set_cached_world_transform(Transform.identity())
        self.parts[self.DATUM_NAME] = datum
        self._children[self.DATUM_NAME] = []

    def add_part(self, part: KinematicPart) -> KinematicPart:
        """Add a part to the kinematic chain.

        :param part: Part to add
        :returns: The added part (for chaining)
        :raises ValueError: If part name already exists or parent not found

        Parts with parent=None are attached to GLOBAL_DATUM.
        """
        if part.name in self.parts:
            raise ValueError(f"Part '{part.name}' already exists in chain")

        # Attach orphan parts to global datum
        if part.parent is None:
            part.parent = self.DATUM_NAME

        # Validate parent exists
        if part.parent not in self.parts:
            raise ValueError(
                f"Parent '{part.parent}' not found for part '{part.name}'"
            )

        # Validate parent frame exists
        parent = self.parts[part.parent]
        if part.parent_frame not in parent.frames:
            logger.warning(
                f"Parent frame '{part.parent_frame}' not found on '{part.parent}', "
                f"using ORIGIN"
            )
            part.parent_frame = "ORIGIN"

        # Add to parts dict and children list
        self.parts[part.name] = part
        if part.parent not in self._children:
            self._children[part.parent] = []
        self._children[part.parent].append(part.name)

        if part.name not in self._children:
            self._children[part.name] = []

        # Invalidate transform cache
        part.invalidate_cache()

        return part

    def set_joint_value(self, part_name: str, value: float, clamp: bool = True) -> None:
        """Set joint value for a part and invalidate descendant caches.

        :param part_name: Part whose joint to modify
        :param value: New joint value (degrees for revolute, mm for prismatic)
        :param clamp: If True, clamp to joint limits
        """
        if part_name not in self.parts:
            raise ValueError(f"Part '{part_name}' not found")

        part = self.parts[part_name]
        part.joint.set_value(value, clamp)

        # Invalidate this part and all descendants
        self._invalidate_subtree(part_name)

    def get_joint_value(self, part_name: str) -> float:
        """Get current joint value for a part.

        :param part_name: Part to query
        :returns: Current joint value
        """
        if part_name not in self.parts:
            raise ValueError(f"Part '{part_name}' not found")
        return self.parts[part_name].joint.value

    def get_world_transform(
        self,
        part_name: str,
        frame_name: str = "ORIGIN",
    ) -> Transform:
        """Compute world transform for a part frame.

        :param part_name: Part to get transform for
        :param frame_name: Frame on the part (default "ORIGIN")
        :returns: World transform (from global origin to part frame)

        This method uses caching for efficiency. Transform is recomputed
        only when a joint value changes in the ancestor chain.
        """
        if part_name not in self.parts:
            raise ValueError(f"Part '{part_name}' not found")

        part = self.parts[part_name]

        # Compute part world transform if not cached
        if not part.is_cache_valid():
            if part.parent is None or part.name == self.DATUM_NAME:
                # Root node - identity
                part.set_cached_world_transform(Transform.identity())
            else:
                # Recursive: get parent world transform at attachment frame
                parent_world = self.get_world_transform(
                    part.parent, part.parent_frame
                )
                # Apply joint transform
                joint_tf = part.joint.get_transform()
                world_tf = parent_world @ joint_tf
                part.set_cached_world_transform(world_tf)

        # Get part world transform and apply frame offset
        part_world = part.get_cached_world_transform()
        if part_world is None:
            part_world = Transform.identity()

        if frame_name == "ORIGIN":
            return part_world

        frame = part.get_frame(frame_name)
        if frame is None:
            logger.warning(f"Frame '{frame_name}' not found on '{part_name}'")
            return part_world

        return part_world @ frame.get_transform()

    def get_relative_transform(
        self,
        from_part: str,
        to_part: str,
        from_frame: str = "ORIGIN",
        to_frame: str = "ORIGIN",
    ) -> Transform:
        """Compute transform from one part/frame to another.

        :param from_part: Source part name
        :param to_part: Target part name
        :param from_frame: Frame on source part
        :param to_frame: Frame on target part
        :returns: Transform from source to target
        """
        from_world = self.get_world_transform(from_part, from_frame)
        to_world = self.get_world_transform(to_part, to_frame)
        return from_world.inverse() @ to_world

    def get_chain_to_root(self, part_name: str) -> List[str]:
        """Get ancestor chain from part to root.

        :param part_name: Starting part
        :returns: List of part names [part, parent, grandparent, ..., GLOBAL_DATUM]
        """
        chain = []
        current = part_name

        while current and current in self.parts:
            chain.append(current)
            current = self.parts[current].parent

        return chain

    def get_children(self, part_name: str) -> List[str]:
        """Get direct children of a part.

        :param part_name: Parent part name
        :returns: List of child part names
        """
        return self._children.get(part_name, [])

    def get_all_descendants(self, part_name: str) -> List[str]:
        """Get all descendants of a part (recursive).

        :param part_name: Root part
        :returns: List of all descendant part names
        """
        descendants = []
        for child in self.get_children(part_name):
            descendants.append(child)
            descendants.extend(self.get_all_descendants(child))
        return descendants

    def get_all_world_transforms(self) -> Dict[str, Transform]:
        """Compute world transforms for all parts.

        :returns: Dict mapping part names to world transforms
        """
        return {
            name: self.get_world_transform(name)
            for name in self.parts
            if name != self.DATUM_NAME
        }

    def validate(self) -> List[str]:
        """Validate the kinematic chain.

        :returns: List of validation issues (empty if valid)

        Checks for:
            - Orphaned parts (parent not found)
            - Missing parent frames
            - Circular dependencies
        """
        issues = []

        for name, part in self.parts.items():
            if name == self.DATUM_NAME:
                continue

            # Check parent exists
            if part.parent and part.parent not in self.parts:
                issues.append(f"Part '{name}' has invalid parent '{part.parent}'")

            # Check parent frame exists
            if part.parent and part.parent in self.parts:
                parent = self.parts[part.parent]
                if part.parent_frame not in parent.frames:
                    issues.append(
                        f"Part '{name}' references missing frame "
                        f"'{part.parent_frame}' on '{part.parent}'"
                    )

        # Check for cycles (simple check via chain length)
        for name in self.parts:
            chain = self.get_chain_to_root(name)
            if len(chain) > len(self.parts):
                issues.append(f"Circular dependency detected involving '{name}'")
                break

        return issues

    def _invalidate_subtree(self, part_name: str) -> None:
        """Recursively invalidate transform caches for a subtree."""
        if part_name in self.parts:
            self.parts[part_name].invalidate_cache()
            for child in self.get_children(part_name):
                self._invalidate_subtree(child)

    # =========================================================================
    # JSON Import/Export
    # =========================================================================

    def export_json(
        self,
        filepath: str,
        include_world_transforms: bool = True,
    ) -> None:
        """Export chain to JSON file.

        :param filepath: Output file path
        :param include_world_transforms: If True, include computed world transforms
        """
        data = self.to_dict(include_world_transforms)
        Path(filepath).parent.mkdir(parents=True, exist_ok=True)
        with open(filepath, 'w') as f:
            json.dump(data, f, indent=2)

    def to_dict(self, include_world_transforms: bool = True) -> dict:
        """Convert chain to dictionary.

        :param include_world_transforms: If True, include computed world transforms
        :returns: JSON-serializable dictionary
        """
        parts_data = {}

        for name, part in self.parts.items():
            if name == self.DATUM_NAME:
                continue

            part_data = part.to_dict()

            if include_world_transforms:
                world_tf = self.get_world_transform(name)
                part_data["world_transform"] = world_tf.to_dict()

            parts_data[name] = part_data

        return {
            "name": self.name,
            "parts": parts_data,
        }

    @classmethod
    def from_json(cls, filepath: str) -> KinematicChain:
        """Load chain from JSON file.

        :param filepath: Input file path
        :returns: Loaded KinematicChain
        """
        with open(filepath, 'r') as f:
            data = json.load(f)
        return cls.from_dict(data)

    @classmethod
    def from_dict(cls, data: dict) -> KinematicChain:
        """Create chain from dictionary.

        :param data: Dictionary with chain data
        :returns: New KinematicChain
        """
        chain = cls(name=data.get("name", "assembly"))

        # Sort parts by dependency (parents first)
        parts_data = data.get("parts", {})
        sorted_names = cls._topological_sort(parts_data)

        for name in sorted_names:
            if name in parts_data:
                part = KinematicPart.from_dict(parts_data[name])
                chain.add_part(part)

        return chain

    @staticmethod
    def _topological_sort(parts_data: dict) -> List[str]:
        """Sort part names so parents come before children."""
        # Build dependency graph
        dependencies: Dict[str, str] = {}
        for name, data in parts_data.items():
            parent = data.get("parent")
            if parent and parent in parts_data:
                dependencies[name] = parent
            else:
                dependencies[name] = None

        # Simple topological sort
        sorted_names = []
        remaining = set(dependencies.keys())

        while remaining:
            # Find parts with no remaining dependencies
            ready = [
                name for name in remaining
                if dependencies[name] is None or dependencies[name] not in remaining
            ]
            if not ready:
                # Cycle or missing parent - add remaining in order
                sorted_names.extend(remaining)
                break
            sorted_names.extend(ready)
            remaining -= set(ready)

        return sorted_names

    def print_tree(self, indent: int = 0) -> None:
        """Print tree structure for debugging."""
        def _print_subtree(name: str, depth: int):
            prefix = "  " * depth
            part = self.parts[name]
            jtype = part.joint.joint_type.value
            print(f"{prefix}{name} [{jtype}]")
            for child in self.get_children(name):
                _print_subtree(child, depth + 1)

        print(f"KinematicChain: {self.name}")
        _print_subtree(self.DATUM_NAME, indent)


__all__ = ["KinematicChain"]
