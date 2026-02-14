"""Interface registry for managing allowed overlap regions.

This module defines the InterfaceRegistry class for tracking and
querying interface volumes across an assembly.

Copyright (c) 2026 yapCAD contributors
License: MIT
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple, Any

from .interface import InterfaceVolume, InterfaceType, CompatibilityResult

# Try to import numpy for overlap calculations
try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None


class InterfaceRegistry:
    """Registry for managing interface volumes in an assembly.

    The registry tracks all interface volumes and provides lookup
    methods for the collision detector to find compatible interfaces
    when overlap is detected between parts.

    Key Features:
        - Register/unregister interface volumes by name
        - Query interfaces by part name or interface type
        - Check overlap compatibility between parts
        - Find all compatible interfaces for a given interface

    Usage:
        >>> registry = InterfaceRegistry()
        >>>
        >>> # Register gear interfaces
        >>> registry.register(GearMeshInterface(
        ...     name="sun_teeth", part_name="SUN_GEAR",
        ...     module=0.75, teeth=18, pressure_angle=20.0, face_width=10.0
        ... ))
        >>> registry.register(GearMeshInterface(
        ...     name="planet_teeth", part_name="PLANET_GEAR",
        ...     module=0.75, teeth=36, pressure_angle=20.0, face_width=10.0
        ... ))
        >>>
        >>> # During collision detection
        >>> is_compatible, results = registry.check_overlap_compatibility(
        ...     "SUN_GEAR", "PLANET_GEAR"
        ... )
        >>> if is_compatible:
        ...     print("Overlap is expected (gear mesh)")

    Notes:
        - Interface names must be unique within the registry
        - Parts can have multiple interfaces (e.g., gear teeth + bearing seat)
        - The registry does not store geometry, only interface definitions
    """

    def __init__(self):
        """Initialize an empty interface registry."""
        self._interfaces: Dict[str, InterfaceVolume] = {}
        self._by_part: Dict[str, List[str]] = {}
        self._by_type: Dict[InterfaceType, List[str]] = {}

    def register(self, interface: InterfaceVolume) -> None:
        """Register an interface volume.

        Args:
            interface: InterfaceVolume to register

        Raises:
            ValueError: If interface with same name already registered
        """
        if interface.name in self._interfaces:
            raise ValueError(f"Interface '{interface.name}' already registered")

        self._interfaces[interface.name] = interface

        # Index by part
        if interface.part_name not in self._by_part:
            self._by_part[interface.part_name] = []
        self._by_part[interface.part_name].append(interface.name)

        # Index by type
        if interface.interface_type not in self._by_type:
            self._by_type[interface.interface_type] = []
        self._by_type[interface.interface_type].append(interface.name)

    def register_many(self, interfaces: List[InterfaceVolume]) -> None:
        """Register multiple interface volumes.

        Args:
            interfaces: List of InterfaceVolume objects to register
        """
        for interface in interfaces:
            self.register(interface)

    def unregister(self, name: str) -> Optional[InterfaceVolume]:
        """Remove an interface from the registry.

        Args:
            name: Interface name to remove

        Returns:
            The removed interface, or None if not found
        """
        if name not in self._interfaces:
            return None

        interface = self._interfaces.pop(name)

        # Remove from part index
        if interface.part_name in self._by_part:
            self._by_part[interface.part_name].remove(name)
            if not self._by_part[interface.part_name]:
                del self._by_part[interface.part_name]

        # Remove from type index
        if interface.interface_type in self._by_type:
            self._by_type[interface.interface_type].remove(name)
            if not self._by_type[interface.interface_type]:
                del self._by_type[interface.interface_type]

        return interface

    def get(self, name: str) -> Optional[InterfaceVolume]:
        """Get interface by name.

        Args:
            name: Interface name to look up

        Returns:
            InterfaceVolume if found, None otherwise
        """
        return self._interfaces.get(name)

    def get_interfaces_for_part(self, part_name: str) -> List[InterfaceVolume]:
        """Get all interfaces belonging to a part.

        Args:
            part_name: Name of the part

        Returns:
            List of InterfaceVolume objects for the part (may be empty)
        """
        names = self._by_part.get(part_name, [])
        return [self._interfaces[n] for n in names]

    def get_interfaces_by_type(self, interface_type: InterfaceType) -> List[InterfaceVolume]:
        """Get all interfaces of a specific type.

        Args:
            interface_type: Type of interface to find

        Returns:
            List of matching InterfaceVolume objects (may be empty)
        """
        names = self._by_type.get(interface_type, [])
        return [self._interfaces[n] for n in names]

    def has_interfaces(self, part_name: str) -> bool:
        """Check if a part has any registered interfaces.

        Args:
            part_name: Name of the part

        Returns:
            True if part has one or more registered interfaces
        """
        return part_name in self._by_part and len(self._by_part[part_name]) > 0

    def find_compatible_interfaces(
        self,
        interface: InterfaceVolume,
        candidates: List[InterfaceVolume] = None
    ) -> List[Tuple[InterfaceVolume, CompatibilityResult]]:
        """Find all interfaces compatible with the given interface.

        Args:
            interface: Interface to find matches for
            candidates: List of candidates to check (default: all registered)

        Returns:
            List of (interface, compatibility_result) tuples for compatible interfaces
        """
        if candidates is None:
            candidates = list(self._interfaces.values())

        compatible = []
        for other in candidates:
            if other.name == interface.name:
                continue

            result = interface.check_compatibility(other)
            if result.is_compatible:
                compatible.append((other, result))

        return compatible

    def check_overlap_compatibility(
        self,
        part_a: str,
        part_b: str,
        overlap_region: Optional[Tuple[Any, Any]] = None
    ) -> Tuple[bool, List[CompatibilityResult]]:
        """Check if overlap between two parts is due to compatible interfaces.

        This is the main entry point for the collision detector to check
        whether detected overlap should be allowed (expected interface)
        or flagged as a collision error.

        Args:
            part_a: First part name
            part_b: Second part name
            overlap_region: Optional (center, size) of overlap bounding box
                for spatial filtering (not yet implemented)

        Returns:
            Tuple of:
                - bool: True if ALL overlapping interfaces are compatible
                - List[CompatibilityResult]: Results for each interface pair checked

        Example:
            >>> is_ok, results = registry.check_overlap_compatibility("SUN", "PLANET")
            >>> if is_ok:
            ...     print("Expected gear mesh overlap")
            ... else:
            ...     for r in results:
            ...         if not r.is_compatible:
            ...             print(f"Incompatible: {r.reason}")
        """
        interfaces_a = self.get_interfaces_for_part(part_a)
        interfaces_b = self.get_interfaces_for_part(part_b)

        if not interfaces_a or not interfaces_b:
            # No interfaces defined - overlap is a collision
            return (False, [])

        results = []
        found_compatible = False
        found_incompatible = False

        for iface_a in interfaces_a:
            for iface_b in interfaces_b:
                # Check spatial overlap between interfaces
                if not iface_a.overlaps_with(iface_b):
                    continue

                # Check compatibility
                result = iface_a.check_compatibility(iface_b)
                results.append(result)

                if result.is_compatible:
                    found_compatible = True
                else:
                    found_incompatible = True

        # If any overlapping interfaces are compatible and none are incompatible,
        # the overlap is allowed
        if found_compatible and not found_incompatible:
            return (True, results)

        # If we found both compatible and incompatible, or only incompatible,
        # treat as collision
        return (False, results)

    def get_compatible_interface_names(
        self,
        part_a: str,
        part_b: str
    ) -> List[str]:
        """Get names of compatible interfaces between two parts.

        Convenience method to get just the interface names for
        parts that have compatible overlap.

        Args:
            part_a: First part name
            part_b: Second part name

        Returns:
            List of compatible interface names (from both parts)
        """
        is_compatible, results = self.check_overlap_compatibility(part_a, part_b)

        if not is_compatible:
            return []

        # Collect interface names from compatible results
        names = set()
        interfaces_a = self.get_interfaces_for_part(part_a)
        interfaces_b = self.get_interfaces_for_part(part_b)

        for iface_a in interfaces_a:
            for iface_b in interfaces_b:
                if iface_a.overlaps_with(iface_b):
                    result = iface_a.check_compatibility(iface_b)
                    if result.is_compatible:
                        names.add(iface_a.name)
                        names.add(iface_b.name)

        return list(names)

    def get_all_interfaces(self) -> List[InterfaceVolume]:
        """Get all registered interfaces.

        Returns:
            List of all InterfaceVolume objects in the registry
        """
        return list(self._interfaces.values())

    def get_all_parts(self) -> List[str]:
        """Get names of all parts with registered interfaces.

        Returns:
            List of part names
        """
        return list(self._by_part.keys())

    def clear(self) -> None:
        """Remove all registered interfaces."""
        self._interfaces.clear()
        self._by_part.clear()
        self._by_type.clear()

    def to_dict(self) -> Dict[str, Any]:
        """Convert registry to dictionary for serialization.

        Returns:
            Dictionary with all interfaces serialized
        """
        return {
            "interfaces": [iface.to_dict() for iface in self._interfaces.values()]
        }

    def __len__(self) -> int:
        """Number of registered interfaces."""
        return len(self._interfaces)

    def __contains__(self, name: str) -> bool:
        """Check if interface name is registered."""
        return name in self._interfaces

    def __iter__(self):
        """Iterate over all registered interfaces."""
        return iter(self._interfaces.values())

    def __repr__(self) -> str:
        """String representation for debugging."""
        return (
            f"InterfaceRegistry({len(self)} interfaces, "
            f"{len(self._by_part)} parts, "
            f"{len(self._by_type)} types)"
        )
