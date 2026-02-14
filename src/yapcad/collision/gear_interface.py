"""Gear mesh interface for involute gear teeth overlap.

This module defines the GearMeshInterface class for representing
gear teeth meshing regions where controlled overlap is expected.

Copyright (c) 2026 yapCAD contributors
License: MIT
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Any

from .interface import InterfaceVolume, InterfaceType, CompatibilityResult


@dataclass
class GearMeshInterface(InterfaceVolume):
    """Interface volume for involute gear teeth meshing.

    Defines the region where gear teeth overlap during meshing.
    Two gears are compatible if they have the same module and
    pressure angle (standard involute tooth profile).

    Attributes:
        module: Gear module in mm (pitch_diameter / teeth)
        teeth: Number of teeth on the gear
        pressure_angle: Pressure angle in degrees (typically 20.0)
        face_width: Width of gear face in mm (axial engagement length)
        helix_angle: Helix angle in degrees for helical gears (0 for spur)
        backlash: Design backlash in mm (0.0 for theoretical)

    Derived Properties:
        pitch_diameter: module * teeth
        addendum: module (tooth height above pitch circle)
        dedendum: 1.25 * module (tooth depth below pitch circle)
        outside_diameter: pitch_diameter + 2 * addendum
        root_diameter: pitch_diameter - 2 * dedendum
        base_diameter: pitch_diameter * cos(pressure_angle)

    Example:
        >>> sun = GearMeshInterface(
        ...     name="sun_teeth", part_name="SUN_GEAR",
        ...     module=0.75, teeth=18, pressure_angle=20.0, face_width=10.0
        ... )
        >>> planet = GearMeshInterface(
        ...     name="planet_teeth", part_name="PLANET_GEAR",
        ...     module=0.75, teeth=36, pressure_angle=20.0, face_width=10.0
        ... )
        >>> result = sun.check_compatibility(planet)
        >>> print(result)  # COMPATIBLE: Gears compatible: m=0.75mm, PA=20.0deg
        >>> print(f"Center distance: {sun.get_center_distance(planet):.2f}mm")

    Notes:
        - Module must match exactly between meshing gears (within 0.1% tolerance)
        - Pressure angle must match (within 0.5 degree tolerance)
        - For helical gears, helix angles must be opposite (LH + RH = 0)
        - Face widths must have sufficient overlap for load transfer
    """
    module: float = 1.0
    teeth: int = 20
    pressure_angle: float = 20.0
    face_width: float = 10.0
    helix_angle: float = 0.0
    backlash: float = 0.0
    interface_type: InterfaceType = field(default=InterfaceType.GEAR_MESH, init=False)

    def __post_init__(self):
        """Set default description if not provided."""
        if not self.description:
            self.description = f"Gear mesh: {self.teeth}T, m={self.module}mm"

    @property
    def pitch_diameter(self) -> float:
        """Pitch diameter = module * teeth."""
        return self.module * self.teeth

    @property
    def addendum(self) -> float:
        """Addendum (tooth height above pitch circle) = module."""
        return self.module

    @property
    def dedendum(self) -> float:
        """Dedendum (tooth depth below pitch circle) = 1.25 * module."""
        return 1.25 * self.module

    @property
    def outside_diameter(self) -> float:
        """Outside diameter (tip of teeth)."""
        return self.pitch_diameter + 2 * self.addendum

    @property
    def root_diameter(self) -> float:
        """Root diameter (bottom of teeth)."""
        return self.pitch_diameter - 2 * self.dedendum

    @property
    def base_diameter(self) -> float:
        """Base diameter for involute profile."""
        return self.pitch_diameter * math.cos(math.radians(self.pressure_angle))

    def get_bounding_cylinder(self) -> Tuple[float, float]:
        """Get cylindrical bounding volume (radius, height).

        Returns:
            Tuple of (radius, height) where radius extends to outside diameter
            and height is the face width.
        """
        radius = self.outside_diameter / 2.0
        return (radius, self.face_width)

    def get_engagement_depth(self) -> float:
        """Face width is the engagement depth for gears.

        Returns:
            Face width in mm
        """
        return self.face_width

    def get_mesh_overlap_depth(self) -> float:
        """Calculate the theoretical tooth mesh overlap depth.

        This is the radial depth where teeth from two meshing gears
        actually overlap, used for collision solver exceptions.

        Returns:
            Depth of tooth engagement zone in mm (approximately 2 * module)
        """
        # Engagement depth is approximately addendum + dedendum - clearance
        # For standard gears, this is about 2 * module
        return 2.0 * self.module

    def get_center_distance(self, other: 'GearMeshInterface') -> float:
        """Calculate theoretical center distance for meshing with another gear.

        Args:
            other: Another gear to mesh with

        Returns:
            Center-to-center distance for proper meshing in mm

        Raises:
            ValueError: If modules don't match
        """
        if abs(self.module - other.module) > self.module * 0.001:
            raise ValueError(
                f"Cannot calculate center distance for different modules: "
                f"{self.module}mm vs {other.module}mm"
            )

        # For external mesh: (D1 + D2) / 2
        # For internal mesh: abs(D1 - D2) / 2 (where D1 > D2)
        return (self.pitch_diameter + other.pitch_diameter) / 2.0

    def get_expected_overlap_volume(self, other: 'GearMeshInterface') -> float:
        """Calculate expected overlap volume when meshing with another gear.

        This provides an estimate of how much volume the gear teeth
        should overlap when properly meshed, useful for validating
        collision detection results.

        Args:
            other: Another gear to mesh with

        Returns:
            Estimated overlap volume in mm^3
        """
        mesh_depth = self.get_mesh_overlap_depth()
        min_face = min(self.face_width, other.face_width)

        # Approximate overlap as arc segment at pitch circle
        # This is a rough approximation - actual overlap is complex
        return math.pi * mesh_depth * (self.pitch_diameter / 2) * min_face * 0.1

    def check_compatibility(self, other: InterfaceVolume) -> CompatibilityResult:
        """Check if this gear can mesh with another interface.

        Compatible conditions:
            1. Other must be a GearMeshInterface
            2. Modules must match (within 0.1% tolerance)
            3. Pressure angles must match (within 0.5 degree tolerance)
            4. Helix angles must be compatible (opposite for external mesh)
            5. Face widths must have sufficient overlap

        Args:
            other: Another interface volume to check against

        Returns:
            CompatibilityResult indicating if gears can properly mesh
        """
        if not isinstance(other, GearMeshInterface):
            return CompatibilityResult(
                is_compatible=False,
                reason=f"Cannot mesh gear with {type(other).__name__}"
            )

        warnings = []

        # Check module match (0.1% tolerance for manufacturing)
        module_diff = abs(self.module - other.module)
        if module_diff > self.module * 0.001:
            return CompatibilityResult(
                is_compatible=False,
                reason=f"Module mismatch: {self.module}mm vs {other.module}mm"
            )

        # Check pressure angle match (0.5 degree tolerance)
        pa_diff = abs(self.pressure_angle - other.pressure_angle)
        if pa_diff > 0.5:
            return CompatibilityResult(
                is_compatible=False,
                reason=f"Pressure angle mismatch: {self.pressure_angle}deg vs {other.pressure_angle}deg"
            )

        # Check helix angle compatibility
        # For helical gears, one must be LH and other RH for external mesh
        if self.helix_angle != 0 or other.helix_angle != 0:
            if abs(self.helix_angle + other.helix_angle) > 0.5:
                warnings.append(
                    f"Helix angle sum is {self.helix_angle + other.helix_angle}deg "
                    f"(expected 0deg for external mesh)"
                )

        # Check face width overlap
        min_face = min(self.face_width, other.face_width)
        if min_face < self.module:
            warnings.append(
                f"Face width {min_face}mm is less than module {self.module}mm"
            )

        # Calculate expected overlap volume
        overlap_volume = self.get_expected_overlap_volume(other)

        return CompatibilityResult(
            is_compatible=True,
            reason=f"Gears compatible: m={self.module}mm, PA={self.pressure_angle}deg",
            warnings=warnings,
            overlap_volume=overlap_volume,
            required_clearance=self.backlash
        )

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization.

        Returns:
            Dictionary with all gear parameters
        """
        d = super().to_dict()
        d.update({
            "module": self.module,
            "teeth": self.teeth,
            "pressure_angle": self.pressure_angle,
            "face_width": self.face_width,
            "helix_angle": self.helix_angle,
            "backlash": self.backlash,
            "pitch_diameter": self.pitch_diameter,
            "outside_diameter": self.outside_diameter,
            "root_diameter": self.root_diameter,
            "base_diameter": self.base_diameter,
        })
        return d


def check_gear_mesh_collision(
    gear_a: GearMeshInterface,
    gear_b: GearMeshInterface,
    actual_center_distance: float,
    tolerance: float = 0.1
) -> CompatibilityResult:
    """Check if two gears are correctly positioned for meshing.

    This function verifies both gear compatibility and proper positioning
    by comparing actual center distance to theoretical.

    Args:
        gear_a: First gear interface
        gear_b: Second gear interface
        actual_center_distance: Measured center-to-center distance in mm
        tolerance: Allowable deviation from theoretical in mm

    Returns:
        CompatibilityResult with mesh status and any warnings

    Example:
        >>> sun = GearMeshInterface(name="sun", part_name="SUN", module=0.75, teeth=18)
        >>> planet = GearMeshInterface(name="planet", part_name="PLANET", module=0.75, teeth=36)
        >>> result = check_gear_mesh_collision(sun, planet, actual_center_distance=20.25)
        >>> if result.is_compatible:
        ...     print("Gears properly meshed")
    """
    # First check basic compatibility
    compat = gear_a.check_compatibility(gear_b)
    if not compat.is_compatible:
        return compat

    # Check center distance
    theoretical_cd = gear_a.get_center_distance(gear_b)
    cd_error = abs(actual_center_distance - theoretical_cd)

    if cd_error > tolerance:
        return CompatibilityResult(
            is_compatible=False,
            reason=f"Center distance error: {cd_error:.3f}mm (max {tolerance}mm)",
            warnings=compat.warnings
        )

    # Add center distance info to warnings
    warnings = compat.warnings.copy()
    warnings.append(
        f"Center distance: {actual_center_distance:.3f}mm "
        f"(theoretical: {theoretical_cd:.3f}mm)"
    )

    # Adjust required clearance based on center distance deviation
    clearance_adjustment = theoretical_cd - actual_center_distance

    return CompatibilityResult(
        is_compatible=True,
        reason=compat.reason,
        warnings=warnings,
        overlap_volume=compat.overlap_volume,
        required_clearance=compat.required_clearance + clearance_adjustment
    )


def create_planetary_gearbox_interfaces(
    gearbox_name: str,
    module: float,
    sun_teeth: int,
    planet_teeth: int,
    ring_teeth: int,
    face_width: float,
    num_planets: int = 3,
    pressure_angle: float = 20.0,
    mesh_plane_z: float = 0.0
) -> List[GearMeshInterface]:
    """Create interface volumes for a complete planetary gearbox.

    This factory function creates GearMeshInterface objects for all gears
    in a planetary gearbox configuration: one sun, one ring, and multiple
    planet gears at the correct orbital positions.

    Args:
        gearbox_name: Base name for the gearbox (e.g., "AXIS1")
        module: Gear module in mm
        sun_teeth: Number of teeth on sun gear
        planet_teeth: Number of teeth on planet gears
        ring_teeth: Number of teeth on ring gear (internal)
        face_width: Gear face width in mm
        num_planets: Number of planet gears (typically 3)
        pressure_angle: Pressure angle in degrees (typically 20.0)
        mesh_plane_z: Z coordinate of the gear mesh plane

    Returns:
        List of GearMeshInterface objects for all gears

    Example:
        >>> interfaces = create_planetary_gearbox_interfaces(
        ...     gearbox_name="AXIS1",
        ...     module=0.75,
        ...     sun_teeth=18,
        ...     planet_teeth=36,
        ...     ring_teeth=90,
        ...     face_width=10.0,
        ...     num_planets=3
        ... )
        >>> for iface in interfaces:
        ...     print(f"{iface.name}: {iface.teeth}T at {iface.center}")

    Notes:
        - Verifies planetary gear geometry: ring_teeth = sun_teeth + 2 * planet_teeth
        - Places planets at equal angular intervals
        - All gears centered at Z = mesh_plane_z
    """
    interfaces = []

    # Verify planetary geometry constraint
    expected_ring = sun_teeth + 2 * planet_teeth
    if ring_teeth != expected_ring:
        import warnings
        warnings.warn(
            f"Ring teeth ({ring_teeth}) != sun_teeth + 2*planet_teeth ({expected_ring}). "
            f"This may not be a valid planetary configuration."
        )

    # Sun gear interface (at center)
    sun_interface = GearMeshInterface(
        name=f"{gearbox_name}_sun_teeth",
        part_name=f"{gearbox_name}_SUN_GEAR",
        module=module,
        teeth=sun_teeth,
        pressure_angle=pressure_angle,
        face_width=face_width,
        center=(0.0, 0.0, mesh_plane_z),
        description=f"Sun gear: {sun_teeth}T, m={module}mm"
    )
    interfaces.append(sun_interface)

    # Ring gear interface (internal gear, same center as sun)
    ring_interface = GearMeshInterface(
        name=f"{gearbox_name}_ring_teeth",
        part_name=f"{gearbox_name}_RING_HOUSING",
        module=module,
        teeth=ring_teeth,
        pressure_angle=pressure_angle,
        face_width=face_width,
        center=(0.0, 0.0, mesh_plane_z),
        description=f"Internal ring gear: {ring_teeth}T, m={module}mm"
    )
    interfaces.append(ring_interface)

    # Planet gear interfaces (at orbital positions)
    orbit_radius = (sun_teeth + planet_teeth) * module / 2.0

    for i in range(num_planets):
        angle_rad = 2 * math.pi * i / num_planets
        cx = orbit_radius * math.cos(angle_rad)
        cy = orbit_radius * math.sin(angle_rad)

        planet_interface = GearMeshInterface(
            name=f"{gearbox_name}_planet_{i+1}_teeth",
            part_name=f"{gearbox_name}_PLANET_GEAR_{i+1}",
            module=module,
            teeth=planet_teeth,
            pressure_angle=pressure_angle,
            face_width=face_width,
            center=(cx, cy, mesh_plane_z),
            description=f"Planet gear {i+1}: {planet_teeth}T at {math.degrees(angle_rad):.0f}deg"
        )
        interfaces.append(planet_interface)

    return interfaces
