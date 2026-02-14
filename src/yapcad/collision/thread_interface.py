"""Thread interface for screw thread engagement overlap.

This module defines the ThreadInterface class for representing
screw thread engagement regions where controlled overlap is expected.

Copyright (c) 2026 yapCAD contributors
License: MIT
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Tuple, Any

from .interface import InterfaceVolume, InterfaceType, CompatibilityResult


class ThreadType(Enum):
    """Standard thread types.

    Attributes:
        ISO_METRIC: ISO metric thread (M2, M3, M4, etc.)
        ISO_METRIC_FINE: ISO metric fine pitch thread
        UNC: Unified National Coarse (imperial)
        UNF: Unified National Fine (imperial)
        ACME: Trapezoidal/ACME thread for power transmission
        CUSTOM: User-defined thread profile
    """
    ISO_METRIC = "ISO"
    ISO_METRIC_FINE = "ISO_F"
    UNC = "UNC"
    UNF = "UNF"
    ACME = "ACME"
    CUSTOM = "CUSTOM"


class ThreadClass(Enum):
    """Thread fit classes per ISO 965-1.

    Defines the tolerance grade for thread fits, affecting
    clearance/interference between mating threads.

    Attributes:
        ISO_6H_6g: Standard fit (normal clearance)
        ISO_6G_6h: Close fit (reduced clearance)
        ISO_5H_4h: Tight fit (precision applications)
    """
    ISO_6H_6g = "6H/6g"
    ISO_6G_6h = "6G/6h"
    ISO_5H_4h = "5H/4h"


# Standard ISO metric coarse thread pitches
_METRIC_COARSE_PITCH: Dict[str, float] = {
    "M1.6": 0.35, "M2": 0.4, "M2.5": 0.45, "M3": 0.5,
    "M3.5": 0.6, "M4": 0.7, "M5": 0.8, "M6": 1.0,
    "M8": 1.25, "M10": 1.5, "M12": 1.75, "M14": 2.0,
    "M16": 2.0, "M20": 2.5, "M24": 3.0, "M30": 3.5,
}


@dataclass
class ThreadInterface(InterfaceVolume):
    """Interface volume for screw thread engagement.

    Defines the region where external threads (bolt/screw) engage
    with internal threads (tapped hole/nut). Compatible threads
    must have matching pitch and diameter.

    Attributes:
        pitch: Thread pitch in mm (distance between threads)
        major_diameter: Major (nominal) diameter in mm
        thread_type: Thread standard (ISO_METRIC, UNC, etc.)
        thread_class: Fit class (tolerance grade)
        engagement_length: Length of thread engagement in mm
        is_internal: True for tapped hole/nut, False for bolt/screw

    Derived Properties:
        minor_diameter: Root diameter of external thread
        pitch_diameter: Diameter where thread flanks meet
        thread_depth: Radial depth of thread profile

    Common ISO Metric Threads:
        M2: pitch=0.4mm, major=2.0mm
        M2.5: pitch=0.45mm, major=2.5mm
        M3: pitch=0.5mm, major=3.0mm
        M4: pitch=0.7mm, major=4.0mm
        M5: pitch=0.8mm, major=5.0mm
        M6: pitch=1.0mm, major=6.0mm

    Example:
        >>> bolt = ThreadInterface.from_metric_size("M3", is_internal=False)
        >>> hole = ThreadInterface.from_metric_size("M3", is_internal=True)
        >>> result = bolt.check_compatibility(hole)
        >>> print(result)  # COMPATIBLE: Threads compatible: M3.0x0.5
    """
    pitch: float = 0.5
    major_diameter: float = 3.0
    thread_type: ThreadType = ThreadType.ISO_METRIC
    thread_class: ThreadClass = ThreadClass.ISO_6H_6g
    engagement_length: float = 6.0
    is_internal: bool = False
    interface_type: InterfaceType = field(default=InterfaceType.THREAD, init=False)

    def __post_init__(self):
        """Set default description if not provided."""
        if not self.description:
            kind = "internal" if self.is_internal else "external"
            self.description = f"Thread {kind}: M{self.major_diameter}x{self.pitch}"

    @property
    def minor_diameter(self) -> float:
        """Minor diameter (root of external thread).

        Approximation for ISO metric: D_minor = D_major - 1.0825 * pitch
        """
        return self.major_diameter - 1.0825 * self.pitch

    @property
    def pitch_diameter(self) -> float:
        """Pitch diameter (where thread flanks meet).

        Approximation for ISO metric: D_pitch = D_major - 0.6495 * pitch
        """
        return self.major_diameter - 0.6495 * self.pitch

    @property
    def thread_depth(self) -> float:
        """Thread depth (radial distance from major to minor diameter)."""
        return (self.major_diameter - self.minor_diameter) / 2.0

    def get_bounding_cylinder(self) -> Tuple[float, float]:
        """Get cylindrical bounding volume (radius, height).

        Returns:
            Tuple of (radius, height) where radius is major_diameter/2
            and height is engagement_length.
        """
        radius = self.major_diameter / 2.0
        return (radius, self.engagement_length)

    def get_engagement_depth(self) -> float:
        """Engagement length is the depth for threads.

        Returns:
            Thread engagement length in mm
        """
        return self.engagement_length

    def check_compatibility(self, other: InterfaceVolume) -> CompatibilityResult:
        """Check if this thread can engage with another interface.

        Compatible conditions:
            1. Other must be a ThreadInterface
            2. One must be internal, one external
            3. Major diameters must match (within 0.1% tolerance)
            4. Pitches must match (within 0.1% tolerance)
            5. Thread types should be compatible

        Args:
            other: Another interface volume to check against

        Returns:
            CompatibilityResult indicating if threads can properly engage
        """
        if not isinstance(other, ThreadInterface):
            return CompatibilityResult(
                is_compatible=False,
                reason=f"Cannot thread with {type(other).__name__}"
            )

        warnings = []

        # Check internal/external pairing
        if self.is_internal == other.is_internal:
            kind = "internal" if self.is_internal else "external"
            return CompatibilityResult(
                is_compatible=False,
                reason=f"Both threads are {kind} - need one internal, one external"
            )

        # Check diameter match (0.1% tolerance)
        dia_diff = abs(self.major_diameter - other.major_diameter)
        if dia_diff > self.major_diameter * 0.001:
            return CompatibilityResult(
                is_compatible=False,
                reason=f"Diameter mismatch: M{self.major_diameter} vs M{other.major_diameter}"
            )

        # Check pitch match (0.1% tolerance)
        pitch_diff = abs(self.pitch - other.pitch)
        if pitch_diff > self.pitch * 0.001:
            return CompatibilityResult(
                is_compatible=False,
                reason=f"Pitch mismatch: {self.pitch}mm vs {other.pitch}mm"
            )

        # Check thread type compatibility
        if self.thread_type != other.thread_type:
            # Some cross-compatibility exists (ISO metric variants)
            compatible_pairs = {
                (ThreadType.ISO_METRIC, ThreadType.ISO_METRIC_FINE),
            }
            pair = (self.thread_type, other.thread_type)
            if pair not in compatible_pairs and tuple(reversed(pair)) not in compatible_pairs:
                warnings.append(
                    f"Thread type mismatch: {self.thread_type.value} vs {other.thread_type.value}"
                )

        # Calculate overlap volume (annular thread region)
        engagement = min(self.engagement_length, other.engagement_length)
        inner_r = self.minor_diameter / 2.0
        outer_r = self.major_diameter / 2.0
        overlap_volume = math.pi * (outer_r**2 - inner_r**2) * engagement

        return CompatibilityResult(
            is_compatible=True,
            reason=f"Threads compatible: M{self.major_diameter}x{self.pitch}",
            warnings=warnings,
            overlap_volume=overlap_volume,
            required_clearance=0.0  # Thread fit class defines clearance
        )

    @classmethod
    def from_metric_size(
        cls,
        size: str,
        engagement_length: float = None,
        is_internal: bool = False,
        name: str = None,
        part_name: str = "",
        **kwargs
    ) -> 'ThreadInterface':
        """Create ThreadInterface from metric size string.

        Args:
            size: Metric size string (e.g., "M3", "M2.5", "M4x0.5")
            engagement_length: Thread engagement in mm (default: 1.5 * diameter)
            is_internal: True for tapped hole
            name: Interface name (default: "thread_{size}")
            part_name: Name of the part this interface belongs to
            **kwargs: Additional arguments passed to constructor

        Returns:
            ThreadInterface configured for the specified size

        Raises:
            ValueError: If size string cannot be parsed

        Example:
            >>> bolt = ThreadInterface.from_metric_size("M3", is_internal=False)
            >>> fine_bolt = ThreadInterface.from_metric_size("M4x0.5")
        """
        # Parse size string
        size_upper = size.upper()

        if "X" in size_upper:
            # Fine pitch specified: M4x0.5
            parts = size_upper.split("X")
            major = float(parts[0].replace("M", ""))
            pitch = float(parts[1])
            thread_type = ThreadType.ISO_METRIC_FINE
        else:
            # Coarse pitch from lookup table
            major = float(size_upper.replace("M", ""))
            pitch = _METRIC_COARSE_PITCH.get(size_upper)
            if pitch is None:
                # Interpolate for non-standard sizes
                pitch = 0.5 if major <= 3 else 0.7 if major <= 4 else major * 0.125
            thread_type = ThreadType.ISO_METRIC

        if engagement_length is None:
            engagement_length = 1.5 * major

        if name is None:
            name = f"thread_{size_upper.replace('.', '_')}"

        return cls(
            name=name,
            part_name=part_name,
            pitch=pitch,
            major_diameter=major,
            thread_type=thread_type,
            engagement_length=engagement_length,
            is_internal=is_internal,
            **kwargs
        )

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization.

        Returns:
            Dictionary with all thread parameters
        """
        d = super().to_dict()
        d.update({
            "pitch": self.pitch,
            "major_diameter": self.major_diameter,
            "minor_diameter": self.minor_diameter,
            "pitch_diameter": self.pitch_diameter,
            "thread_type": self.thread_type.value,
            "thread_class": self.thread_class.value,
            "engagement_length": self.engagement_length,
            "is_internal": self.is_internal,
        })
        return d


def create_bolt_pattern_interfaces(
    pattern_name: str,
    bolt_part_name: str,
    hole_part_name: str,
    thread_size: str,
    positions: List[Tuple[float, float, float]],
    engagement_length: float = None
) -> List[ThreadInterface]:
    """Create interface volumes for a bolt pattern.

    Creates matching pairs of external (bolt) and internal (hole)
    thread interfaces at specified positions.

    Args:
        pattern_name: Base name for the bolt pattern
        bolt_part_name: Part containing the bolts
        hole_part_name: Part containing the tapped holes
        thread_size: Thread size (e.g., "M3", "M2.5")
        positions: List of (x, y, z) positions for each bolt
        engagement_length: Thread engagement in mm (default: from size)

    Returns:
        List of ThreadInterface objects (bolts and holes)

    Example:
        >>> positions = [(10, 0, 0), (-10, 0, 0), (0, 10, 0), (0, -10, 0)]
        >>> interfaces = create_bolt_pattern_interfaces(
        ...     "mounting_bolts",
        ...     "BRACKET", "BASE_PLATE",
        ...     "M3", positions
        ... )
    """
    interfaces = []

    for i, (x, y, z) in enumerate(positions):
        # Bolt (external thread)
        bolt = ThreadInterface.from_metric_size(
            thread_size,
            engagement_length=engagement_length,
            is_internal=False,
            name=f"{pattern_name}_bolt_{i+1}",
            part_name=bolt_part_name,
            center=(x, y, z),
        )
        interfaces.append(bolt)

        # Tapped hole (internal thread)
        hole = ThreadInterface.from_metric_size(
            thread_size,
            engagement_length=engagement_length,
            is_internal=True,
            name=f"{pattern_name}_hole_{i+1}",
            part_name=hole_part_name,
            center=(x, y, z),
        )
        interfaces.append(hole)

    return interfaces
