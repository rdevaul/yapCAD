"""Coordinate frames for kinematic parts.

Copyright (c) 2026 yapCAD contributors
License: MIT
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Any

from .transform import Transform


@dataclass
class CoordinateFrame:
    """Named reference frame on a kinematic part.

    A coordinate frame defines a local coordinate system on a part
    that can be used for attaching child parts or defining features.

    Frames can be defined via:
        - Explicit transform (transform field)
        - External datum lookup (datum_source + datum_name)

    Attributes:
        name: Unique identifier for this frame on the part
        transform: Explicit transform from part origin to frame
        datum_source: Optional external source ID for datum lookup
        datum_name: Optional datum name for external lookup
        description: Human-readable description

    Example::

        # Explicit frame
        mount = CoordinateFrame(
            name="SERVO_MOUNT",
            transform=Transform.from_translation(0, 0, 50),
            description="Servo mounting face"
        )

        # Frame from external datum (lazy lookup)
        from_datum = CoordinateFrame(
            name="HORN_FACE",
            datum_source="cots/servo.json",
            datum_name="output_shaft",
        )
    """

    name: str
    transform: Optional[Transform] = None
    datum_source: Optional[str] = None
    datum_name: Optional[str] = None
    description: str = ""

    # Internal cache
    _cached_transform: Optional[Transform] = field(
        default=None, repr=False, compare=False
    )
    _cache_valid: bool = field(default=False, repr=False, compare=False)

    def get_transform(self) -> Transform:
        """Get frame transform (with caching for external lookups).

        :returns: Transform from part origin to this frame

        If transform is set explicitly, returns it directly.
        If datum_source/datum_name are set, attempts external lookup.
        Falls back to identity transform if lookup fails.
        """
        # Return explicit transform if set
        if self.transform is not None:
            return self.transform

        # Return cached value if valid
        if self._cache_valid and self._cached_transform is not None:
            return self._cached_transform

        # Attempt external datum lookup
        if self.datum_source and self.datum_name:
            tf = self._lookup_datum_transform()
            if tf is not None:
                object.__setattr__(self, '_cached_transform', tf)
                object.__setattr__(self, '_cache_valid', True)
                return tf

        # Fallback to identity
        return Transform.identity()

    def _lookup_datum_transform(self) -> Optional[Transform]:
        """Attempt to lookup transform from external datum registry.

        :returns: Transform if found, None otherwise
        """
        try:
            from yapcad.assembly.datum_registry import DatumRegistry
            from yapcad.assembly.face_mate import datum_to_transform_matrix

            datum = DatumRegistry.get_datum(self.datum_source, self.datum_name)
            if datum is not None:
                matrix = datum_to_transform_matrix(datum)
                return Transform.from_matrix(matrix)
        except ImportError:
            pass  # Assembly module not available
        except Exception:
            pass  # Datum not found or other error

        return None

    def invalidate_cache(self) -> None:
        """Invalidate cached transform (forces re-lookup)."""
        object.__setattr__(self, '_cache_valid', False)
        object.__setattr__(self, '_cached_transform', None)

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        tf = self.get_transform()
        return {
            "name": self.name,
            "transform": tf.to_dict(),
            "datum_source": self.datum_source,
            "datum_name": self.datum_name,
            "description": self.description,
        }

    @classmethod
    def from_dict(cls, data: dict) -> CoordinateFrame:
        """Create frame from dictionary."""
        transform = None
        if "transform" in data:
            transform = Transform.from_dict(data["transform"])

        return cls(
            name=data["name"],
            transform=transform,
            datum_source=data.get("datum_source"),
            datum_name=data.get("datum_name"),
            description=data.get("description", ""),
        )


__all__ = ["CoordinateFrame"]
