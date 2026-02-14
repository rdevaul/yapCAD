"""Central datum registry for cross-DSL datum linking.

This module provides a unified interface to access datums from various sources:
- DSL files (via sidecar JSON or metadata extraction)
- COTS surrogate JSON files (for commercial off-the-shelf parts)
- Programmatically defined datums

The registry enables face-to-face mate computation by providing a single
point of access for datum definitions regardless of their source.

Example usage::

    from yapcad.assembly.datum_registry import DatumRegistry

    # Get a datum from a COTS surrogate file
    stator_face = DatumRegistry.get_datum(
        "cots/xh430_surrogate.json",
        "stator_mounting_face"
    )

    # Get a datum by part name (searches registered sources)
    servo_axis = DatumRegistry.get_datum(
        "AXIS3_SERVO_XH430",
        "rotation_axis"
    )

Copyright (c) 2026 yapCAD contributors
License: MIT
"""

from __future__ import annotations

import json
import os
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Any, Callable

from .datum import Datum, DatumType, PartDefinition

logger = logging.getLogger(__name__)


@dataclass
class DatumSource:
    """Represents a source of datum definitions.

    Attributes:
        source_id: Unique identifier for this source (path or name)
        source_type: Type of source ("json", "dsl", "programmatic")
        file_path: Absolute path to source file (if file-based)
        datums: Dictionary of datum name -> Datum object
        part_names: Set of part names defined in this source
        metadata: Additional source metadata
    """
    source_id: str
    source_type: str
    file_path: Optional[str] = None
    datums: Dict[str, Datum] = field(default_factory=dict)
    part_names: Set[str] = field(default_factory=set)
    metadata: Dict[str, Any] = field(default_factory=dict)


class DatumRegistry:
    """Central registry for datums from DSL files and COTS surrogates.

    The registry provides a unified interface to access datums regardless
    of their source (DSL, STL metadata, surrogate JSON files).

    This is implemented as a singleton to ensure consistent state across
    the application.

    Class Methods:
        get_datum(source, name): Get a datum by source and name
        register_source(source_id, datums): Register a source with datums
        register_datum(source_id, datum): Add a single datum to a source
        list_sources(): List all registered sources
        list_datums(source): List all datums in a source
        clear(): Clear all registered datums (for testing)

    Example::

        # Load from surrogate JSON
        stator = DatumRegistry.get_datum(
            "cots/xh430_surrogate.json",
            "stator_mounting_face"
        )

        # Register programmatic datums
        DatumRegistry.register_datum(
            "LINK_2_3",
            Datum("servo_mount_face", DatumType.PLANE, ...)
        )
    """

    _instance: Optional['DatumRegistry'] = None
    _sources: Dict[str, DatumSource] = {}
    _loaded_sources: Set[str] = set()
    _search_paths: List[str] = []

    @classmethod
    def _get_instance(cls) -> 'DatumRegistry':
        """Get the singleton instance."""
        if cls._instance is None:
            cls._instance = cls()
            cls._instance._sources = {}
            cls._instance._loaded_sources = set()
            cls._instance._search_paths = []
        return cls._instance

    @classmethod
    def get_datum(cls, source: str, name: str) -> Datum:
        """Get a datum by source and name.

        Args:
            source: Source identifier - can be:
                    - DSL file path: "scara_arm/scara_arm.dsl"
                    - Surrogate JSON: "cots/xh430_surrogate.json"
                    - Part name: "AXIS3_SERVO_XH430", "LINK_2_3"
            name: Datum name within the source
                    e.g., "servo_mount_face", "stator_face"

        Returns:
            The Datum object

        Raises:
            KeyError: If datum not found
        """
        instance = cls._get_instance()

        # Lazy load source if needed
        if source not in instance._loaded_sources:
            instance._load_source(source)

        if source not in instance._sources:
            raise KeyError(f"Datum source '{source}' not found")
        if name not in instance._sources[source].datums:
            available = list(instance._sources[source].datums.keys())
            raise KeyError(f"Datum '{name}' not found in '{source}'. "
                          f"Available: {available}")

        return instance._sources[source].datums[name]

    @classmethod
    def get_datums_for_source(cls, source: str) -> Dict[str, Datum]:
        """Get all datums from a source.

        Args:
            source: Source identifier

        Returns:
            Dictionary of datum name -> Datum object
        """
        instance = cls._get_instance()

        if source not in instance._loaded_sources:
            instance._load_source(source)

        if source not in instance._sources:
            raise KeyError(f"Datum source '{source}' not found")

        return instance._sources[source].datums.copy()

    @classmethod
    def register_source(cls, source_id: str,
                       datums: Dict[str, Datum],
                       source_type: str = "programmatic",
                       metadata: Optional[Dict[str, Any]] = None) -> None:
        """Register a source with its datums.

        Args:
            source_id: Unique identifier for this source
            datums: Dictionary of datum name -> Datum object
            source_type: Type of source ("json", "dsl", "programmatic")
            metadata: Optional metadata about the source
        """
        instance = cls._get_instance()

        source = DatumSource(
            source_id=source_id,
            source_type=source_type,
            datums=datums.copy(),
            part_names={source_id},
            metadata=metadata or {}
        )
        instance._sources[source_id] = source
        instance._loaded_sources.add(source_id)

    @classmethod
    def register_datum(cls, source_id: str, datum: Datum) -> None:
        """Add a single datum to a source.

        Creates the source if it doesn't exist.

        Args:
            source_id: Source identifier
            datum: The Datum to register
        """
        instance = cls._get_instance()

        if source_id not in instance._sources:
            instance._sources[source_id] = DatumSource(
                source_id=source_id,
                source_type="programmatic",
                datums={},
                part_names={source_id}
            )
            instance._loaded_sources.add(source_id)

        instance._sources[source_id].datums[datum.name] = datum

    @classmethod
    def register_part(cls, part: PartDefinition) -> None:
        """Register a PartDefinition and all its datums.

        Args:
            part: PartDefinition with datum features
        """
        instance = cls._get_instance()

        datums = {d.name: d for d in part.datums}
        cls.register_source(
            source_id=part.name,
            datums=datums,
            source_type="programmatic",
            metadata={"part_name": part.name}
        )

    @classmethod
    def list_sources(cls) -> List[str]:
        """List all registered source identifiers."""
        instance = cls._get_instance()
        return list(instance._sources.keys())

    @classmethod
    def list_datums(cls, source: str) -> List[str]:
        """List all datum names in a source.

        Args:
            source: Source identifier

        Returns:
            List of datum names
        """
        instance = cls._get_instance()

        if source not in instance._loaded_sources:
            instance._load_source(source)

        if source not in instance._sources:
            return []

        return list(instance._sources[source].datums.keys())

    @classmethod
    def add_search_path(cls, path: str) -> None:
        """Add a directory to search for datum source files.

        Args:
            path: Directory path to add to search paths
        """
        instance = cls._get_instance()
        abs_path = os.path.abspath(path)
        if abs_path not in instance._search_paths:
            instance._search_paths.append(abs_path)

    @classmethod
    def clear(cls) -> None:
        """Clear all registered datums (useful for testing)."""
        instance = cls._get_instance()
        instance._sources.clear()
        instance._loaded_sources.clear()

    # =========================================================================
    # Internal loading methods
    # =========================================================================

    def _load_source(self, source: str) -> None:
        """Load datums from a source file.

        Attempts to load from various source types in order:
        1. JSON surrogate file
        2. DSL file (via sidecar JSON)
        3. Part name lookup in already-loaded sources
        """
        # Try to find the file
        file_path = self._find_source_file(source)

        if file_path:
            if file_path.endswith('.json'):
                self._load_surrogate_json(source, file_path)
            elif file_path.endswith('.dsl'):
                self._load_dsl_datums(source, file_path)
            else:
                logger.warning(f"Unknown source file type: {file_path}")
        else:
            # Try as a part name in already-loaded sources
            self._try_part_name_lookup(source)

        self._loaded_sources.add(source)

    def _find_source_file(self, source: str) -> Optional[str]:
        """Find the source file in search paths.

        Args:
            source: Source identifier (may be relative path or filename)

        Returns:
            Absolute path if found, None otherwise
        """
        # If it's already an absolute path and exists
        if os.path.isabs(source) and os.path.exists(source):
            return source

        # Search in search paths
        for search_path in self._search_paths:
            candidate = os.path.join(search_path, source)
            if os.path.exists(candidate):
                return candidate

        # Try current directory
        if os.path.exists(source):
            return os.path.abspath(source)

        return None

    def _load_surrogate_json(self, source_id: str, json_path: str) -> None:
        """Load datums from a COTS surrogate JSON file.

        JSON structure expected:
        {
            "part_name": "XH430-W350",
            "datums": {
                "datum_name": {
                    "type": "PLANE|AXIS|CIRCLE|POINT|FRAME",
                    "origin": [x, y, z],
                    "normal": [nx, ny, nz],  // for PLANE/CIRCLE
                    "direction": [dx, dy, dz],  // for AXIS
                    "radius": r,  // for CIRCLE
                    "description": "..."
                }
            }
        }
        """
        try:
            with open(json_path, 'r') as f:
                data = json.load(f)
        except (json.JSONDecodeError, IOError) as e:
            logger.error(f"Failed to load surrogate JSON '{json_path}': {e}")
            self._sources[source_id] = DatumSource(
                source_id=source_id,
                source_type="json",
                file_path=json_path,
                datums={},
                metadata={"error": str(e)}
            )
            return

        datums = {}
        for name, datum_data in data.get('datums', {}).items():
            try:
                datum = self._datum_from_dict(name, datum_data)
                datums[name] = datum
            except (KeyError, ValueError) as e:
                logger.warning(f"Failed to parse datum '{name}' from "
                             f"'{json_path}': {e}")

        # Also register by part_name if present
        part_names = {source_id}
        if 'part_name' in data:
            part_names.add(data['part_name'])

        source = DatumSource(
            source_id=source_id,
            source_type="json",
            file_path=json_path,
            datums=datums,
            part_names=part_names,
            metadata={
                'part_name': data.get('part_name'),
                'manufacturer': data.get('manufacturer'),
                'description': data.get('description'),
                'stl_file': data.get('stl_file'),
                'body_dimensions': data.get('body_dimensions')
            }
        )
        self._sources[source_id] = source

        # Also register by part name for easy lookup
        if 'part_name' in data and data['part_name'] != source_id:
            self._sources[data['part_name']] = source

    def _load_dsl_datums(self, source_id: str, dsl_path: str) -> None:
        """Load datums from a DSL file via sidecar JSON.

        DSL files can export datum metadata to a sidecar JSON file
        with the same name but .datums.json extension.
        """
        # Check for sidecar JSON
        sidecar_path = dsl_path.replace('.dsl', '.datums.json')
        if os.path.exists(sidecar_path):
            self._load_surrogate_json(source_id, sidecar_path)
        else:
            # DSL file exists but no datum sidecar
            logger.info(f"DSL file '{dsl_path}' has no datum sidecar "
                       f"(expected: {sidecar_path})")
            self._sources[source_id] = DatumSource(
                source_id=source_id,
                source_type="dsl",
                file_path=dsl_path,
                datums={},
                metadata={"no_datums": True}
            )

    def _try_part_name_lookup(self, source: str) -> None:
        """Try to find datums by part name in already-loaded sources."""
        for src_id, src in self._sources.items():
            if source in src.part_names:
                # Found it - create an alias
                self._sources[source] = src
                return

        # Not found - create empty source
        self._sources[source] = DatumSource(
            source_id=source,
            source_type="not_found",
            datums={}
        )

    def _datum_from_dict(self, name: str, data: Dict[str, Any]) -> Datum:
        """Create a Datum from a dictionary representation.

        Args:
            name: Datum name
            data: Dictionary with datum properties

        Returns:
            Datum object
        """
        type_str = data.get('type', 'POINT').upper()
        try:
            datum_type = DatumType(type_str.lower())
        except ValueError:
            # Try as enum name
            datum_type = DatumType[type_str]

        # Convert 3-element arrays to 4-element homogeneous coords
        origin = data.get('origin', [0, 0, 0])
        if len(origin) == 3:
            origin = list(origin) + [1]  # w=1 for point

        normal = data.get('normal')
        if normal and len(normal) == 3:
            normal = list(normal) + [0]  # w=0 for direction

        direction = data.get('direction')
        if direction and len(direction) == 3:
            direction = list(direction) + [0]  # w=0 for direction

        x_axis = data.get('x_axis')
        if x_axis and len(x_axis) == 3:
            x_axis = list(x_axis) + [0]

        y_axis = data.get('y_axis')
        if y_axis and len(y_axis) == 3:
            y_axis = list(y_axis) + [0]

        return Datum(
            name=name,
            datum_type=datum_type,
            origin=origin,
            normal=normal,
            direction=direction,
            x_axis=x_axis,
            y_axis=y_axis,
            radius=data.get('radius'),
            description=data.get('description', '')
        )


# =============================================================================
# Utility functions
# =============================================================================

def datum_to_transform_matrix(datum: Datum) -> Optional[Any]:
    """Convert a Datum to a 4x4 transformation matrix.

    For FRAME datums, constructs the full transform from axes.
    For PLANE datums, constructs a frame with normal as Z-axis.
    For AXIS datums, constructs a frame with direction as Z-axis.

    Args:
        datum: The datum to convert

    Returns:
        4x4 numpy array (if numpy available), or None
    """
    try:
        import numpy as np
    except ImportError:
        return None

    origin = np.array(datum.origin[:3])

    if datum.datum_type == DatumType.FRAME:
        # Full frame with explicit axes
        x_axis = np.array(datum.x_axis[:3]) if datum.x_axis else np.array([1, 0, 0])
        y_axis = np.array(datum.y_axis[:3]) if datum.y_axis else np.array([0, 1, 0])
        z_axis = np.cross(x_axis, y_axis)
        z_axis = z_axis / np.linalg.norm(z_axis)

    elif datum.datum_type == DatumType.PLANE:
        # Normal becomes Z-axis, derive X and Y
        z_axis = np.array(datum.normal[:3])
        z_axis = z_axis / np.linalg.norm(z_axis)

        # Find a non-parallel reference vector
        if abs(z_axis[0]) < 0.9:
            ref = np.array([1, 0, 0])
        else:
            ref = np.array([0, 1, 0])

        x_axis = np.cross(z_axis, ref)
        x_axis = x_axis / np.linalg.norm(x_axis)
        y_axis = np.cross(z_axis, x_axis)

    elif datum.datum_type == DatumType.AXIS:
        # Direction becomes Z-axis
        z_axis = np.array(datum.direction[:3])
        z_axis = z_axis / np.linalg.norm(z_axis)

        # Find a non-parallel reference vector
        if abs(z_axis[0]) < 0.9:
            ref = np.array([1, 0, 0])
        else:
            ref = np.array([0, 1, 0])

        x_axis = np.cross(z_axis, ref)
        x_axis = x_axis / np.linalg.norm(x_axis)
        y_axis = np.cross(z_axis, x_axis)

    elif datum.datum_type == DatumType.CIRCLE:
        # Normal becomes Z-axis (same as PLANE)
        z_axis = np.array(datum.normal[:3])
        z_axis = z_axis / np.linalg.norm(z_axis)

        if abs(z_axis[0]) < 0.9:
            ref = np.array([1, 0, 0])
        else:
            ref = np.array([0, 1, 0])

        x_axis = np.cross(z_axis, ref)
        x_axis = x_axis / np.linalg.norm(x_axis)
        y_axis = np.cross(z_axis, x_axis)

    else:  # POINT
        # Identity orientation at origin
        x_axis = np.array([1, 0, 0])
        y_axis = np.array([0, 1, 0])
        z_axis = np.array([0, 0, 1])

    # Build 4x4 transform matrix
    T = np.eye(4)
    T[:3, 0] = x_axis
    T[:3, 1] = y_axis
    T[:3, 2] = z_axis
    T[:3, 3] = origin

    return T
