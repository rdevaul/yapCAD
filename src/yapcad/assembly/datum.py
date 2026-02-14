## Datum feature system for yapCAD assembly constraints
## Copyright (c) 2026 yapCAD contributors
## All rights reserved

# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation files
# (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
# BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""Datum feature system for constraint-based assembly in yapCAD.

This module provides datum features (geometric references) and part definitions
that enable declarative assembly constraints. Datums are named geometric features
(points, axes, planes, frames, circles) defined in a part's local coordinate
system that can be used to establish mates and constraints between parts.

Key concepts:
    - **Datum**: A named geometric reference on a part (point, axis, plane, etc.)
    - **PartDefinition**: A part with its named datum features
    - Datums transform with the part when assembled
    - Datums enable declarative mate constraints (FLUSH, CONCENTRIC, etc.)

Example usage::

    from yapcad.assembly.datum import Datum, DatumType, PartDefinition
    from yapcad.geom import point

    # Define a motor with datum features
    motor = PartDefinition(
        name="DDSM115_MOTOR",
        geometry_source="cots/motor.dsl",
        is_printable=False,
        material="Aluminum"
    )

    # Add mounting plane datum
    motor.add_datum(Datum(
        name="stator_face",
        datum_type=DatumType.PLANE,
        origin=point(0, 10, 0),
        normal=point(0, 1, 0, 0),
        description="Stator mounting surface facing +Y"
    ))

    # Add rotation axis datum
    motor.add_datum(Datum(
        name="motor_axis",
        datum_type=DatumType.AXIS,
        origin=point(0, 0, 0),
        direction=point(0, 1, 0, 0),
        description="Motor rotation axis along Y"
    ))

    # Add mounting hole circle datum
    motor.add_datum(Datum(
        name="mounting_holes",
        datum_type=DatumType.CIRCLE,
        origin=point(0, 10, 0),
        normal=point(0, 1, 0, 0),
        radius=15.2,
        description="M2.5 mounting hole pattern at r=15.2mm"
    ))

    # Validate all datums
    issues = motor.validate_datums()
    if issues:
        print("Datum validation issues:", issues)
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Union
from enum import Enum

import yapcad.geom as geom
from yapcad.xform import Matrix


class DatumType(Enum):
    """Types of datum features that can be defined on a part.

    Datum types correspond to common geometric references used in mechanical
    design and assembly constraints:

    - **POINT**: A single point in 3D space
    - **AXIS**: An infinite line defined by origin + direction
    - **PLANE**: An infinite plane defined by origin + normal
    - **FRAME**: A full coordinate frame with origin + 3 orthogonal axes
    - **CIRCLE**: A circle defined by center + normal + radius
    """
    POINT = "point"
    AXIS = "axis"
    PLANE = "plane"
    FRAME = "frame"
    CIRCLE = "circle"


@dataclass
class Datum:
    """A named geometric reference feature on a part.

    Datums are defined in the part's local coordinate system and transform
    with the part when assembled. They provide explicit geometric references
    that can be used to establish assembly mates and constraints.

    All geometric parameters use yapCAD's homogeneous coordinate convention:
        - Points: [x, y, z, 1] (w=1 for positions)
        - Directions/Normals: [x, y, z, 0] (w=0 for direction vectors)

    Attributes:
        name: Unique identifier for this datum within the part
        datum_type: Type of geometric feature (POINT, AXIS, PLANE, etc.)
        origin: Origin point [x, y, z, 1] in part's local coordinates
        direction: Direction vector [x, y, z, 0] for AXIS datum
        normal: Normal vector [x, y, z, 0] for PLANE or CIRCLE datum
        x_axis: X-axis direction [x, y, z, 0] for FRAME datum
        y_axis: Y-axis direction [x, y, z, 0] for FRAME datum
        radius: Radius value (in mm) for CIRCLE datum
        description: Human-readable description of this datum's purpose
        tags: Optional metadata tags for filtering/searching datums

    Example::

        # Mounting plane on motor stator
        stator_face = Datum(
            name="stator_face",
            datum_type=DatumType.PLANE,
            origin=[0, 10, 0, 1],
            normal=[0, 1, 0, 0],
            description="Stator mounting surface facing +Y"
        )

        # Rotation axis through motor shaft
        motor_axis = Datum(
            name="motor_axis",
            datum_type=DatumType.AXIS,
            origin=[0, 0, 0, 1],
            direction=[0, 1, 0, 0],
            description="Motor rotation axis along Y"
        )

        # Mounting hole pattern
        holes = Datum(
            name="mounting_holes",
            datum_type=DatumType.CIRCLE,
            origin=[0, 10, 0, 1],
            normal=[0, 1, 0, 0],
            radius=15.2,
            description="M2.5 mounting holes at 15.2mm radius"
        )
    """
    name: str
    datum_type: DatumType
    origin: List[float] = field(default_factory=lambda: [0.0, 0.0, 0.0, 1.0])

    # For AXIS: direction vector; for PLANE/CIRCLE: normal vector
    direction: Optional[List[float]] = None
    normal: Optional[List[float]] = None

    # For FRAME: additional basis vectors
    x_axis: Optional[List[float]] = None
    y_axis: Optional[List[float]] = None

    # For CIRCLE: radius
    radius: Optional[float] = None

    # Metadata
    description: str = ""
    tags: List[str] = field(default_factory=list)

    def __post_init__(self):
        """Validate datum definition and ensure proper vector format.

        Validates that:
            - Required attributes for each datum type are present
            - Vectors have proper homogeneous coordinates (4 components)
            - Origin is a point (w=1), directions/normals are vectors (w=0)
            - Direction and normal vectors are non-zero

        Raises:
            ValueError: If datum definition is invalid
        """
        # Ensure origin is a proper point vector [x, y, z, 1]
        self.origin = geom.point(
            self.origin[0] if len(self.origin) > 0 else 0.0,
            self.origin[1] if len(self.origin) > 1 else 0.0,
            self.origin[2] if len(self.origin) > 2 else 0.0
        )

        # Validate based on datum type
        if self.datum_type == DatumType.AXIS:
            if self.direction is None:
                raise ValueError(
                    f"Axis datum '{self.name}' requires direction vector"
                )
            # Ensure direction is a proper direction vector [x, y, z, 0]
            self.direction = geom.vect(
                self.direction[0] if len(self.direction) > 0 else 0.0,
                self.direction[1] if len(self.direction) > 1 else 0.0,
                self.direction[2] if len(self.direction) > 2 else 1.0,
                0.0  # w=0 for direction vectors
            )
            # Check for non-zero direction
            if geom.mag(self.direction) < geom.epsilon:
                raise ValueError(
                    f"Axis datum '{self.name}' has zero-length direction"
                )

        elif self.datum_type == DatumType.PLANE:
            if self.normal is None:
                raise ValueError(
                    f"Plane datum '{self.name}' requires normal vector"
                )
            # Ensure normal is a proper direction vector [x, y, z, 0]
            self.normal = geom.vect(
                self.normal[0] if len(self.normal) > 0 else 0.0,
                self.normal[1] if len(self.normal) > 1 else 0.0,
                self.normal[2] if len(self.normal) > 2 else 1.0,
                0.0  # w=0 for direction vectors
            )
            # Check for non-zero normal
            if geom.mag(self.normal) < geom.epsilon:
                raise ValueError(
                    f"Plane datum '{self.name}' has zero-length normal"
                )

        elif self.datum_type == DatumType.CIRCLE:
            if self.normal is None or self.radius is None:
                raise ValueError(
                    f"Circle datum '{self.name}' requires normal and radius"
                )
            # Ensure normal is a proper direction vector [x, y, z, 0]
            self.normal = geom.vect(
                self.normal[0] if len(self.normal) > 0 else 0.0,
                self.normal[1] if len(self.normal) > 1 else 0.0,
                self.normal[2] if len(self.normal) > 2 else 1.0,
                0.0  # w=0 for direction vectors
            )
            # Check for non-zero normal
            if geom.mag(self.normal) < geom.epsilon:
                raise ValueError(
                    f"Circle datum '{self.name}' has zero-length normal"
                )
            # Validate radius
            if not geom.isgoodnum(self.radius) or self.radius <= 0:
                raise ValueError(
                    f"Circle datum '{self.name}' requires positive radius"
                )

        elif self.datum_type == DatumType.FRAME:
            if self.x_axis is None or self.y_axis is None:
                raise ValueError(
                    f"Frame datum '{self.name}' requires x_axis and y_axis"
                )
            # Ensure axes are proper direction vectors [x, y, z, 0]
            self.x_axis = geom.vect(
                self.x_axis[0] if len(self.x_axis) > 0 else 1.0,
                self.x_axis[1] if len(self.x_axis) > 1 else 0.0,
                self.x_axis[2] if len(self.x_axis) > 2 else 0.0,
                0.0  # w=0 for direction vectors
            )
            self.y_axis = geom.vect(
                self.y_axis[0] if len(self.y_axis) > 0 else 0.0,
                self.y_axis[1] if len(self.y_axis) > 1 else 1.0,
                self.y_axis[2] if len(self.y_axis) > 2 else 0.0,
                0.0  # w=0 for direction vectors
            )
            # Check for non-zero axes
            if geom.mag(self.x_axis) < geom.epsilon:
                raise ValueError(
                    f"Frame datum '{self.name}' has zero-length x_axis"
                )
            if geom.mag(self.y_axis) < geom.epsilon:
                raise ValueError(
                    f"Frame datum '{self.name}' has zero-length y_axis"
                )

    def transform(self, matrix: Union[Matrix, List[List[float]]]) -> 'Datum':
        """Return a new Datum transformed by the given 4x4 matrix.

        Transforms the datum's geometric features by applying the given
        transformation matrix. Points are translated, direction vectors
        are rotated (but not translated), and the radius is preserved
        (assuming uniform scaling).

        Args:
            matrix: 4x4 transformation matrix, either as yapCAD Matrix object
                   or as a list of 4 rows (each row is a 4-element list)

        Returns:
            A new Datum with transformed geometry

        Example::

            from yapcad.xform import Translation
            from yapcad.geom import point

            # Original datum at origin
            datum = Datum("mount", DatumType.POINT, origin=point(0, 0, 0))

            # Transform by translation
            T = Translation(point(10, 20, 30))
            new_datum = datum.transform(T)
            # new_datum.origin is now [10, 20, 30, 1]
        """
        # Convert to Matrix if needed
        if not isinstance(matrix, Matrix):
            matrix = Matrix(matrix)

        # Helper to transform a point (w=1)
        def transform_point(p):
            result = matrix.mul(p)
            return result

        # Helper to transform a direction vector (w=0)
        def transform_vector(v):
            if v is None:
                return None
            result = matrix.mul(v)
            # Normalize direction vectors
            mag = geom.mag(result)
            if mag > geom.epsilon:
                result = geom.scale3(result, 1.0 / mag)
            return result

        # Create new datum with transformed geometry
        new_datum = Datum(
            name=self.name,
            datum_type=self.datum_type,
            origin=transform_point(self.origin),
            direction=transform_vector(self.direction),
            normal=transform_vector(self.normal),
            x_axis=transform_vector(self.x_axis),
            y_axis=transform_vector(self.y_axis),
            radius=self.radius,  # Radius preserved (assumes uniform scale)
            description=self.description,
            tags=self.tags.copy()
        )
        return new_datum


@dataclass
class PartDefinition:
    """Definition of a part with its datum features.

    A PartDefinition extends the concept of a part to include explicit
    geometric datum features that can be used for assembly constraints.
    This enables declarative mate definitions (FLUSH, CONCENTRIC, etc.)
    that reference named datums rather than requiring manual transform
    calculations.

    Attributes:
        name: Unique identifier for this part
        datums: Dictionary mapping datum names to Datum objects
        geometry_source: Path to STL file, DSL command, or other geometry
        is_printable: Whether this part is 3D printed (vs. COTS)
        material: Material specification (e.g., "PETG", "Aluminum")
        description: Human-readable description of part function

    Example::

        from yapcad.assembly.datum import PartDefinition, Datum, DatumType
        from yapcad.geom import point

        # Define a wheel bracket part
        bracket = PartDefinition(
            name="WHEEL_BRACKET",
            geometry_source="parts/bracket.dsl",
            is_printable=True,
            material="PETG",
            description="Mounting bracket for drive wheel motor"
        )

        # Add datum for motor mounting face
        bracket.add_datum(Datum(
            name="motor_interface",
            datum_type=DatumType.PLANE,
            origin=point(0, 0, 0),
            normal=point(0, 0, 1, 0),
            description="Motor mounting face normal to +Z"
        ))

        # Add datum for motor bore axis
        bracket.add_datum(Datum(
            name="bore_axis",
            datum_type=DatumType.AXIS,
            origin=point(0, 0, 0),
            direction=point(0, 0, 1, 0),
            description="Motor shaft bore along Z axis"
        ))

        # Validate all datums
        issues = bracket.validate_datums()
        if issues:
            print("Warning:", issues)
    """
    name: str
    datums: Dict[str, Datum] = field(default_factory=dict)

    # Link to geometry (STL file, DSL command, or solid)
    geometry_source: Optional[str] = None

    # Metadata
    is_printable: bool = True
    material: str = "PETG"
    description: str = ""

    def add_datum(self, datum: Datum) -> 'PartDefinition':
        """Add a datum feature to this part.

        Args:
            datum: The Datum object to add

        Returns:
            Self (for method chaining)

        Raises:
            ValueError: If a datum with the same name already exists

        Example::

            bracket = PartDefinition("BRACKET")
            bracket.add_datum(Datum("face", DatumType.PLANE, ...))
            bracket.add_datum(Datum("axis", DatumType.AXIS, ...))
        """
        if datum.name in self.datums:
            raise ValueError(
                f"Datum '{datum.name}' already exists on part '{self.name}'"
            )
        self.datums[datum.name] = datum
        return self

    def get_datum(self, name: str) -> Datum:
        """Get a datum by name.

        Args:
            name: Name of the datum to retrieve

        Returns:
            The Datum object

        Raises:
            KeyError: If no datum with that name exists

        Example::

            bracket = PartDefinition("BRACKET")
            bracket.add_datum(Datum("mount_face", ...))
            face = bracket.get_datum("mount_face")
        """
        if name not in self.datums:
            available = ", ".join(self.datums.keys()) if self.datums else "none"
            raise KeyError(
                f"Datum '{name}' not found on part '{self.name}'. "
                f"Available datums: {available}"
            )
        return self.datums[name]

    def validate_datums(self) -> List[str]:
        """Validate all datum definitions, returning list of issues.

        Performs validation checks on all datums to detect potential errors:
            - Duplicate definitions (same origin and type)
            - Overlapping features that might indicate copy-paste errors

        Returns:
            List of warning/error messages (empty if all valid)

        Example::

            bracket = PartDefinition("BRACKET")
            bracket.add_datum(Datum("face1", DatumType.PLANE, origin=[0,0,0], ...))
            bracket.add_datum(Datum("face2", DatumType.PLANE, origin=[0,0,0], ...))

            issues = bracket.validate_datums()
            # issues will contain warning about duplicate origin
        """
        issues = []

        # Check for duplicate origins (might indicate copy-paste error)
        for name, datum in self.datums.items():
            for other_name, other in self.datums.items():
                if name >= other_name:  # Only check each pair once
                    continue

                # Check if origins are coincident
                if geom.dist(datum.origin, other.origin) < geom.epsilon:
                    # If same type and same origin, likely duplicate
                    if datum.datum_type == other.datum_type:
                        issues.append(
                            f"Datums '{name}' and '{other_name}' have "
                            f"identical type ({datum.datum_type.value}) "
                            f"and origin {datum.origin[:3]}"
                        )

                    # If both are axes/planes, check if parallel/coplanar
                    if datum.datum_type in (DatumType.AXIS, DatumType.PLANE):
                        if other.datum_type in (DatumType.AXIS, DatumType.PLANE):
                            d1 = datum.direction if datum.direction else datum.normal
                            d2 = other.direction if other.direction else other.normal
                            if d1 and d2:
                                # Check if parallel (dot product near 1)
                                dot = abs(geom.dot(d1, d2))
                                if dot > (1.0 - geom.epsilon):
                                    issues.append(
                                        f"Datums '{name}' and '{other_name}' "
                                        f"have same origin and parallel directions"
                                    )

        return issues
