"""
Runtime value wrappers for DSL interpreter.

Values wrap Python/yapCAD objects with DSL type metadata for runtime type
checking and provenance tracking.
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Union
from enum import Enum

from ..types import (
    Type, PrimitiveType, GeometricPrimitiveType, CurveType,
    CompoundCurveType, SurfaceType, SolidType, ListType, DictType,
    OptionalTypeWrapper, FunctionType,
    INT, FLOAT, BOOL, STRING,
    POINT, POINT2D, POINT3D, VECTOR, VECTOR2D, VECTOR3D, TRANSFORM,
    SOLID, REGION2D, SURFACE, SHELL,
)


@dataclass
class Value:
    """
    A runtime value with DSL type information.

    The `data` field holds the actual Python/yapCAD object.
    The `type` field holds the DSL type for runtime validation.
    """
    data: Any
    type: Type

    def __repr__(self) -> str:
        return f"Value({self.data!r}, {self.type})"

    def is_truthy(self) -> bool:
        """Check if this value is truthy in boolean context."""
        if self.type == BOOL:
            return bool(self.data)
        # Lists are truthy if non-empty
        if isinstance(self.type, ListType):
            return len(self.data) > 0
        # Optional is truthy if not None
        if isinstance(self.type, OptionalTypeWrapper):
            return self.data is not None
        # Everything else is truthy (numbers, geometry, etc.)
        return True


# Convenience constructors for primitive values

def int_val(n: int) -> Value:
    """Create an integer value."""
    return Value(int(n), INT)


def float_val(x: float) -> Value:
    """Create a float value."""
    return Value(float(x), FLOAT)


def bool_val(b: bool) -> Value:
    """Create a boolean value."""
    return Value(bool(b), BOOL)


def string_val(s: str) -> Value:
    """Create a string value."""
    return Value(str(s), STRING)


def list_val(items: List[Value], element_type: Type) -> Value:
    """Create a list value from a list of Values."""
    return Value([v.data for v in items], ListType(element_type))


def dict_val(items: Dict[str, Value]) -> Value:
    """Create a dict value."""
    return Value({k: v.data for k, v in items.items()}, DictType())


def none_val(inner_type: Type) -> Value:
    """Create a None value for an optional type."""
    return Value(None, OptionalTypeWrapper(inner_type))


# Geometry value constructors

def point_val(coords: Any, is_2d: bool = False) -> Value:
    """Create a point value (2D or 3D)."""
    return Value(coords, POINT2D if is_2d else POINT3D)


def vector_val(coords: Any, is_2d: bool = False) -> Value:
    """Create a vector value (2D or 3D)."""
    return Value(coords, VECTOR2D if is_2d else VECTOR3D)


def transform_val(matrix: Any) -> Value:
    """Create a transform value."""
    return Value(matrix, TRANSFORM)


def solid_val(solid_data: Any) -> Value:
    """Create a solid value."""
    return Value(solid_data, SOLID)


def region2d_val(region_data: Any) -> Value:
    """Create a region2d value."""
    return Value(region_data, REGION2D)


def surface_val(surface_data: Any) -> Value:
    """Create a surface value."""
    return Value(surface_data, SURFACE)


def shell_val(shell_data: Any) -> Value:
    """Create a shell value."""
    return Value(shell_data, SHELL)


def path3d_val(path_data: Any) -> Value:
    """Create a path3d value (3D curve/spine for sweeps)."""
    from ..types import PATH3D
    return Value(path_data, PATH3D)


# Type checking utilities

def unwrap_value(v: Value) -> Any:
    """Extract the raw Python/yapCAD data from a Value."""
    return v.data


def unwrap_values(values: List[Value]) -> List[Any]:
    """Extract raw data from a list of Values."""
    return [v.data for v in values]


def wrap_value(data: Any, dsl_type: Type) -> Value:
    """Wrap raw data with a DSL type."""
    return Value(data, dsl_type)


def check_type(value: Value, expected: Type) -> bool:
    """Check if a value's type is assignable to expected type."""
    return expected.is_assignable_from(value.type)


def coerce_numeric(value: Value) -> Value:
    """Coerce int to float if needed."""
    if value.type == INT:
        return float_val(float(value.data))
    return value


@dataclass
class EmitResult:
    """
    Result of an emit statement - the final output of a command.

    Contains the emitted geometry plus any metadata attached via `with {...}`.
    """
    value: Value
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def data(self) -> Any:
        """Get the raw geometry data."""
        return self.value.data

    @property
    def type(self) -> Type:
        """Get the DSL type."""
        return self.value.type


@dataclass
class RequireFailure:
    """
    Represents a failed require constraint.

    Used to collect all require failures during execution rather than
    stopping at the first one.
    """
    message: str
    expression_text: Optional[str] = None

    def __str__(self) -> str:
        if self.expression_text:
            return f"require failed: {self.message} (expression: {self.expression_text})"
        return f"require failed: {self.message}"
