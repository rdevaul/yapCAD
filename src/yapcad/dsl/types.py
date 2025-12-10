"""
Type system definitions for the yapCAD DSL.

The type system is organized into five tiers:
    Tier 1: Primitives (int, float, bool, string, point, vector, transform)
    Tier 2: Curve Primitives (line_segment, arc, circle, bezier, nurbs, etc.)
    Tier 3: Compound Curves (path2d, path3d, profile2d, region2d, loop3d)
    Tier 4: Surfaces (surface, shell)
    Tier 5: Solids (solid)

Plus generic types: list<T>, dict
"""

from dataclasses import dataclass, field
from typing import Optional, List, Dict, Set, Union, Tuple
from enum import Enum, auto
from abc import ABC, abstractmethod


class TypeTier(Enum):
    """Tier classification for types in the hierarchy."""
    PRIMITIVE = 1      # int, float, bool, string, point, vector, transform
    CURVE = 2          # line_segment, arc, circle, etc.
    COMPOUND_CURVE = 3 # path2d, path3d, profile2d, region2d, loop3d
    SURFACE = 4        # surface, shell
    SOLID = 5          # solid
    GENERIC = 0        # list<T>, dict


# =============================================================================
# Type Classes
# =============================================================================

@dataclass(frozen=True)
class Type(ABC):
    """Base class for all DSL types."""

    @property
    @abstractmethod
    def name(self) -> str:
        """The type name for display/errors."""
        pass

    @property
    @abstractmethod
    def tier(self) -> TypeTier:
        """The tier this type belongs to."""
        pass

    def is_assignable_from(self, other: "Type") -> bool:
        """Check if this type can accept a value of the other type."""
        return self == other

    def __str__(self) -> str:
        return self.name


@dataclass(frozen=True)
class PrimitiveType(Type):
    """A primitive type (int, float, bool, string)."""
    _name: str

    @property
    def name(self) -> str:
        return self._name

    @property
    def tier(self) -> TypeTier:
        return TypeTier.PRIMITIVE

    def is_assignable_from(self, other: "Type") -> bool:
        if self == other:
            return True
        # int can be promoted to float
        if self._name == "float" and isinstance(other, PrimitiveType) and other._name == "int":
            return True
        return False


@dataclass(frozen=True)
class GeometricPrimitiveType(Type):
    """
    A geometric primitive type (point, vector, transform).

    Points and vectors are dimensionally polymorphic - they can be 2D or 3D.
    point2d/point3d and vector2d/vector3d are specific variants.
    """
    _name: str
    dimension: Optional[int] = None  # None = polymorphic, 2 = 2D, 3 = 3D

    @property
    def name(self) -> str:
        return self._name

    @property
    def tier(self) -> TypeTier:
        return TypeTier.PRIMITIVE

    def is_assignable_from(self, other: "Type") -> bool:
        if self == other:
            return True

        if not isinstance(other, GeometricPrimitiveType):
            return False

        # Polymorphic point/vector can accept specific variants
        if self._name == "point" and other._name in ("point2d", "point3d", "point"):
            return True
        if self._name == "vector" and other._name in ("vector2d", "vector3d", "vector"):
            return True

        # Specific variants only accept same or polymorphic
        if self._name == "point2d" and other._name in ("point2d", "point"):
            return True
        if self._name == "point3d" and other._name in ("point3d", "point"):
            return True
        if self._name == "vector2d" and other._name in ("vector2d", "vector"):
            return True
        if self._name == "vector3d" and other._name in ("vector3d", "vector"):
            return True

        return False


@dataclass(frozen=True)
class CurveType(Type):
    """A Tier 2 curve primitive type."""
    _name: str

    @property
    def name(self) -> str:
        return self._name

    @property
    def tier(self) -> TypeTier:
        return TypeTier.CURVE


@dataclass(frozen=True)
class CompoundCurveType(Type):
    """A Tier 3 compound curve type (paths, profiles, regions)."""
    _name: str

    @property
    def name(self) -> str:
        return self._name

    @property
    def tier(self) -> TypeTier:
        return TypeTier.COMPOUND_CURVE

    def is_assignable_from(self, other: "Type") -> bool:
        if self == other:
            return True

        if not isinstance(other, CompoundCurveType):
            return False

        # region2d can accept profile2d (auto-promotion when closed)
        if self._name == "region2d" and other._name == "profile2d":
            return True

        return False


@dataclass(frozen=True)
class SurfaceType(Type):
    """A Tier 4 surface type."""
    _name: str

    @property
    def name(self) -> str:
        return self._name

    @property
    def tier(self) -> TypeTier:
        return TypeTier.SURFACE


@dataclass(frozen=True)
class SolidType(Type):
    """A Tier 5 solid type."""
    _name: str = "solid"

    @property
    def name(self) -> str:
        return self._name

    @property
    def tier(self) -> TypeTier:
        return TypeTier.SOLID


@dataclass(frozen=True)
class ListType(Type):
    """A generic list type: list<T>."""
    element_type: Type

    @property
    def name(self) -> str:
        return f"list<{self.element_type.name}>"

    @property
    def tier(self) -> TypeTier:
        return TypeTier.GENERIC

    def is_assignable_from(self, other: "Type") -> bool:
        if not isinstance(other, ListType):
            return False
        return self.element_type.is_assignable_from(other.element_type)


@dataclass(frozen=True)
class DictType(Type):
    """A dictionary type (string keys)."""

    @property
    def name(self) -> str:
        return "dict"

    @property
    def tier(self) -> TypeTier:
        return TypeTier.GENERIC


@dataclass(frozen=True)
class OptionalTypeWrapper(Type):
    """An optional type: T?"""
    inner_type: Type

    @property
    def name(self) -> str:
        return f"{self.inner_type.name}?"

    @property
    def tier(self) -> TypeTier:
        return self.inner_type.tier

    def is_assignable_from(self, other: "Type") -> bool:
        # Optional accepts the inner type or another optional of compatible type
        if isinstance(other, OptionalTypeWrapper):
            return self.inner_type.is_assignable_from(other.inner_type)
        return self.inner_type.is_assignable_from(other)


@dataclass(frozen=True)
class NoneType(Type):
    """The none/null type (only valid for optional types)."""

    @property
    def name(self) -> str:
        return "none"

    @property
    def tier(self) -> TypeTier:
        return TypeTier.PRIMITIVE


@dataclass(frozen=True)
class FunctionType(Type):
    """A function type for lambdas and built-in functions."""
    param_types: Tuple[Type, ...]
    return_type: Type

    @property
    def name(self) -> str:
        params = ", ".join(t.name for t in self.param_types)
        return f"({params}) -> {self.return_type.name}"

    @property
    def tier(self) -> TypeTier:
        return TypeTier.GENERIC


@dataclass(frozen=True)
class UnknownType(Type):
    """A placeholder for type inference or error recovery."""

    @property
    def name(self) -> str:
        return "<unknown>"

    @property
    def tier(self) -> TypeTier:
        return TypeTier.PRIMITIVE

    def is_assignable_from(self, other: "Type") -> bool:
        # Unknown accepts anything (for error recovery)
        return True


@dataclass(frozen=True)
class ErrorType(Type):
    """A type representing a type error (prevents cascading errors)."""

    @property
    def name(self) -> str:
        return "<error>"

    @property
    def tier(self) -> TypeTier:
        return TypeTier.PRIMITIVE

    def is_assignable_from(self, other: "Type") -> bool:
        # Error type accepts anything (prevents cascading errors)
        return True


# =============================================================================
# Built-in Type Instances
# =============================================================================

# Tier 1: Primitives
INT = PrimitiveType("int")
FLOAT = PrimitiveType("float")
BOOL = PrimitiveType("bool")
STRING = PrimitiveType("string")

# Tier 1: Geometric Primitives
POINT = GeometricPrimitiveType("point", dimension=None)
POINT2D = GeometricPrimitiveType("point2d", dimension=2)
POINT3D = GeometricPrimitiveType("point3d", dimension=3)
VECTOR = GeometricPrimitiveType("vector", dimension=None)
VECTOR2D = GeometricPrimitiveType("vector2d", dimension=2)
VECTOR3D = GeometricPrimitiveType("vector3d", dimension=3)
TRANSFORM = GeometricPrimitiveType("transform", dimension=None)

# Tier 2: Curve Primitives
LINE_SEGMENT = CurveType("line_segment")
ARC = CurveType("arc")
CIRCLE = CurveType("circle")
ELLIPSE = CurveType("ellipse")
PARABOLA = CurveType("parabola")
HYPERBOLA = CurveType("hyperbola")
CATMULLROM = CurveType("catmullrom")
NURBS = CurveType("nurbs")
BEZIER = CurveType("bezier")

# Tier 3: Compound Curves
PATH2D = CompoundCurveType("path2d")
PATH3D = CompoundCurveType("path3d")
PROFILE2D = CompoundCurveType("profile2d")
REGION2D = CompoundCurveType("region2d")
LOOP3D = CompoundCurveType("loop3d")

# Tier 4: Surfaces
SURFACE = SurfaceType("surface")
SHELL = SurfaceType("shell")

# Tier 5: Solids
SOLID = SolidType("solid")

# Special types
DICT = DictType()
NONE = NoneType()
UNKNOWN = UnknownType()
ERROR = ErrorType()


# =============================================================================
# Type Registry
# =============================================================================

# Map type names to type instances
BUILTIN_TYPES: Dict[str, Type] = {
    # Tier 1: Scalars
    "int": INT,
    "float": FLOAT,
    "bool": BOOL,
    "string": STRING,

    # Tier 1: Geometric
    "point": POINT,
    "point2d": POINT2D,
    "point3d": POINT3D,
    "vector": VECTOR,
    "vector2d": VECTOR2D,
    "vector3d": VECTOR3D,
    "transform": TRANSFORM,

    # Tier 2: Curves
    "line_segment": LINE_SEGMENT,
    "arc": ARC,
    "circle": CIRCLE,
    "ellipse": ELLIPSE,
    "parabola": PARABOLA,
    "hyperbola": HYPERBOLA,
    "catmullrom": CATMULLROM,
    "nurbs": NURBS,
    "bezier": BEZIER,

    # Tier 3: Compound Curves
    "path2d": PATH2D,
    "path3d": PATH3D,
    "profile2d": PROFILE2D,
    "region2d": REGION2D,
    "loop3d": LOOP3D,

    # Tier 4: Surfaces
    "surface": SURFACE,
    "shell": SHELL,

    # Tier 5: Solids
    "solid": SOLID,

    # Generic
    "dict": DICT,
}


def resolve_type_name(name: str) -> Optional[Type]:
    """Look up a type by name."""
    return BUILTIN_TYPES.get(name)


def make_list_type(element_type: Type) -> ListType:
    """Create a list type with the given element type."""
    return ListType(element_type)


def make_optional_type(inner_type: Type) -> OptionalTypeWrapper:
    """Create an optional type wrapping the given type."""
    return OptionalTypeWrapper(inner_type)


# =============================================================================
# Type Compatibility Helpers
# =============================================================================

def is_numeric(t: Type) -> bool:
    """Check if type is numeric (int or float)."""
    return isinstance(t, PrimitiveType) and t.name in ("int", "float")


def is_geometric_primitive(t: Type) -> bool:
    """Check if type is a geometric primitive (point, vector, transform)."""
    return isinstance(t, GeometricPrimitiveType)


def is_curve(t: Type) -> bool:
    """Check if type is a curve (Tier 2)."""
    return isinstance(t, CurveType)


def is_compound_curve(t: Type) -> bool:
    """Check if type is a compound curve (Tier 3)."""
    return isinstance(t, CompoundCurveType)


def is_surface(t: Type) -> bool:
    """Check if type is a surface (Tier 4)."""
    return isinstance(t, SurfaceType)


def is_solid(t: Type) -> bool:
    """Check if type is a solid (Tier 5)."""
    return isinstance(t, SolidType)


def is_geometry(t: Type) -> bool:
    """Check if type is any geometric type (Tier 1 geometric through Tier 5)."""
    return (is_geometric_primitive(t) or is_curve(t) or
            is_compound_curve(t) or is_surface(t) or is_solid(t))


def common_type(t1: Type, t2: Type) -> Optional[Type]:
    """
    Find the common type that both t1 and t2 can be assigned to.

    Returns None if no common type exists.
    """
    if t1.is_assignable_from(t2):
        return t1
    if t2.is_assignable_from(t1):
        return t2

    # Special case: int and float -> float
    if is_numeric(t1) and is_numeric(t2):
        return FLOAT

    # Special case: point variants -> point
    if isinstance(t1, GeometricPrimitiveType) and isinstance(t2, GeometricPrimitiveType):
        if t1.name.startswith("point") and t2.name.startswith("point"):
            return POINT
        if t1.name.startswith("vector") and t2.name.startswith("vector"):
            return VECTOR

    return None


# =============================================================================
# Curve evaluation method types
# =============================================================================

# Method signatures for parametric curves
CURVE_METHODS: Dict[str, Tuple[List[Tuple[str, Type]], Type]] = {
    # method_name: ([(param_name, param_type), ...], return_type)
    "at": ([("t", FLOAT)], POINT),
    "tangent_at": ([("t", FLOAT)], VECTOR),
    "normal_at": ([("t", FLOAT)], VECTOR),
    "curvature_at": ([("t", FLOAT)], FLOAT),
    "length": ([], FLOAT),
    "split_at": ([("t", FLOAT)], UNKNOWN),  # Returns tuple, handle specially
}
