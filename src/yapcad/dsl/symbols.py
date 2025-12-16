"""
Symbol table management for the yapCAD DSL type checker.

Provides scoped symbol tables for tracking variable definitions,
function signatures, and type information during type checking.
"""

from dataclasses import dataclass, field
from typing import Optional, Dict, List, Callable, Any, Tuple
from enum import Enum, auto

from .types import (
    Type, FunctionType, ListType, ERROR, UNKNOWN,
    INT, FLOAT, BOOL, STRING,
    POINT, POINT2D, POINT3D, VECTOR, VECTOR2D, VECTOR3D, TRANSFORM,
    LINE_SEGMENT, ARC, CIRCLE, ELLIPSE, PARABOLA, HYPERBOLA,
    CATMULLROM, NURBS, BEZIER,
    PATH2D, PATH3D, PROFILE2D, REGION2D, LOOP3D,
    SURFACE, SHELL, SOLID, DICT,
    make_list_type, make_optional_type,
)
from .tokens import SourceSpan


class SymbolKind(Enum):
    """The kind of symbol being tracked."""
    VARIABLE = auto()
    PARAMETER = auto()
    FUNCTION = auto()
    COMMAND = auto()
    MODULE = auto()


@dataclass
class Symbol:
    """A symbol in the symbol table."""
    name: str
    kind: SymbolKind
    type: Type
    span: Optional[SourceSpan] = None  # Where it was defined
    is_mutable: bool = False
    has_python_block: bool = False  # Commands with Python blocks


@dataclass
class FunctionSignature:
    """Type signature for a function or built-in."""
    name: str
    params: List[Tuple[str, Type, Optional[Any]]]  # (name, type, default)
    return_type: Type
    is_method: bool = False  # True for methods like curve.at()
    is_variadic: bool = False  # True for functions that accept list<T> as well


@dataclass
class Scope:
    """A single scope in the scope stack."""
    symbols: Dict[str, Symbol] = field(default_factory=dict)
    parent: Optional["Scope"] = None
    name: str = ""  # For debugging: "command MAKE_GEAR", "for loop", etc.

    def define(self, symbol: Symbol) -> None:
        """Define a symbol in this scope."""
        self.symbols[symbol.name] = symbol

    def lookup_local(self, name: str) -> Optional[Symbol]:
        """Look up a symbol in this scope only."""
        return self.symbols.get(name)

    def lookup(self, name: str) -> Optional[Symbol]:
        """Look up a symbol in this scope or any parent scope."""
        symbol = self.symbols.get(name)
        if symbol is not None:
            return symbol
        if self.parent is not None:
            return self.parent.lookup(name)
        return None


class SymbolTable:
    """
    Manages scopes and symbol definitions for type checking.

    Provides:
    - Nested scope management (push/pop)
    - Symbol definition and lookup
    - Built-in function registry
    """

    def __init__(self):
        # Global scope contains built-in functions and types
        self._global_scope = Scope(name="global")
        self._current_scope = self._global_scope

        # Built-in function signatures
        self._builtins: Dict[str, FunctionSignature] = {}

        # Initialize built-ins
        self._init_builtins()

    def _init_builtins(self) -> None:
        """Initialize built-in function signatures."""
        # Tier 1: Point and Vector constructors
        self._register_builtin("point", [
            ("x", FLOAT, None),
            ("y", FLOAT, None),
            ("z", FLOAT, None),  # Optional - presence determines 2D/3D
        ], POINT)

        # Convenience constructor for 2D points (z=0)
        self._register_builtin("point2d", [
            ("x", FLOAT, None),
            ("y", FLOAT, None),
        ], POINT2D)

        self._register_builtin("vector", [
            ("dx", FLOAT, None),
            ("dy", FLOAT, None),
            ("dz", FLOAT, None),
        ], VECTOR)

        # Convenience constructor for 2D vectors (dz=0)
        self._register_builtin("vector2d", [
            ("dx", FLOAT, None),
            ("dy", FLOAT, None),
        ], VECTOR2D)

        # Solid transformation functions (translate, rotate, scale solids directly)
        self._register_builtin("translate", [
            ("s", SOLID, None),
            ("x", FLOAT, None),
            ("y", FLOAT, None),
            ("z", FLOAT, None),
        ], SOLID)
        self._register_builtin("rotate", [
            ("s", SOLID, None),
            ("x", FLOAT, None),
            ("y", FLOAT, None),
            ("z", FLOAT, None),
        ], SOLID)
        self._register_builtin("scale", [
            ("s", SOLID, None),
            ("x", FLOAT, None),
            ("y", FLOAT, None),
            ("z", FLOAT, None),
        ], SOLID)

        # Transform constructors (for advanced use - create transform matrices)
        self._register_builtin("translate_xform", [("v", VECTOR, None)], TRANSFORM)
        self._register_builtin("rotate_xform", [
            ("axis", VECTOR3D, None),
            ("angle", FLOAT, None),
        ], TRANSFORM)
        self._register_builtin("rotate_2d", [("angle", FLOAT, None)], TRANSFORM)
        self._register_builtin("scale_xform", [("factors", VECTOR, None)], TRANSFORM)
        self._register_builtin("scale_uniform", [("factor", FLOAT, None)], TRANSFORM)
        self._register_builtin("mirror", [("plane_normal", VECTOR3D, None)], TRANSFORM)
        self._register_builtin("mirror_2d", [("axis", VECTOR2D, None)], TRANSFORM)
        self._register_builtin("mirror_y", [], TRANSFORM)  # Convenience
        self._register_builtin("identity_transform", [], TRANSFORM)

        # Tier 2: Curve constructors
        self._register_builtin("line", [
            ("start", POINT, None),
            ("end", POINT, None),
        ], LINE_SEGMENT)

        self._register_builtin("arc", [
            ("center", POINT, None),
            ("radius", FLOAT, None),
            ("start_angle", FLOAT, None),
            ("end_angle", FLOAT, None),
        ], ARC)

        self._register_builtin("circle", [
            ("center", POINT, None),
            ("radius", FLOAT, None),
        ], CIRCLE)

        self._register_builtin("ellipse", [
            ("center", POINT, None),
            ("major", FLOAT, None),
            ("minor", FLOAT, None),
            ("rotation", FLOAT, 0.0),
        ], ELLIPSE)

        self._register_builtin("bezier", [
            ("control_points", make_list_type(POINT), None),
        ], BEZIER)

        self._register_builtin("catmullrom", [
            ("control_points", make_list_type(POINT), None),
            ("tension", FLOAT, 0.5),
        ], CATMULLROM)

        self._register_builtin("nurbs", [
            ("control_points", make_list_type(POINT), None),
            ("weights", make_list_type(FLOAT), None),
            ("knots", make_list_type(FLOAT), None),
            ("degree", INT, None),
        ], NURBS)

        # Tier 3: Path/Region operations
        self._register_builtin("path", [
            ("segments", make_list_type(UNKNOWN), None),  # list of curves
        ], PATH2D)  # Returns path2d or path3d based on input

        self._register_builtin("join", [
            ("p1", PATH2D, None),
            ("p2", PATH2D, None),
        ], PATH2D)

        # Path3D constructors for sweep operations
        self._register_builtin("make_path3d", [
            ("segments", make_list_type(PATH3D), None),
        ], PATH3D)

        self._register_builtin("path3d_line", [
            ("start", POINT3D, None),
            ("end", POINT3D, None),
        ], PATH3D)

        self._register_builtin("path3d_arc", [
            ("center", POINT3D, None),
            ("start", POINT3D, None),
            ("end", POINT3D, None),
            ("normal", VECTOR3D, None),
        ], PATH3D)

        self._register_builtin("close", [("p", PROFILE2D, None)], REGION2D)
        self._register_builtin("closeC0", [("p", PROFILE2D, None)], REGION2D)
        self._register_builtin("closeC1", [("p", PROFILE2D, None)], REGION2D)

        self._register_builtin("rectangle", [
            ("width", FLOAT, None),
            ("height", FLOAT, None),
            ("center", POINT2D, None),
        ], REGION2D)

        self._register_builtin("regular_polygon", [
            ("n", INT, None),
            ("radius", FLOAT, None),
            ("center", POINT2D, None),
        ], REGION2D)

        # Tier 4: Surface operations
        self._register_builtin("planar_surface", [
            ("boundary", LOOP3D, None),
        ], SURFACE)

        self._register_builtin("cylindrical_surface", [
            ("axis", VECTOR3D, None),
            ("radius", FLOAT, None),
            ("height", FLOAT, None),
        ], SURFACE)

        self._register_builtin("loft_surface", [
            ("profiles", make_list_type(PATH3D), None),
        ], SURFACE)

        self._register_builtin("shell", [
            ("surfaces", make_list_type(SURFACE), None),
        ], SHELL)

        # Tier 5: Solid operations
        self._register_builtin("extrude", [
            ("profile", REGION2D, None),
            ("height", FLOAT, None),
            ("direction", VECTOR3D, None),
        ], SOLID)

        self._register_builtin("revolve", [
            ("profile", REGION2D, None),
            ("axis", VECTOR3D, None),
            ("angle", FLOAT, None),
        ], SOLID)

        self._register_builtin("sweep", [
            ("profile", REGION2D, None),
            ("spine", PATH3D, None),
        ], SOLID)

        self._register_builtin("sweep_hollow", [
            ("outer_profile", REGION2D, None),
            ("inner_profile", REGION2D, None),
            ("spine", PATH3D, None),
        ], SOLID)

        self._register_builtin("sweep_adaptive", [
            ("profile", REGION2D, None),
            ("spine", PATH3D, None),
            ("threshold", FLOAT, None),  # Angle in degrees
        ], SOLID)

        self._register_builtin("sweep_adaptive_hollow", [
            ("outer_profile", REGION2D, None),
            ("inner_profiles", REGION2D, None),  # Single region2d or list
            ("spine", PATH3D, None),
            ("threshold", FLOAT, None),
        ], SOLID)

        self._register_builtin("loft", [
            ("profiles", make_list_type(REGION2D), None),
        ], SOLID)

        self._register_builtin("box", [
            ("width", FLOAT, None),
            ("depth", FLOAT, None),
            ("height", FLOAT, None),
        ], SOLID)

        self._register_builtin("cylinder", [
            ("radius", FLOAT, None),
            ("height", FLOAT, None),
        ], SOLID)

        self._register_builtin("sphere", [
            ("radius", FLOAT, None),
        ], SOLID)

        self._register_builtin("cone", [
            ("radius1", FLOAT, None),
            ("radius2", FLOAT, None),
            ("height", FLOAT, None),
        ], SOLID)

        # Involute gear - creates a proper gear profile
        self._register_builtin("involute_gear", [
            ("teeth", INT, None),
            ("module_mm", FLOAT, None),
            ("pressure_angle", FLOAT, None),
            ("face_width", FLOAT, None),
        ], SOLID)

        # Boolean operations (variadic)
        self._register_builtin("union", [
            ("a", SOLID, None),
            ("b", SOLID, None),
        ], SOLID, is_variadic=True)

        self._register_builtin("difference", [
            ("a", SOLID, None),
            ("b", SOLID, None),
        ], SOLID, is_variadic=True)

        self._register_builtin("intersection", [
            ("a", SOLID, None),
            ("b", SOLID, None),
        ], SOLID, is_variadic=True)

        # Pattern operations
        self._register_builtin("radial_pattern", [
            ("shape", UNKNOWN, None),  # Any geometry
            ("count", INT, None),
            ("axis", VECTOR3D, None),
            ("center", POINT, None),
        ], UNKNOWN)  # Returns same type as input

        self._register_builtin("linear_pattern", [
            ("shape", UNKNOWN, None),
            ("count", INT, None),
            ("spacing", VECTOR, None),
        ], UNKNOWN)

        # Query operations
        self._register_builtin("volume", [("s", SOLID, None)], FLOAT)
        self._register_builtin("surface_area", [("s", SOLID, None)], FLOAT)
        self._register_builtin("area", [("r", REGION2D, None)], FLOAT)
        self._register_builtin("perimeter", [("r", REGION2D, None)], FLOAT)
        self._register_builtin("centroid", [("s", SOLID, None)], POINT3D)
        self._register_builtin("distance", [
            ("a", POINT, None),
            ("b", POINT, None),
            ("tolerance", FLOAT, None),
        ], FLOAT)

        # Transform application - type-specific versions
        self._register_builtin("apply", [
            ("t", TRANSFORM, None),
            ("shape", SOLID, None),
        ], SOLID)
        self._register_builtin("apply_surface", [
            ("t", TRANSFORM, None),
            ("shape", SURFACE, None),
        ], SURFACE)
        self._register_builtin("apply_point", [
            ("t", TRANSFORM, None),
            ("p", POINT, None),
        ], POINT)
        self._register_builtin("apply_vector", [
            ("t", TRANSFORM, None),
            ("v", VECTOR3D, None),
        ], VECTOR3D)

        # Math functions
        self._register_builtin("sin", [("x", FLOAT, None)], FLOAT)
        self._register_builtin("cos", [("x", FLOAT, None)], FLOAT)
        self._register_builtin("tan", [("x", FLOAT, None)], FLOAT)
        self._register_builtin("asin", [("x", FLOAT, None)], FLOAT)
        self._register_builtin("acos", [("x", FLOAT, None)], FLOAT)
        self._register_builtin("atan", [("x", FLOAT, None)], FLOAT)
        self._register_builtin("atan2", [("y", FLOAT, None), ("x", FLOAT, None)], FLOAT)
        self._register_builtin("sqrt", [("x", FLOAT, None)], FLOAT)
        self._register_builtin("abs", [("x", FLOAT, None)], FLOAT)
        self._register_builtin("min", [("a", FLOAT, None), ("b", FLOAT, None)], FLOAT)
        self._register_builtin("max", [("a", FLOAT, None), ("b", FLOAT, None)], FLOAT)
        self._register_builtin("radians", [("degrees", FLOAT, None)], FLOAT)
        self._register_builtin("degrees", [("radians", FLOAT, None)], FLOAT)
        self._register_builtin("floor", [("x", FLOAT, None)], INT)
        self._register_builtin("ceil", [("x", FLOAT, None)], INT)
        self._register_builtin("round", [("x", FLOAT, None)], INT)
        self._register_builtin("pow", [("base", FLOAT, None), ("exp", FLOAT, None)], FLOAT)
        self._register_builtin("pi", [], FLOAT)  # Constant: pi
        self._register_builtin("tau", [], FLOAT)  # Constant: tau (2*pi)

        # List operations
        self._register_builtin("len", [("list", make_list_type(UNKNOWN), None)], INT)
        self._register_builtin("range", [("end", INT, None)], make_list_type(INT), is_variadic=True)
        self._register_builtin("concat", [
            ("list1", make_list_type(UNKNOWN), None),
            ("list2", make_list_type(UNKNOWN), None),
        ], make_list_type(UNKNOWN))  # Concatenate two lists
        self._register_builtin("reverse", [
            ("list", make_list_type(UNKNOWN), None),
        ], make_list_type(UNKNOWN))  # Reverse a list
        self._register_builtin("flatten", [
            ("list", make_list_type(make_list_type(UNKNOWN)), None),
        ], make_list_type(UNKNOWN))  # Flatten nested list
        self._register_builtin("print", [("value", UNKNOWN, None)], BOOL, is_variadic=True)

        # Utility
        self._register_builtin("is_empty", [("s", SOLID, None)], BOOL)
        self._register_builtin("empty_solid", [], SOLID)
        self._register_builtin("empty_region", [], REGION2D)

    def _register_builtin(
        self,
        name: str,
        params: List[Tuple[str, Type, Optional[Any]]],
        return_type: Type,
        is_variadic: bool = False
    ) -> None:
        """Register a built-in function signature."""
        self._builtins[name] = FunctionSignature(
            name=name,
            params=params,
            return_type=return_type,
            is_variadic=is_variadic
        )

    def push_scope(self, name: str = "") -> None:
        """Push a new scope onto the stack."""
        new_scope = Scope(parent=self._current_scope, name=name)
        self._current_scope = new_scope

    def pop_scope(self) -> None:
        """Pop the current scope."""
        if self._current_scope.parent is not None:
            self._current_scope = self._current_scope.parent

    def define(self, symbol: Symbol) -> bool:
        """
        Define a symbol in the current scope.

        Returns True if successful, False if already defined in current scope.
        """
        if self._current_scope.lookup_local(symbol.name) is not None:
            return False
        self._current_scope.define(symbol)
        return True

    def lookup(self, name: str) -> Optional[Symbol]:
        """Look up a symbol in the current scope chain."""
        return self._current_scope.lookup(name)

    def lookup_builtin(self, name: str) -> Optional[FunctionSignature]:
        """Look up a built-in function signature."""
        return self._builtins.get(name)

    def is_builtin(self, name: str) -> bool:
        """Check if a name is a built-in function."""
        return name in self._builtins

    def current_scope_name(self) -> str:
        """Get the name of the current scope (for debugging)."""
        return self._current_scope.name

    def get_all_builtins(self) -> Dict[str, FunctionSignature]:
        """Get all registered built-in functions."""
        return dict(self._builtins)


# =============================================================================
# Method Signatures for Object Types
# =============================================================================

# Methods available on curve types (Tier 2)
CURVE_METHODS: Dict[str, FunctionSignature] = {
    "at": FunctionSignature(
        name="at",
        params=[("t", FLOAT, None)],
        return_type=POINT,
        is_method=True
    ),
    "tangent_at": FunctionSignature(
        name="tangent_at",
        params=[("t", FLOAT, None)],
        return_type=VECTOR,
        is_method=True
    ),
    "normal_at": FunctionSignature(
        name="normal_at",
        params=[("t", FLOAT, None)],
        return_type=VECTOR,
        is_method=True
    ),
    "curvature_at": FunctionSignature(
        name="curvature_at",
        params=[("t", FLOAT, None)],
        return_type=FLOAT,
        is_method=True
    ),
    "length": FunctionSignature(
        name="length",
        params=[],
        return_type=FLOAT,
        is_method=True
    ),
}

# Methods available on solid types (Tier 5)
SOLID_METHODS: Dict[str, FunctionSignature] = {
    "union": FunctionSignature(
        name="union",
        params=[("other", SOLID, None)],
        return_type=SOLID,
        is_method=True
    ),
    "difference": FunctionSignature(
        name="difference",
        params=[("other", SOLID, None)],
        return_type=SOLID,
        is_method=True
    ),
    "intersection": FunctionSignature(
        name="intersection",
        params=[("other", SOLID, None)],
        return_type=SOLID,
        is_method=True
    ),
    "translate": FunctionSignature(
        name="translate",
        params=[("v", VECTOR, None)],
        return_type=SOLID,
        is_method=True
    ),
    "rotate": FunctionSignature(
        name="rotate",
        params=[("axis", VECTOR3D, None), ("angle", FLOAT, None)],
        return_type=SOLID,
        is_method=True
    ),
    "scale": FunctionSignature(
        name="scale",
        params=[("factors", VECTOR, None)],
        return_type=SOLID,
        is_method=True
    ),
    "apply": FunctionSignature(
        name="apply",
        params=[("t", TRANSFORM, None)],
        return_type=SOLID,
        is_method=True
    ),
}

# Methods available on region2d types (Tier 3)
REGION2D_METHODS: Dict[str, FunctionSignature] = {
    "union": FunctionSignature(
        name="union",
        params=[("other", REGION2D, None)],
        return_type=REGION2D,
        is_method=True
    ),
    "difference": FunctionSignature(
        name="difference",
        params=[("other", REGION2D, None)],
        return_type=REGION2D,
        is_method=True
    ),
    "intersection": FunctionSignature(
        name="intersection",
        params=[("other", REGION2D, None)],
        return_type=REGION2D,
        is_method=True
    ),
}

# Methods available on transform types
TRANSFORM_METHODS: Dict[str, FunctionSignature] = {
    "compose": FunctionSignature(
        name="compose",
        params=[("other", TRANSFORM, None)],
        return_type=TRANSFORM,
        is_method=True
    ),
    "inverse": FunctionSignature(
        name="inverse",
        params=[],
        return_type=TRANSFORM,
        is_method=True
    ),
    "translation": FunctionSignature(
        name="translation",
        params=[],
        return_type=VECTOR3D,
        is_method=True
    ),
    "is_rigid": FunctionSignature(
        name="is_rigid",
        params=[],
        return_type=BOOL,
        is_method=True
    ),
}


def get_method_signature(obj_type: Type, method_name: str) -> Optional[FunctionSignature]:
    """Get the method signature for a type, if it exists."""
    from .types import is_curve, is_solid, is_compound_curve

    if is_curve(obj_type):
        return CURVE_METHODS.get(method_name)

    if is_solid(obj_type):
        return SOLID_METHODS.get(method_name)

    if is_compound_curve(obj_type) and obj_type.name == "region2d":
        return REGION2D_METHODS.get(method_name)

    if obj_type.name == "transform":
        return TRANSFORM_METHODS.get(method_name)

    return None
