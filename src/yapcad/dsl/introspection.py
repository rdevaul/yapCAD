"""
DSL Introspection API for agentic engineering tools.

This module provides programmatic access to the yapCAD DSL's type system,
built-in functions, and methods. It's designed to be queried by AI agents
and other tools that need to understand the DSL's capabilities.

Usage:
    from yapcad.dsl.introspection import (
        get_api_reference,
        get_function_info,
        get_type_info,
        list_functions,
        list_types,
        describe_function,
    )

    # Get complete API reference as a dictionary
    api = get_api_reference()

    # Get info about a specific function
    info = get_function_info("box")
    print(info["signature"])  # "box(width: float, depth: float, height: float) -> solid"

    # List all available functions
    for name in list_functions():
        print(name)
"""

from typing import Dict, List, Optional, Any
from dataclasses import dataclass
import json

from .symbols import (
    SymbolTable,
    FunctionSignature,
    CURVE_METHODS,
    SOLID_METHODS,
    REGION2D_METHODS,
    TRANSFORM_METHODS,
)
from .types import (
    Type, TypeTier,
    INT, FLOAT, BOOL, STRING,
    POINT, POINT2D, POINT3D, VECTOR, VECTOR2D, VECTOR3D, TRANSFORM,
    LINE_SEGMENT, ARC, CIRCLE, ELLIPSE, PARABOLA, HYPERBOLA,
    CATMULLROM, NURBS, BEZIER,
    PATH2D, PATH3D, PROFILE2D, REGION2D, LOOP3D,
    SURFACE, SHELL, SOLID,
)


# =============================================================================
# Function Descriptions (for documentation)
# =============================================================================

FUNCTION_DESCRIPTIONS: Dict[str, str] = {
    # Constructors - Points and Vectors
    "point": "Create a 2D or 3D point from coordinates",
    "vector": "Create a 2D or 3D vector from components",

    # Solid Transformations
    "translate": "Move a solid by offset in x, y, z",
    "rotate": "Rotate a solid by angles around x, y, z axes (degrees)",
    "scale": "Scale a solid by factors in x, y, z",

    # Transform Constructors
    "translate_xform": "Create a translation transform from a vector",
    "rotate_xform": "Create a rotation transform around an axis",
    "rotate_2d": "Create a 2D rotation transform by angle (degrees)",
    "scale_xform": "Create a scaling transform from factor vector",
    "scale_uniform": "Create a uniform scaling transform",
    "mirror": "Create a mirror transform across a plane",
    "mirror_2d": "Create a 2D mirror transform across an axis",
    "mirror_y": "Create a mirror transform across the Y axis",
    "identity_transform": "Create an identity (no-op) transform",

    # Curve Constructors
    "line": "Create a line segment between two points",
    "arc": "Create an arc from center, radius, and start/end angles",
    "circle": "Create a circle from center and radius",
    "ellipse": "Create an ellipse from center, major/minor axes, and rotation",
    "bezier": "Create a Bezier curve from control points",
    "catmullrom": "Create a Catmull-Rom spline from control points",
    "nurbs": "Create a NURBS curve from control points, weights, knots, and degree",

    # Path/Region Operations
    "path": "Create a path from a list of curve segments",
    "join": "Join two paths end-to-end",
    "close": "Close a profile to create a region",
    "closeC0": "Close a profile with C0 (position) continuity",
    "closeC1": "Close a profile with C1 (tangent) continuity",
    "rectangle": "Create a rectangular region from width, height, and center",
    "regular_polygon": "Create a regular polygon region from n sides, radius, and center",

    # Surface Operations
    "planar_surface": "Create a planar surface from a boundary loop",
    "cylindrical_surface": "Create a cylindrical surface from axis, radius, and height",
    "loft_surface": "Create a lofted surface between profile curves",
    "shell": "Create a shell from a list of surfaces",

    # Solid Constructors
    "extrude": "Extrude a 2D region along a direction to create a solid",
    "revolve": "Revolve a 2D region around an axis to create a solid",
    "sweep": "Sweep a 2D profile along a 3D path to create a solid",
    "loft": "Create a solid by lofting between 2D profiles",
    "box": "Create a box solid from width, depth, and height",
    "cylinder": "Create a cylinder solid from radius and height",
    "sphere": "Create a sphere solid from radius",
    "cone": "Create a cone/frustum solid from two radii and height",
    "involute_gear": "Create an involute spur gear solid",

    # Boolean Operations
    "union": "Combine two or more solids into one",
    "difference": "Subtract one solid from another",
    "intersection": "Keep only the overlapping volume of two solids",

    # Pattern Operations
    "radial_pattern": "Create a radial/circular pattern of geometry",
    "linear_pattern": "Create a linear pattern of geometry",

    # Query Operations
    "volume": "Calculate the volume of a solid",
    "surface_area": "Calculate the surface area of a solid",
    "area": "Calculate the area of a 2D region",
    "perimeter": "Calculate the perimeter of a 2D region",
    "centroid": "Calculate the centroid of a solid",
    "distance": "Calculate the distance between two points",

    # Transform Application
    "apply": "Apply a transform to a solid",
    "apply_surface": "Apply a transform to a surface",
    "apply_point": "Apply a transform to a point",
    "apply_vector": "Apply a transform to a vector",

    # Math Functions
    "sin": "Sine of angle (radians)",
    "cos": "Cosine of angle (radians)",
    "tan": "Tangent of angle (radians)",
    "asin": "Arc sine, returns radians",
    "acos": "Arc cosine, returns radians",
    "atan": "Arc tangent, returns radians",
    "atan2": "Two-argument arc tangent, returns radians",
    "sqrt": "Square root",
    "abs": "Absolute value",
    "min": "Minimum of two values",
    "max": "Maximum of two values",
    "radians": "Convert degrees to radians",
    "degrees": "Convert radians to degrees",

    # List Operations
    "len": "Get the length of a list",
    "range": "Generate a list of integers [0, end) or [start, end)",
    "print": "Print a value (for debugging)",

    # Utility
    "is_empty": "Check if a solid is empty",
    "empty_solid": "Create an empty solid (useful for accumulating unions)",
    "empty_region": "Create an empty 2D region",
}

FUNCTION_EXAMPLES: Dict[str, str] = {
    "box": 'box_solid: solid = box(10.0, 20.0, 5.0)',
    "cylinder": 'cyl: solid = cylinder(5.0, 10.0)',
    "sphere": 'ball: solid = sphere(5.0)',
    "union": 'combined: solid = union(box1, box2)',
    "difference": 'with_hole: solid = difference(box, cylinder)',
    "translate": 'moved: solid = translate(box, 10.0, 0.0, 0.0)',
    "rotate": 'rotated: solid = rotate(box, 0.0, 0.0, 45.0)',
    "extrude": 'solid: solid = extrude(rectangle(10.0, 5.0, point(0.0, 0.0)), 3.0, vector(0.0, 0.0, 1.0))',
    "involute_gear": 'gear: solid = involute_gear(24, 2.0, 20.0, 10.0)',
    "point": 'p: point = point(1.0, 2.0, 3.0)',
    "vector": 'v: vector = vector(1.0, 0.0, 0.0)',
    "rectangle": 'rect: region2d = rectangle(10.0, 5.0, point(0.0, 0.0))',
    "circle": 'circ: circle = circle(point(0.0, 0.0), 5.0)',
    "line": 'seg: line_segment = line(point(0.0, 0.0), point(10.0, 0.0))',
    "volume": 'vol: float = volume(box)',
    "radians": 'rad: float = radians(90.0)  # Returns pi/2',
    "range": 'nums: list<int> = range(5)  # [0, 1, 2, 3, 4]',
}


# =============================================================================
# Type Information
# =============================================================================

TYPE_INFO: Dict[str, Dict[str, Any]] = {
    # Primitives
    "int": {"tier": 0, "description": "Integer number", "examples": ["42", "-10", "0"]},
    "float": {"tier": 0, "description": "Floating-point number", "examples": ["3.14", "-1.5", "0.0"]},
    "bool": {"tier": 0, "description": "Boolean true/false", "examples": ["true", "false"]},
    "string": {"tier": 0, "description": "Text string", "examples": ['"hello"', '"gear_1"']},

    # Geometric Primitives (Tier 1)
    "point": {"tier": 1, "description": "2D or 3D point in space", "constructor": "point(x, y, z?)"},
    "point2d": {"tier": 1, "description": "2D point", "constructor": "point(x, y)"},
    "point3d": {"tier": 1, "description": "3D point", "constructor": "point(x, y, z)"},
    "vector": {"tier": 1, "description": "2D or 3D direction vector", "constructor": "vector(dx, dy, dz?)"},
    "vector2d": {"tier": 1, "description": "2D vector", "constructor": "vector(dx, dy)"},
    "vector3d": {"tier": 1, "description": "3D vector", "constructor": "vector(dx, dy, dz)"},
    "transform": {"tier": 1, "description": "Affine transformation matrix", "constructor": "translate_xform(), rotate_xform(), scale_xform()"},

    # Curves (Tier 2)
    "line_segment": {"tier": 2, "description": "Straight line between two points", "constructor": "line(start, end)"},
    "arc": {"tier": 2, "description": "Circular arc", "constructor": "arc(center, radius, start_angle, end_angle)"},
    "circle": {"tier": 2, "description": "Full circle", "constructor": "circle(center, radius)"},
    "ellipse": {"tier": 2, "description": "Ellipse", "constructor": "ellipse(center, major, minor, rotation?)"},
    "bezier": {"tier": 2, "description": "Bezier curve", "constructor": "bezier(control_points)"},
    "catmullrom": {"tier": 2, "description": "Catmull-Rom spline", "constructor": "catmullrom(control_points, tension?)"},
    "nurbs": {"tier": 2, "description": "NURBS curve", "constructor": "nurbs(control_points, weights, knots, degree)"},

    # Compound Curves (Tier 3)
    "path2d": {"tier": 3, "description": "2D path of connected curves", "constructor": "path(segments)"},
    "path3d": {"tier": 3, "description": "3D path of connected curves", "constructor": "path(segments)"},
    "profile2d": {"tier": 3, "description": "Open 2D profile (not closed)", "constructor": "path(segments)"},
    "region2d": {"tier": 3, "description": "Closed 2D region", "constructor": "rectangle(), regular_polygon(), close()"},
    "loop3d": {"tier": 3, "description": "Closed 3D loop", "constructor": "close(path3d)"},

    # Surfaces (Tier 4)
    "surface": {"tier": 4, "description": "3D surface", "constructor": "planar_surface(), cylindrical_surface(), loft_surface()"},
    "shell": {"tier": 4, "description": "Collection of surfaces forming a boundary", "constructor": "shell(surfaces)"},

    # Solids (Tier 5)
    "solid": {"tier": 5, "description": "3D solid volume", "constructor": "box(), cylinder(), sphere(), extrude(), revolve()"},
}


# =============================================================================
# Core API Functions
# =============================================================================

def _format_signature(sig: FunctionSignature) -> str:
    """Format a function signature as a string."""
    params = []
    for name, ptype, default in sig.params:
        type_name = ptype.name if hasattr(ptype, 'name') else str(ptype)
        if default is not None:
            params.append(f"{name}: {type_name} = {default}")
        else:
            params.append(f"{name}: {type_name}")

    ret_type = sig.return_type.name if hasattr(sig.return_type, 'name') else str(sig.return_type)
    return f"{sig.name}({', '.join(params)}) -> {ret_type}"


def _sig_to_dict(sig: FunctionSignature) -> Dict[str, Any]:
    """Convert a FunctionSignature to a dictionary."""
    params = []
    for name, ptype, default in sig.params:
        type_name = ptype.name if hasattr(ptype, 'name') else str(ptype)
        param = {"name": name, "type": type_name}
        if default is not None:
            param["default"] = default
        params.append(param)

    ret_type = sig.return_type.name if hasattr(sig.return_type, 'name') else str(sig.return_type)

    return {
        "name": sig.name,
        "parameters": params,
        "return_type": ret_type,
        "is_variadic": sig.is_variadic,
        "is_method": sig.is_method,
        "signature": _format_signature(sig),
    }


def get_api_reference() -> Dict[str, Any]:
    """
    Get the complete API reference as a dictionary.

    Returns a dictionary with:
    - types: All available types with descriptions
    - functions: All built-in functions with signatures and descriptions
    - methods: Type-specific methods organized by receiver type

    This is the primary entry point for agentic tools to understand
    the DSL's capabilities.
    """
    table = SymbolTable()
    builtins = table.get_all_builtins()

    functions = {}
    for name, sig in builtins.items():
        functions[name] = {
            **_sig_to_dict(sig),
            "description": FUNCTION_DESCRIPTIONS.get(name, ""),
            "example": FUNCTION_EXAMPLES.get(name, ""),
        }

    methods = {
        "curve": {name: _sig_to_dict(sig) for name, sig in CURVE_METHODS.items()},
        "solid": {name: _sig_to_dict(sig) for name, sig in SOLID_METHODS.items()},
        "region2d": {name: _sig_to_dict(sig) for name, sig in REGION2D_METHODS.items()},
        "transform": {name: _sig_to_dict(sig) for name, sig in TRANSFORM_METHODS.items()},
    }

    return {
        "version": "1.0.0",
        "types": TYPE_INFO,
        "functions": functions,
        "methods": methods,
    }


def list_functions(category: Optional[str] = None) -> List[str]:
    """
    List all available built-in function names.

    Args:
        category: Optional filter by category (e.g., "solid", "math", "curve")

    Returns:
        Sorted list of function names
    """
    table = SymbolTable()
    names = list(table.get_all_builtins().keys())

    if category:
        # Filter by rough category based on return type or description
        filtered = []
        for name in names:
            sig = table.lookup_builtin(name)
            ret_type = sig.return_type.name if hasattr(sig.return_type, 'name') else ""
            desc = FUNCTION_DESCRIPTIONS.get(name, "").lower()

            if category == "solid" and ret_type == "solid":
                filtered.append(name)
            elif category == "math" and ret_type == "float" and name in ["sin", "cos", "tan", "asin", "acos", "atan", "atan2", "sqrt", "abs", "min", "max", "radians", "degrees"]:
                filtered.append(name)
            elif category == "curve" and ret_type in ["line_segment", "arc", "circle", "ellipse", "bezier", "catmullrom", "nurbs"]:
                filtered.append(name)
            elif category == "transform" and (ret_type == "transform" or "transform" in desc):
                filtered.append(name)
            elif category == "query" and name in ["volume", "surface_area", "area", "perimeter", "centroid", "distance", "is_empty"]:
                filtered.append(name)

        names = filtered

    return sorted(names)


def list_types(tier: Optional[int] = None) -> List[str]:
    """
    List all available types.

    Args:
        tier: Optional filter by tier (0-5)

    Returns:
        List of type names
    """
    types = []
    for name, info in TYPE_INFO.items():
        if tier is None or info.get("tier") == tier:
            types.append(name)
    return sorted(types)


def get_function_info(name: str) -> Optional[Dict[str, Any]]:
    """
    Get detailed information about a specific function.

    Args:
        name: The function name

    Returns:
        Dictionary with signature, description, example, etc.
        Returns None if function not found.
    """
    table = SymbolTable()
    sig = table.lookup_builtin(name)

    if sig is None:
        return None

    return {
        **_sig_to_dict(sig),
        "description": FUNCTION_DESCRIPTIONS.get(name, ""),
        "example": FUNCTION_EXAMPLES.get(name, ""),
    }


def get_type_info(name: str) -> Optional[Dict[str, Any]]:
    """
    Get information about a type.

    Args:
        name: The type name

    Returns:
        Dictionary with tier, description, constructor info.
        Returns None if type not found.
    """
    return TYPE_INFO.get(name.lower())


def get_methods_for_type(type_name: str) -> Dict[str, Dict[str, Any]]:
    """
    Get all methods available on a given type.

    Args:
        type_name: The type name (e.g., "solid", "curve", "region2d")

    Returns:
        Dictionary mapping method names to their signatures
    """
    type_name = type_name.lower()

    if type_name in ["solid"]:
        return {name: _sig_to_dict(sig) for name, sig in SOLID_METHODS.items()}
    elif type_name in ["curve", "line_segment", "arc", "circle", "ellipse", "bezier", "catmullrom", "nurbs"]:
        return {name: _sig_to_dict(sig) for name, sig in CURVE_METHODS.items()}
    elif type_name in ["region2d"]:
        return {name: _sig_to_dict(sig) for name, sig in REGION2D_METHODS.items()}
    elif type_name in ["transform"]:
        return {name: _sig_to_dict(sig) for name, sig in TRANSFORM_METHODS.items()}

    return {}


def describe_function(name: str) -> str:
    """
    Get a human-readable description of a function.

    Args:
        name: The function name

    Returns:
        Formatted description string
    """
    info = get_function_info(name)
    if info is None:
        return f"Unknown function: {name}"

    lines = [
        f"Function: {name}",
        f"Signature: {info['signature']}",
        f"Description: {info['description'] or 'No description available'}",
    ]

    if info.get("example"):
        lines.append(f"Example: {info['example']}")

    if info.get("is_variadic"):
        lines.append("Note: This function accepts variable arguments or a list")

    return "\n".join(lines)


def get_api_as_json() -> str:
    """
    Get the complete API reference as a JSON string.

    Useful for tools that prefer to parse JSON directly.
    """
    return json.dumps(get_api_reference(), indent=2)


# =============================================================================
# Quick Reference for Common Tasks
# =============================================================================

COMMON_PATTERNS = {
    "create_box": """
# Create a simple box
b: solid = box(width, depth, height)
""",

    "create_cylinder": """
# Create a cylinder
c: solid = cylinder(radius, height)
""",

    "boolean_subtraction": """
# Create a box with a hole through it
block: solid = box(20.0, 20.0, 10.0)
hole: solid = cylinder(5.0, 15.0)
hole_centered: solid = translate(hole, 10.0, 10.0, -2.5)
result: solid = difference(block, hole_centered)
emit result
""",

    "gear_creation": """
# Create an involute spur gear
# Parameters: teeth, module (mm), pressure angle (deg), face width
gear: solid = involute_gear(24, 2.0, 20.0, 10.0)
emit gear
""",

    "extrusion": """
# Extrude a 2D shape into a 3D solid
profile: region2d = rectangle(10.0, 5.0, point(0.0, 0.0))
solid: solid = extrude(profile, 3.0, vector(0.0, 0.0, 1.0))
emit solid
""",

    "transform_chain": """
# Apply multiple transformations to a solid
base: solid = box(10.0, 10.0, 10.0)
rotated: solid = rotate(base, 0.0, 0.0, 45.0)
moved: solid = translate(rotated, 20.0, 0.0, 0.0)
emit moved
""",

    "parametric_design": """
# Parametric design with variables
command MAKE_BRACKET(width: float, height: float, thickness: float, hole_radius: float) -> solid:
    # Create main plate
    plate: solid = box(width, height, thickness)

    # Create hole
    hole: solid = cylinder(hole_radius, thickness + 1.0)
    hole_pos: solid = translate(hole, width/2.0, height/2.0, -0.5)

    # Subtract hole from plate
    result: solid = difference(plate, hole_pos)
    emit result
""",
}


def get_common_pattern(name: str) -> Optional[str]:
    """
    Get a common DSL pattern/example.

    Args:
        name: Pattern name (e.g., "boolean_subtraction", "gear_creation")

    Returns:
        DSL code example or None if pattern not found
    """
    return COMMON_PATTERNS.get(name)


def list_common_patterns() -> List[str]:
    """List all available common pattern names."""
    return list(COMMON_PATTERNS.keys())
