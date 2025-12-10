"""
Built-in function registry for DSL interpreter.

Maps DSL function names to yapCAD implementations.

TODO: Add missing 2D curve primitives:
  - ellipse(center, major, minor, rotation?) -> ellipse
  - parabola(vertex, focus) -> parabola
  - hyperbola(center, a, b) -> hyperbola
  Types are defined in types.py but constructors not yet implemented.
  See yapCAD core: yapcad.geom.Ellipse, yapcad.geom.conic_arc

TODO: Move involute_gear to yapcad.stdlib.gears package
  High-level parametric constructions (gears, fasteners) should be
  DSL modules in a standard library, not hardcoded builtins. This allows:
  - Multiple implementations (figgear vs other algorithms)
  - User customization and extension
  - Proper namespace management (use gears.involute)

TODO: Add fastener primitives via yapcad.stdlib.fasteners package
  - hex_bolt(standard, size, length) using existing fastener catalog
  - hex_nut(standard, size)
  - washer(standard, size)
  - threaded_hole(standard, size, depth, tapped?)
  See: yapcad.contrib.fasteners for existing implementation
"""

from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Optional, Tuple
import math

from .values import (
    Value, int_val, float_val, bool_val, string_val, list_val,
    point_val, vector_val, transform_val, solid_val, region2d_val,
    surface_val, shell_val, wrap_value, unwrap_value, unwrap_values,
)
from ..types import (
    Type, ListType,
    INT, FLOAT, BOOL, STRING, UNKNOWN,
    POINT, POINT2D, POINT3D, VECTOR, VECTOR2D, VECTOR3D, TRANSFORM,
    SOLID, REGION2D, SURFACE, SHELL,
    LINE_SEGMENT, ARC, CIRCLE, ELLIPSE, BEZIER, NURBS,
    PATH2D, PATH3D, PROFILE2D,
)
from ..symbols import FunctionSignature


def _make_sig(name: str, param_types: List[Type], return_type: Type,
              is_variadic: bool = False) -> FunctionSignature:
    """
    Helper to create FunctionSignature from a simple list of parameter types.

    Generates synthetic parameter names (arg0, arg1, ...) for builtins.
    """
    params = [(f"arg{i}", t, None) for i, t in enumerate(param_types)]
    return FunctionSignature(
        name=name,
        params=params,
        return_type=return_type,
        is_variadic=is_variadic,
    )


@dataclass
class BuiltinFunction:
    """
    A built-in function with its implementation and type signature.
    """
    name: str
    signature: FunctionSignature
    implementation: Callable[..., Value]
    doc: str = ""


class BuiltinRegistry:
    """
    Registry of all built-in functions.

    Functions are registered by name and can be looked up for execution.
    """

    def __init__(self):
        self._functions: Dict[str, BuiltinFunction] = {}
        self._methods: Dict[Tuple[str, str], BuiltinFunction] = {}  # (type_name, method_name)
        self._register_all()

    def get_function(self, name: str) -> Optional[BuiltinFunction]:
        """Look up a function by name."""
        return self._functions.get(name)

    def get_method(self, type_name: str, method_name: str) -> Optional[BuiltinFunction]:
        """Look up a method by type and method name."""
        return self._methods.get((type_name, method_name))

    def register(self, func: BuiltinFunction) -> None:
        """Register a function."""
        self._functions[func.name] = func

    def register_method(self, type_name: str, func: BuiltinFunction) -> None:
        """Register a method for a specific type."""
        self._methods[(type_name, func.name)] = func

    def _register_all(self) -> None:
        """Register all built-in functions."""
        self._register_math_functions()
        self._register_point_vector_functions()
        self._register_transform_functions()
        self._register_curve_functions()
        self._register_region_functions()
        self._register_solid_functions()
        self._register_boolean_functions()
        self._register_query_functions()
        self._register_utility_functions()

    # --- Math Functions ---

    def _register_math_functions(self) -> None:
        """Register mathematical functions."""

        def _sin(x: Value) -> Value:
            return float_val(math.sin(x.data))

        def _cos(x: Value) -> Value:
            return float_val(math.cos(x.data))

        def _tan(x: Value) -> Value:
            return float_val(math.tan(x.data))

        def _asin(x: Value) -> Value:
            return float_val(math.asin(x.data))

        def _acos(x: Value) -> Value:
            return float_val(math.acos(x.data))

        def _atan(x: Value) -> Value:
            return float_val(math.atan(x.data))

        def _atan2(y: Value, x: Value) -> Value:
            return float_val(math.atan2(y.data, x.data))

        def _sqrt(x: Value) -> Value:
            return float_val(math.sqrt(x.data))

        def _abs(x: Value) -> Value:
            if x.type == INT:
                return int_val(abs(x.data))
            return float_val(abs(x.data))

        def _floor(x: Value) -> Value:
            return int_val(math.floor(x.data))

        def _ceil(x: Value) -> Value:
            return int_val(math.ceil(x.data))

        def _round(x: Value) -> Value:
            return int_val(round(x.data))

        def _min(*args: Value) -> Value:
            values = [a.data for a in args]
            result = min(values)
            if all(a.type == INT for a in args):
                return int_val(result)
            return float_val(result)

        def _max(*args: Value) -> Value:
            values = [a.data for a in args]
            result = max(values)
            if all(a.type == INT for a in args):
                return int_val(result)
            return float_val(result)

        def _pow(base: Value, exp: Value) -> Value:
            return float_val(math.pow(base.data, exp.data))

        def _radians(deg: Value) -> Value:
            return float_val(math.radians(deg.data))

        def _degrees(rad: Value) -> Value:
            return float_val(math.degrees(rad.data))

        def _pi() -> Value:
            return float_val(math.pi)

        def _tau() -> Value:
            return float_val(math.tau)

        # Register math functions
        math_funcs = [
            ("sin", [FLOAT], FLOAT, _sin),
            ("cos", [FLOAT], FLOAT, _cos),
            ("tan", [FLOAT], FLOAT, _tan),
            ("asin", [FLOAT], FLOAT, _asin),
            ("acos", [FLOAT], FLOAT, _acos),
            ("atan", [FLOAT], FLOAT, _atan),
            ("atan2", [FLOAT, FLOAT], FLOAT, _atan2),
            ("sqrt", [FLOAT], FLOAT, _sqrt),
            ("abs", [FLOAT], FLOAT, _abs),
            ("floor", [FLOAT], INT, _floor),
            ("ceil", [FLOAT], INT, _ceil),
            ("round", [FLOAT], INT, _round),
            ("pow", [FLOAT, FLOAT], FLOAT, _pow),
            ("radians", [FLOAT], FLOAT, _radians),
            ("degrees", [FLOAT], FLOAT, _degrees),
        ]

        for name, param_types, return_type, impl in math_funcs:
            sig = _make_sig(name, param_types, return_type)
            self.register(BuiltinFunction(name, sig, impl))

        # Variadic min/max
        self.register(BuiltinFunction(
            "min",
            _make_sig("min", [FLOAT], FLOAT, is_variadic=True),
            _min,
        ))
        self.register(BuiltinFunction(
            "max",
            _make_sig("max", [FLOAT], FLOAT, is_variadic=True),
            _max,
        ))

        # Constants (zero-arg functions)
        self.register(BuiltinFunction(
            "pi",
            _make_sig("pi", [], FLOAT),
            _pi,
        ))
        self.register(BuiltinFunction(
            "tau",
            _make_sig("tau", [], FLOAT),
            _tau,
        ))

    # --- Point/Vector Functions ---

    def _register_point_vector_functions(self) -> None:
        """Register point and vector construction functions."""

        def _point(*args: Value) -> Value:
            """Create a point (2D or 3D based on argument count)."""
            coords = [a.data for a in args]
            if len(coords) == 2:
                # 2D point - yapCAD uses [x, y, 0, 1] homogeneous coords
                from yapcad.geom import point
                return point_val(point(coords[0], coords[1]), is_2d=True)
            else:
                # 3D point
                from yapcad.geom import point
                return point_val(point(coords[0], coords[1], coords[2]), is_2d=False)

        def _vector(*args: Value) -> Value:
            """Create a vector (2D or 3D based on argument count)."""
            coords = [a.data for a in args]
            if len(coords) == 2:
                from yapcad.geom import vect
                return vector_val(vect(coords[0], coords[1]), is_2d=True)
            else:
                from yapcad.geom import vect
                return vector_val(vect(coords[0], coords[1], coords[2]), is_2d=False)

        def _point2d(x: Value, y: Value) -> Value:
            """Create a 2D point."""
            from yapcad.geom import point
            return point_val(point(x.data, y.data), is_2d=True)

        def _vector2d(dx: Value, dy: Value) -> Value:
            """Create a 2D vector."""
            from yapcad.geom import vect
            return vector_val(vect(dx.data, dy.data), is_2d=True)

        # Point constructors (variadic to allow 2D or 3D)
        self.register(BuiltinFunction(
            "point",
            _make_sig("point", [FLOAT, FLOAT, FLOAT], POINT3D, is_variadic=True),
            _point,
        ))

        # 2D point constructor (convenience)
        self.register(BuiltinFunction(
            "point2d",
            _make_sig("point2d", [FLOAT, FLOAT], POINT2D),
            _point2d,
        ))

        # Vector constructors (variadic to allow 2D or 3D)
        self.register(BuiltinFunction(
            "vector",
            _make_sig("vector", [FLOAT, FLOAT, FLOAT], VECTOR3D, is_variadic=True),
            _vector,
        ))

        # 2D vector constructor (convenience)
        self.register(BuiltinFunction(
            "vector2d",
            _make_sig("vector2d", [FLOAT, FLOAT], VECTOR2D),
            _vector2d,
        ))

    # --- Transform Functions ---

    def _register_transform_functions(self) -> None:
        """Register transform construction and solid transformation functions."""

        def _translate_xform(v: Value) -> Value:
            """Create a translation transform."""
            from yapcad.xform import Translation
            return transform_val(Translation(v.data))

        def _rotate_xform(axis: Value, angle: Value) -> Value:
            """Create a rotation transform."""
            from yapcad.xform import Rotation
            return transform_val(Rotation(axis.data, angle.data))

        def _scale_xform(factors: Value) -> Value:
            """Create a scale transform."""
            from yapcad.xform import Scale
            return transform_val(Scale(factors.data))

        def _identity_transform() -> Value:
            """Create an identity transform."""
            from yapcad.xform import Identity
            return transform_val(Identity())

        def _translate(s: Value, x: Value, y: Value, z: Value) -> Value:
            """Translate a solid by (x, y, z)."""
            from yapcad.geom3d import translatesolid
            from yapcad.geom import point
            delta = point(x.data, y.data, z.data)
            return solid_val(translatesolid(s.data, delta))

        def _rotate(s: Value, x: Value, y: Value, z: Value) -> Value:
            """Rotate a solid by angles (x, y, z) in degrees around respective axes."""
            from yapcad.geom3d import rotatesolid
            from yapcad.geom import point
            result = s.data
            # Apply rotations around each axis in sequence
            # Note: rotatesolid expects angles in degrees (xform.Rotation converts internally)
            if abs(x.data) > 1e-10:
                result = rotatesolid(result, x.data,
                                    cent=point(0, 0, 0), axis=point(1, 0, 0))
            if abs(y.data) > 1e-10:
                result = rotatesolid(result, y.data,
                                    cent=point(0, 0, 0), axis=point(0, 1, 0))
            if abs(z.data) > 1e-10:
                result = rotatesolid(result, z.data,
                                    cent=point(0, 0, 0), axis=point(0, 0, 1))
            return solid_val(result)

        def _scale(s: Value, x: Value, y: Value, z: Value) -> Value:
            """Scale a solid by factors (x, y, z)."""
            from yapcad.geom3d import scalesolid
            return solid_val(scalesolid(s.data, x.data, y.data, z.data))

        # Solid transformation functions (most commonly used)
        self.register(BuiltinFunction(
            "translate",
            _make_sig("translate", [SOLID, FLOAT, FLOAT, FLOAT], SOLID),
            _translate,
        ))
        self.register(BuiltinFunction(
            "rotate",
            _make_sig("rotate", [SOLID, FLOAT, FLOAT, FLOAT], SOLID),
            _rotate,
        ))
        self.register(BuiltinFunction(
            "scale",
            _make_sig("scale", [SOLID, FLOAT, FLOAT, FLOAT], SOLID),
            _scale,
        ))

        # Transform constructors (for advanced use)
        self.register(BuiltinFunction(
            "translate_xform",
            _make_sig("translate_xform", [VECTOR], TRANSFORM),
            _translate_xform,
        ))
        self.register(BuiltinFunction(
            "rotate_xform",
            _make_sig("rotate_xform", [VECTOR3D, FLOAT], TRANSFORM),
            _rotate_xform,
        ))
        self.register(BuiltinFunction(
            "scale_xform",
            _make_sig("scale_xform", [VECTOR], TRANSFORM),
            _scale_xform,
        ))
        self.register(BuiltinFunction(
            "identity_transform",
            _make_sig("identity_transform", [], TRANSFORM),
            _identity_transform,
        ))

        def _apply(t: Value, shape: Value) -> Value:
            """Apply a transform to a shape (solid, surface, point, etc.)."""
            from yapcad.geom3d import issolid, issurface, rotatesolid, rotatesurface
            from yapcad.geom import ispoint, isvect, transform as geom_transform
            from yapcad.geom import point
            from copy import deepcopy

            mat = t.data
            data = shape.data
            shape_type = shape.type

            # Handle solids by applying matrix to each surface
            if shape_type == SOLID or issolid(data):
                # Use rotatesolid with the pre-computed matrix
                # We pass a tiny non-zero angle to avoid early return, but the
                # actual transformation is done entirely by the provided matrix
                result = rotatesolid(data, 0.001, cent=point(0, 0, 0),
                                    axis=point(0, 0, 1.0), mat=mat)
                return solid_val(result)

            # Handle surfaces
            elif shape_type == SURFACE or issurface(data):
                result = rotatesurface(data, 0.001, cent=point(0, 0, 0),
                                      axis=point(0, 0, 1.0), mat=mat)
                return surface_val(result)

            # Handle points/vectors - use matrix multiplication directly
            elif shape_type in (POINT, POINT2D, POINT3D) or ispoint(data):
                result = mat.mul(data)
                return point_val(result)

            elif shape_type in (VECTOR, VECTOR2D, VECTOR3D) or isvect(data):
                result = mat.mul(data)
                return vector_val(result)

            # For other geometry types, use geom.transform
            else:
                try:
                    result = geom_transform(data, mat)
                    return wrap_value(result, shape_type)
                except (ValueError, NotImplementedError):
                    raise ValueError(
                        f"don't know how to apply transform to {shape_type}"
                    )

        # Register type-specific apply functions
        self.register(BuiltinFunction(
            "apply",
            _make_sig("apply", [TRANSFORM, SOLID], SOLID),
            _apply,
        ))
        self.register(BuiltinFunction(
            "apply_surface",
            _make_sig("apply_surface", [TRANSFORM, SURFACE], SURFACE),
            _apply,
        ))
        self.register(BuiltinFunction(
            "apply_point",
            _make_sig("apply_point", [TRANSFORM, POINT], POINT),
            _apply,
        ))
        self.register(BuiltinFunction(
            "apply_vector",
            _make_sig("apply_vector", [TRANSFORM, VECTOR3D], VECTOR3D),
            _apply,
        ))

    # --- Curve Functions ---

    def _register_curve_functions(self) -> None:
        """Register curve construction functions."""

        def _line(start: Value, end: Value) -> Value:
            """Create a line segment."""
            from yapcad.geom import line
            return wrap_value(line(start.data, end.data), LINE_SEGMENT)

        def _arc(center: Value, radius: Value, start_angle: Value, end_angle: Value) -> Value:
            """Create an arc."""
            from yapcad.geom import arc
            return wrap_value(
                arc(center.data, radius.data, start_angle.data, end_angle.data),
                ARC
            )

        def _circle(center: Value, radius: Value) -> Value:
            """Create a circle."""
            from yapcad.geom import arc
            # Full circle is an arc from 0 to 360
            return wrap_value(arc(center.data, radius.data, 0, 360), CIRCLE)

        # Curve constructors
        self.register(BuiltinFunction(
            "line",
            _make_sig("line", [POINT, POINT], LINE_SEGMENT),
            _line,
        ))
        self.register(BuiltinFunction(
            "arc",
            _make_sig("arc", [POINT, FLOAT, FLOAT, FLOAT], ARC),
            _arc,
        ))
        self.register(BuiltinFunction(
            "circle",
            _make_sig("circle", [POINT, FLOAT], CIRCLE),
            _circle,
        ))

    # --- Region Functions ---

    def _register_region_functions(self) -> None:
        """Register 2D region construction functions."""

        def _rectangle(width: Value, height: Value, *args: Value) -> Value:
            """Create a rectangle region2d."""
            from yapcad.geom import geom_util
            w, h = width.data, height.data
            # Optional center point
            cx, cy = 0.0, 0.0
            if args:
                center = args[0].data
                cx, cy = center[0], center[1]

            # Create rectangle as list of line segments
            half_w, half_h = w / 2, h / 2
            corners = [
                (cx - half_w, cy - half_h),
                (cx + half_w, cy - half_h),
                (cx + half_w, cy + half_h),
                (cx - half_w, cy + half_h),
            ]
            from yapcad.geom import point, line
            pts = [point(c[0], c[1]) for c in corners]
            region = [
                line(pts[0], pts[1]),
                line(pts[1], pts[2]),
                line(pts[2], pts[3]),
                line(pts[3], pts[0]),
            ]
            return region2d_val(region)

        def _regular_polygon(n: Value, radius: Value, *args: Value) -> Value:
            """Create a regular polygon region2d."""
            sides = int(n.data)
            r = radius.data
            cx, cy = 0.0, 0.0
            if args:
                center = args[0].data
                cx, cy = center[0], center[1]

            from yapcad.geom import point, line
            import math
            pts = []
            for i in range(sides):
                angle = 2 * math.pi * i / sides
                pts.append(point(cx + r * math.cos(angle), cy + r * math.sin(angle)))

            region = []
            for i in range(sides):
                region.append(line(pts[i], pts[(i + 1) % sides]))
            return region2d_val(region)

        # Region constructors
        self.register(BuiltinFunction(
            "rectangle",
            _make_sig("rectangle", [FLOAT, FLOAT], REGION2D),
            _rectangle,
        ))
        self.register(BuiltinFunction(
            "regular_polygon",
            _make_sig("regular_polygon", [INT, FLOAT], REGION2D),
            _regular_polygon,
        ))

    # --- Solid Functions ---

    def _register_solid_functions(self) -> None:
        """Register 3D solid construction functions."""

        def _box(width: Value, depth: Value, height: Value) -> Value:
            """Create a box solid (rectangular prism)."""
            from yapcad.geom3d_util import prism
            return solid_val(prism(width.data, depth.data, height.data))

        def _cylinder(radius: Value, height: Value) -> Value:
            """Create a solid cylinder using conic with equal radii."""
            from yapcad.geom3d_util import conic
            r = radius.data
            h = height.data
            # conic(base_radius, top_radius, height) - equal radii makes a cylinder
            return solid_val(conic(r, r, h))

        def _sphere(radius: Value) -> Value:
            """Create a sphere solid."""
            from yapcad.geom3d_util import sphere
            return solid_val(sphere(radius.data))

        def _cone(radius1: Value, radius2: Value, height: Value) -> Value:
            """Create a cone/frustum solid using conic."""
            from yapcad.geom3d_util import conic
            # conic(base_radius, top_radius, height)
            # radius2=0 makes a true cone, radius1!=radius2 makes a frustum
            return solid_val(conic(radius1.data, radius2.data, height.data))

        def _extrude(profile: Value, height: Value, *args: Value) -> Value:
            """Extrude a 2D region to create a solid."""
            from yapcad.geom3d_util import extrude
            # extrude takes surface, height
            return solid_val(extrude(profile.data, height.data))

        def _revolve(profile: Value, axis: Value, angle: Value) -> Value:
            """Revolve a 2D region around an axis."""
            from yapcad.geom3d_util import makeRevolutionSolid
            return solid_val(makeRevolutionSolid(profile.data, axis.data, angle.data))

        # Solid constructors
        self.register(BuiltinFunction(
            "box",
            _make_sig("box", [FLOAT, FLOAT, FLOAT], SOLID),
            _box,
        ))
        self.register(BuiltinFunction(
            "cylinder",
            _make_sig("cylinder", [FLOAT, FLOAT], SOLID),
            _cylinder,
        ))
        self.register(BuiltinFunction(
            "sphere",
            _make_sig("sphere", [FLOAT], SOLID),
            _sphere,
        ))
        self.register(BuiltinFunction(
            "cone",
            _make_sig("cone", [FLOAT, FLOAT, FLOAT], SOLID),
            _cone,
        ))
        self.register(BuiltinFunction(
            "extrude",
            _make_sig("extrude", [REGION2D, FLOAT], SOLID),
            _extrude,
        ))
        self.register(BuiltinFunction(
            "revolve",
            _make_sig("revolve", [REGION2D, VECTOR3D, FLOAT], SOLID),
            _revolve,
        ))

        def _involute_gear(teeth: Value, module_mm: Value, pressure_angle: Value,
                          face_width: Value) -> Value:
            """Create an involute spur gear solid.

            Args:
                teeth: Number of teeth (int)
                module_mm: Module in mm (metric gear sizing)
                pressure_angle: Pressure angle in degrees (typically 20)
                face_width: Thickness of the gear (extrusion height)

            Returns:
                A solid representing the gear
            """
            from yapcad.contrib.figgear import make_gear_figure
            from yapcad.geom3d import poly2surfaceXY
            from yapcad.geom3d_util import extrude
            from yapcad.geom import point

            # Generate 2D gear profile
            profile_points, blueprints = make_gear_figure(
                m=module_mm.data,
                z=int(teeth.data),
                alpha_deg=pressure_angle.data,
                bottom_type='spline',
            )

            # Convert to yapCAD points (z=0 for XY plane)
            pts = [point(x, y, 0.0) for x, y in profile_points]

            # Create surface from polygon
            surface, _ = poly2surfaceXY(pts)

            # Extrude to create solid
            gear_solid = extrude(surface, face_width.data)

            return solid_val(gear_solid)

        self.register(BuiltinFunction(
            "involute_gear",
            _make_sig("involute_gear", [INT, FLOAT, FLOAT, FLOAT], SOLID),
            _involute_gear,
        ))

    # --- Boolean Functions ---

    def _register_boolean_functions(self) -> None:
        """Register boolean operations."""

        def _union(*args: Value) -> Value:
            """Union of solids."""
            from yapcad.geom3d import solid_boolean
            if len(args) == 1 and isinstance(args[0].type, ListType):
                operands = args[0].data
            else:
                operands = [a.data for a in args]
            # Perform pairwise unions
            if len(operands) == 0:
                raise ValueError("union requires at least one solid")
            result = operands[0]
            for i in range(1, len(operands)):
                result = solid_boolean(result, operands[i], 'union')
            return solid_val(result)

        def _difference(a: Value, *rest: Value) -> Value:
            """Difference: subtract rest from a."""
            from yapcad.geom3d import solid_boolean
            if len(rest) == 1 and isinstance(rest[0].type, ListType):
                tools = rest[0].data
            else:
                tools = [r.data for r in rest]
            # Perform pairwise differences
            result = a.data
            for tool in tools:
                result = solid_boolean(result, tool, 'difference')
            return solid_val(result)

        def _intersection(*args: Value) -> Value:
            """Intersection of solids."""
            from yapcad.geom3d import solid_boolean
            if len(args) == 1 and isinstance(args[0].type, ListType):
                operands = args[0].data
            else:
                operands = [a.data for a in args]
            # Perform pairwise intersections
            if len(operands) == 0:
                raise ValueError("intersection requires at least one solid")
            result = operands[0]
            for i in range(1, len(operands)):
                result = solid_boolean(result, operands[i], 'intersection')
            return solid_val(result)

        # Boolean operations
        self.register(BuiltinFunction(
            "union",
            _make_sig("union", [SOLID], SOLID, is_variadic=True),
            _union,
        ))
        self.register(BuiltinFunction(
            "difference",
            _make_sig("difference", [SOLID, SOLID], SOLID, is_variadic=True),
            _difference,
        ))
        self.register(BuiltinFunction(
            "intersection",
            _make_sig("intersection", [SOLID], SOLID, is_variadic=True),
            _intersection,
        ))

    # --- Query Functions ---

    def _register_query_functions(self) -> None:
        """Register geometric query functions."""

        def _volume(s: Value) -> Value:
            """Calculate volume of a solid."""
            from yapcad.geom3d import volume
            return float_val(volume(s.data))

        def _surface_area(s: Value) -> Value:
            """Calculate surface area of a solid."""
            from yapcad.geom3d import surface_area
            return float_val(surface_area(s.data))

        def _area(r: Value) -> Value:
            """Calculate area of a region2d."""
            from yapcad.geom import area
            return float_val(area(r.data))

        def _perimeter(r: Value) -> Value:
            """Calculate perimeter of a region2d."""
            from yapcad.geom import perimeter
            return float_val(perimeter(r.data))

        # Query functions
        self.register(BuiltinFunction(
            "volume",
            _make_sig("volume", [SOLID], FLOAT),
            _volume,
        ))
        self.register(BuiltinFunction(
            "surface_area",
            _make_sig("surface_area", [SOLID], FLOAT),
            _surface_area,
        ))
        self.register(BuiltinFunction(
            "area",
            _make_sig("area", [REGION2D], FLOAT),
            _area,
        ))
        self.register(BuiltinFunction(
            "perimeter",
            _make_sig("perimeter", [REGION2D], FLOAT),
            _perimeter,
        ))

    # --- Utility Functions ---

    def _register_utility_functions(self) -> None:
        """Register utility functions."""

        def _len(lst: Value) -> Value:
            """Get length of a list."""
            return int_val(len(lst.data))

        def _range_func(*args: Value) -> Value:
            """Create a range as a list."""
            if len(args) == 1:
                # range(end)
                return list_val([int_val(i) for i in range(int(args[0].data))], INT)
            elif len(args) == 2:
                # range(start, end)
                return list_val([int_val(i) for i in range(int(args[0].data), int(args[1].data))], INT)
            else:
                # range(start, end, step)
                return list_val([int_val(i) for i in range(int(args[0].data), int(args[1].data), int(args[2].data))], INT)

        def _print_val(*args: Value) -> Value:
            """Print values (for debugging)."""
            print(" ".join(str(a.data) for a in args))
            return bool_val(True)

        def _concat(list1: Value, list2: Value) -> Value:
            """Concatenate two lists."""
            combined = list1.data + list2.data
            # Get element type from the list type info, not from data
            elem_type = list1.type.element_type
            # Create Value directly since combined contains raw data
            return Value(combined, ListType(elem_type))

        def _reverse(lst: Value) -> Value:
            """Reverse a list."""
            reversed_data = list(reversed(lst.data))
            # Get element type from the list type info
            elem_type = lst.type.element_type
            return Value(reversed_data, ListType(elem_type))

        def _flatten(lst: Value) -> Value:
            """Flatten a nested list."""
            result = []
            # For list<list<T>>, element type is list<T>, and we want T
            outer_elem_type = lst.type.element_type
            if isinstance(outer_elem_type, ListType):
                inner_elem_type = outer_elem_type.element_type
            else:
                inner_elem_type = UNKNOWN
            for item in lst.data:
                if isinstance(item, list):
                    result.extend(item)
                else:
                    result.append(item)
            return Value(result, ListType(inner_elem_type))

        # Utility functions
        self.register(BuiltinFunction(
            "len",
            _make_sig("len", [ListType(INT)], INT),  # Generic list
            _len,
        ))
        self.register(BuiltinFunction(
            "range",
            _make_sig("range", [INT], ListType(INT), is_variadic=True),
            _range_func,
        ))
        self.register(BuiltinFunction(
            "concat",
            _make_sig("concat", [ListType(UNKNOWN), ListType(UNKNOWN)], ListType(UNKNOWN)),
            _concat,
        ))
        self.register(BuiltinFunction(
            "reverse",
            _make_sig("reverse", [ListType(UNKNOWN)], ListType(UNKNOWN)),
            _reverse,
        ))
        self.register(BuiltinFunction(
            "flatten",
            _make_sig("flatten", [ListType(ListType(UNKNOWN))], ListType(UNKNOWN)),
            _flatten,
        ))
        self.register(BuiltinFunction(
            "print",
            _make_sig("print", [STRING], BOOL, is_variadic=True),
            _print_val,
        ))


# Global singleton registry
_registry: Optional[BuiltinRegistry] = None


def get_builtin_registry() -> BuiltinRegistry:
    """Get the global built-in function registry."""
    global _registry
    if _registry is None:
        _registry = BuiltinRegistry()
    return _registry


def call_builtin(name: str, args: List[Value]) -> Value:
    """
    Call a built-in function by name.

    Raises RuntimeError if function not found.
    """
    registry = get_builtin_registry()
    func = registry.get_function(name)
    if func is None:
        raise RuntimeError(f"Unknown built-in function: {name}")
    return func.implementation(*args)


def call_method(type_name: str, method_name: str, receiver: Value, args: List[Value]) -> Value:
    """
    Call a method on a value.

    Raises RuntimeError if method not found.
    """
    registry = get_builtin_registry()
    method = registry.get_method(type_name, method_name)
    if method is None:
        raise RuntimeError(f"Unknown method: {type_name}.{method_name}")
    return method.implementation(receiver, *args)
