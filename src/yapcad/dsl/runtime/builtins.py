"""
Built-in function registry for DSL interpreter.

Maps DSL function names to yapCAD implementations.

Implementation status (v1.0):
  - 2D curve primitives: DONE (ellipse, parabola, hyperbola implemented)
  - involute_gear: DONE (using yapcad.contrib.figgear)
  - Functional combinators: DONE (union_all, difference_all, sum, etc.)
  - Adaptive sweeps: DONE (sweep_adaptive, sweep_adaptive_hollow)

TODO for 1.0: Add DSL builtins for fasteners (Python API exists):
  - hex_bolt(standard, size, length) -> solid
    Uses yapcad.fasteners.metric_hex_cap_screw / unified_hex_cap_screw
  - hex_nut(standard, size) -> solid
    Uses yapcad.fasteners.metric_hex_nut / unified_hex_nut
  See: yapcad.fasteners, yapcad.threadgen for existing implementation

Future enhancements (v1.1+):
  - Move involute_gear to yapcad.stdlib.gears package for better namespace
  - Add tube/conic_tube/spherical_shell builtins from geom3d_util
  - Add optional 'centered' parameter to cylinder/cone primitives
"""

from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Optional, Tuple
import math

from .values import (
    Value, int_val, float_val, bool_val, string_val, list_val,
    point_val, vector_val, transform_val, solid_val, region2d_val,
    surface_val, shell_val, path3d_val, wrap_value, unwrap_value, unwrap_values,
)
from ..types import (
    Type, ListType,
    INT, FLOAT, BOOL, STRING, UNKNOWN,
    POINT, POINT2D, POINT3D, VECTOR, VECTOR2D, VECTOR3D, TRANSFORM,
    SOLID, REGION2D, SURFACE, SHELL,
    LINE_SEGMENT, ARC, CIRCLE, ELLIPSE, BEZIER, NURBS, CATMULLROM,
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
        self._register_fastener_functions()
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

        def _exp(x: Value) -> Value:
            # Clamp to avoid overflow
            val = x.data
            if val > 700:
                val = 700
            elif val < -700:
                val = -700
            return float_val(math.exp(val))

        def _log(x: Value) -> Value:
            return float_val(math.log(x.data))

        def _log10(x: Value) -> Value:
            return float_val(math.log10(x.data))

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
            ("exp", [FLOAT], FLOAT, _exp),
            ("log", [FLOAT], FLOAT, _log),
            ("log10", [FLOAT], FLOAT, _log10),
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

        # Ellipse constructor
        def _ellipse(center: Value, semi_major: Value, semi_minor: Value,
                     *args: Value) -> Value:
            """Create an ellipse curve.

            Args:
                center: Center point
                semi_major: Semi-major axis length
                semi_minor: Semi-minor axis length
                rotation: (optional) Rotation of major axis in degrees
                start: (optional) Start angle in degrees
                end: (optional) End angle in degrees

            Returns:
                An ellipse curve
            """
            from yapcad.geom import ellipse
            c = center.data
            a = semi_major.data
            b = semi_minor.data
            # Optional parameters
            rotation = args[0].data if len(args) > 0 else 0.0
            start = args[1].data if len(args) > 1 else 0.0
            end = args[2].data if len(args) > 2 else 360.0
            return wrap_value(
                ellipse(c, a, b, rotation=rotation, start=start, end=end),
                ELLIPSE
            )

        self.register(BuiltinFunction(
            "ellipse",
            _make_sig("ellipse", [POINT, FLOAT, FLOAT], ELLIPSE, is_variadic=True),
            _ellipse,
        ))

        # Catmull-Rom spline constructor
        def _catmullrom(points: Value, *args: Value) -> Value:
            """Create a Catmull-Rom spline curve.

            Args:
                points: List of control points
                closed: (optional) Whether the spline is closed (default false)
                alpha: (optional) Parameterization factor 0-1 (default 0.5 centripetal)

            Returns:
                A Catmull-Rom spline curve
            """
            # Extract control points as proper 4-element points [x, y, z, 1]
            pts = []
            for p in points.data:
                if hasattr(p, 'data'):
                    coords = list(p.data)
                else:
                    coords = list(p)
                # Ensure 4-element point with w=1
                if len(coords) == 3:
                    coords.append(1)
                elif len(coords) < 3:
                    coords.extend([0] * (3 - len(coords)))
                    coords.append(1)
                pts.append(coords[:4])

            # Optional parameters
            closed = args[0].data if len(args) > 0 else False
            alpha = args[1].data if len(args) > 1 else 0.5

            # Create Catmull-Rom spline structure
            curve = ['catmullrom', pts, {'closed': closed, 'alpha': alpha}]
            return wrap_value(curve, CATMULLROM)

        self.register(BuiltinFunction(
            "catmullrom",
            _make_sig("catmullrom", [ListType(POINT)], CATMULLROM, is_variadic=True),
            _catmullrom,
        ))

        # NURBS curve constructor
        def _nurbs(points: Value, *args: Value) -> Value:
            """Create a NURBS curve.

            Args:
                points: List of control points
                weights: (optional) List of weights (default all 1.0)
                degree: (optional) Curve degree (default 3)

            Returns:
                A NURBS curve
            """
            # Extract control points as proper 4-element points [x, y, z, 1]
            pts = []
            for p in points.data:
                if hasattr(p, 'data'):
                    coords = list(p.data)
                else:
                    coords = list(p)
                # Ensure 4-element point with w=1
                if len(coords) == 3:
                    coords.append(1)
                elif len(coords) < 3:
                    coords.extend([0] * (3 - len(coords)))
                    coords.append(1)
                pts.append(coords[:4])

            n = len(pts)

            # Optional parameters
            if len(args) > 0 and args[0].data is not None:
                weights = list(args[0].data)
            else:
                weights = [1.0] * n

            degree = int(args[1].data) if len(args) > 1 else 3

            # Generate open uniform knot vector
            # For n control points and degree p, we need n + p + 1 knots
            num_knots = n + degree + 1
            knots = []
            for i in range(num_knots):
                if i <= degree:
                    knots.append(0.0)
                elif i >= num_knots - degree - 1:
                    knots.append(1.0)
                else:
                    knots.append((i - degree) / (n - degree))

            # Create NURBS structure
            curve = ['nurbs', pts, {'degree': degree, 'weights': weights, 'knots': knots}]
            return wrap_value(curve, NURBS)

        self.register(BuiltinFunction(
            "nurbs",
            _make_sig("nurbs", [ListType(POINT)], NURBS, is_variadic=True),
            _nurbs,
        ))

        # Curve sampling functions
        def _sample_curve(curve: Value, t: Value) -> Value:
            """Sample a curve at parameter t in [0, 1].

            Works with line, arc, ellipse, catmullrom, and nurbs curves.
            """
            from yapcad.geom import point, isline, isarc, isellipse
            from yapcad.geom import sample as geom_sample, ellipse_sample
            from yapcad.spline import is_catmullrom, evaluate_catmullrom
            from yapcad.spline import is_nurbs, evaluate_nurbs

            data = curve.data
            u = t.data

            if is_catmullrom(data):
                pt = evaluate_catmullrom(data, u)
                return point_val(pt, is_2d=False)
            elif is_nurbs(data):
                pt = evaluate_nurbs(data, u)
                return point_val(pt, is_2d=False)
            elif isellipse(data):
                pt = ellipse_sample(data, u)
                return point_val(pt, is_2d=False)
            elif isline(data) or isarc(data):
                pt = geom_sample(data, u)
                return point_val(pt, is_2d=False)
            else:
                raise ValueError(f"Cannot sample curve of type {type(data)}")

        def _sample_curve_n(curve: Value, n: Value) -> Value:
            """Sample a curve at n evenly spaced points.

            Returns a list of points along the curve.
            """
            from yapcad.geom import point, isline, isarc, isellipse
            from yapcad.geom import sample as geom_sample, ellipse_sample
            from yapcad.spline import is_catmullrom, evaluate_catmullrom
            from yapcad.spline import is_nurbs, evaluate_nurbs

            data = curve.data
            num = int(n.data)
            if num < 2:
                raise ValueError("n must be at least 2")

            pts = []
            for i in range(num):
                u = i / (num - 1)
                if is_catmullrom(data):
                    pt = evaluate_catmullrom(data, u)
                elif is_nurbs(data):
                    pt = evaluate_nurbs(data, u)
                elif isellipse(data):
                    pt = ellipse_sample(data, u)
                elif isline(data) or isarc(data):
                    pt = geom_sample(data, u)
                else:
                    raise ValueError(f"Cannot sample curve of type {type(data)}")
                pts.append(point_val(pt, is_2d=False))

            return list_val(pts, POINT)

        def _curve_length(curve: Value) -> Value:
            """Compute the approximate length of a curve."""
            from yapcad.geom import length as geom_length, isline, isarc, isellipse
            from yapcad.geom import ellipse_length
            from yapcad.spline import is_catmullrom, sample_catmullrom
            from yapcad.spline import is_nurbs, sample_nurbs

            data = curve.data

            if is_catmullrom(data):
                # Sample and compute polyline length
                pts = sample_catmullrom(data, segments_per_span=20)
                total = 0.0
                for i in range(1, len(pts)):
                    dx = pts[i][0] - pts[i-1][0]
                    dy = pts[i][1] - pts[i-1][1]
                    dz = pts[i][2] - pts[i-1][2]
                    total += math.sqrt(dx*dx + dy*dy + dz*dz)
                return float_val(total)
            elif is_nurbs(data):
                pts = sample_nurbs(data, samples=100)
                total = 0.0
                for i in range(1, len(pts)):
                    dx = pts[i][0] - pts[i-1][0]
                    dy = pts[i][1] - pts[i-1][1]
                    dz = pts[i][2] - pts[i-1][2]
                    total += math.sqrt(dx*dx + dy*dy + dz*dz)
                return float_val(total)
            elif isellipse(data):
                return float_val(ellipse_length(data))
            elif isline(data) or isarc(data):
                return float_val(geom_length(data))
            else:
                raise ValueError(f"Cannot compute length for curve of type {type(data)}")

        self.register(BuiltinFunction(
            "sample_curve",
            _make_sig("sample_curve", [UNKNOWN, FLOAT], POINT),
            _sample_curve,
        ))
        self.register(BuiltinFunction(
            "sample_curve_n",
            _make_sig("sample_curve_n", [UNKNOWN, INT], ListType(POINT)),
            _sample_curve_n,
        ))
        self.register(BuiltinFunction(
            "curve_length",
            _make_sig("curve_length", [UNKNOWN], FLOAT),
            _curve_length,
        ))

        # Path3D constructors for sweep operations
        def _make_path3d(*segments: Value) -> Value:
            """Create a path3d from segments or other path3d objects.

            Each argument can be:
            - A path3d (dict with 'type': 'path3d' and 'segments' list)
            - A segment dict (dict with 'type': 'line' or 'arc')
            - A list of segments
            """
            seg_list = []
            for seg in segments:
                data = seg.data
                if isinstance(data, dict):
                    if data.get('type') == 'path3d':
                        # It's a path3d, extract its segments
                        seg_list.extend(data.get('segments', []))
                    else:
                        # It's a single segment dict
                        seg_list.append(data)
                elif isinstance(data, list):
                    # List of segments
                    seg_list.extend(data)
                else:
                    seg_list.append(data)
            return path3d_val({'type': 'path3d', 'segments': seg_list})

        def _path3d_line(start: Value, end: Value) -> Value:
            """Create a path3d containing a single line segment.

            Args:
                start: Start point [x, y, z]
                end: End point [x, y, z]

            Returns:
                A path3d dict with a single line segment
            """
            s = start.data
            e = end.data
            segment = {
                'type': 'line',
                'start': [s[0], s[1], s[2] if len(s) > 2 else 0],
                'end': [e[0], e[1], e[2] if len(e) > 2 else 0]
            }
            return wrap_value({
                'type': 'path3d',
                'segments': [segment]
            }, PATH3D)

        def _path3d_arc(center: Value, start: Value, end: Value, normal: Value) -> Value:
            """Create a path3d containing a single arc segment.

            Args:
                center: Arc center point [x, y, z]
                start: Arc start point [x, y, z]
                end: Arc end point [x, y, z]
                normal: Arc plane normal [nx, ny, nz]

            Returns:
                A path3d dict with a single arc segment
            """
            c = center.data
            s = start.data
            e = end.data
            n = normal.data
            segment = {
                'type': 'arc',
                'center': [c[0], c[1], c[2] if len(c) > 2 else 0],
                'start': [s[0], s[1], s[2] if len(s) > 2 else 0],
                'end': [e[0], e[1], e[2] if len(e) > 2 else 0],
                'normal': [n[0], n[1], n[2] if len(n) > 2 else 0]
            }
            return wrap_value({
                'type': 'path3d',
                'segments': [segment]
            }, PATH3D)

        self.register(BuiltinFunction(
            "make_path3d",
            _make_sig("make_path3d", [], PATH3D, is_variadic=True),
            _make_path3d,
        ))
        self.register(BuiltinFunction(
            "path3d_line",
            _make_sig("path3d_line", [POINT3D, POINT3D], PATH3D),
            _path3d_line,
        ))
        self.register(BuiltinFunction(
            "path3d_arc",
            _make_sig("path3d_arc", [POINT3D, POINT3D, POINT3D, VECTOR3D], PATH3D),
            _path3d_arc,
        ))

        def _path3d_arc_auto(center: Value, start: Value, end: Value, flip: Value) -> Value:
            """Create a path3d arc with auto-computed normal from geometry.

            The arc plane normal is computed from the cross product of
            (center->start) x (center->end). The flip parameter controls
            which of the two possible arc directions is used.

            Args:
                center: Arc center point [x, y, z]
                start: Arc start point [x, y, z]
                end: Arc end point [x, y, z]
                flip: If true, negate the computed normal (reverses arc direction)

            Returns:
                A path3d dict with a single arc segment
            """
            import math

            c = center.data
            s = start.data
            e = end.data
            do_flip = flip.data

            # Extract coordinates
            cx, cy, cz = c[0], c[1], c[2] if len(c) > 2 else 0
            sx, sy, sz = s[0], s[1], s[2] if len(s) > 2 else 0
            ex, ey, ez = e[0], e[1], e[2] if len(e) > 2 else 0

            # Vectors from center to start and end
            v1 = [sx - cx, sy - cy, sz - cz]
            v2 = [ex - cx, ey - cy, ez - cz]

            # Cross product: v1 x v2
            nx = v1[1] * v2[2] - v1[2] * v2[1]
            ny = v1[2] * v2[0] - v1[0] * v2[2]
            nz = v1[0] * v2[1] - v1[1] * v2[0]

            # Normalize
            mag = math.sqrt(nx*nx + ny*ny + nz*nz)
            if mag < 1e-10:
                # Degenerate case - points are collinear, use Z-up default
                nx, ny, nz = 0.0, 0.0, 1.0
            else:
                nx, ny, nz = nx/mag, ny/mag, nz/mag

            # Apply flip if requested
            if do_flip:
                nx, ny, nz = -nx, -ny, -nz

            segment = {
                'type': 'arc',
                'center': [cx, cy, cz],
                'start': [sx, sy, sz],
                'end': [ex, ey, ez],
                'normal': [nx, ny, nz]
            }
            return wrap_value({
                'type': 'path3d',
                'segments': [segment]
            }, PATH3D)

        self.register(BuiltinFunction(
            "path3d_arc_auto",
            _make_sig("path3d_arc_auto", [POINT3D, POINT3D, POINT3D, BOOL], PATH3D),
            _path3d_arc_auto,
        ))

    # --- Region Functions ---

    def _register_region_functions(self) -> None:
        """Register 2D region construction functions."""

        def _rectangle(width: Value, height: Value, *args: Value) -> Value:
            """Create a rectangle region2d."""
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

        # Polygon from points
        def _polygon(points: Value) -> Value:
            """Create a closed polygon region from a list of points.

            Args:
                points: List of 2D points defining the polygon vertices

            Returns:
                A closed region2d
            """
            from yapcad.geom import point, line

            pts = []
            for p in points.data:
                if hasattr(p, 'data'):
                    pts.append(point(p.data[0], p.data[1]))
                else:
                    pts.append(point(p[0], p[1]))

            if len(pts) < 3:
                raise ValueError("polygon requires at least 3 points")

            # Create closed polygon from line segments
            region = []
            for i in range(len(pts)):
                region.append(line(pts[i], pts[(i + 1) % len(pts)]))
            return region2d_val(region)

        self.register(BuiltinFunction(
            "polygon",
            _make_sig("polygon", [ListType(POINT)], REGION2D),
            _polygon,
        ))

        # Disk (filled circle) region
        def _disk(center: Value, radius: Value, *args: Value) -> Value:
            """Create a filled circular region (disk).

            Args:
                center: Center point
                radius: Radius of the disk
                segments: (optional) Number of sides for polygon approximation (default 64)

            Returns:
                A region2d approximating a disk
            """
            from yapcad.geom import point, line

            c = center.data
            r = radius.data
            cx, cy = c[0], c[1]
            segments = int(args[0].data) if args else 64

            pts = []
            for i in range(segments):
                angle = 2 * math.pi * i / segments
                pts.append(point(cx + r * math.cos(angle), cy + r * math.sin(angle)))

            region = []
            for i in range(segments):
                region.append(line(pts[i], pts[(i + 1) % segments]))
            return region2d_val(region)

        self.register(BuiltinFunction(
            "disk",
            _make_sig("disk", [POINT, FLOAT], REGION2D, is_variadic=True),
            _disk,
        ))

        # 2D Boolean operations
        def _union2d(a: Value, b: Value) -> Value:
            """Boolean union of two 2D regions.

            Args:
                a: First region2d
                b: Second region2d

            Returns:
                Union of the two regions as a region2d
            """
            from yapcad.geom_util import combineglist
            result = combineglist(a.data, b.data, "union")
            return region2d_val(result)

        def _difference2d(a: Value, b: Value) -> Value:
            """Boolean difference of two 2D regions (a minus b).

            Handles chained difference operations where a may already contain holes
            from a previous difference operation. Also works around a known issue
            in combineglist where isinsideXY can fail for points near corners.

            Args:
                a: Region to subtract from (may be region with existing holes)
                b: Region to subtract

            Returns:
                Difference (a - b) as a region2d
            """
            from yapcad.geom import (isgeomlist, isline, isarc, bbox, isinsidebbox,
                                     sample, intersectXY)
            from yapcad.geom_util import combineglist

            a_data = a.data
            b_data = b.data

            # Detect if a_data is a "region with holes" structure:
            # [outer_boundary, hole1, hole2, ...]
            # where each element is itself a geomlist of primitives
            def is_region_with_holes(data):
                """Check if data is a list of geomlists (region with holes)."""
                if not isinstance(data, list) or len(data) < 2:
                    return False
                first = data[0]
                if not isinstance(first, list) or len(first) == 0:
                    return False
                first_elem = first[0]
                if isline(first_elem) or isarc(first_elem):
                    return True
                return False

            def is_b_inside_a(outer, inner):
                """Check if inner is completely inside outer using robust method.

                Uses bounding box check + intersection check to avoid isinsideXY
                corner issues.
                """
                bbox_outer = bbox(outer)
                bbox_inner = bbox(inner)
                if not bbox_outer or not bbox_inner:
                    return False
                # Inner bbox must be inside outer bbox
                if not (isinsidebbox(bbox_outer, bbox_inner[0]) and
                        isinsidebbox(bbox_outer, bbox_inner[1])):
                    return False
                # Check for intersections - if none, inner is either fully inside or outside
                inter = intersectXY(outer, inner, params=True)
                if inter is False or (inter[0] == [] and inter[1] == []):
                    # No intersections - inner is fully inside (already passed bbox check)
                    return True
                return False

            def do_difference(outer, hole):
                """Perform difference, handling the case where hole is inside outer."""
                result = combineglist(outer, hole, "difference")
                # Check if combineglist returned outer unchanged (missed the hole)
                if isgeomlist(result) and not is_region_with_holes(result):
                    # Result is flat geomlist - check if hole was missed
                    if is_b_inside_a(outer, hole):
                        # Hole should have created [outer, hole] but didn't
                        return [outer, hole]
                return result

            if is_region_with_holes(a_data):
                # a_data is [outer, hole1, hole2, ...]
                outer = a_data[0]
                existing_holes = a_data[1:]

                result = do_difference(outer, b_data)

                if is_region_with_holes(result):
                    new_outer = result[0]
                    new_holes = result[1:]
                    return region2d_val([new_outer] + existing_holes + new_holes)
                else:
                    if result:
                        return region2d_val([result] + existing_holes)
                    else:
                        return region2d_val([])
            else:
                result = do_difference(a_data, b_data)
                return region2d_val(result)

        def _intersection2d(a: Value, b: Value) -> Value:
            """Boolean intersection of two 2D regions.

            Args:
                a: First region2d
                b: Second region2d

            Returns:
                Intersection of the two regions as a region2d
            """
            from yapcad.geom_util import combineglist
            result = combineglist(a.data, b.data, "intersection")
            return region2d_val(result)

        self.register(BuiltinFunction(
            "union2d",
            _make_sig("union2d", [REGION2D, REGION2D], REGION2D),
            _union2d,
        ))
        self.register(BuiltinFunction(
            "difference2d",
            _make_sig("difference2d", [REGION2D, REGION2D], REGION2D),
            _difference2d,
        ))
        self.register(BuiltinFunction(
            "intersection2d",
            _make_sig("intersection2d", [REGION2D, REGION2D], REGION2D),
            _intersection2d,
        ))

        # Phase 3: 2D boolean aggregation functions
        def _union2d_all(regions: Value) -> Value:
            """Union all 2D regions in a list."""
            from yapcad.geom_util import combineglist
            if not regions.data:
                raise ValueError("union2d_all requires non-empty list")
            result = regions.data[0]
            if hasattr(result, 'data'):
                result = result.data
            for r in regions.data[1:]:
                r_data = r.data if hasattr(r, 'data') else r
                result = combineglist(result, r_data, "union")
            return region2d_val(result)

        def _difference2d_all(base: Value, tools: Value) -> Value:
            """Subtract all 2D regions in tools list from base."""
            result = base.data
            for tool in tools.data:
                tool_data = tool.data if hasattr(tool, 'data') else tool
                # Use the same robust difference logic as _difference2d
                result_val = _difference2d(
                    Value(result, REGION2D),
                    Value(tool_data, REGION2D)
                )
                result = result_val.data
            return region2d_val(result)

        def _intersection2d_all(regions: Value) -> Value:
            """Intersect all 2D regions in a list."""
            from yapcad.geom_util import combineglist
            if not regions.data:
                raise ValueError("intersection2d_all requires non-empty list")
            result = regions.data[0]
            if hasattr(result, 'data'):
                result = result.data
            for r in regions.data[1:]:
                r_data = r.data if hasattr(r, 'data') else r
                result = combineglist(result, r_data, "intersection")
            return region2d_val(result)

        self.register(BuiltinFunction(
            "union2d_all",
            _make_sig("union2d_all", [ListType(REGION2D)], REGION2D),
            _union2d_all,
        ))
        self.register(BuiltinFunction(
            "difference2d_all",
            _make_sig("difference2d_all", [REGION2D, ListType(REGION2D)], REGION2D),
            _difference2d_all,
        ))
        self.register(BuiltinFunction(
            "intersection2d_all",
            _make_sig("intersection2d_all", [ListType(REGION2D)], REGION2D),
            _intersection2d_all,
        ))

        # Path2D and region from curves
        def _make_path2d(curves: Value) -> Value:
            """Create a 2D path from a list of curves.

            Args:
                curves: List of curves (line, arc, etc.)

            Returns:
                A path2d (open path of curves)
            """
            # Path2D is just a geometry list
            curve_list = []
            for c in curves.data:
                if hasattr(c, 'data'):
                    curve_list.append(c.data)
                else:
                    curve_list.append(c)
            return wrap_value(curve_list, PATH2D)

        def _close_path(path: Value) -> Value:
            """Close an open path to create a region.

            Args:
                path: An open path2d

            Returns:
                A closed region2d
            """
            from yapcad.geom import line, point

            curves = path.data if isinstance(path.data, list) else [path.data]

            if len(curves) == 0:
                raise ValueError("Cannot close empty path")

            # Check if already closed
            first_curve = curves[0]
            last_curve = curves[-1]

            # Get start of first curve and end of last curve
            from yapcad.geom import isline, isarc

            def get_endpoints(curve):
                if isline(curve):
                    return curve[0], curve[1]
                elif isarc(curve):
                    from yapcad.geom import sample
                    return sample(curve, 0), sample(curve, 1)
                else:
                    # For other curves, sample endpoints
                    from yapcad.geom import sample
                    return sample(curve, 0), sample(curve, 1)

            first_start, _ = get_endpoints(first_curve)
            _, last_end = get_endpoints(last_curve)

            # Check if endpoints are close enough
            from yapcad.geom import dist
            if dist(first_start, last_end) > 1e-6:
                # Add closing segment
                closing_line = line(last_end, first_start)
                curves = list(curves) + [closing_line]

            return region2d_val(curves)

        self.register(BuiltinFunction(
            "make_path2d",
            _make_sig("make_path2d", [ListType(UNKNOWN)], PATH2D),
            _make_path2d,
        ))
        self.register(BuiltinFunction(
            "close_path",
            _make_sig("close_path", [PATH2D], REGION2D),
            _close_path,
        ))

        # Region from sampled spline
        def _region_from_spline(spline: Value, *args: Value) -> Value:
            """Convert a spline curve to a closed polygon region.

            Samples the spline and creates a polygon from the sample points.
            Useful for creating regions from Catmull-Rom or NURBS curves.

            Args:
                spline: A catmullrom or nurbs curve (must be closed or will be auto-closed)
                segments: (optional) Number of sample points (default 64)

            Returns:
                A region2d polygon approximating the spline
            """
            from yapcad.geom import point, line
            from yapcad.spline import is_catmullrom, sample_catmullrom
            from yapcad.spline import is_nurbs, sample_nurbs

            data = spline.data
            segments = int(args[0].data) if args else 64

            if is_catmullrom(data):
                pts = sample_catmullrom(data, segments_per_span=segments // max(1, len(data[1]) - 1))
            elif is_nurbs(data):
                pts = sample_nurbs(data, samples=segments)
            else:
                raise ValueError("region_from_spline requires a catmullrom or nurbs curve")

            if len(pts) < 3:
                raise ValueError("Spline produced too few sample points")

            # Create polygon from sample points
            region = []
            for i in range(len(pts)):
                p1 = point(pts[i][0], pts[i][1])
                p2 = point(pts[(i + 1) % len(pts)][0], pts[(i + 1) % len(pts)][1])
                region.append(line(p1, p2))

            return region2d_val(region)

        self.register(BuiltinFunction(
            "region_from_spline",
            _make_sig("region_from_spline", [UNKNOWN], REGION2D, is_variadic=True),
            _region_from_spline,
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

        def _oblate_spheroid(equatorial_diameter: Value, oblateness: Value) -> Value:
            """Create an oblate spheroid (flattened sphere).

            An oblate spheroid has equal X and Y radii (equatorial) and a smaller
            Z radius (polar). The oblateness parameter controls the flattening:
            0 = perfect sphere, higher values = more flattened.

            For reference:
            - Earth's oblateness: ~0.00335
            - Mars' oblateness: ~0.00648

            Args:
                equatorial_diameter: diameter at the equator (X and Y)
                oblateness: geometric flattening (0-1, typically small like 0.006)

            Returns:
                A solid oblate spheroid centered at origin with polar axis along Z
            """
            from yapcad.geom3d_util import oblate_spheroid
            return solid_val(oblate_spheroid(equatorial_diameter.data, oblateness.data))

        def _cone(radius1: Value, radius2: Value, height: Value) -> Value:
            """Create a cone/frustum solid using conic."""
            from yapcad.geom3d_util import conic
            # conic(base_radius, top_radius, height)
            # radius2=0 makes a true cone, radius1!=radius2 makes a frustum
            return solid_val(conic(radius1.data, radius2.data, height.data))

        def _extrude(profile: Value, height: Value, *args: Value) -> Value:
            """Extrude a 2D region to create a solid.

            Args:
                profile: A region2d (closed 2D shape) to extrude
                height: Extrusion height along Z axis

            Returns:
                A solid created by extruding the profile
            """
            from yapcad.geom3d_util import extrude_region2d
            return solid_val(extrude_region2d(profile.data, height.data))

        def _revolve(profile: Value, axis: Value, angle: Value) -> Value:
            """Revolve a 2D region around an axis."""
            from yapcad.geom3d_util import makeRevolutionSolid
            return solid_val(makeRevolutionSolid(profile.data, axis.data, angle.data))

        def _sweep(profile: Value, spine: Value) -> Value:
            """Sweep a 2D profile along a 3D path to create a solid.

            Args:
                profile: A region2d (closed 2D shape) to sweep
                spine: A path3d (3D wire/curve) defining the sweep path

            Returns:
                A solid created by sweeping the profile along the spine
            """
            from yapcad.geom3d_util import sweep_profile_along_path
            return solid_val(sweep_profile_along_path(profile.data, spine.data))

        def _sweep_hollow(outer_profile: Value, inner_profile: Value, spine: Value) -> Value:
            """Sweep a hollow 2D profile along a 3D path to create a solid.

            Creates a hollow tube-like solid by sweeping a profile with a hole.
            The outer_profile defines the outer boundary, inner_profile the hole.

            Args:
                outer_profile: A region2d for the outer boundary
                inner_profile: A region2d for the inner boundary (hole)
                spine: A path3d (3D wire/curve) defining the sweep path

            Returns:
                A hollow solid created by sweeping the profile along the spine
            """
            from yapcad.geom3d_util import sweep_profile_along_path
            return solid_val(sweep_profile_along_path(
                outer_profile.data, spine.data, inner_profile=inner_profile.data
            ))

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
            "oblate_spheroid",
            _make_sig("oblate_spheroid", [FLOAT, FLOAT], SOLID),
            _oblate_spheroid,
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
        self.register(BuiltinFunction(
            "sweep",
            _make_sig("sweep", [REGION2D, PATH3D], SOLID),
            _sweep,
        ))
        self.register(BuiltinFunction(
            "sweep_hollow",
            _make_sig("sweep_hollow", [REGION2D, REGION2D, PATH3D], SOLID),
            _sweep_hollow,
        ))

        def _sweep_adaptive(profile: Value, spine: Value,
                            threshold: Value) -> Value:
            """Sweep a profile along a path with adaptive tangent tracking.

            The profile normal tracks the path tangent. New profile sections
            are generated whenever the tangent direction changes by more than
            the threshold angle.

            Args:
                profile: A region2d (closed 2D shape) to sweep
                spine: A path3d (3D wire/curve) defining the sweep path
                threshold: Angle in degrees that triggers new section (e.g., 5.0)

            Returns:
                A solid created by lofting between adapted sections
            """
            from yapcad.geom3d_util import sweep_adaptive
            return solid_val(sweep_adaptive(
                profile.data, spine.data,
                angle_threshold_deg=threshold.data
            ))

        def _sweep_adaptive_hollow(outer_profile: Value, inner_profiles: Value,
                                   spine: Value, threshold: Value) -> Value:
            """Sweep a hollow profile along a path with adaptive tangent tracking.

            Args:
                outer_profile: A region2d for the outer boundary
                inner_profiles: A region2d or list of region2d for inner void(s)
                spine: A path3d defining the sweep path
                threshold: Angle in degrees that triggers new section

            Returns:
                A hollow solid created by lofting between adapted sections
            """
            from yapcad.geom3d_util import sweep_adaptive
            # Handle single region2d or list
            inners = inner_profiles.data
            return solid_val(sweep_adaptive(
                outer_profile.data, spine.data,
                inner_profiles=inners,
                angle_threshold_deg=threshold.data
            ))

        self.register(BuiltinFunction(
            "sweep_adaptive",
            _make_sig("sweep_adaptive", [REGION2D, PATH3D, FLOAT], SOLID),
            _sweep_adaptive,
        ))
        self.register(BuiltinFunction(
            "sweep_adaptive_hollow",
            _make_sig("sweep_adaptive_hollow", [REGION2D, REGION2D, PATH3D, FLOAT], SOLID),
            _sweep_adaptive_hollow,
        ))

        def _sweep_adaptive_frenet(profile: Value, spine: Value,
                                   threshold: Value) -> Value:
            """Sweep a profile using Frenet frame (natural curvature-following).

            The profile orientation follows the natural Frenet frame of the path,
            where the profile 'up' direction aligns with the curve's normal
            (perpendicular to both tangent and binormal).

            This mode causes the profile to twist naturally with the path curvature,
            which is appropriate for paths like helices where you want the profile
            to follow the curve's intrinsic geometry.

            Args:
                profile: A region2d (closed 2D shape) to sweep
                spine: A path3d (3D wire/curve) defining the sweep path
                threshold: Angle in degrees that triggers new section

            Returns:
                A solid created by lofting with Frenet frame orientation
            """
            from yapcad.geom3d_util import sweep_adaptive
            return solid_val(sweep_adaptive(
                profile.data, spine.data,
                angle_threshold_deg=threshold.data,
                frame_mode='frenet'
            ))

        def _sweep_adaptive_hollow_frenet(outer_profile: Value, inner_profiles: Value,
                                          spine: Value, threshold: Value) -> Value:
            """Sweep a hollow profile using Frenet frame orientation.

            Args:
                outer_profile: A region2d for the outer boundary
                inner_profiles: A region2d or list of region2d for inner void(s)
                spine: A path3d defining the sweep path
                threshold: Angle in degrees that triggers new section

            Returns:
                A hollow solid with Frenet frame orientation
            """
            from yapcad.geom3d_util import sweep_adaptive
            inners = inner_profiles.data
            return solid_val(sweep_adaptive(
                outer_profile.data, spine.data,
                inner_profiles=inners,
                angle_threshold_deg=threshold.data,
                frame_mode='frenet'
            ))

        self.register(BuiltinFunction(
            "sweep_adaptive_frenet",
            _make_sig("sweep_adaptive_frenet", [REGION2D, PATH3D, FLOAT], SOLID),
            _sweep_adaptive_frenet,
        ))
        self.register(BuiltinFunction(
            "sweep_adaptive_hollow_frenet",
            _make_sig("sweep_adaptive_hollow_frenet", [REGION2D, REGION2D, PATH3D, FLOAT], SOLID),
            _sweep_adaptive_hollow_frenet,
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

    # --- Fastener Functions ---

    def _register_fastener_functions(self) -> None:
        """Register fastener generation builtins."""

        def _metric_hex_bolt(size: Value, length: Value) -> Value:
            """Create a metric hex bolt (ISO 4014/4017).

            Args:
                size: Thread size designation (e.g., "M8", "M10")
                length: Total shank length in mm

            Returns:
                yapCAD solid representing the bolt
            """
            from yapcad.fasteners import metric_hex_bolt
            return solid_val(metric_hex_bolt(size.data, length.data))

        def _metric_hex_nut(size: Value) -> Value:
            """Create a metric hex nut (ISO 4032).

            Args:
                size: Thread size designation (e.g., "M8", "M10")

            Returns:
                yapCAD solid representing the nut
            """
            from yapcad.fasteners import metric_hex_nut
            return solid_val(metric_hex_nut(size.data))

        def _unified_hex_bolt(size: Value, length: Value) -> Value:
            """Create a unified (UNC/UNF) hex bolt (ASME B18.2.1).

            Args:
                size: Thread size designation (e.g., "1/4-20", "#10-24")
                length: Total shank length in inches

            Returns:
                yapCAD solid representing the bolt
            """
            from yapcad.fasteners import unified_hex_bolt
            return solid_val(unified_hex_bolt(size.data, length.data))

        def _unified_hex_nut(size: Value) -> Value:
            """Create a unified (UNC/UNF) hex nut (ASME B18.2.2).

            Args:
                size: Thread size designation (e.g., "1/4-20", "#10-24")

            Returns:
                yapCAD solid representing the nut
            """
            from yapcad.fasteners import unified_hex_nut
            return solid_val(unified_hex_nut(size.data))

        self.register(BuiltinFunction(
            "metric_hex_bolt",
            _make_sig("metric_hex_bolt", [STRING, FLOAT], SOLID),
            _metric_hex_bolt,
        ))

        self.register(BuiltinFunction(
            "metric_hex_nut",
            _make_sig("metric_hex_nut", [STRING], SOLID),
            _metric_hex_nut,
        ))

        self.register(BuiltinFunction(
            "unified_hex_bolt",
            _make_sig("unified_hex_bolt", [STRING, FLOAT], SOLID),
            _unified_hex_bolt,
        ))

        self.register(BuiltinFunction(
            "unified_hex_nut",
            _make_sig("unified_hex_nut", [STRING], SOLID),
            _unified_hex_nut,
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

        def _compound(*args: Value) -> Value:
            """Create a compound of multiple solids without boolean operations.

            Unlike union(), this does NOT merge the solids - they remain as
            separate bodies in a compound shape. Useful for assemblies where
            you want to export multiple parts together without merging them.

            Args:
                *args: One or more solids to combine into a compound

            Returns:
                A solid containing all input solids as separate bodies
            """
            from yapcad.brep import occ_available, BrepSolid, attach_brep_to_solid, brep_from_solid
            from yapcad.geom3d import solid

            if len(args) == 1 and isinstance(args[0].type, ListType):
                operands = args[0].data
            else:
                operands = [a.data for a in args]

            if len(operands) == 0:
                raise ValueError("compound requires at least one solid")

            # Collect all surfaces from all solids for the mesh representation
            all_surfaces = []
            for op in operands:
                if len(op) > 1 and op[1]:
                    all_surfaces.extend(op[1])

            # Create result solid with combined surfaces
            result = solid(all_surfaces, [], ['procedure', 'compound'])

            # If OCC available, create a proper compound shape
            if occ_available():
                try:
                    from OCC.Core.TopoDS import TopoDS_Compound
                    from OCC.Core.BRep import BRep_Builder
                    from OCC.Core.TopExp import TopExp_Explorer
                    from OCC.Core.TopAbs import TopAbs_SOLID
                    from OCC.Core.TopoDS import topods

                    compound = TopoDS_Compound()
                    builder = BRep_Builder()
                    builder.MakeCompound(compound)

                    # Add each solid's BREP to the compound
                    for op in operands:
                        brep = brep_from_solid(op)
                        if brep is not None and brep.shape is not None:
                            # Extract all solids from this operand
                            exp = TopExp_Explorer(brep.shape, TopAbs_SOLID)
                            while exp.More():
                                builder.Add(compound, exp.Current())
                                exp.Next()

                    attach_brep_to_solid(result, BrepSolid(compound))
                except Exception:
                    pass

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
        self.register(BuiltinFunction(
            "compound",
            _make_sig("compound", [SOLID], SOLID, is_variadic=True),
            _compound,
        ))

        # Phase 3: List-based aggregation functions
        # These are aliases - the existing functions already handle list arguments
        def _union_all(solids: Value) -> Value:
            """Union all solids in a list."""
            return _union(solids)

        def _difference_all(base: Value, tools: Value) -> Value:
            """Subtract all solids in tools list from base."""
            from yapcad.geom3d import solid_boolean
            result = base.data
            for tool in tools.data:
                tool_data = tool.data if hasattr(tool, 'data') else tool
                result = solid_boolean(result, tool_data, 'difference')
            return solid_val(result)

        def _intersection_all(solids: Value) -> Value:
            """Intersect all solids in a list."""
            return _intersection(solids)

        self.register(BuiltinFunction(
            "union_all",
            _make_sig("union_all", [ListType(SOLID)], SOLID),
            _union_all,
        ))
        self.register(BuiltinFunction(
            "difference_all",
            _make_sig("difference_all", [SOLID, ListType(SOLID)], SOLID),
            _difference_all,
        ))
        self.register(BuiltinFunction(
            "intersection_all",
            _make_sig("intersection_all", [ListType(SOLID)], SOLID),
            _intersection_all,
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
            """Calculate area of a region2d using shoelace formula."""
            from yapcad.geom import isgeomlist, isline, ispoint

            data = r.data
            if not isgeomlist(data):
                return float_val(0.0)

            # Sample all segments into points
            points = []
            for seg in data:
                if isline(seg):
                    p = seg[0]  # Start point of line segment
                    if ispoint(p):
                        points.append(p)
                elif ispoint(seg):
                    points.append(seg)

            if len(points) < 3:
                return float_val(0.0)

            # Shoelace formula for polygon area
            n = len(points)
            area = 0.0
            for i in range(n):
                j = (i + 1) % n
                area += points[i][0] * points[j][1]
                area -= points[j][0] * points[i][1]
            area = abs(area) / 2.0

            return float_val(area)

        def _perimeter(r: Value) -> Value:
            """Calculate perimeter of a region2d as sum of segment lengths."""
            from yapcad.geom import isgeomlist, isline, ispoint, length
            import math

            data = r.data
            if not isgeomlist(data):
                return float_val(0.0)

            total_length = 0.0
            for seg in data:
                if isline(seg):
                    # Line segment - compute length directly
                    p1, p2 = seg[0], seg[1]
                    dx = p2[0] - p1[0]
                    dy = p2[1] - p1[1]
                    total_length += math.sqrt(dx*dx + dy*dy)
                else:
                    # For arcs or other curves, use the length function
                    try:
                        total_length += length(seg)
                    except:
                        pass

            return float_val(total_length)

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

        # Phase 3: Numeric and boolean aggregation functions
        def _sum_list(lst: Value) -> Value:
            """Sum all numeric values in a list."""
            total = 0.0
            for v in lst.data:
                val = v.data if hasattr(v, 'data') else v
                total += float(val)
            return float_val(total)

        def _product_list(lst: Value) -> Value:
            """Multiply all numeric values in a list."""
            result = 1.0
            for v in lst.data:
                val = v.data if hasattr(v, 'data') else v
                result *= float(val)
            return float_val(result)

        def _any_true(lst: Value) -> Value:
            """Return True if any element is truthy."""
            for v in lst.data:
                val = v.data if hasattr(v, 'data') else v
                if val:
                    return bool_val(True)
            return bool_val(False)

        def _all_true(lst: Value) -> Value:
            """Return True if all elements are truthy."""
            for v in lst.data:
                val = v.data if hasattr(v, 'data') else v
                if not val:
                    return bool_val(False)
            return bool_val(True)

        def _min_of(lst: Value) -> Value:
            """Return minimum value in list."""
            if not lst.data:
                raise ValueError("min_of requires non-empty list")
            values = [v.data if hasattr(v, 'data') else v for v in lst.data]
            return float_val(min(values))

        def _max_of(lst: Value) -> Value:
            """Return maximum value in list."""
            if not lst.data:
                raise ValueError("max_of requires non-empty list")
            values = [v.data if hasattr(v, 'data') else v for v in lst.data]
            return float_val(max(values))

        self.register(BuiltinFunction(
            "sum",
            _make_sig("sum", [ListType(FLOAT)], FLOAT),
            _sum_list,
        ))
        self.register(BuiltinFunction(
            "product",
            _make_sig("product", [ListType(FLOAT)], FLOAT),
            _product_list,
        ))
        self.register(BuiltinFunction(
            "any_true",
            _make_sig("any_true", [ListType(BOOL)], BOOL),
            _any_true,
        ))
        self.register(BuiltinFunction(
            "all_true",
            _make_sig("all_true", [ListType(BOOL)], BOOL),
            _all_true,
        ))
        self.register(BuiltinFunction(
            "min_of",
            _make_sig("min_of", [ListType(FLOAT)], FLOAT),
            _min_of,
        ))
        self.register(BuiltinFunction(
            "max_of",
            _make_sig("max_of", [ListType(FLOAT)], FLOAT),
            _max_of,
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
