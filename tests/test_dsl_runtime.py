"""
Tests for DSL runtime (interpreter, values, context, provenance).
"""

import pytest
import math
import textwrap

from yapcad.dsl import (
    tokenize, parse, check,
    Interpreter, ExecutionResult, execute,
    Value, Provenance, create_provenance,
)
from yapcad.dsl.runtime import (
    int_val, float_val, bool_val, string_val, list_val,
    wrap_value, unwrap_value, coerce_numeric,
    ExecutionContext, create_context, Scope,
    get_builtin_registry, call_builtin,
    compute_source_signature, verify_source_signature,
)
from yapcad.dsl.transforms import (
    AstTransform, TreeTransform, TransformPipeline, IdentityTransform,
)
from yapcad.dsl.types import INT, FLOAT, BOOL, STRING, ListType


# --- Value Tests ---

class TestValues:
    """Test runtime value wrappers."""

    def test_int_value(self):
        """Test integer value creation."""
        v = int_val(42)
        assert v.data == 42
        assert v.type == INT

    def test_float_value(self):
        """Test float value creation."""
        v = float_val(3.14)
        assert v.data == 3.14
        assert v.type == FLOAT

    def test_bool_value(self):
        """Test boolean value creation."""
        v_true = bool_val(True)
        v_false = bool_val(False)
        assert v_true.data is True
        assert v_false.data is False
        assert v_true.type == BOOL

    def test_string_value(self):
        """Test string value creation."""
        v = string_val("hello")
        assert v.data == "hello"
        assert v.type == STRING

    def test_list_value(self):
        """Test list value creation."""
        items = [int_val(1), int_val(2), int_val(3)]
        v = list_val(items, INT)
        assert v.data == [1, 2, 3]
        assert isinstance(v.type, ListType)
        assert v.type.element_type == INT

    def test_is_truthy(self):
        """Test truthiness checking."""
        assert bool_val(True).is_truthy() is True
        assert bool_val(False).is_truthy() is False
        assert int_val(0).is_truthy() is True  # Numbers are always truthy
        assert int_val(42).is_truthy() is True
        assert list_val([int_val(1)], INT).is_truthy() is True
        assert list_val([], INT).is_truthy() is False  # Empty list is falsy

    def test_coerce_numeric(self):
        """Test int to float coercion."""
        v_int = int_val(42)
        v_float = coerce_numeric(v_int)
        assert v_float.type == FLOAT
        assert v_float.data == 42.0

        # Float stays float
        v_f = float_val(3.14)
        v_f2 = coerce_numeric(v_f)
        assert v_f2.type == FLOAT


# --- Context Tests ---

class TestContext:
    """Test execution context and scoping."""

    def test_variable_set_get(self):
        """Test setting and getting variables."""
        ctx = ExecutionContext()
        ctx.set_variable("x", int_val(10))
        v = ctx.get_variable("x")
        assert v.data == 10

    def test_undefined_variable(self):
        """Test getting undefined variable returns None."""
        ctx = ExecutionContext()
        assert ctx.get_variable("undefined") is None

    def test_nested_scope(self):
        """Test nested scopes."""
        ctx = ExecutionContext()
        ctx.set_variable("outer", int_val(1))

        with ctx.new_scope("inner"):
            ctx.set_variable("inner_var", int_val(2))
            # Can see outer variable
            assert ctx.get_variable("outer").data == 1
            # Can see inner variable
            assert ctx.get_variable("inner_var").data == 2

        # Outside, inner variable is gone
        assert ctx.get_variable("inner_var") is None
        # Outer still exists
        assert ctx.get_variable("outer").data == 1

    def test_variable_shadowing(self):
        """Test that inner scope can shadow outer variables."""
        ctx = ExecutionContext()
        ctx.set_variable("x", int_val(1))

        with ctx.new_scope("inner"):
            ctx.set_variable("x", int_val(2))
            assert ctx.get_variable("x").data == 2

        # Outside, original value
        assert ctx.get_variable("x").data == 1

    def test_variable_update(self):
        """Test updating existing variables through scope chain."""
        ctx = ExecutionContext()
        ctx.set_variable("x", int_val(1))

        with ctx.new_scope("inner"):
            # Update outer variable
            ctx.update_variable("x", int_val(2))
            assert ctx.get_variable("x").data == 2

        # Change persisted to outer scope
        assert ctx.get_variable("x").data == 2

    def test_require_failure(self):
        """Test recording require failures."""
        ctx = ExecutionContext()
        assert not ctx.has_errors

        ctx.add_require_failure("Value must be positive")
        assert ctx.has_errors
        assert len(ctx.require_failures) == 1
        assert "positive" in ctx.require_failures[0].message

    def test_create_context_with_params(self):
        """Test creating context with parameters."""
        params = {
            "width": float_val(10.0),
            "count": int_val(5),
        }
        ctx = create_context("test_module", "TEST_CMD", params)

        assert ctx.module_name == "test_module"
        assert ctx.command_name == "TEST_CMD"
        assert ctx.get_variable("width").data == 10.0
        assert ctx.get_variable("count").data == 5


# --- Builtin Tests ---

class TestBuiltins:
    """Test built-in function registry."""

    def test_math_functions(self):
        """Test math built-in functions."""
        # sin
        result = call_builtin("sin", [float_val(0.0)])
        assert abs(result.data) < 1e-10

        result = call_builtin("sin", [float_val(math.pi / 2)])
        assert abs(result.data - 1.0) < 1e-10

        # cos
        result = call_builtin("cos", [float_val(0.0)])
        assert abs(result.data - 1.0) < 1e-10

        # sqrt
        result = call_builtin("sqrt", [float_val(4.0)])
        assert result.data == 2.0

    def test_abs_function(self):
        """Test abs function."""
        result = call_builtin("abs", [float_val(-5.0)])
        assert result.data == 5.0

    def test_floor_ceil(self):
        """Test floor and ceil functions."""
        result = call_builtin("floor", [float_val(3.7)])
        assert result.data == 3
        assert result.type == INT

        result = call_builtin("ceil", [float_val(3.2)])
        assert result.data == 4

    def test_radians_degrees(self):
        """Test angle conversion functions."""
        result = call_builtin("radians", [float_val(180.0)])
        assert abs(result.data - math.pi) < 1e-10

        result = call_builtin("degrees", [float_val(math.pi)])
        assert abs(result.data - 180.0) < 1e-10

    def test_range_function(self):
        """Test range function."""
        result = call_builtin("range", [int_val(5)])
        assert result.data == [0, 1, 2, 3, 4]

        result = call_builtin("range", [int_val(2), int_val(5)])
        assert result.data == [2, 3, 4]

    def test_len_function(self):
        """Test len function."""
        lst = list_val([int_val(1), int_val(2), int_val(3)], INT)
        result = call_builtin("len", [lst])
        assert result.data == 3

    def test_unknown_function_raises(self):
        """Test that unknown function raises error."""
        with pytest.raises(RuntimeError, match="Unknown built-in function"):
            call_builtin("nonexistent_function", [])


# --- Provenance Tests ---

class TestProvenance:
    """Test provenance tracking."""

    def test_create_provenance(self):
        """Test creating provenance record."""
        prov = create_provenance(
            module_name="gear_module",
            command_name="MAKE_GEAR",
            parameters={"teeth": 24, "module_mm": 2.0},
            source="command MAKE_GEAR...",
            version="1.0.0",
        )

        assert prov.module_name == "gear_module"
        assert prov.command_name == "MAKE_GEAR"
        assert prov.parameters["teeth"] == 24
        assert prov.version == "1.0.0"
        assert prov.source_signature.startswith("sha256:")

    def test_provenance_to_dict(self):
        """Test provenance serialization."""
        prov = create_provenance(
            module_name="test",
            command_name="CMD",
            parameters={"x": 1},
        )
        d = prov.to_dict()

        assert "invocation" in d
        assert d["invocation"]["module"] == "test"
        assert d["invocation"]["command"] == "CMD"

    def test_source_signature_verification(self):
        """Test source signature verification."""
        source = "let x: int = 42;"
        sig = compute_source_signature(source)

        assert verify_source_signature(source, sig)
        assert not verify_source_signature("let y: int = 0;", sig)

    def test_signature_normalizes_whitespace(self):
        """Test that signature normalizes whitespace."""
        source1 = "let x: int = 42;"
        source2 = "let x: int = 42;  "  # trailing space
        source3 = "  let x: int = 42;"  # leading space

        sig1 = compute_source_signature(source1)
        sig2 = compute_source_signature(source2)
        sig3 = compute_source_signature(source3)

        # All should be equal after normalization
        assert sig1 == sig2 == sig3


# --- Transform Tests ---

class TestTransforms:
    """Test AST transformation framework."""

    def test_identity_transform(self):
        """Test identity transform returns module unchanged."""
        source = """
        module test;
        command TEST() -> int {
            let x: int = 42;
            emit x;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)

        transform = IdentityTransform()
        result = transform.transform(module)

        assert result.name == module.name
        assert len(result.commands) == len(module.commands)

    def test_transform_pipeline(self):
        """Test chaining multiple transforms."""
        source = """
        module test;
        command TEST() -> int {
            emit 42;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)

        pipeline = TransformPipeline()
        pipeline.add(IdentityTransform())
        pipeline.add(IdentityTransform())

        result = pipeline.apply(module)
        assert result.name == module.name


# --- Integration Tests ---

class TestInterpreterBasic:
    """Basic interpreter integration tests (without geometry)."""

    def test_interpreter_simple_arithmetic(self):
        """Test interpreter evaluates arithmetic expressions."""
        source = """
        module test;
        command COMPUTE(a: int, b: int) -> int {
            let sum: int = a + b;
            let product: int = a * b;
            emit sum;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        check_result = check(module)
        assert not check_result.has_errors

        # Execute - note: without actual geometry ops, emit just captures the value
        interpreter = Interpreter()
        result = interpreter.execute(
            module,
            "COMPUTE",
            {"a": 5, "b": 3},
            source=source,
        )

        # For now, execution will fail because we need actual yapCAD for emit
        # This is expected - the interpreter framework is in place
        # Full integration requires connecting to yapCAD geometry

    def test_interpreter_command_not_found(self):
        """Test interpreter returns error for unknown command."""
        source = """
        module test;
        command EXISTING() -> int {
            emit 1;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)

        interpreter = Interpreter()
        result = interpreter.execute(
            module,
            "NONEXISTENT",
            {},
        )

        assert not result.success
        assert "not found" in result.error_message


# --- End-to-End Integration Tests ---

class TestDSLIntegration:
    """Full end-to-end integration tests with yapCAD geometry."""

    def test_box_creation(self):
        """Test DSL can create a box solid with correct volume."""
        from yapcad.geom3d import issolid, volumeof

        source = """
        module test_design;

        command MAKE_BOX(width: float, height: float, depth: float) -> solid {
            let b: solid = box(width, height, depth);
            emit b;
        }
        """

        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(
            module,
            "MAKE_BOX",
            {"width": 10.0, "height": 20.0, "depth": 5.0},
            source=source,
        )

        assert exec_result.success
        assert issolid(exec_result.geometry)
        assert abs(volumeof(exec_result.geometry) - 1000.0) < 0.01

    def test_arithmetic_in_dsl(self):
        """Test arithmetic expressions in DSL work correctly."""
        source = """
        module calc;

        command CALC(a: float, b: float) -> float {
            let result: float = a * 2.0 + b;
            emit result;
        }
        """

        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(
            module,
            "CALC",
            {"a": 5.0, "b": 3.0},
            source=source,
        )

        assert exec_result.success
        assert exec_result.geometry == 13.0  # 5*2 + 3

    def test_provenance_tracking(self):
        """Test that provenance is correctly captured."""
        source = """
        module provenance_test;

        command MY_CMD(x: float) -> float {
            emit x;
        }
        """

        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)

        interpreter = Interpreter()
        exec_result = interpreter.execute(
            module,
            "MY_CMD",
            {"x": 42.0},
            source=source,
        )

        assert exec_result.success
        assert exec_result.provenance is not None
        assert exec_result.provenance.module_name == "provenance_test"
        assert exec_result.provenance.command_name == "MY_CMD"
        assert exec_result.provenance.parameters["x"] == 42.0


class TestDSLPackaging:
    """Test DSL-to-Package integration."""

    def test_package_from_dsl(self, tmp_path):
        """Test creating a package from DSL source."""
        from yapcad.dsl import package_from_dsl
        from yapcad.geom3d import issolid, volumeof

        source = textwrap.dedent("""
        module pkg_test;

        command MAKE_BOX(w: float, h: float, d: float) -> solid {
            let result: solid = box(w, h, d);
            emit result;
        }
        """)

        pkg_dir = tmp_path / "test_pkg"
        result = package_from_dsl(
            source,
            "MAKE_BOX",
            {"w": 10.0, "h": 5.0, "d": 2.0},
            pkg_dir,
            name="test_box",
            version="1.0.0",
            description="Test box package",
        )

        assert result.success
        assert result.manifest is not None
        assert result.manifest.data["name"] == "test_box"
        assert result.manifest.data["version"] == "1.0.0"

        # Check generator metadata includes DSL info
        generator = result.manifest.data.get("generator", {})
        assert generator.get("tool") == "yapCAD-DSL"
        assert "dsl" in generator
        assert generator["dsl"]["module"] == "pkg_test"
        assert generator["dsl"]["command"] == "MAKE_BOX"

        # Check DSL source attachment
        attachments = result.manifest.data.get("attachments", [])
        dsl_attachment = next((a for a in attachments if a["id"] == "dsl-source"), None)
        assert dsl_attachment is not None
        assert dsl_attachment["format"] == "dsl"

        # Verify source file exists
        source_file = result.manifest.root / "attachments" / "source.dsl"
        assert source_file.exists()

    def test_package_from_dsl_error_handling(self, tmp_path):
        """Test package_from_dsl handles DSL errors gracefully."""
        from yapcad.dsl import package_from_dsl

        source = """
        module bad;
        command BAD() -> int {
            let x: int = "not an int";  // Type error
            emit x;
        }
        """

        pkg_dir = tmp_path / "bad_pkg"
        result = package_from_dsl(
            source,
            "BAD",
            {},
            pkg_dir,
            name="bad_pkg",
            version="1.0.0",
        )

        assert not result.success
        assert "Type error" in result.error_message or "DSL execution failed" in result.error_message


class TestNew2DFeatures:
    """Test new 2D geometry DSL features."""

    def test_ellipse_constructor(self):
        """Test ellipse curve construction."""
        source = """
        module test_ellipse;

        command MAKE_ELLIPSE() -> ellipse {
            let center: point = point(0.0, 0.0);
            let e: ellipse = ellipse(center, 50.0, 30.0);
            emit e;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(module, "MAKE_ELLIPSE", {}, source=source)
        assert exec_result.success
        # Ellipse is ['ellipse', center, metadata]
        assert exec_result.geometry[0] == 'ellipse'

    def test_catmullrom_constructor(self):
        """Test Catmull-Rom spline construction."""
        source = """
        module test_spline;

        command MAKE_SPLINE() -> catmullrom {
            let p1: point = point(0.0, 0.0);
            let p2: point = point(10.0, 5.0);
            let p3: point = point(20.0, -3.0);
            let p4: point = point(30.0, 0.0);
            let pts: list<point> = [p1, p2, p3, p4];
            let spline: catmullrom = catmullrom(pts);
            emit spline;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(module, "MAKE_SPLINE", {}, source=source)
        assert exec_result.success
        # Catmullrom is ['catmullrom', points, metadata]
        assert exec_result.geometry[0] == 'catmullrom'
        assert len(exec_result.geometry[1]) == 4  # 4 control points

    def test_curve_sampling(self):
        """Test curve sampling functions."""
        source = """
        module test_sample;

        command SAMPLE_LINE() -> point {
            let start: point = point(0.0, 0.0);
            let end: point = point(10.0, 0.0);
            let l: line_segment = line(start, end);
            let mid: point = sample_curve(l, 0.5);
            emit mid;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(module, "SAMPLE_LINE", {}, source=source)
        assert exec_result.success
        # Midpoint should be at (5, 0)
        assert abs(exec_result.geometry[0] - 5.0) < 0.01
        assert abs(exec_result.geometry[1] - 0.0) < 0.01

    def test_curve_length(self):
        """Test curve length computation."""
        source = """
        module test_length;

        command LINE_LENGTH() -> float {
            let start: point = point(0.0, 0.0);
            let end: point = point(10.0, 0.0);
            let l: line_segment = line(start, end);
            let len: float = curve_length(l);
            emit len;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(module, "LINE_LENGTH", {}, source=source)
        assert exec_result.success
        assert abs(exec_result.geometry - 10.0) < 0.01

    def test_polygon_from_points(self):
        """Test polygon region construction from points."""
        source = """
        module test_polygon;

        command MAKE_TRIANGLE() -> region2d {
            let p1: point = point(0.0, 0.0);
            let p2: point = point(10.0, 0.0);
            let p3: point = point(5.0, 8.66);
            let pts: list<point> = [p1, p2, p3];
            let tri: region2d = polygon(pts);
            emit tri;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(module, "MAKE_TRIANGLE", {}, source=source)
        assert exec_result.success
        # Triangle region should have 3 line segments
        assert len(exec_result.geometry) == 3

    def test_disk_region(self):
        """Test disk (filled circle) region construction."""
        source = """
        module test_disk;

        command MAKE_DISK() -> region2d {
            let center: point = point(0.0, 0.0);
            let d: region2d = disk(center, 10.0, 32);
            emit d;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(module, "MAKE_DISK", {}, source=source)
        assert exec_result.success
        # Disk with 32 segments should have 32 line segments
        assert len(exec_result.geometry) == 32

    def test_2d_boolean_difference(self):
        """Test 2D boolean difference operation."""
        source = """
        module test_bool2d;

        command MAKE_PLATE_WITH_HOLE() -> region2d {
            let plate: region2d = rectangle(100.0, 50.0);
            let hole: region2d = disk(point(0.0, 0.0), 15.0, 32);
            let result: region2d = difference2d(plate, hole);
            emit result;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(module, "MAKE_PLATE_WITH_HOLE", {}, source=source)
        assert exec_result.success
        # Result should be a non-empty geometry list
        assert len(exec_result.geometry) > 0

    def test_2d_boolean_union(self):
        """Test 2D boolean union operation."""
        source = """
        module test_union2d;

        command UNION_RECTS() -> region2d {
            let r1: region2d = rectangle(20.0, 10.0, point(-5.0, 0.0));
            let r2: region2d = rectangle(20.0, 10.0, point(5.0, 0.0));
            let result: region2d = union2d(r1, r2);
            emit result;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(module, "UNION_RECTS", {}, source=source)
        assert exec_result.success
        # Result should be a non-empty geometry list
        assert len(exec_result.geometry) > 0

    def test_make_path2d_and_close(self):
        """Test path2d construction and closing."""
        source = """
        module test_path2d;

        command MAKE_PATH() -> region2d {
            let l1: line_segment = line(point(0.0, 0.0), point(10.0, 0.0));
            let l2: line_segment = line(point(10.0, 0.0), point(10.0, 10.0));
            let l3: line_segment = line(point(10.0, 10.0), point(0.0, 10.0));
            let path: path2d = make_path2d([l1, l2, l3]);
            let region: region2d = close_path(path);
            emit region;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(module, "MAKE_PATH", {}, source=source)
        assert exec_result.success
        # Closed path should have 4 segments (3 + closing)
        assert len(exec_result.geometry) == 4

    def test_nurbs_constructor(self):
        """Test NURBS curve construction."""
        source = """
        module test_nurbs;

        command MAKE_NURBS() -> nurbs {
            let p1: point = point(0.0, 0.0);
            let p2: point = point(5.0, 10.0);
            let p3: point = point(15.0, 10.0);
            let p4: point = point(20.0, 0.0);
            let pts: list<point> = [p1, p2, p3, p4];
            let curve: nurbs = nurbs(pts);
            emit curve;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(module, "MAKE_NURBS", {}, source=source)
        assert exec_result.success
        # NURBS is ['nurbs', points, metadata]
        assert exec_result.geometry[0] == 'nurbs'
        assert len(exec_result.geometry[1]) == 4

    def test_spline_sampling(self):
        """Test sampling a catmullrom spline for multiple points."""
        source = """
        module test_spline_sample;

        command SAMPLE_SPLINE() -> list<point> {
            let p1: point = point(0.0, 0.0);
            let p2: point = point(10.0, 5.0);
            let p3: point = point(20.0, -3.0);
            let p4: point = point(30.0, 0.0);
            let pts: list<point> = [p1, p2, p3, p4];
            let spline: catmullrom = catmullrom(pts);
            let samples: list<point> = sample_curve_n(spline, 10);
            emit samples;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(module, "SAMPLE_SPLINE", {}, source=source)
        assert exec_result.success
        # Should return 10 sample points
        assert len(exec_result.geometry) == 10
        # First point should be close to p1 (0,0)
        assert abs(exec_result.geometry[0][0] - 0.0) < 1.0
        # Last point should be close to p4 (30,0)
        assert abs(exec_result.geometry[-1][0] - 30.0) < 1.0

    def test_ellipse_sampling(self):
        """Test sampling an ellipse curve."""
        source = """
        module test_ellipse_sample;

        command SAMPLE_ELLIPSE() -> point {
            let center: point = point(0.0, 0.0);
            let e: ellipse = ellipse(center, 50.0, 30.0);
            let top: point = sample_curve(e, 0.25);
            emit top;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(module, "SAMPLE_ELLIPSE", {}, source=source)
        assert exec_result.success
        # At t=0.25 (quarter turn), should be at top of ellipse: (0, semi_minor)
        # For ellipse with semi_major=50, semi_minor=30, top is at (0, 30)
        assert abs(exec_result.geometry[0]) < 1.0  # x near 0
        assert abs(exec_result.geometry[1] - 30.0) < 1.0  # y near 30

    def test_chained_difference2d(self):
        """Test chained 2D difference operations preserve all holes."""
        source = """
        module test_chained_diff;

        command PLATE_WITH_THREE_HOLES() -> region2d {
            let plate: region2d = rectangle(100.0, 60.0);
            let hole1: region2d = disk(point(0.0, 0.0), 10.0, 16);
            let hole2: region2d = disk(point(-30.0, 0.0), 5.0, 12);
            let hole3: region2d = disk(point(30.0, 0.0), 5.0, 12);

            let step1: region2d = difference2d(plate, hole1);
            let step2: region2d = difference2d(step1, hole2);
            let result: region2d = difference2d(step2, hole3);
            emit result;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(module, "PLATE_WITH_THREE_HOLES", {}, source=source)
        assert exec_result.success
        # Result should be a nested list: [outer, hole1, hole2, hole3]
        # 4 geometry lists total
        assert isinstance(exec_result.geometry, list)
        assert len(exec_result.geometry) == 4  # outer + 3 holes
        # Each element should be a geomlist (list of lines)
        from yapcad.geom import isline
        for poly in exec_result.geometry:
            assert isinstance(poly, list)
            assert len(poly) > 0
            assert isline(poly[0])


class TestDXFExport:
    """Test DXF export functionality."""

    def test_dxf_export_basic(self):
        """Test basic DXF export of 2D geometry."""
        import tempfile
        import os

        source = """
        module test_dxf;

        command MAKE_RECT() -> region2d {
            let r: region2d = rectangle(50.0, 30.0);
            emit r;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(module, "MAKE_RECT", {}, source=source)
        assert exec_result.success

        # Export to DXF
        from yapcad.ezdxf_exporter import write_dxf, is_2d_geometry

        assert is_2d_geometry(exec_result.geometry)

        with tempfile.TemporaryDirectory() as tmpdir:
            dxf_path = os.path.join(tmpdir, "test.dxf")
            success = write_dxf(exec_result.geometry, dxf_path)
            assert success
            assert os.path.exists(dxf_path)

            # Verify file content
            with open(dxf_path, 'r') as f:
                content = f.read()
                assert 'LINE' in content  # Should contain LINE entities
                assert 'SOLID' not in content.split('\n')  # No SOLID entities

    def test_dxf_export_with_holes(self):
        """Test DXF export of region with holes."""
        import tempfile
        import os

        source = """
        module test_dxf_holes;

        command PLATE_WITH_HOLE() -> region2d {
            let plate: region2d = rectangle(80.0, 50.0);
            let hole: region2d = disk(point(0.0, 0.0), 15.0, 24);
            let result: region2d = difference2d(plate, hole);
            emit result;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(module, "PLATE_WITH_HOLE", {}, source=source)
        assert exec_result.success

        from yapcad.ezdxf_exporter import write_dxf

        with tempfile.TemporaryDirectory() as tmpdir:
            dxf_path = os.path.join(tmpdir, "plate_hole.dxf")
            success = write_dxf(exec_result.geometry, dxf_path)
            assert success

            # Count LINE entities - should have 4 (rect) + 24 (hole) = 28
            import ezdxf
            doc = ezdxf.readfile(dxf_path)
            msp = doc.modelspace()
            line_count = len(list(msp.query('LINE')))
            assert line_count == 28  # 4 + 24


# --- Conditional Expression Tests ---

class TestConditionalExpressions:
    """Test ternary conditional expression support."""

    def test_conditional_true_branch(self):
        """Test conditional expression selects true branch."""
        source = """
        module test;
        command COND_TEST(flag: bool) -> float {
            let result: float = 10.0 if flag else 20.0;
            emit result;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(
            module,
            "COND_TEST",
            {"flag": True},
            source=source,
        )

        assert exec_result.success
        assert exec_result.geometry == 10.0

    def test_conditional_false_branch(self):
        """Test conditional expression selects false branch."""
        source = """
        module test;
        command COND_TEST(flag: bool) -> float {
            let result: float = 10.0 if flag else 20.0;
            emit result;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)

        interpreter = Interpreter()
        exec_result = interpreter.execute(
            module,
            "COND_TEST",
            {"flag": False},
            source=source,
        )

        assert exec_result.success
        assert exec_result.geometry == 20.0

    def test_conditional_with_comparison(self):
        """Test conditional expression with comparison condition."""
        source = """
        module test;
        command THRESHOLD(value: float) -> float {
            let result: float = 100.0 if value > 50.0 else 0.0;
            emit result;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()

        # Above threshold
        exec_result = interpreter.execute(
            module,
            "THRESHOLD",
            {"value": 75.0},
            source=source,
        )
        assert exec_result.success
        assert exec_result.geometry == 100.0

        # Below threshold
        exec_result = interpreter.execute(
            module,
            "THRESHOLD",
            {"value": 25.0},
            source=source,
        )
        assert exec_result.success
        assert exec_result.geometry == 0.0

    def test_conditional_nested(self):
        """Test nested (chained) conditional expressions."""
        source = """
        module test;
        command GRADE(score: float) -> string {
            let grade: string = "A" if score >= 90.0 else ("B" if score >= 80.0 else "C");
            emit grade;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()

        # A grade
        exec_result = interpreter.execute(
            module,
            "GRADE",
            {"score": 95.0},
            source=source,
        )
        assert exec_result.success
        assert exec_result.geometry == "A"

        # B grade
        exec_result = interpreter.execute(
            module,
            "GRADE",
            {"score": 85.0},
            source=source,
        )
        assert exec_result.success
        assert exec_result.geometry == "B"

        # C grade
        exec_result = interpreter.execute(
            module,
            "GRADE",
            {"score": 70.0},
            source=source,
        )
        assert exec_result.success
        assert exec_result.geometry == "C"

    def test_conditional_in_expression(self):
        """Test conditional expression used inline in arithmetic."""
        source = """
        module test;
        command CALC(value: float, use_metric: bool) -> float {
            let factor: float = 1.0 if use_metric else 25.4;
            let result: float = value * factor;
            emit result;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()

        # Metric
        exec_result = interpreter.execute(
            module,
            "CALC",
            {"value": 10.0, "use_metric": True},
            source=source,
        )
        assert exec_result.success
        assert abs(exec_result.geometry - 10.0) < 0.001

        # Imperial
        exec_result = interpreter.execute(
            module,
            "CALC",
            {"value": 10.0, "use_metric": False},
            source=source,
        )
        assert exec_result.success
        assert abs(exec_result.geometry - 254.0) < 0.001

    def test_conditional_with_geometry(self):
        """Test conditional expression selecting geometry."""
        from yapcad.geom3d import issolid, volumeof

        source = """
        module test;
        command MAKE_SHAPE(use_cube: bool, size: float) -> solid {
            let shape: solid = box(size, size, size) if use_cube else cylinder(size/2.0, size);
            emit shape;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()

        # Cube: volume = size^3 = 1000
        exec_result = interpreter.execute(
            module,
            "MAKE_SHAPE",
            {"use_cube": True, "size": 10.0},
            source=source,
        )
        assert exec_result.success
        assert issolid(exec_result.geometry)
        assert abs(volumeof(exec_result.geometry) - 1000.0) < 1.0

        # Cylinder: volume = pi * r^2 * h = pi * 25 * 10 ~ 785.4
        exec_result = interpreter.execute(
            module,
            "MAKE_SHAPE",
            {"use_cube": False, "size": 10.0},
            source=source,
        )
        assert exec_result.success
        assert issolid(exec_result.geometry)
        expected_vol = math.pi * 25.0 * 10.0
        # Mesh approximation can have ~1% error
        assert abs(volumeof(exec_result.geometry) - expected_vol) < expected_vol * 0.01

    def test_conditional_type_error_condition(self):
        """Test that non-boolean condition raises type error."""
        source = """
        module test;
        command BAD_COND(x: int) -> float {
            let result: float = 1.0 if x else 2.0;
            emit result;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        # Should have an error about condition type
        assert result.has_errors
        error_msgs = [d.message for d in result.diagnostics]
        assert any("boolean" in msg.lower() for msg in error_msgs)

    def test_conditional_type_error_branches(self):
        """Test that incompatible branch types raise type error."""
        source = """
        module test;
        command BAD_TYPES(flag: bool) -> float {
            let result: float = 1.0 if flag else "nope";
            emit result;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        # Should have an error about incompatible types
        assert result.has_errors


class TestFunctionalCombinators:
    """Test Phase 3 functional combinator functions."""

    def test_sum_function(self):
        """Test sum() aggregation function."""
        source = """
        module test;
        command CALC_SUM() -> float {
            let values: list<float> = [1.0, 2.0, 3.0, 4.0, 5.0];
            let total: float = sum(values);
            emit total;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(module, "CALC_SUM", {}, source=source)
        assert exec_result.success
        assert abs(exec_result.geometry - 15.0) < 0.001

    def test_product_function(self):
        """Test product() aggregation function."""
        source = """
        module test;
        command CALC_PRODUCT() -> float {
            let values: list<float> = [1.0, 2.0, 3.0, 4.0, 5.0];
            let result: float = product(values);
            emit result;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(module, "CALC_PRODUCT", {}, source=source)
        assert exec_result.success
        assert abs(exec_result.geometry - 120.0) < 0.001

    def test_min_of_function(self):
        """Test min_of() aggregation function."""
        source = """
        module test;
        command FIND_MIN() -> float {
            let values: list<float> = [5.0, 2.0, 8.0, 1.0, 9.0];
            let minimum: float = min_of(values);
            emit minimum;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(module, "FIND_MIN", {}, source=source)
        assert exec_result.success
        assert abs(exec_result.geometry - 1.0) < 0.001

    def test_max_of_function(self):
        """Test max_of() aggregation function."""
        source = """
        module test;
        command FIND_MAX() -> float {
            let values: list<float> = [5.0, 2.0, 8.0, 1.0, 9.0];
            let maximum: float = max_of(values);
            emit maximum;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(module, "FIND_MAX", {}, source=source)
        assert exec_result.success
        assert abs(exec_result.geometry - 9.0) < 0.001

    def test_any_true_function(self):
        """Test any_true() boolean aggregation."""
        source = """
        module test;
        command CHECK_ANY() -> bool {
            let values: list<bool> = [false, false, true, false];
            let result: bool = any_true(values);
            emit result;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(module, "CHECK_ANY", {}, source=source)
        assert exec_result.success
        assert exec_result.geometry is True

    def test_any_true_all_false(self):
        """Test any_true() with all false values."""
        source = """
        module test;
        command CHECK_ANY_FALSE() -> bool {
            let values: list<bool> = [false, false, false];
            let result: bool = any_true(values);
            emit result;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(module, "CHECK_ANY_FALSE", {}, source=source)
        assert exec_result.success
        assert exec_result.geometry is False

    def test_all_true_function(self):
        """Test all_true() boolean aggregation."""
        source = """
        module test;
        command CHECK_ALL() -> bool {
            let values: list<bool> = [true, true, true];
            let result: bool = all_true(values);
            emit result;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(module, "CHECK_ALL", {}, source=source)
        assert exec_result.success
        assert exec_result.geometry is True

    def test_all_true_with_false(self):
        """Test all_true() with one false value."""
        source = """
        module test;
        command CHECK_ALL_FALSE() -> bool {
            let values: list<bool> = [true, false, true];
            let result: bool = all_true(values);
            emit result;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(module, "CHECK_ALL_FALSE", {}, source=source)
        assert exec_result.success
        assert exec_result.geometry is False

    def test_union_all_function(self):
        """Test union_all() solid aggregation."""
        from yapcad.geom3d import issolid, volumeof

        source = """
        module test;
        command MAKE_UNION() -> solid {
            let boxes: list<solid> = [
                box(10.0, 10.0, 10.0),
                translate(box(10.0, 10.0, 10.0), 20.0, 0.0, 0.0),
                translate(box(10.0, 10.0, 10.0), 40.0, 0.0, 0.0)
            ];
            let result: solid = union_all(boxes);
            emit result;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(module, "MAKE_UNION", {}, source=source)
        assert exec_result.success
        assert issolid(exec_result.geometry)
        # 3 boxes of 1000 each = 3000 total volume
        assert abs(volumeof(exec_result.geometry) - 3000.0) < 10.0

    def test_difference_all_function(self):
        """Test difference_all() solid aggregation."""
        from yapcad.geom3d import issolid, volumeof

        source = """
        module test;
        command MAKE_DIFF() -> solid {
            let plate: solid = box(100.0, 100.0, 10.0);
            let holes: list<solid> = [
                translate(cylinder(5.0, 12.0), -20.0, 0.0, -1.0),
                translate(cylinder(5.0, 12.0), 0.0, 0.0, -1.0),
                translate(cylinder(5.0, 12.0), 20.0, 0.0, -1.0)
            ];
            let result: solid = difference_all(plate, holes);
            emit result;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(module, "MAKE_DIFF", {}, source=source)
        assert exec_result.success
        assert issolid(exec_result.geometry)
        # 100*100*10 = 100000, minus 3 cylinders of pi*25*10 each
        plate_vol = 100000.0
        hole_vol = math.pi * 25.0 * 10.0 * 3
        expected = plate_vol - hole_vol
        # Allow some tolerance for mesh approximation
        assert abs(volumeof(exec_result.geometry) - expected) < expected * 0.02

    def test_intersection_all_function(self):
        """Test intersection_all() solid aggregation."""
        from yapcad.geom3d import issolid, volumeof

        source = """
        module test;
        command MAKE_INTERSECT() -> solid {
            let solids: list<solid> = [
                box(20.0, 20.0, 20.0),
                translate(box(20.0, 20.0, 20.0), 5.0, 0.0, 0.0)
            ];
            let result: solid = intersection_all(solids);
            emit result;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(module, "MAKE_INTERSECT", {}, source=source)
        assert exec_result.success
        assert issolid(exec_result.geometry)
        # Two 20x20x20 boxes offset by 5 in X: overlap is 15x20x20 = 6000
        expected = 15.0 * 20.0 * 20.0
        assert abs(volumeof(exec_result.geometry) - expected) < expected * 0.02

    def test_union2d_all_function(self):
        """Test union2d_all() 2D region aggregation."""
        source = """
        module test;
        command MAKE_2D_UNION() -> region2d {
            let shapes: list<region2d> = [
                rectangle(10.0, 10.0),
                disk(point(20.0, 0.0), 5.0, 32)
            ];
            let result: region2d = union2d_all(shapes);
            emit result;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(module, "MAKE_2D_UNION", {}, source=source)
        assert exec_result.success
        # Result should be a valid 2D geometry
        assert exec_result.geometry is not None

    def test_difference2d_all_function(self):
        """Test difference2d_all() 2D region aggregation."""
        source = """
        module test;
        command MAKE_2D_DIFF() -> region2d {
            let plate: region2d = rectangle(100.0, 50.0);
            let holes: list<region2d> = [
                disk(point(-30.0, 0.0), 8.0, 32),
                disk(point(0.0, 0.0), 8.0, 32),
                disk(point(30.0, 0.0), 8.0, 32)
            ];
            let result: region2d = difference2d_all(plate, holes);
            emit result;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(module, "MAKE_2D_DIFF", {}, source=source)
        assert exec_result.success
        assert exec_result.geometry is not None

    def test_combinators_with_comprehension(self):
        """Test functional combinators combined with list comprehension."""
        source = """
        module test;
        command CALC_SQUARES_SUM() -> float {
            let values: list<float> = [1.0, 2.0, 3.0, 4.0];
            let squares: list<float> = [x * x for x in values];
            let total: float = sum(squares);
            emit total;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(module, "CALC_SQUARES_SUM", {}, source=source)
        assert exec_result.success
        # 1 + 4 + 9 + 16 = 30
        assert abs(exec_result.geometry - 30.0) < 0.001

    def test_difference_all_with_comprehension(self):
        """Test difference_all with list comprehension for hole array."""
        from yapcad.geom3d import issolid, volumeof

        source = """
        module test;
        command MAKE_HOLE_PATTERN() -> solid {
            let plate: solid = cylinder(50.0, 10.0);
            let angles: list<float> = [0.0, 60.0, 120.0, 180.0, 240.0, 300.0];
            let holes: list<solid> = [
                translate(cylinder(5.0, 12.0), 30.0 * cos(radians(a)), 30.0 * sin(radians(a)), -1.0)
                for a in angles
            ];
            let result: solid = difference_all(plate, holes);
            emit result;
        }
        """
        tokens = tokenize(textwrap.dedent(source))
        module = parse(tokens)
        result = check(module)
        assert not result.has_errors

        interpreter = Interpreter()
        exec_result = interpreter.execute(module, "MAKE_HOLE_PATTERN", {}, source=source)
        assert exec_result.success
        assert issolid(exec_result.geometry)
        # Should be plate volume minus 6 holes
        plate_vol = math.pi * 50.0 * 50.0 * 10.0
        holes_vol = 6 * math.pi * 25.0 * 10.0
        expected = plate_vol - holes_vol
        assert abs(volumeof(exec_result.geometry) - expected) < expected * 0.05
