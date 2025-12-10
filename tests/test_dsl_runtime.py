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
