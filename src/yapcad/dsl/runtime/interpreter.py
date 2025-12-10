"""
Tree-walking interpreter for DSL execution.

Evaluates AST nodes to produce geometry and other values.
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Callable, Union

from .values import (
    Value, EmitResult, RequireFailure,
    int_val, float_val, bool_val, string_val, list_val, dict_val, none_val,
    wrap_value, unwrap_value, coerce_numeric,
)
from .context import ExecutionContext, create_context
from .builtins import call_builtin, call_method, get_builtin_registry
from .provenance import Provenance, create_provenance

from ..ast import (
    AstNode, Module, Command, FunctionDef, Parameter,
    Statement, LetStatement, VarDecl, AssignmentStatement,
    RequireStatement, AssertStatement,
    EmitStatement, ForStatement, WhileStatement, IfStatement,
    ExpressionStatement, ReturnStatement, PassStatement,
    Block, PythonBlock, NativeBlock, NativeFunction, ElifBranch,
    Expression, Literal, Identifier, BinaryOp, UnaryOp,
    FunctionCall, MethodCall, MemberAccess, IndexAccess,
    ListLiteral, ListComprehension, RangeExpr, DictLiteral,
    IfExpr, MatchExpr, MatchArm, LambdaExpr, PythonExpr,
    Pattern, LiteralPattern, IdentifierPattern, WildcardPattern,
)
from ..types import (
    Type, ListType, DictType, resolve_type_name,
    INT, FLOAT, BOOL, STRING,
)
from ..tokens import SourceSpan, TokenType


@dataclass
class ExecutionResult:
    """Result of executing a command."""
    success: bool
    emit_result: Optional[EmitResult] = None
    require_failures: List[RequireFailure] = field(default_factory=list)
    provenance: Optional[Provenance] = None
    error_message: Optional[str] = None

    @property
    def geometry(self) -> Any:
        """Get the emitted geometry data."""
        if self.emit_result:
            return self.emit_result.data
        return None

    @property
    def metadata(self) -> Dict[str, Any]:
        """Get the emit metadata."""
        if self.emit_result:
            return self.emit_result.metadata
        return {}


class Interpreter:
    """
    Tree-walking interpreter for DSL execution.

    Evaluates AST nodes by dispatching to type-specific methods.
    """

    def __init__(self, transforms: List[Callable[[Module], Module]] = None):
        """
        Initialize the interpreter.

        Args:
            transforms: Optional list of AST transformations to apply before execution
        """
        self.transforms = transforms or []
        self.native_functions: Dict[str, Callable] = {}  # Functions from native blocks

    def execute(
        self,
        module: Module,
        command_name: str,
        parameters: Dict[str, Any],
        source: str = "",
    ) -> ExecutionResult:
        """
        Execute a command from a module.

        Args:
            module: The parsed and type-checked module
            command_name: Name of the command to execute
            parameters: Parameter values (raw Python values, not wrapped)
            source: Original source code for error messages

        Returns:
            ExecutionResult with geometry and provenance
        """
        # Apply optional AST transforms
        transformed_module = module
        for transform in self.transforms:
            transformed_module = transform(transformed_module)

        # Execute native blocks to register their exported functions
        for native_block in transformed_module.native_blocks:
            try:
                self._execute_native_block(native_block)
            except Exception as e:
                return ExecutionResult(
                    success=False,
                    error_message=f"Error in native block: {e}"
                )

        # Find the command
        command = None
        for cmd in transformed_module.commands:
            if cmd.name == command_name:
                command = cmd
                break

        if command is None:
            return ExecutionResult(
                success=False,
                error_message=f"Command '{command_name}' not found in module '{module.name}'"
            )

        # Wrap parameters as Values
        wrapped_params = self._wrap_parameters(command, parameters)

        # Create execution context
        ctx = create_context(
            module_name=module.name,
            command_name=command_name,
            parameters=wrapped_params,
            source=source,
        )

        # Execute the command body
        try:
            self._execute_command(command, ctx)
        except RuntimeError as e:
            return ExecutionResult(
                success=False,
                error_message=str(e),
                require_failures=ctx.require_failures,
            )

        # Build provenance
        provenance = create_provenance(
            module_name=module.name,
            command_name=command_name,
            parameters=parameters,
            source=source,
        )

        # Check for require failures
        if ctx.require_failures:
            return ExecutionResult(
                success=False,
                emit_result=ctx.emit_result,
                require_failures=ctx.require_failures,
                provenance=provenance,
            )

        # Check for emit
        if ctx.emit_result is None:
            return ExecutionResult(
                success=False,
                error_message="Command did not emit any geometry",
                provenance=provenance,
            )

        return ExecutionResult(
            success=True,
            emit_result=ctx.emit_result,
            provenance=provenance,
        )

    def _wrap_parameters(self, command: Command, parameters: Dict[str, Any]) -> Dict[str, Value]:
        """Wrap raw parameter values as DSL Values."""
        wrapped = {}
        for param in command.parameters:
            if param.name in parameters:
                raw_value = parameters[param.name]
                # Handle optional type annotations (new Pythonic syntax)
                if param.type_annotation is not None and hasattr(param.type_annotation, 'name'):
                    param_type = resolve_type_name(param.type_annotation.name)
                else:
                    # Infer type from the raw value
                    if isinstance(raw_value, bool):
                        param_type = BOOL
                    elif isinstance(raw_value, int):
                        param_type = INT
                    elif isinstance(raw_value, float):
                        param_type = FLOAT
                    elif isinstance(raw_value, str):
                        param_type = STRING
                    else:
                        param_type = FLOAT  # Default fallback
                wrapped[param.name] = wrap_value(raw_value, param_type)
            elif param.default is not None:
                # Use default value - we'd need to evaluate it
                # For now, this is a limitation
                pass
        return wrapped

    def _execute_command(self, command: Command, ctx: ExecutionContext) -> None:
        """Execute a command's body statements."""
        for stmt in command.body.statements:
            self._execute_statement(stmt, ctx)
            if ctx.should_return:
                break

    def _execute_statement(self, stmt: Statement, ctx: ExecutionContext) -> None:
        """Execute a statement."""
        if isinstance(stmt, LetStatement):  # Also handles VarDecl (alias)
            self._execute_let(stmt, ctx)
        elif isinstance(stmt, AssignmentStatement):
            self._execute_assignment(stmt, ctx)
        elif isinstance(stmt, RequireStatement):  # Also handles AssertStatement (alias)
            self._execute_require(stmt, ctx)
        elif isinstance(stmt, EmitStatement):
            self._execute_emit(stmt, ctx)
        elif isinstance(stmt, ForStatement):
            self._execute_for(stmt, ctx)
        elif isinstance(stmt, WhileStatement):
            self._execute_while(stmt, ctx)
        elif isinstance(stmt, IfStatement):
            self._execute_if_statement(stmt, ctx)
        elif isinstance(stmt, PassStatement):
            pass  # PassStatement does nothing
        elif isinstance(stmt, ExpressionStatement):
            self._evaluate(stmt.expression, ctx)
        elif isinstance(stmt, ReturnStatement):
            self._execute_return(stmt, ctx)
        elif isinstance(stmt, Block):
            self._execute_block(stmt, ctx)
        elif isinstance(stmt, PythonBlock):
            self._execute_python_block(stmt, ctx)
        else:
            raise RuntimeError(f"Unknown statement type: {type(stmt).__name__}")

    def _execute_let(self, stmt: LetStatement, ctx: ExecutionContext) -> None:
        """Execute a let statement."""
        value = self._evaluate(stmt.initializer, ctx)
        ctx.set_variable(stmt.name, value)

    def _execute_assignment(self, stmt: AssignmentStatement, ctx: ExecutionContext) -> None:
        """Execute an assignment statement."""
        value = self._evaluate(stmt.value, ctx)
        if not ctx.update_variable(stmt.target.name, value):
            ctx.add_error(f"Cannot assign to undefined variable '{stmt.target.name}'", stmt.span)

    def _execute_require(self, stmt: RequireStatement, ctx: ExecutionContext) -> None:
        """Execute a require statement."""
        condition = self._evaluate(stmt.condition, ctx)
        if not condition.is_truthy():
            message = stmt.message if stmt.message else "Constraint violated"
            ctx.add_require_failure(message)

    def _execute_emit(self, stmt: EmitStatement, ctx: ExecutionContext) -> None:
        """Execute an emit statement."""
        value = self._evaluate(stmt.value, ctx)
        metadata = {}
        if stmt.metadata:
            # Handle both old DictLiteral syntax and new dict metadata
            if isinstance(stmt.metadata, DictLiteral):
                for key, expr in stmt.metadata.entries.items():
                    metadata[key] = self._evaluate(expr, ctx).data
            elif isinstance(stmt.metadata, dict):
                # New syntax: emit value, name="x", material="y"
                for key, expr in stmt.metadata.items():
                    metadata[key] = self._evaluate(expr, ctx).data
        ctx.set_emit(value, metadata)

    def _execute_for(self, stmt: ForStatement, ctx: ExecutionContext) -> None:
        """Execute a for loop."""
        iterable = self._evaluate(stmt.iterable, ctx)

        # Handle different iterable types
        if isinstance(iterable.type, ListType):
            items = iterable.data
        elif hasattr(iterable.data, '__iter__'):
            items = list(iterable.data)
        else:
            raise RuntimeError(f"Cannot iterate over {iterable.type}")

        # Determine element type
        if isinstance(iterable.type, ListType):
            elem_type = iterable.type.element_type
        else:
            elem_type = INT  # Default for ranges

        with ctx.new_scope("for-loop"):
            for item in items:
                ctx.set_variable(stmt.variable, wrap_value(item, elem_type))
                for body_stmt in stmt.body.statements:
                    self._execute_statement(body_stmt, ctx)
                    if ctx.should_return:
                        return

    def _execute_while(self, stmt: WhileStatement, ctx: ExecutionContext) -> None:
        """Execute a while loop."""
        with ctx.new_scope("while-loop"):
            while True:
                condition = self._evaluate(stmt.condition, ctx)
                if not condition.is_truthy():
                    break
                for body_stmt in stmt.body.statements:
                    self._execute_statement(body_stmt, ctx)
                    if ctx.should_return:
                        return

    def _execute_if_statement(self, stmt: IfStatement, ctx: ExecutionContext) -> None:
        """Execute a block-level if statement."""
        condition = self._evaluate(stmt.condition, ctx)

        if condition.is_truthy():
            # Execute then branch
            with ctx.new_scope("if-then"):
                for body_stmt in stmt.then_branch.statements:
                    self._execute_statement(body_stmt, ctx)
                    if ctx.should_return:
                        return
        else:
            # Check elif branches
            executed = False
            for elif_branch in stmt.elif_branches:
                elif_cond = self._evaluate(elif_branch.condition, ctx)
                if elif_cond.is_truthy():
                    with ctx.new_scope("elif"):
                        for body_stmt in elif_branch.body.statements:
                            self._execute_statement(body_stmt, ctx)
                            if ctx.should_return:
                                return
                    executed = True
                    break

            # Execute else branch if no elif matched
            if not executed and stmt.else_branch is not None:
                with ctx.new_scope("else"):
                    for body_stmt in stmt.else_branch.statements:
                        self._execute_statement(body_stmt, ctx)
                        if ctx.should_return:
                            return

    def _execute_return(self, stmt: ReturnStatement, ctx: ExecutionContext) -> None:
        """Execute a return statement."""
        if stmt.value:
            value = self._evaluate(stmt.value, ctx)
            ctx.signal_return(value)
        else:
            ctx.signal_return(none_val(INT))  # Return unit/none

    def _execute_block(self, block: Block, ctx: ExecutionContext) -> None:
        """Execute a block of statements."""
        with ctx.new_scope("block"):
            for stmt in block.statements:
                self._execute_statement(stmt, ctx)
                if ctx.should_return:
                    break

    def _execute_python_block(self, block: PythonBlock, ctx: ExecutionContext) -> None:
        """Execute an embedded Python block."""
        # Build the execution namespace with current variables
        namespace = {}
        scope = ctx.current_scope
        while scope:
            for name, value in scope.variables.items():
                namespace[name] = value.data
            scope = scope.parent

        # Add yapCAD imports
        try:
            exec("from yapcad.geom import *", namespace)
            exec("from yapcad.geom3d import *", namespace)
        except ImportError:
            pass

        # Execute the Python code
        try:
            exec(block.code, namespace)
        except Exception as e:
            raise RuntimeError(f"Python block error: {e}")

        # Check for return value (indicated by 'return' in code)
        # The DSL requires explicit `return <value> as <type>` syntax
        # For now, we extract any variable named '_return_value' and '_return_type'
        if '_return_value' in namespace and '_return_type' in namespace:
            return_type = resolve_type_name(namespace['_return_type'])
            ctx.set_variable(block.result_var or '_python_result',
                           wrap_value(namespace['_return_value'], return_type))

    def _execute_native_block(self, block: NativeBlock) -> None:
        """Execute a native Python block and register exported functions.

        The native block's Python code is executed in a namespace with yapCAD
        imports available. Functions declared in the exports section are then
        extracted from the namespace and stored for later use by DSL commands.
        """
        # Create namespace with yapCAD imports
        namespace = {}
        try:
            exec("from yapcad.geom import *", namespace)
            exec("from yapcad.geom3d import *", namespace)
        except ImportError:
            pass

        # Also import commonly used modules
        try:
            exec("import math", namespace)
            exec("from math import pi, sin, cos, tan, sqrt, atan2", namespace)
        except ImportError:
            pass

        # Execute the Python code
        try:
            exec(block.code, namespace)
        except Exception as e:
            raise RuntimeError(f"Native block execution error: {e}")

        # Extract exported functions from namespace
        for func_decl in block.exports:
            func_name = func_decl.name
            if func_name not in namespace:
                raise RuntimeError(
                    f"Native block exports '{func_name}' but it was not defined in the Python code"
                )
            func = namespace[func_name]
            if not callable(func):
                raise RuntimeError(
                    f"Native block exports '{func_name}' but it is not callable"
                )
            # Store the function with its declared return type
            self.native_functions[func_name] = func

    def _evaluate(self, expr: Expression, ctx: ExecutionContext) -> Value:
        """Evaluate an expression to produce a Value."""
        if isinstance(expr, Literal):
            return self._eval_literal(expr)
        elif isinstance(expr, Identifier):
            return self._eval_identifier(expr, ctx)
        elif isinstance(expr, BinaryOp):
            return self._eval_binary_op(expr, ctx)
        elif isinstance(expr, UnaryOp):
            return self._eval_unary_op(expr, ctx)
        elif isinstance(expr, FunctionCall):
            return self._eval_function_call(expr, ctx)
        elif isinstance(expr, MethodCall):
            return self._eval_method_call(expr, ctx)
        elif isinstance(expr, MemberAccess):
            return self._eval_member_access(expr, ctx)
        elif isinstance(expr, IndexAccess):
            return self._eval_index_access(expr, ctx)
        elif isinstance(expr, ListLiteral):
            return self._eval_list_literal(expr, ctx)
        elif isinstance(expr, ListComprehension):
            return self._eval_list_comprehension(expr, ctx)
        elif isinstance(expr, RangeExpr):
            return self._eval_range(expr, ctx)
        elif isinstance(expr, DictLiteral):
            return self._eval_dict_literal(expr, ctx)
        elif isinstance(expr, IfExpr):
            return self._eval_if_expr(expr, ctx)
        elif isinstance(expr, MatchExpr):
            return self._eval_match_expr(expr, ctx)
        elif isinstance(expr, LambdaExpr):
            return self._eval_lambda(expr, ctx)
        elif isinstance(expr, PythonExpr):
            return self._eval_python_expr(expr, ctx)
        else:
            raise RuntimeError(f"Unknown expression type: {type(expr).__name__}")

    def _eval_literal(self, lit: Literal) -> Value:
        """Evaluate a literal value."""
        if lit.literal_type == TokenType.INT_LITERAL:
            return int_val(int(lit.value))
        elif lit.literal_type == TokenType.FLOAT_LITERAL:
            return float_val(float(lit.value))
        elif lit.literal_type == TokenType.BOOL_LITERAL:
            # Value is already a bool from the parser
            if isinstance(lit.value, bool):
                return bool_val(lit.value)
            return bool_val(str(lit.value).lower() == "true")
        elif lit.literal_type == TokenType.STRING_LITERAL:
            return string_val(lit.value)
        else:
            raise RuntimeError(f"Unknown literal type: {lit.literal_type}")

    def _eval_identifier(self, ident: Identifier, ctx: ExecutionContext) -> Value:
        """Evaluate an identifier (variable lookup)."""
        value = ctx.get_variable(ident.name)
        if value is None:
            raise RuntimeError(f"Undefined variable: {ident.name}")
        return value

    def _eval_binary_op(self, op: BinaryOp, ctx: ExecutionContext) -> Value:
        """Evaluate a binary operation."""
        left = self._evaluate(op.left, ctx)
        right = self._evaluate(op.right, ctx)

        # Short-circuit for logical operators
        if op.operator == TokenType.AND:
            if not left.is_truthy():
                return bool_val(False)
            return bool_val(right.is_truthy())
        elif op.operator == TokenType.OR:
            if left.is_truthy():
                return bool_val(True)
            return bool_val(right.is_truthy())

        # Arithmetic operators
        if op.operator == TokenType.PLUS:
            if left.type == STRING or right.type == STRING:
                return string_val(str(left.data) + str(right.data))
            # List concatenation
            if isinstance(left.type, ListType) and isinstance(right.type, ListType):
                combined = left.data + right.data
                # Use the element type from the ListType, not from the data
                # Create Value directly since combined contains raw data (not Value objects)
                elem_type = left.type.element_type
                return Value(combined, ListType(elem_type))
            result = left.data + right.data
            if left.type == INT and right.type == INT:
                return int_val(result)
            return float_val(result)
        elif op.operator == TokenType.MINUS:
            result = left.data - right.data
            if left.type == INT and right.type == INT:
                return int_val(result)
            return float_val(result)
        elif op.operator == TokenType.STAR:
            result = left.data * right.data
            if left.type == INT and right.type == INT:
                return int_val(result)
            return float_val(result)
        elif op.operator == TokenType.SLASH:
            result = left.data / right.data
            return float_val(result)
        elif op.operator == TokenType.PERCENT:
            result = left.data % right.data
            if left.type == INT and right.type == INT:
                return int_val(result)
            return float_val(result)
        elif op.operator == TokenType.DOUBLE_SLASH:
            # Integer division
            result = left.data // right.data
            return int_val(int(result))
        elif op.operator == TokenType.DOUBLE_STAR:
            # Power operator
            result = left.data ** right.data
            if left.type == INT and right.type == INT and right.data >= 0:
                return int_val(int(result))
            return float_val(result)

        # Comparison operators
        elif op.operator == TokenType.EQ:
            return bool_val(left.data == right.data)
        elif op.operator == TokenType.NE:
            return bool_val(left.data != right.data)
        elif op.operator == TokenType.LT:
            return bool_val(left.data < right.data)
        elif op.operator == TokenType.GT:
            return bool_val(left.data > right.data)
        elif op.operator == TokenType.LE:
            return bool_val(left.data <= right.data)
        elif op.operator == TokenType.GE:
            return bool_val(left.data >= right.data)

        else:
            raise RuntimeError(f"Unknown binary operator: {op.operator}")

    def _eval_unary_op(self, op: UnaryOp, ctx: ExecutionContext) -> Value:
        """Evaluate a unary operation."""
        operand = self._evaluate(op.operand, ctx)

        if op.operator == TokenType.MINUS:
            if operand.type == INT:
                return int_val(-operand.data)
            return float_val(-operand.data)
        elif op.operator == TokenType.NOT:
            return bool_val(not operand.is_truthy())
        else:
            raise RuntimeError(f"Unknown unary operator: {op.operator}")

    def _eval_function_call(self, call: FunctionCall, ctx: ExecutionContext) -> Value:
        """Evaluate a function call."""
        # Evaluate arguments
        args = [self._evaluate(arg, ctx) for arg in call.arguments]

        # Extract function name from callee (which is an Identifier)
        if isinstance(call.callee, Identifier):
            func_name = call.callee.name
        else:
            raise RuntimeError(f"Unsupported callee type: {type(call.callee).__name__}")

        # Check for native function first
        if func_name in self.native_functions:
            return self._call_native_function(func_name, args)

        # Call the built-in function
        return call_builtin(func_name, args)

    def _call_native_function(self, func_name: str, args: List[Value]) -> Value:
        """Call a native function with the given arguments."""
        func = self.native_functions[func_name]

        # Unwrap Value objects to raw Python values
        raw_args = [arg.data for arg in args]

        try:
            result = func(*raw_args)
        except Exception as e:
            raise RuntimeError(f"Error calling native function '{func_name}': {e}")

        # Wrap the result - infer type from the result
        # For solids, we use SOLID type
        from ..types import SOLID
        try:
            from yapcad.geom3d import issolid
            if issolid(result):
                return wrap_value(result, SOLID)
        except ImportError:
            pass

        # Default to FLOAT for numeric, or wrap as-is
        if isinstance(result, (int, float)):
            return float_val(float(result))
        elif isinstance(result, bool):
            return bool_val(result)
        elif isinstance(result, str):
            return string_val(result)
        elif isinstance(result, list):
            return list_val([wrap_value(r, FLOAT) for r in result], FLOAT)
        else:
            # Default: wrap as SOLID since most native functions return geometry
            return wrap_value(result, SOLID)

    def _eval_method_call(self, call: MethodCall, ctx: ExecutionContext) -> Value:
        """Evaluate a method call."""
        receiver = self._evaluate(call.object, ctx)
        args = [self._evaluate(arg, ctx) for arg in call.arguments]

        # Get the type name for method lookup
        type_name = receiver.type.name if hasattr(receiver.type, 'name') else str(receiver.type)

        return call_method(type_name, call.method, receiver, args)

    def _eval_member_access(self, access: MemberAccess, ctx: ExecutionContext) -> Value:
        """Evaluate member access (e.g., point.x)."""
        obj = self._evaluate(access.object, ctx)

        # Handle dict access
        if isinstance(obj.type, DictType):
            return wrap_value(obj.data.get(access.member), STRING)

        # Handle named tuple fields (common in yapCAD)
        if hasattr(obj.data, access.member):
            return wrap_value(getattr(obj.data, access.member), FLOAT)

        # Handle list-based points (index by x=0, y=1, z=2)
        if access.member == 'x' and hasattr(obj.data, '__getitem__'):
            return float_val(obj.data[0])
        elif access.member == 'y' and hasattr(obj.data, '__getitem__'):
            return float_val(obj.data[1])
        elif access.member == 'z' and hasattr(obj.data, '__getitem__'):
            return float_val(obj.data[2])

        raise RuntimeError(f"Unknown member: {access.member}")

    def _eval_index_access(self, access: IndexAccess, ctx: ExecutionContext) -> Value:
        """Evaluate index access (e.g., list[0])."""
        obj = self._evaluate(access.object, ctx)
        index = self._evaluate(access.index, ctx)

        if isinstance(obj.type, ListType):
            elem_type = obj.type.element_type
            return wrap_value(obj.data[int(index.data)], elem_type)
        else:
            return wrap_value(obj.data[int(index.data)], FLOAT)

    def _eval_list_literal(self, lst: ListLiteral, ctx: ExecutionContext) -> Value:
        """Evaluate a list literal."""
        if not lst.elements:
            return list_val([], INT)  # Empty list, default to int

        values = [self._evaluate(elem, ctx) for elem in lst.elements]
        elem_type = values[0].type if values else INT
        return list_val(values, elem_type)

    def _eval_list_comprehension(self, comp: ListComprehension, ctx: ExecutionContext) -> Value:
        """Evaluate a list comprehension."""
        iterable = self._evaluate(comp.iterable, ctx)

        results = []
        with ctx.new_scope("comprehension"):
            for item in iterable.data:
                if isinstance(iterable.type, ListType):
                    elem_type = iterable.type.element_type
                else:
                    elem_type = INT
                ctx.set_variable(comp.variable, wrap_value(item, elem_type))

                # Check condition if present
                if comp.condition:
                    cond = self._evaluate(comp.condition, ctx)
                    if not cond.is_truthy():
                        continue

                # Evaluate the expression
                value = self._evaluate(comp.element_expr, ctx)
                results.append(value)

        if results:
            return list_val(results, results[0].type)
        return list_val([], INT)

    def _eval_range(self, range_expr: RangeExpr, ctx: ExecutionContext) -> Value:
        """Evaluate a range expression."""
        start = self._evaluate(range_expr.start, ctx)
        end = self._evaluate(range_expr.end, ctx)

        # Create a list of integers
        values = [int_val(i) for i in range(int(start.data), int(end.data))]
        return list_val(values, INT)

    def _eval_dict_literal(self, dct: DictLiteral, ctx: ExecutionContext) -> Value:
        """Evaluate a dict literal."""
        result = {}
        for key, value_expr in dct.entries.items():
            result[key] = self._evaluate(value_expr, ctx)
        return dict_val(result)

    def _eval_if_expr(self, if_expr: IfExpr, ctx: ExecutionContext) -> Value:
        """Evaluate an if expression."""
        condition = self._evaluate(if_expr.condition, ctx)
        if condition.is_truthy():
            return self._evaluate(if_expr.then_branch, ctx)
        elif if_expr.else_branch:
            return self._evaluate(if_expr.else_branch, ctx)
        else:
            return none_val(INT)

    def _eval_match_expr(self, match: MatchExpr, ctx: ExecutionContext) -> Value:
        """Evaluate a match expression."""
        subject = self._evaluate(match.subject, ctx)

        for arm in match.arms:
            if self._pattern_matches(arm.pattern, subject, ctx):
                return self._evaluate(arm.body, ctx)

        raise RuntimeError("No match arm matched")

    def _pattern_matches(self, pattern: Pattern, value: Value, ctx: ExecutionContext) -> bool:
        """Check if a pattern matches a value, binding variables if needed."""
        if isinstance(pattern, WildcardPattern):
            return True
        elif isinstance(pattern, LiteralPattern):
            # pattern.value is a Literal node, so we need pattern.value.value
            return value.data == pattern.value.value
        elif isinstance(pattern, IdentifierPattern):
            ctx.set_variable(pattern.name, value)
            return True
        else:
            return False

    def _eval_lambda(self, lam: LambdaExpr, ctx: ExecutionContext) -> Value:
        """Evaluate a lambda expression (returns a callable)."""
        # Capture the current scope for closure
        captured_scope = ctx.current_scope

        def lambda_impl(*args: Value) -> Value:
            with ctx.new_scope("lambda"):
                # Bind parameters
                for i, param in enumerate(lam.parameters):
                    if i < len(args):
                        ctx.set_variable(param, args[i])
                # Evaluate body
                return self._evaluate(lam.body, ctx)

        # Return as a special callable value
        from ..types import FunctionType
        return wrap_value(lambda_impl, FunctionType([], FLOAT))  # Simplified type

    def _eval_python_expr(self, expr: PythonExpr, ctx: ExecutionContext) -> Value:
        """Evaluate an inline Python expression."""
        # Build namespace from context
        namespace = {}
        scope = ctx.current_scope
        while scope:
            for name, value in scope.variables.items():
                namespace[name] = value.data
            scope = scope.parent

        # Evaluate the expression
        result = eval(expr.code, namespace)

        # Wrap result - need type annotation
        # expr.return_type is a TypeNode, get its name
        if expr.return_type and hasattr(expr.return_type, 'name'):
            result_type = resolve_type_name(expr.return_type.name)
        else:
            result_type = FLOAT
        return wrap_value(result, result_type)


# Convenience function for simple execution
def execute(
    module: Module,
    command_name: str,
    parameters: Dict[str, Any],
    source: str = "",
) -> ExecutionResult:
    """
    Execute a command from a module.

    This is a convenience wrapper around Interpreter.execute().
    """
    interpreter = Interpreter()
    return interpreter.execute(module, command_name, parameters, source)


def compile_and_run(
    source: str,
    command_name: str,
    parameters: Dict[str, Any],
) -> ExecutionResult:
    """
    High-level API to compile and run DSL source code in one call.

    This is the simplest way to execute DSL code:

        from yapcad.dsl import compile_and_run

        result = compile_and_run('''
            module my_design;
            command MAKE_BOX(w: float, h: float, d: float) -> solid {
                emit box(w, h, d);
            }
        ''', "MAKE_BOX", {"w": 10.0, "h": 20.0, "d": 5.0})

        if result.success:
            geometry = result.geometry
        else:
            print(f"Error: {result.error_message}")

    Args:
        source: DSL source code as a string
        command_name: Name of the command to execute
        parameters: Parameter values (raw Python values)

    Returns:
        ExecutionResult with geometry, provenance, and any errors
    """
    # Import lexer, parser, checker
    from ..lexer import tokenize
    from ..parser import parse
    from ..checker import check

    # Tokenize
    try:
        tokens = tokenize(source)
    except Exception as e:
        return ExecutionResult(
            success=False,
            error_message=f"Lexer error: {e}",
        )

    # Parse (pass source for native block extraction)
    try:
        module = parse(tokens, source=source)
    except Exception as e:
        return ExecutionResult(
            success=False,
            error_message=f"Parser error: {e}",
        )

    # Type check
    check_result = check(module)
    if check_result.has_errors:
        error_msgs = [str(d) for d in check_result.diagnostics if d.severity.name == 'ERROR']
        return ExecutionResult(
            success=False,
            error_message=f"Type errors: {'; '.join(error_msgs)}",
        )

    # Execute
    interpreter = Interpreter()
    return interpreter.execute(module, command_name, parameters, source)
