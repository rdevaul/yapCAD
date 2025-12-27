"""
Type checker for the yapCAD DSL.

Traverses the AST and validates types, collecting diagnostics
for type errors, undefined identifiers, and other semantic issues.
"""

from typing import List, Optional, Dict, Any, Union
from dataclasses import dataclass

from .ast import (
    AstNode, Module, Command, FunctionDef, Parameter, Block,
    Statement, LetStatement, VarDecl, AssignmentStatement,
    RequireStatement, AssertStatement,
    EmitStatement, ForStatement, WhileStatement, IfStatement,
    ExpressionStatement, ReturnStatement, PassStatement,
    PythonBlock, NativeBlock, NativeFunctionDecl, NativeFunction, Decorator,
    ElifBranch,
    Expression, Literal, Identifier, BinaryOp, UnaryOp,
    FunctionCall, MethodCall, MemberAccess, IndexAccess,
    ListLiteral, ListComprehension, RangeExpr, ConditionalExpr, IfExpr, MatchExpr,
    MatchArm, Pattern, LiteralPattern, IdentifierPattern, WildcardPattern,
    LambdaExpr, PythonExpr, DictLiteral,
    TypeNode, SimpleType, GenericType, OptionalType as AstOptionalType,
    UseStatement, ExportUseStatement,
)
from .tokens import TokenType, SourceSpan
from .types import (
    Type, PrimitiveType, GeometricPrimitiveType, CurveType,
    CompoundCurveType, SurfaceType, SolidType,
    ListType, DictType, OptionalTypeWrapper, FunctionType,
    ERROR, UNKNOWN, NONE,
    INT, FLOAT, BOOL, STRING,
    POINT, POINT2D, POINT3D, VECTOR, VECTOR2D, VECTOR3D, TRANSFORM,
    SOLID, REGION2D, DICT,
    resolve_type_name, make_list_type, make_optional_type,
    is_numeric, is_curve, is_geometry, common_type,
)
from .symbols import (
    SymbolTable, Symbol, SymbolKind, FunctionSignature,
    get_method_signature,
)
from .errors import (
    Diagnostic, DiagnosticCollector, ErrorSeverity,
    DslError, TypeError as DslTypeError,
)


@dataclass
class CheckResult:
    """Result of type checking a module."""
    diagnostics: List[Diagnostic]
    has_errors: bool
    has_warnings: bool
    has_python_blocks: bool  # Module contains Python blocks (requires review)


class TypeChecker:
    """
    Type checker for the yapCAD DSL.

    Traverses the AST and validates:
    - Type compatibility in assignments and function calls
    - Return type matching
    - Require expression boolean constraint
    - Emit target type matching command return type
    - Undefined identifier detection
    - Python block flagging
    """

    def __init__(self, max_errors: int = 20):
        self.symbols = SymbolTable()
        self.diagnostics = DiagnosticCollector()
        self.max_errors = max_errors
        self._current_command: Optional[Command] = None
        self._has_python_blocks = False

    def check(self, module: Module) -> CheckResult:
        """Type check a complete module."""
        self._check_module(module)

        return CheckResult(
            diagnostics=self.diagnostics.diagnostics,
            has_errors=self.diagnostics.has_errors,
            has_warnings=self.diagnostics.has_warnings,
            has_python_blocks=self._has_python_blocks
        )

    # =========================================================================
    # Module/Command Checking
    # =========================================================================

    def _check_module(self, module: Module) -> None:
        """Check a complete module."""
        # First pass: register native block exports (they provide functions to commands)
        for native_block in module.native_blocks:
            self._register_native_block(native_block)

        # Register @native decorated functions
        for native_func in getattr(module, 'native_functions', []):
            self._register_native_function(native_func)

        # Second pass: register all commands/functions
        for command in module.commands:
            self._register_command(command)

        # Third pass: check command/function bodies
        for command in module.commands:
            self._check_command(command)

    def _register_native_block(self, native_block: NativeBlock) -> None:
        """Register exported functions from a native block."""
        self._has_python_blocks = True  # Native blocks contain Python

        for func_decl in native_block.exports:
            return_type = self._resolve_type_node(func_decl.return_type)

            # Build function type from parameters
            param_types = []
            for param in func_decl.parameters:
                param_type = self._resolve_type_node(param.type_annotation)
                param_types.append(param_type)

            func_type = FunctionType(tuple(param_types), return_type)

            symbol = Symbol(
                name=func_decl.name,
                kind=SymbolKind.FUNCTION,
                type=func_type,
                span=func_decl.span
            )

            if not self.symbols.define(symbol):
                self._error(
                    f"Function '{func_decl.name}' is already defined",
                    func_decl.span,
                    "E205"
                )

        self._warning(
            "Native Python block requires manual review for type safety",
            native_block.span,
            "W212"
        )

    def _register_native_function(self, native_func: NativeFunction) -> None:
        """Register a @native decorated function."""
        self._has_python_blocks = True  # Native functions contain Python

        return_type = self._resolve_type_node(native_func.return_type)

        # Build function type from parameters
        param_types = []
        for param in native_func.parameters:
            if param.type_annotation is not None:
                param_type = self._resolve_type_node(param.type_annotation)
            else:
                param_type = UNKNOWN
            param_types.append(param_type)

        func_type = FunctionType(tuple(param_types), return_type)

        symbol = Symbol(
            name=native_func.name,
            kind=SymbolKind.FUNCTION,
            type=func_type,
            span=native_func.span
        )

        if not self.symbols.define(symbol):
            self._error(
                f"Function '{native_func.name}' is already defined",
                native_func.span,
                "E205"
            )

        self._warning(
            "Native function requires manual review for type safety",
            native_func.span,
            "W213"
        )

    def _register_command(self, command: Command) -> None:
        """Register a command/function in the symbol table."""
        # Handle optional return type (new Pythonic syntax)
        if command.return_type is not None:
            return_type = self._resolve_type_node(command.return_type)
        else:
            return_type = UNKNOWN  # Will be inferred

        # Build function type from parameters
        param_types = []
        for param in command.parameters:
            if param.type_annotation is not None:
                param_type = self._resolve_type_node(param.type_annotation)
            else:
                param_type = UNKNOWN  # Type inference
            param_types.append(param_type)

        func_type = FunctionType(tuple(param_types), return_type)

        symbol = Symbol(
            name=command.name,
            kind=SymbolKind.COMMAND,
            type=func_type,
            span=command.span
        )

        if not self.symbols.define(symbol):
            self._error(
                f"Command '{command.name}' is already defined",
                command.span,
                "E201"
            )

    def _check_command(self, command: Command) -> None:
        """Type check a command/function definition."""
        self._current_command = command
        self.symbols.push_scope(f"function {command.name}")

        # Register parameters
        for param in command.parameters:
            if param.type_annotation is not None:
                param_type = self._resolve_type_node(param.type_annotation)
            else:
                param_type = UNKNOWN
            symbol = Symbol(
                name=param.name,
                kind=SymbolKind.PARAMETER,
                type=param_type,
                span=param.span
            )
            if not self.symbols.define(symbol):
                self._error(
                    f"Duplicate parameter '{param.name}'",
                    param.span,
                    "E202"
                )

            # Check default value type
            if param.default_value is not None:
                default_type = self._check_expression(param.default_value)
                if param_type != UNKNOWN and not param_type.is_assignable_from(default_type):
                    self._error(
                        f"Default value type '{default_type}' is not assignable to "
                        f"parameter type '{param_type}'",
                        param.default_value.span,
                        "E203"
                    )

        # Check body
        self._check_block(command.body)

        self.symbols.pop_scope()
        self._current_command = None

    # =========================================================================
    # Statement Checking
    # =========================================================================

    def _check_block(self, block: Block) -> Optional[Type]:
        """Check a block of statements, returning the final expression type if any."""
        for stmt in block.statements:
            self._check_statement(stmt)

        if block.final_expression is not None:
            return self._check_expression(block.final_expression)

        return None

    def _check_statement(self, stmt: Statement) -> None:
        """Check a statement."""
        if self.diagnostics.error_count >= self.max_errors:
            return

        if isinstance(stmt, LetStatement):  # Also handles VarDecl (alias)
            self._check_let_statement(stmt)
        elif isinstance(stmt, AssignmentStatement):
            self._check_assignment_statement(stmt)
        elif isinstance(stmt, RequireStatement):  # Also handles AssertStatement (alias)
            self._check_require_statement(stmt)
        elif isinstance(stmt, EmitStatement):
            self._check_emit_statement(stmt)
        elif isinstance(stmt, ForStatement):
            self._check_for_statement(stmt)
        elif isinstance(stmt, WhileStatement):
            self._check_while_statement(stmt)
        elif isinstance(stmt, IfStatement):
            self._check_if_statement(stmt)
        elif isinstance(stmt, PassStatement):
            pass  # PassStatement has no semantic content to check
        elif isinstance(stmt, ExpressionStatement):
            self._check_expression(stmt.expression)
        elif isinstance(stmt, ReturnStatement):
            self._check_return_statement(stmt)
        elif isinstance(stmt, PythonBlock):
            self._check_python_block(stmt)
        else:
            self._warning(
                f"Unknown statement type: {type(stmt).__name__}",
                stmt.span,
                "W201"
            )

    def _check_let_statement(self, stmt: LetStatement) -> None:
        """Check a let statement."""
        init_type = self._check_expression(stmt.initializer)

        if stmt.type_annotation is not None:
            declared_type = self._resolve_type_node(stmt.type_annotation)
            if not declared_type.is_assignable_from(init_type):
                self._error(
                    f"Cannot assign '{init_type}' to variable of type '{declared_type}'",
                    stmt.initializer.span,
                    "E210"
                )
            var_type = declared_type
        else:
            var_type = init_type

        symbol = Symbol(
            name=stmt.name,
            kind=SymbolKind.VARIABLE,
            type=var_type,
            span=stmt.span,
            is_mutable=False
        )

        if not self.symbols.define(symbol):
            self._error(
                f"Variable '{stmt.name}' is already defined in this scope",
                stmt.span,
                "E211"
            )

    def _check_assignment_statement(self, stmt: AssignmentStatement) -> None:
        """Check an assignment statement."""
        target_type = self._check_expression(stmt.target)
        value_type = self._check_expression(stmt.value)

        if not target_type.is_assignable_from(value_type):
            self._error(
                f"Cannot assign '{value_type}' to target of type '{target_type}'",
                stmt.value.span,
                "E212"
            )

        # Check that target is an l-value (identifier or member/index access)
        if isinstance(stmt.target, Identifier):
            symbol = self.symbols.lookup(stmt.target.name)
            if symbol is not None and not symbol.is_mutable:
                # DSL variables are immutable by default, but reassignment is allowed
                # Update: actually, let's allow reassignment for now
                pass

    def _check_require_statement(self, stmt: RequireStatement) -> None:
        """Check a require statement."""
        cond_type = self._check_expression(stmt.condition)

        if cond_type != BOOL and cond_type != ERROR:
            self._error(
                f"Require condition must be boolean, got '{cond_type}'",
                stmt.condition.span,
                "E220"
            )

        if stmt.message is not None:
            msg_type = self._check_expression(stmt.message)
            if msg_type != STRING and msg_type != ERROR:
                self._error(
                    f"Require message must be string, got '{msg_type}'",
                    stmt.message.span,
                    "E221"
                )

    def _check_emit_statement(self, stmt: EmitStatement) -> None:
        """Check an emit statement."""
        value_type = self._check_expression(stmt.value)

        # Check that emit type matches command return type
        if self._current_command is not None and self._current_command.return_type is not None:
            return_type = self._resolve_type_node(self._current_command.return_type)
            if not return_type.is_assignable_from(value_type):
                self._error(
                    f"Emit value type '{value_type}' does not match command "
                    f"return type '{return_type}'",
                    stmt.value.span,
                    "E230"
                )

        # Check metadata if present
        # Metadata can be a DictLiteral (old syntax) or a plain dict (new kwargs syntax)
        if stmt.metadata is not None:
            if isinstance(stmt.metadata, DictLiteral):
                self._check_expression(stmt.metadata)
            elif isinstance(stmt.metadata, dict):
                # New syntax: emit value, name="x", material="y"
                # metadata is stored as a dict of key -> expression
                for key, expr in stmt.metadata.items():
                    self._check_expression(expr)

    def _check_for_statement(self, stmt: ForStatement) -> None:
        """Check a for statement."""
        iterable_type = self._check_expression(stmt.iterable)

        # Determine element type from iterable
        if isinstance(iterable_type, ListType):
            elem_type = iterable_type.element_type
        elif isinstance(iterable_type, RangeExpr):
            elem_type = INT
        else:
            # Range expressions parsed as RangeExpr, check for int range
            elem_type = INT  # Assume numeric iteration

        # Create new scope for loop body
        self.symbols.push_scope("for loop")

        symbol = Symbol(
            name=stmt.variable,
            kind=SymbolKind.VARIABLE,
            type=elem_type,
            span=stmt.span,
            is_mutable=False
        )
        self.symbols.define(symbol)

        self._check_block(stmt.body)

        self.symbols.pop_scope()

    def _check_while_statement(self, stmt: WhileStatement) -> None:
        """Check a while statement."""
        cond_type = self._check_expression(stmt.condition)

        if cond_type != BOOL and cond_type != ERROR:
            self._error(
                f"While condition must be boolean, got '{cond_type}'",
                stmt.condition.span,
                "E225"
            )

        # Create new scope for loop body
        self.symbols.push_scope("while loop")
        self._check_block(stmt.body)
        self.symbols.pop_scope()

    def _check_if_statement(self, stmt: IfStatement) -> None:
        """Check a block-level if statement."""
        cond_type = self._check_expression(stmt.condition)

        if cond_type != BOOL and cond_type != ERROR:
            self._error(
                f"If condition must be boolean, got '{cond_type}'",
                stmt.condition.span,
                "E226"
            )

        # Check then branch
        self.symbols.push_scope("if then")
        self._check_block(stmt.then_branch)
        self.symbols.pop_scope()

        # Check elif branches
        for elif_branch in stmt.elif_branches:
            elif_cond_type = self._check_expression(elif_branch.condition)
            if elif_cond_type != BOOL and elif_cond_type != ERROR:
                self._error(
                    f"Elif condition must be boolean, got '{elif_cond_type}'",
                    elif_branch.condition.span,
                    "E227"
                )
            self.symbols.push_scope("elif")
            self._check_block(elif_branch.body)
            self.symbols.pop_scope()

        # Check else branch if present
        if stmt.else_branch is not None:
            self.symbols.push_scope("else")
            self._check_block(stmt.else_branch)
            self.symbols.pop_scope()

    def _check_return_statement(self, stmt: ReturnStatement) -> None:
        """Check a return statement (for Python blocks)."""
        value_type = self._check_expression(stmt.value)
        declared_type = self._resolve_type_node(stmt.return_type)

        # The declared type is the bridge type - we trust the user's annotation
        # but flag it as requiring review
        self._has_python_blocks = True

    def _check_python_block(self, stmt: PythonBlock) -> None:
        """Check a Python block (just flag it for review)."""
        self._has_python_blocks = True
        self._warning(
            "Python block requires manual review for type safety",
            stmt.span,
            "W210"
        )

    # =========================================================================
    # Expression Checking
    # =========================================================================

    def _check_expression(self, expr: Expression) -> Type:
        """Check an expression and return its type."""
        if self.diagnostics.error_count >= self.max_errors:
            return ERROR

        if isinstance(expr, Literal):
            return self._check_literal(expr)
        elif isinstance(expr, Identifier):
            return self._check_identifier(expr)
        elif isinstance(expr, BinaryOp):
            return self._check_binary_op(expr)
        elif isinstance(expr, UnaryOp):
            return self._check_unary_op(expr)
        elif isinstance(expr, FunctionCall):
            return self._check_function_call(expr)
        elif isinstance(expr, MethodCall):
            return self._check_method_call(expr)
        elif isinstance(expr, MemberAccess):
            return self._check_member_access(expr)
        elif isinstance(expr, IndexAccess):
            return self._check_index_access(expr)
        elif isinstance(expr, ListLiteral):
            return self._check_list_literal(expr)
        elif isinstance(expr, ListComprehension):
            return self._check_list_comprehension(expr)
        elif isinstance(expr, RangeExpr):
            return self._check_range_expr(expr)
        elif isinstance(expr, ConditionalExpr):
            return self._check_conditional_expr(expr)
        elif isinstance(expr, IfExpr):
            return self._check_if_expr(expr)
        elif isinstance(expr, MatchExpr):
            return self._check_match_expr(expr)
        elif isinstance(expr, DictLiteral):
            return self._check_dict_literal(expr)
        elif isinstance(expr, LambdaExpr):
            return self._check_lambda_expr(expr)
        elif isinstance(expr, PythonExpr):
            return self._check_python_expr(expr)
        else:
            self._warning(
                f"Unknown expression type: {type(expr).__name__}",
                expr.span,
                "W220"
            )
            return ERROR

    def _check_literal(self, expr: Literal) -> Type:
        """Check a literal and return its type."""
        if expr.literal_type == TokenType.INT_LITERAL:
            return INT
        elif expr.literal_type == TokenType.FLOAT_LITERAL:
            return FLOAT
        elif expr.literal_type == TokenType.STRING_LITERAL:
            return STRING
        elif expr.literal_type == TokenType.BOOL_LITERAL:
            return BOOL
        else:
            return ERROR

    def _check_identifier(self, expr: Identifier) -> Type:
        """Check an identifier reference."""
        # First check local symbols
        symbol = self.symbols.lookup(expr.name)
        if symbol is not None:
            return symbol.type

        # Check if it's a type used as constructor (e.g., point(...))
        type_val = resolve_type_name(expr.name)
        if type_val is not None:
            # Types used as identifiers are constructors
            return UNKNOWN  # Will be resolved in function call

        # Check if it's a built-in function
        if self.symbols.is_builtin(expr.name):
            return UNKNOWN  # Function type, resolved in call

        self._error(
            f"Undefined identifier '{expr.name}'",
            expr.span,
            "E240"
        )
        return ERROR

    def _check_binary_op(self, expr: BinaryOp) -> Type:
        """Check a binary operation."""
        left_type = self._check_expression(expr.left)
        right_type = self._check_expression(expr.right)

        op = expr.operator

        # Arithmetic operators
        if op in (TokenType.PLUS, TokenType.MINUS, TokenType.STAR,
                  TokenType.SLASH, TokenType.PERCENT, TokenType.DOUBLE_SLASH,
                  TokenType.DOUBLE_STAR):
            if is_numeric(left_type) and is_numeric(right_type):
                # int op int -> int (except for division), int op float -> float
                # // (integer division) always returns int
                # ** (power) follows standard type promotion
                if op == TokenType.DOUBLE_SLASH:
                    return INT  # Integer division always returns int
                if left_type == FLOAT or right_type == FLOAT or op == TokenType.SLASH:
                    return FLOAT
                return INT
            # Vector/point arithmetic
            if isinstance(left_type, GeometricPrimitiveType) and is_numeric(right_type):
                return left_type  # scalar multiplication
            if isinstance(left_type, GeometricPrimitiveType) and isinstance(right_type, GeometricPrimitiveType):
                if op in (TokenType.PLUS, TokenType.MINUS):
                    return common_type(left_type, right_type) or left_type

            # List concatenation with +
            if op == TokenType.PLUS and isinstance(left_type, ListType) and isinstance(right_type, ListType):
                # Lists of compatible element types can be concatenated
                elem_common = common_type(left_type.element_type, right_type.element_type)
                if elem_common is not None:
                    return make_list_type(elem_common)
                # If no common type, use left element type
                return left_type

            if left_type != ERROR and right_type != ERROR:
                self._error(
                    f"Cannot apply '{op.name}' to '{left_type}' and '{right_type}'",
                    expr.span,
                    "E250"
                )
            return ERROR

        # Comparison operators
        if op in (TokenType.LT, TokenType.GT, TokenType.LE, TokenType.GE):
            if is_numeric(left_type) and is_numeric(right_type):
                return BOOL
            if left_type != ERROR and right_type != ERROR:
                self._error(
                    f"Cannot compare '{left_type}' and '{right_type}' with '{op.name}'",
                    expr.span,
                    "E251"
                )
            return ERROR

        # Equality operators
        if op in (TokenType.EQ, TokenType.NE):
            # Most types can be compared for equality
            return BOOL

        # Logical operators
        if op in (TokenType.AND, TokenType.OR):
            if left_type != BOOL and left_type != ERROR:
                self._error(
                    f"Left operand of '{op.name}' must be boolean, got '{left_type}'",
                    expr.left.span,
                    "E252"
                )
            if right_type != BOOL and right_type != ERROR:
                self._error(
                    f"Right operand of '{op.name}' must be boolean, got '{right_type}'",
                    expr.right.span,
                    "E253"
                )
            return BOOL

        return ERROR

    def _check_unary_op(self, expr: UnaryOp) -> Type:
        """Check a unary operation."""
        operand_type = self._check_expression(expr.operand)

        if expr.operator == TokenType.NOT:
            if operand_type != BOOL and operand_type != ERROR:
                self._error(
                    f"Operand of '!' must be boolean, got '{operand_type}'",
                    expr.operand.span,
                    "E260"
                )
            return BOOL

        if expr.operator == TokenType.MINUS:
            if is_numeric(operand_type):
                return operand_type
            if isinstance(operand_type, GeometricPrimitiveType):
                return operand_type  # Negating vectors/points
            if operand_type != ERROR:
                self._error(
                    f"Cannot negate '{operand_type}'",
                    expr.operand.span,
                    "E261"
                )
            return ERROR

        return ERROR

    def _check_function_call(self, expr: FunctionCall) -> Type:
        """Check a function call."""
        # Get the callee name
        callee_name = None
        if isinstance(expr.callee, Identifier):
            callee_name = expr.callee.name

        if callee_name is None:
            # Complex callee expression
            callee_type = self._check_expression(expr.callee)
            return UNKNOWN

        # Check for built-in function
        builtin = self.symbols.lookup_builtin(callee_name)
        if builtin is not None:
            return self._check_builtin_call(builtin, expr)

        # Check for type constructor
        type_val = resolve_type_name(callee_name)
        if type_val is not None:
            # Type constructor - check argument count loosely
            for arg in expr.arguments:
                self._check_expression(arg)
            for arg in expr.named_arguments.values():
                self._check_expression(arg)
            return type_val

        # Check for user-defined command or native function
        symbol = self.symbols.lookup(callee_name)
        if symbol is not None and symbol.kind in (SymbolKind.COMMAND, SymbolKind.FUNCTION):
            if isinstance(symbol.type, FunctionType):
                # Check argument count and types for native functions
                func_type = symbol.type
                for i, arg in enumerate(expr.arguments):
                    arg_type = self._check_expression(arg)
                    if i < len(func_type.param_types):
                        param_type = func_type.param_types[i]
                        if param_type != UNKNOWN and not param_type.is_assignable_from(arg_type):
                            if arg_type != ERROR:
                                self._error(
                                    f"Argument {i+1} expects '{param_type}', got '{arg_type}'",
                                    arg.span,
                                    "E274"
                                )
                return func_type.return_type
            return UNKNOWN

        self._error(
            f"Unknown function '{callee_name}'",
            expr.callee.span,
            "E270"
        )
        return ERROR

    def _check_builtin_call(self, sig: FunctionSignature, expr: FunctionCall) -> Type:
        """Check a call to a built-in function."""
        # Check argument count
        required_params = [p for p in sig.params if p[2] is None]
        if len(expr.arguments) < len(required_params):
            self._error(
                f"Function '{sig.name}' requires at least {len(required_params)} "
                f"arguments, got {len(expr.arguments)}",
                expr.span,
                "E271"
            )

        # Check argument types and collect them for type inference
        arg_types = []
        for i, arg in enumerate(expr.arguments):
            arg_type = self._check_expression(arg)
            arg_types.append(arg_type)
            if i < len(sig.params):
                param_name, param_type, _ = sig.params[i]
                if param_type != UNKNOWN and not param_type.is_assignable_from(arg_type):
                    if arg_type != ERROR:
                        self._error(
                            f"Argument '{param_name}' expects '{param_type}', "
                            f"got '{arg_type}'",
                            arg.span,
                            "E272"
                        )

        # Check named arguments
        param_names = {p[0] for p in sig.params}
        for name, arg in expr.named_arguments.items():
            arg_type = self._check_expression(arg)
            if name not in param_names:
                self._error(
                    f"Unknown parameter '{name}' for function '{sig.name}'",
                    arg.span,
                    "E273"
                )

        # Infer return type for list functions based on argument types
        return_type = sig.return_type
        if isinstance(return_type, ListType) and return_type.element_type == UNKNOWN:
            # Functions that return list<unknown> need type inference
            if sig.name in ('concat', 'reverse', 'flatten'):
                if arg_types and isinstance(arg_types[0], ListType):
                    if sig.name == 'flatten' and isinstance(arg_types[0].element_type, ListType):
                        # flatten: list<list<T>> -> list<T>
                        return_type = arg_types[0].element_type
                    else:
                        # concat, reverse: list<T> -> list<T>
                        return_type = arg_types[0]

        return return_type

    def _check_method_call(self, expr: MethodCall) -> Type:
        """Check a method call."""
        obj_type = self._check_expression(expr.object)

        # Get method signature
        method_sig = get_method_signature(obj_type, expr.method)
        if method_sig is None:
            if obj_type != ERROR:
                self._error(
                    f"Type '{obj_type}' has no method '{expr.method}'",
                    expr.span,
                    "E280"
                )
            return ERROR

        # Check arguments
        for i, arg in enumerate(expr.arguments):
            arg_type = self._check_expression(arg)
            if i < len(method_sig.params):
                param_name, param_type, _ = method_sig.params[i]
                if not param_type.is_assignable_from(arg_type):
                    if arg_type != ERROR:
                        self._error(
                            f"Method argument '{param_name}' expects '{param_type}', "
                            f"got '{arg_type}'",
                            arg.span,
                            "E281"
                        )

        return method_sig.return_type

    def _check_member_access(self, expr: MemberAccess) -> Type:
        """Check member access (e.g., point.x)."""
        obj_type = self._check_expression(expr.object)

        # Known member accesses
        if isinstance(obj_type, GeometricPrimitiveType):
            if obj_type.name.startswith("point") or obj_type.name.startswith("vector"):
                if expr.member in ("x", "y", "z"):
                    return FLOAT

        if obj_type != ERROR:
            self._error(
                f"Type '{obj_type}' has no member '{expr.member}'",
                expr.span,
                "E290"
            )
        return ERROR

    def _check_index_access(self, expr: IndexAccess) -> Type:
        """Check index access (e.g., list[0])."""
        obj_type = self._check_expression(expr.object)
        idx_type = self._check_expression(expr.index)

        if idx_type != INT and idx_type != ERROR:
            self._error(
                f"Index must be integer, got '{idx_type}'",
                expr.index.span,
                "E291"
            )

        if isinstance(obj_type, ListType):
            return obj_type.element_type

        if obj_type != ERROR:
            self._error(
                f"Cannot index into '{obj_type}'",
                expr.object.span,
                "E292"
            )
        return ERROR

    def _check_list_literal(self, expr: ListLiteral) -> Type:
        """Check a list literal."""
        if not expr.elements:
            return make_list_type(UNKNOWN)

        # Infer element type from first element
        elem_type = self._check_expression(expr.elements[0])

        for i, elem in enumerate(expr.elements[1:], start=1):
            t = self._check_expression(elem)
            if not elem_type.is_assignable_from(t):
                ct = common_type(elem_type, t)
                if ct is not None:
                    elem_type = ct
                else:
                    self._error(
                        f"List element {i} has type '{t}', expected '{elem_type}'",
                        elem.span,
                        "E293"
                    )

        return make_list_type(elem_type)

    def _check_list_comprehension(self, expr: ListComprehension) -> Type:
        """Check a list comprehension."""
        iterable_type = self._check_expression(expr.iterable)

        # Determine element type
        if isinstance(iterable_type, ListType):
            iter_elem_type = iterable_type.element_type
        else:
            iter_elem_type = INT  # Assume range

        # Create scope for comprehension variable
        self.symbols.push_scope("list comprehension")

        symbol = Symbol(
            name=expr.variable,
            kind=SymbolKind.VARIABLE,
            type=iter_elem_type,
            span=expr.span
        )
        self.symbols.define(symbol)

        # Check condition if present
        if expr.condition is not None:
            cond_type = self._check_expression(expr.condition)
            if cond_type != BOOL and cond_type != ERROR:
                self._error(
                    f"Comprehension condition must be boolean, got '{cond_type}'",
                    expr.condition.span,
                    "E294"
                )

        # Check element expression
        elem_type = self._check_expression(expr.element_expr)

        self.symbols.pop_scope()

        return make_list_type(elem_type)

    def _check_range_expr(self, expr: RangeExpr) -> Type:
        """Check a range expression (e.g., 0..10)."""
        start_type = self._check_expression(expr.start)
        end_type = self._check_expression(expr.end)

        if start_type != INT and start_type != ERROR:
            self._error(
                f"Range start must be integer, got '{start_type}'",
                expr.start.span,
                "E295"
            )
        if end_type != INT and end_type != ERROR:
            self._error(
                f"Range end must be integer, got '{end_type}'",
                expr.end.span,
                "E296"
            )

        return make_list_type(INT)

    def _check_conditional_expr(self, expr: ConditionalExpr) -> Type:
        """Check a ternary conditional expression (e.g., x if cond else y)."""
        cond_type = self._check_expression(expr.condition)
        if cond_type != BOOL and cond_type != ERROR:
            self._error(
                f"Conditional expression requires boolean condition, got '{cond_type}'",
                expr.condition.span,
                "E310"
            )

        true_type = self._check_expression(expr.true_branch)
        false_type = self._check_expression(expr.false_branch)

        # Both branches must have compatible types
        if true_type != ERROR and false_type != ERROR:
            ct = common_type(true_type, false_type)
            if ct is None and true_type != false_type:
                self._error(
                    f"Conditional branches have incompatible types: '{true_type}' and '{false_type}'",
                    expr.span,
                    "E311"
                )
                return true_type  # Return first branch type as fallback
            return ct if ct is not None else true_type

        return true_type if true_type != ERROR else false_type

    def _check_if_expr(self, expr: IfExpr) -> Type:
        """Check an if expression."""
        cond_type = self._check_expression(expr.condition)
        if cond_type != BOOL and cond_type != ERROR:
            self._error(
                f"If condition must be boolean, got '{cond_type}'",
                expr.condition.span,
                "E297"
            )

        then_type = self._check_block(expr.then_branch)

        if expr.else_branch is not None:
            if isinstance(expr.else_branch, Block):
                else_type = self._check_block(expr.else_branch)
            else:  # IfExpr (else if)
                else_type = self._check_if_expr(expr.else_branch)

            # Types should match
            if then_type and else_type:
                ct = common_type(then_type, else_type)
                if ct is None and then_type != ERROR and else_type != ERROR:
                    self._warning(
                        f"If branches have different types: '{then_type}' and '{else_type}'",
                        expr.span,
                        "W230"
                    )
                return ct or then_type

        return then_type or UNKNOWN

    def _check_match_expr(self, expr: MatchExpr) -> Type:
        """Check a match expression."""
        subject_type = self._check_expression(expr.subject)

        result_type: Optional[Type] = None

        for arm in expr.arms:
            # Check pattern (just validate it)
            self._check_pattern(arm.pattern, subject_type)

            # Check body
            arm_type = self._check_expression(arm.body)

            if result_type is None:
                result_type = arm_type
            else:
                ct = common_type(result_type, arm_type)
                if ct is not None:
                    result_type = ct

        return result_type or UNKNOWN

    def _check_pattern(self, pattern: Pattern, expected_type: Type) -> None:
        """Check a match pattern."""
        if isinstance(pattern, LiteralPattern):
            lit_type = self._check_literal(pattern.value)
            if not expected_type.is_assignable_from(lit_type):
                self._warning(
                    f"Pattern type '{lit_type}' may not match subject type '{expected_type}'",
                    pattern.span,
                    "W231"
                )
        elif isinstance(pattern, IdentifierPattern):
            # Binding pattern - introduces variable in arm scope
            # For now, just note it
            pass
        elif isinstance(pattern, WildcardPattern):
            # Always matches
            pass

    def _check_dict_literal(self, expr: DictLiteral) -> Type:
        """Check a dictionary literal."""
        for key, value in expr.entries.items():
            self._check_expression(value)
        return DICT

    def _check_lambda_expr(self, expr: LambdaExpr) -> Type:
        """Check a lambda expression."""
        # Create scope for lambda parameters
        self.symbols.push_scope("lambda")

        for param in expr.parameters:
            symbol = Symbol(
                name=param,
                kind=SymbolKind.PARAMETER,
                type=UNKNOWN,  # Type inference would happen here
                span=expr.span
            )
            self.symbols.define(symbol)

        body_type = self._check_expression(expr.body)

        self.symbols.pop_scope()

        param_types = tuple(UNKNOWN for _ in expr.parameters)
        return FunctionType(param_types, body_type)

    def _check_python_expr(self, expr: PythonExpr) -> Type:
        """Check a Python expression."""
        self._has_python_blocks = True
        self._warning(
            "Python expression requires manual review for type safety",
            expr.span,
            "W211"
        )
        return self._resolve_type_node(expr.return_type)

    # =========================================================================
    # Type Resolution
    # =========================================================================

    def _resolve_type_node(self, node: Optional[TypeNode]) -> Type:
        """Resolve an AST type node to a Type instance."""
        if node is None:
            return UNKNOWN  # Type will be inferred

        if isinstance(node, SimpleType):
            resolved = resolve_type_name(node.name)
            if resolved is None:
                self._error(
                    f"Unknown type '{node.name}'",
                    node.span,
                    "E299"
                )
                return ERROR
            return resolved

        elif isinstance(node, GenericType):
            if node.name == "list":
                if len(node.type_args) != 1:
                    self._error(
                        f"list<T> requires exactly one type argument",
                        node.span,
                        "E298"
                    )
                    return ERROR
                elem_type = self._resolve_type_node(node.type_args[0])
                return make_list_type(elem_type)
            else:
                self._error(
                    f"Unknown generic type '{node.name}'",
                    node.span,
                    "E297"
                )
                return ERROR

        elif isinstance(node, AstOptionalType):
            inner = self._resolve_type_node(node.inner)
            return make_optional_type(inner)

        return ERROR

    # =========================================================================
    # Diagnostics
    # =========================================================================

    def _error(self, message: str, span: SourceSpan, code: str) -> None:
        """Record an error diagnostic."""
        self.diagnostics.add(Diagnostic(
            code=code,
            message=message,
            severity=ErrorSeverity.ERROR,
            span=span
        ))

    def _warning(self, message: str, span: SourceSpan, code: str) -> None:
        """Record a warning diagnostic."""
        self.diagnostics.add(Diagnostic(
            code=code,
            message=message,
            severity=ErrorSeverity.WARNING,
            span=span
        ))


def check(module: Module, max_errors: int = 20) -> CheckResult:
    """
    Convenience function to type check a module.

    Args:
        module: The parsed module AST
        max_errors: Maximum errors before stopping (default 20)

    Returns:
        CheckResult with diagnostics
    """
    checker = TypeChecker(max_errors=max_errors)
    return checker.check(module)
