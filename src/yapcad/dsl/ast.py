"""
Abstract Syntax Tree (AST) node definitions for the yapCAD DSL v2 (Pythonic Syntax).

The AST represents the structure of a parsed DSL program, which can then
be type-checked and compiled/interpreted.

Changes from v1:
- FunctionDef replaces Command (uses 'def' keyword)
- VarDecl replaces LetStatement (no 'let' keyword, optional type annotation)
- AssertStatement replaces RequireStatement
- NativeFunction for @native decorated functions
- Added PassStatement, WhileStatement, ElifBranch
"""

from dataclasses import dataclass, field
from typing import Optional, List, Union, Any
from abc import ABC, abstractmethod
from .tokens import SourceSpan, TokenType


# =============================================================================
# Base Classes
# =============================================================================

@dataclass
class AstNode(ABC):
    """Base class for all AST nodes."""
    span: SourceSpan  # Source location for error reporting

    def accept(self, visitor: "AstVisitor") -> Any:
        """Accept a visitor for traversal."""
        method_name = f"visit_{self.__class__.__name__}"
        method = getattr(visitor, method_name, visitor.generic_visit)
        return method(self)


class AstVisitor(ABC):
    """Base class for AST visitors."""

    def generic_visit(self, node: AstNode) -> Any:
        """Default visit method."""
        raise NotImplementedError(f"No visitor for {node.__class__.__name__}")


# =============================================================================
# Type Nodes
# =============================================================================

@dataclass
class TypeNode(AstNode):
    """Base class for type annotations."""
    pass


@dataclass
class SimpleType(TypeNode):
    """A simple type like 'int', 'float', 'solid', etc."""
    name: str  # The type name


@dataclass
class GenericType(TypeNode):
    """A generic type like 'list[point3d]' or 'dict[str, int]'."""
    name: str  # 'list' or 'dict'
    type_args: List[TypeNode]  # Type arguments


@dataclass
class OptionalType(TypeNode):
    """An optional type, e.g., 'point3d?'."""
    inner: TypeNode


# =============================================================================
# Expression Nodes
# =============================================================================

@dataclass
class Expression(AstNode):
    """Base class for all expressions."""
    pass


@dataclass
class Literal(Expression):
    """A literal value (int, float, string, bool)."""
    value: Union[int, float, str, bool]
    literal_type: TokenType  # INT_LITERAL, FLOAT_LITERAL, STRING_LITERAL, BOOL_LITERAL


@dataclass
class Identifier(Expression):
    """A variable or function name reference."""
    name: str


@dataclass
class BinaryOp(Expression):
    """A binary operation (e.g., a + b, x and y)."""
    left: Expression
    operator: TokenType  # Includes AND, OR for logical operators
    right: Expression


@dataclass
class UnaryOp(Expression):
    """A unary operation (e.g., not x, -n)."""
    operator: TokenType  # Includes NOT for logical negation
    operand: Expression


@dataclass
class FunctionCall(Expression):
    """A function or constructor call (e.g., point(1, 2, 3))."""
    callee: Expression  # Identifier or member access
    arguments: List[Expression]
    named_arguments: dict[str, Expression] = field(default_factory=dict)


@dataclass
class MethodCall(Expression):
    """A method call (e.g., curve.at(0.5))."""
    object: Expression
    method: str
    arguments: List[Expression]
    named_arguments: dict[str, Expression] = field(default_factory=dict)


@dataclass
class MemberAccess(Expression):
    """Member access (e.g., point.x)."""
    object: Expression
    member: str


@dataclass
class IndexAccess(Expression):
    """Index access (e.g., list[0])."""
    object: Expression
    index: Expression


@dataclass
class ListLiteral(Expression):
    """A list literal (e.g., [1, 2, 3])."""
    elements: List[Expression]


@dataclass
class ComprehensionClause(AstNode):
    """A single for clause in a list comprehension.

    Represents: for variable in iterable [if condition1] [if condition2] ...

    Multiple if conditions are allowed per clause, all must be true.
    """
    variable: str
    iterable: Expression
    conditions: List[Expression] = field(default_factory=list)


@dataclass
class ListComprehension(Expression):
    """A list comprehension with one or more for clauses.

    Single clause (original):
        [f(x) for x in items if cond]

    Multiple clauses (nested):
        [f(x, y) for x in xs for y in ys]
        [f(x, y) for x in xs for y in ys if x < y]
        [f(x, y) for x in xs if x > 0 for y in ys if y < 10]

    The clauses are evaluated left-to-right as nested loops.
    """
    element_expr: Expression
    clauses: List[ComprehensionClause] = field(default_factory=list)

    # Legacy compatibility properties for single-clause comprehensions
    @property
    def variable(self) -> str:
        """Get variable name (for single-clause compatibility)."""
        if self.clauses:
            return self.clauses[0].variable
        return ""

    @property
    def iterable(self) -> Expression:
        """Get iterable (for single-clause compatibility)."""
        if self.clauses:
            return self.clauses[0].iterable
        # Return a dummy - shouldn't happen in valid AST
        from .tokens import SourceLocation
        return Identifier(
            span=self.span,
            name=""
        )

    @property
    def condition(self) -> Optional[Expression]:
        """Get first condition (for single-clause compatibility)."""
        if self.clauses and self.clauses[0].conditions:
            return self.clauses[0].conditions[0]
        return None

    @classmethod
    def from_single(cls, span: SourceSpan, element_expr: Expression,
                    variable: str, iterable: Expression,
                    condition: Optional[Expression] = None) -> "ListComprehension":
        """Create a single-clause comprehension (backward compatibility)."""
        clause = ComprehensionClause(
            span=span,
            variable=variable,
            iterable=iterable,
            conditions=[condition] if condition else []
        )
        return cls(span=span, element_expr=element_expr, clauses=[clause])


@dataclass
class RangeExpr(Expression):
    """A range expression (e.g., 0..10 or range(10))."""
    start: Expression
    end: Expression
    step: Optional[Expression] = None  # For range(start, end, step)


@dataclass
class ConditionalExpr(Expression):
    """A ternary conditional expression (e.g., x if condition else y).

    Unlike IfExpr (which uses blocks), this is for inline expressions:
        value: float = a if condition else b

    Both branches must be expressions and the else is required.
    """
    condition: Expression
    true_branch: Expression
    false_branch: Expression


@dataclass
class IfExpr(Expression):
    """An if-else expression (returns a value)."""
    condition: Expression
    then_branch: "Block"
    elif_branches: List["ElifBranch"] = field(default_factory=list)
    else_branch: Optional["Block"] = None


@dataclass
class ElifBranch(AstNode):
    """An elif branch in an if expression/statement."""
    condition: Expression
    body: "Block"


@dataclass
class MatchExpr(Expression):
    """A match expression."""
    subject: Expression
    arms: List["MatchArm"]


@dataclass
class MatchArm(AstNode):
    """A single arm of a match expression."""
    pattern: "Pattern"
    body: Expression


@dataclass
class Pattern(AstNode):
    """Base class for match patterns."""
    pass


@dataclass
class LiteralPattern(Pattern):
    """A literal pattern (e.g., 'match x { case 42: ... }')."""
    value: Literal


@dataclass
class IdentifierPattern(Pattern):
    """A binding pattern (e.g., 'match x { case n: ... }')."""
    name: str


@dataclass
class WildcardPattern(Pattern):
    """The wildcard pattern '_'."""
    pass


@dataclass
class LambdaExpr(Expression):
    """A lambda/anonymous function (e.g., (x) => x * 2)."""
    parameters: List[str]
    body: Expression


@dataclass
class PythonExpr(Expression):
    """A python block that returns a value (legacy support)."""
    code: str
    return_type: "TypeNode"  # The 'as <type>' annotation


# =============================================================================
# Statement Nodes
# =============================================================================

@dataclass
class Statement(AstNode):
    """Base class for all statements."""
    pass


@dataclass
class VarDecl(Statement):
    """A variable declaration (Pythonic style, no 'let' keyword).

    Syntax options:
        x = 42              # Type inferred
        x: int = 42         # Explicit type
        x: int              # Declaration without initialization (rare)
    """
    name: str
    type_annotation: Optional[TypeNode]
    initializer: Optional[Expression]


# Keep LetStatement as alias for backward compatibility during transition
LetStatement = VarDecl


@dataclass
class AssignmentStatement(Statement):
    """An assignment to an existing variable (e.g., x = 5)."""
    target: Expression  # Identifier or member/index access
    value: Expression


@dataclass
class AssertStatement(Statement):
    """An assert statement (e.g., assert x > 0, "x must be positive")."""
    condition: Expression
    message: Optional[Expression] = None  # String literal message


# Keep RequireStatement as alias for backward compatibility
RequireStatement = AssertStatement


@dataclass
class EmitStatement(Statement):
    """An emit statement with optional metadata kwargs.

    Syntax:
        emit gear                              # Simple emit
        emit gear, name="spur", material="steel"  # With metadata
    """
    value: Expression
    metadata: dict[str, Expression] = field(default_factory=dict)


@dataclass
class ReturnStatement(Statement):
    """A return statement."""
    value: Optional[Expression] = None


@dataclass
class PassStatement(Statement):
    """A pass statement (placeholder for empty blocks)."""
    pass


@dataclass
class DictLiteral(Expression):
    """A dictionary literal (e.g., {"key": value, ...})."""
    entries: dict[str, Expression]


@dataclass
class ForStatement(Statement):
    """A for loop (e.g., for i in range(n):)."""
    variable: str
    iterable: Expression
    body: "Block"


@dataclass
class WhileStatement(Statement):
    """A while loop (e.g., while condition:)."""
    condition: Expression
    body: "Block"


@dataclass
class IfStatement(Statement):
    """An if statement (doesn't return a value).

    Syntax:
        if condition:
            ...
        elif condition:
            ...
        else:
            ...
    """
    condition: Expression
    then_branch: "Block"
    elif_branches: List[ElifBranch] = field(default_factory=list)
    else_branch: Optional["Block"] = None


@dataclass
class ExpressionStatement(Statement):
    """An expression used as a statement."""
    expression: Expression


@dataclass
class Block(AstNode):
    """A block of statements (indented block).

    In Pythonic syntax, blocks are delimited by indentation (INDENT/DEDENT)
    rather than braces.
    """
    statements: List[Statement]
    # The last statement may be an expression that produces a value
    final_expression: Optional[Expression] = None


@dataclass
class PythonBlock(Statement):
    """An inline Python block (legacy support)."""
    code: str


# =============================================================================
# Decorators and Function Definitions
# =============================================================================

@dataclass
class Decorator(AstNode):
    """A decorator (e.g., @native)."""
    name: str
    arguments: List[Expression] = field(default_factory=list)


@dataclass
class Parameter(AstNode):
    """A function parameter."""
    name: str
    type_annotation: Optional[TypeNode] = None  # Optional for type inference
    default_value: Optional[Expression] = None


@dataclass
class FunctionDef(AstNode):
    """A function definition (uses 'def' keyword).

    Syntax:
        def function_name(param1: type1, param2: type2) -> return_type:
            ...

    Replaces the old 'command' syntax.
    """
    name: str
    parameters: List[Parameter]
    return_type: Optional[TypeNode]  # Optional for type inference
    body: Block
    decorators: List[Decorator] = field(default_factory=list)


# Keep Command as alias for backward compatibility
Command = FunctionDef


@dataclass
class NativeFunction(AstNode):
    """A native function (Python code with DSL type signature).

    Syntax:
        @native
        def function_name(param1: type1, param2: type2) -> return_type:
            '''Python code here'''
            ...

    The function body contains Python code that is executed directly.
    """
    name: str
    parameters: List[Parameter]
    return_type: TypeNode
    python_code: str  # The Python code in the function body


# =============================================================================
# Legacy Native Block Support (for backward compatibility)
# =============================================================================

@dataclass
class NativeFunctionDecl(AstNode):
    """A function declaration in a native block's exports section (legacy).

    Represents: fn name(param1: type1, param2: type2) -> return_type;
    """
    name: str
    parameters: List["Parameter"]
    return_type: TypeNode


@dataclass
class NativeBlock(AstNode):
    """A native Python block with exported function declarations (legacy).

    Syntax:
        native python {
            # Python code here
        } exports {
            fn func_name(param: type) -> return_type;
        }

    The Python code is executed to define functions, which are then
    made available to the DSL with the declared type signatures.

    Note: This is the legacy syntax. New code should use @native decorator.
    """
    code: str
    exports: List[NativeFunctionDecl]


# =============================================================================
# Module-Level Declarations
# =============================================================================

@dataclass
class UseStatement(AstNode):
    """A use/import statement (e.g., use yapcad.stdlib.transforms).

    Syntax:
        use module.path
        use module.path as alias
        use module.path.{item1, item2}  # (future: selective imports)
    """
    module_path: List[str]  # ['yapcad', 'stdlib', 'transforms']
    alias: Optional[str] = None  # 'as' alias


@dataclass
class ExportStatement(AstNode):
    """An export statement (e.g., export function_name)."""
    name: str


@dataclass
class ExportUseStatement(AstNode):
    """An export use statement (e.g., export use other.module)."""
    module_path: List[str]


@dataclass
class Module(AstNode):
    """A complete DSL module.

    Syntax:
        module module_name

        use other.module

        @native
        def native_func(...):
            ...

        def my_function(...):
            ...
    """
    name: Optional[str]  # module name, or None for scripts
    uses: List[Union[UseStatement, ExportUseStatement]] = field(default_factory=list)
    native_blocks: List[NativeBlock] = field(default_factory=list)  # Legacy support
    native_functions: List[NativeFunction] = field(default_factory=list)  # New style
    functions: List[FunctionDef] = field(default_factory=list)  # All non-native functions
    exports: List[ExportStatement] = field(default_factory=list)

    @property
    def commands(self) -> List[FunctionDef]:
        """Backward compatibility alias for functions."""
        return self.functions


# =============================================================================
# Visitor Helpers
# =============================================================================

class PrintVisitor(AstVisitor):
    """Debug visitor that prints the AST structure."""

    def __init__(self, indent: int = 0):
        self.indent = indent

    def _print(self, text: str) -> None:
        print("  " * self.indent + text)

    def generic_visit(self, node: AstNode) -> None:
        self._print(f"{node.__class__.__name__}")
        for name, value in node.__dict__.items():
            if name == "span":
                continue
            if isinstance(value, AstNode):
                self._print(f"  {name}:")
                PrintVisitor(self.indent + 2).generic_visit(value)
            elif isinstance(value, list):
                self._print(f"  {name}: [")
                for item in value:
                    if isinstance(item, AstNode):
                        PrintVisitor(self.indent + 2).generic_visit(item)
                    else:
                        self._print(f"    {item!r}")
                self._print("  ]")
            else:
                self._print(f"  {name}: {value!r}")


def print_ast(node: AstNode) -> None:
    """Print an AST node for debugging."""
    PrintVisitor().generic_visit(node)
