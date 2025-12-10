"""
yapCAD Domain-Specific Language (DSL) compiler.

This module provides:
- Lexer: Tokenizes DSL source code
- Parser: Builds AST from tokens
- Type checker: Validates types and constraints
- Interpreter: Executes DSL commands to generate geometry
- Transforms: AST transformations for optimization

Usage:
    from yapcad.dsl import tokenize, parse, check, Lexer, Parser, TypeChecker

    # Simple parsing
    tokens = tokenize('let x: int = 42;')

    # Or parse and type check a complete module
    source = '''
    module my_design;

    command MAKE_BOX(width: float, height: float) -> solid {
        let box: solid = box(width, height, 10.0);
        emit box;
    }
    '''
    tokens = tokenize(source)
    module = parse(tokens)
    result = check(module)
    if result.has_errors:
        for diag in result.diagnostics:
            print(diag)
"""

from .tokens import (
    Token,
    TokenType,
    SourceLocation,
    SourceSpan,
    KEYWORDS,
    is_type_token,
)

from .lexer import (
    Lexer,
    tokenize,
)

from .parser import (
    Parser,
    parse,
)

from .ast import (
    # Base
    AstNode,
    AstVisitor,
    # Types
    TypeNode,
    SimpleType,
    GenericType,
    OptionalType,
    # Expressions
    Expression,
    Literal,
    Identifier,
    BinaryOp,
    UnaryOp,
    FunctionCall,
    MethodCall,
    MemberAccess,
    IndexAccess,
    ListLiteral,
    ListComprehension,
    RangeExpr,
    IfExpr,
    MatchExpr,
    MatchArm,
    Pattern,
    LiteralPattern,
    IdentifierPattern,
    WildcardPattern,
    LambdaExpr,
    PythonExpr,
    DictLiteral,
    # Statements
    Statement,
    LetStatement,
    AssignmentStatement,
    RequireStatement,
    EmitStatement,
    ForStatement,
    WhileStatement,
    IfStatement,
    ExpressionStatement,
    ReturnStatement,
    Block,
    PythonBlock,
    ElifBranch,
    # Declarations
    Parameter,
    Command,
    UseStatement,
    ExportUseStatement,
    Module,
    # Helpers
    print_ast,
)

from .errors import (
    DslError,
    LexerError,
    ParserError,
    TypeError,
    SemanticError,
    Diagnostic,
    DiagnosticCollector,
    ErrorSeverity,
)

from .types import (
    Type,
    TypeTier,
    PrimitiveType,
    GeometricPrimitiveType,
    CurveType,
    CompoundCurveType,
    SurfaceType,
    SolidType,
    ListType,
    DictType,
    OptionalTypeWrapper,
    FunctionType,
    # Built-in type instances
    INT, FLOAT, BOOL, STRING,
    POINT, POINT2D, POINT3D,
    VECTOR, VECTOR2D, VECTOR3D,
    TRANSFORM,
    SOLID, REGION2D, SURFACE, SHELL,
    resolve_type_name,
    make_list_type,
    make_optional_type,
)

from .symbols import (
    SymbolTable,
    Symbol,
    SymbolKind,
    FunctionSignature,
    Scope,
    get_method_signature,
)

from .checker import (
    TypeChecker,
    CheckResult,
    check,
)

from .runtime import (
    # Interpreter
    Interpreter,
    ExecutionResult,
    execute,
    compile_and_run,
    # Values
    Value,
    EmitResult,
    # Context
    ExecutionContext,
    # Provenance
    Provenance,
    create_provenance,
)

from .transforms import (
    AstTransform,
    TreeTransform,
    TransformPipeline,
)

from .packaging import (
    package_from_dsl,
    PackageResult,
)

from .introspection import (
    get_api_reference,
    get_function_info,
    get_type_info,
    get_methods_for_type,
    list_functions,
    list_types,
    describe_function,
    get_api_as_json,
    get_common_pattern,
    list_common_patterns,
)

__all__ = [
    # Tokens
    'Token',
    'TokenType',
    'SourceLocation',
    'SourceSpan',
    'KEYWORDS',
    'is_type_token',

    # Lexer
    'Lexer',
    'tokenize',

    # Parser
    'Parser',
    'parse',

    # AST nodes
    'AstNode',
    'AstVisitor',
    'TypeNode',
    'SimpleType',
    'GenericType',
    'OptionalType',
    'Expression',
    'Literal',
    'Identifier',
    'BinaryOp',
    'UnaryOp',
    'FunctionCall',
    'MethodCall',
    'MemberAccess',
    'IndexAccess',
    'ListLiteral',
    'ListComprehension',
    'RangeExpr',
    'IfExpr',
    'MatchExpr',
    'MatchArm',
    'Pattern',
    'LiteralPattern',
    'IdentifierPattern',
    'WildcardPattern',
    'LambdaExpr',
    'PythonExpr',
    'DictLiteral',
    'Statement',
    'LetStatement',
    'AssignmentStatement',
    'RequireStatement',
    'EmitStatement',
    'ForStatement',
    'ExpressionStatement',
    'ReturnStatement',
    'Block',
    'PythonBlock',
    'Parameter',
    'Command',
    'UseStatement',
    'ExportUseStatement',
    'Module',
    'print_ast',

    # Errors
    'DslError',
    'LexerError',
    'ParserError',
    'TypeError',
    'SemanticError',
    'Diagnostic',
    'DiagnosticCollector',
    'ErrorSeverity',

    # Type system
    'Type',
    'TypeTier',
    'PrimitiveType',
    'GeometricPrimitiveType',
    'CurveType',
    'CompoundCurveType',
    'SurfaceType',
    'SolidType',
    'ListType',
    'DictType',
    'OptionalTypeWrapper',
    'FunctionType',
    'INT', 'FLOAT', 'BOOL', 'STRING',
    'POINT', 'POINT2D', 'POINT3D',
    'VECTOR', 'VECTOR2D', 'VECTOR3D',
    'TRANSFORM',
    'SOLID', 'REGION2D', 'SURFACE', 'SHELL',
    'resolve_type_name',
    'make_list_type',
    'make_optional_type',

    # Symbol table
    'SymbolTable',
    'Symbol',
    'SymbolKind',
    'FunctionSignature',
    'Scope',
    'get_method_signature',

    # Type checker
    'TypeChecker',
    'CheckResult',
    'check',

    # Runtime/Interpreter
    'Interpreter',
    'ExecutionResult',
    'execute',
    'compile_and_run',
    'Value',
    'EmitResult',
    'ExecutionContext',
    'Provenance',
    'create_provenance',

    # Transforms
    'AstTransform',
    'TreeTransform',
    'TransformPipeline',

    # Packaging
    'package_from_dsl',
    'PackageResult',

    # Introspection (for agentic tools)
    'get_api_reference',
    'get_function_info',
    'get_type_info',
    'get_methods_for_type',
    'list_functions',
    'list_types',
    'describe_function',
    'get_api_as_json',
    'get_common_pattern',
    'list_common_patterns',
]
