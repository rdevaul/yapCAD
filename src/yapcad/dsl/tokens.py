"""
Token types for the yapCAD DSL lexer.

yapCAD DSL v2 - Pythonic Syntax

Token type categories follow error code ranges from the roadmap:
- E0xx: Lexer errors
- E1xx: Parser errors
- E2xx: Type errors
- E3xx: Semantic errors
"""

from enum import Enum, auto
from dataclasses import dataclass
from typing import Any, Optional


class TokenType(Enum):
    """All token types recognized by the DSL lexer."""

    # --- Literals ---
    INT_LITERAL = auto()        # 42, 0xff, 0b1010
    FLOAT_LITERAL = auto()      # 3.14, 1e-9, 2.5E+10
    STRING_LITERAL = auto()     # "hello", """multi\nline"""
    BOOL_LITERAL = auto()       # True, False

    # --- Identifiers ---
    IDENTIFIER = auto()         # user-defined names

    # --- Keywords (Pythonic) ---
    MODULE = auto()             # module
    USE = auto()                # use (import-like)
    DEF = auto()                # def (function definition)
    RETURN = auto()             # return
    EMIT = auto()               # emit (with metadata)
    IF = auto()                 # if
    ELIF = auto()               # elif
    ELSE = auto()               # else
    FOR = auto()                # for
    IN = auto()                 # in
    WHILE = auto()              # while (future)
    ASSERT = auto()             # assert
    PASS = auto()               # pass
    AS = auto()                 # as
    MATCH = auto()              # match (pattern matching)
    EXPORT = auto()             # export
    NATIVE = auto()             # native (for @native decorator)

    # --- Deprecated keywords (kept for error messages) ---
    COMMAND = auto()            # command (use 'def' instead)
    LET = auto()                # let (no longer needed)
    REQUIRE = auto()            # require (use 'assert' instead)
    WITH = auto()               # with (use emit kwargs instead)
    PYTHON = auto()             # python (use @native instead)
    FN = auto()                 # fn (use def instead)
    EXPORTS = auto()            # exports (use @native instead)

    # --- Closure keywords ---
    CLOSE = auto()              # close
    CLOSE_C0 = auto()           # closeC0
    CLOSE_C1 = auto()           # closeC1

    # --- Type keywords: Tier 1 - Primitives ---
    TYPE_INT = auto()           # int
    TYPE_FLOAT = auto()         # float
    TYPE_STRING = auto()        # string / str
    TYPE_BOOL = auto()          # bool
    TYPE_POINT = auto()         # point
    TYPE_POINT2D = auto()       # point2d
    TYPE_POINT3D = auto()       # point3d
    TYPE_VECTOR = auto()        # vector
    TYPE_VECTOR2D = auto()      # vector2d
    TYPE_VECTOR3D = auto()      # vector3d
    TYPE_TRANSFORM = auto()     # transform

    # --- Type keywords: Tier 2 - Curves ---
    TYPE_LINE_SEGMENT = auto()  # line_segment
    TYPE_ARC = auto()           # arc
    TYPE_CIRCLE = auto()        # circle
    TYPE_ELLIPSE = auto()       # ellipse
    TYPE_PARABOLA = auto()      # parabola
    TYPE_HYPERBOLA = auto()     # hyperbola
    TYPE_CATMULLROM = auto()    # catmullrom
    TYPE_NURBS = auto()         # nurbs
    TYPE_BEZIER = auto()        # bezier

    # --- Type keywords: Tier 3 - Compound curves ---
    TYPE_PATH2D = auto()        # path2d
    TYPE_PATH3D = auto()        # path3d
    TYPE_PROFILE2D = auto()     # profile2d
    TYPE_REGION2D = auto()      # region2d
    TYPE_LOOP3D = auto()        # loop3d

    # --- Type keywords: Tier 4 - Surfaces ---
    TYPE_SURFACE = auto()       # surface
    TYPE_SHELL = auto()         # shell

    # --- Type keywords: Tier 5 - Solids ---
    TYPE_SOLID = auto()         # solid

    # --- Generic type keywords ---
    TYPE_LIST = auto()          # list
    TYPE_DICT = auto()          # dict

    # --- Arithmetic operators ---
    PLUS = auto()               # +
    MINUS = auto()              # -
    STAR = auto()               # *
    SLASH = auto()              # /
    DOUBLE_SLASH = auto()       # // (integer division)
    PERCENT = auto()            # %
    DOUBLE_STAR = auto()        # ** (power)

    # --- Comparison operators ---
    LT = auto()                 # <
    GT = auto()                 # >
    LE = auto()                 # <=
    GE = auto()                 # >=
    EQ = auto()                 # ==
    NE = auto()                 # !=

    # --- Logical operators (keyword-based, Python style) ---
    AND = auto()                # and
    OR = auto()                 # or
    NOT = auto()                # not

    # --- Assignment ---
    ASSIGN = auto()             # =
    PLUS_ASSIGN = auto()        # += (future)
    MINUS_ASSIGN = auto()       # -= (future)

    # --- Delimiters ---
    LBRACE = auto()             # { (for dict literals only)
    RBRACE = auto()             # } (for dict literals only)
    LPAREN = auto()             # (
    RPAREN = auto()             # )
    LBRACKET = auto()           # [
    RBRACKET = auto()           # ]
    COLON = auto()              # :
    SEMICOLON = auto()          # ; (deprecated, kept for error messages)
    COMMA = auto()              # ,
    DOT = auto()                # .
    ARROW = auto()              # ->
    DOUBLE_ARROW = auto()       # => (for lambdas)
    RANGE = auto()              # .. (range literal syntax)
    QUESTION = auto()           # ? (optional type)
    UNDERSCORE = auto()         # _ (wildcard in match)
    AT = auto()                 # @ (decorator)
    HASH = auto()               # # (comment start - usually skipped)

    # --- Indentation tokens (Python-style blocks) ---
    INDENT = auto()             # Increase in indentation level
    DEDENT = auto()             # Decrease in indentation level
    NEWLINE = auto()            # Significant newline (end of statement)

    # --- Special ---
    EOF = auto()                # end of file

    # --- Native/Python block content ---
    NATIVE_BLOCK = auto()       # Content of @native decorated function


@dataclass(frozen=True)
class SourceLocation:
    """Represents a position in source code."""
    line: int           # 1-indexed line number
    column: int         # 1-indexed column number
    offset: int         # 0-indexed character offset from start
    filename: Optional[str] = None

    def __str__(self) -> str:
        if self.filename:
            return f"{self.filename}:{self.line}:{self.column}"
        return f"{self.line}:{self.column}"


@dataclass(frozen=True)
class SourceSpan:
    """Represents a range in source code."""
    start: SourceLocation
    end: SourceLocation

    def __str__(self) -> str:
        if self.start.filename:
            return f"{self.start.filename}:{self.start.line}:{self.start.column}-{self.end.line}:{self.end.column}"
        return f"{self.start.line}:{self.start.column}-{self.end.line}:{self.end.column}"


@dataclass(frozen=True)
class Token:
    """A single token from the lexer."""
    type: TokenType
    value: Any              # The actual value (int, float, str, etc.)
    lexeme: str             # The original source text
    span: SourceSpan        # Location in source

    def __str__(self) -> str:
        if self.type in (TokenType.INT_LITERAL, TokenType.FLOAT_LITERAL,
                         TokenType.STRING_LITERAL, TokenType.IDENTIFIER):
            return f"{self.type.name}({self.value!r})"
        return self.type.name


# Keyword mapping - maps string to token type
KEYWORDS: dict[str, TokenType] = {
    # Core keywords (Pythonic)
    "module": TokenType.MODULE,
    "use": TokenType.USE,
    "def": TokenType.DEF,
    "return": TokenType.RETURN,
    "emit": TokenType.EMIT,
    "if": TokenType.IF,
    "elif": TokenType.ELIF,
    "else": TokenType.ELSE,
    "for": TokenType.FOR,
    "in": TokenType.IN,
    "while": TokenType.WHILE,
    "assert": TokenType.ASSERT,
    "pass": TokenType.PASS,
    "as": TokenType.AS,
    "match": TokenType.MATCH,
    "export": TokenType.EXPORT,
    "native": TokenType.NATIVE,

    # Logical operators (Python style)
    "and": TokenType.AND,
    "or": TokenType.OR,
    "not": TokenType.NOT,

    # Deprecated keywords (kept for helpful error messages)
    "command": TokenType.COMMAND,
    "let": TokenType.LET,
    "require": TokenType.REQUIRE,
    "with": TokenType.WITH,
    "python": TokenType.PYTHON,
    "fn": TokenType.FN,
    "exports": TokenType.EXPORTS,

    # Closure
    "close": TokenType.CLOSE,
    "closeC0": TokenType.CLOSE_C0,
    "closeC1": TokenType.CLOSE_C1,

    # Boolean literals (Python style: True/False)
    "True": TokenType.BOOL_LITERAL,
    "False": TokenType.BOOL_LITERAL,
    # Also accept lowercase for compatibility
    "true": TokenType.BOOL_LITERAL,
    "false": TokenType.BOOL_LITERAL,

    # Tier 1 types
    "int": TokenType.TYPE_INT,
    "float": TokenType.TYPE_FLOAT,
    "string": TokenType.TYPE_STRING,
    "str": TokenType.TYPE_STRING,  # Python alias
    "bool": TokenType.TYPE_BOOL,
    "point": TokenType.TYPE_POINT,
    "point2d": TokenType.TYPE_POINT2D,
    "point3d": TokenType.TYPE_POINT3D,
    "vector": TokenType.TYPE_VECTOR,
    "vector2d": TokenType.TYPE_VECTOR2D,
    "vector3d": TokenType.TYPE_VECTOR3D,
    "transform": TokenType.TYPE_TRANSFORM,

    # Tier 2 types
    "line_segment": TokenType.TYPE_LINE_SEGMENT,
    "arc": TokenType.TYPE_ARC,
    "circle": TokenType.TYPE_CIRCLE,
    "ellipse": TokenType.TYPE_ELLIPSE,
    "parabola": TokenType.TYPE_PARABOLA,
    "hyperbola": TokenType.TYPE_HYPERBOLA,
    "catmullrom": TokenType.TYPE_CATMULLROM,
    "nurbs": TokenType.TYPE_NURBS,
    "bezier": TokenType.TYPE_BEZIER,

    # Tier 3 types
    "path2d": TokenType.TYPE_PATH2D,
    "path3d": TokenType.TYPE_PATH3D,
    "profile2d": TokenType.TYPE_PROFILE2D,
    "region2d": TokenType.TYPE_REGION2D,
    "loop3d": TokenType.TYPE_LOOP3D,

    # Tier 4 types
    "surface": TokenType.TYPE_SURFACE,
    "shell": TokenType.TYPE_SHELL,

    # Tier 5 types
    "solid": TokenType.TYPE_SOLID,

    # Generic types
    "list": TokenType.TYPE_LIST,
    "dict": TokenType.TYPE_DICT,
}


# Set of deprecated keywords for helpful error messages
DEPRECATED_KEYWORDS: set[str] = {
    "command",  # Use 'def' instead
    "let",      # No longer needed, just use 'name: type = value' or 'name = value'
    "require",  # Use 'assert' instead
    "fn",       # Use 'def' instead
    "exports",  # Use @native decorator instead
}

# Mapping of deprecated keyword to suggestion
DEPRECATED_SUGGESTIONS: dict[str, str] = {
    "command": "Use 'def' to define functions",
    "let": "Variable declarations no longer need 'let'. Use 'name: type = value' or 'name = value'",
    "require": "Use 'assert condition, \"message\"' instead",
    "fn": "Use 'def' to define functions",
    "exports": "Use '@native' decorator instead of 'native python { } exports { }'",
    "with": "Use keyword arguments with emit: 'emit value, name=\"x\", material=\"y\"'",
}


def is_type_token(token_type: TokenType) -> bool:
    """Check if a token type represents a type keyword."""
    return token_type.name.startswith("TYPE_")


def is_deprecated_keyword(keyword: str) -> bool:
    """Check if a keyword is deprecated."""
    return keyword in DEPRECATED_KEYWORDS


def get_deprecation_message(keyword: str) -> Optional[str]:
    """Get the deprecation message for a keyword, if any."""
    return DEPRECATED_SUGGESTIONS.get(keyword)
