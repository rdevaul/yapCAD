"""
DSL-specific exceptions and error handling.

Error code ranges (from Phase 3 roadmap):
- E0xx: Lexer errors
- E1xx: Parser errors
- E2xx: Type errors
- E3xx: Semantic errors
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional, List
from .tokens import SourceSpan, SourceLocation


class ErrorSeverity(Enum):
    """Severity levels for diagnostics."""
    ERROR = "error"
    WARNING = "warning"
    INFO = "info"
    HINT = "hint"


@dataclass
class Diagnostic:
    """A single diagnostic message (error, warning, etc.)."""
    code: str                       # E001, E101, etc.
    message: str                    # Human-readable message
    severity: ErrorSeverity
    span: SourceSpan
    source_line: Optional[str] = None   # The actual line of source code
    hints: List[str] = field(default_factory=list)
    related: List["Diagnostic"] = field(default_factory=list)

    def format(self, show_source: bool = True) -> str:
        """Format the diagnostic for display."""
        parts = []

        # Header: location: severity[code]: message
        loc = f"{self.span.start}"
        parts.append(f"{loc}: {self.severity.value}[{self.code}]: {self.message}")

        # Source line with caret
        if show_source and self.source_line is not None:
            parts.append(f"  |")
            line_num = str(self.span.start.line)
            parts.append(f"{line_num:>3} | {self.source_line}")

            # Caret pointing to the error
            col = self.span.start.column
            end_col = self.span.end.column if self.span.start.line == self.span.end.line else len(self.source_line) + 1
            underline_len = max(1, end_col - col)
            parts.append(f"    | {' ' * (col - 1)}{'^' * underline_len}")

        # Hints
        for hint in self.hints:
            parts.append(f"    = hint: {hint}")

        # Related diagnostics
        for related in self.related:
            parts.append(f"    --> {related.span.start}: {related.message}")

        return "\n".join(parts)

    def to_json(self) -> dict:
        """Convert to JSON-serializable dict for tooling integration."""
        return {
            "code": self.code,
            "message": self.message,
            "severity": self.severity.value,
            "range": {
                "start": {
                    "line": self.span.start.line,
                    "column": self.span.start.column,
                    "offset": self.span.start.offset,
                },
                "end": {
                    "line": self.span.end.line,
                    "column": self.span.end.column,
                    "offset": self.span.end.offset,
                },
            },
            "hints": self.hints,
            "related": [r.to_json() for r in self.related],
        }


class DslError(Exception):
    """Base exception for DSL errors."""

    def __init__(self, diagnostic: Diagnostic):
        self.diagnostic = diagnostic
        super().__init__(diagnostic.message)

    def __str__(self) -> str:
        return self.diagnostic.format()


class LexerError(DslError):
    """Error during lexical analysis (E0xx)."""
    pass


class ParserError(DslError):
    """Error during parsing (E1xx)."""
    pass


class TypeError(DslError):
    """Error during type checking (E2xx)."""
    pass


class SemanticError(DslError):
    """Error during semantic analysis (E3xx)."""
    pass


# --- Lexer error codes ---

def error_unexpected_character(char: str, span: SourceSpan, source_line: str = None) -> LexerError:
    """E001: Unexpected character."""
    diag = Diagnostic(
        code="E001",
        message=f"unexpected character '{char}'",
        severity=ErrorSeverity.ERROR,
        span=span,
        source_line=source_line,
    )
    return LexerError(diag)


def error_unterminated_string(span: SourceSpan, source_line: str = None) -> LexerError:
    """E002: Unterminated string literal."""
    diag = Diagnostic(
        code="E002",
        message="unterminated string literal",
        severity=ErrorSeverity.ERROR,
        span=span,
        source_line=source_line,
        hints=["string literals must be closed with matching quotes"],
    )
    return LexerError(diag)


def error_unterminated_multiline_string(span: SourceSpan, source_line: str = None) -> LexerError:
    """E003: Unterminated multi-line string."""
    diag = Diagnostic(
        code="E003",
        message='unterminated multi-line string (expected closing """)',
        severity=ErrorSeverity.ERROR,
        span=span,
        source_line=source_line,
    )
    return LexerError(diag)


def error_unterminated_comment(span: SourceSpan, source_line: str = None) -> LexerError:
    """E004: Unterminated multi-line comment."""
    diag = Diagnostic(
        code="E004",
        message="unterminated multi-line comment (expected closing */)",
        severity=ErrorSeverity.ERROR,
        span=span,
        source_line=source_line,
    )
    return LexerError(diag)


def error_invalid_escape_sequence(seq: str, span: SourceSpan, source_line: str = None) -> LexerError:
    """E005: Invalid escape sequence in string."""
    diag = Diagnostic(
        code="E005",
        message=f"invalid escape sequence '\\{seq}'",
        severity=ErrorSeverity.ERROR,
        span=span,
        source_line=source_line,
        hints=["valid escape sequences: \\n, \\t, \\r, \\\", \\\\, \\0, \\x##, \\u####"],
    )
    return LexerError(diag)


def error_invalid_number_literal(text: str, span: SourceSpan, source_line: str = None) -> LexerError:
    """E006: Invalid number literal."""
    diag = Diagnostic(
        code="E006",
        message=f"invalid number literal '{text}'",
        severity=ErrorSeverity.ERROR,
        span=span,
        source_line=source_line,
    )
    return LexerError(diag)


def error_invalid_hex_literal(text: str, span: SourceSpan, source_line: str = None) -> LexerError:
    """E007: Invalid hexadecimal literal."""
    diag = Diagnostic(
        code="E007",
        message=f"invalid hexadecimal literal '{text}'",
        severity=ErrorSeverity.ERROR,
        span=span,
        source_line=source_line,
        hints=["hex literals must contain at least one hex digit: 0x1, 0xFF, etc."],
    )
    return LexerError(diag)


def error_invalid_binary_literal(text: str, span: SourceSpan, source_line: str = None) -> LexerError:
    """E008: Invalid binary literal."""
    diag = Diagnostic(
        code="E008",
        message=f"invalid binary literal '{text}'",
        severity=ErrorSeverity.ERROR,
        span=span,
        source_line=source_line,
        hints=["binary literals must contain only 0 and 1: 0b101, 0b1111, etc."],
    )
    return LexerError(diag)


# --- Parser error codes ---

def error_unexpected_token(expected: str, found: str, span: SourceSpan,
                           source_line: str = None) -> ParserError:
    """E101: Unexpected token."""
    diag = Diagnostic(
        code="E101",
        message=f"expected {expected}, found {found}",
        severity=ErrorSeverity.ERROR,
        span=span,
        source_line=source_line,
    )
    return ParserError(diag)


def error_unexpected_eof(expected: str, span: SourceSpan) -> ParserError:
    """E102: Unexpected end of file."""
    diag = Diagnostic(
        code="E102",
        message=f"unexpected end of file, expected {expected}",
        severity=ErrorSeverity.ERROR,
        span=span,
    )
    return ParserError(diag)


def error_invalid_expression(span: SourceSpan, source_line: str = None) -> ParserError:
    """E103: Invalid expression."""
    diag = Diagnostic(
        code="E103",
        message="invalid expression",
        severity=ErrorSeverity.ERROR,
        span=span,
        source_line=source_line,
    )
    return ParserError(diag)


# --- Type error codes ---

def error_type_mismatch(expected: str, found: str, span: SourceSpan,
                        source_line: str = None) -> TypeError:
    """E201: Type mismatch."""
    diag = Diagnostic(
        code="E201",
        message=f"type mismatch: expected '{expected}', found '{found}'",
        severity=ErrorSeverity.ERROR,
        span=span,
        source_line=source_line,
    )
    return TypeError(diag)


def error_undefined_identifier(name: str, span: SourceSpan,
                               source_line: str = None) -> TypeError:
    """E202: Undefined identifier."""
    diag = Diagnostic(
        code="E202",
        message=f"undefined identifier '{name}'",
        severity=ErrorSeverity.ERROR,
        span=span,
        source_line=source_line,
    )
    return TypeError(diag)


# --- Semantic error codes ---

def error_require_failed(message: str, span: SourceSpan,
                         source_line: str = None) -> SemanticError:
    """E301: Require constraint failed."""
    diag = Diagnostic(
        code="E301",
        message=f"require constraint failed: {message}",
        severity=ErrorSeverity.ERROR,
        span=span,
        source_line=source_line,
    )
    return SemanticError(diag)


# --- Warnings ---

def warning_python_block(span: SourceSpan, source_line: str = None) -> Diagnostic:
    """W001: Python block requires manual approval."""
    return Diagnostic(
        code="W001",
        message="python block requires manual approval before execution",
        severity=ErrorSeverity.WARNING,
        span=span,
        source_line=source_line,
        hints=["packages containing python blocks are marked as 'requires-review'"],
    )


class DiagnosticCollector:
    """Collects diagnostics during compilation."""

    def __init__(self, max_errors: int = 20):
        self.diagnostics: List[Diagnostic] = []
        self.max_errors = max_errors
        self._error_count = 0

    def add(self, diagnostic: Diagnostic) -> None:
        """Add a diagnostic."""
        self.diagnostics.append(diagnostic)
        if diagnostic.severity == ErrorSeverity.ERROR:
            self._error_count += 1

    def add_error(self, error: DslError) -> None:
        """Add an error exception as a diagnostic."""
        self.add(error.diagnostic)

    @property
    def error_count(self) -> int:
        return self._error_count

    @property
    def warning_count(self) -> int:
        return sum(1 for d in self.diagnostics if d.severity == ErrorSeverity.WARNING)

    @property
    def has_errors(self) -> bool:
        return self._error_count > 0

    @property
    def has_warnings(self) -> bool:
        return self.warning_count > 0

    @property
    def should_stop(self) -> bool:
        """Check if we've hit the max error limit."""
        return self._error_count >= self.max_errors

    def format_all(self, show_source: bool = True) -> str:
        """Format all diagnostics for display."""
        parts = [d.format(show_source) for d in self.diagnostics]
        if self._error_count > 0:
            parts.append(f"\n{self._error_count} error(s), {self.warning_count} warning(s)")
        elif self.warning_count > 0:
            parts.append(f"\n{self.warning_count} warning(s)")
        return "\n\n".join(parts)

    def to_json(self) -> dict:
        """Convert all diagnostics to JSON format."""
        return {
            "diagnostics": [d.to_json() for d in self.diagnostics],
            "error_count": self._error_count,
            "warning_count": self.warning_count,
        }
