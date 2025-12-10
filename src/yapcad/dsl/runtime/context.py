"""
Execution context for DSL interpreter.

Manages variable scopes, tracks execution state, and collects errors/warnings.
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Callable
from contextlib import contextmanager

from .values import Value, RequireFailure, EmitResult
from ..types import Type
from ..errors import Diagnostic, DiagnosticCollector, ErrorSeverity
from ..tokens import SourceSpan


@dataclass
class Scope:
    """
    A single scope containing variable bindings.

    Scopes form a chain via the `parent` field for lexical scoping.
    """
    variables: Dict[str, Value] = field(default_factory=dict)
    parent: Optional["Scope"] = None
    name: str = "anonymous"  # For debugging

    def get(self, name: str) -> Optional[Value]:
        """Look up a variable in this scope or parent scopes."""
        if name in self.variables:
            return self.variables[name]
        if self.parent:
            return self.parent.get(name)
        return None

    def set(self, name: str, value: Value) -> None:
        """Set a variable in this scope (shadowing parent if exists)."""
        self.variables[name] = value

    def update(self, name: str, value: Value) -> bool:
        """
        Update an existing variable (mutable assignment).

        Searches up the scope chain to find where the variable is defined.
        Returns True if found and updated, False if not found.
        """
        if name in self.variables:
            self.variables[name] = value
            return True
        if self.parent:
            return self.parent.update(name, value)
        return False

    def contains(self, name: str) -> bool:
        """Check if a variable exists in this scope or parents."""
        return self.get(name) is not None


@dataclass
class ExecutionContext:
    """
    The full execution context for interpreting DSL code.

    Tracks:
    - Variable scopes
    - Require failures
    - Emit result
    - Diagnostics (errors/warnings)
    - Module/command being executed
    """
    # Scope chain
    current_scope: Scope = field(default_factory=lambda: Scope(name="global"))

    # Execution state
    module_name: str = ""
    command_name: str = ""
    parameters: Dict[str, Any] = field(default_factory=dict)

    # Results
    emit_result: Optional[EmitResult] = None
    require_failures: List[RequireFailure] = field(default_factory=list)

    # Diagnostics
    diagnostics: DiagnosticCollector = field(default_factory=DiagnosticCollector)

    # Source tracking for error messages
    source_lines: List[str] = field(default_factory=list)

    # Control flow flags
    _should_return: bool = False
    _return_value: Optional[Value] = None

    def get_variable(self, name: str) -> Optional[Value]:
        """Look up a variable in the current scope chain."""
        return self.current_scope.get(name)

    def set_variable(self, name: str, value: Value) -> None:
        """Define a new variable in the current scope."""
        self.current_scope.set(name, value)

    def update_variable(self, name: str, value: Value) -> bool:
        """Update an existing variable (for assignment statements)."""
        return self.current_scope.update(name, value)

    @contextmanager
    def new_scope(self, name: str = "block"):
        """
        Context manager to create a new nested scope.

        Usage:
            with ctx.new_scope("for-loop"):
                # variables defined here are local to this scope
                ctx.set_variable("i", int_val(0))
        """
        old_scope = self.current_scope
        self.current_scope = Scope(parent=old_scope, name=name)
        try:
            yield self.current_scope
        finally:
            self.current_scope = old_scope

    def add_require_failure(self, message: str, expression_text: str = None) -> None:
        """Record a require constraint failure."""
        self.require_failures.append(RequireFailure(message, expression_text))

    def set_emit(self, value: Value, metadata: Dict[str, Any] = None) -> None:
        """Set the emit result for this command."""
        self.emit_result = EmitResult(value, metadata or {})

    def add_error(self, message: str, span: SourceSpan) -> None:
        """Add an error diagnostic."""
        source_line = self._get_source_line(span.start.line)
        diag = Diagnostic(
            code="E400",  # Runtime errors use E4xx
            message=message,
            severity=ErrorSeverity.ERROR,
            span=span,
            source_line=source_line,
        )
        self.diagnostics.add(diag)

    def add_warning(self, message: str, span: SourceSpan) -> None:
        """Add a warning diagnostic."""
        source_line = self._get_source_line(span.start.line)
        diag = Diagnostic(
            code="W400",
            message=message,
            severity=ErrorSeverity.WARNING,
            span=span,
            source_line=source_line,
        )
        self.diagnostics.add(diag)

    def _get_source_line(self, line_num: int) -> Optional[str]:
        """Get a source line for error messages."""
        if 1 <= line_num <= len(self.source_lines):
            return self.source_lines[line_num - 1]
        return None

    @property
    def has_errors(self) -> bool:
        """Check if any errors occurred."""
        return self.diagnostics.has_errors or len(self.require_failures) > 0

    @property
    def has_warnings(self) -> bool:
        """Check if any warnings occurred."""
        return self.diagnostics.has_warnings

    def signal_return(self, value: Value) -> None:
        """Signal an early return from a command."""
        self._should_return = True
        self._return_value = value

    @property
    def should_return(self) -> bool:
        """Check if early return was signaled."""
        return self._should_return

    @property
    def return_value(self) -> Optional[Value]:
        """Get the return value if early return was signaled."""
        return self._return_value

    def clear_return(self) -> None:
        """Clear the return signal (used after handling return)."""
        self._should_return = False
        self._return_value = None


def create_context(
    module_name: str,
    command_name: str,
    parameters: Dict[str, Value],
    source: str = "",
) -> ExecutionContext:
    """
    Create a new execution context for a command invocation.

    Args:
        module_name: The module being executed
        command_name: The command being executed
        parameters: Parameter values passed to the command
        source: The source code (for error messages)

    Returns:
        A fresh ExecutionContext with parameters bound in scope
    """
    ctx = ExecutionContext(
        module_name=module_name,
        command_name=command_name,
        parameters={k: v.data for k, v in parameters.items()},
        source_lines=source.split('\n') if source else [],
    )

    # Bind parameters to the global scope
    for name, value in parameters.items():
        ctx.set_variable(name, value)

    return ctx
