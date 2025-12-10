"""
DSL Runtime - Tree-walking interpreter for DSL execution.

This module provides:
- Interpreter: Executes DSL commands to generate geometry
- Value: Runtime value wrappers with type metadata
- ExecutionContext: Variable scope management
- Provenance: Invocation metadata for audit trails
- BuiltinRegistry: Built-in function implementations
"""

from .values import (
    Value,
    EmitResult,
    RequireFailure,
    int_val,
    float_val,
    bool_val,
    string_val,
    list_val,
    dict_val,
    none_val,
    point_val,
    vector_val,
    transform_val,
    solid_val,
    region2d_val,
    surface_val,
    shell_val,
    wrap_value,
    unwrap_value,
    unwrap_values,
    check_type,
    coerce_numeric,
)

from .context import (
    Scope,
    ExecutionContext,
    create_context,
)

from .builtins import (
    BuiltinFunction,
    BuiltinRegistry,
    get_builtin_registry,
    call_builtin,
    call_method,
)

from .interpreter import (
    Interpreter,
    ExecutionResult,
    execute,
    compile_and_run,
)

from .provenance import (
    Provenance,
    create_provenance,
    compute_source_signature,
    verify_source_signature,
)

__all__ = [
    # Values
    'Value',
    'EmitResult',
    'RequireFailure',
    'int_val',
    'float_val',
    'bool_val',
    'string_val',
    'list_val',
    'dict_val',
    'none_val',
    'point_val',
    'vector_val',
    'transform_val',
    'solid_val',
    'region2d_val',
    'surface_val',
    'shell_val',
    'wrap_value',
    'unwrap_value',
    'unwrap_values',
    'check_type',
    'coerce_numeric',

    # Context
    'Scope',
    'ExecutionContext',
    'create_context',

    # Builtins
    'BuiltinFunction',
    'BuiltinRegistry',
    'get_builtin_registry',
    'call_builtin',
    'call_method',

    # Interpreter
    'Interpreter',
    'ExecutionResult',
    'execute',
    'compile_and_run',

    # Provenance
    'Provenance',
    'create_provenance',
    'compute_source_signature',
    'verify_source_signature',
]
