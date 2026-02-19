from __future__ import annotations

from typing import Any, Dict, Optional, Tuple

import anyio

from yapcad.dsl.runtime.interpreter import compile_and_run
from yapcad.dsl.lexer import tokenize
from yapcad.dsl.parser import parse

from .ast_serialize import ast_to_dict


class DslTimeoutError(RuntimeError):
    pass


async def eval_dsl(
    *,
    source: str,
    command: str,
    parameters: Dict[str, Any],
    timeout_s: float = 30.0,
    recursion_limit: Optional[int] = None,
):
    """Evaluate DSL source with a hard timeout.

    Runs compilation + execution in a worker thread.
    """

    def _run_sync():
        return compile_and_run(
            source=source,
            command_name=command,
            parameters=parameters,
            recursion_limit=recursion_limit,
        )

    try:
        with anyio.fail_after(timeout_s):
            return await anyio.to_thread.run_sync(_run_sync)
    except TimeoutError as exc:  # anyio raises built-in TimeoutError
        raise DslTimeoutError(f"DSL evaluation exceeded timeout ({timeout_s}s)") from exc


def parse_dsl(*, source: str) -> Tuple[dict, Any]:
    """Parse DSL source into AST (module) and return serialized dict + module."""
    tokens = tokenize(source)
    module = parse(tokens, source=source)
    return ast_to_dict(module), module
