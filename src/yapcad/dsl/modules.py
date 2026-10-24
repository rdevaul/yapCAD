"""Resolution of the DSL's ``use`` statements.

Both the type checker and the interpreter need to know what a module's
``use`` statements bring into scope, so resolution lives here rather than in
either of them.

Imported commands are reached through the importing module's own namespace:

- ``use shapes`` merges ``shapes``'s commands under their own names. A name
  that is already taken is an error rather than a silent first-wins, because
  silently binding a call to another module's command of the same name
  produces wrong geometry with no diagnostic.
- ``use shapes as s`` puts them behind a dot instead: ``s.CUBE(...)``. Calls
  *between* the imported module's own commands are rewritten to match, so a
  library whose commands call each other keeps working when it is aliased.

A command's defining file, not the file that imported it, decides whether two
imports of the same name are a conflict or the same command arriving twice by
different routes (the diamond case: two modules that both ``use base``).
"""

from __future__ import annotations

import copy
import dataclasses
import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple

from .ast import (
    ExportUseStatement,
    FunctionCall,
    Identifier,
    Module,
    UseStatement,
)

__all__ = [
    "UseMergeResult",
    "UseError",
    "resolve_uses",
    "clear_module_cache",
]


class UseError(Exception):
    """A ``use`` statement that cannot be honoured."""


@dataclass
class UseMergeResult:
    success: bool
    merged_module: Optional[Module] = None
    error_message: Optional[str] = None


# realpath -> (stamp, module, origins). The stamp is (mtime_ns, size): a
# cache that never noticed an edit had every later build in a long-lived
# process -- the viewer, a pytest session, tools/build_ycpkg.py -- reusing
# stale geometry.
_Stamp = Tuple[int, int]
_Origins = Dict[str, Optional[str]]
_CACHE: Dict[str, Tuple[_Stamp, Module, _Origins]] = {}


def clear_module_cache() -> None:
    """Forget every cached module.

    Cached modules are invalidated automatically when the file changes on
    disk; this is for hosts that need to force a reload anyway.
    """
    _CACHE.clear()


def _stamp(real_path: str) -> Tuple[int, int]:
    info = os.stat(real_path)
    return (info.st_mtime_ns, info.st_size)


def _cached(real_path: str):
    hit = _CACHE.get(real_path)
    if hit is None:
        return None
    stamp, module, origins = hit
    try:
        if stamp != _stamp(real_path):
            _CACHE.pop(real_path, None)
            return None
    except OSError:
        _CACHE.pop(real_path, None)
        return None
    return module, origins


def _qualify(name: str, namespace: str) -> str:
    """Put ``name`` behind ``namespace``.

    A name that is already qualified is left alone: a module's own aliases
    are its business, and re-qualifying them would produce ``a.b.NAME``,
    which the call syntax cannot express.
    """
    return name if "." in name else f"{namespace}.{name}"


def _rewrite_calls(node, mapping: Dict[str, str]) -> None:
    """Point unqualified calls at their namespaced names, in place.

    A generic walk over dataclass fields rather than a per-node-type
    visitor, so that an expression kind added later cannot silently escape
    the rewrite. Callers pass a deep copy: the originals may be shared with
    the module cache.
    """
    if isinstance(node, FunctionCall) and isinstance(node.callee, Identifier):
        replacement = mapping.get(node.callee.name)
        if replacement is not None:
            node.callee.name = replacement

    if dataclasses.is_dataclass(node):
        for field_ in dataclasses.fields(node):
            _rewrite_calls(getattr(node, field_.name, None), mapping)
    elif isinstance(node, (list, tuple)):
        for item in node:
            _rewrite_calls(item, mapping)
    elif isinstance(node, dict):
        for item in node.values():
            _rewrite_calls(item, mapping)


def _namespaced_commands(module: Module, namespace: str) -> List:
    """``module``'s commands, copied under ``namespace``.

    Calls between them are rewritten too, so that an imported command that
    calls a sibling still reaches the sibling and not some same-named
    command in the importing module.
    """
    mapping = {
        cmd.name: f"{namespace}.{cmd.name}"
        for cmd in module.commands
        if "." not in cmd.name
    }
    out = []
    for cmd in module.commands:
        renamed = copy.deepcopy(cmd)
        renamed.name = mapping.get(cmd.name, cmd.name)
        _rewrite_calls(renamed, mapping)
        out.append(renamed)
    return out


def _locate(module_path: List[str], source_dir: Optional[str]) -> str:
    """The file ``use module_path`` refers to, or raise ``UseError``."""
    path_text = ".".join(module_path)
    if not source_dir:
        raise UseError(
            f"Cannot resolve 'use {path_text}': the source directory is not "
            "known. Pass source_path= (or use the CLI, which does) so that "
            "imports can be located relative to the importing file."
        )

    relative = os.path.join(*module_path) + ".dsl"
    root = os.path.realpath(source_dir)
    candidate = os.path.realpath(os.path.join(root, relative))

    # Keep imports inside the project. 'use ..foo' resolves to a '..'
    # path component, which would otherwise read any .dsl file on the box.
    if os.path.commonpath([root, candidate]) != root:
        raise UseError(
            f"Cannot resolve 'use {path_text}': it points outside the source "
            f"directory. Imports must live under '{root}', but this resolves "
            f"to '{candidate}'."
        )
    if not os.path.isfile(candidate):
        raise UseError(
            f"Cannot resolve 'use {path_text}': no file '{relative}' "
            f"under '{root}'."
        )
    return candidate


def _load(real_path: str, loading: Set[str]):
    """Parse and check the module at ``real_path``, with its imports merged.

    Returns ``(module, origins)``, where ``origins`` maps each command name
    to the file that defines it, or None when this module defines it.
    """
    from .checker import check
    from .lexer import tokenize
    from .parser import parse

    hit = _cached(real_path)
    if hit is not None:
        return hit

    source_text = open(real_path, encoding="utf-8").read()
    module = parse(tokenize(source_text), source=source_text)

    source_dir = os.path.dirname(real_path)
    module, origins = _resolve(module, source_dir, loading)

    result = check(module)
    if result.has_errors:
        messages = "; ".join(
            diag.message
            for diag in result.diagnostics
            if diag.severity.name == "ERROR"
        )
        raise UseError(f"{os.path.basename(real_path)}: {messages}")

    _CACHE[real_path] = (_stamp(real_path), module, origins)
    return module, origins


def _resolve(
    module: Module, source_dir: Optional[str], loading: Set[str]
) -> Tuple[Module, Dict[str, Optional[str]]]:
    """Merge ``module``'s imports, or raise ``UseError``."""
    if not module.uses:
        return module, {cmd.name: None for cmd in module.commands}

    commands = list(module.commands)
    origins: Dict[str, Optional[str]] = {cmd.name: None for cmd in commands}
    native_blocks = list(module.native_blocks)
    native_functions = list(module.native_functions)
    exports = list(module.exports)
    aliases: Dict[str, str] = {}

    for statement in module.uses:
        if isinstance(statement, ExportUseStatement):
            raise UseError(
                f"'export use {'.'.join(statement.module_path)}' is not "
                "supported. Re-exporting another module's commands has no "
                "implementation yet; use a plain 'use' statement."
            )
        if not isinstance(statement, UseStatement):
            raise UseError(
                f"Unsupported import statement: {type(statement).__name__}"
            )

        path_text = ".".join(statement.module_path)
        real_path = _locate(statement.module_path, source_dir)

        cached = _cached(real_path)
        if cached is not None:
            used, used_origins = cached
        elif real_path in loading:
            raise UseError(
                f"Circular import detected: 'use {path_text}' is already "
                "being loaded further up the chain."
            )
        else:
            loading.add(real_path)
            try:
                used, used_origins = _load(real_path, loading)
            except UseError as exc:
                raise UseError(
                    f"Failed to load module '{path_text}': {exc}"
                ) from exc
            except Exception as exc:
                raise UseError(
                    f"Failed to load module '{path_text}' from "
                    f"{real_path}: {exc}"
                ) from exc
            finally:
                # Leaving the path behind would make a later, unrelated
                # import of it look like a cycle.
                loading.discard(real_path)

        namespace = statement.alias
        if namespace:
            if aliases.get(namespace, real_path) != real_path:
                raise UseError(
                    f"Import alias '{namespace}' is already used for a "
                    f"different module ({aliases[namespace]}). Give "
                    f"'use {path_text}' a different alias."
                )
            aliases[namespace] = real_path
            incoming = _namespaced_commands(used, namespace)
        else:
            incoming = list(used.commands)

        for command in incoming:
            original = (
                _unqualify(command.name, namespace)
                if namespace else command.name
            )
            defined_in = used_origins.get(original) or real_path

            if command.name not in origins:
                origins[command.name] = defined_in
                commands.append(command)
                continue
            if origins[command.name] == defined_in:
                continue  # the same command reached by two routes

            previous = origins[command.name]
            where = (
                "defined in this module"
                if previous is None
                else f"already imported from {previous}"
            )
            raise UseError(
                f"'use {path_text}' would redefine command "
                f"'{command.name}', which is {where}. Alias one of the "
                f"imports ('use {path_text} as name') so that its commands "
                "sit behind their own namespace."
            )

        native_blocks.extend(used.native_blocks)
        native_functions.extend(used.native_functions)

    merged = Module(
        span=module.span,
        name=module.name,
        uses=[],
        functions=commands,
        native_blocks=native_blocks,
        native_functions=native_functions,
        exports=exports,
    )
    return merged, origins


def _unqualify(name: str, namespace: str) -> str:
    prefix = namespace + "."
    return name[len(prefix):] if name.startswith(prefix) else name


def resolve_uses(
    module: Module, source_dir: Optional[str]
) -> UseMergeResult:
    """Merge ``module``'s ``use`` statements into one module.

    Never raises for a bad import: the failure is reported in the result so
    that callers can turn it into a diagnostic instead of a traceback.
    """
    try:
        merged, _origins = _resolve(module, source_dir, set())
    except UseError as exc:
        return UseMergeResult(success=False, error_message=str(exc))
    return UseMergeResult(success=True, merged_module=merged)
