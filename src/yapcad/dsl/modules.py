"""Resolution of the DSL's ``use`` statements.

Both the type checker and the interpreter need to know what a module's
``use`` statements bring into scope, so resolution lives here rather than in
either of them.

Every imported command keeps a namespace of its own, taken from the alias
when there is one and otherwise from the last component of the module path:
``use shapes`` and ``use parts.shapes as shapes`` both offer
``shapes.CUBE(...)``. Calls *between* an imported module's own commands are
rewritten to match, so a library whose commands call each other keeps working
however it is imported, and never reaches a same-named command somewhere else.

An unaliased ``use`` additionally offers each of its commands under the bare
name, which is what makes ``use shapes`` followed by ``CUBE(...)`` work. A
bare name that two imports both offer is left unbound rather than being
resolved first-come or rejected outright: an ambiguity is reported where the
name is *called*, so importing two modules that happen to share one command
name costs nothing until you call that name without saying which you mean. A
name the importing module defines itself always wins, and the import remains
reachable through its namespace.

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
    ImportScope,
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

_MISSING = object()


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

    # Name -> command, in the order the names were bound. A name that two
    # different modules both define is unbound again and recorded in
    # ``conflicts`` instead; whether that matters is up to the call sites.
    bound: Dict[str, object] = {}
    origins: _Origins = {}
    for command in module.commands:
        bound[command.name] = command
        origins[command.name] = None

    native_blocks = list(module.native_blocks)
    native_functions = list(module.native_functions)
    exports = list(module.exports)
    aliases: Dict[str, str] = {}
    conflicts: Dict[str, Tuple[Tuple[str, str], ...]] = {}
    # Bare name -> [(qualified name, defining file)], collected across the
    # use statements and bound afterwards, once it is known how many
    # modules offer each name.
    offered: Dict[str, List[Tuple[str, str]]] = {}

    def bind(name: str, command, defined_in: str) -> None:
        previous = origins.get(name, _MISSING)
        if previous is _MISSING and name not in conflicts:
            bound[name] = command
            origins[name] = defined_in
            return
        if previous == defined_in:
            return  # the same command reached by two routes
        bound.pop(name, None)
        origins.pop(name, None)
        providers = list(conflicts.get(name, ()))
        if isinstance(previous, str):
            _append_provider(providers, (name, previous))
        _append_provider(providers, (name, defined_in))
        conflicts[name] = tuple(providers)

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

        namespace = statement.alias or statement.module_path[-1]
        if statement.alias:
            # An alias is something the author chose and can change, so
            # pointing two of them at different modules is a mistake worth
            # refusing. Two module paths that merely end in the same
            # component are not: their namespaces share, and only a command
            # name they both define becomes ambiguous.
            if aliases.get(namespace, real_path) != real_path:
                raise UseError(
                    f"Import alias '{namespace}' is already used for a "
                    f"different module ({aliases[namespace]}). Give "
                    f"'use {path_text}' a different alias."
                )
            aliases[namespace] = real_path

        incoming = _namespaced_commands(used, namespace)
        for command, original in zip(incoming, used.commands):
            # The file that defines the command, which may be neither this
            # module nor the one being imported.
            defined_in = used_origins.get(original.name) or real_path
            bind(command.name, command, defined_in)
            if statement.alias is None and used_origins.get(
                original.name
            ) is None:
                _append_provider(
                    offered.setdefault(original.name, []),
                    (command.name, defined_in),
                )

        native_blocks.extend(used.native_blocks)
        native_functions.extend(used.native_functions)

    local_names = {command.name for command in module.commands}
    bindings: Dict[str, str] = {}
    shadowed: Dict[str, Tuple[str, ...]] = {}
    for bare, providers in offered.items():
        if bare in local_names:
            shadowed[bare] = tuple(
                qualified for qualified, _path in providers
            )
        elif len(providers) > 1:
            conflicts[bare] = tuple(providers)
        elif providers[0][0] in bound:
            bindings[bare] = providers[0][0]
        else:
            # The one module offering the bare name lost its namespace to a
            # conflict, so there is nothing to bind it to.
            conflicts[bare] = (
                conflicts.get(providers[0][0]) or tuple(providers)
            )

    # Point this module's own calls at the namespaced commands. Resolving
    # them here, while it is known which module each name came from, is what
    # lets the merged module be one flat namespace: a command imported from
    # here on carries calls that already say what they mean.
    if bindings:
        for command in module.commands:
            _rewrite_calls(command, bindings)

    merged = Module(
        span=module.span,
        name=module.name,
        uses=[],
        functions=list(bound.values()),
        native_blocks=native_blocks,
        native_functions=native_functions,
        exports=exports,
        imports=ImportScope(
            bindings=bindings, conflicts=conflicts, shadowed=shadowed
        ),
    )
    return merged, origins


def _append_provider(
    providers: List[Tuple[str, str]], provider: Tuple[str, str]
) -> None:
    """Add ``provider`` unless that file is already accounted for."""
    if all(path != provider[1] for _qualified, path in providers):
        providers.append(provider)


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
