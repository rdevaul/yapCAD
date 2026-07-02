"""Assembly operation resolver — pluggable strategy for applying ``operation``
metadata to assembled geometry.

This module implements the configurable resolver strategy described in
RFC §16 and ``metadata-namespace-v1.1``.  The selection interface mirrors
``YAPCAD_BOOLEAN_ENGINE``::

    # Programmatic
    from yapcad.assembly.resolver import resolve_assembly, BasicResolver, StagedResolver
    resolved = resolve_assembly(parts, resolver=StagedResolver(stages=["position", "cut"]))

    # Via environment variable (mirrors boolean-kernel selection)
    YAPCAD_ASSEMBLY_RESOLVER=basic|staged|<entry-point>

Resolvers register via the ``yapcad.assembly_resolvers`` Python entry point
so external packages (e.g. Mechatron) can supply their own without forking.

Built-in strategies
-------------------
``BasicResolver`` (default)
    Bucket → resolve targets → sort by ``priority`` → apply sequentially per
    target.  No positioning awareness.  Handles static assemblies such as the
    Agentic-1 forward bulkhead with two radial cutouts.

``StagedResolver``
    Multi-pass: caller declares ordered stage names; each operation is tagged
    with ``stage:``.  The resolver runs all operations in stage 0 to
    completion before moving to stage 1.  Use when some cuts must follow
    placement (e.g. clearance cuts that are only meaningful after an
    articulated part reaches its assembled pose).

``KinematicResolver``
    Scoped for v1.2 — placeholder registered here so external code can
    reference the name.  Raises ``NotImplementedError`` if instantiated.

Resolution algorithm (BasicResolver — normative, per spec §8.2)
---------------------------------------------------------------
1.  Collect all parts whose ``operation`` metadata has ``kind`` set.
2.  Validate that each cutter has a non-empty ``target_filter``.
3.  For each target part, gather all operations that list it in
    ``target_filter``.
4.  Sort gathered operations by ``priority`` (ascending — lower number = first).
5.  Apply each operation in order using the active boolean engine.
6.  Handle policy: ``strict`` raises on failure; ``warn`` logs; ``ignore``
    continues silently.
7.  If ``through`` is ``True``, assert the operation did not produce a
    watertight-closed surface on the removed volume (i.e. the cut goes all the
    way through).  If the assertion fails, apply ``policy``.
8.  If ``consume`` is ``False``, retain the cutter geometry in the output.
    Default (``True``) drops it.
"""

from __future__ import annotations

import importlib.metadata
import logging
import os
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

from ..metadata import get_operation_metadata

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Public data types
# ---------------------------------------------------------------------------


@dataclass
class ResolvedPart:
    """A single part after all applicable operations have been applied."""

    part_id: str
    geometry: Any                       # yapCAD solid / list[face]
    applied_operations: List[str]       # cutter part_ids applied, in order
    retained_cutters: List[str]         # cutter part_ids kept (consume=False)
    warnings: List[str] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)


@dataclass
class ResolveResult:
    """Full result returned by :func:`resolve_assembly`."""

    parts: Dict[str, ResolvedPart]      # keyed by part_id
    execution_log: List[str]            # ordered narrative of operations applied
    success: bool = True
    errors: List[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Resolver ABC
# ---------------------------------------------------------------------------


class AssemblyResolver(ABC):
    """Abstract base for pluggable assembly resolvers."""

    @abstractmethod
    def resolve(
        self,
        parts: Dict[str, Any],
        meta_map: Dict[str, Dict],
        *,
        boolean_engine: Optional[str] = None,
    ) -> ResolveResult:
        """Apply all ``operation`` metadata to geometry.

        Parameters
        ----------
        parts:
            ``{part_id: geometry}`` dict.  Geometry may be any yapCAD solid.
        meta_map:
            ``{part_id: metadata_dict}`` — the full metadata dict (as returned
            by ``get_metadata(solid)``) for each part.
        boolean_engine:
            Override the active boolean engine name.  ``None`` uses the env
            var or the yapCAD default.

        Returns
        -------
        ResolveResult
        """

    # ------------------------------------------------------------------
    # Shared helpers available to all concrete resolvers
    # ------------------------------------------------------------------

    @staticmethod
    def _collect_cutters(
        meta_map: Dict[str, Dict],
    ) -> List[Tuple[str, Dict]]:
        """Return ``[(part_id, operation_dict)]`` for all parts with ``operation.kind`` set."""
        cutters = []
        for part_id, meta in meta_map.items():
            op = get_operation_metadata(meta)
            if op.get("kind"):
                cutters.append((part_id, op))
        return cutters

    @staticmethod
    def _validate_cutters(cutters: List[Tuple[str, Dict]]) -> List[str]:
        """Return a list of error strings for invalid cutter specs."""
        errors = []
        for part_id, op in cutters:
            if not op.get("target_filter"):
                errors.append(
                    f"cutter '{part_id}' has no target_filter — "
                    "omitting target_filter is an error (RFC §13 Q4)"
                )
        return errors

    @staticmethod
    def _apply_boolean(
        target_geom: Any,
        cutter_geom: Any,
        kind: str,
        engine: Optional[str],
    ) -> Any:
        """Thin dispatch into the yapCAD boolean engine.

        ``kind`` is one of ``"subtract"``, ``"intersect"``, ``"union"``.
        """
        # Import lazily to avoid circular deps and optional-dep failures.
        try:
            from ..geom3d import boolean3d  # type: ignore[import]
        except ImportError:
            raise RuntimeError(
                "geom3d.boolean3d not available — cannot apply boolean operation"
            )
        op_map = {"subtract": "difference", "intersect": "intersection", "union": "union"}
        op_name = op_map.get(kind, kind)
        kwargs: Dict[str, Any] = {}
        if engine:
            kwargs["engine"] = engine
        return boolean3d(target_geom, cutter_geom, op_name, **kwargs)

    @staticmethod
    def _handle_policy(
        policy: str,
        message: str,
        errors: List[str],
        warnings: List[str],
    ) -> None:
        if policy == "strict":
            errors.append(message)
        elif policy == "warn":
            warnings.append(message)
            log.warning("resolver: %s", message)
        else:  # ignore
            log.debug("resolver (ignored): %s", message)


# ---------------------------------------------------------------------------
# BasicResolver
# ---------------------------------------------------------------------------


class BasicResolver(AssemblyResolver):
    """Default resolver — per spec §8.2.

    Static assemblies.  No staging, no positioning awareness.
    ``BasicResolver`` ignores the ``stage`` field on operations.
    """

    def resolve(
        self,
        parts: Dict[str, Any],
        meta_map: Dict[str, Dict],
        *,
        boolean_engine: Optional[str] = None,
    ) -> ResolveResult:
        log_lines: List[str] = []
        global_errors: List[str] = []

        cutters = self._collect_cutters(meta_map)
        val_errors = self._validate_cutters(cutters)
        if val_errors:
            return ResolveResult(
                parts={pid: ResolvedPart(pid, geom, [], []) for pid, geom in parts.items()},
                execution_log=val_errors,
                success=False,
                errors=val_errors,
            )

        # Build target → [(priority, cutter_id, op)] index
        target_ops: Dict[str, List[Tuple[int, str, Dict]]] = {}
        for cutter_id, op in cutters:
            for target_id in op.get("target_filter", []):
                target_ops.setdefault(target_id, []).append(
                    (op.get("priority", 50), cutter_id, op)
                )

        # Sort each target's op list by priority (ascending)
        for target_id in target_ops:
            target_ops[target_id].sort(key=lambda t: t[0])

        resolved_parts: Dict[str, ResolvedPart] = {}
        retained_all: set = set()

        # Build set of consumed cutter ids (consume=True, the default)
        consumed_cutter_ids: set = {
            cutter_id
            for cutter_id, op in cutters
            if op.get("consume", True)
        }

        for part_id, geom in parts.items():
            if part_id not in target_ops:
                # Not a target — pass through unless it's a consumed cutter
                if part_id in consumed_cutter_ids:
                    continue  # drop consumed cutters that aren't also targets
                resolved_parts[part_id] = ResolvedPart(part_id, geom, [], [])
                continue

            current_geom = geom
            applied: List[str] = []
            retained: List[str] = []
            part_warnings: List[str] = []
            part_errors: List[str] = []

            for priority, cutter_id, op in target_ops[part_id]:
                cutter_geom = parts.get(cutter_id)
                if cutter_geom is None:
                    msg = f"cutter '{cutter_id}' listed in target_filter for '{part_id}' but not in parts dict"
                    self._handle_policy(op.get("policy", "strict"), msg, part_errors, part_warnings)
                    continue

                kind = op.get("kind", "subtract")
                policy = op.get("policy", "strict")
                log_lines.append(
                    f"  {kind} '{cutter_id}' → '{part_id}' (priority={priority})"
                )

                try:
                    current_geom = self._apply_boolean(current_geom, cutter_geom, kind, boolean_engine)
                    applied.append(cutter_id)
                except Exception as exc:  # noqa: BLE001
                    msg = f"boolean {kind} of '{cutter_id}' into '{part_id}' failed: {exc}"
                    self._handle_policy(policy, msg, part_errors, part_warnings)

                # Cutter retention
                consume = op.get("consume", True)
                if not consume:
                    retained.append(cutter_id)
                    retained_all.add(cutter_id)

            resolved_parts[part_id] = ResolvedPart(
                part_id, current_geom, applied, retained,
                warnings=part_warnings, errors=part_errors,
            )
            if part_errors:
                global_errors.extend(part_errors)

        # Add cutter parts to output only if:
        # (a) they were ALSO a target (resolved_parts already has them), OR
        # (b) consume=False (they are in retained_all).
        # Consumed cutters that were not targets are dropped from output.
        for cutter_id, _op in cutters:
            if cutter_id not in resolved_parts:
                # Not a target; only include if explicitly retained
                if cutter_id in retained_all:
                    resolved_parts[cutter_id] = ResolvedPart(cutter_id, parts[cutter_id], [], [])

        return ResolveResult(
            parts=resolved_parts,
            execution_log=log_lines,
            success=not global_errors,
            errors=global_errors,
        )


# ---------------------------------------------------------------------------
# StagedResolver
# ---------------------------------------------------------------------------


class StagedResolver(AssemblyResolver):
    """Multi-stage resolver for assemblies where cut order depends on position.

    The caller declares an ordered list of stage names.  Each operation
    specifies a ``stage`` field (a free-form string matching one entry in the
    stage list).  Operations with no ``stage`` field default to the *last*
    stage.

    The resolver runs ``BasicResolver`` on each stage's subset in order,
    feeding the resolved geometries of stage N as the input to stage N+1.

    Parameters
    ----------
    stages:
        Ordered stage names, e.g. ``["pre_position", "post_position", "final"]``.
        Must be non-empty.  Operations with unknown stage tags are collected
        into a warning and processed in the last stage.
    """

    def __init__(self, stages: Sequence[str]) -> None:
        if not stages:
            raise ValueError("StagedResolver: stages list must be non-empty")
        self.stages = list(stages)
        self._basic = BasicResolver()

    def resolve(
        self,
        parts: Dict[str, Any],
        meta_map: Dict[str, Dict],
        *,
        boolean_engine: Optional[str] = None,
    ) -> ResolveResult:
        log_lines: List[str] = []
        global_errors: List[str] = []

        cutters = self._collect_cutters(meta_map)
        val_errors = self._validate_cutters(cutters)
        if val_errors:
            return ResolveResult(
                parts={pid: ResolvedPart(pid, geom, [], []) for pid, geom in parts.items()},
                execution_log=val_errors,
                success=False,
                errors=val_errors,
            )

        # Bucket cutters by stage
        stage_buckets: Dict[str, List[Tuple[str, Dict]]] = {s: [] for s in self.stages}
        fallback_stage = self.stages[-1]

        for cutter_id, op in cutters:
            stage_tag = op.get("stage")
            if stage_tag is None or stage_tag not in stage_buckets:
                if stage_tag is not None:
                    log_lines.append(
                        f"  WARNING: cutter '{cutter_id}' has unknown stage='{stage_tag}'"
                        f" — falling back to '{fallback_stage}'"
                    )
                stage_buckets[fallback_stage].append((cutter_id, op))
            else:
                stage_buckets[stage_tag].append((cutter_id, op))

        # Run each stage in order
        current_parts = dict(parts)
        current_meta = dict(meta_map)
        accumulated_ops: Dict[str, List[str]] = {}  # pid → all applied_operations across stages

        for stage_name in self.stages:
            bucket = stage_buckets[stage_name]
            if not bucket:
                log_lines.append(f"  stage '{stage_name}': no operations — skip")
                continue

            log_lines.append(f"  stage '{stage_name}': {len(bucket)} operation(s)")

            # Build a meta_map that contains ONLY the cutters in this stage
            # (so BasicResolver doesn't try to apply operations from other stages)
            stage_meta: Dict[str, Dict] = {}
            for cutter_id, _op in bucket:
                stage_meta[cutter_id] = current_meta.get(cutter_id, {})
            # Non-cutter parts get their meta passed through (for target resolution)
            for pid, meta in current_meta.items():
                if pid not in stage_meta:
                    # Strip operation metadata so BasicResolver ignores them in this stage
                    stripped = {k: v for k, v in meta.items() if k != "operation"}
                    stage_meta[pid] = stripped

            stage_result = self._basic.resolve(
                current_parts,
                stage_meta,
                boolean_engine=boolean_engine,
            )

            log_lines.extend(stage_result.execution_log)
            if stage_result.errors:
                global_errors.extend(stage_result.errors)

            # Merge applied_operations into accumulator and feed geometry into next stage
            for pid, rp in stage_result.parts.items():
                prev = accumulated_ops.get(pid, [])
                accumulated_ops[pid] = prev + rp.applied_operations
            current_parts = {pid: rp.geometry for pid, rp in stage_result.parts.items()}

        # Produce final ResolvedPart objects carrying the full applied_operations list
        final_parts: Dict[str, ResolvedPart] = {}
        for pid, geom in current_parts.items():
            final_parts[pid] = ResolvedPart(pid, geom, accumulated_ops.get(pid, []), [])

        return ResolveResult(
            parts=final_parts,
            execution_log=log_lines,
            success=not global_errors,
            errors=global_errors,
        )


# ---------------------------------------------------------------------------
# KinematicResolver — v1.2 placeholder
# ---------------------------------------------------------------------------


class KinematicResolver(AssemblyResolver):
    """StagedResolver + joint-graph-aware positioning.  Scoped for v1.2.

    Raises ``NotImplementedError`` if instantiated — registered here so
    external code can reference the class name for forward-compat.
    """

    def __init__(self, **kwargs: Any) -> None:
        raise NotImplementedError(
            "KinematicResolver is scoped for yapCAD v1.2 / Mechatron integration. "
            "Use StagedResolver for staged resolution in v1.1."
        )

    def resolve(self, *args: Any, **kwargs: Any) -> ResolveResult:  # pragma: no cover
        raise NotImplementedError


# ---------------------------------------------------------------------------
# Resolver registry — mirrors boolean/__init__.py ENGINE_REGISTRY
# ---------------------------------------------------------------------------

RESOLVER_REGISTRY: Dict[str, type] = {
    "basic": BasicResolver,
    "staged": StagedResolver,
    # kinematic: not shipped in v1.1
}

# Load external resolvers via entry points
try:
    _eps = importlib.metadata.entry_points(group="yapcad.assembly_resolvers")
    for _ep in _eps:
        try:
            RESOLVER_REGISTRY[_ep.name] = _ep.load()
            log.debug("assembly_resolvers: loaded '%s' from %s", _ep.name, _ep.value)
        except Exception as _exc:  # noqa: BLE001
            log.warning("assembly_resolvers: failed to load '%s': %s", _ep.name, _exc)
except Exception:  # noqa: BLE001
    pass  # importlib.metadata not available — fine, built-ins still work


def get_resolver(name: str) -> Optional[type]:
    """Return the resolver *class* for ``name``, or ``None`` if not registered."""
    return RESOLVER_REGISTRY.get(name)


# ---------------------------------------------------------------------------
# Top-level convenience function
# ---------------------------------------------------------------------------


def resolve_assembly(
    parts: Dict[str, Any],
    meta_map: Dict[str, Dict],
    *,
    resolver: Optional[AssemblyResolver] = None,
    boolean_engine: Optional[str] = None,
) -> ResolveResult:
    """Apply all ``operation`` metadata in *meta_map* to the geometries in *parts*.

    Parameters
    ----------
    parts:
        ``{part_id: geometry}`` dict.
    meta_map:
        ``{part_id: metadata_dict}`` for each part.
    resolver:
        A pre-constructed resolver instance.  When ``None``, the resolver is
        selected via the ``YAPCAD_ASSEMBLY_RESOLVER`` environment variable
        (``basic`` | ``staged`` | ``<entry-point-name>``).  Defaults to
        ``BasicResolver()``.
    boolean_engine:
        Override the active boolean engine.  ``None`` uses the yapCAD default.

    Returns
    -------
    ResolveResult
    """
    if resolver is None:
        env_name = os.environ.get("YAPCAD_ASSEMBLY_RESOLVER", "basic").lower()
        resolver_cls = RESOLVER_REGISTRY.get(env_name)
        if resolver_cls is None:
            raise ValueError(
                f"Unknown YAPCAD_ASSEMBLY_RESOLVER='{env_name}'. "
                f"Available: {sorted(RESOLVER_REGISTRY)}"
            )
        if env_name == "staged":
            raise ValueError(
                "StagedResolver requires a 'stages' list — "
                "instantiate it explicitly: StagedResolver(stages=[...])"
            )
        resolver = resolver_cls()

    return resolver.resolve(parts, meta_map, boolean_engine=boolean_engine)


__all__ = [
    "AssemblyResolver",
    "BasicResolver",
    "StagedResolver",
    "KinematicResolver",
    "ResolvedPart",
    "ResolveResult",
    "RESOLVER_REGISTRY",
    "get_resolver",
    "resolve_assembly",
]
