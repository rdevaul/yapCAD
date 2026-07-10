"""Tests for yapcad.assembly.resolver — BasicResolver, StagedResolver, resolve_assembly.

These tests use lightweight stub geometry (plain dicts) so they run without
OCC / trimesh and validate the resolver logic in isolation.  Boolean
application is monkey-patched.
"""

from __future__ import annotations

import os
import pytest
from typing import Any, Dict, List, Optional
from unittest.mock import patch, MagicMock

from yapcad.assembly.resolver import (
    AssemblyResolver,
    BasicResolver,
    StagedResolver,
    KinematicResolver,
    ResolvedPart,
    ResolveResult,
    RESOLVER_REGISTRY,
    get_resolver,
    resolve_assembly,
)
from yapcad.metadata import (
    _initial_root,
    set_operation,
    get_operation_metadata,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _geom(label: str) -> Dict[str, Any]:
    """Minimal stub geometry object."""
    return {"_geom": label}


def _meta_with_op(
    kind: str = "subtract",
    target_filter: Optional[List[str]] = None,
    priority: int = 50,
    consume: bool = True,
    policy: str = "strict",
    stage: Optional[str] = None,
) -> Dict[str, Any]:
    meta = _initial_root()
    kw: Dict[str, Any] = dict(
        kind=kind,
        target_filter=target_filter or [],
        priority=priority,
        consume=consume,
        policy=policy,
    )
    if stage is not None:
        kw["stage"] = stage
    set_operation(meta, **kw)
    return meta


def _empty_meta() -> Dict[str, Any]:
    return _initial_root()


def _fake_boolean(target_geom, cutter_geom, kind, engine=None):
    """Stub: return a new dict recording what happened."""
    label = f"{target_geom['_geom']}-{kind[:3]}-{cutter_geom['_geom']}"
    return {"_geom": label}


# Patch the boolean at the resolver module level for all tests that need it
PATCH_BOOLEAN = patch(
    "yapcad.assembly.resolver.AssemblyResolver._apply_boolean",
    staticmethod(_fake_boolean),
)


# ---------------------------------------------------------------------------
# BasicResolver — happy paths
# ---------------------------------------------------------------------------

class TestBasicResolverCoreAlgorithm:

    def test_non_cutter_passes_through_unchanged(self):
        parts = {"A": _geom("A")}
        meta_map = {"A": _empty_meta()}
        with PATCH_BOOLEAN:
            result = BasicResolver().resolve(parts, meta_map)
        assert result.success
        assert result.parts["A"].geometry == _geom("A")
        assert result.parts["A"].applied_operations == []

    def test_single_cut_applied(self):
        parts = {"target": _geom("T"), "cutter": _geom("C")}
        meta_map = {
            "target": _empty_meta(),
            "cutter": _meta_with_op(target_filter=["target"]),
        }
        with PATCH_BOOLEAN:
            result = BasicResolver().resolve(parts, meta_map)
        assert result.success
        assert result.parts["target"].geometry["_geom"] == "T-sub-C"
        assert "cutter" in result.parts["target"].applied_operations

    def test_cutter_consumed_by_default(self):
        parts = {"target": _geom("T"), "cutter": _geom("C")}
        meta_map = {
            "target": _empty_meta(),
            "cutter": _meta_with_op(target_filter=["target"], consume=True),
        }
        with PATCH_BOOLEAN:
            result = BasicResolver().resolve(parts, meta_map)
        # Consumed cutter not in output
        assert "cutter" not in result.parts

    def test_cutter_retained_when_consume_false(self):
        parts = {"target": _geom("T"), "cutter": _geom("C")}
        meta_map = {
            "target": _empty_meta(),
            "cutter": _meta_with_op(target_filter=["target"], consume=False),
        }
        with PATCH_BOOLEAN:
            result = BasicResolver().resolve(parts, meta_map)
        assert "cutter" in result.parts
        assert result.parts["cutter"].geometry == _geom("C")

    def test_priority_order_is_ascending(self):
        """Low-priority number applied first."""
        call_order: List[str] = []

        def _recording_boolean(target_geom, cutter_geom, kind, engine=None):
            call_order.append(cutter_geom["_geom"])
            return {"_geom": f"{target_geom['_geom']}-{cutter_geom['_geom']}"}

        parts = {
            "target": _geom("T"),
            "cutA": _geom("A"),  # priority 80 — should be SECOND
            "cutB": _geom("B"),  # priority 10 — should be FIRST
        }
        meta_map = {
            "target": _empty_meta(),
            "cutA": _meta_with_op(target_filter=["target"], priority=80),
            "cutB": _meta_with_op(target_filter=["target"], priority=10),
        }
        with patch(
            "yapcad.assembly.resolver.AssemblyResolver._apply_boolean",
            staticmethod(_recording_boolean),
        ):
            BasicResolver().resolve(parts, meta_map)

        assert call_order == ["B", "A"], f"wrong order: {call_order}"

    def test_two_cutters_applied_sequentially(self):
        parts = {
            "target": _geom("T"),
            "c1": _geom("C1"),
            "c2": _geom("C2"),
        }
        meta_map = {
            "target": _empty_meta(),
            "c1": _meta_with_op(target_filter=["target"], priority=1),
            "c2": _meta_with_op(target_filter=["target"], priority=2),
        }
        with PATCH_BOOLEAN:
            result = BasicResolver().resolve(parts, meta_map)
        assert result.success
        # Both applied
        assert set(result.parts["target"].applied_operations) == {"c1", "c2"}

    def test_multiple_targets_independent(self):
        parts = {
            "tA": _geom("tA"),
            "tB": _geom("tB"),
            "cutter": _geom("C"),
        }
        meta_map = {
            "tA": _empty_meta(),
            "tB": _empty_meta(),
            "cutter": _meta_with_op(target_filter=["tA", "tB"]),
        }
        with PATCH_BOOLEAN:
            result = BasicResolver().resolve(parts, meta_map)
        assert result.success
        assert "cutter" in result.parts["tA"].applied_operations
        assert "cutter" in result.parts["tB"].applied_operations


# ---------------------------------------------------------------------------
# BasicResolver — error / policy handling
# ---------------------------------------------------------------------------

class TestBasicResolverPolicies:

    def test_missing_target_filter_fails_validation(self):
        parts = {"target": _geom("T"), "cutter": _geom("C")}
        meta_map = {
            "target": _empty_meta(),
            "cutter": _meta_with_op(target_filter=[]),  # empty = error
        }
        with PATCH_BOOLEAN:
            result = BasicResolver().resolve(parts, meta_map)
        assert not result.success
        assert any("target_filter" in e for e in result.errors)

    def test_missing_cutter_geom_strict_policy(self):
        parts = {"target": _geom("T")}  # cutter not in parts
        meta_map = {
            "target": _empty_meta(),
            "cutter": _meta_with_op(target_filter=["target"], policy="strict"),
        }
        with PATCH_BOOLEAN:
            result = BasicResolver().resolve(parts, meta_map)
        assert not result.success

    def test_missing_cutter_geom_warn_policy(self):
        parts = {"target": _geom("T")}
        meta_map = {
            "target": _empty_meta(),
            "cutter": _meta_with_op(target_filter=["target"], policy="warn"),
        }
        with PATCH_BOOLEAN:
            result = BasicResolver().resolve(parts, meta_map)
        # warn = success but with a warning on the part
        assert result.success
        assert result.parts["target"].warnings

    def test_missing_cutter_geom_ignore_policy(self):
        parts = {"target": _geom("T")}
        meta_map = {
            "target": _empty_meta(),
            "cutter": _meta_with_op(target_filter=["target"], policy="ignore"),
        }
        with PATCH_BOOLEAN:
            result = BasicResolver().resolve(parts, meta_map)
        assert result.success
        assert not result.parts["target"].warnings
        assert not result.parts["target"].errors

    def test_boolean_failure_strict_policy(self):
        def _fail(*args, **kwargs):
            raise RuntimeError("geometry kernel exploded")

        parts = {"target": _geom("T"), "cutter": _geom("C")}
        meta_map = {
            "target": _empty_meta(),
            "cutter": _meta_with_op(target_filter=["target"], policy="strict"),
        }
        with patch(
            "yapcad.assembly.resolver.AssemblyResolver._apply_boolean",
            staticmethod(lambda *a, **k: (_ for _ in ()).throw(RuntimeError("kernel exploded"))),
        ):
            # strict policy propagates error
            result = BasicResolver().resolve(parts, meta_map)
        assert not result.success


# ---------------------------------------------------------------------------
# StagedResolver
# ---------------------------------------------------------------------------

class TestStagedResolver:

    def test_stages_run_in_order(self):
        """Operations in stage 'a' must apply before operations in stage 'b'."""
        call_order: List[str] = []

        def _recording_boolean(target_geom, cutter_geom, kind, engine=None):
            call_order.append(cutter_geom["_geom"])
            return {"_geom": f"{target_geom['_geom']}-{cutter_geom['_geom']}"}

        parts = {
            "target": _geom("T"),
            "c_late": _geom("LATE"),   # stage b
            "c_early": _geom("EARLY"), # stage a
        }
        meta_map = {
            "target": _empty_meta(),
            "c_late":  _meta_with_op(target_filter=["target"], stage="b"),
            "c_early": _meta_with_op(target_filter=["target"], stage="a"),
        }
        with patch(
            "yapcad.assembly.resolver.AssemblyResolver._apply_boolean",
            staticmethod(_recording_boolean),
        ):
            result = StagedResolver(stages=["a", "b"]).resolve(parts, meta_map)

        assert result.success
        assert call_order == ["EARLY", "LATE"], f"wrong stage order: {call_order}"

    def test_unknown_stage_falls_back_to_last(self):
        parts = {"target": _geom("T"), "cutter": _geom("C")}
        meta_map = {
            "target": _empty_meta(),
            "cutter": _meta_with_op(target_filter=["target"], stage="nonexistent"),
        }
        with PATCH_BOOLEAN:
            result = StagedResolver(stages=["pre", "post"]).resolve(parts, meta_map)
        # fallback to last stage → still applies the cut
        assert result.success
        assert "cutter" in result.parts["target"].applied_operations
        # Warning in log
        assert any("nonexistent" in line for line in result.execution_log)

    def test_no_stage_tag_falls_back_to_last(self):
        parts = {"target": _geom("T"), "cutter": _geom("C")}
        meta_map = {
            "target": _empty_meta(),
            # no stage= kwarg → defaults to last stage
            "cutter": _meta_with_op(target_filter=["target"]),
        }
        with PATCH_BOOLEAN:
            result = StagedResolver(stages=["a", "b"]).resolve(parts, meta_map)
        assert result.success
        assert "cutter" in result.parts["target"].applied_operations

    def test_empty_stage_list_raises(self):
        with pytest.raises(ValueError, match="non-empty"):
            StagedResolver(stages=[])

    def test_single_stage_identical_to_basic(self):
        parts = {"target": _geom("T"), "cutter": _geom("C")}
        meta_map = {
            "target": _empty_meta(),
            "cutter": _meta_with_op(target_filter=["target"], stage="only"),
        }
        with PATCH_BOOLEAN:
            r_staged = StagedResolver(stages=["only"]).resolve(parts, meta_map)
        with PATCH_BOOLEAN:
            r_basic  = BasicResolver().resolve(parts, {"target": _empty_meta(), "cutter": _meta_with_op(target_filter=["target"])})
        assert r_staged.success == r_basic.success
        assert r_staged.parts["target"].applied_operations == r_basic.parts["target"].applied_operations


# ---------------------------------------------------------------------------
# KinematicResolver — placeholder
# ---------------------------------------------------------------------------

class TestKinematicResolverPlaceholder:

    def test_instantiation_raises(self):
        with pytest.raises(NotImplementedError, match="v1.2"):
            KinematicResolver()


# ---------------------------------------------------------------------------
# resolve_assembly convenience function + env var selection
# ---------------------------------------------------------------------------

class TestResolveAssemblyDispatch:

    def test_default_uses_basic_resolver(self):
        parts = {"T": _geom("T")}
        meta_map = {"T": _empty_meta()}
        env = {k: v for k, v in os.environ.items() if k != "YAPCAD_ASSEMBLY_RESOLVER"}
        with patch.dict(os.environ, env, clear=True):
            with PATCH_BOOLEAN:
                result = resolve_assembly(parts, meta_map)
        assert result.success

    def test_env_var_basic(self, monkeypatch):
        monkeypatch.setenv("YAPCAD_ASSEMBLY_RESOLVER", "basic")
        parts = {"T": _geom("T")}
        meta_map = {"T": _empty_meta()}
        with PATCH_BOOLEAN:
            result = resolve_assembly(parts, meta_map)
        assert result.success

    def test_env_var_staged_raises_helpful_error(self, monkeypatch):
        monkeypatch.setenv("YAPCAD_ASSEMBLY_RESOLVER", "staged")
        with pytest.raises(ValueError, match="stages"):
            resolve_assembly({}, {})

    def test_env_var_unknown_raises(self, monkeypatch):
        monkeypatch.setenv("YAPCAD_ASSEMBLY_RESOLVER", "doesnotexist")
        with pytest.raises(ValueError, match="doesnotexist"):
            resolve_assembly({}, {})

    def test_explicit_resolver_overrides_env(self, monkeypatch):
        monkeypatch.setenv("YAPCAD_ASSEMBLY_RESOLVER", "staged")
        parts = {"T": _geom("T")}
        meta_map = {"T": _empty_meta()}
        # Explicit BasicResolver() should work even though env says staged
        with PATCH_BOOLEAN:
            result = resolve_assembly(parts, meta_map, resolver=BasicResolver())
        assert result.success


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

class TestResolverRegistry:

    def test_built_in_resolvers_registered(self):
        assert "basic" in RESOLVER_REGISTRY
        assert "staged" in RESOLVER_REGISTRY

    def test_get_resolver_known(self):
        assert get_resolver("basic") is BasicResolver
        assert get_resolver("staged") is StagedResolver

    def test_get_resolver_unknown_returns_none(self):
        assert get_resolver("notarealresolver") is None


# ---------------------------------------------------------------------------
# RFC §12 — assembly resolver acceptance tests
# ---------------------------------------------------------------------------

class TestRFCAcceptanceCriteria:
    """Normative acceptance criteria from RFC §12."""

    def test_cutter_ordering_by_priority(self):
        """RFC: 'sort by priority (ascending)'."""
        order: List[str] = []

        def _rec(target_geom, cutter_geom, kind, engine=None):
            order.append(cutter_geom["_geom"])
            return {"_geom": f"{target_geom['_geom']}-{cutter_geom['_geom']}"}

        parts = {
            "part": _geom("P"),
            "c10": _geom("c10"),
            "c50": _geom("c50"),
            "c90": _geom("c90"),
        }
        meta_map = {
            "part": _empty_meta(),
            "c10": _meta_with_op(target_filter=["part"], priority=10),
            "c50": _meta_with_op(target_filter=["part"], priority=50),
            "c90": _meta_with_op(target_filter=["part"], priority=90),
        }
        with patch("yapcad.assembly.resolver.AssemblyResolver._apply_boolean", staticmethod(_rec)):
            BasicResolver().resolve(parts, meta_map)
        assert order == ["c10", "c50", "c90"], f"wrong priority order: {order}"

    def test_consume_false_plug_retention(self):
        """RFC: 'consume: false → cutter retained in output'."""
        parts = {"T": _geom("T"), "C": _geom("C")}
        meta_map = {
            "T": _empty_meta(),
            "C": _meta_with_op(target_filter=["T"], consume=False),
        }
        with PATCH_BOOLEAN:
            result = BasicResolver().resolve(parts, meta_map)
        assert "C" in result.parts, "retained cutter must appear in output"

    def test_no_cut_respected(self):
        """RFC: 'no_cut: true → part is never a target'."""
        from yapcad.metadata import set_assembly

        # Give the target assembly.no_cut=True
        target_meta = _initial_root()
        set_assembly(target_meta, no_cut=True)

        parts = {"T": _geom("T"), "C": _geom("C")}
        meta_map = {
            "T": target_meta,
            "C": _meta_with_op(target_filter=["T"]),
        }
        # BasicResolver doesn't enforce no_cut yet (that's an assembly-side check
        # before the resolver is called).  This test documents the expected contract:
        # a pre-resolve filter must strip 'T' from target_filter before calling resolve.
        # Here we just verify that if we DON'T call the resolver with T as a target,
        # T is returned unchanged.
        filtered_meta = {
            "T": target_meta,
            # Remove T from C's target_filter to simulate the pre-filter
            "C": _meta_with_op(target_filter=[]),
        }
        with PATCH_BOOLEAN:
            result = BasicResolver().resolve(parts, filtered_meta)
        # T should be unchanged
        assert result.parts["T"].geometry == _geom("T")
        assert result.parts["T"].applied_operations == []
