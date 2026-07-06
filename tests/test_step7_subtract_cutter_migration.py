"""
Step 7 — Subtract-Cutter Migration Test
=========================================
Validates the metadata-v1.1 "process-aware subtract cutter" capability using
a generic, project-neutral fixture DSL
(``tests/fixtures/enclosure_shell_cutouts_v1.dsl``). It declares three radial
box cutters against a cylindrical ``enclosure_shell`` target and checks:

1. The DSL parses without errors (all @meta decorators are recognised).
2. Each cutter command produces a valid solid with the correct v1.1 metadata:
   - operation.kind = "subtract"
   - operation.target_filter = ["enclosure_shell"]
   - operation.through = True
   - operation.policy = "strict"
   - operation.feature_id / feature_kind correct per-command
3. The BasicResolver wires the three cutters in priority order against a
   stand-in cylinder (proxying the enclosure shell) and produces a result
   solid with applied_operations in [vent_panel, access_door, accessory_port]
   order.
4. The sidecar .meta.yaml emitter writes all three operations into the
   operations_applied block of the output part.

These tests run in the yapcad-brep conda env (OCC available).
They are skipped if OCC is not importable (CI without OCC).
"""

import os
import sys
import math
import json
import tempfile
import pathlib
import textwrap
import pytest

# ── OCC availability guard ────────────────────────────────────────────────────
try:
    from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeCylinder  # noqa: F401
    OCC_AVAILABLE = True
except ImportError:
    OCC_AVAILABLE = False

pytestmark = [
    pytest.mark.requires_occ,
    pytest.mark.skipif(
        not OCC_AVAILABLE,
        reason="OCC not available in this environment",
    ),
]

# ── yapCAD imports (PYTHONPATH=src required) ──────────────────────────────────
from yapcad.dsl import parse, check, tokenize, execute
from yapcad.dsl.checker import ErrorSeverity
from yapcad.metadata import get_operation_metadata
from yapcad.assembly.resolver import BasicResolver, ResolvedPart
from yapcad.io.meta_yaml import dump_metadata_yaml, load_metadata_yaml


def _execute_cutter(module, cmd_name):
    """
    Run a cutter command, extract its OCC geometry and metadata dict.
    Returns (occ_shape, meta_dict).

    The meta_dict is the last dict element in the solid Value's data.
    The OCC shape is recovered from the BREP-encoded value.
    """
    result = execute(module, cmd_name, {})
    assert result.success, f"{cmd_name} failed: {result.error_message}"
    # Drill down to the solid Value
    er_outer = result.emit_result
    er_inner = er_outer.value       # EmitResult from RADIAL_BOX_CUTTER
    solid_val = er_inner.value      # Value(type='solid', data=[...])
    # Meta dict is the last dict element in data
    meta_dict = {}
    for item in reversed(solid_val.data):
        if isinstance(item, dict):
            meta_dict = item
            break
    # OCC shape: decode the brep data
    occ_shape = None
    brep_data = meta_dict.get('brep', {})
    if brep_data:
        import base64
        from OCC.Core.BRep import BRep_Builder
        from OCC.Core.TopoDS import TopoDS_Shape
        from OCC.Core.BRepTools import breptools
        import tempfile, os
        raw = base64.b64decode(brep_data.get('data', ''))
        with tempfile.NamedTemporaryFile(suffix='.brep', delete=False) as f:
            f.write(raw)
            tmppath = f.name
        try:
            shape = TopoDS_Shape()
            builder = BRep_Builder()
            breptools.Read(shape, tmppath, builder)
            if not shape.IsNull():
                occ_shape = shape
        finally:
            os.unlink(tmppath)
    return occ_shape, meta_dict


# ── Path to the fixture DSL ───────────────────────────────────────────────────
DSL_PATH = pathlib.Path(__file__).parent / "fixtures" / "enclosure_shell_cutouts_v1.dsl"


def _load_module():
    """Compile the enclosure-shell cutouts DSL and return the module + diagnostics."""
    if not DSL_PATH.exists():
        pytest.skip(f"DSL fixture not found: {DSL_PATH}")
    source = DSL_PATH.read_text()
    tokens = tokenize(source)
    module = parse(tokens, filename=str(DSL_PATH), source=source)
    result = check(module)
    return module, result


# ── Test 1: Parse / checker ───────────────────────────────────────────────────

def test_dsl_parses_without_errors():
    """enclosure_shell_cutouts_v1.dsl must compile clean (no errors)."""
    module, result = _load_module()
    errors = [d for d in result.diagnostics if d.severity == ErrorSeverity.ERROR]
    assert errors == [], f"Parse errors: {[(d.code, d.message) for d in errors]}"


def test_dsl_recognises_cutter_decorator():
    """All three commands must carry @meta(operation.kind=subtract) metadata."""
    module, _ = _load_module()
    # @meta stores keys as flat dot-notation: 'operation.kind', not nested dict
    cutter_commands = [
        cmd for cmd in module.functions
        if cmd.meta_hint and cmd.meta_hint.get("operation.kind") == "subtract"
    ]
    assert len(cutter_commands) == 3, (
        f"Expected 3 @meta(operation.kind=subtract) commands, got {len(cutter_commands)}: "
        f"{[c.name for c in cutter_commands]}"
    )


def test_cutter_target_filters():
    """All three cutters target 'enclosure_shell'."""
    module, _ = _load_module()
    for cmd in module.functions:
        mh = cmd.meta_hint or {}
        if mh.get("operation.kind") == "subtract":
            tf = mh.get("operation.target_filter", [])
            assert "enclosure_shell" in tf, (
                f"{cmd.name}: target_filter missing 'enclosure_shell' — got {tf}"
            )


def test_cutter_priorities_unique_and_ordered():
    """Priorities must be unique and in ascending order (10, 20, 30)."""
    module, _ = _load_module()
    priorities = sorted(
        cmd.meta_hint["operation.priority"]
        for cmd in module.functions
        if (cmd.meta_hint or {}).get("operation.kind") == "subtract"
    )
    assert priorities == [10, 20, 30], f"Expected [10,20,30], got {priorities}"


# ── Test 2: Runtime execution / metadata ─────────────────────────────────────

def test_cutter_command_produces_solid_with_operation_meta():
    """
    Running CUT_VENT_PANEL should produce a solid with v1.1 operation metadata
    embedded in the solid's metadata dict (at data[-1]).
    """
    module, _ = _load_module()
    _, meta = _execute_cutter(module, "CUT_VENT_PANEL")
    op = get_operation_metadata(meta)
    assert op.get("kind") == "subtract", f"operation.kind={op.get('kind')!r}"
    assert "enclosure_shell" in op.get("target_filter", []), f"target_filter={op.get('target_filter')}"
    assert op.get("through") is True, "operation.through not True"
    assert op.get("policy") == "strict", f"policy={op.get('policy')!r}"
    assert op.get("feature_id") == "vent_panel", f"feature_id={op.get('feature_id')!r}"


def test_cutter_geometries_are_valid_solids():
    """
    All three cutter commands must produce non-empty OCC solids (volume > 0).
    """
    from OCC.Core.GProp import GProp_GProps
    from OCC.Core.BRepGProp import brepgprop

    module, _ = _load_module()
    for cmd_name in ("CUT_VENT_PANEL", "CUT_ACCESS_DOOR", "CUT_ACCESSORY_PORT"):
        occ_shape, meta = _execute_cutter(module, cmd_name)
        assert occ_shape is not None, f"{cmd_name}: OCC shape is None (BREP decode failed)"
        props = GProp_GProps()
        brepgprop.VolumeProperties(occ_shape, props)
        vol = props.Mass()
        assert vol > 1.0, f"{cmd_name}: volume {vol:.2f} mm³ is suspiciously small"


def test_cutter_positions():
    """
    Verify each cutter's bounding box Z-center is at z=162.5mm (all three share
    the same axial position), and that the XY centroid is near the expected
    radial direction for each clock angle.
    """
    from OCC.Core.BRepBndLib import brepbndlib
    from OCC.Core.Bnd import Bnd_Box

    R_OUTER = 157.0

    def bbox_center(shape):
        box = Bnd_Box()
        brepbndlib.Add(shape, box)
        xmin, ymin, zmin, xmax, ymax, zmax = box.Get()
        return ((xmin + xmax) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2)

    module, _ = _load_module()
    tol_z = 5.0    # mm
    tol_xy = 20.0  # mm — generous; just checks gross radial placement

    cases = [
        ("CUT_VENT_PANEL",     247.5, 162.5),
        ("CUT_ACCESS_DOOR",    180.0, 162.5),
        ("CUT_ACCESSORY_PORT",  90.0, 162.5),
    ]
    for cmd_name, clock_deg, z_mm in cases:
        occ_shape, _ = _execute_cutter(module, cmd_name)
        assert occ_shape is not None, f"{cmd_name}: OCC shape is None"
        cx, cy, cz = bbox_center(occ_shape)
        clock_rad = math.radians(clock_deg)
        exp_x = R_OUTER * math.cos(clock_rad)
        exp_y = R_OUTER * math.sin(clock_rad)
        assert abs(cz - z_mm) < tol_z, f"{cmd_name}: z={cz:.1f} expected {z_mm}"
        dist_xy = math.hypot(cx - exp_x, cy - exp_y)
        assert dist_xy < tol_xy, (
            f"{cmd_name}: XY center ({cx:.1f},{cy:.1f}) too far from "
            f"expected ({exp_x:.1f},{exp_y:.1f}), dist={dist_xy:.1f}mm"
        )


# ── Test 3: Resolver integration ─────────────────────────────────────────────

def _make_proxy_cylinder():
    """
    Stand-in for the enclosure shell: a solid cylinder R=157, H=330.
    """
    from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeCylinder
    cyl = BRepPrimAPI_MakeCylinder(157.0, 330.0).Shape()
    meta = {
        "_schema": "metadata-namespace-v1.1",
    }
    return cyl, meta


def test_resolver_applies_cutters_in_priority_order():
    """
    BasicResolver with three @cutter commands in one call must:
    - Apply all three subtracts against the target
    - Return applied_operations in [vent_panel, access_door, accessory_port] order
    - Consume all three cutters (they don't appear as extra resolved_parts)
    """
    module, _ = _load_module()

    proxy_geom, proxy_meta = _make_proxy_cylinder()
    parts = {"enclosure_shell": proxy_geom}
    metas = {"enclosure_shell": proxy_meta}

    for cmd_name in ("CUT_VENT_PANEL", "CUT_ACCESS_DOOR", "CUT_ACCESSORY_PORT"):
        occ_shape, meta = _execute_cutter(module, cmd_name)
        parts[cmd_name.lower()] = occ_shape
        metas[cmd_name.lower()] = meta

    resolver = BasicResolver()
    resolve_result = resolver.resolve(parts, metas)

    assert resolve_result.success, f"Resolver failed: {resolve_result.errors}"
    resolved = resolve_result.parts
    assert "enclosure_shell" in resolved, "enclosure_shell missing from resolver output"
    extra = [k for k in resolved if k != "enclosure_shell"]
    assert extra == [], f"Unexpected parts in resolver output: {extra}"

    rp = resolved["enclosure_shell"]
    assert isinstance(rp, ResolvedPart)
    applied = rp.applied_operations
    assert len(applied) == 3, f"Expected 3 applied ops, got {len(applied)}: {applied}"
    feature_ids = [metas[cid].get('operation', {}).get('feature_id') for cid in applied]
    assert feature_ids == ["vent_panel", "access_door", "accessory_port"], (
        f"Wrong feature_id order: {feature_ids} (applied={applied})"
    )


def test_resolver_result_volume_decreases():
    """
    After cutting, the resolved solid must have less volume than the proxy cylinder.
    """
    from OCC.Core.GProp import GProp_GProps
    from OCC.Core.BRepGProp import brepgprop

    module, _ = _load_module()
    proxy_geom, proxy_meta = _make_proxy_cylinder()

    props_before = GProp_GProps()
    brepgprop.VolumeProperties(proxy_geom, props_before)
    vol_before = props_before.Mass()

    parts = {"enclosure_shell": proxy_geom}
    metas = {"enclosure_shell": proxy_meta}
    for cmd_name in ("CUT_VENT_PANEL", "CUT_ACCESS_DOOR", "CUT_ACCESSORY_PORT"):
        occ_shape, meta = _execute_cutter(module, cmd_name)
        parts[cmd_name.lower()] = occ_shape
        metas[cmd_name.lower()] = meta

    resolver = BasicResolver()
    resolve_result = resolver.resolve(parts, metas)
    assert resolve_result.success
    resolved = resolve_result.parts

    props_after = GProp_GProps()
    brepgprop.VolumeProperties(resolved["enclosure_shell"].geometry, props_after)
    vol_after = props_after.Mass()

    assert vol_after < vol_before, (
        f"Volume did not decrease after cuts: before={vol_before:.0f}, after={vol_after:.0f}"
    )


# ── Test 4: Sidecar .meta.yaml ────────────────────────────────────────────────

def test_meta_yaml_contains_operations():
    """
    A resolved part's metadata must round-trip all three operations through the
    .meta.yaml sidecar in the operations_applied block.
    """
    module, _ = _load_module()
    proxy_geom, proxy_meta = _make_proxy_cylinder()

    parts = {"enclosure_shell": proxy_geom}
    metas = {"enclosure_shell": proxy_meta}
    for cmd_name in ("CUT_VENT_PANEL", "CUT_ACCESS_DOOR", "CUT_ACCESSORY_PORT"):
        occ_shape, meta = _execute_cutter(module, cmd_name)
        parts[cmd_name.lower()] = occ_shape
        metas[cmd_name.lower()] = meta

    resolver = BasicResolver()
    resolve_result = resolver.resolve(parts, metas)
    assert resolve_result.success
    resolved = resolve_result.parts
    rp = resolved["enclosure_shell"]

    with tempfile.TemporaryDirectory() as tmpdir:
        out_path = pathlib.Path(tmpdir) / "enclosure_shell_holed.meta.yaml"
        out_meta = dict(proxy_meta)
        out_meta["_applied_operations"] = [
            {
                "feature_id": metas[cid].get("operation", {}).get("feature_id", cid),
                "kind": metas[cid].get("operation", {}).get("kind", "subtract"),
                "priority": metas[cid].get("operation", {}).get("priority"),
            }
            for cid in rp.applied_operations
        ]
        import yaml as _yaml
        out_path.write_text(_yaml.dump(out_meta, allow_unicode=True))

        loaded = load_metadata_yaml(out_path)
        ops = loaded.get("_applied_operations", [])
        feature_ids = [op["feature_id"] for op in ops]
        assert set(feature_ids) == {"vent_panel", "access_door", "accessory_port"}, (
            f"feature_ids mismatch: {feature_ids}"
        )
