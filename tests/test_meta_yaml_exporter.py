"""Tests for the .meta.yaml sidecar emitter (Step 5 — metadata-namespace-v1.1).

Covers:
- dump_metadata_yaml round-trip (write → load → assert equality)
- sidecar_path_for path computation
- Parts with no metadata produce no sidecar (no-op)
- _provenance block is present and contains required keys
- load_metadata_yaml round-trip
- v1.0 sidecar loads cleanly as v1.1 (backwards compat)
- write_step_with_meta / write_stl_with_meta emit the sidecar
- explicit yaml_path override
- source_dsl sha256 provenance
"""

from __future__ import annotations

import hashlib
import os
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from yapcad.io.meta_yaml import dump_metadata_yaml, load_metadata_yaml, sidecar_path_for
from yapcad.metadata import (
    set_material,
    add_tags,
    get_assembly_metadata,
    get_operation_metadata,
    add_bolt_pattern,
    set_operation,
    _initial_root,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_part_with_meta(**kwargs) -> MagicMock:
    """Return a mock part that carries a minimal v1.1 metadata dict."""
    part = MagicMock()
    meta_dict = _initial_root()
    meta_dict.update(kwargs)
    part._metadata = meta_dict
    return part


def _part_with_real_meta() -> MagicMock:
    """Build a mock whose metadata dict has assembly + operation namespaces."""
    part = MagicMock()
    meta = {
        "schema": "metadata-namespace-v1.1",
        "entityId": "abc-123",
        "tags": ["agentic-1"],
        "layer": "default",
        "material": {"name": "PETG-CF", "grade": "DML-spec"},
        "assembly": {
            "joint_kind": "radial",
            "bolt_patterns": [
                {
                    "id": "bot_ring",
                    "ring": "radial",
                    "R_mm": 155.5,
                    "z_mm": 10.0,
                    "count": 8,
                    "fastener": {"designation": "M8x1.25", "thread": "heat-set", "head": "SHCS"},
                }
            ],
        },
        "operation": {
            "kind": "subtract",
            "target_filter": ["forward_bulkhead"],
            "priority": 10,
            "through": True,
            "consume": False,
            "policy": "strict",
        },
    }
    part._metadata = meta
    return part


# ---------------------------------------------------------------------------
# sidecar_path_for
# ---------------------------------------------------------------------------

class TestSidecarPathFor:
    def test_step_extension(self):
        p = sidecar_path_for("/some/dir/my_part.step")
        assert p.name == "my_part.meta.yaml"
        assert str(p.parent) == "/some/dir"

    def test_stl_extension(self):
        p = sidecar_path_for("/out/rocket_CURRENT.stl")
        assert p.name == "rocket_CURRENT.meta.yaml"

    def test_multi_suffix(self):
        # Agentic-1 convention uses underscores, not dots: boattail_CURRENT.step
        # All dot-suffixes are stripped; "_CURRENT" is part of the stem.
        p = sidecar_path_for("/out/boattail_CURRENT.step")
        assert p.name == "boattail_CURRENT.meta.yaml"

    def test_relative_path(self):
        p = sidecar_path_for("output/part.step")
        assert p.name == "part.meta.yaml"


# ---------------------------------------------------------------------------
# dump_metadata_yaml — basic round-trip
# ---------------------------------------------------------------------------

class TestDumpMetadataYaml:
    def test_no_metadata_returns_none(self, tmp_path):
        """A part with no metadata dict should not write a sidecar."""
        part = MagicMock()
        # Ensure _extract_metadata finds nothing.
        part._metadata = None
        part.metadata = None
        part.meta = None

        with patch("yapcad.io.meta_yaml._extract_metadata", return_value=None):
            geo = tmp_path / "part.step"
            geo.write_text("ISO-10303-21;")
            result = dump_metadata_yaml(part, geo)

        assert result is None
        assert not (tmp_path / "part.meta.yaml").exists()

    def test_basic_round_trip(self, tmp_path):
        """Write sidecar, reload it, assert all top-level keys survive."""
        geo = tmp_path / "my_part.step"
        geo.write_text("ISO-10303-21;")

        meta = {
            "schema": "metadata-namespace-v1.1",
            "entityId": "eid-001",
            "tags": ["test", "round-trip"],
            "layer": "default",
            "material": {"name": "Al6061", "grade": "T6"},
        }

        with patch("yapcad.io.meta_yaml._extract_metadata", return_value=meta):
            part = MagicMock()
            out = dump_metadata_yaml(part, geo)

        assert out is not None
        sidecar = Path(out)
        assert sidecar.exists()
        assert sidecar.name == "my_part.meta.yaml"

        loaded = load_metadata_yaml(sidecar)
        assert loaded["schema"] == "metadata-namespace-v1.1"
        assert loaded["entityId"] == "eid-001"
        assert loaded["tags"] == ["test", "round-trip"]
        assert loaded["material"]["name"] == "Al6061"
        assert "_provenance" in loaded

    def test_provenance_keys_present(self, tmp_path):
        geo = tmp_path / "part.step"
        geo.write_text("ISO-10303-21;")

        meta = {"schema": "metadata-namespace-v1.1", "entityId": "x"}

        with patch("yapcad.io.meta_yaml._extract_metadata", return_value=meta):
            out = dump_metadata_yaml(MagicMock(), geo)

        loaded = load_metadata_yaml(out)
        prov = loaded["_provenance"]
        assert "yapcad_version" in prov
        assert "exported_at" in prov
        # exported_at should end with Z (UTC)
        assert prov["exported_at"].endswith("Z")

    def test_explicit_yaml_path(self, tmp_path):
        geo = tmp_path / "part.step"
        geo.write_text("ISO-10303-21;")
        explicit = tmp_path / "custom_sidecar.yaml"

        meta = {"schema": "metadata-namespace-v1.1", "entityId": "y"}

        with patch("yapcad.io.meta_yaml._extract_metadata", return_value=meta):
            out = dump_metadata_yaml(MagicMock(), geo, yaml_path=explicit)

        assert Path(out).resolve() == explicit.resolve()
        assert explicit.exists()

    def test_source_dsl_provenance(self, tmp_path):
        geo = tmp_path / "part.step"
        geo.write_text("ISO-10303-21;")

        dsl_file = tmp_path / "forward_bulkhead.dsl"
        dsl_content = b"# fake DSL content for sha256 test\nFORWARD_BULKHEAD()\n"
        dsl_file.write_bytes(dsl_content)
        expected_sha = hashlib.sha256(dsl_content).hexdigest()

        meta = {"schema": "metadata-namespace-v1.1", "entityId": "z"}

        with patch("yapcad.io.meta_yaml._extract_metadata", return_value=meta):
            out = dump_metadata_yaml(
                MagicMock(),
                geo,
                source_dsl=str(dsl_file),
                source_command="FORWARD_BULKHEAD",
            )

        loaded = load_metadata_yaml(out)
        prov = loaded["_provenance"]
        assert prov["source_dsl"] == str(dsl_file)
        assert prov["source_dsl_sha256"] == expected_sha
        assert prov["source_command"] == "FORWARD_BULKHEAD"

    def test_full_assembly_operation_round_trip(self, tmp_path):
        """Full namespace round-trip with assembly + operation sections."""
        geo = tmp_path / "boattail.step"
        geo.write_text("ISO-10303-21;")

        meta = {
            "schema": "metadata-namespace-v1.1",
            "entityId": "eid-bt-001",
            "tags": ["agentic-1", "boattail"],
            "layer": "default",
            "assembly": {
                "joint_kind": "radial",
                "bolt_patterns": [
                    {
                        "id": "ring_a",
                        "ring": "radial",
                        "R_mm": 155.5,
                        "z_mm": 10.0,
                        "count": 8,
                    }
                ],
            },
            "operation": {
                "kind": "subtract",
                "target_filter": ["boattail"],
                "priority": 10,
                "through": True,
                "consume": False,
                "policy": "strict",
            },
        }

        with patch("yapcad.io.meta_yaml._extract_metadata", return_value=meta):
            out = dump_metadata_yaml(MagicMock(), geo)

        loaded = load_metadata_yaml(out)

        assert loaded["assembly"]["joint_kind"] == "radial"
        assert loaded["assembly"]["bolt_patterns"][0]["R_mm"] == 155.5
        assert loaded["operation"]["kind"] == "subtract"
        assert loaded["operation"]["priority"] == 10
        assert loaded["operation"]["target_filter"] == ["boattail"]


# ---------------------------------------------------------------------------
# Backwards compatibility: v1.0 sidecar loads cleanly
# ---------------------------------------------------------------------------

class TestBackwardsCompat:
    def test_v1_0_sidecar_loads(self, tmp_path):
        """A v1.0 metadata dict (no assembly/operation keys) should load fine."""
        sidecar = tmp_path / "legacy.meta.yaml"
        import yaml
        legacy = {
            "schema": "metadata-namespace-v1.0",
            "entityId": "legacy-001",
            "tags": ["legacy"],
            "layer": "default",
            "material": {"name": "PLA"},
        }
        sidecar.write_text(yaml.dump(legacy))

        loaded = load_metadata_yaml(sidecar)
        assert loaded["schema"] == "metadata-namespace-v1.0"
        assert "assembly" not in loaded
        assert "operation" not in loaded
        assert loaded["material"]["name"] == "PLA"


# ---------------------------------------------------------------------------
# write_step_with_meta / write_stl_with_meta integration
# ---------------------------------------------------------------------------

class TestWithMetaExporters:
    def test_write_step_with_meta_emits_sidecar(self, tmp_path):
        from yapcad.io.step import write_step_with_meta

        geo = tmp_path / "out.step"
        meta = {"schema": "metadata-namespace-v1.1", "entityId": "step-001"}
        part = MagicMock()

        with (
            patch("yapcad.io.step.write_step") as mock_write,
            patch("yapcad.io.meta_yaml._extract_metadata", return_value=meta),
        ):
            result = write_step_with_meta(part, str(geo), name="test_part")

        mock_write.assert_called_once_with(part, str(geo), name="test_part")
        assert result is not None
        assert result.endswith(".meta.yaml")

    def test_write_stl_with_meta_emits_sidecar(self, tmp_path):
        from yapcad.io.stl import write_stl_with_meta

        geo = tmp_path / "out.stl"
        meta = {"schema": "metadata-namespace-v1.1", "entityId": "stl-001"}
        part = MagicMock()

        with (
            patch("yapcad.io.stl.write_stl") as mock_write,
            patch("yapcad.io.meta_yaml._extract_metadata", return_value=meta),
        ):
            result = write_stl_with_meta(part, str(geo))

        mock_write.assert_called_once()
        assert result is not None
        assert result.endswith(".meta.yaml")

    def test_emit_sidecar_false_skips_yaml(self, tmp_path):
        from yapcad.io.step import write_step_with_meta

        geo = tmp_path / "out.step"
        part = MagicMock()

        with patch("yapcad.io.step.write_step"):
            result = write_step_with_meta(part, str(geo), emit_sidecar=False)

        assert result is None
        assert not (tmp_path / "out.meta.yaml").exists()

    def test_no_metadata_no_sidecar(self, tmp_path):
        from yapcad.io.step import write_step_with_meta

        geo = tmp_path / "out.step"
        part = MagicMock()

        with (
            patch("yapcad.io.step.write_step"),
            patch("yapcad.io.meta_yaml._extract_metadata", return_value=None),
        ):
            result = write_step_with_meta(part, str(geo))

        assert result is None
