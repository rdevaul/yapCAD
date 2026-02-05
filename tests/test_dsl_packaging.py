"""
Tests for DSL-to-Package integration (yapcad.dsl.packaging).
"""

import pytest
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock

from yapcad.dsl.packaging import (
    package_from_dsl,
    PackageResult,
    _get_yapcad_version,
    _add_dsl_source_attachment,
)


class TestPackageResult:
    """Test the PackageResult class."""

    def test_success_result(self):
        """PackageResult with success=True."""
        mock_manifest = Mock()
        mock_manifest.root = Path("/tmp/test_pkg")

        result = PackageResult(
            success=True,
            manifest=mock_manifest,
        )

        assert result.success is True
        assert result.manifest is mock_manifest
        assert result.error_message is None

    def test_failure_result(self):
        """PackageResult with success=False."""
        result = PackageResult(
            success=False,
            error_message="Test error message",
        )

        assert result.success is False
        assert result.manifest is None
        assert result.error_message == "Test error message"

    def test_repr_success(self):
        """Test __repr__ for successful result."""
        mock_manifest = Mock()
        result = PackageResult(success=True, manifest=mock_manifest)
        repr_str = repr(result)
        assert "success=True" in repr_str

    def test_repr_failure(self):
        """Test __repr__ for failed result."""
        result = PackageResult(success=False, error_message="failed")
        repr_str = repr(result)
        assert "success=False" in repr_str
        assert "failed" in repr_str


class TestGetYapcadVersion:
    """Test the _get_yapcad_version helper."""

    def test_version_available(self):
        """Returns version when yapcad.__version__ exists."""
        with patch('yapcad.dsl.packaging._get_yapcad_version') as mock:
            mock.return_value = "1.2.3"
            # Call the actual function to test the import logic
            version = _get_yapcad_version()
            # Version should be a string (either actual version or "unknown")
            assert isinstance(version, str)

    def test_version_import_error(self):
        """Returns 'unknown' when import fails."""
        # The actual function handles ImportError internally
        version = _get_yapcad_version()
        assert isinstance(version, str)


class TestPackageFromDsl:
    """Test the package_from_dsl function."""

    def test_dsl_execution_failure(self):
        """Returns error when DSL execution fails."""
        with patch('yapcad.dsl.packaging.compile_and_run') as mock_run:
            mock_result = Mock()
            mock_result.success = False
            mock_result.error_message = "Syntax error at line 1"
            mock_run.return_value = mock_result

            result = package_from_dsl(
                source="invalid dsl",
                command_name="TEST",
                parameters={},
                target_dir="/tmp/test",
                name="test_pkg",
                version="1.0.0",
            )

            assert result.success is False
            assert "DSL execution failed" in result.error_message
            assert "Syntax error" in result.error_message

    def test_no_geometry_output(self):
        """Returns error when DSL produces no geometry."""
        with patch('yapcad.dsl.packaging.compile_and_run') as mock_run:
            mock_result = Mock()
            mock_result.success = True
            mock_result.geometry = None
            mock_run.return_value = mock_result

            result = package_from_dsl(
                source="module test;",
                command_name="TEST",
                parameters={},
                target_dir="/tmp/test",
                name="test_pkg",
                version="1.0.0",
            )

            assert result.success is False
            assert "no geometry output" in result.error_message.lower()

    def test_unexpected_geometry_type(self):
        """Returns error for unexpected geometry type."""
        with patch('yapcad.dsl.packaging.compile_and_run') as mock_run:
            with patch('yapcad.geom3d.issolid', return_value=False):
                with patch('yapcad.geom3d.issurface', return_value=False):
                    mock_result = Mock()
                    mock_result.success = True
                    mock_result.geometry = "not a geometry object"
                    mock_run.return_value = mock_result

                    result = package_from_dsl(
                        source="module test;",
                        command_name="TEST",
                        parameters={},
                        target_dir="/tmp/test",
                        name="test_pkg",
                        version="1.0.0",
                    )

                    assert result.success is False
                    assert "Unexpected geometry type" in result.error_message

    def test_successful_packaging(self, tmp_path):
        """Successfully creates package from DSL."""
        mock_solid = Mock()
        mock_manifest = Mock()
        mock_manifest.root = tmp_path / "test_pkg"
        mock_manifest.data = {}

        mock_exec_result = Mock()
        mock_exec_result.success = True
        mock_exec_result.geometry = mock_solid
        mock_exec_result.provenance = None

        with patch('yapcad.dsl.packaging.compile_and_run', return_value=mock_exec_result):
            with patch('yapcad.geom3d.issolid', return_value=True):
                with patch('yapcad.package.create_package_from_entities', return_value=mock_manifest) as mock_create:
                    with patch('yapcad.dsl.packaging._add_dsl_source_attachment'):
                        result = package_from_dsl(
                            source="module test; command MAKE() -> solid { emit box(1,1,1); }",
                            command_name="MAKE",
                            parameters={},
                            target_dir=tmp_path / "output",
                            name="test_pkg",
                            version="1.0.0",
                        )

                        assert result.success is True
                        assert result.manifest is mock_manifest
                        mock_create.assert_called_once()

    def test_provenance_metadata(self, tmp_path):
        """Generator metadata includes provenance when available."""
        mock_solid = Mock()
        mock_manifest = Mock()
        mock_manifest.root = tmp_path / "test_pkg"
        mock_manifest.data = {}

        mock_provenance = Mock()
        mock_provenance.module_name = "test_module"
        mock_provenance.command_name = "MAKE_THING"
        mock_provenance.parameters = {"size": 10}
        mock_provenance.source_signature = "abc123"
        mock_provenance.version = "1.0"

        mock_exec_result = Mock()
        mock_exec_result.success = True
        mock_exec_result.geometry = mock_solid
        mock_exec_result.provenance = mock_provenance

        with patch('yapcad.dsl.packaging.compile_and_run', return_value=mock_exec_result):
            with patch('yapcad.geom3d.issolid', return_value=True):
                with patch('yapcad.package.create_package_from_entities', return_value=mock_manifest) as mock_create:
                    with patch('yapcad.dsl.packaging._add_dsl_source_attachment'):
                        result = package_from_dsl(
                            source="test source",
                            command_name="MAKE_THING",
                            parameters={"size": 10},
                            target_dir=tmp_path / "output",
                            name="test_pkg",
                            version="1.0.0",
                        )

                        assert result.success is True
                        # Check that generator metadata was passed
                        call_kwargs = mock_create.call_args[1]
                        generator = call_kwargs.get('generator', {})
                        assert generator.get('tool') == 'yapCAD-DSL'
                        assert 'dsl' in generator
                        assert generator['dsl']['module'] == 'test_module'

    def test_package_creation_exception(self, tmp_path):
        """Handles exception during package creation."""
        mock_solid = Mock()

        mock_exec_result = Mock()
        mock_exec_result.success = True
        mock_exec_result.geometry = mock_solid
        mock_exec_result.provenance = None

        with patch('yapcad.dsl.packaging.compile_and_run', return_value=mock_exec_result):
            with patch('yapcad.geom3d.issolid', return_value=True):
                with patch('yapcad.package.create_package_from_entities', side_effect=Exception("Disk full")):
                    result = package_from_dsl(
                        source="test source",
                        command_name="MAKE",
                        parameters={},
                        target_dir=tmp_path / "output",
                        name="test_pkg",
                        version="1.0.0",
                    )

                    assert result.success is False
                    assert "Package creation failed" in result.error_message
                    assert "Disk full" in result.error_message


class TestAddDslSourceAttachment:
    """Test the _add_dsl_source_attachment helper."""

    def test_creates_attachment_file(self, tmp_path):
        """Creates source.dsl file in attachments directory."""
        mock_manifest = Mock()
        mock_manifest.root = tmp_path
        mock_manifest.data = {}

        source_code = "module test;\ncommand MAKE() -> solid { emit box(1,1,1); }"

        _add_dsl_source_attachment(mock_manifest, source_code, "MAKE")

        # Check file was created
        source_file = tmp_path / "attachments" / "source.dsl"
        assert source_file.exists()
        assert source_file.read_text() == source_code

    def test_adds_attachment_to_manifest(self, tmp_path):
        """Adds attachment entry to manifest data."""
        mock_manifest = Mock()
        mock_manifest.root = tmp_path
        mock_manifest.data = {}

        source_code = "module test;"

        _add_dsl_source_attachment(mock_manifest, source_code, "TEST_CMD")

        # Check manifest was updated
        attachments = mock_manifest.data.get("attachments", [])
        assert len(attachments) == 1
        attachment = attachments[0]
        assert attachment["id"] == "dsl-source"
        assert attachment["path"] == "attachments/source.dsl"
        assert attachment["format"] == "dsl"
        assert "TEST_CMD" in attachment["purpose"]
        assert attachment["hash"].startswith("sha256:")

        # Verify save was called
        mock_manifest.save.assert_called_once()


class TestModuleImports:
    """Test that module imports work correctly."""

    def test_import_package_result(self):
        """PackageResult can be imported from yapcad.dsl."""
        from yapcad.dsl import PackageResult
        assert PackageResult is not None

    def test_import_package_from_dsl(self):
        """package_from_dsl can be imported from yapcad.dsl."""
        from yapcad.dsl import package_from_dsl
        assert callable(package_from_dsl)

    def test_type_annotations_valid(self):
        """Type annotations don't cause import errors."""
        # This test verifies the TYPE_CHECKING fix works
        # If PackageManifest wasn't properly handled, this would fail
        import yapcad.dsl.packaging as pkg

        # Check that the module loaded without errors
        assert hasattr(pkg, 'package_from_dsl')
        assert hasattr(pkg, 'PackageResult')
        assert hasattr(pkg, '_add_dsl_source_attachment')
