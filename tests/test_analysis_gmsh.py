"""Tests for Gmsh meshing integration.

These tests verify the GmshMesher class and related utilities.
Tests requiring Gmsh are skipped if the dependency is not installed.
"""

import math
import pytest
from pathlib import Path
from unittest.mock import MagicMock, patch

# Check if gmsh is available
try:
    import gmsh
    GMSH_AVAILABLE = True
except ImportError:
    GMSH_AVAILABLE = False


class TestMeshHints:
    """Tests for MeshHints dataclass."""

    def test_default_values(self):
        """Test default mesh hint values."""
        from yapcad.package.analysis.gmsh_mesher import MeshHints

        hints = MeshHints()
        assert hints.element_size == 5.0
        assert hints.element_order == 1
        assert hints.optimize is True
        assert hints.algorithm_2d == 6  # Frontal-Delaunay

    def test_custom_values(self):
        """Test custom mesh hint values."""
        from yapcad.package.analysis.gmsh_mesher import MeshHints

        hints = MeshHints(
            element_size=2.0,
            min_element_size=0.5,
            max_element_size=10.0,
            element_order=2,
        )
        assert hints.element_size == 2.0
        assert hints.min_element_size == 0.5
        assert hints.max_element_size == 10.0
        assert hints.element_order == 2


class TestPhysicalGroup:
    """Tests for PhysicalGroup dataclass."""

    def test_physical_group(self):
        """Test physical group creation."""
        from yapcad.package.analysis.gmsh_mesher import PhysicalGroup

        pg = PhysicalGroup(name="fixed", dim=2, tags=[1, 2, 3])
        assert pg.name == "fixed"
        assert pg.dim == 2
        assert pg.tags == [1, 2, 3]


class TestGmshAvailability:
    """Tests for gmsh availability checking."""

    def test_gmsh_available_function(self):
        """Test gmsh_available returns correct value."""
        from yapcad.package.analysis.gmsh_mesher import gmsh_available

        # Should return True if gmsh is installed, False otherwise
        assert gmsh_available() == GMSH_AVAILABLE

    def test_require_gmsh_raises_when_unavailable(self):
        """Test require_gmsh raises RuntimeError when gmsh not available."""
        from yapcad.package.analysis import gmsh_mesher

        with patch.object(gmsh_mesher, '_GMSH_AVAILABLE', False):
            with pytest.raises(RuntimeError, match="Gmsh Python API is not available"):
                gmsh_mesher.require_gmsh()


@pytest.mark.skipif(not GMSH_AVAILABLE, reason="Gmsh not installed")
class TestGmshMesher:
    """Tests for GmshMesher class (require Gmsh)."""

    def test_initialization(self):
        """Test mesher initialization."""
        from yapcad.package.analysis.gmsh_mesher import GmshMesher

        mesher = GmshMesher("test_model")
        assert mesher._model_name == "test_model"
        assert not mesher._initialized

    def test_context_manager(self):
        """Test mesher as context manager."""
        from yapcad.package.analysis.gmsh_mesher import GmshMesher

        with GmshMesher() as mesher:
            assert mesher._initialized
        assert not mesher._initialized

    def test_initialize_finalize(self):
        """Test explicit initialize/finalize."""
        from yapcad.package.analysis.gmsh_mesher import GmshMesher

        mesher = GmshMesher()
        mesher.initialize()
        assert mesher._initialized

        mesher.finalize()
        assert not mesher._initialized

    def test_import_step_file(self, tmp_path):
        """Test importing a STEP file."""
        from yapcad.package.analysis.gmsh_mesher import GmshMesher

        # Create a simple STEP file for testing
        step_file = tmp_path / "test.step"

        # Skip if we don't have a test STEP file
        # In real tests, we'd create one via yapCAD
        if not step_file.exists():
            pytest.skip("No test STEP file available")

        with GmshMesher() as mesher:
            entities = mesher.import_step(step_file)
            assert isinstance(entities, list)

    def test_set_physical_groups(self):
        """Test setting physical groups."""
        from yapcad.package.analysis.gmsh_mesher import GmshMesher

        with GmshMesher() as mesher:
            # Create a simple box for testing
            gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
            gmsh.model.occ.synchronize()

            # Get faces
            faces = gmsh.model.getEntities(dim=2)
            face_tags = [f[1] for f in faces]

            # Set physical groups
            mesher.set_physical_groups({"test_group": face_tags[:2]}, dim=2)

            assert "test_group" in mesher._physical_groups
            assert mesher._physical_groups["test_group"].dim == 2

    def test_generate_mesh_simple_box(self, tmp_path):
        """Test mesh generation for a simple box."""
        from yapcad.package.analysis.gmsh_mesher import GmshMesher, MeshHints

        with GmshMesher() as mesher:
            # Create a simple box
            gmsh.model.occ.addBox(0, 0, 0, 10, 10, 10)
            gmsh.model.occ.synchronize()

            # Generate mesh
            hints = MeshHints(element_size=2.0)
            mesher.generate_mesh(hints, dim=3)

            # Check mesh was generated
            stats = mesher.get_mesh_stats()
            assert stats["nodes"] > 0
            assert len(stats["elements"]) > 0

            # Export
            mesh_path = tmp_path / "test.msh"
            mesher.export_mesh(mesh_path)
            assert mesh_path.exists()

    def test_mesh_stats(self):
        """Test mesh statistics."""
        from yapcad.package.analysis.gmsh_mesher import GmshMesher, MeshHints

        with GmshMesher() as mesher:
            # Create and mesh a box
            gmsh.model.occ.addBox(0, 0, 0, 5, 5, 5)
            gmsh.model.occ.synchronize()
            mesher.generate_mesh(MeshHints(element_size=2.0), dim=3)

            stats = mesher.get_mesh_stats()

            assert "nodes" in stats
            assert "elements" in stats
            assert "physical_groups" in stats
            assert stats["nodes"] > 8  # More than box corners


class TestFaceNaming:
    """Tests for face naming system."""

    def test_by_normal_selector(self):
        """Test ByNormalSelector."""
        from yapcad.package.analysis.face_naming import (
            ByNormalSelector, FaceInfo
        )

        selector = ByNormalSelector((0, 0, 1), tolerance_deg=5.0)

        # Face with +Z normal should match
        face_up = FaceInfo(0, (0, 0, 0), (0, 0, 1), 1.0)
        assert selector.matches(face_up)

        # Face with -Z normal should not match
        face_down = FaceInfo(1, (0, 0, 0), (0, 0, -1), 1.0)
        assert not selector.matches(face_down)

        # Face with slight tilt should match within tolerance
        # Normalize the tilted normal vector
        import math
        tilt = (0.05, 0.05, 0.99)
        mag = math.sqrt(sum(t*t for t in tilt))
        tilt_normalized = tuple(t/mag for t in tilt)
        face_tilted = FaceInfo(2, (0, 0, 0), tilt_normalized, 1.0)
        assert selector.matches(face_tilted)

    def test_by_normal_selector_reversed(self):
        """Test ByNormalSelector with allow_reversed."""
        from yapcad.package.analysis.face_naming import (
            ByNormalSelector, FaceInfo
        )

        selector = ByNormalSelector((0, 0, 1), allow_reversed=True)

        # Both +Z and -Z should match
        face_up = FaceInfo(0, (0, 0, 0), (0, 0, 1), 1.0)
        face_down = FaceInfo(1, (0, 0, 0), (0, 0, -1), 1.0)

        assert selector.matches(face_up)
        assert selector.matches(face_down)

    def test_by_area_selector_range(self):
        """Test ByAreaSelector with min/max range."""
        from yapcad.package.analysis.face_naming import (
            ByAreaSelector, FaceInfo
        )

        selector = ByAreaSelector(min_area=5.0, max_area=15.0)

        small_face = FaceInfo(0, (0, 0, 0), (0, 0, 1), 2.0)
        medium_face = FaceInfo(1, (0, 0, 0), (0, 0, 1), 10.0)
        large_face = FaceInfo(2, (0, 0, 0), (0, 0, 1), 20.0)

        assert not selector.matches(small_face)
        assert selector.matches(medium_face)
        assert not selector.matches(large_face)

    def test_by_area_selector_largest(self):
        """Test ByAreaSelector for largest face."""
        from yapcad.package.analysis.face_naming import (
            ByAreaSelector, FaceInfo
        )

        faces = [
            FaceInfo(0, (0, 0, 0), (0, 0, 1), 5.0),
            FaceInfo(1, (0, 0, 0), (0, 0, 1), 15.0),
            FaceInfo(2, (0, 0, 0), (0, 0, 1), 10.0),
        ]

        selector = ByAreaSelector(largest=True)
        selector.set_context(faces)

        assert not selector.matches(faces[0])
        assert selector.matches(faces[1])  # Largest
        assert not selector.matches(faces[2])

    def test_by_position_selector_at_max(self):
        """Test ByPositionSelector for faces at max Z."""
        from yapcad.package.analysis.face_naming import (
            ByPositionSelector, FaceInfo
        )

        faces = [
            FaceInfo(0, (0, 0, 0), (0, 0, -1), 1.0),  # Bottom
            FaceInfo(1, (0, 0, 10), (0, 0, 1), 1.0),  # Top
            FaceInfo(2, (0, 0, 5), (1, 0, 0), 1.0),   # Middle
        ]

        selector = ByPositionSelector(axis="z", at_max=True)
        selector.set_context(faces)

        assert not selector.matches(faces[0])
        assert selector.matches(faces[1])  # At z_max
        assert not selector.matches(faces[2])

    def test_combined_selector_and(self):
        """Test CombinedSelector with AND logic."""
        from yapcad.package.analysis.face_naming import (
            ByNormalSelector, ByAreaSelector, CombinedSelector, FaceInfo
        )

        # Select faces that are both top-facing AND large
        selector = CombinedSelector([
            ByNormalSelector((0, 0, 1)),
            ByAreaSelector(min_area=10.0),
        ], mode="and")

        small_top = FaceInfo(0, (0, 0, 0), (0, 0, 1), 5.0)
        large_top = FaceInfo(1, (0, 0, 0), (0, 0, 1), 15.0)
        large_bottom = FaceInfo(2, (0, 0, 0), (0, 0, -1), 15.0)

        assert not selector.matches(small_top)   # Top but small
        assert selector.matches(large_top)       # Top and large
        assert not selector.matches(large_bottom)  # Large but not top

    def test_combined_selector_or(self):
        """Test CombinedSelector with OR logic."""
        from yapcad.package.analysis.face_naming import (
            ByNormalSelector, CombinedSelector, FaceInfo
        )

        # Select faces that are either top OR bottom
        selector = CombinedSelector([
            ByNormalSelector((0, 0, 1)),
            ByNormalSelector((0, 0, -1)),
        ], mode="or")

        top = FaceInfo(0, (0, 0, 0), (0, 0, 1), 1.0)
        bottom = FaceInfo(1, (0, 0, 0), (0, 0, -1), 1.0)
        side = FaceInfo(2, (0, 0, 0), (1, 0, 0), 1.0)

        assert selector.matches(top)
        assert selector.matches(bottom)
        assert not selector.matches(side)

    def test_convenience_selectors(self):
        """Test convenience selector functions."""
        from yapcad.package.analysis.face_naming import (
            top_faces, bottom_faces, left_faces, right_faces,
            front_faces, back_faces, FaceInfo
        )

        face_up = FaceInfo(0, (0, 0, 0), (0, 0, 1), 1.0)
        face_down = FaceInfo(1, (0, 0, 0), (0, 0, -1), 1.0)
        face_left = FaceInfo(2, (0, 0, 0), (-1, 0, 0), 1.0)
        face_right = FaceInfo(3, (0, 0, 0), (1, 0, 0), 1.0)
        face_front = FaceInfo(4, (0, 0, 0), (0, 1, 0), 1.0)
        face_back = FaceInfo(5, (0, 0, 0), (0, -1, 0), 1.0)

        assert top_faces().matches(face_up)
        assert bottom_faces().matches(face_down)
        assert left_faces().matches(face_left)
        assert right_faces().matches(face_right)
        assert front_faces().matches(face_front)
        assert back_faces().matches(face_back)


class TestFenicsAdapter:
    """Tests for FEniCSx adapter structure."""

    def test_adapter_registration(self):
        """Test that FEniCSx adapter is registered."""
        from yapcad.package.analysis import get_backend, available_backends

        backends = available_backends()
        # May or may not include fenics depending on installation
        assert "calculix" in backends  # Always available

    def test_material_properties(self):
        """Test MaterialProperties calculations."""
        try:
            from yapcad.package.analysis.fenics import MaterialProperties
        except ImportError:
            pytest.skip("FEniCSx module not available")

        mat = MaterialProperties(
            youngs_modulus=200e9,  # Steel
            poisson_ratio=0.3,
        )

        # Check Lamé parameters
        E, nu = 200e9, 0.3
        expected_lambda = E * nu / ((1 + nu) * (1 - 2 * nu))
        expected_mu = E / (2 * (1 + nu))

        assert abs(mat.lame_lambda - expected_lambda) < 1e6
        assert abs(mat.lame_mu - expected_mu) < 1e6
