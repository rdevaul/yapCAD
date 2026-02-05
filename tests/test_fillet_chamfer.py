"""
Test script for fillet and chamfer operations.

This script tests both the Python API and DSL builtins for fillet/chamfer.
"""

import math
import pytest

from yapcad.geom import point
from yapcad.geom3d import solid
from yapcad.geom3d_util import prism, conic
from yapcad.brep import (
    BrepSolid,
    attach_brep_to_solid,
    brep_from_solid,
    has_brep_data,
    occ_available,
    fillet_all_edges,
    chamfer_all_edges,
)

# Skip all tests if OCC not available
pytestmark = pytest.mark.skipif(
    not occ_available(), reason="pythonocc-core not available"
)


class TestFilletPythonAPI:
    """Test fillet operations via Python API."""

    def test_fillet_box_basic(self):
        """Basic fillet test on a box."""
        from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox

        shape = BRepPrimAPI_MakeBox(20.0, 10.0, 5.0).Shape()
        brep = BrepSolid(shape)

        filleted = fillet_all_edges(brep, 1.0)

        assert filleted is not None
        assert isinstance(filleted, BrepSolid)

    def test_fillet_tessellation(self):
        """Verify filleted shape tessellates correctly."""
        from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox

        shape = BRepPrimAPI_MakeBox(10.0, 10.0, 10.0).Shape()
        brep = BrepSolid(shape)
        filleted = fillet_all_edges(brep, 0.5)

        surface = filleted.tessellate()

        assert surface[0] == 'surface'
        assert len(surface[1]) > 0  # vertices
        assert len(surface[2]) > 0  # normals
        assert len(surface[3]) > 0  # triangles

    def test_fillet_on_yapCAD_prism(self):
        """Test fillet on a prism created via yapCAD."""
        sld = prism(30.0, 20.0, 10.0)
        assert has_brep_data(sld)

        brep = brep_from_solid(sld)
        filleted = fillet_all_edges(brep, 2.0)

        surface = filleted.tessellate()
        assert surface[0] == 'surface'

    def test_fillet_small_radius(self):
        """Test fillet with very small radius."""
        from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox

        shape = BRepPrimAPI_MakeBox(50.0, 50.0, 50.0).Shape()
        brep = BrepSolid(shape)
        filleted = fillet_all_edges(brep, 0.1)

        assert filleted is not None

    def test_fillet_large_radius(self):
        """Test fillet with radius approaching edge limit."""
        from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox

        # 10x10x10 box, max fillet radius is ~5 (half the smallest dimension)
        shape = BRepPrimAPI_MakeBox(10.0, 10.0, 10.0).Shape()
        brep = BrepSolid(shape)

        # Use radius of 3 (safe margin)
        filleted = fillet_all_edges(brep, 3.0)
        assert filleted is not None


class TestChamferPythonAPI:
    """Test chamfer operations via Python API."""

    def test_chamfer_box_basic(self):
        """Basic chamfer test on a box."""
        from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox

        shape = BRepPrimAPI_MakeBox(20.0, 10.0, 5.0).Shape()
        brep = BrepSolid(shape)

        chamfered = chamfer_all_edges(brep, 0.5)

        assert chamfered is not None
        assert isinstance(chamfered, BrepSolid)

    def test_chamfer_tessellation(self):
        """Verify chamfered shape tessellates correctly."""
        from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox

        shape = BRepPrimAPI_MakeBox(10.0, 10.0, 10.0).Shape()
        brep = BrepSolid(shape)
        chamfered = chamfer_all_edges(brep, 0.5)

        surface = chamfered.tessellate()

        assert surface[0] == 'surface'
        assert len(surface[1]) > 0  # vertices
        assert len(surface[2]) > 0  # normals
        assert len(surface[3]) > 0  # triangles

    def test_chamfer_on_yapCAD_prism(self):
        """Test chamfer on a prism created via yapCAD."""
        sld = prism(30.0, 20.0, 10.0)
        assert has_brep_data(sld)

        brep = brep_from_solid(sld)
        chamfered = chamfer_all_edges(brep, 1.0)

        surface = chamfered.tessellate()
        assert surface[0] == 'surface'


class TestFilletChamferIntegration:
    """Integration tests for fillet/chamfer with yapCAD workflow."""

    def test_fillet_then_attach_brep(self):
        """Test workflow: create BREP, fillet, attach to yapCAD solid."""
        from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox

        # Create and fillet
        shape = BRepPrimAPI_MakeBox(15.0, 15.0, 15.0).Shape()
        brep = BrepSolid(shape)
        filleted = fillet_all_edges(brep, 1.5)

        # Create yapCAD solid
        surface = filleted.tessellate()
        sld = solid([surface], [], ['procedure', 'test_fillet'])
        attach_brep_to_solid(sld, filleted)

        # Verify
        assert has_brep_data(sld)
        restored = brep_from_solid(sld)
        assert restored is not None

    def test_chamfer_then_attach_brep(self):
        """Test workflow: create BREP, chamfer, attach to yapCAD solid."""
        from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox

        # Create and chamfer
        shape = BRepPrimAPI_MakeBox(15.0, 15.0, 15.0).Shape()
        brep = BrepSolid(shape)
        chamfered = chamfer_all_edges(brep, 0.8)

        # Create yapCAD solid
        surface = chamfered.tessellate()
        sld = solid([surface], [], ['procedure', 'test_chamfer'])
        attach_brep_to_solid(sld, chamfered)

        # Verify
        assert has_brep_data(sld)
        restored = brep_from_solid(sld)
        assert restored is not None


class TestDSLBuiltins:
    """Test DSL builtin functions for fillet and chamfer."""

    def test_dsl_fillet_builtin_registered(self):
        """Verify fillet is registered as a DSL builtin."""
        from yapcad.dsl.runtime.builtins import get_builtin_registry

        registry = get_builtin_registry()
        fillet_func = registry.get_function("fillet")

        assert fillet_func is not None
        assert fillet_func.name == "fillet"

    def test_dsl_chamfer_builtin_registered(self):
        """Verify chamfer is registered as a DSL builtin."""
        from yapcad.dsl.runtime.builtins import get_builtin_registry

        registry = get_builtin_registry()
        chamfer_func = registry.get_function("chamfer")

        assert chamfer_func is not None
        assert chamfer_func.name == "chamfer"

    def test_dsl_fillet_type_signature(self):
        """Verify fillet has correct type signature in symbols."""
        from yapcad.dsl.symbols import SymbolTable

        table = SymbolTable()
        sig = table.lookup_builtin("fillet")

        assert sig is not None
        assert sig.name == "fillet"
        assert len(sig.params) == 2
        assert sig.params[0][0] == "s"  # solid parameter
        assert sig.params[1][0] == "radius"  # radius parameter
        assert sig.return_type.name == "solid"

    def test_dsl_chamfer_type_signature(self):
        """Verify chamfer has correct type signature in symbols."""
        from yapcad.dsl.symbols import SymbolTable

        table = SymbolTable()
        sig = table.lookup_builtin("chamfer")

        assert sig is not None
        assert sig.name == "chamfer"
        assert len(sig.params) == 2
        assert sig.params[0][0] == "s"  # solid parameter
        assert sig.params[1][0] == "distance"  # distance parameter
        assert sig.return_type.name == "solid"


if __name__ == "__main__":
    # Run basic tests when executed directly
    print("Testing fillet and chamfer operations...")

    if not occ_available():
        print("ERROR: pythonocc-core not available. Cannot run tests.")
        exit(1)

    print("1. Testing fillet on box...")
    from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox

    shape = BRepPrimAPI_MakeBox(20.0, 10.0, 5.0).Shape()
    brep = BrepSolid(shape)
    filleted = fillet_all_edges(brep, 1.0)
    surface = filleted.tessellate()
    print(f"   Filleted box: {len(surface[1])} vertices, {len(surface[3])} triangles")

    print("2. Testing chamfer on box...")
    shape = BRepPrimAPI_MakeBox(20.0, 10.0, 5.0).Shape()
    brep = BrepSolid(shape)
    chamfered = chamfer_all_edges(brep, 0.5)
    surface = chamfered.tessellate()
    print(f"   Chamfered box: {len(surface[1])} vertices, {len(surface[3])} triangles")

    print("3. Testing fillet on yapCAD prism...")
    sld = prism(30.0, 20.0, 10.0)
    brep = brep_from_solid(sld)
    filleted = fillet_all_edges(brep, 2.0)
    surface = filleted.tessellate()
    print(f"   Filleted prism: {len(surface[1])} vertices, {len(surface[3])} triangles")

    print("4. Testing DSL builtins registered...")
    from yapcad.dsl.runtime.builtins import get_builtin_registry
    registry = get_builtin_registry()
    assert registry.get_function("fillet") is not None
    assert registry.get_function("chamfer") is not None
    print("   fillet and chamfer builtins registered correctly")

    print("\nAll tests passed!")
