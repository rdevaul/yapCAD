"""Tests for Phase 2 DSL builtins: 5 new geometry primitives from geom3d_util."""

import pytest
from yapcad.dsl.runtime.builtins import get_builtin_registry, call_builtin
from yapcad.dsl.runtime.values import Value, float_val
from yapcad.dsl.types import FLOAT, SOLID, REGION2D

try:
    from OCP.GProp import GProp_GProps
    from OCP.BRepGProp import BRepGProp
    HAS_OCC = True
except ImportError:
    HAS_OCC = False

needs_occ = pytest.mark.skipif(not HAS_OCC, reason="OCC not available")


@pytest.fixture
def registry():
    return get_builtin_registry()


class TestDodecahedron:
    def test_exists(self, registry):
        assert registry.get_function("dodecahedron") is not None

    def test_returns_solid(self):
        result = call_builtin("dodecahedron", [float_val(20.0)])
        assert result.type == SOLID
        assert result.data is not None


class TestTube:
    def test_exists(self, registry):
        assert registry.get_function("tube") is not None

    def test_returns_solid(self):
        result = call_builtin("tube", [float_val(20.0), float_val(2.0), float_val(50.0)])
        assert result.type == SOLID
        assert result.data is not None

    @needs_occ
    def test_tube_less_volume_than_cylinder(self):
        """A tube should have less volume than a solid cylinder of same OD and length."""
        from yapcad.geom3d_util import tube, conic
        tube_solid = tube(20.0, 2.0, 50.0)
        cyl_solid = conic(10.0, 10.0, 50.0)
        props_tube = GProp_GProps()
        BRepGProp.VolumeProperties_s(tube_solid, props_tube)
        props_cyl = GProp_GProps()
        BRepGProp.VolumeProperties_s(cyl_solid, props_cyl)
        assert props_tube.Mass() < props_cyl.Mass()


class TestConicTube:
    def test_exists(self, registry):
        assert registry.get_function("conic_tube") is not None

    def test_returns_solid(self):
        result = call_builtin("conic_tube", [
            float_val(30.0), float_val(20.0), float_val(2.0), float_val(50.0)
        ])
        assert result.type == SOLID
        assert result.data is not None


class TestSphericalShell:
    def test_exists(self, registry):
        assert registry.get_function("spherical_shell") is not None

    def test_returns_solid(self):
        result = call_builtin("spherical_shell", [float_val(30.0), float_val(2.0)])
        assert result.type == SOLID
        assert result.data is not None

    @needs_occ
    def test_shell_less_volume_than_sphere(self):
        from yapcad.geom3d_util import spherical_shell, sphere
        shell = spherical_shell(30.0, 2.0)
        solid_sphere = sphere(15.0)
        props_shell = GProp_GProps()
        BRepGProp.VolumeProperties_s(shell, props_shell)
        props_sphere = GProp_GProps()
        BRepGProp.VolumeProperties_s(solid_sphere, props_sphere)
        assert props_shell.Mass() < props_sphere.Mass()


class TestHelicalExtrude:
    def test_exists(self, registry):
        assert registry.get_function("helical_extrude") is not None

    @needs_occ
    def test_returns_solid(self):
        """Create a simple rectangle region and helically extrude it."""
        from yapcad.geom import point, line
        p1 = point(5, -1, 0)
        p2 = point(10, -1, 0)
        p3 = point(10, 1, 0)
        p4 = point(5, 1, 0)
        profile = [line(p1, p2), line(p2, p3), line(p3, p4), line(p4, p1)]
        result = call_builtin("helical_extrude", [
            Value(profile, REGION2D),
            float_val(20.0),
            float_val(360.0),
        ])
        assert result.type == SOLID
        assert result.data is not None
