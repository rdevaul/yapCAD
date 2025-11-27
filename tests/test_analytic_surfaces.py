"""Tests for the analytic surfaces module."""

import pytest
from math import pi, sqrt

from yapcad.analytic_surfaces import (
    plane_surface, is_plane_surface, evaluate_plane_surface, plane_surface_normal,
    sphere_surface, is_sphere_surface, evaluate_sphere_surface, sphere_surface_normal,
    cylinder_surface, is_cylinder_surface, evaluate_cylinder_surface, cylinder_surface_normal,
    cone_surface, is_cone_surface, evaluate_cone_surface, cone_surface_normal,
    torus_surface, is_torus_surface, evaluate_torus_surface, torus_surface_normal,
    is_analytic_surface, evaluate_surface, surface_normal, tessellate_surface
)


class TestPlaneSurface:
    """Test plane surface operations."""

    def test_create_plane(self):
        """Test plane surface creation."""
        p = plane_surface([0, 0, 0], [0, 0, 1])
        assert is_plane_surface(p)
        assert p[0] == 'plane_surface'

    def test_plane_evaluation(self):
        """Test point evaluation on plane."""
        p = plane_surface([0, 0, 0], [0, 0, 1])
        pt = evaluate_plane_surface(p, 0.5, 0.5)
        # Point should be at some offset in the XY plane at z=0
        assert abs(pt[2]) < 1e-10

    def test_plane_normal(self):
        """Test plane normal is constant."""
        p = plane_surface([0, 0, 0], [0, 0, 1])
        n1 = plane_surface_normal(p, 0, 0)
        n2 = plane_surface_normal(p, 0.5, 0.5)
        assert abs(n1[2] - 1.0) < 1e-10
        assert abs(n2[2] - 1.0) < 1e-10


class TestSphereSurface:
    """Test sphere surface operations."""

    def test_create_sphere(self):
        """Test sphere surface creation."""
        s = sphere_surface([0, 0, 0], 1.0)
        assert is_sphere_surface(s)
        assert s[2]['radius'] == 1.0

    def test_sphere_evaluation(self):
        """Test point evaluation on sphere."""
        s = sphere_surface([0, 0, 0], 1.0)

        # At (u=0, v=0), should be at (1, 0, 0)
        pt = evaluate_sphere_surface(s, 0, 0)
        assert abs(pt[0] - 1.0) < 1e-10
        assert abs(pt[1]) < 1e-10
        assert abs(pt[2]) < 1e-10

        # At (u=pi/2, v=0), should be at (0, 1, 0)
        pt2 = evaluate_sphere_surface(s, pi/2, 0)
        assert abs(pt2[0]) < 1e-10
        assert abs(pt2[1] - 1.0) < 1e-10

        # At (u=0, v=pi/2), should be at (0, 0, 1)
        pt3 = evaluate_sphere_surface(s, 0, pi/2)
        assert abs(pt3[0]) < 1e-10
        assert abs(pt3[2] - 1.0) < 1e-10

    def test_sphere_normal(self):
        """Test sphere normal is outward."""
        s = sphere_surface([0, 0, 0], 1.0)
        n = sphere_surface_normal(s, 0, 0)
        # Normal at (1,0,0) should point in +x direction
        assert abs(n[0] - 1.0) < 1e-10
        assert abs(n[1]) < 1e-10


class TestCylinderSurface:
    """Test cylinder surface operations."""

    def test_create_cylinder(self):
        """Test cylinder surface creation."""
        c = cylinder_surface([0, 0, 0], [0, 0, 1], 1.0)
        assert is_cylinder_surface(c)
        assert c[2]['radius'] == 1.0

    def test_cylinder_evaluation(self):
        """Test point evaluation on cylinder."""
        c = cylinder_surface([0, 0, 0], [0, 0, 1], 1.0, v_range=(0, 2))
        pt = evaluate_cylinder_surface(c, 0, 1)
        # Should be at radius 1, height 1
        r = sqrt(pt[0]**2 + pt[1]**2)
        assert abs(r - 1.0) < 1e-10
        assert abs(pt[2] - 1.0) < 1e-10


class TestConeSurface:
    """Test cone surface operations."""

    def test_create_cone(self):
        """Test cone surface creation."""
        c = cone_surface([0, 0, 0], [0, 0, 1], pi/4)
        assert is_cone_surface(c)
        assert abs(c[2]['half_angle'] - pi/4) < 1e-10

    def test_cone_evaluation(self):
        """Test point evaluation on cone."""
        c = cone_surface([0, 0, 0], [0, 0, 1], pi/4, v_range=(0, 2))
        pt = evaluate_cone_surface(c, 0, 1)
        # At distance 1 from apex with 45 degree half-angle
        # radius should equal height
        r = sqrt(pt[0]**2 + pt[1]**2)
        assert abs(r - pt[2]) < 0.1


class TestTorusSurface:
    """Test torus surface operations."""

    def test_create_torus(self):
        """Test torus surface creation."""
        t = torus_surface([0, 0, 0], [0, 0, 1], 2.0, 0.5)
        assert is_torus_surface(t)
        assert t[2]['major_radius'] == 2.0
        assert t[2]['minor_radius'] == 0.5

    def test_torus_evaluation(self):
        """Test point evaluation on torus."""
        t = torus_surface([0, 0, 0], [0, 0, 1], 2.0, 0.5)
        pt = evaluate_torus_surface(t, 0, 0)
        # At (u=0, v=0), should be at (0, 2.5, 0) (outer edge)
        r_xy = sqrt(pt[0]**2 + pt[1]**2)
        assert abs(r_xy - 2.5) < 1e-10
        assert abs(pt[2]) < 1e-10


class TestGenericOperations:
    """Test generic surface operations."""

    def test_is_analytic_surface(self):
        """Test generic type checking."""
        p = plane_surface([0, 0, 0], [0, 0, 1])
        s = sphere_surface([0, 0, 0], 1.0)
        assert is_analytic_surface(p)
        assert is_analytic_surface(s)
        assert not is_analytic_surface([1, 2, 3])

    def test_generic_evaluate_surface(self):
        """Test generic evaluate_surface function."""
        s = sphere_surface([0, 0, 0], 1.0)
        pt = evaluate_surface(s, 0, 0)
        assert abs(pt[0] - 1.0) < 1e-10

    def test_generic_surface_normal(self):
        """Test generic surface_normal function."""
        s = sphere_surface([0, 0, 0], 1.0)
        n = surface_normal(s, 0, 0)
        assert abs(n[0] - 1.0) < 1e-10


class TestTessellation:
    """Test surface tessellation."""

    def test_tessellate_sphere(self):
        """Test sphere tessellation."""
        s = sphere_surface([0, 0, 0], 1.0)
        mesh = tessellate_surface(s, u_divisions=8, v_divisions=4)

        assert mesh[0] == 'surface'
        assert len(mesh[1]) > 0  # Has vertices
        assert len(mesh[2]) > 0  # Has normals
        assert len(mesh[3]) > 0  # Has faces

    def test_tessellate_plane(self):
        """Test plane tessellation."""
        p = plane_surface([0, 0, 0], [0, 0, 1])
        mesh = tessellate_surface(p, u_divisions=4, v_divisions=4)

        assert mesh[0] == 'surface'
        assert len(mesh[1]) == 25  # (4+1)*(4+1) vertices

    def test_tessellation_is_valid_surface(self):
        """Test that tessellation produces valid yapcad surface."""
        from yapcad.geom3d import issurface
        s = sphere_surface([0, 0, 0], 1.0)
        mesh = tessellate_surface(s, u_divisions=8, v_divisions=4)
        assert issurface(mesh)
