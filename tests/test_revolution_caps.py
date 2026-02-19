"""Tests for makeRevolutionSolid disc caps (watertight geometry).

Verifies that revolution solids produce closed, manifold meshes when
the contour has non-zero radius at the endpoints (flat disc caps needed).
"""
import pytest
import math
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from yapcad.geom3d_util import makeRevolutionSolid, makeRevolutionSurface


def _mesh_from_surface(surf):
    """Extract vertices and faces from a yapCAD surface for analysis."""
    verts = surf[1]  # list of [x, y, z, 1] points
    faces = surf[2]  # list of [i, j, k] index triples (normals)
    tris = surf[3]   # list of [i, j, k] face index triples
    return verts, tris


def _check_edge_manifold(faces):
    """Check that every edge in the mesh has exactly 2 adjacent faces.
    Returns (is_manifold, boundary_edges, non_manifold_edges).
    """
    edge_count = {}
    for face in faces:
        n = len(face)
        for i in range(n):
            a, b = face[i], face[(i + 1) % n]
            edge = (min(a, b), max(a, b))
            edge_count[edge] = edge_count.get(edge, 0) + 1

    boundary = [e for e, c in edge_count.items() if c == 1]
    non_manifold = [e for e, c in edge_count.items() if c > 2]
    is_manifold = len(boundary) == 0 and len(non_manifold) == 0
    return is_manifold, boundary, non_manifold


def _check_watertight_trimesh(solid):
    """Use trimesh to verify watertightness if available."""
    try:
        import trimesh
        import numpy as np
    except ImportError:
        pytest.skip("trimesh not available")

    surf = solid[1][0]  # first surface
    verts = surf[1]
    faces = surf[3]

    vertices = np.array([[v[0], v[1], v[2]] for v in verts])
    triangles = np.array(faces)

    mesh = trimesh.Trimesh(vertices=vertices, faces=triangles)
    return mesh.is_watertight, mesh.is_volume


# ── Test Cases ────────────────────────────────────────────────


class TestRevolutionCaps:
    """Test that makeRevolutionSolid produces closed geometry."""

    def test_cylinder_revolution_is_capped(self):
        """A cylinder via revolution (constant r) should have top+bottom caps."""
        radius = 10.0
        height = 20.0
        contour = lambda z: radius
        solid = makeRevolutionSolid(contour, 0, height, 10, arcSamples=24)

        surf = solid[1][0]
        verts, faces = _mesh_from_surface(surf)

        # Should have disc caps, so mesh should be manifold
        is_manifold, boundary, non_manifold = _check_edge_manifold(faces)
        assert is_manifold, (
            f"Cylinder revolution not manifold: "
            f"{len(boundary)} boundary edges, {len(non_manifold)} non-manifold edges"
        )

    def test_truncated_cone_is_capped(self):
        """A truncated cone (r varies linearly, r>0 at both ends) should be capped."""
        def contour(z):
            return 10.0 + (z / 20.0) * (-5.0)  # r goes from 10 to 5
        solid = makeRevolutionSolid(contour, 0, 20, 10, arcSamples=24)

        surf = solid[1][0]
        _, faces = _mesh_from_surface(surf)

        is_manifold, boundary, non_manifold = _check_edge_manifold(faces)
        assert is_manifold, (
            f"Truncated cone not manifold: "
            f"{len(boundary)} boundary edges, {len(non_manifold)} non-manifold edges"
        )

    def test_dome_no_bottom_cap_needed(self):
        """A hemisphere (r=0 at top) should only cap the bottom."""
        radius = 10.0
        def contour(z):
            r_sq = radius**2 - z**2
            return math.sqrt(max(0, r_sq))
        solid = makeRevolutionSolid(contour, 0, radius, 20, arcSamples=24)

        surf = solid[1][0]
        _, faces = _mesh_from_surface(surf)

        is_manifold, boundary, non_manifold = _check_edge_manifold(faces)
        assert is_manifold, (
            f"Dome not manifold: "
            f"{len(boundary)} boundary edges, {len(non_manifold)} non-manifold edges"
        )

    def test_sphere_already_closed(self):
        """Full sphere contour (r=0 at both ends) needs no disc caps."""
        radius = 10.0
        def contour(z):
            r_sq = radius**2 - z**2
            return math.sqrt(max(0, r_sq))
        solid = makeRevolutionSolid(contour, -radius, radius, 40, arcSamples=24)

        surf = solid[1][0]
        _, faces = _mesh_from_surface(surf)

        is_manifold, boundary, non_manifold = _check_edge_manifold(faces)
        assert is_manifold, (
            f"Sphere not manifold: "
            f"{len(boundary)} boundary edges, {len(non_manifold)} non-manifold edges"
        )

    def test_cylinder_watertight_trimesh(self):
        """Verify cylinder revolution is watertight via trimesh."""
        radius = 10.0
        contour = lambda z: radius
        solid = makeRevolutionSolid(contour, 0, 20, 10, arcSamples=24)

        is_watertight, is_volume = _check_watertight_trimesh(solid)
        assert is_watertight, "Cylinder revolution mesh not watertight"
        assert is_volume, "Cylinder revolution mesh not a valid volume"

    def test_truncated_cone_watertight_trimesh(self):
        """Verify truncated cone is watertight via trimesh."""
        def contour(z):
            return 10.0 + (z / 20.0) * (-5.0)
        solid = makeRevolutionSolid(contour, 0, 20, 10, arcSamples=24)

        is_watertight, is_volume = _check_watertight_trimesh(solid)
        assert is_watertight, "Truncated cone mesh not watertight"
        assert is_volume, "Truncated cone mesh not a valid volume"

    def test_dome_watertight_trimesh(self):
        """Verify hemisphere with bottom cap is watertight via trimesh."""
        radius = 10.0
        def contour(z):
            r_sq = radius**2 - z**2
            return math.sqrt(max(0, r_sq))
        solid = makeRevolutionSolid(contour, 0, radius, 20, arcSamples=24)

        is_watertight, is_volume = _check_watertight_trimesh(solid)
        assert is_watertight, "Dome mesh not watertight"
        assert is_volume, "Dome mesh not a valid volume"

    def test_bottle_shape_capped(self):
        """A bottle shape (wider middle, narrow top, flat bottom) should be capped."""
        def contour(z):
            if z < 5:
                return 8.0  # flat bottom cylinder
            elif z < 15:
                return 8.0 + 4.0 * math.sin((z - 5) / 10.0 * math.pi)  # bulge
            else:
                return 3.0  # narrow neck
        solid = makeRevolutionSolid(contour, 0, 20, 40, arcSamples=24)

        surf = solid[1][0]
        _, faces = _mesh_from_surface(surf)

        is_manifold, boundary, non_manifold = _check_edge_manifold(faces)
        assert is_manifold, (
            f"Bottle shape not manifold: "
            f"{len(boundary)} boundary edges, {len(non_manifold)} non-manifold edges"
        )

    def test_tank_contour_capped(self):
        """Simulate the tank v3 contour: dome top, cylindrical body, open bottom."""
        radius = 25.0
        dome_height = 15.0
        body_height = 50.0
        total = dome_height + body_height

        def contour(z):
            if z <= body_height:
                return radius  # cylindrical body
            else:
                # ellipsoidal dome
                t = (z - body_height) / dome_height
                r_sq = 1.0 - t**2
                return radius * math.sqrt(max(0, r_sq))

        solid = makeRevolutionSolid(contour, 0, total, 40, arcSamples=32)

        surf = solid[1][0]
        _, faces = _mesh_from_surface(surf)

        is_manifold, boundary, non_manifold = _check_edge_manifold(faces)
        assert is_manifold, (
            f"Tank contour not manifold: "
            f"{len(boundary)} boundary edges, {len(non_manifold)} non-manifold edges"
        )


class TestRevolutionCapsRegression:
    """Regression tests: ensure existing behavior isn't broken."""

    def test_cone_to_point(self):
        """A cone (r goes to 0) should still work with pole cap."""
        def contour(z):
            return 10.0 * (1.0 - z / 20.0)
        solid = makeRevolutionSolid(contour, 0, 20, 10, arcSamples=24)
        surf = solid[1][0]
        verts = surf[1]
        faces = surf[3]
        assert len(verts) > 0
        assert len(faces) > 0

    def test_surface_only_mode(self):
        """makeRevolutionSurface without return_brep should still work."""
        contour = lambda z: 5.0
        surf = makeRevolutionSurface(contour, 0, 10, 5, arcSamples=12)
        assert surf[0] == 'surface'
        assert len(surf[1]) > 0  # vertices
        assert len(surf[3]) > 0  # faces
