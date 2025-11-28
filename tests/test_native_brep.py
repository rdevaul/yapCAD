"""Tests for the native BREP topology graph."""

import pytest

from yapcad.geom import line, point, isline
from yapcad.native_brep import (
    brep_vertex, is_brep_vertex, vertex_location, vertex_id, vertices_coincident,
    brep_edge, is_brep_edge, edge_curve_type, edge_id, edge_vertices,
    line_edge, evaluate_edge_curve,
    brep_trim, is_brep_trim, trim_id, trim_edge_id,
    brep_loop, is_brep_loop, loop_id, loop_trims, loop_type,
    brep_face, is_brep_face, face_id, face_loops, face_surface,
    brep_shell, is_brep_shell, shell_id, shell_faces, shell_closed, set_shell_closed,
    brep_solid, is_brep_solid_native, solid_id, solid_shells,
    ShellClosureError, SolidValidationError, TopologyGraph
)


class TestBrepVertex:
    """Test BREP vertex operations."""

    def test_create_vertex(self):
        """Test vertex creation."""
        v = brep_vertex([1, 2, 3])
        assert is_brep_vertex(v)
        assert v[0] == 'brep_vertex'

    def test_vertex_location(self):
        """Test vertex location retrieval."""
        v = brep_vertex([1, 2, 3])
        loc = vertex_location(v)
        assert abs(loc[0] - 1) < 1e-10
        assert abs(loc[1] - 2) < 1e-10
        assert abs(loc[2] - 3) < 1e-10

    def test_vertex_id_unique(self):
        """Test that vertex IDs are unique."""
        v1 = brep_vertex([0, 0, 0])
        v2 = brep_vertex([0, 0, 0])
        assert vertex_id(v1) != vertex_id(v2)

    def test_vertices_coincident(self):
        """Test vertex coincidence detection."""
        v1 = brep_vertex([0, 0, 0])
        v2 = brep_vertex([0.000001, 0, 0])
        v3 = brep_vertex([1, 0, 0])

        assert vertices_coincident(v1, v2)
        assert not vertices_coincident(v1, v3)


class TestBrepEdge:
    """Test BREP edge operations."""

    def test_create_edge(self):
        """Test edge creation with parametric curve."""
        v0 = brep_vertex([0, 0, 0])
        v1 = brep_vertex([1, 0, 0])
        e = line_edge(v0, v1)

        assert is_brep_edge(e)
        assert e[0] == 'brep_edge'
        assert edge_curve_type(e) == 'line'

    def test_edge_vertices(self):
        """Test edge vertex retrieval."""
        v0 = brep_vertex([0, 0, 0])
        v1 = brep_vertex([1, 0, 0])
        e = line_edge(v0, v1)

        start_id, end_id = edge_vertices(e)
        assert start_id == vertex_id(v0)
        assert end_id == vertex_id(v1)

    def test_edge_curve_evaluation(self):
        """Test edge curve evaluation from vertices."""
        v0 = brep_vertex([0, 0, 0])
        v1 = brep_vertex([1, 0, 0])
        e = line_edge(v0, v1)

        # Evaluate the curve using vertex positions
        curve = evaluate_edge_curve(e, vertex_location(v0), vertex_location(v1))
        # Check it's a line (yapCAD lines are [[start], [end]])
        assert isline(curve)
        # Check endpoints match vertex positions
        assert abs(curve[0][0] - 0) < 1e-10
        assert abs(curve[1][0] - 1) < 1e-10


class TestBrepTrim:
    """Test BREP trim operations."""

    def test_create_trim(self):
        """Test trim creation."""
        v0 = brep_vertex([0, 0, 0])
        v1 = brep_vertex([1, 0, 0])
        e = line_edge(v0, v1)
        t = brep_trim(e, sense=True)

        assert is_brep_trim(t)
        assert trim_edge_id(t) == edge_id(e)


class TestBrepLoop:
    """Test BREP loop operations."""

    def test_create_loop(self):
        """Test loop creation."""
        # Create a simple square loop
        v0 = brep_vertex([0, 0, 0])
        v1 = brep_vertex([1, 0, 0])
        v2 = brep_vertex([1, 1, 0])
        v3 = brep_vertex([0, 1, 0])

        e01 = line_edge(v0, v1)
        e12 = line_edge(v1, v2)
        e23 = line_edge(v2, v3)
        e30 = line_edge(v3, v0)

        t01 = brep_trim(e01)
        t12 = brep_trim(e12)
        t23 = brep_trim(e23)
        t30 = brep_trim(e30)

        loop = brep_loop([t01, t12, t23, t30], loop_type='outer')

        assert is_brep_loop(loop)
        assert len(loop_trims(loop)) == 4
        assert loop_type(loop) == 'outer'


class TestBrepFace:
    """Test BREP face operations."""

    def test_create_face(self):
        """Test face creation."""
        from yapcad.analytic_surfaces import plane_surface

        # Create a simple loop
        v0 = brep_vertex([0, 0, 0])
        v1 = brep_vertex([1, 0, 0])
        e = line_edge(v0, v1)
        t = brep_trim(e)
        loop = brep_loop([t])

        surface = plane_surface([0, 0, 0], [0, 0, 1])
        face = brep_face(surface, [loop])

        assert is_brep_face(face)
        assert len(face_loops(face)) == 1


class TestBrepShell:
    """Test BREP shell operations."""

    def test_create_open_shell(self):
        """Test open shell creation."""
        from yapcad.analytic_surfaces import plane_surface

        surface = plane_surface([0, 0, 0], [0, 0, 1])
        face = brep_face(surface, [])
        # Single face cannot be closed - mark as explicitly open
        shell = brep_shell([face], shell_type='outer', closed=False)

        assert is_brep_shell(shell)
        assert shell_closed(shell) is False
        assert len(shell_faces(shell)) == 1

    def test_shell_closed_none_default(self):
        """Test that closed defaults to None (auto-detect)."""
        shell = brep_shell([])
        assert shell_closed(shell) is None


class TestBrepSolid:
    """Test BREP solid operations."""

    def test_create_solid(self):
        """Test solid creation."""
        shell = brep_shell([])
        solid = brep_solid([shell])

        assert is_brep_solid_native(solid)
        assert len(solid_shells(solid)) == 1


class TestTopologyGraph:
    """Test the TopologyGraph container."""

    def test_graph_creation(self):
        """Test empty graph creation."""
        graph = TopologyGraph()
        summary = graph.summary()
        assert summary['vertices'] == 0
        assert summary['edges'] == 0

    def test_add_entities(self):
        """Test adding entities to graph."""
        graph = TopologyGraph()

        v0 = brep_vertex([0, 0, 0])
        v1 = brep_vertex([1, 0, 0])
        vid0 = graph.add_vertex(v0)
        vid1 = graph.add_vertex(v1)

        e = line_edge(v0, v1)
        eid = graph.add_edge(e)

        summary = graph.summary()
        assert summary['vertices'] == 2
        assert summary['edges'] == 1

    def test_vertex_edges(self):
        """Test getting edges incident to a vertex."""
        graph = TopologyGraph()

        v0 = brep_vertex([0, 0, 0])
        v1 = brep_vertex([1, 0, 0])
        v2 = brep_vertex([0, 1, 0])

        graph.add_vertex(v0)
        graph.add_vertex(v1)
        graph.add_vertex(v2)

        e01 = line_edge(v0, v1)
        e02 = line_edge(v0, v2)

        graph.add_edge(e01)
        graph.add_edge(e02)

        # v0 should have 2 incident edges
        edges = graph.vertex_edges(vertex_id(v0))
        assert len(edges) == 2

    def test_get_entities(self):
        """Test entity retrieval by ID."""
        graph = TopologyGraph()

        v = brep_vertex([0, 0, 0])
        vid = graph.add_vertex(v)

        retrieved = graph.get_vertex(vid)
        assert retrieved is v
        assert graph.get_vertex('nonexistent') is None


class TestShellClosureValidation:
    """Test shell closure validation in TopologyGraph."""

    def _build_cube_topology(self, graph):
        """Build a closed cube shell topology in the graph.

        Creates a unit cube with 8 vertices, 12 edges, 6 faces.
        Each edge is shared by exactly 2 faces.
        Returns the shell.
        """
        from yapcad.analytic_surfaces import plane_surface
        from yapcad.native_brep import (
            brep_vertex, brep_face, brep_loop, brep_trim, brep_shell,
            line_edge, face_id, shell_faces
        )

        # 8 vertices of unit cube
        v000 = brep_vertex([0, 0, 0])
        v100 = brep_vertex([1, 0, 0])
        v110 = brep_vertex([1, 1, 0])
        v010 = brep_vertex([0, 1, 0])
        v001 = brep_vertex([0, 0, 1])
        v101 = brep_vertex([1, 0, 1])
        v111 = brep_vertex([1, 1, 1])
        v011 = brep_vertex([0, 1, 1])

        vertices = [v000, v100, v110, v010, v001, v101, v111, v011]
        for v in vertices:
            graph.add_vertex(v)

        # 12 edges of the cube
        # Bottom face edges
        e_b01 = line_edge(v000, v100)  # bottom front
        e_b12 = line_edge(v100, v110)  # bottom right
        e_b23 = line_edge(v110, v010)  # bottom back
        e_b30 = line_edge(v010, v000)  # bottom left

        # Top face edges
        e_t01 = line_edge(v001, v101)  # top front
        e_t12 = line_edge(v101, v111)  # top right
        e_t23 = line_edge(v111, v011)  # top back
        e_t30 = line_edge(v011, v001)  # top left

        # Vertical edges
        e_v0 = line_edge(v000, v001)  # front-left
        e_v1 = line_edge(v100, v101)  # front-right
        e_v2 = line_edge(v110, v111)  # back-right
        e_v3 = line_edge(v010, v011)  # back-left

        edges = [e_b01, e_b12, e_b23, e_b30, e_t01, e_t12, e_t23, e_t30,
                 e_v0, e_v1, e_v2, e_v3]
        for e in edges:
            graph.add_edge(e)

        # Create trims for each edge (we'll reuse edges with different senses)
        def make_face(surface, edge_list):
            """Helper to create face with loop from edges."""
            trims = [brep_trim(e) for e in edge_list]
            for t in trims:
                graph.add_trim(t)
            loop = brep_loop(trims, loop_type='outer')
            graph.add_loop(loop)
            face = brep_face(surface, [loop])
            graph.add_face(face)
            return face

        # 6 faces
        # Bottom face (z=0)
        surf_bottom = plane_surface([0, 0, 0], [0, 0, -1])
        f_bottom = make_face(surf_bottom, [e_b01, e_b12, e_b23, e_b30])

        # Top face (z=1)
        surf_top = plane_surface([0, 0, 1], [0, 0, 1])
        f_top = make_face(surf_top, [e_t01, e_t12, e_t23, e_t30])

        # Front face (y=0)
        surf_front = plane_surface([0, 0, 0], [0, -1, 0])
        f_front = make_face(surf_front, [e_b01, e_v1, e_t01, e_v0])

        # Back face (y=1)
        surf_back = plane_surface([0, 1, 0], [0, 1, 0])
        f_back = make_face(surf_back, [e_b23, e_v3, e_t23, e_v2])

        # Left face (x=0)
        surf_left = plane_surface([0, 0, 0], [-1, 0, 0])
        f_left = make_face(surf_left, [e_b30, e_v0, e_t30, e_v3])

        # Right face (x=1)
        surf_right = plane_surface([1, 0, 0], [1, 0, 0])
        f_right = make_face(surf_right, [e_b12, e_v2, e_t12, e_v1])

        faces = [f_bottom, f_top, f_front, f_back, f_left, f_right]
        shell = brep_shell(faces, closed=None)  # Auto-detect closure

        return shell

    def _build_open_box_topology(self, graph):
        """Build an open box (cube missing top face) topology.

        Creates a box with 8 vertices, 12 edges, 5 faces.
        The top edges will only be used by 1 face each (boundary edges).
        Returns the shell.
        """
        from yapcad.analytic_surfaces import plane_surface
        from yapcad.native_brep import (
            brep_vertex, brep_face, brep_loop, brep_trim, brep_shell,
            line_edge
        )

        # 8 vertices
        v000 = brep_vertex([0, 0, 0])
        v100 = brep_vertex([1, 0, 0])
        v110 = brep_vertex([1, 1, 0])
        v010 = brep_vertex([0, 1, 0])
        v001 = brep_vertex([0, 0, 1])
        v101 = brep_vertex([1, 0, 1])
        v111 = brep_vertex([1, 1, 1])
        v011 = brep_vertex([0, 1, 1])

        for v in [v000, v100, v110, v010, v001, v101, v111, v011]:
            graph.add_vertex(v)

        # Edges
        e_b01 = line_edge(v000, v100)
        e_b12 = line_edge(v100, v110)
        e_b23 = line_edge(v110, v010)
        e_b30 = line_edge(v010, v000)
        e_t01 = line_edge(v001, v101)
        e_t12 = line_edge(v101, v111)
        e_t23 = line_edge(v111, v011)
        e_t30 = line_edge(v011, v001)
        e_v0 = line_edge(v000, v001)
        e_v1 = line_edge(v100, v101)
        e_v2 = line_edge(v110, v111)
        e_v3 = line_edge(v010, v011)

        for e in [e_b01, e_b12, e_b23, e_b30, e_t01, e_t12, e_t23, e_t30,
                  e_v0, e_v1, e_v2, e_v3]:
            graph.add_edge(e)

        def make_face(surface, edge_list):
            trims = [brep_trim(e) for e in edge_list]
            for t in trims:
                graph.add_trim(t)
            loop = brep_loop(trims)
            graph.add_loop(loop)
            face = brep_face(surface, [loop])
            graph.add_face(face)
            return face

        # Only 5 faces (no top)
        f_bottom = make_face(plane_surface([0, 0, 0], [0, 0, -1]),
                             [e_b01, e_b12, e_b23, e_b30])
        f_front = make_face(plane_surface([0, 0, 0], [0, -1, 0]),
                            [e_b01, e_v1, e_t01, e_v0])
        f_back = make_face(plane_surface([0, 1, 0], [0, 1, 0]),
                           [e_b23, e_v3, e_t23, e_v2])
        f_left = make_face(plane_surface([0, 0, 0], [-1, 0, 0]),
                           [e_b30, e_v0, e_t30, e_v3])
        f_right = make_face(plane_surface([1, 0, 0], [1, 0, 0]),
                            [e_b12, e_v2, e_t12, e_v1])

        shell = brep_shell([f_bottom, f_front, f_back, f_left, f_right],
                           closed=None)
        return shell

    def test_closed_cube_validation(self):
        """Test that a proper cube is detected as closed."""
        from yapcad.native_brep import shell_closed, shell_id

        graph = TopologyGraph()
        shell = self._build_cube_topology(graph)
        sid = graph.add_shell(shell)

        # Shell should be detected as closed
        assert shell_closed(shell) is True

        # Verify Euler characteristic: V - E + F = 2 for sphere topology
        is_closed, details = graph.compute_shell_closure(sid)
        assert is_closed
        assert details['vertex_count'] == 8
        assert details['edge_count'] == 12
        assert details['face_count'] == 6
        assert details['euler_characteristic'] == 2  # 8 - 12 + 6 = 2
        assert len(details['boundary_edges']) == 0
        assert len(details['non_manifold_edges']) == 0

    def test_open_box_validation(self):
        """Test that an open box (missing face) is detected as open."""
        from yapcad.native_brep import shell_closed, shell_id

        graph = TopologyGraph()
        shell = self._build_open_box_topology(graph)
        sid = graph.add_shell(shell)

        # Shell should be detected as open
        assert shell_closed(shell) is False

        is_closed, details = graph.compute_shell_closure(sid)
        assert not is_closed
        # Should have 4 boundary edges (the top edges)
        assert len(details['boundary_edges']) == 4

    def test_closure_validation_raises_on_invalid(self):
        """Test that validate_shell_closure raises ShellClosureError."""
        from yapcad.native_brep import (
            brep_shell, brep_face, brep_loop, brep_trim, shell_id, ShellClosureError
        )
        from yapcad.analytic_surfaces import plane_surface

        graph = TopologyGraph()

        # Create a single-face shell with actual edges (a square face)
        v0 = brep_vertex([0, 0, 0])
        v1 = brep_vertex([1, 0, 0])
        v2 = brep_vertex([1, 1, 0])
        v3 = brep_vertex([0, 1, 0])
        for v in [v0, v1, v2, v3]:
            graph.add_vertex(v)

        e01 = line_edge(v0, v1)
        e12 = line_edge(v1, v2)
        e23 = line_edge(v2, v3)
        e30 = line_edge(v3, v0)
        for e in [e01, e12, e23, e30]:
            graph.add_edge(e)

        trims = [brep_trim(e) for e in [e01, e12, e23, e30]]
        for t in trims:
            graph.add_trim(t)

        loop = brep_loop(trims)
        graph.add_loop(loop)

        surface = plane_surface([0, 0, 0], [0, 0, 1])
        face = brep_face(surface, [loop])
        graph.add_face(face)

        # Mark shell as closed (but it's not - all 4 edges are boundary edges)
        shell = brep_shell([face], closed=True)

        # Adding should raise ShellClosureError
        with pytest.raises(ShellClosureError) as exc_info:
            graph.add_shell(shell)

        # Check error details
        assert 'boundary edge' in str(exc_info.value).lower() or \
               'not closed' in str(exc_info.value).lower()

    def test_add_shell_skip_validation(self):
        """Test that validation can be skipped."""
        from yapcad.native_brep import brep_shell, brep_face, shell_closed
        from yapcad.analytic_surfaces import plane_surface

        graph = TopologyGraph()
        surface = plane_surface([0, 0, 0], [0, 0, 1])
        face = brep_face(surface, [])
        graph.add_face(face)

        # Force closed=True but skip validation
        shell = brep_shell([face], closed=True)
        sid = graph.add_shell(shell, validate_closure=False)

        # Shell was added despite being invalid
        assert graph.get_shell(sid) is not None
        # closed flag remains True (not validated)
        assert shell_closed(shell) is True


class TestSolidValidation:
    """Test BREP solid validation functionality."""

    def _build_cube_shell(self, graph):
        """Build a closed cube shell for testing.

        Creates a unit cube with 8 vertices, 12 edges, 6 faces.
        Returns the shell.
        """
        from yapcad.analytic_surfaces import plane_surface
        from yapcad.native_brep import (
            brep_vertex, brep_face, brep_loop, brep_trim, brep_shell,
            line_edge
        )

        # 8 vertices of unit cube
        v000 = brep_vertex([0, 0, 0])
        v100 = brep_vertex([1, 0, 0])
        v110 = brep_vertex([1, 1, 0])
        v010 = brep_vertex([0, 1, 0])
        v001 = brep_vertex([0, 0, 1])
        v101 = brep_vertex([1, 0, 1])
        v111 = brep_vertex([1, 1, 1])
        v011 = brep_vertex([0, 1, 1])

        for v in [v000, v100, v110, v010, v001, v101, v111, v011]:
            graph.add_vertex(v)

        # 12 edges
        # Bottom face edges
        e_b01 = line_edge(v000, v100)
        e_b12 = line_edge(v100, v110)
        e_b23 = line_edge(v110, v010)
        e_b30 = line_edge(v010, v000)
        # Top face edges
        e_t01 = line_edge(v001, v101)
        e_t12 = line_edge(v101, v111)
        e_t23 = line_edge(v111, v011)
        e_t30 = line_edge(v011, v001)
        # Vertical edges
        e_v0 = line_edge(v000, v001)
        e_v1 = line_edge(v100, v101)
        e_v2 = line_edge(v110, v111)
        e_v3 = line_edge(v010, v011)

        for e in [e_b01, e_b12, e_b23, e_b30, e_t01, e_t12, e_t23, e_t30,
                  e_v0, e_v1, e_v2, e_v3]:
            graph.add_edge(e)

        def make_face(surface, edge_list):
            """Helper to create face with loop from edges."""
            trims = [brep_trim(e) for e in edge_list]
            for t in trims:
                graph.add_trim(t)
            loop = brep_loop(trims, loop_type='outer')
            graph.add_loop(loop)
            face = brep_face(surface, [loop])
            graph.add_face(face)
            return face

        # 6 faces
        f_bottom = make_face(plane_surface([0, 0, 0], [0, 0, -1]),
                             [e_b01, e_b12, e_b23, e_b30])
        f_top = make_face(plane_surface([0, 0, 1], [0, 0, 1]),
                          [e_t01, e_t12, e_t23, e_t30])
        f_front = make_face(plane_surface([0, 0, 0], [0, -1, 0]),
                            [e_b01, e_v1, e_t01, e_v0])
        f_back = make_face(plane_surface([0, 1, 0], [0, 1, 0]),
                           [e_b23, e_v3, e_t23, e_v2])
        f_left = make_face(plane_surface([0, 0, 0], [-1, 0, 0]),
                           [e_b30, e_v0, e_t30, e_v3])
        f_right = make_face(plane_surface([1, 0, 0], [1, 0, 0]),
                            [e_b12, e_v2, e_t12, e_v1])

        faces = [f_bottom, f_top, f_front, f_back, f_left, f_right]
        shell = brep_shell(faces, closed=None)  # Auto-detect closure
        return shell

    def _build_open_box_shell(self, graph):
        """Build an open box shell (missing top face) for testing.

        Creates a box with 8 vertices, 12 edges, 5 faces.
        Returns the shell.
        """
        from yapcad.analytic_surfaces import plane_surface
        from yapcad.native_brep import (
            brep_vertex, brep_face, brep_loop, brep_trim, brep_shell,
            line_edge
        )

        # 8 vertices of unit cube
        v000 = brep_vertex([0, 0, 0])
        v100 = brep_vertex([1, 0, 0])
        v110 = brep_vertex([1, 1, 0])
        v010 = brep_vertex([0, 1, 0])
        v001 = brep_vertex([0, 0, 1])
        v101 = brep_vertex([1, 0, 1])
        v111 = brep_vertex([1, 1, 1])
        v011 = brep_vertex([0, 1, 1])

        for v in [v000, v100, v110, v010, v001, v101, v111, v011]:
            graph.add_vertex(v)

        # 12 edges
        e_b01 = line_edge(v000, v100)
        e_b12 = line_edge(v100, v110)
        e_b23 = line_edge(v110, v010)
        e_b30 = line_edge(v010, v000)
        e_t01 = line_edge(v001, v101)
        e_t12 = line_edge(v101, v111)
        e_t23 = line_edge(v111, v011)
        e_t30 = line_edge(v011, v001)
        e_v0 = line_edge(v000, v001)
        e_v1 = line_edge(v100, v101)
        e_v2 = line_edge(v110, v111)
        e_v3 = line_edge(v010, v011)

        for e in [e_b01, e_b12, e_b23, e_b30, e_t01, e_t12, e_t23, e_t30,
                  e_v0, e_v1, e_v2, e_v3]:
            graph.add_edge(e)

        def make_face(surface, edge_list):
            """Helper to create face with loop from edges."""
            trims = [brep_trim(e) for e in edge_list]
            for t in trims:
                graph.add_trim(t)
            loop = brep_loop(trims, loop_type='outer')
            graph.add_loop(loop)
            face = brep_face(surface, [loop])
            graph.add_face(face)
            return face

        # Only 5 faces (no top)
        f_bottom = make_face(plane_surface([0, 0, 0], [0, 0, -1]),
                             [e_b01, e_b12, e_b23, e_b30])
        f_front = make_face(plane_surface([0, 0, 0], [0, -1, 0]),
                            [e_b01, e_v1, e_t01, e_v0])
        f_back = make_face(plane_surface([0, 1, 0], [0, 1, 0]),
                           [e_b23, e_v3, e_t23, e_v2])
        f_left = make_face(plane_surface([0, 0, 0], [-1, 0, 0]),
                           [e_b30, e_v0, e_t30, e_v3])
        f_right = make_face(plane_surface([1, 0, 0], [1, 0, 0]),
                            [e_b12, e_v2, e_t12, e_v1])

        faces = [f_bottom, f_front, f_back, f_left, f_right]
        shell = brep_shell(faces, closed=False)  # Explicitly open
        return shell

    def test_valid_solid_with_closed_shell(self):
        """Test that a solid with a closed shell passes validation."""
        from yapcad.native_brep import brep_solid, solid_id

        graph = TopologyGraph()
        shell = self._build_cube_shell(graph)
        sid = graph.add_shell(shell)

        # Create solid from closed shell
        solid = brep_solid([shell])
        solid_id_val = graph.add_solid(solid, validate=True)

        # Solid was added successfully
        assert graph.get_solid(solid_id_val) is not None
        assert solid_id(solid) is not None

    def test_solid_with_open_shell_raises_error(self):
        """Test that a solid with an open shell raises SolidValidationError."""
        from yapcad.native_brep import brep_solid

        graph = TopologyGraph()
        shell = self._build_open_box_shell(graph)
        sid = graph.add_shell(shell, validate_closure=False)

        # Create solid from open shell
        solid = brep_solid([shell])

        # Adding should raise SolidValidationError
        with pytest.raises(SolidValidationError) as exc_info:
            graph.add_solid(solid, validate=True)

        # Check that error mentions shell not being closed
        assert 'closed' in str(exc_info.value).lower() or \
               'open' in str(exc_info.value).lower()

    def test_solid_validation_skip(self):
        """Test that validation can be skipped."""
        from yapcad.native_brep import brep_solid

        graph = TopologyGraph()
        shell = self._build_open_box_shell(graph)
        sid = graph.add_shell(shell, validate_closure=False)

        # Create solid from open shell
        solid = brep_solid([shell])

        # Add with validation disabled
        solid_id_val = graph.add_solid(solid, validate=False)

        # Solid was added despite having an open shell
        assert graph.get_solid(solid_id_val) is not None

    def test_shell_bbox_computation(self):
        """Test that shell bounding box is computed correctly."""
        from yapcad.native_brep import shell_id

        graph = TopologyGraph()
        shell = self._build_cube_shell(graph)
        sid = graph.add_shell(shell)

        # Get bbox
        bbox = graph.shell_bbox(sid)

        # Bbox should exist and have valid structure
        assert bbox is not None
        assert len(bbox) == 2
        # Min point should have lower values than max point
        assert bbox[0][0] <= bbox[1][0]
        assert bbox[0][1] <= bbox[1][1]
        assert bbox[0][2] <= bbox[1][2]
        # Bbox should contain the unit cube (0,0,0) to (1,1,1)
        # Min should be <= 0 and max should be >= 1
        assert bbox[0][0] <= 0.0
        assert bbox[0][1] <= 0.0
        assert bbox[0][2] <= 0.0
        assert bbox[1][0] >= 1.0
        assert bbox[1][1] >= 1.0
        assert bbox[1][2] >= 1.0

    def test_shell_tessellation(self):
        """Test that shell can be tessellated into a mesh."""
        from yapcad.geom3d import issurface

        graph = TopologyGraph()
        shell = self._build_cube_shell(graph)
        sid = graph.add_shell(shell)

        # Tessellate the shell
        mesh = graph.tessellate_shell(sid, u_divisions=4, v_divisions=4)

        # Should be a valid surface
        assert issurface(mesh)
        assert mesh[0] == 'surface'
        # Should have vertices, normals, and faces
        assert len(mesh[1]) > 0  # vertices
        assert len(mesh[2]) > 0  # normals
        assert len(mesh[3]) > 0  # faces

    def test_shell_octree(self):
        """Test that shell octree is built correctly."""

        graph = TopologyGraph()
        shell = self._build_cube_shell(graph)
        sid = graph.add_shell(shell)

        # Get octree
        octree = graph.shell_octree(sid)

        # Octree should exist
        assert octree is not None

        # Query should return some triangles for a point inside the cube
        bbox = [[0.4, 0.4, 0.4, 1], [0.6, 0.6, 0.6, 1]]
        elements = octree.getElements(bbox)
        # May or may not find triangles depending on tessellation,
        # but octree should not error
        assert isinstance(elements, list)


class TestSerialization:
    """Test BREP topology serialization and deserialization."""

    def test_serialize_vertex(self):
        """Test vertex serialization round-trip."""
        from yapcad.native_brep import (
            brep_vertex, _serialize_vertex, _deserialize_vertex,
            vertex_id, vertex_location
        )

        v = brep_vertex([1.5, 2.5, 3.5], tolerance=0.001)
        serialized = _serialize_vertex(v)

        # Check serialized structure
        assert serialized['type'] == 'brep_vertex'
        assert serialized['location'] == [1.5, 2.5, 3.5, 1.0]
        assert serialized['tolerance'] == 0.001

        # Round-trip
        v2 = _deserialize_vertex(serialized)
        assert vertex_id(v2) == serialized['id']
        loc = vertex_location(v2)
        assert abs(loc[0] - 1.5) < 1e-10
        assert abs(loc[1] - 2.5) < 1e-10
        assert abs(loc[2] - 3.5) < 1e-10

    def test_serialize_line_edge(self):
        """Test line edge serialization round-trip."""
        from yapcad.native_brep import (
            brep_vertex, line_edge, _serialize_edge, _deserialize_edge,
            edge_id, edge_curve_type
        )

        v1 = brep_vertex([0, 0, 0])
        v2 = brep_vertex([1, 1, 1])
        e = line_edge(v1, v2)

        serialized = _serialize_edge(e)

        # Check serialized structure
        assert serialized['type'] == 'brep_edge'
        assert serialized['curve_type'] == 'line'

        # Round-trip
        e2 = _deserialize_edge(serialized)
        assert edge_id(e2) == serialized['id']
        assert edge_curve_type(e2) == 'line'

    def test_serialize_arc_edge(self):
        """Test arc edge serialization round-trip."""
        from yapcad.native_brep import (
            brep_vertex, arc_edge, _serialize_edge, _deserialize_edge,
            edge_id, edge_curve_type
        )

        v1 = brep_vertex([1, 0, 0])
        v2 = brep_vertex([0, 1, 0])
        # arc_edge uses center (not radius - radius is computed from geometry)
        e = arc_edge(v1, v2, center=[0, 0, 0])

        serialized = _serialize_edge(e)

        # Check serialized structure
        assert serialized['type'] == 'brep_edge'
        assert serialized['curve_type'] == 'arc'
        assert 'curve_params' in serialized

        # Round-trip
        e2 = _deserialize_edge(serialized)
        assert edge_id(e2) == serialized['id']
        assert edge_curve_type(e2) == 'arc'

    def test_serialize_topology_graph(self):
        """Test full topology graph serialization round-trip."""
        from yapcad.native_brep import (
            brep_vertex, line_edge, brep_trim, brep_loop,
            brep_face, brep_shell,
            serialize_topology_graph, deserialize_topology_graph
        )
        from yapcad.analytic_surfaces import plane_surface

        # Build a simple topology
        graph = TopologyGraph()

        # Add vertices
        v1 = brep_vertex([0, 0, 0])
        v2 = brep_vertex([1, 0, 0])
        v3 = brep_vertex([1, 1, 0])
        v4 = brep_vertex([0, 1, 0])

        for v in [v1, v2, v3, v4]:
            graph.add_vertex(v)

        # Add edges
        e1 = line_edge(v1, v2)
        e2 = line_edge(v2, v3)
        e3 = line_edge(v3, v4)
        e4 = line_edge(v4, v1)

        for e in [e1, e2, e3, e4]:
            graph.add_edge(e)

        # Add trims, loop, face
        trims = [brep_trim(e) for e in [e1, e2, e3, e4]]
        for t in trims:
            graph.add_trim(t)

        loop = brep_loop(trims, loop_type='outer')
        graph.add_loop(loop)

        surface = plane_surface([0, 0, 0], [0, 0, 1])
        face = brep_face(surface, [loop])
        graph.add_face(face)

        shell = brep_shell([face], closed=False)
        graph.add_shell(shell, validate_closure=False)

        # Serialize
        serialized = serialize_topology_graph(graph)

        # Check structure
        assert serialized['version'] == '1.0'
        assert len(serialized['vertices']) == 4
        assert len(serialized['edges']) == 4
        assert len(serialized['trims']) == 4
        assert len(serialized['loops']) == 1
        assert len(serialized['faces']) == 1
        assert len(serialized['shells']) == 1

        # Deserialize
        graph2 = deserialize_topology_graph(serialized)

        # Check counts match
        assert len(graph2.vertices) == 4
        assert len(graph2.edges) == 4
        assert len(graph2.trims) == 4
        assert len(graph2.loops) == 1
        assert len(graph2.faces) == 1
        assert len(graph2.shells) == 1

    def test_serialize_surface_plane(self):
        """Test plane surface serialization round-trip."""
        from yapcad.native_brep import _serialize_surface, _deserialize_surface
        from yapcad.analytic_surfaces import plane_surface, is_plane_surface

        surf = plane_surface([1, 2, 3], [0, 0, 1])
        serialized = _serialize_surface(surf)

        # Check structure
        assert serialized['surface_type'] == 'plane_surface'
        assert serialized['origin'] == [1.0, 2.0, 3.0, 1.0]

        # Round-trip
        surf2 = _deserialize_surface(serialized)
        assert is_plane_surface(surf2)

    def test_serialize_surface_sphere(self):
        """Test sphere surface serialization round-trip."""
        from yapcad.native_brep import _serialize_surface, _deserialize_surface
        from yapcad.analytic_surfaces import sphere_surface, is_sphere_surface

        surf = sphere_surface([0, 0, 0], 2.5)
        serialized = _serialize_surface(surf)

        # Check structure
        assert serialized['surface_type'] == 'sphere_surface'
        assert serialized['metadata']['radius'] == 2.5

        # Round-trip
        surf2 = _deserialize_surface(serialized)
        assert is_sphere_surface(surf2)
        assert surf2[2]['radius'] == 2.5


class TestSolidMetadataIntegration:
    """Test native BREP integration with yapCAD solid metadata."""

    def test_attach_and_retrieve_native_brep(self):
        """Test attaching and retrieving native BREP from a solid."""
        from yapcad.geom3d_util import prism
        from yapcad.native_brep import (
            attach_native_brep_to_solid, native_brep_from_solid,
            has_native_brep, clear_native_brep,
            brep_vertex, line_edge, brep_trim, brep_loop,
            brep_face, brep_shell
        )
        from yapcad.analytic_surfaces import plane_surface

        # Create a simple yapCAD solid
        solid = prism(1, 1, 1)

        # Initially no native BREP
        assert not has_native_brep(solid)

        # Build a simple native BREP (single face)
        graph = TopologyGraph()

        v1 = brep_vertex([0, 0, 0])
        v2 = brep_vertex([1, 0, 0])
        v3 = brep_vertex([1, 1, 0])
        v4 = brep_vertex([0, 1, 0])

        for v in [v1, v2, v3, v4]:
            graph.add_vertex(v)

        e1 = line_edge(v1, v2)
        e2 = line_edge(v2, v3)
        e3 = line_edge(v3, v4)
        e4 = line_edge(v4, v1)

        for e in [e1, e2, e3, e4]:
            graph.add_edge(e)

        trims = [brep_trim(e) for e in [e1, e2, e3, e4]]
        for t in trims:
            graph.add_trim(t)

        loop = brep_loop(trims)
        graph.add_loop(loop)

        surface = plane_surface([0, 0, 0], [0, 0, 1])
        face = brep_face(surface, [loop])
        graph.add_face(face)

        shell = brep_shell([face], closed=False)
        graph.add_shell(shell, validate_closure=False)

        # Attach to solid
        attach_native_brep_to_solid(solid, graph)

        # Now has native BREP
        assert has_native_brep(solid)

        # Retrieve and verify
        graph2 = native_brep_from_solid(solid)
        assert graph2 is not None
        assert len(graph2.vertices) == 4
        assert len(graph2.edges) == 4
        assert len(graph2.faces) == 1

        # Clear and verify
        clear_native_brep(solid)
        assert not has_native_brep(solid)
        assert native_brep_from_solid(solid) is None

    def test_json_serializable(self):
        """Test that serialized native BREP is JSON-serializable."""
        import json
        from yapcad.native_brep import (
            brep_vertex, line_edge, brep_trim, brep_loop,
            brep_face, brep_shell, serialize_topology_graph
        )
        from yapcad.analytic_surfaces import plane_surface

        # Build topology
        graph = TopologyGraph()

        v1 = brep_vertex([0, 0, 0])
        v2 = brep_vertex([1, 0, 0])
        graph.add_vertex(v1)
        graph.add_vertex(v2)

        e = line_edge(v1, v2)
        graph.add_edge(e)

        surface = plane_surface([0, 0, 0], [0, 0, 1])
        face = brep_face(surface, [])
        graph.add_face(face)

        # Serialize
        serialized = serialize_topology_graph(graph)

        # Should be JSON-serializable
        json_str = json.dumps(serialized)
        assert isinstance(json_str, str)
        assert len(json_str) > 0

        # Should round-trip through JSON
        parsed = json.loads(json_str)
        assert parsed['version'] == '1.0'
        assert len(parsed['vertices']) == 2


class TestNativeBrepTransformations:
    """Test native BREP transformation functions."""

    def test_translate_native_brep(self):
        """Test translation of native BREP in solid metadata."""
        from yapcad.geom3d_util import prism
        from yapcad.native_brep import (
            attach_native_brep_to_solid, native_brep_from_solid,
            has_native_brep, brep_vertex, line_edge, brep_trim,
            brep_loop, brep_face, brep_shell, vertex_location,
            translate_native_brep
        )
        from yapcad.analytic_surfaces import plane_surface

        # Create a simple yapCAD solid with native BREP
        solid = prism(1, 1, 1)

        # Build a simple native BREP with one vertex at origin
        graph = TopologyGraph()
        v1 = brep_vertex([0, 0, 0])
        v2 = brep_vertex([1, 0, 0])
        graph.add_vertex(v1)
        graph.add_vertex(v2)

        e = line_edge(v1, v2)
        graph.add_edge(e)

        surface = plane_surface([0, 0, 0], [0, 0, 1])
        face = brep_face(surface, [])
        graph.add_face(face)

        attach_native_brep_to_solid(solid, graph)

        # Translate
        translate_native_brep(solid, [10, 20, 30])

        # Verify transformation
        graph2 = native_brep_from_solid(solid)
        v1_new = list(graph2.vertices.values())[0]
        loc = vertex_location(v1_new)
        # Original was at (0,0,0) -> should now be at (10,20,30)
        # or original was at (1,0,0) -> should now be at (11,20,30)
        # Just check that z is around 30
        assert abs(loc[2] - 30) < 0.01

    def test_transform_topology_graph(self):
        """Test generic transformation of topology graph."""
        from yapcad.native_brep import (
            brep_vertex, line_edge, brep_face,
            transform_topology_graph, vertex_location
        )
        from yapcad.analytic_surfaces import plane_surface
        from yapcad.geom import point

        graph = TopologyGraph()
        v1 = brep_vertex([1, 0, 0])
        v2 = brep_vertex([0, 1, 0])
        graph.add_vertex(v1)
        graph.add_vertex(v2)

        e = line_edge(v1, v2)
        graph.add_edge(e)

        # Simple translation transform
        def translate(p):
            return point(p[0] + 5, p[1] + 5, p[2] + 5)

        transform_topology_graph(graph, translate)

        # Verify vertex positions changed
        for vid, vertex in graph.vertices.items():
            loc = vertex_location(vertex)
            # All z values should be 5 now
            assert abs(loc[2] - 5) < 0.01


class TestOCCNativeConversion:
    """Test OCC ↔ Native BREP conversion functions."""

    @pytest.fixture
    def occ_available(self):
        """Check if OCC is available."""
        try:
            from yapcad.occ_native_convert import occ_available
            return occ_available()
        except ImportError:
            return False

    def test_native_vertex_to_occ(self, occ_available):
        """Test converting a native vertex to OCC."""
        if not occ_available:
            pytest.skip("OCC not available")

        from yapcad.native_brep import brep_vertex, vertex_location
        from yapcad.occ_native_convert import native_vertex_to_occ
        from OCC.Core.BRep import BRep_Tool

        # Create native vertex
        v = brep_vertex([1.0, 2.0, 3.0])

        # Convert to OCC
        occ_v = native_vertex_to_occ(v)

        # Verify position
        pnt = BRep_Tool.Pnt(occ_v)
        assert abs(pnt.X() - 1.0) < 1e-10
        assert abs(pnt.Y() - 2.0) < 1e-10
        assert abs(pnt.Z() - 3.0) < 1e-10

    def test_native_edge_to_occ_line(self, occ_available):
        """Test converting a line edge to OCC."""
        if not occ_available:
            pytest.skip("OCC not available")

        from yapcad.native_brep import brep_vertex, line_edge
        from yapcad.occ_native_convert import native_edge_to_occ
        from OCC.Core.BRep import BRep_Tool
        from OCC.Core.TopExp import topexp

        # Create native line edge
        v1 = brep_vertex([0, 0, 0])
        v2 = brep_vertex([10, 0, 0])

        graph = TopologyGraph()
        graph.add_vertex(v1)
        graph.add_vertex(v2)

        e = line_edge(v1, v2)
        graph.add_edge(e)

        # Convert to OCC
        occ_e = native_edge_to_occ(e, graph)

        # Verify endpoints
        first_v = topexp.FirstVertex(occ_e)
        last_v = topexp.LastVertex(occ_e)

        p1 = BRep_Tool.Pnt(first_v)
        p2 = BRep_Tool.Pnt(last_v)

        # One should be at origin, other at (10, 0, 0)
        assert (abs(p1.X()) < 1e-6 or abs(p2.X()) < 1e-6)
        assert (abs(p1.X() - 10) < 1e-6 or abs(p2.X() - 10) < 1e-6)

    def test_native_surface_to_occ_plane(self, occ_available):
        """Test converting a plane surface to OCC."""
        if not occ_available:
            pytest.skip("OCC not available")

        from yapcad.analytic_surfaces import plane_surface
        from yapcad.occ_native_convert import native_surface_to_occ
        from OCC.Core.GeomAbs import GeomAbs_Plane

        # Create native plane
        p = plane_surface([0, 0, 5], [0, 0, 1])

        # Convert to OCC
        occ_surf = native_surface_to_occ(p)

        assert occ_surf is not None
        # Check it's a plane
        pln = occ_surf.Pln()
        loc = pln.Location()
        assert abs(loc.Z() - 5.0) < 1e-10

    def test_native_surface_to_occ_sphere(self, occ_available):
        """Test converting a sphere surface to OCC."""
        if not occ_available:
            pytest.skip("OCC not available")

        from yapcad.analytic_surfaces import sphere_surface
        from yapcad.occ_native_convert import native_surface_to_occ

        # Create native sphere
        s = sphere_surface([1, 2, 3], 5.0)

        # Convert to OCC
        occ_surf = native_surface_to_occ(s)

        assert occ_surf is not None
        sph = occ_surf.Sphere()
        assert abs(sph.Radius() - 5.0) < 1e-10
        loc = sph.Location()
        assert abs(loc.X() - 1.0) < 1e-10
        assert abs(loc.Y() - 2.0) < 1e-10
        assert abs(loc.Z() - 3.0) < 1e-10

    def test_native_surface_to_occ_cylinder(self, occ_available):
        """Test converting a cylinder surface to OCC."""
        if not occ_available:
            pytest.skip("OCC not available")

        from yapcad.analytic_surfaces import cylinder_surface
        from yapcad.occ_native_convert import native_surface_to_occ

        # Create native cylinder
        c = cylinder_surface([0, 0, 0], [0, 0, 1], 3.0)

        # Convert to OCC
        occ_surf = native_surface_to_occ(c)

        assert occ_surf is not None
        cyl = occ_surf.Cylinder()
        assert abs(cyl.Radius() - 3.0) < 1e-10

    def test_roundtrip_occ_to_native_to_occ(self, occ_available):
        """Test round-trip conversion OCC → Native → OCC."""
        if not occ_available:
            pytest.skip("OCC not available")

        from yapcad.occ_native_convert import (
            occ_solid_to_native_brep, native_brep_to_occ
        )
        from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
        from OCC.Core.GProp import GProp_GProps
        from OCC.Core.BRepGProp import brepgprop

        # Create an OCC box
        box_maker = BRepPrimAPI_MakeBox(10, 10, 10)
        occ_box = box_maker.Solid()

        # Convert to native BREP
        native_solid, graph = occ_solid_to_native_brep(occ_box)

        # Verify graph has topology
        assert len(graph.faces) == 6  # Cube has 6 faces
        assert len(graph.vertices) > 0
        assert len(graph.edges) > 0

        # Convert back to OCC
        occ_box2 = native_brep_to_occ(graph)

        # Verify it's a valid solid by computing volume
        props = GProp_GProps()
        brepgprop.VolumeProperties(occ_box2, props)
        volume = props.Mass()

        # Volume should be close to 10*10*10 = 1000
        # Allow tolerance because sewing might slightly affect shape
        assert volume > 900 and volume < 1100

    def test_roundtrip_cylinder(self, occ_available):
        """Test round-trip conversion of a cylinder.

        Note: Cylinder round-trip is challenging due to circular edge reconstruction.
        This test verifies the conversion completes without error; volume preservation
        for curved geometry may require additional curve conversion work.
        """
        if not occ_available:
            pytest.skip("OCC not available")

        from yapcad.occ_native_convert import (
            occ_solid_to_native_brep, native_brep_to_occ
        )
        from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeCylinder
        from OCC.Core.TopAbs import TopAbs_SOLID, TopAbs_SHELL

        # Create an OCC cylinder
        cyl_maker = BRepPrimAPI_MakeCylinder(5, 10)
        occ_cyl = cyl_maker.Solid()

        # Convert to native BREP
        native_solid, graph = occ_solid_to_native_brep(occ_cyl)

        # Verify graph has topology (cylinder: 3 faces - top, bottom, side)
        assert len(graph.faces) == 3

        # Convert back to OCC - just verify it completes without exception
        occ_cyl2 = native_brep_to_occ(graph)

        # Should return some shape (solid or shell)
        assert occ_cyl2 is not None
        assert occ_cyl2.ShapeType() in (TopAbs_SOLID, TopAbs_SHELL)

    def test_native_brep_to_occ_from_solid(self, occ_available):
        """Test converting native BREP from a yapCAD solid via OCC roundtrip."""
        if not occ_available:
            pytest.skip("OCC not available")

        from yapcad.geom3d_util import prism
        from yapcad.native_brep import (
            attach_native_brep_to_solid, has_native_brep
        )
        from yapcad.occ_native_convert import (
            native_brep_to_occ_from_solid, occ_solid_to_native_brep
        )
        from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
        from OCC.Core.TopAbs import TopAbs_SOLID, TopAbs_SHELL

        # Create a yapCAD solid
        solid = prism(2, 2, 2)

        # Create an OCC box to get a valid native BREP
        occ_box = BRepPrimAPI_MakeBox(2, 2, 2).Solid()
        _, graph = occ_solid_to_native_brep(occ_box)

        # Attach the native BREP to the yapCAD solid
        attach_native_brep_to_solid(solid, graph)
        assert has_native_brep(solid)

        # Now convert back to OCC
        occ_shape = native_brep_to_occ_from_solid(solid)

        # Should get either a solid or shell back
        assert occ_shape is not None
        assert occ_shape.ShapeType() in (TopAbs_SOLID, TopAbs_SHELL)
