"""Tests for Phase 4 DSL builtins: path3d_eval, path3d_length, split_solid."""
import pytest
import math
from yapcad.dsl.runtime.builtins import get_builtin_registry, call_builtin
from yapcad.dsl.runtime.values import (
    Value, float_val, point_val, vector_val, solid_val, path3d_val, wrap_value,
)
from yapcad.dsl.types import FLOAT, SOLID, POINT3D, PATH3D
from yapcad.geom import point, vect
from yapcad.geom3d import issolid
from yapcad.geom3d_util import prism


@pytest.fixture
def registry():
    return get_builtin_registry()


def make_line_path(start, end):
    """Helper: create a path3d Value with a single line segment."""
    return path3d_val({
        'type': 'path3d',
        'segments': [{
            'type': 'line',
            'start': [start[0], start[1], start[2]],
            'end': [end[0], end[1], end[2]],
        }]
    })


def make_two_segment_path():
    """Helper: path from (0,0,0) -> (10,0,0) -> (10,10,0)."""
    return path3d_val({
        'type': 'path3d',
        'segments': [
            {'type': 'line', 'start': [0, 0, 0], 'end': [10, 0, 0]},
            {'type': 'line', 'start': [10, 0, 0], 'end': [10, 10, 0]},
        ]
    })


class TestPath3dEval:
    def test_exists(self, registry):
        assert registry.get_function("path3d_eval") is not None

    def test_eval_start(self):
        path = make_line_path([0, 0, 0], [10, 0, 0])
        result = call_builtin("path3d_eval", [path, float_val(0.0)])
        assert result.type == POINT3D
        p = result.data
        assert abs(p[0] - 0.0) < 0.01
        assert abs(p[1] - 0.0) < 0.01

    def test_eval_end(self):
        path = make_line_path([0, 0, 0], [10, 0, 0])
        result = call_builtin("path3d_eval", [path, float_val(1.0)])
        p = result.data
        assert abs(p[0] - 10.0) < 0.01

    def test_eval_midpoint(self):
        path = make_line_path([0, 0, 0], [10, 0, 0])
        result = call_builtin("path3d_eval", [path, float_val(0.5)])
        p = result.data
        assert abs(p[0] - 5.0) < 0.01

    def test_eval_two_segments(self):
        path = make_two_segment_path()
        # t=0.25 should be midpoint of first segment
        result = call_builtin("path3d_eval", [path, float_val(0.25)])
        p = result.data
        assert abs(p[0] - 5.0) < 0.01
        assert abs(p[1] - 0.0) < 0.01

    def test_eval_second_segment(self):
        path = make_two_segment_path()
        # t=0.75 should be midpoint of second segment
        result = call_builtin("path3d_eval", [path, float_val(0.75)])
        p = result.data
        assert abs(p[0] - 10.0) < 0.01
        assert abs(p[1] - 5.0) < 0.01


class TestPath3dLength:
    def test_exists(self, registry):
        assert registry.get_function("path3d_length") is not None

    def test_simple_length(self):
        path = make_line_path([0, 0, 0], [10, 0, 0])
        result = call_builtin("path3d_length", [path])
        assert result.type == FLOAT
        assert abs(result.data - 10.0) < 0.01

    def test_diagonal_length(self):
        path = make_line_path([0, 0, 0], [3, 4, 0])
        result = call_builtin("path3d_length", [path])
        assert abs(result.data - 5.0) < 0.01

    def test_two_segment_length(self):
        path = make_two_segment_path()
        result = call_builtin("path3d_length", [path])
        assert abs(result.data - 20.0) < 0.01

    def test_3d_length(self):
        path = make_line_path([0, 0, 0], [1, 1, 1])
        result = call_builtin("path3d_length", [path])
        assert abs(result.data - math.sqrt(3)) < 0.01


class TestSplitSolid:
    def test_exists(self, registry):
        assert registry.get_function("split_solid") is not None

    @pytest.mark.skipif(True, reason="Requires OCC (pythonocc-core)")
    def test_split_box(self):
        """Split a box in half — requires OCC."""
        box = prism(20, 10, 10)
        result = call_builtin("split_solid", [
            solid_val(box),
            point_val(point(10, 0, 0)),
            vector_val(vect(1, 0, 0, 0)),
        ])
        # Should return a list of two solids
        assert result.data is not None
