"""Tests for Phase 3 DSL builtins: text_solid, engrave_text, text_width."""
import pytest
from yapcad.dsl.runtime.builtins import get_builtin_registry, call_builtin
from yapcad.dsl.runtime.values import Value, float_val, string_val, solid_val, vector_val
from yapcad.dsl.types import FLOAT, SOLID, STRING
from yapcad.geom3d import issolid, solidbbox
from yapcad.geom3d_util import prism
from yapcad.geom import point, vect


@pytest.fixture
def registry():
    return get_builtin_registry()


class TestTextSolid:
    def test_exists(self, registry):
        assert registry.get_function("text_solid") is not None

    def test_basic_text(self):
        result = call_builtin("text_solid", [string_val("HI"), float_val(5.0), float_val(1.0), float_val(1.0)])
        assert result.type == SOLID
        assert issolid(result.data)
        assert len(result.data[1]) > 0

    def test_text_with_height(self):
        result = call_builtin("text_solid", [string_val("AB"), float_val(10.0), float_val(2.0), float_val(1.5)])
        assert issolid(result.data)
        bbox = solidbbox(result.data)
        assert bbox is not None
        y_extent = bbox[1][1] - bbox[0][1]
        assert 8.0 < y_extent < 12.0, f"Expected ~10mm height, got {y_extent}"

    def test_single_char(self):
        result = call_builtin("text_solid", [string_val("A"), float_val(5.0), float_val(1.0), float_val(1.0)])
        assert issolid(result.data)
        assert len(result.data[1]) > 0

    def test_text_depth(self):
        thin = call_builtin("text_solid", [string_val("X"), float_val(5.0), float_val(1.0), float_val(1.0)])
        thick = call_builtin("text_solid", [string_val("X"), float_val(5.0), float_val(5.0), float_val(1.0)])
        thin_bbox = solidbbox(thin.data)
        thick_bbox = solidbbox(thick.data)
        thin_z = thin_bbox[1][2] - thin_bbox[0][2]
        thick_z = thick_bbox[1][2] - thick_bbox[0][2]
        assert thick_z > thin_z

    def test_empty_text(self):
        result = call_builtin("text_solid", [string_val(""), float_val(5.0), float_val(1.0), float_val(1.0)])
        assert result.type == SOLID


class TestTextWidth:
    def test_exists(self, registry):
        assert registry.get_function("text_width") is not None

    def test_basic_width(self):
        result = call_builtin("text_width", [string_val("HELLO"), float_val(5.0), float_val(1.0)])
        assert result.type == FLOAT
        assert result.data > 0

    def test_longer_text_wider(self):
        w1 = call_builtin("text_width", [string_val("A"), float_val(5.0), float_val(1.0)])
        w5 = call_builtin("text_width", [string_val("ABCDE"), float_val(5.0), float_val(1.0)])
        assert w5.data > w1.data

    def test_height_scales_width(self):
        w_small = call_builtin("text_width", [string_val("AB"), float_val(5.0), float_val(1.0)])
        w_large = call_builtin("text_width", [string_val("AB"), float_val(10.0), float_val(1.0)])
        ratio = w_large.data / w_small.data
        assert 1.8 < ratio < 2.2, f"Expected ~2x ratio, got {ratio}"

    def test_empty_width(self):
        result = call_builtin("text_width", [string_val(""), float_val(5.0), float_val(1.0)])
        assert result.data == 0.0


class TestEngraveText:
    def test_exists(self, registry):
        assert registry.get_function("engrave_text") is not None

    def test_basic_engrave(self):
        plate = prism(50, 30, 5)
        pos = point(5, 10, 2.5)
        normal = vect(0, 0, 1, 0)
        result = call_builtin("engrave_text", [
            solid_val(plate), string_val("V1"),
            vector_val(pos), vector_val(normal),
            float_val(3.0), float_val(0.5), float_val(1.0)
        ])
        assert result.type == SOLID
        assert issolid(result.data)

    def test_engrave_preserves_structure(self):
        plate = prism(50, 30, 5)
        result = call_builtin("engrave_text", [
            solid_val(plate), string_val("AB"),
            vector_val(point(0, 0, 2.5)), vector_val(vect(0, 0, 1, 0)),
            float_val(3.0), float_val(0.5), float_val(1.0)
        ])
        assert len(result.data[1]) > 0
