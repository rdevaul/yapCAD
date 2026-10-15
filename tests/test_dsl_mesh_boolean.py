"""Retained booleans agree with the eager mesh path at public boundaries."""

import copy
import json
import os
import unittest
from unittest.mock import patch

import numpy as np
import pytest

from yapcad.dsl import compile_and_run
from yapcad.dsl.runtime.builtins import get_builtin_registry
from yapcad.dsl.runtime.values import solid_val
from yapcad.geom3d import solid_boolean
from yapcad.dsl.runtime.mesh_boolean import ManifoldSolid, native, try_boolean


pytest.importorskip("manifold3d")


def run(body, result_type="solid"):
    result = compile_and_run(
        f"module test\ncommand TEST() -> {result_type}:\n" + body, "TEST", {}
    )
    assert result.success, result.error_message
    return result.geometry


class RuntimeTests(unittest.TestCase):
    def setUp(self):
        # These tests exercise the mesh path even when OCC is installed.
        engine = patch.dict(os.environ, {"YAPCAD_BOOLEAN_ENGINE": "manifold"})
        engine.start()
        self.addCleanup(engine.stop)

    def test_brep_inputs_keep_auto_dispatch(self):
        from yapcad.brep import occ_available, has_brep_data

        if not occ_available():
            self.skipTest("requires OCC")
        with patch.dict(os.environ):
            os.environ.pop("YAPCAD_BOOLEAN_ENGINE", None)
            a = run("    emit box(2.0, 2.0, 2.0)\n")
            self.assertTrue(has_brep_data(a))
            self.assertIsNone(try_boolean([a, a], "union"))
            result = (
                get_builtin_registry()
                .get_function("union")
                .implementation(solid_val(a), solid_val(a))
                .data
            )
            self.assertNotIsInstance(result, ManifoldSolid)
            self.assertTrue(has_brep_data(result))

    def test_retained_and_eager_booleans_agree(self):
        a = run("    emit box(10.0, 10.0, 4.0)\n")
        tools = [
            run(f"    emit translate(cylinder(1.0, 8.0), {x}.0, 0.0, -4.0)\n")
            for x in (-2, 2)
        ]
        registry = get_builtin_registry()
        for operation in ("union", "difference", "intersection"):
            result = (
                registry.get_function(operation)
                .implementation(*map(solid_val, [a, *tools]))
                .data
            )
            self.assertIsInstance(result, ManifoldSolid)
            eager = solid_boolean(
                solid_boolean(a, tools[0], operation, engine="manifold"),
                tools[1],
                operation,
                engine="manifold",
            )
            self.assertAlmostEqual(
                native(result).volume(), native(eager).volume(), places=7
            )
            np.testing.assert_allclose(
                native(result).bounding_box(),
                native(eager).bounding_box(),
                atol=1e-7,
            )

    def test_list_comprehensions_and_accumulation(self):
        expected = run(
            (
                "    emit difference(box(10.0, 10.0, 4.0), cylinder(1.0,"
                " 10.0), translate(cylinder(1.0, 10.0), 3.0, 0.0, 0.0))\n"
            )
        )
        batched = run(
            (
                "    tools: list<solid> = [translate(cylinder(1.0, 10.0)"
                ", 3.0 * i, 0.0, 0.0) for i in range(2)]\n"
                "    emit difference_all(box(10.0, 10.0, 4.0), tools)\n"
            )
        )
        self.assertAlmostEqual(
            native(expected).volume(), native(batched).volume(), places=7
        )
        empty = run(
            (
                "    s: solid = empty_solid()\n"
                "    for i in range(3):\n"
                "        s = union(s, translate(box(1.0, 1.0, 1.0), 2.0 "
                "* i, 0.0, 0.0))\n"
                "    emit s\n"
            )
        )
        self.assertAlmostEqual(native(empty).volume(), 3.0)

    def test_reassignment_and_early_emit(self):
        self.assertEqual(
            run(
                (
                    "    x: float = 1.0\n"
                    "    x = x + 1.0\n"
                    "    for i in range(3):\n"
                    "        x = x + 2.0\n"
                    "        x = x + 3.0\n"
                    "        if i == 1:\n"
                    "            emit x\n"
                    "    emit 999.0\n"
                ),
                "float",
            ),
            12.0,
        )
        wrong = compile_and_run(
            (
                "module t\n"
                "command TEST() -> float:\n"
                "    x: float = 1.0\n"
                '    x = "bad"\n'
                "    emit x\n"
            ),
            "TEST",
            {},
        )
        self.assertFalse(wrong.success)

    def test_mesh_access_and_copy_do_not_leave_stale_geometry(self):
        result = run(
            "    emit union(box(2.0, 2.0, 2.0), box(1.0, 1.0, 1.0))\n"
        )
        result = try_boolean([result], "union")
        cloned = copy.deepcopy(result)
        self.assertIsNotNone(result.body)
        self.assertIsNotNone(cloned.body)
        for surf in result[1]:
            for vertex in surf[1]:
                vertex[0] += 10
        self.assertIsNone(result.body)
        self.assertAlmostEqual(native(result).bounding_box()[0], 9.0)
        self.assertAlmostEqual(native(cloned).bounding_box()[0], -1.0)
        encoded = json.loads(json.dumps(cloned))
        self.assertTrue(encoded[1])
        self.assertAlmostEqual(native(encoded).volume(), 8.0)

    def test_structural_list_mutation_materializes_first(self):
        result = run("    emit box(2.0, 2.0, 2.0)\n")
        retained = try_boolean([result], "union")
        surfaces = retained.pop(1)
        self.assertTrue(surfaces)
        self.assertIsNone(retained.body)
        retained.insert(1, surfaces)
        self.assertAlmostEqual(native(retained).volume(), 8.0)

    def test_native_transforms_match_mesh_transforms(self):
        result = run(
            (
                "    emit rotate(translate(union(box(2.0, 3.0, 4.0), box"
                "(1.0, 1.0, 1.0)), 1.0, 2.0, 3.0), 12.0, 23.0, 34.0)\n"
            )
        )
        eager = run(
            (
                "    emit rotate(translate(box(2.0, 3.0, 4.0), 1.0, 2.0,"
                " 3.0), 12.0, 23.0, 34.0)\n"
            )
        )
        np.testing.assert_allclose(
            native(result).bounding_box(),
            native(eager).bounding_box(),
            atol=1e-7,
        )
        self.assertAlmostEqual(
            native(result).volume(), native(eager).volume(), places=7
        )

    def test_engine_override_falls_back(self):
        with patch.dict(os.environ, {"YAPCAD_BOOLEAN_ENGINE": "native"}):
            result = run(
                "    emit union(box(2.0, 2.0, 2.0), box(1.0, 1.0, 1.0))\n"
            )
        self.assertNotIsInstance(result, ManifoldSolid)

    def test_boolean_helper_condition(self):
        result = compile_and_run(
            (
                "module t\n"
                "command yes() -> bool:\n"
                "    emit true\n"
                "command TEST() -> float:\n"
                "    if true and not not yes():\n"
                "        emit 2.0\n"
                "    emit 1.0\n"
            ),
            "TEST",
            {},
        )
        self.assertTrue(result.success, result.error_message)
        self.assertEqual(result.geometry, 2.0)

    def test_mesh_conversion_happens_at_the_public_boundary(self):
        from yapcad.boolean import manifold_engine

        with patch.object(
            manifold_engine,
            "_from_manifold",
            wraps=manifold_engine._from_manifold,
        ) as convert:
            result = run(
                (
                    "    s: solid = box(10.0, 10.0, 1.0)\n"
                    "    for i in range(12):\n"
                    "        s = difference(s, translate(box(0.3, 0.3, 2.0),"
                    " 0.5 * i - 3.0, 0.0, 0.0))\n"
                    "    emit s\n"
                )
            )
        self.assertEqual(convert.call_count, 1)
        self.assertIs(type(result), list)
        self.assertTrue(result[1])

    def test_three_argument_text(self):
        result = run('    emit text_solid("ABC", 5.0, 1.0)\n')
        self.assertGreater(native(result).volume(), 0)


if __name__ == "__main__":
    unittest.main()
