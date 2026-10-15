"""Simplification must not emit faces its own library calls degenerate.

A collapse used to be accepted on any non-zero area, so it could leave a
sliver of order 1e-13 -- three distinct but collinear vertices.  That costs
nothing in field error, so no other test here objected.  Mesh transport
carries such a face since small faces stopped being dropped, so it no longer
opens a hole; but ``triangle_normal`` still has no normal for it at the
default tolerance, and mesh tools, slicers and STL importers downstream drop
or choke on zero-area facets.

These tests look at what the simplifier emits, not at what a consumer makes
of it, so they stay honest whatever transport tolerance is in force.
"""
import textwrap

import numpy as np
import pytest

from yapcad.dsl import compile_and_run
from yapcad.geom import epsilon
from yapcad.geometry_utils import triangle_normal

PLATE = textwrap.dedent("""
    module plate
    command PLATE() -> solid:
        let plate: solid = box(60.0, 30.0, 6.0)
        let hole: solid = translate(cylinder(5.0, 8.0), 0.0, 0.0, -4.0)
        let left: solid = translate(cylinder(2.25, 8.0), -21.0, 0.0, -4.0)
        let right: solid = translate(cylinder(2.25, 8.0), 21.0, 0.0, -4.0)
        emit difference(plate, hole, left, right)
""")


def _faces(solid):
    """Every emitted face as an array of its three vertices."""
    out = []
    for surf in solid[1]:
        verts, faces = surf[1], surf[3]
        for face in faces:
            out.append(np.array([verts[i][:3] for i in face], dtype=float))
    return out


def _cross_magnitudes(solid):
    tris = _faces(solid)
    assert tris, "the solid carries no faces"
    return np.array([
        float(np.linalg.norm(np.cross(t[1] - t[0], t[2] - t[0]))) for t in tris
    ])


def _build(cell, simplify=None):
    result = compile_and_run(PLATE, "PLATE", {}, representation="sdf",
                             sdf_cell_mm=cell, sdf_simplify=simplify)
    assert result.success, result.error_message
    return result.geometry


@pytest.mark.sdf_simplify
@pytest.mark.parametrize("cell", (1.0, 0.5))
def test_no_emitted_face_is_degenerate_to_the_library(cell):
    # Nothing at or under the tolerance where yapCAD's predicates give up.
    assert _cross_magnitudes(_build(cell)).min() > epsilon


@pytest.mark.sdf_simplify
@pytest.mark.parametrize("cell", (1.0, 0.5))
def test_every_emitted_face_has_a_normal_at_the_default_tolerance(cell):
    # The library's own predicate must agree that these faces are real.
    missing = [t for t in _faces(_build(cell))
               if triangle_normal(t[0], t[1], t[2]) is None]
    assert missing == []


@pytest.mark.sdf_simplify
def test_the_floor_does_not_cost_the_simplification():
    simplified = len(_faces(_build(0.5)))
    uniform = len(_faces(_build(0.5, simplify=False)))
    assert simplified < 0.1 * uniform
