"""Ground-truth tests for the native mesh boolean engine.

The engine's existing tests check point containment and outward normals,
never whether a result is closed or has the right volume, so open and
double-walled results passed.  These check both, against closed-form
volumes, with ``engine="native"`` forced -- otherwise, with OCC installed,
``prism()`` carries a BREP and the booleans quietly go through OCC.

Cases the engine gets right are regression tests.  Cases it gets wrong are
``xfail(strict=True)`` with the defect named: they document the current
state, and a fix that makes one pass will fail the run until its marker is
removed.
"""

import pytest

from yapcad.brep import _clear_brep_data
from yapcad.geom import point
from yapcad.geom3d import issolidclosed, solid_boolean, translatesolid
from yapcad.geom3d_util import conic, prism, sphere
from yapcad.mesh import mesh_view


def mesh_only(solid):
    _clear_brep_data(solid)
    return solid


def moved(solid, delta):
    return mesh_only(translatesolid(solid, point(*delta)))


def signed_volume(solid):
    total = 0.0
    for _, a, b, c in mesh_view(solid):
        total += (a[0] * (b[1] * c[2] - b[2] * c[1])
                  - a[1] * (b[0] * c[2] - b[2] * c[0])
                  + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6.0
    return total


def native(a, b, op):
    return solid_boolean(a, b, op, engine="native")


def assert_correct(result, volume):
    assert issolidclosed(result), "result is not a closed solid"
    assert signed_volume(result) == pytest.approx(volume, rel=1e-6,
                                                  abs=1e-9)


# ---------------------------------------------------------------------------
# What the engine gets right
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("op,volume", [
    ("union", 64.0), ("intersection", 8.0), ("difference", 56.0),
])
def test_nested_boxes(op, volume):
    """No triangle crosses the other surface, so every one is classified by
    a single containment test -- the path the octree fallback used to
    replace with a clip against every target plane."""
    assert_correct(native(mesh_only(prism(4, 4, 4)),
                          mesh_only(prism(2, 2, 2)), op), volume)


@pytest.mark.parametrize("op,volume", [
    ("union", 8.0), ("difference", 0.0),
])
def test_a_sphere_inside_a_box(op, volume):
    # geom3d_util.sphere takes a diameter: this one sits wholly inside.
    assert_correct(native(mesh_only(prism(2, 2, 2)),
                          mesh_only(sphere(1.2)), op)
                   if op == "union" else
                   native(mesh_only(sphere(1.2)),
                          mesh_only(prism(2, 2, 2)), op), volume)


def test_disjoint_union_keeps_both():
    assert_correct(native(mesh_only(prism(2, 2, 2)),
                          moved(prism(2, 2, 2), (5.0, 0.0, 0.0)), "union"),
                   16.0)


@pytest.mark.parametrize("op,volume", [
    ("intersection", 0.0), ("difference", 8.0),
])
def test_face_sharing_boxes_except_union(op, volume):
    assert_correct(native(mesh_only(prism(2, 2, 2)),
                          moved(prism(2, 2, 2), (2.0, 0.0, 0.0)), op),
                   volume)


# ---------------------------------------------------------------------------
# Known defects
# ---------------------------------------------------------------------------

CRACKS = ("each triangle is cut independently, so neighbours' cut points "
          "do not coincide and the result has gaps and T-junctions")
COPLANAR = ("coincident faces are not classified as same/opposite, so they "
            "are kept as internal double walls or dropped")


@pytest.mark.xfail(strict=True, reason=CRACKS)
@pytest.mark.parametrize("op,volume", [
    ("union", 12.175), ("intersection", 3.825), ("difference", 4.175),
])
def test_overlapping_boxes_in_general_position(op, volume):
    """The volume comes out exact; the mesh is not closed."""
    assert_correct(native(mesh_only(prism(2, 2, 2)),
                          moved(prism(2, 2, 2), (0.75, 0.3, 0.2)), op),
                   volume)


@pytest.mark.xfail(strict=True, reason=COPLANAR)
def test_face_sharing_union():
    assert_correct(native(mesh_only(prism(2, 2, 2)),
                          moved(prism(2, 2, 2), (2.0, 0.0, 0.0)), "union"),
                   16.0)


@pytest.mark.xfail(strict=True, reason=COPLANAR)
@pytest.mark.parametrize("op,volume", [
    ("union", 12.0), ("intersection", 4.0), ("difference", 4.0),
])
def test_boxes_with_coplanar_side_faces(op, volume):
    """Four side faces coincide over the overlap; volumes come out 67%
    wrong."""
    assert_correct(native(mesh_only(prism(2, 2, 2)),
                          moved(prism(2, 2, 2), (1.0, 0.0, 0.0)), op),
                   volume)


@pytest.mark.xfail(strict=True, reason=CRACKS)
@pytest.mark.parametrize("op", ["union", "intersection", "difference"])
def test_box_and_through_cylinder(op):
    result = native(mesh_only(prism(2, 2, 2)),
                    mesh_only(conic(0.6, 0.6, 3.0,
                                    center=point(0, 0, -1.5))), op)
    assert issolidclosed(result)


@pytest.mark.xfail(strict=True, reason=CRACKS)
@pytest.mark.parametrize("op", ["union", "intersection", "difference"])
def test_overlapping_spheres(op):
    result = native(mesh_only(sphere(1.0)),
                    moved(sphere(1.0), (0.8, 0.1, 0.0)), op)
    assert issolidclosed(result)
