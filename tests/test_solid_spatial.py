import copy
import pytest

from yapcad.geom import point
from yapcad.geom3d import (
    reversesurface,
    solid_boolean,
    solid_contains_point,
    solids_intersect,
    translatesolid,
    _iter_triangles_from_solid,
    tri2p0n,
)
from yapcad.geom3d_util import prism, sphere, conic, tube


def _offset_point(base, normal, scale):
    return point(base[0] + normal[0] * scale,
                 base[1] + normal[1] * scale,
                 base[2] + normal[2] * scale)


def _assert_normals_outward(sld, *, label='solid', step=1e-3):
    triangles = list(_iter_triangles_from_solid(sld))
    assert triangles, f"{label} produced no faces"

    max_checks = min(32, len(triangles))
    stride = max(1, len(triangles) // max_checks)
    checked = 0

    for idx, tri in enumerate(triangles):
        if idx % stride != 0:
            continue
        try:
            center, normal = tri2p0n(tri)
        except ValueError:
            continue

        scale = step
        detected = False
        flipped = False
        for _ in range(6):
            inside_pt = _offset_point(center, normal, -scale)
            outside_pt = _offset_point(center, normal, scale)
            inside = solid_contains_point(sld, inside_pt)
            outside = solid_contains_point(sld, outside_pt)
            if inside and not outside:
                detected = True
                break
            if outside and not inside:
                flipped = True
                break
            scale *= 2.0
        if flipped:
            raise AssertionError(
                f"normal orientation check failed for {label} at center {center[:3]}"
            )
        if detected:
            checked += 1

    assert checked > 0, f"no triangles sampled for {label}"


def _translate(sld, delta):
    return translatesolid(sld, point(delta[0], delta[1], delta[2]))


def _make_cavity_solid():
    outer = prism(2, 2, 2)
    inner = prism(1, 1, 1)
    cavity_surfaces = outer[1] + [reversesurface(s) for s in inner[1]]
    return ['solid', cavity_surfaces, [], []]


def test_solid_contains_point_basic():
    cube = prism(2, 2, 2)
    assert solid_contains_point(cube, point(0.1, 0.1, 0.1))
    assert solid_contains_point(cube, point(1.0, 0.0, 0.0))
    assert not solid_contains_point(cube, point(2.5, 0.0, 0.0))


def test_solid_contains_point_cavity():
    cavity = _make_cavity_solid()
    assert not solid_contains_point(cavity, point(0.1, 0.0, 0.0))
    assert solid_contains_point(cavity, point(0.9, 0.0, 0.0))


@pytest.mark.slow
def test_solid_boolean_cubes():
    a = prism(2, 2, 2)
    b = _translate(prism(2, 2, 2), (0.75, 0.0, 0.0))

    union = solid_boolean(a, b, 'union')
    assert solid_contains_point(union, point(-0.9, 0.0, 0.0))
    assert solid_contains_point(union, point(1.6, 0.0, 0.0))

    intersection = solid_boolean(a, b, 'intersection')
    assert solid_contains_point(intersection, point(0.6, 0.0, 0.0))
    assert not solid_contains_point(intersection, point(-0.9, 0.0, 0.0))

    difference = solid_boolean(a, b, 'difference')
    assert solid_contains_point(difference, point(-0.9, 0.0, 0.0))
    assert not solid_contains_point(difference, point(0.6, 0.0, 0.0))


@pytest.mark.slow
def test_solid_boolean_spheres():
    a = sphere(2.0)
    b = copy.deepcopy(a)
    b = translatesolid(b, point(1.0, 0.0, 0.0))

    union = solid_boolean(a, b, 'union')
    assert solid_contains_point(union, point(-0.8, 0.0, 0.0))
    assert solid_contains_point(union, point(1.8, 0.0, 0.0))

    intersection = solid_boolean(a, b, 'intersection')
    assert solid_contains_point(intersection, point(0.5, 0.0, 0.0))
    assert not solid_contains_point(intersection, point(-1.8, 0.0, 0.0))

    difference = solid_boolean(a, b, 'difference')
    assert solid_contains_point(difference, point(-0.8, 0.0, 0.0))
    assert not solid_contains_point(difference, point(0.5, 0.0, 0.0))


@pytest.mark.slow
def test_solids_intersect_detection():
    base = prism(2, 2, 2)
    shifted = _translate(prism(2, 2, 2), (0.5, 0.0, 0.0))
    assert solids_intersect(base, shifted)

    far = _translate(prism(2, 2, 2), (10.0, 0.0, 0.0))
    assert not solids_intersect(base, far)


@pytest.mark.slow
def test_solid_boolean_sphere_normals():
    a = sphere(2.0)
    b = copy.deepcopy(a)
    b = translatesolid(b, point(1.0, 0.0, 0.0))

    union = solid_boolean(a, b, 'union')
    _assert_normals_outward(union, label='sphere union')

    difference = solid_boolean(a, b, 'difference')
    _assert_normals_outward(difference, label='sphere difference')


@pytest.mark.slow
def test_solid_boolean_box_minus_cylinder():
    box = prism(2, 2, 2)
    drill = conic(0.6, 0.6, 3.0, center=point(0.0, 0.0, -1.5))

    result = solid_boolean(box, drill, 'difference')

    assert not solid_contains_point(result, point(0.0, 0.0, 0.0))
    assert solid_contains_point(result, point(0.9, 0.0, 0.0))
    _assert_normals_outward(result, label='box minus cylinder')


@pytest.mark.slow
def test_solid_boolean_tube_side_hole():
    shell = tube(outer_diameter=3.0, wall_thickness=0.5, length=4.0,
                 base_point=point(0.0, 0.0, -2.0))
    drill = conic(0.6, 0.6, 4.2, center=point(1.2, 0.0, -2.1))

    result = solid_boolean(shell, drill, 'difference')

    assert not solid_contains_point(result, point(1.2, 0.0, 0.0))
    assert solid_contains_point(result, point(0.0, 1.3, 0.0))
    _assert_normals_outward(result, label='tube side hole')
