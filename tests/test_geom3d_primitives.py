import math

from yapcad.combine import Boolean
from yapcad.geom import arc, point
from yapcad.geom3d import solidbbox
from yapcad.geom3d_util import tube, conic_tube, spherical_shell, stack_solids
from yapcad.geom_util import geomlist2poly_components
from yapcad.poly import Polygon


def _radius(point_xyz, center_xy):
    dx = point_xyz[0] - center_xy[0]
    dy = point_xyz[1] - center_xy[1]
    return math.hypot(dx, dy)


def _surface_vertices(surface):
    return surface[1]


def _surface_normals(surface):
    return surface[2]


def test_boolean_difference_handles_arcs():
    outer = Polygon()
    outer.addArc(arc(point(0, 0), 10.0))
    hole = Polygon()
    hole.addArc(arc(point(0, 0), 2.0))

    shape = Boolean('difference', [outer, hole])
    outer_poly, holes = geomlist2poly_components(shape.geom)[0]
    assert len(outer_poly) > 3
    assert len(holes) == 1


def test_tube_basic_properties():
    sld = tube(outer_diameter=10.0, wall_thickness=2.0, length=8.0)
    assert sld[0] == 'solid'
    assert len(sld[1]) == 4

    bbox = solidbbox(sld)
    outer_radius = 5.0
    assert math.isclose(bbox[0][2], 0.0, abs_tol=1e-6)
    assert math.isclose(bbox[1][2], 8.0, abs_tol=1e-6)
    assert math.isclose(bbox[1][0], outer_radius, abs_tol=1e-4)

    outer_side = sld[1][1]
    verts = _surface_vertices(outer_side)
    norms = _surface_normals(outer_side)
    sample_idx = 0
    radius_vec = (verts[sample_idx][0], verts[sample_idx][1])
    normal = norms[sample_idx]
    assert radius_vec[0] * normal[0] + radius_vec[1] * normal[1] > 0

    inner_side = sld[1][3]
    verts_inner = _surface_vertices(inner_side)
    norms_inner = _surface_normals(inner_side)
    center = (0.0, 0.0)
    vec_to_center = (center[0] - verts_inner[0][0], center[1] - verts_inner[0][1])
    dot_inner = vec_to_center[0] * norms_inner[0][0] + vec_to_center[1] * norms_inner[0][1]
    assert dot_inner > 0


def test_conic_tube_tapers():
    sld = conic_tube(bottom_outer_diameter=12.0,
                     top_outer_diameter=8.0,
                     wall_thickness=1.5,
                     length=6.0)
    bbox = solidbbox(sld)
    assert math.isclose(bbox[0][2], 0.0, abs_tol=1e-6)
    assert math.isclose(bbox[1][2], 6.0, abs_tol=1e-6)

    top_surface = sld[1][2]
    center_xy = (0.0, 0.0)
    radii = [_radius(v[:3], center_xy) for v in _surface_vertices(top_surface)]
    assert math.isclose(max(radii), 4.0, rel_tol=1e-2)


def test_spherical_shell_full():
    sld = spherical_shell(outer_diameter=10.0, wall_thickness=1.0)
    assert len(sld[1]) == 2
    outer_surface = sld[1][0]
    inner_surface = sld[1][1]
    center = (0.0, 0.0, 0.0)
    sample_outer = _surface_vertices(outer_surface)[0]
    normal_outer = _surface_normals(outer_surface)[0]
    vec_outer = (sample_outer[0] - center[0], sample_outer[1] - center[1], sample_outer[2] - center[2])
    assert vec_outer[0] * normal_outer[0] + vec_outer[1] * normal_outer[1] + vec_outer[2] * normal_outer[2] > 0

    sample_inner = _surface_vertices(inner_surface)[0]
    normal_inner = _surface_normals(inner_surface)[0]
    vec_inner = (center[0] - sample_inner[0], center[1] - sample_inner[1], center[2] - sample_inner[2])
    assert vec_inner[0] * normal_inner[0] + vec_inner[1] * normal_inner[1] + vec_inner[2] * normal_inner[2] > 0


def test_spherical_shell_cap_has_conic_side():
    sld = spherical_shell(outer_diameter=10.0,
                          wall_thickness=1.0,
                          solid_angle=2 * math.pi)
    assert len(sld[1]) == 3
    conic_surface = sld[1][2]
    radii = {_radius(v[:3], (0.0, 0.0)) for v in _surface_vertices(conic_surface)}
    assert min(radii) < max(radii)


def test_stack_solids_spacing_directives():
    a = tube(outer_diameter=4.0, wall_thickness=0.5, length=2.0)
    b = tube(outer_diameter=4.0, wall_thickness=0.5, length=2.0)
    c, d = stack_solids(['space:3.0', a, 'space:1.5', b], axis='z', start=0.0, gap=0.0)

    bbox_c = solidbbox(c)
    bbox_d = solidbbox(d)

    assert math.isclose(bbox_c[0][2], 3.0, abs_tol=1e-6)
    assert math.isclose(bbox_d[0][2] - bbox_c[1][2], 1.5, abs_tol=1e-6)
