import math
import os
import shutil

import pytest

from yapcad.combine import Boolean
from yapcad.geom import arc, point
from yapcad.geom3d import poly2surfaceXY, reversesurface
from yapcad.geom3d_util import extrude
from yapcad.geom_util import geomlist2poly, geomlist2poly_components
from yapcad.io import write_stl
from yapcad.poly import Polygon, Rect


@pytest.mark.slow
@pytest.mark.parametrize('corner_radius', [0.8])
def test_extrude_pentagon_plate(tmp_path, corner_radius):
    outer_radius = 12.0
    outer = Polygon()
    for i in range(5):
        ang = (2 * math.pi * i / 5.0) + math.radians(18.0)
        corner = point(outer_radius * math.cos(ang),
                       outer_radius * math.sin(ang))
        outer.addArc(arc(corner, corner_radius))

    hole_centers = []
    inner_ring_radius = 4.5
    hole_radius = 1.2
    for i in range(5):
        ang = 2 * math.pi * i / 5.0
        hole_centers.append(point(inner_ring_radius * math.cos(ang),
                                  inner_ring_radius * math.sin(ang)))
    hole_centers.append(point(0, 0))

    holes = []
    for center in hole_centers:
        hole = Polygon()
        hole.addArc(arc(center, hole_radius))
        holes.append(hole)

    shape = outer
    for hole in holes:
        shape = Boolean('difference', [shape, hole])

    components = geomlist2poly_components(shape.geom)
    assert components, 'expected at least one polygon component'
    outer_loop, hole_loops = components[0]
    assert len(hole_loops) == len(holes)

    surface, _ = poly2surfaceXY(outer_loop, holepolys=hole_loops)
    solid_plate = extrude(surface, distance=3.0)

    assert solid_plate[0] == 'solid'
    assert len(solid_plate[1]) == 3  # caps + side strip

    strip = solid_plate[1][1]
    base_vertices = len(reversesurface(surface)[1])
    strip_vertices = strip[1]
    strip_normals = strip[2]

    loops = []
    reversed_surface = reversesurface(surface)
    if reversed_surface[4]:
        loops.append(reversed_surface[4])
    loops.extend(loop for loop in reversed_surface[5] if loop)

    assert len(loops) == 1 + len(holes)

    outer_loop_indices = loops[0]
    for vidx in outer_loop_indices:
        if vidx >= base_vertices:
            continue
        vertex = strip_vertices[vidx]
        normal = strip_normals[vidx]
        outward = (vertex[0], vertex[1])
        assert outward[0] * normal[0] + outward[1] * normal[1] > 0

    for hole_idx, loop_indices in enumerate(loops[1:]):
        center = hole_centers[hole_idx]
        for vidx in loop_indices:
            if vidx >= base_vertices:
                continue
            vertex = strip_vertices[vidx]
            normal = strip_normals[vidx]
            toward_center = (center[0] - vertex[0], center[1] - vertex[1])
            assert toward_center[0] * normal[0] + toward_center[1] * normal[1] > 0

    stl_path = tmp_path / 'extrusion_test.stl'
    write_stl(solid_plate, stl_path, binary=True, name='extrusion_test')
    assert stl_path.exists()

    export_path = os.environ.get('YAPCAD_STL_OUTPUT')
    if export_path:
        target_path = export_path
        if os.path.isdir(export_path):
            target_path = os.path.join(export_path, 'extrusion_test.stl')
        else:
            os.makedirs(os.path.dirname(target_path) or '.', exist_ok=True)
        shutil.copyfile(stl_path, target_path)
