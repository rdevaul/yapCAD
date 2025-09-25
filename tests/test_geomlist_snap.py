from yapcad.geom import line, point, dist, epsilon
from yapcad.geom_util import geomlist2poly

def test_geomlist2poly_snaps_small_gap():
    line_one = line(point(0, 0), point(1, 0))
    line_two = line(point(1.001, 0), point(1.001, 1))

    poly = geomlist2poly([line_one, line_two])

    assert dist(poly[1], point(1, 0)) <= epsilon
