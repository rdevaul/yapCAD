from yapcad.poly import Polyline, Polygon
from yapcad.geom import point, vclose


def test_polyline_center_and_sample():
    pts = [point(0,0), point(2,0), point(2,2)]
    pl = Polyline(pts)
    assert vclose(pl.getCenter(), point(4/3, 2/3))
    assert vclose(pl.sample(0.75), point(2, 1))


def test_polygon_wrap_sampling():
    square = Polygon([point(0,0), point(1,0), point(1,1), point(0,1)])
    assert vclose(square.getCenter(), point(0.5, 0.5))
    assert vclose(square.sample(1.25), point(1, 0))
