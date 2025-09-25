from yapcad.geom import point, epsilon
from yapcad.geom3d import poly2surfaceXY


def test_poly2surfacexy_generates_ccw_triangles():
    poly = [
        point(0, 0),
        point(4, 0),
        point(4, 4),
        point(2, 4),
        point(2, 2),
        point(1, 2),
        point(1, 4),
        point(0, 4),
        point(0, 0),
    ]

    surf, _ = poly2surfaceXY(poly)

    verts = surf[1]
    for face in surf[3]:
        p0, p1, p2 = [verts[i] for i in face]
        area = ((p1[0] - p0[0]) * (p2[1] - p0[1]) -
                (p2[0] - p0[0]) * (p1[1] - p0[1])) / 2.0
        assert area >= -epsilon
