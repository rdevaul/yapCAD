"""FaceNamer against real solids (issue #107).

``FaceNamer`` imported four ``geom3d`` helpers that did not exist, so it had
never constructed; the selector tests in ``test_analysis_gmsh.py`` only ever
fed it hand-built ``FaceInfo`` records.  These build it from solids, and in
particular from solids whose surfaces do not line up with their faces.
"""

import math

import pytest

from yapcad import sdf
from yapcad.geom import point
from yapcad.geom3d import solid, surface
from yapcad.geom3d_util import conic, prism
from yapcad.package.analysis import FaceNamer
from yapcad.package.analysis.face_naming import (
    ByAreaSelector,
    CombinedSelector,
    bottom_faces,
    faces_at_z_max,
    largest_face,
    right_faces,
    top_faces,
)
from yapcad.sdf.deferred import DeferredMeshError


def one_surface(sld):
    """The same solid with all its surfaces merged into one, as a mesh
    boolean or an SDF part would deliver it."""
    verts, norms, faces = [], [], []
    for s in sld[1]:
        base = len(verts)
        verts += s[1]
        norms += s[2]
        faces += [[base + i for i in f] for f in s[3]]
    return solid([surface(verts, norms, faces)])


BOX = prism(4, 2, 6, center=point(1, 0, 3))   # x -1..3, y -1..1, z 0..6


@pytest.mark.parametrize("box", [BOX, one_surface(BOX)], ids=["six", "one"])
class TestBox:
    def test_six_faces_however_the_surfaces_are_divided(self, box):
        faces = FaceNamer(box).faces
        assert len(faces) == 6
        assert sum(f.area for f in faces) == pytest.approx(2 * (8 + 24 + 12))

    def test_each_face_is_measured(self, box):
        for f in FaceNamer(box).faces:
            axis = max(range(3), key=lambda k: abs(f.normal[k]))
            assert abs(f.normal[axis]) == pytest.approx(1.0)
            half = (2.0, 1.0, 3.0)[axis]
            centre = (1.0, 0.0, 3.0)
            assert f.centroid[axis] - centre[axis] == pytest.approx(
                math.copysign(half, f.normal[axis]))
            # Area of the face perpendicular to ``axis``.
            dims = [4.0, 2.0, 6.0]
            dims.pop(axis)
            assert f.area == pytest.approx(dims[0] * dims[1])

    def test_selectors_pick_the_right_faces(self, box):
        namer = FaceNamer(box)
        result = namer.name_faces({
            "top": top_faces(),
            "bottom": bottom_faces(),
            "big_right": CombinedSelector([right_faces(),
                                           ByAreaSelector(min_area=10.0)]),
        })
        faces = namer.faces
        assert [faces[i].normal[2] for i in result["top"]] == [pytest.approx(1.0)]
        assert [faces[i].normal[2] for i in result["bottom"]] == [pytest.approx(-1.0)]
        assert len(result["big_right"]) == 1      # the 2x6 face, area 12
        assert namer.get_named_faces() == result

    def test_largest_and_highest(self, box):
        namer = FaceNamer(box)
        result = namer.name_faces({"largest": largest_face(),
                                   "highest": faces_at_z_max()})
        largest = [namer.faces[i] for i in result["largest"]]
        # Two 4x6 faces tie for largest.
        assert len(largest) == 2 and all(f.area == pytest.approx(24.0) for f in largest)
        assert [namer.faces[i].centroid[2] for i in result["highest"]] == [pytest.approx(6.0)]

    def test_triangles_partition_the_mesh(self, box):
        faces = FaceNamer(box).faces
        located = [t for f in faces for t in f.triangles]
        assert len(located) == len(set(located)) == sum(len(s[3]) for s in box[1])


def test_sdf_part_is_one_surface_and_six_faces():
    part = sdf.to_solid(sdf.box((20.0, 16.0, 10.0)), resolution=16)
    assert len(part[1]) == 1
    areas = sorted(f.area for f in FaceNamer(part).faces)
    assert areas == pytest.approx([160, 160, 200, 200, 320, 320])


class TestCurved:
    def test_cylinder_is_two_caps_and_a_wall(self):
        cyl = conic(2.0, 2.0, 5.0, max_chord=0.05)
        faces = FaceNamer(cyl).faces
        assert len(faces) == 3
        wall = max(faces, key=lambda f: f.area)
        assert wall.area == pytest.approx(2 * math.pi * 2.0 * 5.0, rel=1e-3)
        # A whole wall faces every way at once, so it has no mean normal and
        # no direction selector picks it up.
        assert wall.normal == (0.0, 0.0, 0.0)
        assert not top_faces().matches(wall)

    def test_crease_angle_controls_splitting(self):
        # An octagonal prism: 45 degrees between neighbouring sides.
        octagon = conic(2.0, 2.0, 5.0, angr=45.0)
        assert len(FaceNamer(octagon).faces) == 10
        assert len(FaceNamer(octagon, crease_deg=50.0).faces) == 3


class TestRefusals:
    def test_not_a_solid(self):
        with pytest.raises(ValueError):
            FaceNamer(["not", "a", "solid"])

    def test_deferred_solid(self):
        with sdf.deferring():
            part = sdf.to_solid(sdf.box((2.0, 2.0, 2.0)))
        with pytest.raises(DeferredMeshError):
            FaceNamer(part)

    def test_metadata_round_trip(self):
        namer = FaceNamer(BOX)
        namer.name_faces({"top": top_faces()})
        meta = namer.to_metadata()
        assert meta["face_names"] == namer.get_named_faces()
        assert len(meta["faces"]) == 6
