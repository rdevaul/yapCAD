"""Regression tests for centred scaling and the Geometry solid transforms.

Every representation of a solid must follow a transform the same way.
These pin five defects that let them disagree:

* ``geom.scale`` and ``geom3d.scalesurface`` composed a centred scale as
  ``T(-c) S T(c)`` -- ``S(p + c) - c``, a scale about ``-c`` -- while both
  BREP hooks scale about ``c``.  A centred scale of a BREP solid left its
  mesh and its BREP in different places.
* ``geom3d.scalesolid`` left the BREP attached, unscaled, after a
  non-uniform scale.  An attached BREP is read as authoritative, so a saved
  document and an analytic STEP export silently described the unscaled part.
* ``Geometry.mirror`` on a mesh-only solid did nothing.
* ``Geometry.scale`` on a mesh-only solid produced an invalid solid, its
  surfaces missing their ``'surface'`` tag.
* ``Geometry.scale`` and ``Geometry.mirror`` on a BREP solid rebuilt the
  mesh from the BREP and dropped the solid's metadata.
"""

import math

import pytest

from yapcad.brep import (
    _clear_brep_data,
    brep_from_solid,
    has_brep_data,
    occ_available,
)
from yapcad.geom import arc, line, point, scale
from yapcad.geom3d import (
    issolid,
    issolidclosed,
    mirrorsolid,
    scalesolid,
    scalesurface,
    solidbbox,
    translatesolid,
    volumeof,
)
from yapcad.geom3d_util import prism
from yapcad.geometry import Geometry
from yapcad.metadata import get_solid_metadata

needs_occ = pytest.mark.skipif(not occ_available(),
                               reason="pythonocc-core is not available")

CENTRE = point(20.0, 0.0, 0.0)


def centred(p, s, c=CENTRE):
    """Where a centred scale must send ``p``: ``c + S (p - c)``."""
    sx, sy, sz = s if isinstance(s, tuple) else (s, s, s)
    return (c[0] + sx * (p[0] - c[0]),
            c[1] + sy * (p[1] - c[1]),
            c[2] + sz * (p[2] - c[2]))


def mesh_only_box(offset=0.0):
    """A 10-unit cube with its BREP stripped, so only the mesh remains."""
    solid = prism(10.0, 10.0, 10.0)
    _clear_brep_data(solid)
    if offset:
        solid = translatesolid(solid, point(offset, 0.0, 0.0))
    return solid


def x_range(solid):
    lo, hi = solidbbox(solid)[:2]
    return lo[0], hi[0]


def signed_volume(solid):
    """Positive for outward winding; volumeof takes the absolute value."""
    total = 0.0
    for surf in solid[1]:
        verts = surf[1]
        for a, b, c in surf[3]:
            p, q, r = verts[a], verts[b], verts[c]
            total += (p[0] * (q[1] * r[2] - q[2] * r[1])
                      - p[1] * (q[0] * r[2] - q[2] * r[0])
                      + p[2] * (q[0] * r[1] - q[1] * r[0])) / 6.0
    return total


def brep_x_range(solid):
    from OCC.Core.Bnd import Bnd_Box
    from OCC.Core.BRepBndLib import brepbndlib
    box = Bnd_Box()
    brepbndlib.AddOptimal(brep_from_solid(solid).shape, box, False, False)
    x0, _, _, x1, _, _ = box.Get()
    return x0, x1


# ---------------------------------------------------------------------------
# Centred scale composition
# ---------------------------------------------------------------------------


class TestCentredScale:

    def test_a_point_scales_about_the_centre(self):
        assert scale(point(5.0, 0.0, 0.0), 2.0, cent=CENTRE)[:3] == \
            pytest.approx(centred(point(5.0, 0.0, 0.0), 2.0))

    def test_the_centre_is_a_fixed_point(self):
        assert scale(CENTRE, 3.0, cent=CENTRE)[:3] == \
            pytest.approx(CENTRE[:3])

    def test_scaling_about_the_origin_is_unchanged(self):
        assert scale(point(5.0, 1.0, 2.0), 2.0)[:3] == \
            pytest.approx((10.0, 2.0, 4.0))

    def test_a_line_scales_about_the_centre(self):
        result = scale(line(point(0, 0, 0), point(10, 0, 0)), 2.0,
                       cent=CENTRE)
        assert result[0][:3] == pytest.approx(centred(point(0, 0, 0), 2.0))
        assert result[1][:3] == pytest.approx(centred(point(10, 0, 0), 2.0))

    def test_an_arc_keeps_its_centre_relationship(self):
        result = scale(arc(point(25.0, 0.0), 3.0), 2.0, cent=CENTRE)
        assert result[0][:2] == pytest.approx((30.0, 0.0))
        assert result[1][0] == pytest.approx(6.0)

    def test_a_surface_scales_about_the_centre_non_uniformly(self):
        surf = mesh_only_box()[1][0]
        result = scalesurface(surf, 2.0, 3.0, 0.5, CENTRE)
        for before, after in zip(surf[1], result[1]):
            assert after[:3] == pytest.approx(
                centred(before, (2.0, 3.0, 0.5)))

    def test_a_solid_scales_about_the_centre(self):
        assert x_range(scalesolid(mesh_only_box(), 2.0, cent=CENTRE)) == \
            pytest.approx((-30.0, -10.0))

    @needs_occ
    def test_mesh_and_brep_agree_after_a_centred_scale(self):
        """The reason this matters: both BREP hooks already scaled about
        the centre, so the mesh was the representation that disagreed."""
        solid = scalesolid(prism(10.0, 10.0, 10.0), 2.0, cent=CENTRE)
        assert x_range(solid) == pytest.approx((-30.0, -10.0))
        assert brep_x_range(solid) == pytest.approx((-30.0, -10.0),
                                                    abs=1e-6)


# ---------------------------------------------------------------------------
# Non-uniform scale and the BREP
# ---------------------------------------------------------------------------


@needs_occ
class TestNonUniformScaleDropsTheBrep:

    def test_a_non_uniform_scale_removes_the_stale_brep(self):
        solid = scalesolid(prism(10.0, 10.0, 10.0), 3.0, 1.0, 1.0)
        assert not has_brep_data(solid)
        assert x_range(solid) == pytest.approx((-15.0, 15.0))

    def test_a_uniform_scale_keeps_and_scales_the_brep(self):
        solid = scalesolid(prism(10.0, 10.0, 10.0), 3.0)
        assert has_brep_data(solid)
        assert brep_x_range(solid) == pytest.approx((-15.0, 15.0), abs=1e-6)

    def test_the_original_keeps_its_brep(self):
        original = prism(10.0, 10.0, 10.0)
        scalesolid(original, 3.0, 1.0, 1.0)
        assert has_brep_data(original)

    def test_a_saved_document_no_longer_claims_the_unscaled_part(self):
        """The consequence that made the stale BREP a real bug: it was
        written out as the authoritative definition of the solid."""
        from yapcad.io.geometry_json import geometry_to_json
        solid = scalesolid(prism(10.0, 10.0, 10.0), 3.0, 1.0, 1.0)
        document = geometry_to_json([solid], units="mm")
        entry = next(e for e in document["entities"]
                     if e["type"] == "solid")
        assert entry["representations"]["authoritative"] == "mesh"
        assert "brep" not in entry["representations"]


# ---------------------------------------------------------------------------
# The Geometry wrapper's solid transforms
# ---------------------------------------------------------------------------


class TestGeometryMirror:

    def test_a_mesh_only_solid_is_mirrored(self):
        g = Geometry(mesh_only_box(offset=20.0))
        g.mirror("yz")
        assert x_range(g.elem) == pytest.approx((-25.0, -15.0))

    def test_the_mirrored_solid_is_closed_and_wound_outward(self):
        g = Geometry(mesh_only_box(offset=20.0))
        g.mirror("yz")
        assert x_range(g.elem)[1] < 0.0   # it really moved
        assert issolidclosed(g.elem)
        assert signed_volume(g.elem) == pytest.approx(1000.0)

    @pytest.mark.parametrize("plane,axis", [("yz", 0), ("xz", 1), ("xy", 2)])
    def test_every_plane_mirrors_its_own_axis(self, plane, axis):
        offset = point(*[20.0 if i == axis else 0.0 for i in range(3)])
        g = Geometry(translatesolid(mesh_only_box(), offset))
        g.mirror(plane)
        lo, hi = solidbbox(g.elem)[:2]
        assert (lo[axis], hi[axis]) == pytest.approx((-25.0, -15.0))

    def test_it_matches_mirrorsolid(self):
        # Offset along y, so that mirroring across xz actually moves it.
        solid = translatesolid(mesh_only_box(), point(0.0, 20.0, 0.0))
        g = Geometry(solid)
        g.mirror("xz")
        assert solidbbox(g.elem) == solidbbox(mirrorsolid(solid, "xz"))
        assert solidbbox(g.elem)[1][1] < 0.0

    @needs_occ
    def test_a_brep_solid_keeps_mesh_brep_and_metadata_together(self):
        solid = translatesolid(prism(10.0, 10.0, 10.0),
                               point(20.0, 0.0, 0.0))
        get_solid_metadata(solid)["name"] = "bracket"
        g = Geometry(solid)
        g.mirror("yz")
        assert x_range(g.elem) == pytest.approx((-25.0, -15.0))
        assert brep_x_range(g.elem) == pytest.approx((-25.0, -15.0),
                                                     abs=1e-6)
        assert get_solid_metadata(g.elem)["name"] == "bracket"


class TestGeometryScale:

    def test_a_mesh_only_solid_stays_a_valid_solid(self):
        g = Geometry(mesh_only_box())
        g.scale(2.0)
        assert issolid(g.elem, fast=False)
        assert issolidclosed(g.elem)
        assert volumeof(g.elem) == pytest.approx(8000.0)

    def test_it_scales_about_the_given_centre(self):
        g = Geometry(mesh_only_box())
        g.scale(2.0, cent=CENTRE)
        assert x_range(g.elem) == pytest.approx((-30.0, -10.0))

    def test_the_material_and_construction_slots_survive(self):
        solid = mesh_only_box()
        solid[2] = ["steel"]
        solid[3] = ["procedure", "prism(10,10,10)"]
        g = Geometry(solid)
        g.scale(2.0)
        assert g.elem[2] == ["steel"]
        assert g.elem[3] == ["procedure", "prism(10,10,10)"]

    def test_anisotropic_scaling_is_still_refused(self):
        """The Geometry wrapper's documented contract is unchanged."""
        with pytest.raises(NotImplementedError):
            Geometry(mesh_only_box()).scale(2.0, sy=1.0)

    @needs_occ
    def test_a_brep_solid_keeps_mesh_brep_and_metadata_together(self):
        solid = prism(10.0, 10.0, 10.0)
        get_solid_metadata(solid)["name"] = "bracket"
        g = Geometry(solid)
        g.scale(2.0, cent=CENTRE)
        assert x_range(g.elem) == pytest.approx((-30.0, -10.0))
        assert brep_x_range(g.elem) == pytest.approx((-30.0, -10.0),
                                                     abs=1e-6)
        assert get_solid_metadata(g.elem)["name"] == "bracket"
        assert math.isclose(volumeof(g.elem), 8000.0, rel_tol=1e-9)
