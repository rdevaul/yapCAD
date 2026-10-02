"""Gears as fields.

A ``straight_bevel_gear`` node stores the gear's specification -- a dozen
numbers -- rather than its tooth outline, so the tree stays a few hundred
bytes and re-evaluates from the same spec :mod:`yapcad.gears.bevel` builds
its BREP from.  Its field is the tooth section at the outer plane, as an
exact polygon, carried to the apex (:func:`~yapcad.sdf.planar.apex_extrude`)
between the inner and outer planes, minus the bore.  That is the same
construction as the BREP generator's ruled loft, whose two sections are the
same outline scaled toward the apex.

It replays through OCC by calling that generator, so it is ``csg_exact``.
The two differ only in how a flank runs between its sample points: the
field joins them with straight segments, the BREP with an interpolating
B-spline.  With the default 17 samples per flank the field lies within a
few micrometres of the BREP; with the BREP generator's own default of 5 the
gap is larger, so the field samples the involute more densely.
"""

import math
from functools import lru_cache

from yapcad.sdf.node import FieldProps, NodeSpec, SdfError, analyze, \
    make_node, register, register_backend
from yapcad.sdf.ops import subtract, translate
from yapcad.sdf.planar import apex_extrude, polygon
from yapcad.sdf.primitives import cylinder

__all__ = ["straight_bevel_gear", "DEFAULT_FLANK_SAMPLES"]

#: Samples per tooth flank in the field's polygon.
DEFAULT_FLANK_SAMPLES = 17

_SPEC_FIELDS = ("teeth", "mate_teeth", "outer_module_mm", "face_width_mm",
                "shaft_angle_deg", "pressure_angle_deg", "backlash_mm",
                "bore_diameter_mm")


def _spec(params):
    from yapcad.gears.bevel import StraightBevelGearSpec
    return StraightBevelGearSpec(**{k: params[k] for k in _SPEC_FIELDS})


@lru_cache(maxsize=32)
def _tree_for(items):
    """The field tree the node stands for, built once per specification."""
    from yapcad.gears.bevel import derive_straight_bevel_geometry, \
        tooth_section_points
    params = dict(items)
    spec = _spec(params)
    geometry = derive_straight_bevel_geometry(spec)
    outer = tooth_section_points(spec, geometry, cone_distance="outer",
                                 flank_samples=params["flank_samples"])[:-1]
    profile = polygon([(x, y) for x, y, _ in outer], symmetry=spec.teeth)
    body = apex_extrude(profile, geometry.outer_plane_z_mm,
                        geometry.inner_plane_z_mm, geometry.outer_plane_z_mm)
    if spec.bore_diameter_mm > 0.0:
        # As the BREP generator cuts it: from below the gear to above it.
        height = geometry.outer_plane_z_mm + 2.0
        bore = translate(cylinder(spec.bore_diameter_mm / 2.0, height),
                         (0.0, 0.0, height / 2.0 - 1.0))
        body = subtract(body, bore)
    return body


def _tree(params):
    return _tree_for(tuple(sorted(params.items())))


def _analyze(params, _children):
    props = analyze(_tree(params))
    return FieldProps(exact=props.exact, lipschitz=props.lipschitz,
                      bounds=props.bounds, csg_exact=True)


def _eval(params, p, _children, ev):
    return ev(_tree(params), p)


def _occ(params, _children, _build, _region):
    from yapcad.brep import brep_from_solid
    from yapcad.gears.bevel import make_straight_bevel_gear
    solid = make_straight_bevel_gear(_spec(params),
                                     flank_samples=params["flank_samples"])
    return brep_from_solid(solid).shape


register(NodeSpec(
    kind="straight_bevel_gear",
    min_children=0,
    max_children=0,
    analyze=_analyze,
    backends={"numpy": _eval},
))
register_backend("straight_bevel_gear", "occ", _occ)


def straight_bevel_gear(spec, flank_samples=DEFAULT_FLANK_SAMPLES):
    """A straight bevel gear as a field, from a
    :class:`yapcad.gears.bevel.StraightBevelGearSpec`.

    The frame is the BREP generator's: pitch apex at the origin, axis +z,
    a tooth centred on +x.
    """
    from yapcad.gears.bevel import derive_straight_bevel_geometry
    # Validated by the BREP generator's own rules, raising its own errors:
    # an invalid gear fails the same way with or without OCC.
    derive_straight_bevel_geometry(spec)
    samples = int(flank_samples)
    if samples < 3:
        raise SdfError("straight_bevel_gear: flank_samples must be >= 3")
    params = {k: getattr(spec, k) for k in _SPEC_FIELDS}
    params = {k: (int(v) if k in ("teeth", "mate_teeth") else float(v))
              for k, v in params.items()}
    params["flank_samples"] = samples
    return make_node("straight_bevel_gear", params)


def tip_land(spec, flank_samples=DEFAULT_FLANK_SAMPLES):
    """Width of a tooth's top land at the inner (small) end: the thinnest
    feature of the gear, and so what sets the cell size it can be meshed
    at."""
    from yapcad.gears.bevel import derive_straight_bevel_geometry, \
        tooth_section_points
    geometry = derive_straight_bevel_geometry(spec)
    inner = tooth_section_points(spec, geometry, cone_distance="inner",
                                 flank_samples=flank_samples)
    return math.dist(inner[flank_samples - 1][:2], inner[flank_samples][:2])
