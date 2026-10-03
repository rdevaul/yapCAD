"""Gears as fields.

``spur_gear`` covers spur, helical and herringbone gears: the involute
profile :mod:`yapcad.contrib.figgear` generates -- with the same parameters
the DSL builtins pass it -- as an exact symmetric polygon, extruded from
``z = 0`` to the face width, and twisted for a helix.

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
from yapcad.sdf.planar import apex_extrude, extrude, polygon, twist
from yapcad.sdf.primitives import cylinder

__all__ = ["straight_bevel_gear", "spur_gear", "DEFAULT_FLANK_SAMPLES"]

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


# ---------------------------------------------------------------------------
# spur, helical and herringbone gears
# ---------------------------------------------------------------------------


@lru_cache(maxsize=32)
def _spur_tree(items):
    from yapcad.contrib.figgear import make_gear_figure
    p = dict(items)
    options = {"bottom_type": "spline"}
    if p["involute_step"] is not None:
        options["involute_step"] = p["involute_step"]
    if p["spline_division_num"] is not None:
        options["spline_division_num"] = p["spline_division_num"]
    points, _ = make_gear_figure(m=p["module_mm"], z=p["teeth"],
                                 alpha_deg=p["pressure_angle_deg"],
                                 **options)
    profile = polygon(points, symmetry=p["teeth"])
    width = p["face_width_mm"]
    body = translate(extrude(profile, width), (0.0, 0.0, 0.5 * width))
    if p["helix_angle_deg"]:
        rate = math.tan(math.radians(p["helix_angle_deg"])) / \
            (0.5 * p["teeth"] * p["module_mm"])
        body = twist(body, rate,
                     fold_z=0.5 * width if p["herringbone"] else None)
    return body


def _spur_items(params):
    return tuple(sorted(params.items()))


def _analyze_spur(params, _children):
    return analyze(_spur_tree(_spur_items(params)))


def _eval_spur(params, p, _children, ev):
    return ev(_spur_tree(_spur_items(params)), p)


register(NodeSpec(
    kind="spur_gear",
    min_children=0,
    max_children=0,
    analyze=_analyze_spur,
    backends={"numpy": _eval_spur},
))


def spur_gear(teeth, module_mm, face_width_mm, pressure_angle_deg=20.0,
              helix_angle_deg=0.0, herringbone=False, involute_step=None,
              spline_division_num=None):
    """An involute gear as a field: spur, or helical/herringbone with a
    nonzero ``helix_angle_deg``.  Axis +z, faces at ``z = 0`` and ``z =
    face_width_mm``; a helix turns counter-clockwise going up for a
    positive angle, matching :func:`yapcad.geom3d_util.helical_extrude`.
    """
    try:
        n = int(teeth)
        m, w, a, h = (float(v) for v in (module_mm, face_width_mm,
                                         pressure_angle_deg, helix_angle_deg))
    except (TypeError, ValueError):
        raise SdfError("spur_gear: dimensions must be numbers") from None
    if n < 6 or m <= 0 or w <= 0 or not 0 < a < 45 or abs(h) >= 60:
        raise SdfError("spur_gear: needs teeth >= 6, positive module and "
                       "width, a pressure angle in (0, 45) and a helix "
                       "angle under 60 degrees")
    return make_node("spur_gear", {
        "teeth": n, "module_mm": m, "face_width_mm": w,
        "pressure_angle_deg": a, "helix_angle_deg": h,
        "herringbone": bool(herringbone),
        "involute_step": None if involute_step is None
        else float(involute_step),
        "spline_division_num": None if spline_division_num is None
        else int(spline_division_num),
    })
