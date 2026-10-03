"""Involute spur, helical and herringbone gears as SDF-authored solids."""

import math


def make_involute_gear_sdf(teeth, module_mm, face_width_mm, *,
                           pressure_angle_deg=20.0, helix_angle_deg=0.0,
                           herringbone=False, involute_step=None,
                           spline_division_num=None, cell_mm=None):
    """Mesh :func:`yapcad.sdf.spur_gear` into a closed solid; no OCC needed.

    The profile is :func:`yapcad.contrib.figgear.make_gear_figure`'s, with
    the same options the DSL builtins pass it; axis +z, faces at ``z = 0``
    and ``z = face_width_mm``.  ``cell_mm`` defaults to an eighth of the
    module, fine enough for the tooth tips.
    """
    from yapcad.sdf.convert import to_solid
    from yapcad.sdf.gears import spur_gear

    node = spur_gear(teeth, module_mm, face_width_mm,
                     pressure_angle_deg=pressure_angle_deg,
                     helix_angle_deg=helix_angle_deg, herringbone=herringbone,
                     involute_step=involute_step,
                     spline_division_num=spline_division_num)
    cell = cell_mm or module_mm / 8.0
    lo, hi = node.bounds
    longest = max(hi[i] - lo[i] for i in range(3))
    kind = "herringbone" if herringbone else (
        "helical" if helix_angle_deg else "spur")
    return to_solid(node, resolution=max(16, math.ceil(longest / cell)),
                    metadata={"gear": {"kind": kind, "teeth": int(teeth),
                                       "module_mm": float(module_mm),
                                       "helix_angle_deg":
                                           float(helix_angle_deg)}})
