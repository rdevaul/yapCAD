module mounting_plate

# =============================================================================
# Circular mounting plate — configurable OD
#
# 2D layout:  mounting_plate_2d()  — exact circle sketch in the 2D pane
#             All circles are native arc primitives tagged with their param name.
#             Drag handles appear at the 3-o'clock position of each circle.
#             Drag to resize plate, bolt circle, or holes live.
#
# 3D solid:   mounting_plate_3d()  — plate with through-bores
#
# Metric hole sizes:  "M3" "M4" "M5" "M6" "M8" "M10"
# Unified center:     "1/4-20"  "5/16-18"  "3/8-16"  "1/2-13"  "1/2-20"
# Pentagon starts at 12 o'clock (90°), one hole per 72°
# =============================================================================

native python {
import math

# Tap-drill diameters (mm) — hole is drilled to this size before tapping
TAP_DRILL = {
    "M3": 2.5, "M4": 3.3, "M5": 4.2,
    "M6": 5.0, "M8": 6.75, "M10": 8.5,
}
# Clearance bores (mm) for unified fasteners
CLEARANCE = {
    "1/4-20": 6.8, "5/16-18": 8.4, "3/8-16": 10.0,
    "1/2-13": 13.5, "1/2-20": 13.5,
}

def plate_circles(plate_dia_mm, bolt_circle_mm, metric_size, unified_size):
    """Return a geomlist of tagged native circle primitives."""
    from yapcad.geom import arc, point as mkpoint

    plate_r   = plate_dia_mm / 2.0
    bc_r      = bolt_circle_mm / 2.0
    tap_r     = TAP_DRILL[metric_size] / 2.0
    center_r  = CLEARANCE[unified_size] / 2.0

    # Plate outline
    circles = [arc(mkpoint(0, 0), plate_r,   meta={"param": "plate_dia_mm",   "label": "Plate OD"})]

    # Pentagon of metric tapped holes — 12 o'clock, 72° steps
    for i in range(5):
        theta = math.radians(90.0 + i * 72.0)
        cx = bc_r * math.cos(theta)
        cy = bc_r * math.sin(theta)
        circles.append(arc(mkpoint(cx, cy), tap_r,
                           meta={"param": "metric_hole_r", "label": "M-hole {}".format(i+1)}))

    # Bolt circle (dashed in 2D — same primitive, viewer can distinguish by param name)
    circles.append(arc(mkpoint(0, 0), bc_r,
                       meta={"param": "bolt_circle_mm", "label": "Bolt Circle"}))

    # Center clearance hole
    circles.append(arc(mkpoint(0, 0), center_r,
                       meta={"param": "center_hole_r", "label": "Center Hole"}))

    return circles

def plate_solid(plate_dia_mm, bolt_circle_mm, thickness_mm, metric_size, unified_size):
    """Return a 3D solid mounting plate with through-bores."""
    from yapcad.geom import point as mkpoint
    from yapcad.geom3d import poly2surfaceXY
    from yapcad.geom3d_util import extrude

    plate_r   = plate_dia_mm / 2.0
    bc_r      = bolt_circle_mm / 2.0
    tap_r     = TAP_DRILL[metric_size] / 2.0
    center_r  = CLEARANCE[unified_size] / 2.0

    def circle_poly(cx, cy, r, n=64):
        pts = []
        for i in range(n):
            a = 2 * math.pi * i / n
            pts.append(mkpoint(cx + r*math.cos(a), cy + r*math.sin(a)))
        pts.append(pts[0])
        return pts

    # Build plate as extrusion minus bores
    plate_poly = circle_poly(0, 0, plate_r)
    surf = poly2surfaceXY(plate_poly)
    base = extrude(surf, thickness_mm)

    # Subtract metric holes
    from yapcad.boolean import difference
    for i in range(5):
        theta = math.radians(90.0 + i * 72.0)
        cx = bc_r * math.cos(theta)
        cy = bc_r * math.sin(theta)
        hole_poly = circle_poly(cx, cy, tap_r)
        hs = poly2surfaceXY(hole_poly)
        hole = extrude(hs, thickness_mm + 2)
        base = difference(base, hole)

    # Subtract center hole
    ch_poly = circle_poly(0, 0, center_r)
    cs = poly2surfaceXY(ch_poly)
    ch = extrude(cs, thickness_mm + 2)
    base = difference(base, ch)

    return base

} exports {
    fn plate_circles(plate_dia_mm: float, bolt_circle_mm: float, metric_size: string, unified_size: string) -> list;
    fn plate_solid(plate_dia_mm: float, bolt_circle_mm: float, thickness_mm: float, metric_size: string, unified_size: string) -> solid;
}


command mounting_plate_2d(
    plate_dia_mm:   float @ui(widget="circle_r", snap="mm", unit="diameter", label="Plate OD (mm)")    = 350.0,
    bolt_circle_mm: float @ui(widget="circle_r", snap="mm", unit="diameter", label="Bolt Circle (mm)") = 130.0,
    metric_size:    string                                                                         = "M4",
    unified_size:   string                                                                         = "1/2-13",
) -> list:
    result: list = plate_circles(plate_dia_mm, bolt_circle_mm, metric_size, unified_size)
    emit result


command mounting_plate_3d(
    plate_dia_mm:    float = 350.0,
    bolt_circle_mm:  float = 130.0,
    thickness_mm:    float = 10.0,
    metric_size:     string = "M4",
    unified_size:    string = "1/2-13",
) -> solid:
    result: solid = plate_solid(plate_dia_mm, bolt_circle_mm, thickness_mm, metric_size, unified_size)
    emit result
