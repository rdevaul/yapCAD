# Circular Mounting Plate — 2D Profile
# Interactive circle-radius drag widget demo.
#
# Drag the orange diamond handles on each circle to resize:
#   - Plate radius (outer boundary)
#   - Tap drill radius (5 bolt-circle holes, shared)
#   - Center hole radius
#
# The snap modes constrain drag values to standard hardware dimensions.

module mounting_plate_2d

command mounting_plate_2d(
    plate_dia_mm:   float @ui(widget="circle_r", label="Plate diameter", snap="mm") = 350.0,
    tap_dia_mm:     float @ui(widget="circle_r", label="Tap drill dia", snap="metric_tap") = 3.3,
    center_dia_mm:  float @ui(widget="circle_r", label="Center hole dia", snap="unified_clearance") = 13.5,
    bolt_circle_mm: float = 130.0
) -> region2d:
    plate_r: float = plate_dia_mm / 2.0
    tap_r:   float = tap_dia_mm / 2.0
    clr_r:   float = center_dia_mm / 2.0
    bc:      float = bolt_circle_mm / 2.0

    # Two-pi divided by 5 bolt holes = 72 degrees per hole
    two_pi: float = pi() * 2.0
    step:   float = two_pi / 5.0

    # Outer plate disk
    plate: region2d = disk(point(0.0, 0.0, 0.0), plate_r)

    # 5 evenly-spaced bolt holes on the bolt circle
    h0: region2d = disk(point(bc * cos(0.0 * step), bc * sin(0.0 * step), 0.0), tap_r)
    h1: region2d = disk(point(bc * cos(1.0 * step), bc * sin(1.0 * step), 0.0), tap_r)
    h2: region2d = disk(point(bc * cos(2.0 * step), bc * sin(2.0 * step), 0.0), tap_r)
    h3: region2d = disk(point(bc * cos(3.0 * step), bc * sin(3.0 * step), 0.0), tap_r)
    h4: region2d = disk(point(bc * cos(4.0 * step), bc * sin(4.0 * step), 0.0), tap_r)

    # Center clearance hole
    hc: region2d = disk(point(0.0, 0.0, 0.0), clr_r)

    # Subtract all holes at once
    result: region2d = difference2d_all(plate, [h0, h1, h2, h3, h4, hc])
    emit result
