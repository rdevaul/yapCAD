module spline_demo

# Interactive lathe demo.
# Edit the spline profile by dragging control points in the 2D view.
# Click "Set as defaults" to bake the current positions into the DSL source.
# Then evaluate lathe_solid — it calls spline_profile() with no arguments,
# so it always uses the current DSL defaults.
#
# Coordinate convention: X = radius from rotation axis (keep >= 0), Y = height.

command spline_profile(
    p0: point2d @ui(widget="point2d", label="P0 (base)") = point2d(3.0, 0.0),
    p1: point2d @ui(widget="point2d", label="P1") = point2d(8.0, 10.0),
    p2: point2d @ui(widget="point2d", label="P2 (waist)") = point2d(5.0, 20.0),
    p3: point2d @ui(widget="point2d", label="P3") = point2d(9.0, 30.0),
    p4: point2d @ui(widget="point2d", label="P4 (top)") = point2d(3.0, 40.0)
) -> catmullrom:
    let pts: list[point] = [p0, p1, p2, p3, p4]
    let curve: catmullrom = catmullrom(pts)
    emit curve

command lathe_solid() -> solid:
    let profile: catmullrom = spline_profile()
    let region: region2d = region_from_spline(profile)
    let axis: vector = vector(0.0, 0.0, 1.0)
    let result: solid = revolve(region, axis, 360.0)
    emit result
