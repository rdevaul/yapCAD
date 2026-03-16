module spline_demo

# Interactive Catmull-Rom spline / lathe demo.
#
# The 2D profile is in the (radius, height) plane:
#   X axis = radius from the rotation axis (must stay >= 0)
#   Y axis = height along the rotation axis
#
# Editing p0..p4 in spline_profile updates the 2D view live.
# Switching to lathe_solid and evaluating uses the same parameters
# to generate the solid of revolution.

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

command lathe_solid(
    p0: point2d @ui(widget="point2d", label="P0 (base)") = point2d(3.0, 0.0),
    p1: point2d @ui(widget="point2d", label="P1") = point2d(8.0, 10.0),
    p2: point2d @ui(widget="point2d", label="P2 (waist)") = point2d(5.0, 20.0),
    p3: point2d @ui(widget="point2d", label="P3") = point2d(9.0, 30.0),
    p4: point2d @ui(widget="point2d", label="P4 (top)") = point2d(3.0, 40.0)
) -> solid:
    let pts: list[point] = [p0, p1, p2, p3, p4]
    let profile: catmullrom = catmullrom(pts)
    let region: region2d = region_from_spline(profile)
    let axis: vector = vector(0.0, 0.0, 1.0)
    let result: solid = revolve(region, axis, 360.0)
    emit result
