module spline_demo

# Shared control point defaults — same values in both commands.
# Edit p0..p4 in either command's parameter panel; both commands
# use the same geometry. When you update params and evaluate
# lathe_solid, it uses whatever p0..p4 you set.

command spline_profile(
    p0: point2d @ui(widget="point2d", label="P0") = point2d(5.0, 0.0),
    p1: point2d @ui(widget="point2d", label="P1") = point2d(8.0, 10.0),
    p2: point2d @ui(widget="point2d", label="P2") = point2d(6.0, 20.0),
    p3: point2d @ui(widget="point2d", label="P3") = point2d(9.0, 30.0),
    p4: point2d @ui(widget="point2d", label="P4") = point2d(5.0, 40.0)
) -> catmullrom:
    let pts: list[point] = [p0, p1, p2, p3, p4]
    let curve: catmullrom = catmullrom(pts)
    emit curve

# lathe_solid takes the same p0..p4 parameters as spline_profile.
# The workbench shares parameter values across commands in the same file,
# so updating p0..p4 in one panel updates both.
command lathe_solid(
    p0: point2d @ui(widget="point2d", label="P0") = point2d(5.0, 0.0),
    p1: point2d @ui(widget="point2d", label="P1") = point2d(8.0, 10.0),
    p2: point2d @ui(widget="point2d", label="P2") = point2d(6.0, 20.0),
    p3: point2d @ui(widget="point2d", label="P3") = point2d(9.0, 30.0),
    p4: point2d @ui(widget="point2d", label="P4") = point2d(5.0, 40.0)
) -> solid:
    let pts: list[point] = [p0, p1, p2, p3, p4]
    let profile: catmullrom = catmullrom(pts)
    let region: region2d = region_from_spline(profile)
    let axis: vector = vector(0.0, 0.0, 1.0)
    let result: solid = revolve(region, axis, 360.0)
    emit result
