module spline_demo

# Interactive Catmull-Rom spline demo.
# spline_profile: drag control points on the 2D canvas or edit in the parameter panel.
# lathe_solid: revolves the same profile around the Y axis to produce a solid of revolution.

command spline_profile(
    p0: point2d @ui(widget="point2d", label="P0") = point2d(0.0, 0.0),
    p1: point2d @ui(widget="point2d", label="P1") = point2d(10.0, 15.0),
    p2: point2d @ui(widget="point2d", label="P2") = point2d(20.0, 5.0),
    p3: point2d @ui(widget="point2d", label="P3") = point2d(30.0, 18.0),
    p4: point2d @ui(widget="point2d", label="P4") = point2d(40.0, 0.0)
) -> catmullrom:
    let pts: list[point] = [p0, p1, p2, p3, p4]
    let curve: catmullrom = catmullrom(pts)
    emit curve

# Solid of revolution: control points define the lathe profile.
# p0.x is radius at the base, p4.x is radius at the top.
# The profile is revolved 360 degrees around the vertical axis.
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

command demo_box(
    width: float  = 40.0,
    depth: float  = 20.0,
    height: float = 10.0
) -> solid:
    let b: solid = box(width, depth, height)
    emit b
