module spline_demo

# Interactive Catmull-Rom spline demo.
# Each control point is a point2d parameter with widget="point2d".
# Drag handles on the 2D canvas, or edit values in the parameter panel.

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

command demo_box(
    width: float  = 40.0,
    depth: float  = 20.0,
    height: float = 10.0
) -> solid:
    let b: solid = box(width, depth, height)
    emit b
