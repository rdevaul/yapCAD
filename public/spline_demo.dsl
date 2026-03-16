module spline_demo

# Interactive Catmull-Rom spline demo
# Control points are parameters — adjust them in the parameter panel
# or drag them directly on the 2D canvas

command spline_profile(
    x0: float = 0.0,
    y0: float = 0.0,
    x1: float = 10.0,
    y1: float = 15.0,
    x2: float = 20.0,
    y2: float = 5.0,
    x3: float = 30.0,
    y3: float = 18.0,
    x4: float = 40.0,
    y4: float = 0.0
) -> catmullrom:
    let pts: list[point] = [
        point(x0, y0, 0.0),
        point(x1, y1, 0.0),
        point(x2, y2, 0.0),
        point(x3, y3, 0.0),
        point(x4, y4, 0.0)
    ]
    let curve: catmullrom = catmullrom(pts)
    emit curve

command demo_box(
    width: float  = 40.0,
    depth: float  = 20.0,
    height: float = 10.0
) -> solid:
    let b: solid = box(width, depth, height)
    emit b
