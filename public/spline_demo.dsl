module spline_demo

# Catmull-Rom spline demo
# catmullrom() takes a list of 3D points as its first argument.
# Additional args: closed (bool), alpha (float, 0-1)

command spline_profile() -> catmullrom:
    let pts: list[point] = [
        point(0.0,  0.0, 0.0),
        point(10.0, 5.0, 0.0),
        point(20.0, 2.0, 0.0),
        point(30.0, 8.0, 0.0),
        point(40.0, 0.0, 0.0)
    ]
    let curve: catmullrom = catmullrom(pts)
    emit curve

command demo_box() -> solid:
    let b: solid = box(40.0, 20.0, 10.0)
    emit b
