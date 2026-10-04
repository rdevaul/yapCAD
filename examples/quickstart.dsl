module quickstart

# A mounting bracket: a filleted plate with a bore and two bolt holes.
# Build it as a mesh/BREP solid (the default) or as a signed distance field:
#   python -m yapcad.dsl run examples/quickstart.dsl BRACKET -o bracket.stl
#   python -m yapcad.dsl run examples/quickstart.dsl BRACKET -o bracket.stl \
#       --representation sdf
command BRACKET(width: float = 40.0, depth: float = 30.0,
                thickness: float = 6.0, bore: float = 8.0) -> solid:
    let plate: solid = fillet(box(width, depth, thickness), 1.0)
    let hole: solid = translate(cylinder(bore / 2.0, thickness + 2.0),
                                0.0, 0.0, -thickness / 2.0 - 1.0)
    let bolt: solid = cylinder(2.0, thickness + 2.0)
    let left: solid = translate(bolt, -width / 2.0 + 6.0, 0.0,
                                -thickness / 2.0 - 1.0)
    let right: solid = translate(bolt, width / 2.0 - 6.0, 0.0,
                                 -thickness / 2.0 - 1.0)
    emit difference(plate, hole, left, right)
