# Fillet and Chamfer Demo
# Demonstrates the new fillet() and chamfer() DSL operations
#
# Run with: python -m yapcad.dsl.interpreter examples/fillet_chamfer_demo.dsl
# Or view in OpenGL: python -m yapcad.dsl.interpreter examples/fillet_chamfer_demo.dsl --view

# Create a simple box
base_box: solid = box(30.0, 20.0, 10.0)

# Apply fillet to round all edges
# The fillet radius determines how rounded the edges become
fillet_radius: float = 2.0
filleted_box: solid = fillet(base_box, fillet_radius)

# Move it for display
filleted_display: solid = translate(filleted_box, -20.0, 0.0, 0.0)

# Create another box for chamfer demo
chamfer_box: solid = box(30.0, 20.0, 10.0)

# Apply chamfer to bevel all edges
# The chamfer distance determines the size of the bevel
chamfer_distance: float = 1.5
chamfered_box: solid = chamfer(chamfer_box, chamfer_distance)

# Move it for display
chamfered_display: solid = translate(chamfered_box, 20.0, 0.0, 0.0)

# Combine for output
# Left: filleted box, Right: chamfered box
result: solid = union(filleted_display, chamfered_display)

# Export
export result
