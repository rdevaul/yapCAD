# Sample yapCAD DSL Code
# This is a simple example that creates a cube

# Create a basic cube
let my_cube = cube(10, 10, 10);

# Translate it to a different position  
let translated_cube = translate(my_cube, 5, 5, 5);

# Create a sphere
let my_sphere = sphere(5);

# Combine them using union
let combined = union(translated_cube, my_sphere);

# Emit the final result
emit combined;