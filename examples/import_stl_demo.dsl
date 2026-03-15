module import_stl_demo;

# =============================================================================
# STL Import — native Python block wrapper
#
# Imports an external STL file into a yapCAD solid.
# Both binary and ASCII STL formats are supported.
#
# Commands:
#   IMPORT_STL(path)              — import at original scale/position
#   IMPORT_STL_SCALED(path, scale) — apply uniform scale (e.g. 25.4 for in→mm)
#
# path must be an absolute filesystem path on the server.
#
# Author: Jarvis / Dark Matter Lab  2026-03-06
# =============================================================================

native python {
from yapcad.io.stl import import_stl
from yapcad.geom3d import scale as _scale_solid

def load_stl_solid(path):
    """Load an STL file and return a yapCAD solid."""
    sld = import_stl(str(path))
    if not sld or not sld[1]:
        raise ValueError(f"No geometry found in STL file: {path}")
    return sld

def load_stl_solid_scaled(path, scale):
    """Load an STL file and apply a uniform scale factor."""
    sld = load_stl_solid(path)
    if scale != 1.0:
        sld = _scale_solid(sld, float(scale))
    return sld

} exports {
    fn load_stl_solid(path: string) -> solid;
    fn load_stl_solid_scaled(path: string, scale: float) -> solid;
}


command IMPORT_STL(path: string) -> solid:
    result: solid = load_stl_solid(path)
    emit result


command IMPORT_STL_SCALED(
    path:  string,
    scale: float = 1.0
) -> solid:
    result: solid = load_stl_solid_scaled(path, scale)
    emit result
