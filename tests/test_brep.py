from yapcad.geometry import Geometry
from yapcad.brep import BrepSolid
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox

def test_brep_tessellate():
    """
    Test that a BrepSolid can be tessellated into a yapCAD surface.
    """
    # Create a simple BREP box
    box_shape = BRepPrimAPI_MakeBox(10, 10, 10).Shape()
    brep_solid = BrepSolid(box_shape)

    # Wrap it in a Geometry object
    geo = Geometry(brep_solid)

    # Trigger tessellation
    surface = geo.surface()

    # Assert that the returned surface is a valid yapCAD surface
    assert isinstance(surface, list)
    assert len(surface) == 4
    assert surface[0] == 'surface'
    assert isinstance(surface[1], list) # vertices
    assert isinstance(surface[2], list) # normals
    assert isinstance(surface[3], list) # faces
