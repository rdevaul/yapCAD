## yapCAD BREP representation
## =====================================

## Copyright (c) 2025 Richard W. DeVaul
## Copyright (c) 2025 yapCAD contributors
## All rights reserved

# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation files
# (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
# BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""
This module provides the core classes for the Boundary Representation (BREP)
model in yapCAD. These classes wrap the underlying `pythonocc-core` objects
to provide a more Pythonic and yapCAD-integrated interface.
"""

from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Vertex, TopoDS_Edge, TopoDS_Face, TopoDS_Shell, TopoDS_Solid
from OCC.Core.BRep import BRep_Tool
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.BRepMesh import BRepMesh_IncrementalMesh

from yapcad.geom import point

class BrepVertex:
    """A wrapper for a TopoDS_Vertex."""
    def __init__(self, shape: TopoDS_Vertex):
        self._shape = shape

    @property
    def shape(self):
        return self._shape

class BrepEdge:
    """A wrapper for a TopoDS_Edge."""
    def __init__(self, shape: TopoDS_Edge):
        self._shape = shape

    @property
    def shape(self):
        return self._shape

class BrepFace:
    """A wrapper for a TopoDS_Face."""
    def __init__(self, shape: TopoDS_Face):
        self._shape = shape

    @property
    def shape(self):
        return self._shape

class BrepSolid:
    """A wrapper for a TopoDS_Solid."""
    def __init__(self, shape: TopoDS_Solid):
        self._shape = shape

    @property
    def shape(self):
        return self._shape

    def tessellate(self, deflection=0.5):
        """
        Generate a faceted representation of the BREP model.

        This method will use `pythonocc-core`'s meshing capabilities to
        generate a triangular mesh of the BREP model. The resulting
        faceted representation will be returned in the same format as
        `yapcad.geom3d.poly2surface`.
        """
        
        mesh = BRepMesh_IncrementalMesh(self._shape, deflection)
        
        all_vertices = []
        all_triangles = []
        all_normals = []
        vertex_offset = 0

        explorer = TopExp_Explorer(self._shape, TopAbs_FACE)
        while explorer.More():
            face = explorer.Current()
            loc = TopLoc_Location()
            
            triangulation = BRep_Tool.Triangulation(face, loc)

            if triangulation:
                trsf = loc.Transformation()
                
                for i in range(1, triangulation.NbNodes() + 1):
                    pnt = triangulation.Node(i).Transformed(trsf)
                    all_vertices.append(point(pnt.X(), pnt.Y(), pnt.Z()))
                
                if triangulation.HasNormals():
                    for i in range(1, triangulation.NbNodes() + 1):
                        n = triangulation.Normal(i)
                        all_normals.append(point(n.X(), n.Y(), n.Z()))

                for i in range(1, triangulation.NbTriangles() + 1):
                    t = triangulation.Triangle(i)
                    n1, n2, n3 = t.Get()
                    all_triangles.append((n1 - 1 + vertex_offset, n2 - 1 + vertex_offset, n3 - 1 + vertex_offset))
                
                vertex_offset += triangulation.NbNodes()
            
            explorer.Next()
        
        return ['surface', all_vertices, all_normals, all_triangles]


def is_brep(obj):
    """Check if an object is a yapCAD BREP object."""
    return isinstance(obj, (BrepVertex, BrepEdge, BrepFace, BrepSolid))