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
Core Boundary Representation (BREP) helpers for yapCAD.

These classes wrap `pythonocc-core` objects when that dependency is available.
Importing this module on systems without pythonocc-core should not explode;
instead, we raise a clear runtime error the first time a BREP feature is
requested so users know to activate the conda environment.
"""

import math
from typing import Any, Optional

try:  # pragma: no cover - exercised indirectly in environments with OCC
    from OCC.Core.TopoDS import (
        TopoDS_Shape,
        TopoDS_Vertex,
        TopoDS_Edge,
        TopoDS_Face,
        TopoDS_Shell,
        TopoDS_Solid,
    )
    from OCC.Core.BRep import BRep_Tool
    from OCC.Core.TopExp import TopExp_Explorer
    from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_REVERSED
    from OCC.Core.TopLoc import TopLoc_Location
    from OCC.Core.BRepMesh import BRepMesh_IncrementalMesh

    _OCC_IMPORT_ERROR: Optional[Exception] = None
    _HAVE_OCC = True
except ImportError as exc:  # pragma: no cover - handled during runtime detection
    TopoDS_Shape = TopoDS_Vertex = TopoDS_Edge = TopoDS_Face = TopoDS_Shell = TopoDS_Solid = Any
    TopExp_Explorer = TopLoc_Location = BRep_Tool = BRepMesh_IncrementalMesh = Any
    TopAbs_FACE = 0
    TopAbs_REVERSED = -1
    _OCC_IMPORT_ERROR = exc
    _HAVE_OCC = False

from yapcad.geom import point


def occ_available() -> bool:
    """Return True when pythonocc-core imports succeeded."""
    return _HAVE_OCC


def require_occ() -> None:
    """
    Raise a descriptive error if pythonocc-core is not installed/activated.
    """
    if _HAVE_OCC:
        return
    raise RuntimeError(
        "pythonocc-core is not available. Activate the yapcad-brep conda "
        "environment (see docs/BREP_integration_strategy.md) before using "
        "yapCAD's BREP features."
    ) from _OCC_IMPORT_ERROR

class BrepVertex:
    """A wrapper for a TopoDS_Vertex."""
    def __init__(self, shape: TopoDS_Vertex):
        require_occ()
        self._shape = shape

    @property
    def shape(self):
        return self._shape

class BrepEdge:
    """A wrapper for a TopoDS_Edge."""
    def __init__(self, shape: TopoDS_Edge):
        require_occ()
        self._shape = shape

    @property
    def shape(self):
        return self._shape

class BrepFace:
    """A wrapper for a TopoDS_Face."""
    def __init__(self, shape: TopoDS_Face):
        require_occ()
        self._shape = shape

    @property
    def shape(self):
        return self._shape

class BrepSolid:
    """A wrapper for a TopoDS_Solid."""
    def __init__(self, shape: TopoDS_Solid):
        require_occ()
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
        
        require_occ()
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
                orientation = getattr(face, "Orientation", lambda: None)()
                reverse = orientation == TopAbs_REVERSED
                local_vertices = []

                for i in range(1, triangulation.NbNodes() + 1):
                    pnt = triangulation.Node(i).Transformed(trsf)
                    vertex_point = point(pnt.X(), pnt.Y(), pnt.Z())
                    local_vertices.append(vertex_point)
                    all_vertices.append(vertex_point)
                
                has_normals = triangulation.HasNormals()
                local_normals = []
                if has_normals:
                    for i in range(1, triangulation.NbNodes() + 1):
                        n = triangulation.Normal(i)
                        vec = gp_Vec(n.X(), n.Y(), n.Z())
                        vec.Transform(trsf)
                        if reverse:
                            vec = gp_Vec(-vec.X(), -vec.Y(), -vec.Z())
                        local_normals.append(point(vec.X(), vec.Y(), vec.Z()))
                else:
                    accum = [[0.0, 0.0, 0.0] for _ in range(triangulation.NbNodes())]

                for i in range(1, triangulation.NbTriangles() + 1):
                    t = triangulation.Triangle(i)
                    n1, n2, n3 = t.Get()
                    if reverse:
                        n2, n3 = n3, n2
                    all_triangles.append((n1 - 1 + vertex_offset,
                                          n2 - 1 + vertex_offset,
                                          n3 - 1 + vertex_offset))
                    if not has_normals:
                        v0 = local_vertices[n1 - 1]
                        v1 = local_vertices[n2 - 1]
                        v2 = local_vertices[n3 - 1]
                        nx = ((v1[1] - v0[1]) * (v2[2] - v0[2]) -
                              (v1[2] - v0[2]) * (v2[1] - v0[1]))
                        ny = ((v1[2] - v0[2]) * (v2[0] - v0[0]) -
                              (v1[0] - v0[0]) * (v2[2] - v0[2]))
                        nz = ((v1[0] - v0[0]) * (v2[1] - v0[1]) -
                              (v1[1] - v0[1]) * (v2[0] - v0[0]))
                        accum[n1 - 1][0] += nx
                        accum[n1 - 1][1] += ny
                        accum[n1 - 1][2] += nz
                        accum[n2 - 1][0] += nx
                        accum[n2 - 1][1] += ny
                        accum[n2 - 1][2] += nz
                        accum[n3 - 1][0] += nx
                        accum[n3 - 1][1] += ny
                        accum[n3 - 1][2] += nz

                if has_normals:
                    all_normals.extend(local_normals)
                else:
                    for vec in accum:
                        length = math.sqrt(vec[0] ** 2 + vec[1] ** 2 + vec[2] ** 2)
                        if length <= 1e-12:
                            nx, ny, nz = 0.0, 0.0, 1.0
                        else:
                            nx, ny, nz = vec[0] / length, vec[1] / length, vec[2] / length
                        all_normals.append(point(nx, ny, nz))
                
                vertex_offset += triangulation.NbNodes()
            
            explorer.Next()
        
        if not all_normals or len(all_normals) != len(all_vertices):
            default_normal = point(0, 0, 1)
            all_normals = [default_normal for _ in all_vertices]

        triangle_lists = [[tri[0], tri[1], tri[2]] for tri in all_triangles]

        # Boundary/holes placeholders keep the structure compatible with issurface.
        return ['surface', all_vertices, all_normals, triangle_lists, [], []]


def is_brep(obj):
    """Check if an object is a yapCAD BREP object."""
    return isinstance(obj, (BrepVertex, BrepEdge, BrepFace, BrepSolid))


__all__ = [
    "BrepVertex",
    "BrepEdge",
    "BrepFace",
    "BrepSolid",
    "is_brep",
    "occ_available",
    "require_occ",
]
