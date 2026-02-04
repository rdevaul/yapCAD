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

import base64
import math
import os
import tempfile
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
    from OCC.Core.BRep import BRep_Tool, BRep_Builder
    from OCC.Core.TopExp import TopExp_Explorer
    from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_REVERSED
    from OCC.Core.TopLoc import TopLoc_Location
    from OCC.Core.BRepMesh import BRepMesh_IncrementalMesh
    from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_Transform
    from OCC.Core.gp import gp_Trsf, gp_Vec, gp_Dir, gp_Pnt, gp_Ax1, gp_Ax2
    from OCC.Core.BRepTools import breptools
    from OCC.Core.TopoDS import topods

    _OCC_IMPORT_ERROR: Optional[Exception] = None
    _HAVE_OCC = True
except ImportError as exc:  # pragma: no cover - handled during runtime detection
    TopoDS_Shape = TopoDS_Vertex = TopoDS_Edge = TopoDS_Face = TopoDS_Shell = TopoDS_Solid = Any
    TopExp_Explorer = TopLoc_Location = BRep_Tool = BRepMesh_IncrementalMesh = Any
    TopAbs_FACE = 0
    TopAbs_REVERSED = -1
    BRep_Builder = Any
    breptools = None
    gp_Trsf = gp_Vec = gp_Dir = gp_Pnt = gp_Ax1 = gp_Ax2 = None
    topods = None
    _OCC_IMPORT_ERROR = exc
    _HAVE_OCC = False

from yapcad.geom import point
from yapcad.metadata import ensure_solid_id, get_solid_metadata

_BREP_SOLID_CACHE: dict[str, "BrepSolid"] = {}


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


def _shape_to_bytes(shape) -> bytes:
    require_occ()
    if breptools is None:
        raise RuntimeError("BRepTools write support is unavailable")
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".brep")
    tmp.close()
    breptools.Write(shape, tmp.name)
    with open(tmp.name, "rb") as handle:
        data = handle.read()
    os.remove(tmp.name)
    return data


def _shape_from_bytes(data: bytes) -> TopoDS_Shape:
    require_occ()
    if breptools is None:
        raise RuntimeError("BRepTools read support is unavailable")
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".brep")
    tmp.close()
    with open(tmp.name, "wb") as handle:
        handle.write(data)
    builder = BRep_Builder()
    shape = TopoDS_Shape()
    breptools.Read(shape, tmp.name, builder)
    os.remove(tmp.name)
    return shape


def attach_brep_to_solid(solid: list, brep: "BrepSolid") -> None:
    """Embed serialized BREP data into the solid metadata and cache the shape."""
    require_occ()
    meta = get_solid_metadata(solid, create=True)
    solid_id = ensure_solid_id(solid)
    encoded = base64.b64encode(_shape_to_bytes(brep.shape)).decode("ascii")
    meta["brep"] = {
        "encoding": "brep-ascii-base64",
        "data": encoded,
    }
    _BREP_SOLID_CACHE[solid_id] = brep


def brep_from_solid(solid: list) -> Optional["BrepSolid"]:
    """Return the cached BrepSolid for ``solid`` if metadata is present."""
    if not occ_available():
        return None
    meta = get_solid_metadata(solid, create=False)
    if not meta:
        return None
    solid_id = meta.get("entityId")
    if not solid_id:
        return None
    cached = _BREP_SOLID_CACHE.get(solid_id)
    if cached:
        return cached
    brep_info = meta.get("brep")
    if not brep_info:
        return None
    data = brep_info.get("data")
    if not data:
        return None
    try:
        decoded = base64.b64decode(data)
    except Exception:
        return None
    try:
        shape = _shape_from_bytes(decoded)
    except Exception:
        return None
    # Don't force-cast to Solid - shape may be a Compound (from disconnected unions)
    # BrepSolid.tessellate() uses TopExp_Explorer which handles both
    brep = BrepSolid(shape)
    _BREP_SOLID_CACHE[solid_id] = brep
    return brep


def has_brep_data(solid: list) -> bool:
    meta = get_solid_metadata(solid, create=False)
    return bool(meta and meta.get("brep"))


def translate_brep_solid(solid: list, delta) -> None:
    """Apply a translation to the stored BREP shape if present."""
    if not occ_available() or gp_Trsf is None or gp_Vec is None:
        return
    dx = float(delta[0]) if len(delta) > 0 else 0.0
    dy = float(delta[1]) if len(delta) > 1 else 0.0
    dz = float(delta[2]) if len(delta) > 2 else 0.0
    trsf = gp_Trsf()
    trsf.SetTranslation(gp_Vec(dx, dy, dz))
    _apply_trsf_to_brep(solid, trsf)


def rotate_brep_solid(solid: list, ang: float, center, axis) -> None:
    if not occ_available() or gp_Trsf is None or gp_Pnt is None or gp_Ax1 is None or gp_Dir is None:
        return
    ax = axis or point(0, 0, 1)
    magnitude = math.sqrt(ax[0] ** 2 + ax[1] ** 2 + ax[2] ** 2)
    if magnitude <= 1e-12:
        return
    trsf = gp_Trsf()
    direction = gp_Dir(ax[0], ax[1], ax[2])
    trsf.SetRotation(gp_Ax1(gp_Pnt(center[0], center[1], center[2]), direction), math.radians(ang))
    _apply_trsf_to_brep(solid, trsf)


def mirror_brep_solid(solid: list, plane: str) -> None:
    if not occ_available() or gp_Trsf is None or gp_Ax2 is None or gp_Dir is None or gp_Pnt is None:
        return
    mapping = {
        'xy': (0.0, 0.0, 1.0),
        'yz': (1.0, 0.0, 0.0),
        'xz': (0.0, 1.0, 0.0),
    }
    normal = mapping.get(plane)
    if normal is None:
        return
    trsf = gp_Trsf()
    trsf.SetMirror(gp_Ax2(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(*normal)))
    _apply_trsf_to_brep(solid, trsf)


def scale_brep_solid(solid: list, factor: float, center=None) -> None:
    if not occ_available() or gp_Trsf is None or gp_Pnt is None:
        return
    if center is None:
        center = point(0, 0, 0)
    trsf = gp_Trsf()
    trsf.SetScale(gp_Pnt(float(center[0]), float(center[1]), float(center[2])), float(factor))
    _apply_trsf_to_brep(solid, trsf)


def _apply_trsf_to_brep(solid: list, trsf) -> None:
    if not occ_available() or BRepBuilderAPI_Transform is None:
        return
    brep = brep_from_solid(solid)
    if brep is None:
        return
    builder = BRepBuilderAPI_Transform(brep.shape, trsf, True)
    shape = builder.Shape()
    if topods is not None:
        try:
            shape = topods.Solid(shape)
        except Exception:
            pass
    # Regenerate entity ID to avoid cache collision with original solid
    # (deepcopy preserves entity ID, but transformed solid should have its own)
    _regenerate_solid_id(solid)
    attach_brep_to_solid(solid, BrepSolid(shape))


def _regenerate_solid_id(solid: list) -> str:
    """Force a new entity ID for a solid, clearing its BREP cache entry."""
    import uuid
    meta = get_solid_metadata(solid, create=True)
    old_id = meta.get('entityId')
    if old_id and old_id in _BREP_SOLID_CACHE:
        del _BREP_SOLID_CACHE[old_id]
    new_id = str(uuid.uuid4())
    meta['entityId'] = new_id
    meta['id'] = new_id
    return new_id


def transform_brep_shape(brep: "BrepSolid", trsf) -> Optional["BrepSolid"]:
    """Return a transformed BrepSolid (does not touch metadata)."""
    if not occ_available() or BRepBuilderAPI_Transform is None:
        return None
    builder = BRepBuilderAPI_Transform(brep.shape, trsf, True)
    shape = builder.Shape()
    if topods is not None:
        try:
            shape = topods.Solid(shape)
        except Exception:
            pass
    return BrepSolid(shape)


def apply_matrix_to_brep_solid(solid: list, matrix) -> None:
    """Apply a yapCAD transformation matrix to the stored BREP shape.

    This function converts a yapCAD Matrix (4x4 homogeneous transformation)
    to an OCC gp_Trsf and applies it to the BREP data attached to the solid.

    Args:
        solid: A yapCAD solid with BREP metadata attached
        matrix: A yapCAD Matrix instance (from yapcad.xform)

    Note:
        The matrix must represent an affine transformation (rotation, translation,
        uniform scale, or composition thereof). Non-uniform scaling will cause
        the BREP to be dropped.
    """
    if not occ_available() or gp_Trsf is None:
        return

    brep = brep_from_solid(solid)
    if brep is None:
        return

    # Extract transformation parameters from the matrix
    # yapCAD matrices are 4x4 homogeneous transformation matrices
    # Row-major format: m[row][col]
    try:
        m = matrix.m if hasattr(matrix, 'm') else matrix

        # Check for non-uniform scale (approximate check)
        # Extract the 3x3 rotation/scale submatrix
        sx = math.sqrt(m[0][0]**2 + m[1][0]**2 + m[2][0]**2)
        sy = math.sqrt(m[0][1]**2 + m[1][1]**2 + m[2][1]**2)
        sz = math.sqrt(m[0][2]**2 + m[1][2]**2 + m[2][2]**2)

        # If non-uniform scale, we cannot preserve BREP with gp_Trsf
        # (would need gp_GTrsf which changes topology)
        eps = 1e-6
        if abs(sx - sy) > eps or abs(sy - sz) > eps or abs(sx - sz) > eps:
            # Non-uniform scale detected - clear BREP data to avoid stale/inconsistent state
            _clear_brep_data(solid)
            return

        # Build OCC transformation matrix
        # gp_Trsf uses column vectors, so we need to transpose
        trsf = gp_Trsf()

        # Set the transformation values (1-indexed in OCC)
        # OCC gp_Trsf expects: row, col (1-indexed), value
        # The matrix format is:
        #   [R11 R12 R13 Tx]
        #   [R21 R22 R23 Ty]
        #   [R31 R32 R33 Tz]
        #   [0   0   0   1 ]
        trsf.SetValues(
            m[0][0], m[0][1], m[0][2], m[0][3],
            m[1][0], m[1][1], m[1][2], m[1][3],
            m[2][0], m[2][1], m[2][2], m[2][3]
        )

        _apply_trsf_to_brep(solid, trsf)
    except Exception:
        # If matrix conversion fails, BREP is not updated
        pass


def _clear_brep_data(solid: list) -> None:
    """Remove BREP data from a solid to avoid stale/inconsistent state."""
    meta = get_solid_metadata(solid, create=False)
    if meta is None:
        return
    solid_id = meta.get('entityId')
    if solid_id and solid_id in _BREP_SOLID_CACHE:
        del _BREP_SOLID_CACHE[solid_id]
    if 'brep' in meta:
        del meta['brep']


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
    """A wrapper for a TopoDS_Shape (Solid or Compound of Solids)."""
    def __init__(self, shape):
        require_occ()
        self._shape = shape

    @property
    def shape(self):
        return self._shape

    def tessellate(self, deflection=0.5, _debug=False):
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

        if _debug:
            print(f"[tessellate] Shape type: {self._shape.ShapeType()}")

        explorer = TopExp_Explorer(self._shape, TopAbs_FACE)
        _face_count = 0
        while explorer.More():
            _face_count += 1
            # Cast to TopoDS_Face for OCP compatibility (pythonocc may not need this)
            face = topods.Face(explorer.Current())
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
                        local_normals.append([vec.X(), vec.Y(), vec.Z(), 0.0])
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
                        all_normals.append([nx, ny, nz, 0.0])
                
                vertex_offset += triangulation.NbNodes()
            
            explorer.Next()

        if _debug:
            print(f"[tessellate] Faces: {_face_count}, Vertices: {len(all_vertices)}")

        if not all_normals or len(all_normals) != len(all_vertices):
            default_normal = [0.0, 0.0, 1.0, 0.0]
            all_normals = [default_normal for _ in all_vertices]

        triangle_lists = [[tri[0], tri[1], tri[2]] for tri in all_triangles]

        # Boundary/holes placeholders keep the structure compatible with issurface.
        return ['surface', all_vertices, all_normals, triangle_lists, [], []]


# =============================================================================
# Fillet and Chamfer Operations
# =============================================================================

def _get_all_edges(shape) -> list:
    """Extract all edges from a TopoDS_Shape.

    Returns a list of TopoDS_Edge objects.
    """
    require_occ()
    from OCC.Core.TopAbs import TopAbs_EDGE

    edges = []
    explorer = TopExp_Explorer(shape, TopAbs_EDGE)
    while explorer.More():
        edge = explorer.Current()
        edges.append(edge)
        explorer.Next()
    return edges


def fillet_all_edges(brep_solid: BrepSolid, radius: float) -> BrepSolid:
    """Apply fillet (rounded edge) to all edges of a BREP solid.

    Args:
        brep_solid: A BrepSolid object containing the shape to fillet
        radius: Fillet radius in model units

    Returns:
        A new BrepSolid with filleted edges

    Raises:
        RuntimeError: If OCC is not available or fillet operation fails
    """
    require_occ()

    from OCC.Core.BRepFilletAPI import BRepFilletAPI_MakeFillet
    from OCC.Core.TopAbs import TopAbs_EDGE

    shape = brep_solid.shape
    fillet_maker = BRepFilletAPI_MakeFillet(shape)

    # Add all edges with the specified radius
    edges = _get_all_edges(shape)
    if not edges:
        # No edges to fillet, return original
        return brep_solid

    for edge in edges:
        # Cast to TopoDS_Edge if needed
        if topods is not None:
            try:
                edge = topods.Edge(edge)
            except Exception:
                continue
        try:
            fillet_maker.Add(radius, edge)
        except Exception:
            # Skip edges that can't be filleted (e.g., too small)
            continue

    # Build the filleted shape
    try:
        fillet_maker.Build()
        if not fillet_maker.IsDone():
            raise RuntimeError("Fillet operation failed - IsDone() returned False")
        result_shape = fillet_maker.Shape()
    except Exception as e:
        raise RuntimeError(f"Fillet operation failed: {e}")

    return BrepSolid(result_shape)


def chamfer_all_edges(brep_solid: BrepSolid, distance: float) -> BrepSolid:
    """Apply chamfer (beveled edge) to all edges of a BREP solid.

    Creates symmetric 45-degree chamfers with equal distances on both faces.

    Args:
        brep_solid: A BrepSolid object containing the shape to chamfer
        distance: Chamfer distance from the edge in model units

    Returns:
        A new BrepSolid with chamfered edges

    Raises:
        RuntimeError: If OCC is not available or chamfer operation fails
    """
    require_occ()

    from OCC.Core.BRepFilletAPI import BRepFilletAPI_MakeChamfer

    shape = brep_solid.shape
    chamfer_maker = BRepFilletAPI_MakeChamfer(shape)

    edges = _get_all_edges(shape)
    if not edges:
        # No edges to chamfer, return original
        return brep_solid

    for edge in edges:
        # Cast to TopoDS_Edge if needed
        if topods is not None:
            try:
                edge = topods.Edge(edge)
            except Exception:
                continue

        try:
            # Add symmetric chamfer (single distance applies to both faces equally)
            chamfer_maker.Add(distance, edge)
        except Exception:
            # Skip edges that can't be chamfered
            continue

    # Build the chamfered shape
    try:
        chamfer_maker.Build()
        if not chamfer_maker.IsDone():
            raise RuntimeError("Chamfer operation failed - IsDone() returned False")
        result_shape = chamfer_maker.Shape()
    except Exception as e:
        raise RuntimeError(f"Chamfer operation failed: {e}")

    return BrepSolid(result_shape)


def fillet_edges(brep_solid: BrepSolid, edges: list, radius: float) -> BrepSolid:
    """Apply fillet to specific edges of a BREP solid.

    Args:
        brep_solid: A BrepSolid object containing the shape to fillet
        edges: List of BrepEdge objects to fillet
        radius: Fillet radius in model units

    Returns:
        A new BrepSolid with filleted edges

    Raises:
        RuntimeError: If OCC is not available or fillet operation fails
    """
    require_occ()

    from OCC.Core.BRepFilletAPI import BRepFilletAPI_MakeFillet

    shape = brep_solid.shape
    fillet_maker = BRepFilletAPI_MakeFillet(shape)

    if not edges:
        return brep_solid

    for brep_edge in edges:
        edge = brep_edge.shape if isinstance(brep_edge, BrepEdge) else brep_edge
        try:
            fillet_maker.Add(radius, edge)
        except Exception:
            continue

    try:
        fillet_maker.Build()
        if not fillet_maker.IsDone():
            raise RuntimeError("Fillet operation failed - IsDone() returned False")
        result_shape = fillet_maker.Shape()
    except Exception as e:
        raise RuntimeError(f"Fillet operation failed: {e}")

    return BrepSolid(result_shape)


def chamfer_edges(brep_solid: BrepSolid, edges: list, distance: float) -> BrepSolid:
    """Apply chamfer to specific edges of a BREP solid.

    Args:
        brep_solid: A BrepSolid object containing the shape to chamfer
        edges: List of BrepEdge objects to chamfer
        distance: Chamfer distance from the edge in model units

    Returns:
        A new BrepSolid with chamfered edges

    Raises:
        RuntimeError: If OCC is not available or chamfer operation fails
    """
    require_occ()

    from OCC.Core.BRepFilletAPI import BRepFilletAPI_MakeChamfer

    shape = brep_solid.shape
    chamfer_maker = BRepFilletAPI_MakeChamfer(shape)

    if not edges:
        return brep_solid

    for brep_edge in edges:
        edge = brep_edge.shape if isinstance(brep_edge, BrepEdge) else brep_edge
        try:
            # Add symmetric chamfer (single distance)
            chamfer_maker.Add(distance, edge)
        except Exception:
            continue

    try:
        chamfer_maker.Build()
        if not chamfer_maker.IsDone():
            raise RuntimeError("Chamfer operation failed - IsDone() returned False")
        result_shape = chamfer_maker.Shape()
    except Exception as e:
        raise RuntimeError(f"Chamfer operation failed: {e}")

    return BrepSolid(result_shape)


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
    "attach_brep_to_solid",
    "brep_from_solid",
    "has_brep_data",
    "apply_matrix_to_brep_solid",
    "fillet_all_edges",
    "chamfer_all_edges",
    "fillet_edges",
    "chamfer_edges",
]
