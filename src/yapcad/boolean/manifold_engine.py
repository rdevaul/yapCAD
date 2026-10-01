"""Mesh booleans through the manifold3d library, called directly.

`manifold <https://github.com/elalish/manifold>`_ is a robust, exact-
topology mesh boolean library: given closed, consistently oriented meshes it
returns closed, consistently oriented meshes, coplanar and touching faces
included.  This engine calls its Python bindings directly rather than going
through trimesh, so the only addition is the single ``manifold3d`` wheel.

It works in float64 throughout via ``Mesh64``; the older float32 ``Mesh``
would quietly round coordinates to about seven significant digits, so a
manifold3d without ``Mesh64`` is treated as unavailable rather than used at
reduced precision.

Manifold refuses inputs that are not closed 2-manifolds -- an STL with a
hole, a self-touching edge.  :func:`solid_boolean` raises
:class:`NotManifoldError` for those, naming the operand, instead of guessing.
"""

import numpy as np

from yapcad.geom import epsilon

try:  # pragma: no cover - optional dependency
    import manifold3d as _manifold
except ImportError:  # pragma: no cover - optional dependency
    _manifold = None

ENGINE_NAME = "manifold"

_OPERATORS = {
    "union": lambda a, b: a + b,
    "difference": lambda a, b: a - b,
    "intersection": lambda a, b: a ^ b,
}


class NotManifoldError(ValueError):
    """An operand is not a closed, consistently oriented 2-manifold."""


def is_available():
    """True when manifold3d is installed with float64 (``Mesh64``) support."""
    return _manifold is not None and hasattr(_manifold, "Mesh64")


def _require():
    if _manifold is None:
        raise RuntimeError(
            "the manifold boolean engine needs manifold3d "
            "(pip install manifold3d, or the yapCAD 'manifold' extra)"
        )
    if not hasattr(_manifold, "Mesh64"):
        raise RuntimeError(
            "the installed manifold3d has no float64 Mesh64; upgrade to "
            "manifold3d>=3.0 rather than run booleans at float32 precision"
        )


def _raw_triangles(solid):
    """Every triangle of a solid, unfiltered.

    yapcad.mesh.mesh_view skips faces with area <= epsilon, which can open a
    closed input, so surfaces are read directly.  A BREP-only solid, which
    has no surfaces, is tessellated through mesh_view since that is the only
    mesh it has.
    """
    tris = []
    for surf in solid[1]:
        verts = surf[1]
        for face in surf[3]:
            if len(face) == 3:
                tris.append((verts[face[0]], verts[face[1]], verts[face[2]]))
    if tris or solid[1]:
        return tris
    from yapcad.mesh import mesh_view
    return [(v0, v1, v2) for _, v0, v1, v2 in mesh_view(solid)]


def _to_manifold(solid, role):
    """Weld a yapCAD solid into an indexed mesh and build a Manifold.

    Vertices are merged with yapCAD's own position key -- the quantisation
    geom3d.issolidclosed uses -- so "the same vertex" means the same thing
    here as everywhere else in yapCAD.  The first occurrence supplies the
    coordinates, which keeps the result deterministic.
    """
    tris = _raw_triangles(solid)
    if not tris:
        return _manifold.Manifold()
    pts = np.asarray([[p[0], p[1], p[2]] for tri in tris for p in tri],
                     dtype=np.float64)
    keys = np.round(pts / epsilon).astype(np.int64)
    _, first, inverse = np.unique(keys, axis=0, return_index=True,
                                  return_inverse=True)
    inverse = inverse.ravel()
    order = np.argsort(first, kind="stable")      # keep first-seen order
    remap = np.empty_like(order)
    remap[order] = np.arange(len(order))
    verts = pts[first[order]]
    faces = remap[inverse].reshape(-1, 3)
    # Triangles that collapse when welded are degenerate, not geometry.
    keep = ((faces[:, 0] != faces[:, 1]) & (faces[:, 1] != faces[:, 2])
            & (faces[:, 0] != faces[:, 2]))
    faces = faces[keep]

    mesh = _manifold.Mesh64(
        vert_properties=np.ascontiguousarray(verts),
        tri_verts=np.ascontiguousarray(faces.astype(np.uint64)),
    )
    result = _manifold.Manifold(mesh)
    status = result.status()
    if status != _manifold.Error.NoError:
        raise NotManifoldError(
            f"the {role} operand is not a closed, consistently oriented "
            f"mesh (manifold3d: {status.name}); repair it first, e.g. with "
            "yapcad.io.stl's mesh checks or pymeshfix"
        )
    return result


def _from_manifold(result, operation):
    """A yapCAD solid from a Manifold: flat-shaded, nothing dropped.

    Built as a triangle soup with one face normal per triangle, the same
    shape the other engines return, so hard edges shade as edges.  Unlike
    native._surface_from_triangles it drops no triangle, small or not: the
    output is manifold by construction and dropping any face would open it.
    """
    from yapcad.geom import point
    from yapcad.geom3d import solid, surface

    record = ["boolean", f"{ENGINE_NAME}:{operation}"]
    if result.is_empty():
        return solid([], [], record)
    mesh = result.to_mesh64()
    verts = np.asarray(mesh.vert_properties, dtype=np.float64)[:, :3]
    faces = np.asarray(mesh.tri_verts, dtype=np.int64)
    corners = verts[faces]                                  # (T, 3, 3)
    normal = np.cross(corners[:, 1] - corners[:, 0],
                      corners[:, 2] - corners[:, 0])
    length = np.linalg.norm(normal, axis=1, keepdims=True)
    normal = np.divide(normal, length, out=np.zeros_like(normal),
                       where=length > 0)

    out_verts, out_norms, out_faces = [], [], []
    for t in range(len(faces)):
        n = [float(normal[t, 0]), float(normal[t, 1]), float(normal[t, 2]),
             0.0]
        base = len(out_verts)
        for k in range(3):
            c = corners[t, k]
            out_verts.append(point(float(c[0]), float(c[1]), float(c[2])))
            out_norms.append(list(n))
        out_faces.append([base, base + 1, base + 2])
    return solid([surface(out_verts, out_norms, out_faces)], [], record)


def solid_boolean(a, b, operation):
    """Union, intersection or difference of two solids through manifold3d.

    :raises NotManifoldError: when an operand is not a closed 2-manifold.
    """
    _require()
    op = operation.lower()
    if op not in _OPERATORS:
        raise ValueError(f"unsupported solid boolean operation {operation!r}")
    ma = _to_manifold(a, "first")
    mb = _to_manifold(b, "second")
    return _from_manifold(_OPERATORS[op](ma, mb), op)


__all__ = ["ENGINE_NAME", "NotManifoldError", "is_available",
           "solid_boolean"]
