"""Turn an SDF node DAG into an ordinary yapCAD solid.

This is what makes the rest of yapCAD work on SDF parts: the output is a
plain ``['solid', [surface], material, construction]``, so drawables,
exporters, the collision system and the package writer need no knowledge of
fields at all.

The solid is SDF-authoritative.  Its mesh is a *preview* regenerated from
the field, and the field itself rides in the construction slot that Phase 0
made round-trip, alongside the meshing parameters needed to reproduce the
preview exactly — which matters because yapCAD signs packages.
"""

from yapcad.geom3d import solid, surface
from yapcad.sdf.contour import dual_contour, manifold_defects
from yapcad.sdf.evaluate import DEFAULT_BACKEND
from yapcad.sdf.node import SdfError, to_construction


def mesh_to_surface(mesh):
    """Build a yapCAD surface from a :class:`~yapcad.sdf.contour.Mesh`.

    Vertices become homogeneous points and normals become direction vectors
    with ``w = 0``, matching what the rest of ``geom3d`` produces.  The
    boundary and hole loops are left empty: a dual-contoured surface is
    closed, so it has no boundary.
    """
    vertices = [[float(x), float(y), float(z), 1.0]
                for x, y, z in mesh.vertices]
    normals = [[float(x), float(y), float(z), 0.0] for x, y, z in mesh.normals]
    faces = [[int(a), int(b), int(c)] for a, b, c in mesh.triangles]
    return surface(vertices, normals, faces)


def to_solid(node, resolution=64, bounds=None, padding=None, chunk=32,
             prune=True, backend=DEFAULT_BACKEND, metadata=None, check=True):
    """Mesh ``node`` and return it as a yapCAD solid.

    :param resolution: cells across the longest axis; the dominant quality
        and cost knob.
    :param bounds: region to mesh, defaulting to the node's own bounds.  An
        unbounded field — a bare gyroid, a half space — must be given one,
        or intersected with a solid first.
    :param metadata: optional dictionary merged into the solid's metadata.
    :param check: verify the mesh is a closed manifold before wrapping it in
        a solid.  Leave it on: the downstream — ``volumeof``, the exporters,
        the collision system — assumes closed geometry, and a silently
        broken solid is worse than a refusal.
    :returns: a solid whose construction slot holds the generating field.

    The meshing parameters are recorded in the construction record so the
    preview can be regenerated bit-for-bit; see
    :func:`yapcad.sdf.node.to_construction`.
    """
    mesh = dual_contour(node, resolution=resolution, bounds=bounds,
                        padding=padding, chunk=chunk, prune=prune,
                        backend=backend)
    if check:
        defects = manifold_defects(mesh)
        if any(defects.values()):
            raise SdfError(
                f"dual contouring produced an unusable mesh at resolution "
                f"{resolution}: {defects['nonmanifold']} non-manifold edges, "
                f"{defects['boundary']} boundary edges, "
                f"{defects['coincident']} coincident vertices. The surface "
                f"most likely passes through one cell twice -- a wall or a "
                f"gap thinner than a cell. Raise the resolution, or pass "
                f"check=False to accept the mesh as-is."
            )
    surf = mesh_to_surface(mesh)
    meshing = {
        "method": "dual-contouring",
        "resolution": int(resolution),
    }
    # Only record what was actually specified; a null in the document would
    # have to round-trip as a null, and an absent key says the same thing.
    if padding is not None:
        meshing["padding"] = float(padding)
    if bounds is not None:
        meshing["bounds"] = [float(v) for v in bounds[0]] + \
                            [float(v) for v in bounds[1]]
    record = to_construction(node, meshing=meshing)
    result = solid([surf], [], record)
    if metadata:
        from yapcad.metadata import set_solid_metadata
        set_solid_metadata(result, dict(metadata))
    return result


__all__ = ["mesh_to_surface", "to_solid"]
