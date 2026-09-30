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
             prune=True, backend=DEFAULT_BACKEND, metadata=None, check=True,
             brep=False):
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
    :param brep: also replay the tree through OCC and attach the exact BREP
        as a *derived* representation, so that
        :func:`yapcad.io.step.write_step_analytic` exports analytic surfaces
        instead of facets.  ``True`` requires a ``csg_exact`` tree and
        pythonocc-core, and raises otherwise; ``"auto"`` attaches one when
        both hold and quietly skips it when not.  Off by default, because a
        BREP payload is megabytes where the tree is kilobytes.
    :returns: a solid whose construction slot holds the generating field.

    The meshing parameters are recorded in the construction record so the
    preview can be regenerated bit-for-bit; see
    :func:`yapcad.sdf.node.to_construction`.
    """
    if brep not in (False, True, "auto"):
        raise SdfError(f"brep must be True, False or 'auto', got {brep!r}")
    replayed = _replay(node, brep)

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
    if replayed is not None:
        _attach_derived_brep(result, replayed, node)
    return result


def _replay(node, mode):
    """Replay ``node`` through OCC per the ``brep`` argument, or ``None``."""
    if mode is False:
        return None
    from yapcad.sdf import occ
    if mode == "auto" and not (occ.occ_available() and node.csg_exact):
        return None
    # Fail before meshing, which is the expensive part.
    return occ.to_brep(node)


def _attach_derived_brep(result, replayed, node):
    """Attach ``replayed`` and tag it with the tree it came from.

    The tag is what lets the serialiser tell a BREP derived from this field
    from a contradictory one.  Anything that later re-attaches a BREP --
    an OCC boolean, a transform -- writes a fresh record without it, and the
    solid then correctly refuses to serialise as SDF-authoritative: the
    BREP no longer derives from the tree.
    """
    from yapcad.brep import attach_brep_to_solid
    from yapcad.metadata import get_solid_metadata
    attach_brep_to_solid(result, replayed)
    get_solid_metadata(result)["brep"]["derivedFrom"] = {"sdf": node.digest}


__all__ = ["mesh_to_surface", "to_solid"]
