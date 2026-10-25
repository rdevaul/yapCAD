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

import math
import warnings

from yapcad.geom3d import solid, surface
from yapcad.sdf.contour import (
    dual_contour,
    manifold_defects,
    proves_empty,
)
from yapcad.sdf.evaluate import DEFAULT_BACKEND
from yapcad.sdf.node import (
    SdfError,
    bounds_intersect,
    bounds_is_empty,
    make_node,
    to_construction,
)


class NonManifoldMeshError(SdfError):
    """Raised by :func:`to_solid` when the mesh fails its manifold check.

    A subclass of :class:`~yapcad.sdf.node.SdfError`, so existing handlers
    still catch it; separate so that a caller able to respond -- by meshing
    again at a higher resolution -- can catch exactly this and nothing else.
    """


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


#: Cells across the longest axis that automatic refinement will not pass.
#: The corner lattice is ``(n+1) ** 3`` float64 values, so 256 is about
#: 135 MB and 512 would be over a gigabyte.  :mod:`yapcad.sdf.booleans`
#: starts no finer than this for the same reason.
MAX_RESOLUTION = 256


def to_solid(node, resolution=64, bounds=None, padding=None, chunk=32,
             prune=True, backend=DEFAULT_BACKEND, metadata=None, check=True,
             brep=False, cell_mm=None, refine=True,
             refine_limit=MAX_RESOLUTION):
    """Mesh ``node`` and return it as a yapCAD solid.

    :param resolution: cells across the longest axis; the dominant quality
        and cost knob.
    :param bounds: region to mesh, defaulting to the node's own bounds.  An
        unbounded field — a bare gyroid, a half space — must be given one,
        or intersected with a solid first.
    :param cell_mm: the cell size to mesh at, in model units, instead of
        deriving one from ``resolution``.  This is the hint the refinement
        warning hands back: pass it and the part meshes in one pass.  It is
        rounded down to a whole number of cells across the longest axis, so
        the cell used is never coarser than asked.
    :param refine: when the cell asked for is too coarse to find the
        surface — an empty mesh, or one dual contouring cannot close —
        halve it and mesh again, up to ``refine_limit`` cells across, and
        warn with the cell it settled on.  That cell is the coarsest one
        that *found* the surface, which is not the same as one that
        measures it well: a remnant a cell and a bit thick meshes, and its
        volume is still rough.  The warning names the cell so that a caller
        who needs better can ask for better.  A region the Lipschitz bound
        proves holds nothing is not refined: it meshes to an empty solid,
        which is what an empty CSG result is.
    :param refine_limit: the most cells across the longest axis refinement
        will ask for.
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
    :func:`yapcad.sdf.node.to_construction`.  What is recorded is the
    resolution the mesh was *made* at, refinement included, which is also
    what later booleans read to keep from meshing the part coarser than it
    already is.
    """
    if brep not in (False, True, "auto"):
        raise SdfError(f"brep must be True, False or 'auto', got {brep!r}")
    replayed = _replay(node, brep)

    bodies = _bodies(node)
    region = node.bounds if bounds is None else bounds
    start = _cells_across(region, resolution, cell_mm)
    limit = max(int(refine_limit), start)

    def plan_for(cells):
        """What to mesh, where, and at what resolution, for ``cells``
        across the whole region.

        A compound's bodies are meshed one at a time at the cell size the
        whole region would have had, so the result neither fuses touching
        bodies nor depends on how they were grouped.
        """
        if len(bodies) == 1:
            return [(node, bounds, cells)]
        cell = _longest(region) / float(cells)
        plan = []
        for body in bodies:
            body_bounds = body.bounds if bounds is None else \
                bounds_intersect(body.bounds, bounds)
            if bounds_is_empty(body_bounds):
                continue
            plan.append((
                body,
                None if bounds is None else body_bounds,
                max(2, math.ceil(_longest(body_bounds) / cell)),
            ))
        return plan

    def attempt(cells):
        """Mesh everything at ``cells``.

        Returns the surfaces, and the bodies whose lattice came back without
        a single triangle.  Those are not an answer on their own: empty
        space and a part thinner than a cell look the same here, which is
        what :func:`~yapcad.sdf.contour.proves_empty` is for.
        """
        surfaces, absent = [], []
        for body, body_bounds, body_cells in plan_for(cells):
            mesh = dual_contour(body, resolution=body_cells,
                                bounds=body_bounds, padding=padding,
                                chunk=chunk, prune=prune, backend=backend)
            if len(mesh.triangles) == 0:
                absent.append((body, body_bounds, body_cells))
                continue
            if check:
                defects = manifold_defects(mesh)
                if any(defects.values()):
                    raise NonManifoldMeshError(
                        f"dual contouring produced an unusable mesh at "
                        f"resolution {body_cells}: "
                        f"{defects['nonmanifold']} non-manifold edges, "
                        f"{defects['boundary']} boundary edges, "
                        f"{defects['coincident']} coincident vertices. The "
                        f"surface most likely passes through one cell twice "
                        f"-- a wall or a gap thinner than a cell. Raise the "
                        f"resolution, or pass check=False to accept the "
                        f"mesh as-is."
                    )
            surfaces.append(mesh_to_surface(mesh))
        return surfaces, absent

    cells = start
    reason = None
    while True:
        try:
            surfaces, absent = attempt(cells)
        except NonManifoldMeshError as exc:
            failure = exc
        else:
            failure = None
            missed = [item for item in absent
                      if not proves_empty(item[0], resolution=item[2],
                                          bounds=item[1], padding=padding,
                                          chunk=chunk, backend=backend)]
            if missed:
                failure = missed
        if failure is None:
            break
        if reason is None:
            reason = failure
        nxt = min(2 * cells, limit)
        if not refine or nxt == cells:
            raise _exhausted(failure, region, start, cells, refine)
        cells = nxt

    if cells != start:
        warnings.warn(_refined(reason, region, start, cells, surfaces),
                      RuntimeWarning, stacklevel=2)
    meshing = {
        "method": "dual-contouring",
        "resolution": int(cells),
    }
    # Only record what was actually specified; a null in the document would
    # have to round-trip as a null, and an absent key says the same thing.
    if padding is not None:
        meshing["padding"] = float(padding)
    if bounds is not None:
        meshing["bounds"] = [float(v) for v in bounds[0]] + \
                            [float(v) for v in bounds[1]]
    record = to_construction(node, meshing=meshing)
    result = solid(surfaces, [], record)
    if metadata:
        from yapcad.metadata import set_solid_metadata
        set_solid_metadata(result, dict(metadata))
    if replayed is not None:
        _attach_derived_brep(result, replayed, node)
    return result


def _cells_across(region, resolution, cell_mm):
    """The lattice the caller asked for, as cells across the longest axis:
    ``resolution`` itself, or what ``cell_mm`` works out to."""
    if cell_mm is None:
        return max(2, int(resolution))
    cell = float(cell_mm)
    if not math.isfinite(cell) or cell <= 0.0:
        raise SdfError("to_solid: cell_mm must be positive and finite")
    longest = _longest(region)
    if not math.isfinite(longest) or longest <= 0.0:
        raise SdfError(
            "to_solid: cell_mm needs a region of finite extent to measure "
            "against; pass bounds= for an unbounded field"
        )
    return max(2, math.ceil(longest / cell))


def _cell_of(region, cells):
    """The cell size a lattice of ``cells`` across ``region`` works out to."""
    return _longest(region) / float(cells)


def _span(region, start, final):
    """How the lattice was tried, for a message."""
    if start == final:
        return f"cell {_cell_of(region, start):.3g} mm (resolution {start})"
    return (f"cell {_cell_of(region, start):.3g} mm down to "
            f"{_cell_of(region, final):.3g} mm (resolution {start} to "
            f"{final})")


def _refined(reason, region, start, final, surfaces):
    """The warning for a mesh that only came out after refinement."""
    opening = ("could not close the mesh"
               if isinstance(reason, NonManifoldMeshError)
               else "found no surface")
    vertices = sum(len(surf[1]) for surf in surfaces)
    cell = _cell_of(region, final)
    return (
        f"to_solid: dual contouring {opening} at "
        f"cell {_cell_of(region, start):.3g} mm (resolution {start}), so it "
        f"refined to cell {cell:.3g} mm (resolution {final}), where the "
        f"mesh came out closed with {vertices} vertices. Pass "
        f"cell_mm={cell:.3g} to mesh it in one pass."
    )


def _exhausted(failure, region, start, final, refine):
    """The refusal when even the finest lattice tried found nothing."""
    tried = _span(region, start, final)
    if not refine and start == final:
        tried += ", and refine=False"
    if isinstance(failure, NonManifoldMeshError):
        return NonManifoldMeshError(f"{failure} Tried at {tried}.")
    where = ""
    if len(failure) == 1 and failure[0][1] is not None:
        lo, hi = failure[0][1]
        where = (" The region in question spans "
                 f"{tuple(round(v, 4) for v in lo)} to "
                 f"{tuple(round(v, 4) for v in hi)}.")
    elif len(failure) > 1:
        where = f" {len(failure)} bodies came back empty."
    return SdfError(
        f"to_solid: dual contouring found no surface at {tried}, and the "
        f"Lipschitz bound could not prove the region holds nothing either, "
        f"so this is either a part thinner than the finest cell tried or an "
        f"empty result that cannot be told apart from one. Mesh it finer "
        f"with cell_mm=, or raise refine_limit above {final}." + where
    )


def _longest(region):
    lo, hi = region
    return max(hi[i] - lo[i] for i in range(3))


def _bodies(node):
    """The separate bodies of ``node``: its compound's members, with any
    transforms above the compound pushed down onto them, or ``[node]``."""
    if node.kind == "compound":
        return [b for child in node.children for b in _bodies(child)]
    if node.kind == "transform":
        inner = _bodies(node.children[0])
        if len(inner) > 1:
            return [make_node("transform", node.p, (b,)) for b in inner]
    return [node]


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


__all__ = ["MAX_RESOLUTION", "NonManifoldMeshError",
           "mesh_to_surface", "to_solid"]
