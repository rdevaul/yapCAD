"""Building the mesh a deferred solid has been carrying a plan for.

:func:`yapcad.sdf.convert.deferring` lets a design build without meshing
anything, because most of what a design does with a part never needs
triangles.  This module is the other half: the test for whether a solid is
still a plan, the function that turns it into a mesh, and the refusal that
everything needing triangles raises instead of reading an empty surface list
and believing it.

That refusal is the point.  An empty surface list is already legal in yapCAD
-- it is what an empty CSG result looks like -- and the code that consumes
solids was written when the only way to get one was to mean it.  So
``issolidclosed`` calls an empty solid closed, ``volumeof`` calls it zero and
``solid_contains_point`` puts every point outside: perfectly reasonable
answers about nothing at all, and all wrong about a part whose mesh has not
been built yet.  :func:`require_mesh` is how those sites tell the two apart.
"""

import warnings

from yapcad.construction import get_construction, set_construction
from yapcad.sdf.convert import DEFERRED, to_solid
from yapcad.sdf.node import (
    SdfError,
    from_construction,
    is_sdf_construction,
    meshing_from_construction,
)

__all__ = ["DeferredMeshError", "is_deferred", "deferred_plan",
           "plan_from_sources", "realize", "realize_all", "require_mesh",
           "realized"]


class DeferredMeshError(SdfError):
    """Raised when an operation needs triangles a solid has not built yet."""


def _record(target):
    try:
        record = get_construction(target)
    except Exception:
        return None
    return record if is_sdf_construction(record) else None


def is_deferred(target):
    """True if ``target`` is a solid carrying a field and no mesh.

    False for a field-authored solid that has its mesh, and false for a
    field that was *meshed* to nothing -- an empty CSG result is a real
    answer rather than a postponed one, and it says so by carrying the
    meshing record that produced it.

    A field-authored solid with no surfaces and no meshing record at all is
    read as deferred.  That is the backstop: every site that composes a new
    tree has to carry the plan with it, and one that forgets would otherwise
    produce exactly the silent empty solid this module exists to prevent.
    Being wrong this way costs a loud refusal; being wrong the other way
    costs a wrong answer.
    """
    record = _record(target)
    if record is None:
        return False
    meshing = meshing_from_construction(record)
    if meshing is None:
        return not target[1]
    return meshing.get("method") == DEFERRED


def deferred_plan(target):
    """The meshing parameters :func:`realize` would use, or ``None``."""
    if not is_deferred(target):
        return None
    return dict(meshing_from_construction(_record(target)) or {})


def plan_from_sources(tree, sources):
    """The meshing plan a composed field inherits from the solids it came
    from, or ``None`` if none of them were deferred.

    Resolutions are per-solid -- cells across *that* solid's longest axis --
    so they cannot simply be carried over to a tree that spans more space.
    They are converted to a cell size, the finest is kept, and that is
    turned back into a resolution for the new tree, which is how
    :mod:`yapcad.sdf.booleans` sizes a boolean against its operands.
    """
    cells = []
    for source in sources:
        plan = deferred_plan(source)
        if not plan or not plan.get("resolution"):
            continue
        node = from_construction(_record(source))
        if node is None:
            continue
        span = _longest(node.bounds)
        if span is not None and span > 0.0:
            cells.append(span / float(plan["resolution"]))
    if not cells:
        return None
    span = _longest(tree.bounds)
    if span is None or span <= 0.0:
        return {"method": DEFERRED, "resolution": 64}
    import math
    return {"method": DEFERRED,
            "resolution": max(2, int(math.ceil(span / min(cells))))}


def _longest(bounds):
    import math
    lo, hi = bounds
    spans = [hi[i] - lo[i] for i in range(3)]
    if not all(math.isfinite(v) for v in spans):
        return None
    return max(spans)


def realize(target, simplify=False, **overrides):
    """Mesh a deferred solid, in place, and return it.

    In place because an assembly, a compound and the caller all hold the
    same list: handing back a new solid would leave every other reference
    still holding the plan.  A solid that already has its mesh is returned
    untouched, so this is safe to call defensively.

    :param simplify: simplify the new mesh against its field.  Off by
        default, because this builds exactly the mesh the deferred
        :func:`~yapcad.sdf.to_solid` would have built, and yapCAD
        simplifies *finished* parts rather than intermediate ones --
        simplifying costs about twenty times the meshing.  Use
        :func:`realize_all`, which follows ``YAPCAD_SDF_SIMPLIFY``, to
        finish a whole design the way the eager path would.
    :param overrides: passed to :func:`yapcad.sdf.to_solid` over the
        recorded plan -- ``resolution``, ``cell_mm``, ``brep`` and the rest.
    """
    if not is_deferred(target):
        return target
    record = _record(target)
    node = from_construction(record)
    if node is None:
        raise DeferredMeshError(
            "this solid is marked as deferred but carries no field to mesh")
    plan = dict(meshing_from_construction(record) or {})
    kwargs = {"resolution": int(plan.get("resolution", 64))}
    if "padding" in plan:
        kwargs["padding"] = float(plan["padding"])
    if "bounds" in plan:
        flat = plan["bounds"]
        kwargs["bounds"] = (tuple(flat[:3]), tuple(flat[3:]))
    kwargs.update(overrides)
    kwargs["defer"] = False

    built = to_solid(node, **kwargs)
    from yapcad.sdf.simplify import simplify_finished
    built = simplify_finished(built, simplify)
    target[1] = built[1]
    set_construction(target, get_construction(built))
    return target


def realize_all(value, simplify=None, **overrides):
    """Mesh every deferred solid reachable from ``value``, in place.

    Handles a solid, a nested list or tuple of them, a mapping of them, and
    an object carrying a ``geometry`` mapping -- which is what an
    :class:`~yapcad.assembly.assembly.Assembly` is.  Each distinct solid is
    meshed once however many times it appears.

    This is the "finish the design" call, so unlike :func:`realize` it
    simplifies by default, which is what the eager path does to a finished
    part.  ``simplify=False`` turns that off.

    :returns: the number of solids meshed.
    """
    if simplify is None:
        from yapcad.sdf.simplify import _enabled
        simplify = _enabled()
    seen = set()

    def walk(item):
        count = 0
        if isinstance(item, list) and item and item[0] == "solid":
            if id(item) not in seen:
                seen.add(id(item))
                if is_deferred(item):
                    realize(item, simplify=simplify, **overrides)
                    count += 1
            return count
        if isinstance(item, (list, tuple)):
            for child in item:
                count += walk(child)
            return count
        if hasattr(item, "items"):
            for child in item.values():
                count += walk(child)
            return count
        geometry = getattr(item, "geometry", None)
        if geometry is not None and hasattr(geometry, "items"):
            count += walk(geometry)
        return count

    return walk(value)


def require_mesh(target, operation, hint=None):
    """Refuse ``operation`` on a solid whose mesh has not been built.

    Call this wherever an empty surface list would otherwise be taken at
    face value.  A no-op for every solid that is not deferred, including a
    genuinely empty one.
    """
    if not is_deferred(target):
        return target
    raise DeferredMeshError(
        f"{operation} needs triangles, and this solid carries a field whose "
        f"mesh is deferred. Build it with yapcad.sdf.realize(solid)"
        + (f" {hint}" if hint else "")
        + ", or yapcad.sdf.realize_all(result) for a whole design, first."
    )


def realized(target, operation, **overrides):
    """Mesh ``target`` if it is deferred, warning that nobody asked.

    For the places where refusing would be worse than paying -- a viewer
    that would otherwise show nothing.  Prefer :func:`require_mesh` wherever
    the caller can reasonably be told to realize first.
    """
    if not is_deferred(target):
        return target
    plan = deferred_plan(target) or {}
    warnings.warn(
        f"{operation} needed triangles, so this solid's deferred mesh was "
        f"built at resolution {plan.get('resolution')}. Call "
        f"yapcad.sdf.realize(solid) where you want to pay for that.",
        RuntimeWarning, stacklevel=3)
    return realize(target, **overrides)
