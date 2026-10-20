"""Retain native DSL booleans until a consumer requests mutable mesh surfaces.

No global geometry cache. Reading surfaces retires the native copy so edits
through the public yapCAD list representation cannot leave a stale cache.
"""

from __future__ import annotations
import copy


class ManifoldSolid(list):
    def __init__(self, body, operation):
        super().__init__(
            ["solid", [], [], ["boolean", f"manifold:{operation}"], {}]
        )
        self.body = body

    def materialize(self):
        if self.body is not None:
            from yapcad.boolean.manifold_engine import _from_manifold

            mesh = _from_manifold(self.body, "retained")
            list.__setitem__(self, 1, mesh[1])
            self.body = None
        return self

    def __getitem__(self, key):
        if isinstance(key, slice) or key in (1, -4):
            self.materialize()
        return super().__getitem__(key)

    def __setitem__(self, key, value):
        if isinstance(key, slice) or key in (0, 1, -5, -4):
            self.materialize()
        return super().__setitem__(key, value)

    def __repr__(self):
        self.materialize()
        return super().__repr__()

    def __ne__(self, other):
        return not self == other

    def __delitem__(self, key):
        self.materialize()
        return super().__delitem__(key)

    def pop(self, index=-1):
        self.materialize()
        return super().pop(index)

    def insert(self, index, value):
        self.materialize()
        return super().insert(index, value)

    def append(self, value):
        self.materialize()
        return super().append(value)

    def extend(self, values):
        self.materialize()
        return super().extend(values)

    def clear(self):
        self.materialize()
        return super().clear()

    def remove(self, value):
        self.materialize()
        return super().remove(value)

    def reverse(self):
        self.materialize()
        return super().reverse()

    def sort(self, *args, **kwargs):
        self.materialize()
        return super().sort(*args, **kwargs)

    def __iadd__(self, values):
        self.materialize()
        return super().__iadd__(values)

    def __imul__(self, count):
        self.materialize()
        return super().__imul__(count)

    def copy(self):
        self.materialize()
        return list(self)

    def __eq__(self, other):
        self.materialize()
        if isinstance(other, ManifoldSolid):
            other.materialize()
        return super().__eq__(other)

    def __iter__(self):
        self.materialize()
        return super().__iter__()

    def __deepcopy__(self, memo):
        if self.body is None:
            return copy.deepcopy(list(self), memo)
        out = ManifoldSolid(self.body, "copy")
        memo[id(self)] = out
        for i in (2, 3, 4):
            list.__setitem__(
                out, i, copy.deepcopy(list.__getitem__(self, i), memo)
            )
        return out


def native(solid):
    from yapcad.boolean.manifold_engine import _to_manifold

    if isinstance(solid, ManifoldSolid) and solid.body is not None:
        return solid.body
    return _to_manifold(solid, "DSL")


def eligible(operands):
    """True when the batched manifold3d path is the right thing to run.

    The engine decision belongs to yapcad.boolean.selection rather than
    being repeated here: this path is an optimisation of whatever
    solid_boolean would have done, so it must not run when solid_boolean
    would have chosen a different engine.
    """
    from yapcad.boolean.selection import forced_engine, preferred_mesh_engine
    from yapcad.brep import has_brep_data
    from yapcad.sdf.booleans import is_sdf_solid

    # An outright override bypasses the SDF and OCC paths in solid_boolean,
    # so it bypasses the deference to them here too.
    forced = forced_engine()
    if forced:
        return forced == "manifold"

    engine, _named = preferred_mesh_engine()
    if engine != "manifold":
        return False
    # Otherwise solid_boolean would hand these operands to the field or OCC
    # engine before any mesh engine saw them.
    return not any(is_sdf_solid(s) or has_brep_data(s) for s in operands)


def try_boolean(operands, operation):
    """Return a retained mesh result, or None for upstream engine dispatch."""
    from yapcad.boolean import manifold_engine

    if (
        not operands
        or not manifold_engine.is_available()
        or not eligible(operands)
    ):
        return None
    from manifold3d import Manifold, OpType

    try:
        bodies = [native(s) for s in operands]
    except manifold_engine.NotManifoldError:
        # Preserve the existing validation/fallback policy for invalid inputs.
        return None
    op = {
        "union": OpType.Add,
        "difference": OpType.Subtract,
        "intersection": OpType.Intersect,
    }[operation]
    return ManifoldSolid(Manifold.batch_boolean(bodies, op), operation)


def try_transform(solid, operation, xyz):
    if not isinstance(solid, ManifoldSolid) or solid.body is None:
        return None
    return ManifoldSolid(getattr(solid.body, operation)(xyz), operation)


def materialize(value):
    """Do not leak deferred native objects across the public execution API."""
    if isinstance(value, ManifoldSolid):
        return list(value)
    if isinstance(value, list):
        # Ordinary solids are already materialized; only traverse DSL lists.
        if value and value[0] == "solid":
            return value
        return [materialize(item) for item in value]
    return value
