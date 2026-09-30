"""Keep an SDF solid's authoritative tree in step with transforms of it.

The solid transforms in :mod:`yapcad.geom3d` move a solid's surfaces and
then give each representation a chance to follow -- ``translate_brep_solid``,
``translate_native_brep`` and so on.  Without an equivalent for the SDF tree,
translating an SDF-authoritative solid moved its mesh and left the tree that
*defines* it where it was: the document then described a sphere at the
origin carrying a mesh a hundred units away.  This module is that missing
hook.

Three things happen to a transformed SDF solid:

* the tree is wrapped in a ``transform`` node for the same matrix, so the
  authoritative definition follows the geometry;
* the recorded meshing parameters are dropped, because a mesh transformed
  after meshing is not what those parameters would regenerate, and the
  document should not claim a reproducibility it cannot honour;
* a derived BREP is re-tagged against the new tree when the BREP hook
  transformed it too -- replaying ``transform(tree, M)`` is literally the
  replayed shape carried through ``M``, so the tag stays truthful -- and is
  dropped when it did not (a non-uniform scale leaves the BREP behind).  A
  derived BREP carries nothing the tree does not, so dropping one is safe.
"""

import numpy as np

from yapcad.construction import get_construction, set_construction
from yapcad.sdf.node import from_construction, is_sdf_construction, \
    to_construction
from yapcad.sdf.ops import transform

__all__ = ["follow_transform", "is_similarity"]


def is_similarity(matrix, tolerance=1e-9):
    """True if the linear part of ``matrix`` is a scaled orthogonal map."""
    linear = np.asarray([row[:3] for row in matrix[:3]], dtype=float)
    sv = np.linalg.svd(linear, compute_uv=False)
    return float(sv.max() - sv.min()) <= tolerance * max(1.0, float(sv.max()))


def follow_transform(original, transformed, matrix):
    """Carry ``transformed``'s SDF tree through ``matrix``.

    :param original: the solid before the transform, read for its tree and
        for whether its BREP was derived from that tree.
    :param transformed: the transformed copy, updated in place.
    :param matrix: the row-major affine 4x4 the surfaces were moved by.

    A solid without an SDF construction record is left alone.
    """
    record = get_construction(original)
    if not is_sdf_construction(record):
        return transformed
    before = from_construction(record)
    after = transform(before, matrix)
    set_construction(transformed, to_construction(after))
    _follow_derived_brep(original, transformed, before, after, matrix)
    return transformed


def _follow_derived_brep(original, transformed, before, after, matrix):
    from yapcad.metadata import get_solid_metadata

    old_meta = get_solid_metadata(original) or {}
    old_brep = old_meta.get("brep") or {}
    if (old_brep.get("derivedFrom") or {}).get("sdf") != before.digest:
        return  # not derived from this tree; not ours to manage

    meta = get_solid_metadata(transformed)
    brep = meta.get("brep")
    if brep is None:
        return
    if is_similarity(matrix):
        brep["derivedFrom"] = {"sdf": after.digest}
        return
    # The BREP hooks only follow similarity transforms, so this BREP was
    # left behind.  geom3d.scalesolid already drops a BREP on a non-uniform
    # scale; this is the safety net for any caller that does not.  It is
    # derived, so dropping it loses nothing the tree does not define.
    from yapcad.brep import _clear_brep_data
    _clear_brep_data(transformed)
