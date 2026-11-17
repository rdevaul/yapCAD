"""The watertightness check shared by the CLIs and the mesh tools.

One place decides what "printable" means: triangle count, whether the
mesh is watertight, how many separate bodies it holds, and how big it
is.  The deferred-meshing STL path in ``python -m yapcad.dsl run``, the
tooling in ``tools/`` and any scripted check all read the same
:class:`MeshCheck`, so their verdicts cannot drift apart.

The mesh argument may be a trimesh mesh, a path to a mesh file, or a
yapCAD solid.  Reading a mesh file needs trimesh.  A yapCAD solid does
not: without trimesh its triangles are welded and checked with numpy,
so a plain ``pip install yapcad`` can still answer.  The component count
never depends on trimesh's optional graph engines (scipy or networkx);
when they are missing it is counted here.  ``None`` in a report means
unverified, and callers that require watertightness must not read it as
passed.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

__all__ = ["MeshCheck", "check_mesh", "trimesh_available"]


@dataclass(frozen=True)
class MeshCheck:
    """What came back from checking one mesh."""

    triangles: int
    watertight: Optional[bool]
    components: Optional[int]
    bounds: Optional[Tuple[Tuple[float, float, float],
                           Tuple[float, float, float]]]

    def summary(self) -> str:
        watertight = ("yes" if self.watertight else
                      "no" if self.watertight is not None else
                      "unknown (trimesh not installed)")
        components = (str(self.components) if self.components is not None
                      else "?")
        bounds = ("-" if self.bounds is None else
                  f"({self.bounds[0][0]:.3f}, {self.bounds[0][1]:.3f}, "
                  f"{self.bounds[0][2]:.3f}) to "
                  f"({self.bounds[1][0]:.3f}, {self.bounds[1][1]:.3f}, "
                  f"{self.bounds[1][2]:.3f})")
        return (f"mesh: {self.triangles} triangles, "
                f"watertight: {watertight}, "
                f"components: {components}, bbox: {bounds}")


def trimesh_available() -> bool:
    try:
        import trimesh  # noqa: F401
    except ImportError:
        return False
    return True


def _as_trimesh(target):
    """Adapt a trimesh mesh, a mesh file path, or a yapCAD solid."""
    import trimesh

    if isinstance(target, (str, bytes)) or hasattr(target, "__fspath__"):
        return trimesh.load(target, force="mesh")
    if hasattr(target, "is_watertight") and hasattr(target, "faces"):
        return target
    # A yapCAD solid: reuse the mesh engine's converter, which dedupes
    # vertices and repairs the winding so the checks below see one
    # closed manifold rather than per-facet noise.
    from yapcad.boolean.trimesh_engine import _solid_to_mesh
    return _solid_to_mesh(target)


def _welded(target):
    """``(vertices, faces)`` arrays of a yapCAD solid, vertices welded.

    ``mesh_view`` hands back loose triangles, so corners that coincide
    are merged here.  The weld tolerance is relative to the mesh's own
    size, so the verdict does not depend on the model's units.
    """
    import numpy as np
    from yapcad.mesh import mesh_view

    corners = np.array([[v0, v1, v2] for _, v0, v1, v2 in mesh_view(target)],
                       dtype=float).reshape(-1, 3)
    if not len(corners):
        return corners, np.zeros((0, 3), dtype=np.int64)
    span = float(np.ptp(corners, axis=0).max()) or 1.0
    keys = np.round(corners / (span * 1e-9)).astype(np.int64)
    _, first, inverse = np.unique(keys, axis=0, return_index=True,
                                  return_inverse=True)
    faces = inverse.reshape(-1, 3)
    # A triangle whose corners weld together has no area and no edges.
    distinct = ((faces[:, 0] != faces[:, 1]) & (faces[:, 1] != faces[:, 2])
                & (faces[:, 0] != faces[:, 2]))
    return corners[first], faces[distinct]


def _edge_groups(faces):
    """Each face edge's group id, with the groups' sizes.

    Edges are grouped without regard to direction; row ``3 * f + k`` is
    edge ``k`` of face ``f``.
    """
    import numpy as np

    edges = np.stack([faces, np.roll(faces, -1, axis=1)],
                     axis=2).reshape(-1, 2)
    _, group, counts = np.unique(np.sort(edges, axis=1), axis=0,
                                 return_inverse=True, return_counts=True)
    return edges, group.reshape(-1), counts


def _is_watertight(faces) -> bool:
    """True when every edge is shared by exactly two faces, which is the
    test trimesh's ``is_watertight`` makes."""
    import numpy as np

    if not len(faces):
        return False
    _, _, counts = _edge_groups(faces)
    return bool(np.all(counts == 2))


def _count_components(faces) -> int:
    """Bodies in the mesh: sets of faces connected through shared edges."""
    import numpy as np

    if not len(faces):
        return 0
    _, group, _ = _edge_groups(faces)
    face_of = np.repeat(np.arange(len(faces)), 3)
    label = np.arange(len(faces))
    while True:
        # Every face takes the smallest label on any of its edges, and
        # then follows that label to its own, which halves the distance
        # left to travel each pass.
        lowest = np.full(int(group.max()) + 1, len(faces), dtype=np.int64)
        np.minimum.at(lowest, group, label[face_of])
        new = label.copy()
        np.minimum.at(new, face_of, lowest[group])
        new = new[new]
        if np.array_equal(new, label):
            break
        label = new
    return int(len(np.unique(label)))


def _check_solid(target) -> MeshCheck:
    """The check for a yapCAD solid when trimesh is not installed."""
    vertices, faces = _welded(target)
    bounds = None
    if len(vertices):
        lo, hi = vertices.min(axis=0), vertices.max(axis=0)
        bounds = (tuple(float(v) for v in lo), tuple(float(v) for v in hi))
    return MeshCheck(triangles=int(len(faces)),
                     watertight=_is_watertight(faces),
                     components=_count_components(faces), bounds=bounds)


def check_mesh(target) -> MeshCheck:
    """Check a mesh for the things STL fabrication cares about.

    :param target: a trimesh mesh, a path to a mesh file, or a yapCAD
        solid.  Deferred solids are refused by the conversion layer with
        a named reason; this check never meshes anything itself.
    """
    try:
        import trimesh  # noqa: F401
    except ImportError:
        if isinstance(target, (str, bytes)) or hasattr(target, "__fspath__"):
            raise ImportError(
                "trimesh is required to read mesh files; install "
                "trimesh[easy] or pass a yapCAD solid") from None
        return _check_solid(target)

    mesh = _as_trimesh(target)
    triangles = int(len(mesh.faces))
    watertight = bool(mesh.is_watertight)

    try:
        components = len(mesh.split(only_watertight=False))
    except Exception:
        # trimesh splits with scipy or networkx and has neither here.
        import numpy as np
        components = _count_components(np.asarray(mesh.faces))

    bounds = None
    if len(mesh.vertices):
        lo, hi = mesh.bounds
        bounds = (tuple(float(v) for v in lo), tuple(float(v) for v in hi))

    return MeshCheck(triangles=triangles, watertight=watertight,
                     components=components, bounds=bounds)
