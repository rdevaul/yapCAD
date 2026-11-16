"""The watertightness check shared by the CLIs and the mesh tools.

One place decides what "printable" means: triangle count, whether the
mesh is watertight, how many separate bodies it holds, and how big it
is.  The deferred-meshing STL path in ``python -m yapcad.dsl run``, the
tooling in ``tools/`` and any scripted check all read the same
:class:`MeshCheck`, so their verdicts cannot drift apart.

The mesh argument may be a trimesh mesh, a path to a mesh file, or a
yapCAD solid.  watertight and components are computed with trimesh when
it is installed; without it the numbers are reported as ``None`` and the
triangle count (which needs no trimesh) still is.  Callers that require
watertightness must treat ``None`` as unverified, not as passed.
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


def _count_triangles(target) -> int:
    """Triangle count without trimesh, for the degraded report."""
    from yapcad.geometry_utils import triangles_from_mesh
    from yapcad.mesh import mesh_view

    return len(triangles_from_mesh(mesh_view(target)))


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
        return MeshCheck(triangles=_count_triangles(target),
                         watertight=None, components=None, bounds=None)

    mesh = _as_trimesh(target)
    triangles = int(len(mesh.faces))
    watertight = bool(mesh.is_watertight)

    components = None
    try:
        components = len(mesh.split(only_watertight=False))
    except Exception:
        try:
            import trimesh.graph
            labels = trimesh.graph.connected_component_labels(mesh.edges)
            components = int(labels.max(initial=0) + 1) if len(labels) \
                else 0
        except Exception:
            components = None

    bounds = None
    if len(mesh.vertices):
        lo, hi = mesh.bounds
        bounds = (tuple(float(v) for v in lo), tuple(float(v) for v in hi))

    return MeshCheck(triangles=triangles, watertight=watertight,
                     components=components, bounds=bounds)
