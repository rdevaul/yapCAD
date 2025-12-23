"""Utilities for working with triangulated views of yapCAD surfaces."""

from __future__ import annotations

from typing import Iterable, Iterator, Sequence, Tuple

from yapcad.geom3d import issolid, issurface
from yapcad.geometry_utils import to_vec3, triangle_normal

Vec3 = Tuple[float, float, float]
TriTuple = Tuple[Vec3, Vec3, Vec3, Vec3]


def mesh_view(obj: Sequence) -> Iterator[TriTuple]:
    """Yield triangles for a surface or solid as ``(normal, v0, v1, v2)``.

    Normals are unit vectors. Vertices are returned as ``(x, y, z)`` tuples.
    Faces with degenerate geometry (zero area) are skipped silently.

    For BREP-only solids (those with empty surfaces list), tessellation is
    performed via pythonocc-core if available.
    """

    surfaces: Iterable[Sequence]
    if issurface(obj):
        surfaces = [obj]
    elif issolid(obj):
        surfaces = obj[1]
        # Handle BREP-only solids with empty surfaces list
        if not surfaces:
            yield from _mesh_view_brep(obj)
            return
    else:
        raise ValueError("mesh_view expects a surface or solid")

    for surf in surfaces:
        verts = surf[1]
        faces = surf[3] if len(surf) > 3 else []

        for face in faces:
            if len(face) != 3:
                continue
            idx0, idx1, idx2 = face
            v0 = to_vec3(verts[idx0])
            v1 = to_vec3(verts[idx1])
            v2 = to_vec3(verts[idx2])

            calc_normal = triangle_normal(v0, v1, v2)
            if calc_normal is None:
                continue

            yield calc_normal, v0, v1, v2


def _mesh_view_brep(solid: Sequence) -> Iterator[TriTuple]:
    """Yield triangles from a BREP-only solid via tessellation."""
    try:
        from yapcad.brep import brep_from_solid, occ_available
        if not occ_available():
            return
        brep = brep_from_solid(solid)
        if brep is None:
            return
        # Tessellate the BREP shape
        surface = brep.tessellate(deflection=0.5)
        # surface format: ['surface', vertices, normals, triangles, [], []]
        verts = surface[1]
        faces = surface[3]
        for face in faces:
            if len(face) != 3:
                continue
            idx0, idx1, idx2 = face
            v0 = to_vec3(verts[idx0])
            v1 = to_vec3(verts[idx1])
            v2 = to_vec3(verts[idx2])
            calc_normal = triangle_normal(v0, v1, v2)
            if calc_normal is None:
                continue
            yield calc_normal, v0, v1, v2
    except ImportError:
        return

