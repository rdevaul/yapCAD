"""Trimesh-backed boolean engine for yapCAD solids.

This engine is optional. It converts yapCAD solids to ``trimesh.Trimesh``
instances, dispatches boolean operations via :mod:`trimesh.boolean`, and
converts the resulting mesh back into a yapCAD solid.

Availability depends on both the ``trimesh`` package and at least one
boolean backend supported by ``trimesh`` (e.g. Blender, OpenSCAD, Cork).
"""

from __future__ import annotations

from typing import Iterable

import os
import numpy as np

try:
    import trimesh
except ImportError:  # pragma: no cover - optional dependency
    trimesh = None  # type: ignore[assignment]

from . import native as _native

ENGINE_NAME = "trimesh"


def engines_available() -> set[str]:
    """Return the set of trimesh boolean backends that are operational."""

    if trimesh is None:  # pragma: no cover - optional dependency
        return set()
    return set(trimesh.boolean.engines_available)


def is_available(backend: str | None = None) -> bool:
    """Check whether the engine can run (trimesh + backend present)."""

    available = engines_available()
    if not available:
        return False
    if backend is None:
        return True
    return backend in available


def _solid_to_mesh(sld) -> "trimesh.Trimesh":
    if trimesh is None:  # pragma: no cover - optional dependency
        raise RuntimeError("trimesh is not installed")

    triangles = list(_native._iter_triangles_from_solid(sld))
    if not triangles:
        return trimesh.Trimesh(vertices=np.zeros((0, 3)), faces=np.zeros((0, 3), dtype=np.int64), process=False)

    vertex_map: dict[tuple[float, float, float], int] = {}
    verts = []
    faces = []
    for tri in triangles:
        face_inds = []
        for pt in tri:
            key = (round(pt[0], 9), round(pt[1], 9), round(pt[2], 9))
            idx = vertex_map.get(key)
            if idx is None:
                idx = len(verts)
                vertex_map[key] = idx
                verts.append([pt[0], pt[1], pt[2]])
            face_inds.append(idx)
        faces.append(face_inds)

    verts_np = np.asarray(verts, dtype=float)
    faces_np = np.asarray(faces, dtype=np.int64)
    mesh = trimesh.Trimesh(vertices=verts_np, faces=faces_np, process=False)
    mesh.remove_unreferenced_vertices()
    if not mesh.is_watertight:
        mesh.merge_vertices()
        mesh.remove_duplicate_faces()
        mesh.remove_degenerate_faces()
        mesh.process(validate=True)
    try:
        mesh.fix_normals()
    except Exception:
        pass
    return mesh


def _mesh_to_solid(mesh: "trimesh.Trimesh", operation: str) -> list:
    triangles = np.asarray(mesh.triangles)
    if triangles.size == 0:
        return _native._geom3d().solid([], [], ['boolean', f'{ENGINE_NAME}:{operation}'])

    tri_points = []
    for tri in triangles:
        tri_points.append([
            _native.point(tri[0][0], tri[0][1], tri[0][2]),
            _native.point(tri[1][0], tri[1][1], tri[1][2]),
            _native.point(tri[2][0], tri[2][1], tri[2][2]),
        ])

    surface = _native._surface_from_triangles(tri_points)
    if surface is None:
        return _native._geom3d().solid([], [], ['boolean', f'{ENGINE_NAME}:{operation}'])
    return _native._geom3d().solid([surface], [], ['boolean', f'{ENGINE_NAME}:{operation}'])


def solid_boolean(a, b, operation: str, tol=_native._DEFAULT_RAY_TOL, *, stitch: bool = False,
                  backend: str | None = None):
    """Perform a boolean between ``a`` and ``b`` using trimesh."""

    if trimesh is None:  # pragma: no cover - optional dependency
        raise RuntimeError("trimesh is not installed; install trimesh to enable this engine")

    available = engines_available()
    if backend is not None and backend not in available:
        raise RuntimeError(
            f"trimesh backend '{backend}' is not available; install the appropriate binary (available: {available})"
        )
    if backend is None and not available:
        raise RuntimeError(
            "no trimesh boolean backends are available; install Blender, OpenSCAD, Cork, or another supported engine"
        )

    mesh_a = _solid_to_mesh(a)
    mesh_b = _solid_to_mesh(b)

    backup_cache = None
    if backend == 'blender':
        backup_cache = os.environ.get('ARCH_CACHE_LINE_SIZE')
        os.environ['ARCH_CACHE_LINE_SIZE'] = '64'

    op = operation.lower()
    try:
        if op == 'union':
            result = trimesh.boolean.union([mesh_a, mesh_b], engine=backend, check_volume=False)
        elif op == 'intersection':
            result = trimesh.boolean.intersection([mesh_a, mesh_b], engine=backend, check_volume=False)
        elif op == 'difference':
            result = trimesh.boolean.difference([mesh_a, mesh_b], engine=backend, check_volume=False)
        else:
            raise RuntimeError(f"unsupported boolean operation '{operation}' for trimesh engine")
    except Exception as exc:  # pragma: no cover - depends on external binaries
        raise RuntimeError(f"trimesh boolean operation failed: {exc}") from exc
    finally:
        if backend == 'blender':
            if backup_cache is None:
                os.environ.pop('ARCH_CACHE_LINE_SIZE', None)
            else:
                os.environ['ARCH_CACHE_LINE_SIZE'] = backup_cache

    if result is None or result.faces.size == 0:
        return _native._geom3d().solid([], [], ['boolean', f'{ENGINE_NAME}:{operation}'])

    return _mesh_to_solid(result, operation)


__all__ = ['ENGINE_NAME', 'is_available', 'solid_boolean', 'engines_available']
