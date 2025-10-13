"""Shared mesh-diagnostics helpers for yapCAD solids."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional

try:
    import numpy as np
except ImportError as exc:  # pragma: no cover - optional dependency
    raise SystemExit("numpy is required for mesh diagnostics (pip install numpy)") from exc

try:
    import trimesh
except ImportError as exc:  # pragma: no cover - optional dependency
    raise SystemExit("trimesh is required for mesh diagnostics (pip install trimesh[easy])") from exc


@dataclass
class MeshDiagnostics:
    triangles: int
    vertices: int
    is_watertight: bool
    is_volume: bool
    euler_number: int
    boundary_edges: int
    self_intersections: int
    min_edge_length: Optional[float]
    max_edge_length: Optional[float]
    orphan_vertices: int

    def to_dict(self) -> dict:
        return asdict(self)


def _boundary_mask(mesh: "trimesh.Trimesh") -> Optional[np.ndarray]:
    try:
        unique = mesh.edges_unique
        inverse = mesh.edges_unique_inverse
        if unique is None or inverse is None:
            return None
        counts = np.bincount(inverse, minlength=len(unique))
        if counts.size == 0:
            return np.zeros(0, dtype=bool)
        return counts == 1
    except Exception:
        return None


def diagnose_mesh(mesh: "trimesh.Trimesh") -> MeshDiagnostics:
    boundary_mask = _boundary_mask(mesh)
    boundary_count = int(boundary_mask.sum()) if boundary_mask is not None else 0

    self_intersections = 0
    try:
        faces = mesh.self_intersecting_faces
        if faces is not None:
            self_intersections = int(len(faces))
    except Exception:
        pass

    min_len = max_len = None
    try:
        lengths = mesh.edges_unique_length
        if lengths is not None and len(lengths) > 0:
            min_len = float(np.min(lengths))
            max_len = float(np.max(lengths))
    except Exception:
        pass

    orphan_vertices = 0
    try:
        vf = mesh.vertex_faces
        if vf is not None:
            orphan_vertices = int((vf < 0).sum())
    except Exception:
        orphan_vertices = 0

    is_volume = bool(getattr(mesh, "is_volume", mesh.is_watertight))

    return MeshDiagnostics(
        triangles=int(len(mesh.faces)) if mesh.faces is not None else 0,
        vertices=int(len(mesh.vertices)) if mesh.vertices is not None else 0,
        is_watertight=bool(mesh.is_watertight),
        is_volume=is_volume,
        euler_number=int(getattr(mesh, "euler_number", 0)),
        boundary_edges=boundary_count,
        self_intersections=int(self_intersections),
        min_edge_length=min_len,
        max_edge_length=max_len,
        orphan_vertices=int(orphan_vertices),
    )


def export_boundary(mesh: "trimesh.Trimesh", path: Path) -> Optional[Path]:
    mask = _boundary_mask(mesh)
    if mask is None or not mask.any():
        return None
    segments = mesh.edges_unique[mask]
    if segments.size == 0:
        return None
    line = trimesh.load_path(mesh.vertices[segments])
    path.parent.mkdir(parents=True, exist_ok=True)
    line.export(str(path))
    return path

