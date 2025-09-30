"""Common geometric helpers shared across exporters and validators."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence, Tuple

from yapcad.geom import cross, epsilon, mag

Vec3 = Tuple[float, float, float]


@dataclass(frozen=True)
class Triangle:
    """Immutable triangle representation in XYZ space."""

    normal: Vec3
    v0: Vec3
    v1: Vec3
    v2: Vec3


def to_vec3(point_like: Sequence[float]) -> Vec3:
    """Return the XYZ components of a yapCAD point/vector as a tuple."""

    if len(point_like) < 3:
        raise ValueError("value must have at least three components")
    return float(point_like[0]), float(point_like[1]), float(point_like[2])


def to_point(vec: Sequence[float]) -> Tuple[float, float, float, float]:
    """Lift an XYZ tuple into yapCAD homogeneous point form."""

    if len(vec) < 3:
        raise ValueError("vector must have three components")
    return float(vec[0]), float(vec[1]), float(vec[2]), 1.0


def to_vector(vec: Sequence[float]) -> Tuple[float, float, float, float]:
    """Lift an XYZ tuple into yapCAD homogeneous vector form (w=0)."""

    if len(vec) < 3:
        raise ValueError("vector must have three components")
    return float(vec[0]), float(vec[1]), float(vec[2]), 0.0


def triangle_normal(v0: Vec3, v1: Vec3, v2: Vec3) -> Vec3 | None:
    """Return the unit normal of a triangle or ``None`` if degenerate."""

    ax, ay, az = v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]
    bx, by, bz = v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]
    n = cross([ax, ay, az, 1.0], [bx, by, bz, 1.0])
    length = mag(n)
    if length <= epsilon:
        return None
    return (n[0] / length, n[1] / length, n[2] / length)


def triangle_area(v0: Vec3, v1: Vec3, v2: Vec3) -> float:
    """Return the area of a triangle."""

    ax, ay, az = v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]
    bx, by, bz = v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]
    n = cross([ax, ay, az, 1.0], [bx, by, bz, 1.0])
    return 0.5 * mag(n)


def triangle_is_degenerate(v0: Vec3, v1: Vec3, v2: Vec3, tol: float = epsilon) -> bool:
    """Return ``True`` if a triangle collapses under the given tolerance."""

    return triangle_area(v0, v1, v2) <= tol


def orient_triangle(v0: Vec3, v1: Vec3, v2: Vec3, preferred_normal: Vec3) -> Tuple[Vec3, Vec3, Vec3]:
    """Ensure triangle winding aligns with ``preferred_normal``."""

    current = triangle_normal(v0, v1, v2)
    if current is None:
        return v0, v1, v2
    dot = current[0] * preferred_normal[0] + current[1] * preferred_normal[1] + current[2] * preferred_normal[2]
    if dot < 0:
        return v0, v2, v1
    return v0, v1, v2


def triangle_centroid(v0: Vec3, v1: Vec3, v2: Vec3) -> Vec3:
    """Return the centroid of a triangle."""

    return (
        (v0[0] + v1[0] + v2[0]) / 3.0,
        (v0[1] + v1[1] + v2[1]) / 3.0,
        (v0[2] + v1[2] + v2[2]) / 3.0,
    )


def triangles_from_mesh(mesh: Iterable[Tuple[Vec3, Vec3, Vec3, Vec3]]) -> Iterable[Triangle]:
    """Convert ``mesh_view`` output into ``Triangle`` instances."""

    for normal, v0, v1, v2 in mesh:
        yield Triangle(normal=normal, v0=v0, v1=v1, v2=v2)


__all__ = [
    "Triangle",
    "Vec3",
    "to_vec3",
    "to_point",
    "to_vector",
    "triangle_normal",
    "triangle_area",
    "triangle_is_degenerate",
    "orient_triangle",
    "triangle_centroid",
    "triangles_from_mesh",
]
