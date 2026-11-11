"""Face naming system for boundary condition assignment.

This module provides utilities for naming faces of yapCAD solids, which can
then be used to assign boundary conditions in analysis plans.

Face names can be assigned:
1. At creation time via DSL `with { face_names: {...} }` syntax
2. Post-hoc via selectors (by normal, by area, by position)
3. Interactively in the viewer (future)

The face names are stored in solid metadata and propagate through to
Gmsh physical groups when meshing.

Copyright (c) 2025 yapCAD contributors
MIT License
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

from yapcad.geom import point, vect


@dataclass
class FaceInfo:
    """Information about a face for selection/naming.

    Attributes:
        index: Face index in the solid
        centroid: Area-weighted face centroid
        normal: Area-weighted mean outward normal, normalised.  A patch that
            faces every way at once -- a whole cylinder wall or sphere --
            has no mean direction and gets ``(0, 0, 0)``, which no
            ``ByNormalSelector`` matches.
        area: Face area
        name: Optional assigned name
        tags: Optional list of tags
        triangles: ``(surface_index, face_index)`` pairs locating the face's
            triangles in the solid it was extracted from
    """
    index: int
    centroid: Tuple[float, float, float]
    normal: Tuple[float, float, float]
    area: float
    name: Optional[str] = None
    tags: List[str] = field(default_factory=list)
    triangles: List[Tuple[int, int]] = field(default_factory=list)


class FaceSelector:
    """Base class for face selection predicates."""

    def matches(self, face: FaceInfo) -> bool:
        """Return True if face matches this selector."""
        raise NotImplementedError


class ByNormalSelector(FaceSelector):
    """Select faces by normal direction."""

    def __init__(
        self,
        direction: Tuple[float, float, float],
        tolerance_deg: float = 5.0,
        allow_reversed: bool = False,
    ):
        """
        Args:
            direction: Target normal direction (will be normalized)
            tolerance_deg: Angular tolerance in degrees
            allow_reversed: Also match faces with reversed normal
        """
        mag = math.sqrt(sum(d*d for d in direction))
        self.direction = tuple(d/mag for d in direction) if mag > 1e-10 else (0, 0, 1)
        self.tolerance_rad = math.radians(tolerance_deg)
        self.allow_reversed = allow_reversed

    def matches(self, face: FaceInfo) -> bool:
        dot = sum(a*b for a, b in zip(self.direction, face.normal))
        angle = math.acos(max(-1, min(1, abs(dot) if self.allow_reversed else dot)))
        return angle <= self.tolerance_rad


class ByAreaSelector(FaceSelector):
    """Select faces by area criteria."""

    def __init__(
        self,
        min_area: Optional[float] = None,
        max_area: Optional[float] = None,
        largest: bool = False,
        smallest: bool = False,
    ):
        """
        Args:
            min_area: Minimum area threshold
            max_area: Maximum area threshold
            largest: Select only the largest face
            smallest: Select only the smallest face
        """
        self.min_area = min_area
        self.max_area = max_area
        self.largest = largest
        self.smallest = smallest
        self._all_faces: List[FaceInfo] = []

    def set_context(self, faces: List[FaceInfo]) -> None:
        """Set all faces for largest/smallest comparison."""
        self._all_faces = faces

    def matches(self, face: FaceInfo) -> bool:
        if self.min_area is not None and face.area < self.min_area:
            return False
        if self.max_area is not None and face.area > self.max_area:
            return False

        if self.largest and self._all_faces:
            max_area = max(f.area for f in self._all_faces)
            return math.isclose(face.area, max_area, rel_tol=1e-6)

        if self.smallest and self._all_faces:
            min_area = min(f.area for f in self._all_faces)
            return math.isclose(face.area, min_area, rel_tol=1e-6)

        return True


class ByPositionSelector(FaceSelector):
    """Select faces by centroid position."""

    def __init__(
        self,
        axis: str = "z",
        at_min: bool = False,
        at_max: bool = False,
        above: Optional[float] = None,
        below: Optional[float] = None,
        tolerance: float = 1e-6,
    ):
        """
        Args:
            axis: Coordinate axis ("x", "y", or "z")
            at_min: Select face(s) at minimum coordinate
            at_max: Select face(s) at maximum coordinate
            above: Select faces with coordinate above this value
            below: Select faces with coordinate below this value
            tolerance: Tolerance for at_min/at_max comparison
        """
        self.axis_idx = {"x": 0, "y": 1, "z": 2}.get(axis.lower(), 2)
        self.at_min = at_min
        self.at_max = at_max
        self.above = above
        self.below = below
        self.tolerance = tolerance
        self._all_faces: List[FaceInfo] = []

    def set_context(self, faces: List[FaceInfo]) -> None:
        """Set all faces for min/max comparison."""
        self._all_faces = faces

    def matches(self, face: FaceInfo) -> bool:
        coord = face.centroid[self.axis_idx]

        if self.above is not None and coord <= self.above:
            return False
        if self.below is not None and coord >= self.below:
            return False

        if self.at_min and self._all_faces:
            min_coord = min(f.centroid[self.axis_idx] for f in self._all_faces)
            return abs(coord - min_coord) <= self.tolerance

        if self.at_max and self._all_faces:
            max_coord = max(f.centroid[self.axis_idx] for f in self._all_faces)
            return abs(coord - max_coord) <= self.tolerance

        return True


class CombinedSelector(FaceSelector):
    """Combine multiple selectors with AND/OR logic."""

    def __init__(self, selectors: List[FaceSelector], mode: str = "and"):
        """
        Args:
            selectors: List of selectors to combine
            mode: "and" (all must match) or "or" (any must match)
        """
        self.selectors = selectors
        self.mode = mode.lower()

    def set_context(self, faces: List[FaceInfo]) -> None:
        """Propagate context to child selectors."""
        for sel in self.selectors:
            if hasattr(sel, 'set_context'):
                sel.set_context(faces)

    def matches(self, face: FaceInfo) -> bool:
        if self.mode == "or":
            return any(sel.matches(face) for sel in self.selectors)
        return all(sel.matches(face) for sel in self.selectors)


class FaceNamer:
    """Utility for naming faces of yapCAD solids.

    This class extracts face information from solids and applies names
    based on selectors or explicit assignments.

    A yapCAD surface is a triangle mesh with no notion of B-rep faces, and
    its boundaries are an accident of construction: ``prism`` makes one
    surface per side, but a mesh boolean or an SDF part is a single surface
    covering every side.  So faces are recovered from the triangles: a face
    is a maximal patch of triangles joined across edges whose dihedral angle
    is at most ``crease_deg``, found across the whole solid regardless of
    how its surfaces are divided.  A box gives six faces however it was
    built; a smooth cylinder wall stays one face, as in a B-rep, provided
    it is tessellated finer than ``crease_deg`` per segment.

    Face indices are in order of each face's first triangle, so they are
    deterministic for a given solid.  They are *mesh* faces: they match the
    B-rep face order Gmsh uses only where the two decompositions agree.
    """

    def __init__(self, solid: Any, crease_deg: float = 30.0):
        """
        Args:
            solid: yapCAD solid geometry
            crease_deg: Largest dihedral angle, in degrees, across which two
                triangles still belong to the same face

        Raises:
            ValueError: if ``solid`` is not a solid
            DeferredMeshError: if its mesh is deferred (realize it first)
        """
        self.solid = solid
        self.crease_deg = crease_deg
        self._faces: List[FaceInfo] = []
        self._extract_faces()

    def _extract_faces(self) -> None:
        """Group the solid's triangles into faces and measure each one."""
        from yapcad.geom3d import (
            issolid, _require_mesh, _triangle_frame, _canonical_edge_key,
        )

        if not issolid(self.solid):
            raise ValueError("FaceNamer requires a yapCAD solid")
        _require_mesh(self.solid, "FaceNamer")

        # Per triangle: id, points, unit normal (None if degenerate), area.
        tris = []
        for si, surf in enumerate(self.solid[1]):
            vertices = surf[1]
            for fi, face in enumerate(surf[3]):
                p0, p1, p2, cp = _triangle_frame(vertices, face)
                length = math.sqrt(cp[0]**2 + cp[1]**2 + cp[2]**2)
                unit = (cp[0]/length, cp[1]/length, cp[2]/length) if length > 0.0 else None
                tris.append(((si, fi), (p0, p1, p2), unit, length / 2.0, cp))

        # Edges by position, so faces join across surface boundaries and
        # duplicated vertices alike.
        edges: Dict[Any, List[int]] = {}
        for t, (_, pts, _, _, _) in enumerate(tris):
            for a, b in ((0, 1), (1, 2), (2, 0)):
                edges.setdefault(_canonical_edge_key(pts[a], pts[b]), []).append(t)

        cos_crease = math.cos(math.radians(self.crease_deg))

        def joined(t, u):
            # A degenerate triangle has no normal to compare, and letting it
            # join either side would bridge two faces across a crease.
            nt, nu = tris[t][2], tris[u][2]
            if nt is None or nu is None:
                return False
            return nt[0]*nu[0] + nt[1]*nu[1] + nt[2]*nu[2] >= cos_crease

        neighbours: List[List[int]] = [[] for _ in tris]
        for members in edges.values():
            for i, t in enumerate(members):
                for u in members[i + 1:]:
                    if t != u and joined(t, u):
                        neighbours[t].append(u)
                        neighbours[u].append(t)

        seen = [False] * len(tris)
        for seed in range(len(tris)):
            if seen[seed] or tris[seed][2] is None:
                continue
            seen[seed] = True
            group, stack = [], [seed]
            while stack:
                t = stack.pop()
                group.append(t)
                for u in neighbours[t]:
                    if not seen[u]:
                        seen[u] = True
                        stack.append(u)
            group.sort()
            self._faces.append(self._measure(len(self._faces), [tris[t] for t in group]))

    @staticmethod
    def _measure(index: int, members: List[Any]) -> FaceInfo:
        """Area, area-weighted centroid and mean normal of a triangle group."""
        area = 0.0
        centroid = [0.0, 0.0, 0.0]
        net = [0.0, 0.0, 0.0]
        twice_area = 0.0
        for _, pts, _, a, cp in members:
            area += a
            twice_area += 2.0 * a
            for k in range(3):
                centroid[k] += a * (pts[0][k] + pts[1][k] + pts[2][k]) / 3.0
                net[k] += cp[k]
        centroid = tuple(c / area for c in centroid)
        length = math.sqrt(net[0]**2 + net[1]**2 + net[2]**2)
        # The same threshold geom3d.surface_normal uses for "no direction".
        if length <= 1e-9 * twice_area:
            normal = (0.0, 0.0, 0.0)
        else:
            normal = (net[0]/length, net[1]/length, net[2]/length)
        return FaceInfo(
            index=index,
            centroid=centroid,
            normal=normal,
            area=area,
            triangles=[m[0] for m in members],
        )

    @property
    def faces(self) -> List[FaceInfo]:
        """Get list of face information."""
        return self._faces

    def name_faces(self, assignments: Dict[str, Union[FaceSelector, List[int]]]) -> Dict[str, List[int]]:
        """Apply face name assignments.

        Args:
            assignments: Mapping of names to selectors or explicit face indices

        Returns:
            Mapping of names to matched face indices
        """
        result: Dict[str, List[int]] = {}

        for name, selector in assignments.items():
            if isinstance(selector, list):
                # Explicit indices
                result[name] = selector
                for idx in selector:
                    if idx < len(self._faces):
                        self._faces[idx].name = name
            else:
                # Selector-based
                if hasattr(selector, 'set_context'):
                    selector.set_context(self._faces)

                matched = []
                for face in self._faces:
                    if selector.matches(face):
                        face.name = name
                        matched.append(face.index)

                result[name] = matched

        return result

    def get_named_faces(self) -> Dict[str, List[int]]:
        """Get all named faces.

        Returns:
            Mapping of names to face indices
        """
        named: Dict[str, List[int]] = {}
        for face in self._faces:
            if face.name:
                if face.name not in named:
                    named[face.name] = []
                named[face.name].append(face.index)
        return named

    def to_metadata(self) -> Dict[str, Any]:
        """Convert face naming to metadata format.

        Returns:
            Dictionary suitable for solid metadata
        """
        return {
            "face_names": self.get_named_faces(),
            "faces": [
                {
                    "index": f.index,
                    "centroid": f.centroid,
                    "normal": f.normal,
                    "area": f.area,
                    "name": f.name,
                    "tags": f.tags,
                    "triangles": f.triangles,
                }
                for f in self._faces
            ]
        }


# Convenience functions for common selections

def top_faces(tolerance_deg: float = 5.0) -> ByNormalSelector:
    """Select faces with +Z normal (top faces)."""
    return ByNormalSelector((0, 0, 1), tolerance_deg)


def bottom_faces(tolerance_deg: float = 5.0) -> ByNormalSelector:
    """Select faces with -Z normal (bottom faces)."""
    return ByNormalSelector((0, 0, -1), tolerance_deg)


def front_faces(tolerance_deg: float = 5.0) -> ByNormalSelector:
    """Select faces with +Y normal (front faces)."""
    return ByNormalSelector((0, 1, 0), tolerance_deg)


def back_faces(tolerance_deg: float = 5.0) -> ByNormalSelector:
    """Select faces with -Y normal (back faces)."""
    return ByNormalSelector((0, -1, 0), tolerance_deg)


def left_faces(tolerance_deg: float = 5.0) -> ByNormalSelector:
    """Select faces with -X normal (left faces)."""
    return ByNormalSelector((-1, 0, 0), tolerance_deg)


def right_faces(tolerance_deg: float = 5.0) -> ByNormalSelector:
    """Select faces with +X normal (right faces)."""
    return ByNormalSelector((1, 0, 0), tolerance_deg)


def largest_face() -> ByAreaSelector:
    """Select the largest face by area."""
    return ByAreaSelector(largest=True)


def smallest_face() -> ByAreaSelector:
    """Select the smallest face by area."""
    return ByAreaSelector(smallest=True)


def faces_at_z_min(tolerance: float = 1e-6) -> ByPositionSelector:
    """Select faces at minimum Z coordinate."""
    return ByPositionSelector(axis="z", at_min=True, tolerance=tolerance)


def faces_at_z_max(tolerance: float = 1e-6) -> ByPositionSelector:
    """Select faces at maximum Z coordinate."""
    return ByPositionSelector(axis="z", at_max=True, tolerance=tolerance)


__all__ = [
    "FaceInfo",
    "FaceSelector",
    "ByNormalSelector",
    "ByAreaSelector",
    "ByPositionSelector",
    "CombinedSelector",
    "FaceNamer",
    # Convenience functions
    "top_faces",
    "bottom_faces",
    "front_faces",
    "back_faces",
    "left_faces",
    "right_faces",
    "largest_face",
    "smallest_face",
    "faces_at_z_min",
    "faces_at_z_max",
]
