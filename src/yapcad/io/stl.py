"""STL import and export utilities for yapCAD surfaces and solids."""

from __future__ import annotations

import io
import os
import re
import struct
from typing import Iterable, List, Sequence, TextIO, Tuple, Union

from yapcad.geom import point, vect
from yapcad.geom3d import surface, solid
from yapcad.geometry_utils import Triangle, triangles_from_mesh, triangle_normal
from yapcad.mesh import mesh_view

_HEADER_SIZE = 80
_STRUCT_TRIANGLE = struct.Struct('<12fH')
_VERTEX_TOL = 1e-9  # Tolerance for vertex deduplication


def write_stl(obj: Sequence, path_or_file, *, binary: bool = True, name: str = 'yapCAD') -> None:
    """Write ``obj`` (surface or solid) to STL.

    ``path_or_file`` can be a filesystem path or an open binary/text stream.
    """

    triangles = list(triangles_from_mesh(mesh_view(obj)))

    if binary:
        _write_binary(triangles, path_or_file, name)
    else:
        _write_ascii(triangles, path_or_file, name)


def _write_binary(triangles: Iterable[Triangle], path_or_file, name: str) -> None:
    close_when_done = False
    if hasattr(path_or_file, 'write'):
        stream = path_or_file
    else:
        stream = open(path_or_file, 'wb')
        close_when_done = True

    try:
        header = (name[:_HEADER_SIZE]).encode('ascii', errors='replace')
        header = header.ljust(_HEADER_SIZE, b' ')
        stream.write(header)
        stream.write(struct.pack('<I', len(triangles)))

        for tri in triangles:
            data = _STRUCT_TRIANGLE.pack(
                *tri.normal,
                *tri.v0,
                *tri.v1,
                *tri.v2,
                0,
            )
            stream.write(data)
    finally:
        if close_when_done:
            stream.close()


def _write_ascii(triangles: Iterable[Triangle], path_or_file, name: str) -> None:
    close_when_done = False
    if hasattr(path_or_file, 'write'):
        stream = path_or_file
    else:
        stream = open(path_or_file, 'w', encoding='ascii')
        close_when_done = True

    try:
        print(f"solid {name}", file=stream)
        for tri in triangles:
            print(f"  facet normal {tri.normal[0]:.6e} {tri.normal[1]:.6e} {tri.normal[2]:.6e}", file=stream)
            print("    outer loop", file=stream)
            print(f"      vertex {tri.v0[0]:.6e} {tri.v0[1]:.6e} {tri.v0[2]:.6e}", file=stream)
            print(f"      vertex {tri.v1[0]:.6e} {tri.v1[1]:.6e} {tri.v1[2]:.6e}", file=stream)
            print(f"      vertex {tri.v2[0]:.6e} {tri.v2[1]:.6e} {tri.v2[2]:.6e}", file=stream)
            print("    endloop", file=stream)
            print("  endfacet", file=stream)
        print(f"endsolid {name}", file=stream)
    finally:
        if close_when_done:
            stream.close()


# ---------------------------------------------------------------------------
# STL Import
# ---------------------------------------------------------------------------


def _is_binary_stl(data: bytes) -> bool:
    """Determine if STL data is binary format.

    Binary STL has 80-byte header + 4-byte count, then 50 bytes per triangle.
    ASCII STL starts with 'solid' keyword.
    """
    if len(data) < 84:
        return False

    # Check if it looks like ASCII (starts with 'solid')
    try:
        header = data[:80].decode('ascii', errors='ignore').strip().lower()
        if header.startswith('solid'):
            # Could still be binary if 'solid' is just in header
            # Check if the file size matches binary format expectation
            tri_count = struct.unpack('<I', data[80:84])[0]
            expected_size = 84 + (tri_count * 50)
            if len(data) == expected_size:
                # Size matches binary format - but double check for ASCII markers
                rest = data[84:min(200, len(data))]
                if b'facet' in rest or b'vertex' in rest:
                    return False
                return True
            # Size doesn't match, likely ASCII
            return False
        return True
    except (UnicodeDecodeError, struct.error):
        return True


def _parse_binary_stl(data: bytes) -> List[Triangle]:
    """Parse binary STL data into triangles."""
    if len(data) < 84:
        raise ValueError("Invalid binary STL: file too small")

    tri_count = struct.unpack('<I', data[80:84])[0]
    triangles = []
    offset = 84

    for _ in range(tri_count):
        if offset + 50 > len(data):
            break
        values = _STRUCT_TRIANGLE.unpack(data[offset:offset + 50])
        normal = (values[0], values[1], values[2])
        v0 = (values[3], values[4], values[5])
        v1 = (values[6], values[7], values[8])
        v2 = (values[9], values[10], values[11])
        triangles.append(Triangle(normal=normal, v0=v0, v1=v1, v2=v2))
        offset += 50

    return triangles


def _parse_ascii_stl(text: str) -> List[Triangle]:
    """Parse ASCII STL text into triangles."""
    triangles = []

    # Pattern for facet block
    facet_pattern = re.compile(
        r'facet\s+normal\s+([eE\d.+-]+)\s+([eE\d.+-]+)\s+([eE\d.+-]+)\s+'
        r'outer\s+loop\s+'
        r'vertex\s+([eE\d.+-]+)\s+([eE\d.+-]+)\s+([eE\d.+-]+)\s+'
        r'vertex\s+([eE\d.+-]+)\s+([eE\d.+-]+)\s+([eE\d.+-]+)\s+'
        r'vertex\s+([eE\d.+-]+)\s+([eE\d.+-]+)\s+([eE\d.+-]+)\s+'
        r'endloop\s+endfacet',
        re.IGNORECASE
    )

    for match in facet_pattern.finditer(text):
        groups = match.groups()
        normal = (float(groups[0]), float(groups[1]), float(groups[2]))
        v0 = (float(groups[3]), float(groups[4]), float(groups[5]))
        v1 = (float(groups[6]), float(groups[7]), float(groups[8]))
        v2 = (float(groups[9]), float(groups[10]), float(groups[11]))
        triangles.append(Triangle(normal=normal, v0=v0, v1=v1, v2=v2))

    return triangles


def _vertex_key(v: Tuple[float, float, float], tol: float = _VERTEX_TOL) -> Tuple[int, int, int]:
    """Create a hashable key for vertex deduplication."""
    scale = 1.0 / tol
    return (int(round(v[0] * scale)), int(round(v[1] * scale)), int(round(v[2] * scale)))


def _triangles_to_surface(triangles: List[Triangle], deduplicate: bool = True) -> list:
    """Convert triangles to a yapCAD surface.

    Parameters
    ----------
    triangles : list of Triangle
        Input triangles with normals and vertices.
    deduplicate : bool
        If True, merge coincident vertices to create a proper indexed mesh.
        If False, each triangle gets its own vertices (3 * len(triangles) vertices).

    Returns
    -------
    surface
        yapCAD surface: ['surface', vertices, normals, faces, boundary, holes]
    """
    if not triangles:
        return surface()

    vertices = []
    normals = []
    faces = []

    if deduplicate:
        # Build vertex index with deduplication
        vertex_map = {}  # key -> vertex index
        vertex_normals = {}  # key -> list of normals for averaging

        for tri in triangles:
            face_indices = []
            for v in (tri.v0, tri.v1, tri.v2):
                key = _vertex_key(v)
                if key not in vertex_map:
                    vertex_map[key] = len(vertices)
                    vertices.append(point(v[0], v[1], v[2]))
                    vertex_normals[key] = []
                vertex_normals[key].append(tri.normal)
                face_indices.append(vertex_map[key])
            faces.append(face_indices)

        # Average normals at each vertex
        for key in sorted(vertex_map.keys(), key=lambda k: vertex_map[k]):
            norms = vertex_normals[key]
            avg_n = [
                sum(n[0] for n in norms) / len(norms),
                sum(n[1] for n in norms) / len(norms),
                sum(n[2] for n in norms) / len(norms),
            ]
            # Normalize
            length = (avg_n[0]**2 + avg_n[1]**2 + avg_n[2]**2) ** 0.5
            if length > 1e-10:
                avg_n = [avg_n[0]/length, avg_n[1]/length, avg_n[2]/length]
            normals.append(vect(avg_n[0], avg_n[1], avg_n[2], 0))
    else:
        # No deduplication - each triangle gets separate vertices
        for i, tri in enumerate(triangles):
            base = i * 3
            vertices.append(point(tri.v0[0], tri.v0[1], tri.v0[2]))
            vertices.append(point(tri.v1[0], tri.v1[1], tri.v1[2]))
            vertices.append(point(tri.v2[0], tri.v2[1], tri.v2[2]))
            normals.append(vect(tri.normal[0], tri.normal[1], tri.normal[2], 0))
            normals.append(vect(tri.normal[0], tri.normal[1], tri.normal[2], 0))
            normals.append(vect(tri.normal[0], tri.normal[1], tri.normal[2], 0))
            faces.append([base, base + 1, base + 2])

    return surface(vertices, normals, faces, [], [])


def read_stl(path_or_file, *, deduplicate: bool = True) -> list:
    """Read an STL file and return a yapCAD solid.

    Parameters
    ----------
    path_or_file : str or path-like or file-like
        Path to STL file, or an open binary file object.
    deduplicate : bool, optional
        If True (default), merge coincident vertices to create a proper
        indexed mesh. If False, each triangle gets its own vertices.

    Returns
    -------
    solid
        yapCAD solid containing a single surface with the imported mesh.
        The solid structure is: ['solid', [surface], boundary, holes]

    Examples
    --------
    >>> from yapcad.io.stl import read_stl, write_stl
    >>> my_solid = read_stl('model.stl')
    >>> # Round-trip: export and re-import
    >>> write_stl(my_solid, 'copy.stl')
    >>> copy_solid = read_stl('copy.stl')
    """
    close_when_done = False

    if hasattr(path_or_file, 'read'):
        data = path_or_file.read()
        if isinstance(data, str):
            data = data.encode('utf-8')
    else:
        with open(path_or_file, 'rb') as f:
            data = f.read()

    if _is_binary_stl(data):
        triangles = _parse_binary_stl(data)
    else:
        text = data.decode('utf-8', errors='replace')
        triangles = _parse_ascii_stl(text)

    if not triangles:
        # Return an empty solid (empty surfaces list)
        return ['solid', [], [], []]

    surf = _triangles_to_surface(triangles, deduplicate=deduplicate)
    return solid([surf], [], [])


def import_stl(path: str, *, deduplicate: bool = True) -> list:
    """Import an STL file as a yapCAD solid.

    This is an alias for :func:`read_stl` for API consistency with
    :func:`yapcad.io.step_importer.import_step`.

    Parameters
    ----------
    path : str
        Path to STL file.
    deduplicate : bool, optional
        If True (default), merge coincident vertices.

    Returns
    -------
    solid
        yapCAD solid containing the imported mesh.
    """
    return read_stl(path, deduplicate=deduplicate)


__all__ = ['write_stl', 'read_stl', 'import_stl']
