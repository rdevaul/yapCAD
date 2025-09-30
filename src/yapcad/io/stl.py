"""STL export utilities for yapCAD surfaces and solids."""

from __future__ import annotations

import io
import os
import struct
from typing import Iterable, Sequence, TextIO

from yapcad.geometry_utils import Triangle, triangles_from_mesh
from yapcad.mesh import mesh_view

_HEADER_SIZE = 80
_STRUCT_TRIANGLE = struct.Struct('<12fH')


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


__all__ = ['write_stl']
