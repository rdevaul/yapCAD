from __future__ import annotations

from typing import Any, Dict, List, Tuple

from yapcad.geom3d import issolid, issurface


def _collect_surfaces(entity: Any) -> List[list]:
    if issurface(entity):
        return [entity]
    if issolid(entity):
        # solid structure: ['solid', [surfaces], voids, ...]
        return [s for s in (entity[1] or []) if issurface(s)]
    return []


def solidish_to_mesh(entities: List[Any]) -> Dict[str, List]:
    """Convert yapCAD surface/solid entities into a single indexed triangle mesh.

    Output format:
      - vertices: flat list [x,y,z,...]
      - normals: flat list [nx,ny,nz,...]
      - indices: flat list of triangle indices [i0,i1,i2,...]

    Notes:
      - If normals are missing, zeros are emitted.
      - Sketch/curve entities are ignored (mesh will be empty).
    """

    vertices: List[float] = []
    normals: List[float] = []
    indices: List[int] = []

    vert_offset = 0

    for ent in entities:
        for surf in _collect_surfaces(ent):
            verts = surf[1] or []
            norms = surf[2] or []
            faces = surf[3] or []

            # vertices/normals
            for i, v in enumerate(verts):
                vertices.extend([float(v[0]), float(v[1]), float(v[2])])
                if i < len(norms):
                    n = norms[i]
                    normals.extend([float(n[0]), float(n[1]), float(n[2])])
                else:
                    normals.extend([0.0, 0.0, 0.0])

            # indices
            for f in faces:
                if len(f) != 3:
                    continue
                indices.extend([vert_offset + int(f[0]), vert_offset + int(f[1]), vert_offset + int(f[2])])

            vert_offset += len(verts)

    return {"vertices": vertices, "normals": normals, "indices": indices}
