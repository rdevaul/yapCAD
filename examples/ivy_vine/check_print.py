"""Check an STL of the vine for the things a supportless print needs.

    python check_print.py out/ivy_vine.stl

Reports, from the triangles alone: how many separate bodies there are and
whether each is watertight; whether everything sits on the bed; and how much
surface overhangs more steeply than 45 degrees without being on the bed.
"""

import sys

import numpy as np
import trimesh

BED = 0.06            # anything this close to z = 0 counts as on the bed
LIMIT_DEG = 45.0      # steepest overhang printed without support
SLACK_DEG = 3.0       # meshing noise on a nominal 45 degree face


def body_labels(mesh):
    """Label each face with the connected body it belongs to."""
    parent = np.arange(len(mesh.vertices))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for a, b in np.vstack([mesh.faces[:, [0, 1]], mesh.faces[:, [1, 2]]]):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    roots = np.array([find(i) for i in mesh.faces[:, 0]])
    return np.unique(roots, return_inverse=True)[1]


def main(path):
    mesh = trimesh.load(path, force="mesh")
    lo, hi = mesh.bounds
    print(f"{path}: {len(mesh.faces)} triangles, "
          f"{hi[0]-lo[0]:.1f} x {hi[1]-lo[1]:.1f} x {hi[2]-lo[2]:.1f} mm, "
          f"z {lo[2]:+.3f} .. {hi[2]:+.3f}")

    labels = body_labels(mesh)
    count = int(labels.max()) + 1
    watertight = floating = 0
    volume = 0.0
    for body in range(count):
        faces = mesh.faces[labels == body]
        tris = mesh.vertices[faces]
        edges = np.sort(faces[:, [0, 1, 1, 2, 2, 0]].reshape(-1, 2), axis=1)
        _, uses = np.unique(edges, axis=0, return_counts=True)
        watertight += bool((uses == 2).all())
        floating += bool(tris[:, :, 2].min() > BED)
        volume += float(np.einsum("ij,ij->i", tris[:, 0],
                                  np.cross(tris[:, 1], tris[:, 2])).sum()) / 6.0
    print(f"bodies: {count}, watertight: {watertight}")
    print(f"bodies not touching the bed: {floating}")
    print(f"volume: {volume / 1000.0:.1f} cm^3")

    normals = mesh.face_normals
    areas = mesh.area_faces
    centres = mesh.triangles_center
    # Angle of a downward face from vertical: 90 is a flat ceiling.
    down = normals[:, 2] < 0.0
    overhang = np.degrees(np.arcsin(np.clip(-normals[:, 2], 0.0, 1.0)))
    on_bed = centres[:, 2] < BED
    steep = down & ~on_bed & (overhang > LIMIT_DEG + SLACK_DEG)
    print(f"bed contact area: {areas[down & on_bed & (overhang > 80)].sum():.0f} mm^2")
    print(f"overhang steeper than {LIMIT_DEG:.0f} deg, off the bed: "
          f"{areas[steep].sum():.2f} mm^2 in {int(steep.sum())} triangles")
    if steep.any():
        # Group by height to say where it is.
        z = centres[steep, 2]
        for z0 in np.arange(0.0, hi[2] + 0.5, 0.5):
            band = (z >= z0) & (z < z0 + 0.5)
            if band.any():
                a = areas[steep][band].sum()
                worst = overhang[steep][band].max()
                print(f"   z {z0:.1f}-{z0+0.5:.1f}: {a:7.2f} mm^2, "
                      f"up to {worst:.0f} deg")
    return mesh, steep


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "out/ivy_vine.stl")
