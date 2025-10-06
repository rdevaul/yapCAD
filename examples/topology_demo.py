#!/usr/bin/env python3
"""
Demonstration of yapCAD solid topology analysis functions.

This script demonstrates the new issolidclosed() and volumeof() functions
for analyzing 3D solid geometry.
"""

from yapcad.geom import *
from yapcad.geom3d import *
from yapcad.geom3d_util import *

def main():
    print("yapCAD Solid Topology Analysis Demo")
    print("=" * 50)
    print()

    # Example 1: Simple cube
    print("Example 1: Unit Cube")
    print("-" * 30)
    cube = prism(1, 1, 1)
    is_closed = issolidclosed(cube)
    print(f"  Is closed: {is_closed}")
    if is_closed:
        vol = volumeof(cube)
        print(f"  Volume: {vol:.6f}")
        print(f"  Expected: 1.0")
    print()

    # Example 2: 2x3x4 Rectangular prism
    print("Example 2: Rectangular Prism (2×3×4)")
    print("-" * 30)
    box = prism(2, 3, 4)
    is_closed = issolidclosed(box)
    print(f"  Is closed: {is_closed}")
    if is_closed:
        vol = volumeof(box)
        print(f"  Volume: {vol:.6f}")
        print(f"  Expected: 24.0")
    print()

    # Example 3: Sphere
    print("Example 3: Sphere (radius 2.0)")
    print("-" * 30)
    import math
    radius = 2.0
    sph = sphere(2 * radius, depth=3)
    is_closed = issolidclosed(sph)
    print(f"  Is closed: {is_closed}")
    if is_closed:
        vol = volumeof(sph)
        expected = (4.0/3.0) * math.pi * (radius ** 3)
        error = abs(vol - expected) / expected * 100
        print(f"  Volume: {vol:.6f}")
        print(f"  Expected: {expected:.6f}")
        print(f"  Error: {error:.2f}%")
    print()

    # Example 4: Extruded square with hole
    print("Example 4: Extruded Square with Hole")
    print("-" * 30)
    outer = [point(0, 0), point(4, 0), point(4, 4), point(0, 4), point(0, 0)]
    hole = [point(1, 1), point(3, 1), point(3, 3), point(1, 3), point(1, 1)]
    surf, _ = poly2surfaceXY(outer, holepolys=[hole])
    extruded = extrude(surf, 2.0)
    is_closed = issolidclosed(extruded)
    print(f"  Is closed: {is_closed}")
    if is_closed:
        vol = volumeof(extruded)
        # (4*4 - 2*2) * 2 = 24
        expected = (4 * 4 - 2 * 2) * 2
        print(f"  Volume: {vol:.6f}")
        print(f"  Expected: {expected:.1f}")
    print()

    # Example 5: Cylinder
    print("Example 5: Cylinder (radius 2, height 5)")
    print("-" * 30)
    radius = 2.0
    height = 5.0
    cylinder = conic(radius, radius, height)
    is_closed = issolidclosed(cylinder)
    print(f"  Is closed: {is_closed}")
    if is_closed:
        vol = volumeof(cylinder)
        expected = math.pi * (radius ** 2) * height
        error = abs(vol - expected) / expected * 100
        print(f"  Volume: {vol:.6f}")
        print(f"  Expected: {expected:.6f}")
        print(f"  Error: {error:.2f}%")
    print()

    # Example 6: Tetrahedron
    print("Example 6: Regular Tetrahedron")
    print("-" * 30)
    # Create tetrahedron vertices
    tet_verts = [
        point(1, 1, 1),
        point(-1, -1, 1),
        point(-1, 1, -1),
        point(1, -1, -1)
    ]

    # Create faces
    faces = [
        [0, 2, 1],
        [0, 1, 3],
        [1, 2, 3],
        [2, 0, 3]
    ]

    # Create surface
    all_verts = []
    all_norms = []
    for f_idx, f in enumerate(faces):
        tri = [tet_verts[f[0]], tet_verts[f[1]], tet_verts[f[2]]]
        p0, n = tri2p0n(tri)
        for v_idx in f:
            all_verts.append(tet_verts[v_idx])
            all_norms.append(n)

    new_faces = []
    for i in range(len(faces)):
        new_faces.append([i * 3, i * 3 + 1, i * 3 + 2])

    surf_tet = surface(all_verts, all_norms, new_faces)
    tet = solid([surf_tet])

    is_closed = issolidclosed(tet)
    print(f"  Is closed: {is_closed}")
    if is_closed:
        vol = volumeof(tet)
        # Calculate edge length and expected volume: V = a³/(6√2)
        edge = line(tet_verts[0], tet_verts[1])
        a = length(edge)
        expected = (a ** 3) / (6 * math.sqrt(2))
        error = abs(vol - expected) / expected * 100
        print(f"  Edge length: {a:.6f}")
        print(f"  Volume: {vol:.6f}")
        print(f"  Expected: {expected:.6f}")
        print(f"  Error: {error:.4f}%")
    print()

    # Example 7: Open box (missing top) - should NOT be closed
    print("Example 7: Open Box (missing top)")
    print("-" * 30)
    bottom = rectangularPlane(2, 2, point(0, 0, 0))

    side1 = rectangularPlane(2, 1, point(0, 0, 0.5))
    side1 = rotatesurface(side1, 90, axis=point(1, 0, 0))
    side1 = translatesurface(side1, point(0, -1, 0))

    side2 = rectangularPlane(2, 1, point(0, 0, 0.5))
    side2 = rotatesurface(side2, 90, axis=point(1, 0, 0))
    side2 = translatesurface(side2, point(0, 1, 0))

    side3 = rectangularPlane(2, 1, point(0, 0, 0.5))
    side3 = rotatesurface(side3, 90, axis=point(0, 1, 0))
    side3 = translatesurface(side3, point(-1, 0, 0))

    side4 = rectangularPlane(2, 1, point(0, 0, 0.5))
    side4 = rotatesurface(side4, 90, axis=point(0, 1, 0))
    side4 = translatesurface(side4, point(1, 0, 0))

    open_box = solid([bottom, side1, side2, side3, side4])
    is_closed = issolidclosed(open_box)
    print(f"  Is closed: {is_closed}")
    print(f"  (Expected: False - box has open top)")
    if is_closed:
        vol = volumeof(open_box)
        print(f"  Volume: {vol:.6f}")
    else:
        print(f"  Cannot compute volume - solid not closed")
    print()

    print("Demo complete!")
    print("=" * 50)

if __name__ == "__main__":
    main()
