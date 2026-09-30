"""Minimal software renderer for triangle meshes, used by ``sdf_demo.py``.

yapCAD's interactive viewers need a window; this produces a PNG headlessly
and with no dependencies beyond numpy, so the SDF gallery can be regenerated
in CI or over ssh.  It is a z-buffered rasteriser with Gouraud shading, not a
serious renderer -- it exists so that a change to the mesher can be *looked
at* as well as asserted about.

Shading uses the per-vertex normals supplied by the caller, which for an SDF
mesh come from the analytic field gradient rather than from the triangles.
Where a triangle's own plane disagrees sharply with a vertex normal the
corner is shaded flat instead: a vertex normal cannot represent a crease,
because the gradient of a CSG field is genuinely discontinuous there, and
interpolating across it saws along every sharp edge.
"""

import struct
import zlib

import numpy as np

__all__ = ["render", "write_png"]


def write_png(path, rgb):
    """Write an ``(H, W, 3)`` uint8 array to ``path`` as a PNG."""
    height, width, _ = rgb.shape
    raw = b"".join(b"\x00" + rgb[y].tobytes() for y in range(height))

    def chunk(tag, data):
        body = tag + data
        return (struct.pack(">I", len(data)) + body
                + struct.pack(">I", zlib.crc32(body) & 0xFFFFFFFF))

    blob = (b"\x89PNG\r\n\x1a\n"
            + chunk(b"IHDR",
                    struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0))
            + chunk(b"IDAT", zlib.compress(raw, 6))
            + chunk(b"IEND", b""))
    with open(path, "wb") as handle:
        handle.write(blob)


def _camera(vertices, azimuth, elevation, zoom):
    """Frame the mesh: return the eye point, a view basis and its radius."""
    low = vertices.min(axis=0)
    high = vertices.max(axis=0)
    centre = 0.5 * (low + high)
    radius = float(np.linalg.norm(high - low)) * 0.5 or 1.0

    a = np.radians(azimuth)
    e = np.radians(elevation)
    direction = np.array([np.cos(e) * np.cos(a),
                          np.cos(e) * np.sin(a),
                          np.sin(e)])
    eye = centre + direction * radius * zoom

    up = np.array([0.0, 0.0, 1.0])
    forward = centre - eye
    forward /= np.linalg.norm(forward)
    if abs(np.dot(forward, up)) > 0.999:
        up = np.array([0.0, 1.0, 0.0])
    right = np.cross(forward, up)
    right /= np.linalg.norm(right)
    return eye, np.stack([right, np.cross(right, forward), forward]), radius


def render(vertices, normals, triangles, size=900, supersample=2,
           azimuth=35.0, elevation=26.0, zoom=2.9, fov=32.0,
           base=(0.62, 0.68, 0.78), background=(0.109, 0.121, 0.145),
           crease_angle=35.0):
    """Render a triangle mesh, returning an ``(size, size, 3)`` uint8 image.

    :param vertices: ``(V, 3)`` positions.
    :param normals: ``(V, 3)`` unit normals.
    :param triangles: ``(T, 3)`` vertex indices.
    :param crease_angle: degrees; beyond this a corner is shaded flat.
    """
    background = np.asarray(background, dtype=np.float32)
    if len(triangles) == 0:
        canvas = np.empty((size, size, 3), dtype=np.uint8)
        canvas[:] = (background * 255.0).astype(np.uint8)
        return canvas

    dim = size * supersample
    eye, basis, radius = _camera(vertices, azimuth, elevation, zoom)

    local = (vertices - eye) @ basis.T
    depth = local[:, 2]
    focal = (dim * 0.5) / np.tan(np.radians(fov) * 0.5)
    near = radius * 1e-3

    safe = np.maximum(depth, near)
    screen_x = local[:, 0] / safe * focal + dim * 0.5
    screen_y = dim * 0.5 - local[:, 1] / safe * focal

    colour = np.empty((dim, dim, 3), dtype=np.float32)
    colour[:] = background
    zbuffer = np.full((dim, dim), np.inf, dtype=np.float32)

    tri = triangles[np.all(depth[triangles] > near, axis=1)]
    x = screen_x[tri]
    y = screen_y[tri]
    z = depth[tri]

    area = ((x[:, 1] - x[:, 0]) * (y[:, 2] - y[:, 0])
            - (x[:, 2] - x[:, 0]) * (y[:, 1] - y[:, 0]))
    keep = np.abs(area) > 1e-9
    tri, x, y, z, area = tri[keep], x[keep], y[keep], z[keep], area[keep]
    if len(tri) == 0:
        return (np.clip(colour, 0, 1) ** (1 / 2.2)
                ).reshape(size, supersample, size, supersample,
                          3).mean((1, 3)).astype(np.uint8)

    edge_a = vertices[tri[:, 1]] - vertices[tri[:, 0]]
    edge_b = vertices[tri[:, 2]] - vertices[tri[:, 0]]
    face_n = np.cross(edge_a, edge_b)
    face_len = np.linalg.norm(face_n, axis=1, keepdims=True)
    face_n = np.divide(face_n, face_len, out=np.zeros_like(face_n),
                       where=face_len > 0)
    corner = normals[tri]
    sharp = np.einsum("tij,tj->ti", corner, face_n) < np.cos(
        np.radians(crease_angle))
    corner = np.where(sharp[:, :, None], face_n[:, None, :], corner)

    lo_x = np.clip(np.floor(x.min(axis=1)).astype(int), 0, dim - 1)
    hi_x = np.clip(np.ceil(x.max(axis=1)).astype(int), 0, dim - 1)
    lo_y = np.clip(np.floor(y.min(axis=1)).astype(int), 0, dim - 1)
    hi_y = np.clip(np.ceil(y.max(axis=1)).astype(int), 0, dim - 1)

    lights = []
    for vec, tint, strength in (
        ((0.45, 0.35, 0.82), (1.00, 0.97, 0.92), 0.85),
        ((-0.62, -0.30, 0.25), (0.45, 0.55, 0.78), 0.40),
    ):
        direction = np.asarray(vec, dtype=np.float32)
        lights.append((direction / np.linalg.norm(direction),
                       np.asarray(tint, dtype=np.float32), strength))
    base = np.asarray(base, dtype=np.float32)

    for t in range(len(tri)):
        x0, x1, y0, y1 = lo_x[t], hi_x[t], lo_y[t], hi_y[t]
        if x1 < x0 or y1 < y0:
            continue
        gx, gy = np.meshgrid(np.arange(x0, x1 + 1) + 0.5,
                             np.arange(y0, y1 + 1) + 0.5)

        ax, bx, cx = x[t]
        ay, by, cy = y[t]
        inv = 1.0 / area[t]
        w0 = ((bx - gx) * (cy - gy) - (cx - gx) * (by - gy)) * inv
        w1 = ((cx - gx) * (ay - gy) - (ax - gx) * (cy - gy)) * inv
        w2 = 1.0 - w0 - w1
        inside = (w0 >= 0) & (w1 >= 0) & (w2 >= 0)
        if not inside.any():
            continue

        za, zb, zc = z[t]
        frag_z = w0 * za + w1 * zb + w2 * zc
        sub_z = zbuffer[y0:y1 + 1, x0:x1 + 1]
        win = inside & (frag_z < sub_z)
        if not win.any():
            continue

        n = (w0[win, None] * corner[t, 0]
             + w1[win, None] * corner[t, 1]
             + w2[win, None] * corner[t, 2])
        length = np.linalg.norm(n, axis=1, keepdims=True)
        n = np.divide(n, length, out=np.zeros_like(n), where=length > 0)

        shade = np.full((int(win.sum()), 3), 0.16, dtype=np.float32) * base
        for direction, tint, strength in lights:
            lambert = np.clip(n @ direction, 0.0, 1.0)[:, None]
            shade = shade + base * tint * lambert * strength
        rim = np.clip(1.0 - np.abs(n @ basis[2]), 0.0, 1.0)[:, None] ** 3
        shade = shade + rim * 0.22

        colour[y0:y1 + 1, x0:x1 + 1][win] = np.clip(shade, 0.0, 1.0)
        sub_z[win] = frag_z[win]

    image = np.clip(colour, 0.0, 1.0) ** (1.0 / 2.2)
    image = image.reshape(size, supersample, size, supersample, 3).mean((1, 3))
    return (image * 255.0 + 0.5).astype(np.uint8)
