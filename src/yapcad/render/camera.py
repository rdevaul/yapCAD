"""The camera, and the one thing both render backends must agree about.

A :class:`Camera` turns a pixel grid into rays.  Every backend -- the sphere
tracer here, a triangle rasteriser later -- consumes the *same* camera object,
so the two cannot disagree about where the eye is or how millimetres map to
pixels.  That is the same discipline ``docs/SDF-DESIGN.md`` section 6 asks of
the shader emitter, applied one tier down.

Orthographic is the default because it is the projection engineering drawings
are read in: a millimetre is the same number of pixels everywhere, which is
what makes :func:`yapcad.render.raster.scale_bar` honest.  Perspective is
available for a presentation view, and ``mm_per_px`` is then only nominal --
the value at the target plane.

Screen axes follow the drawing conventions: the front view puts ``+X`` right
and ``+Z`` up; the top view puts ``+X`` right and ``+Y`` up.
"""

import math
from dataclasses import dataclass

import numpy as np

#: ``name -> (view direction, up hint)``.  The direction points *from the eye
#: toward the model*, so the ``front`` view looks along ``+Y`` from ``-Y``.
VIEWS = {
    "iso": ((-1.0, -1.0, -1.0), (0.0, 0.0, 1.0)),
    "front": ((0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
    "back": ((0.0, -1.0, 0.0), (0.0, 0.0, 1.0)),
    "left": ((1.0, 0.0, 0.0), (0.0, 0.0, 1.0)),
    "right": ((-1.0, 0.0, 0.0), (0.0, 0.0, 1.0)),
    "top": ((0.0, 0.0, -1.0), (0.0, 1.0, 0.0)),
    "bottom": ((0.0, 0.0, 1.0), (0.0, 1.0, 0.0)),
}

ORTHOGRAPHIC = "orthographic"
PERSPECTIVE = "perspective"


def _unit(v, what="vector"):
    a = np.asarray(v, dtype=float).reshape(3)
    n = float(np.linalg.norm(a))
    if not (n > 0.0) or not math.isfinite(n):
        raise ValueError(f"{what} must have finite non-zero length, got {v!r}")
    return a / n


def _basis(direction, up_hint):
    """Right-handed screen basis: ``right x up == -direction``."""
    d = _unit(direction, "view direction")
    up = _unit(up_hint, "up hint")
    if abs(float(np.dot(up, d))) > 0.999:
        up = np.array([0.0, 1.0, 0.0]) if abs(d[1]) < 0.9 \
            else np.array([1.0, 0.0, 0.0])
    right = np.cross(d, up)
    right /= np.linalg.norm(right)
    up = np.cross(right, d)
    return d, right, up


@dataclass(frozen=True)
class Camera:
    """Where the view is taken from, and how big a pixel is.

    :param target: the point the view is centred on.
    :param direction: unit vector from the eye toward ``target``.
    :param up: unit screen-up vector, perpendicular to ``direction``.
    :param extent: the view height in millimetres at the target plane.
    :param pixels: ``(width, height)`` in pixels.
    :param distance: eye-to-target distance.  For an orthographic view this
        only says where the rays start, and must put them clear of the model.
    :param far: how far a ray may travel before it counts as a miss.
    :param projection: ``"orthographic"`` or ``"perspective"``.
    :param fov: vertical field of view in degrees, perspective only.
    """

    target: tuple
    direction: tuple
    up: tuple
    extent: float
    pixels: tuple = (640, 480)
    distance: float = 1000.0
    far: float = 4000.0
    projection: str = ORTHOGRAPHIC
    fov: float = 32.0

    def __post_init__(self):
        if self.projection not in (ORTHOGRAPHIC, PERSPECTIVE):
            raise ValueError(
                f"projection must be {ORTHOGRAPHIC!r} or {PERSPECTIVE!r}, "
                f"got {self.projection!r}")
        if not (self.extent > 0.0) or not math.isfinite(self.extent):
            raise ValueError(f"extent must be positive, got {self.extent}")
        if not math.isfinite(self.far) or self.far <= 0.0:
            raise ValueError(
                f"far must be a positive finite distance, got {self.far}; an "
                "unbounded field needs explicit bounds so the march can stop")
        w, h = self.pixels
        if int(w) < 2 or int(h) < 2:
            raise ValueError(f"pixels must be at least 2x2, got {self.pixels}")

    @property
    def width(self):
        """Pixel width."""
        return int(self.pixels[0])

    @property
    def height(self):
        """Pixel height."""
        return int(self.pixels[1])

    @property
    def aspect(self):
        return self.width / self.height

    @property
    def mm_per_px(self):
        """Millimetres per pixel.  Exact for orthographic; at the target
        plane for perspective."""
        return self.extent / self.height

    @property
    def is_orthographic(self):
        return self.projection == ORTHOGRAPHIC

    @property
    def basis(self):
        """``(direction, right, up)`` as an orthonormal right-handed set.

        The stored ``up`` is treated as a hint and re-orthogonalised, so a
        hand-built camera cannot produce a skewed image.
        """
        d = _unit(self.direction, "view direction")
        right = np.cross(d, _unit(self.up, "up vector"))
        right /= np.linalg.norm(right)
        return d, right, np.cross(right, d)

    @property
    def eye(self):
        """The eye position."""
        return (np.asarray(self.target, dtype=float)
                - np.asarray(self.direction, dtype=float) * self.distance)

    def rays(self):
        """Ray origins and directions, one per pixel, row-major from the top.

        :return: ``(origins, directions)``, both ``(width * height, 3)``.
        """
        d, right, up = self.basis
        half_h = 0.5 * self.extent
        half_w = half_h * self.aspect
        us = (np.arange(self.width) + 0.5) / self.width * 2.0 - 1.0
        vs = 1.0 - (np.arange(self.height) + 0.5) / self.height * 2.0
        u, v = np.meshgrid(us * half_w, vs * half_h, indexing="xy")
        offset = (u.ravel()[:, None] * right + v.ravel()[:, None] * up)
        if self.is_orthographic:
            origins = self.eye + offset
            dirs = np.broadcast_to(d, origins.shape).copy()
            return origins, dirs
        dirs = d * self.distance + offset
        dirs /= np.linalg.norm(dirs, axis=1)[:, None]
        origins = np.broadcast_to(self.eye, dirs.shape).copy()
        return origins, dirs

    def project(self, points):
        """World points to ``(column, row)`` pixel coordinates.

        Used to place the axis triad, not to rasterise anything.
        """
        pts = np.atleast_2d(np.asarray(points, dtype=float))
        d, right, up = self.basis
        rel = pts - np.asarray(self.target, dtype=float)
        u = rel @ right
        v = rel @ up
        if not self.is_orthographic:
            depth = rel @ d + self.distance
            scale = self.distance / np.where(np.abs(depth) < 1e-9, 1e-9, depth)
            u, v = u * scale, v * scale
        half_h = 0.5 * self.extent
        half_w = half_h * self.aspect
        col = (u / half_w * 0.5 + 0.5) * self.width - 0.5
        row = (0.5 - v / half_h * 0.5) * self.height - 0.5
        return np.stack([col, row], axis=1)

    @classmethod
    def fit(cls, bounds, view="iso", pixels=(640, 480), margin=1.08,
            projection=ORTHOGRAPHIC, fov=32.0, direction=None, up=None,
            extent=None):
        """Frame an axis-aligned box.

        The eight corners are projected onto the screen axes and the extent
        taken from that, so the model fills the frame rather than a sphere
        around it.

        :param bounds: ``((x0, y0, z0), (x1, y1, z1))``; must be finite.
        :param view: a key of :data:`VIEWS`, ignored if ``direction`` is given.
        :param margin: fraction of slack around the model.
        :param extent: force the view height in millimetres instead of
            fitting it.  Several views of one model are only comparable if
            they share a scale, and each view's natural fit differs -- an
            isometric projection of a box reaches further across the screen
            than a face-on one -- so a multi-view sheet fits every view,
            takes the largest, and passes it back here.
        """
        lo = np.asarray(bounds[0], dtype=float)
        hi = np.asarray(bounds[1], dtype=float)
        if not (np.all(np.isfinite(lo)) and np.all(np.isfinite(hi))):
            raise ValueError(
                "cannot frame an unbounded field; pass explicit finite bounds")
        if np.any(hi < lo):
            raise ValueError(f"bounds are inverted: {bounds!r}")
        if direction is None:
            try:
                direction, hint = VIEWS[str(view)]
            except KeyError:
                raise ValueError(
                    f"unknown view {view!r}; expected one of "
                    f"{', '.join(VIEWS)}") from None
        else:
            hint = (0.0, 0.0, 1.0)
        d, right, up_vec = _basis(direction, up if up is not None else hint)

        centre = 0.5 * (lo + hi)
        half = 0.5 * (hi - lo)
        # Support of the box along each screen axis.
        reach_u = float(np.abs(right) @ half)
        reach_v = float(np.abs(up_vec) @ half)
        reach_d = float(np.abs(d) @ half)
        w, h = int(pixels[0]), int(pixels[1])
        aspect = w / h
        forced = extent is not None
        if forced:
            extent = float(extent)
            if not (extent > 0.0) or not math.isfinite(extent):
                raise ValueError(f"extent must be positive, got {extent}")
        else:
            extent = 2.0 * max(reach_v, reach_u / aspect) * float(margin)
            if extent <= 0.0:                  # a degenerate (flat) box
                extent = max(2.0 * reach_d * float(margin), 1.0)
        depth = max(reach_d, 1e-6)
        if projection == PERSPECTIVE:
            distance = depth + 0.5 * extent / math.tan(math.radians(fov) * 0.5)
            if not forced:
                # At that distance the frame is `extent` tall only at the
                # target plane, so re-derive what the fov actually gives.
                extent = 2.0 * distance * math.tan(math.radians(fov) * 0.5)
        else:
            distance = depth * 1.05 + 0.01 * extent
        return cls(target=tuple(centre), direction=tuple(d), up=tuple(up_vec),
                   extent=extent, pixels=(w, h), distance=distance,
                   far=distance + 2.0 * depth + extent,
                   projection=projection, fov=fov)
