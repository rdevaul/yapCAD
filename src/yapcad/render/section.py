"""A planar section, which is not ray-marched at all.

Cutting a field is a single evaluation of the tree over a plane of points.
There is no marching, no meshing and no boolean: the slice *is* the field
restricted to a plane, so a section costs one ``evaluate`` call and comes out
cheaper than any shaded view of the same part.

What you get is more than a mesh renderer can give you.  Outside the material
the value is the clearance to the nearest surface, so the exterior of the
slice is a clearance map; inside, its magnitude is the local half-thickness.
The zero contour is the cut outline.
"""

import math
import warnings
from dataclasses import dataclass

import numpy as np

from yapcad.sdf.evaluate import DEFAULT_BACKEND, evaluate

from yapcad.render.palette import get_theme
from yapcad.render.raster import round_number, to_image

#: Which world axes become screen right and up for a cut normal to each axis.
#: Chosen to match the matching orthographic view, so a section and a plan
#: read the same way round.
PLANES = {
    "x": (1, 2, "+Y", "+Z"),
    "y": (0, 2, "+X", "+Z"),
    "z": (0, 1, "+X", "+Y"),
}


@dataclass
class Section:
    """A sampled slice through a field.

    :param distance: ``(h, w)`` signed distance; negative inside material.
    :param axis: ``"x"``, ``"y"`` or ``"z"`` -- the cut plane's normal.
    :param where: the plane's coordinate along ``axis``, in millimetres.
    :param mm_per_px: the slice scale.
    :param span: ``((u0, v0), (u1, v1))`` in-plane extent in millimetres.
    :param labels: the screen-right and screen-up axis names.
    """

    distance: np.ndarray
    axis: str
    where: float
    mm_per_px: float
    span: tuple
    labels: tuple

    @property
    def shape(self):
        return self.distance.shape

    @property
    def area_mm2(self):
        """Cut area: the material the plane passes through."""
        return float((self.distance < 0.0).sum()) * self.mm_per_px ** 2

    def is_empty(self):
        """True if the plane misses the model entirely."""
        return not bool((self.distance < 0.0).any())


def section(node, axis="z", where=None, bounds=None, pixels=(640, 480),
            margin=1.08, backend=DEFAULT_BACKEND, extent=None):
    """Sample ``node`` on a plane.

    :param axis: the cut plane's normal.
    :param where: the plane's coordinate; defaults to the middle of the
        model along ``axis``.
    :param extent: force the slice height in millimetres, so a section can
        be put on a sheet at the same scale as the views beside it.
    """
    axis = str(axis).lower()
    if axis not in PLANES:
        raise ValueError(f"axis must be one of x, y, z; got {axis!r}")
    box = node.bounds if bounds is None else bounds
    lo = np.asarray(box[0], dtype=float)
    hi = np.asarray(box[1], dtype=float)
    if not (np.all(np.isfinite(lo)) and np.all(np.isfinite(hi))):
        raise ValueError(
            "cannot section an unbounded field; pass explicit finite bounds")
    normal_axis = "xyz".index(axis)
    if where is None:
        where = 0.5 * (lo[normal_axis] + hi[normal_axis])
    where = float(where)
    if not lo[normal_axis] <= where <= hi[normal_axis]:
        warnings.warn(
            f"the {axis} = {where:g} plane lies outside the model's "
            f"{axis} range [{lo[normal_axis]:g}, {hi[normal_axis]:g}]; the "
            "section will be empty",
            RuntimeWarning, stacklevel=2)

    u_axis, v_axis, u_name, v_name = PLANES[axis]
    width, height = int(pixels[0]), int(pixels[1])
    aspect = width / height
    centre_u = 0.5 * (lo[u_axis] + hi[u_axis])
    centre_v = 0.5 * (lo[v_axis] + hi[v_axis])
    reach_u = 0.5 * (hi[u_axis] - lo[u_axis])
    reach_v = 0.5 * (hi[v_axis] - lo[v_axis])
    if extent is None:
        extent = 2.0 * max(reach_v, reach_u / aspect) * float(margin)
        if extent <= 0.0:
            extent = 1.0
    extent = float(extent)
    if not (extent > 0.0) or not math.isfinite(extent):
        raise ValueError(f"extent must be positive, got {extent}")
    mm_per_px = extent / height
    half_v = 0.5 * extent
    half_u = half_v * aspect

    us = centre_u + np.linspace(-half_u, half_u, width, endpoint=False) \
        + 0.5 * mm_per_px
    vs = centre_v - (np.linspace(-half_v, half_v, height, endpoint=False)
                     + 0.5 * mm_per_px)
    grid_u, grid_v = np.meshgrid(us, vs, indexing="xy")
    points = np.empty((grid_u.size, 3))
    points[:, normal_axis] = where
    points[:, u_axis] = grid_u.ravel()
    points[:, v_axis] = grid_v.ravel()

    values = evaluate(node, points, backend=backend)
    return Section(distance=np.asarray(values).reshape(height, width),
                   axis=axis, where=where, mm_per_px=mm_per_px,
                   span=((centre_u - half_u, centre_v - half_v),
                         (centre_u + half_u, centre_v + half_v)),
                   labels=(u_name, v_name))


def section_colours(slice_, theme="light", limit=None):
    """Colour a :class:`Section` into a float ``(h, w, 3)`` image.

    Material is flat; the clearance outside it is a one-hue sequential ramp,
    darkest at the wall.  Returns ``(rgb, limit)`` so the caller can label
    the key with the distance the ramp tops out at.
    """
    theme = get_theme(theme)
    distance = slice_.distance
    outside = distance > 0.0
    if limit is None:
        if outside.any():
            limit = round_number(max(
                float(np.percentile(distance[outside], 90)),
                4.0 * slice_.mm_per_px))
        else:
            limit = round_number(max(4.0 * slice_.mm_per_px, 1.0))
    limit = float(limit)

    steps = np.asarray(theme.ramp)
    fraction = np.clip(distance / limit, 0.0, 1.0)
    index = np.rint((1.0 - fraction) * (len(steps) - 1)).astype(int)
    rgb = steps[index]
    rgb[~outside] = theme.material
    rgb[np.abs(distance) < 0.75 * slice_.mm_per_px] = theme.line
    return rgb, limit


def section_image(slice_, theme="light", limit=None):
    """A uint8 image of a :class:`Section`, without annotation."""
    rgb, _ = section_colours(slice_, theme=theme, limit=limit)
    return to_image(rgb)
