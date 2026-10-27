"""Turning a G-buffer into a readable picture.

The line work matters more than the lighting.  A shaded surface tells a reader
roughly which way a face points; a drawn crease tells them exactly where one
face stops and the next begins, and that is what makes a chamfer, a gusset or
a bolt hole legible at 512 pixels.

Every edge here comes from the composited buffer rather than from any one
part, which is the point of :mod:`yapcad.render.buffers`: the silhouette
between two parts is found the same way as the silhouette between a part and
the background, so it does not matter which backend drew either side.
"""

import numpy as np

from yapcad.render.palette import assign_groups, get_theme

#: Normal change, summed over the three components, that counts as a crease.
#: About 35 degrees, which keeps a 45-degree chamfer and loses tessellation
#: noise.
CREASE = 0.6

#: Depth change, in pixels' worth of millimetres, that counts as a step.
STEP = 4.0


def edges(gbuf, mm_per_px, crease=CREASE, step=STEP):
    """Boolean mask of the lines to draw.

    Four kinds, found together: the silhouette against the background, the
    boundary between two parts, a depth step within one part, and a crease
    where the normal turns sharply.
    """
    hit = gbuf.hit
    depth = np.where(hit, gbuf.depth, 0.0)
    mask = np.zeros(hit.shape, dtype=bool)
    for axis in (0, 1):
        lead = [slice(None), slice(None)]
        trail = [slice(None), slice(None)]
        lead[axis] = slice(0, -1)
        trail[axis] = slice(1, None)
        lead, trail = tuple(lead), tuple(trail)
        both = hit[lead] & hit[trail]
        turn = np.abs(gbuf.normal[lead] - gbuf.normal[trail]).sum(axis=-1)
        jump = np.abs(depth[lead] - depth[trail])
        inner = both & ((turn > crease)
                        | (jump > step * mm_per_px)
                        | (gbuf.part[lead] != gbuf.part[trail]))
        silhouette = hit[lead] ^ hit[trail]
        mask[lead] |= inner | (silhouette & hit[lead])
        mask[trail] |= silhouette & hit[trail]
    return mask


def group_slots(groups):
    """The palette slot each part's group key maps to, in part order."""
    assigned = assign_groups(groups)
    return [assigned.get(key, -1) for key in groups]


def legend_entries(gbuf, theme):
    """``[(label, rgb)]`` for the distinct groups in the scene.

    Returns nothing when every part shares the default material: a legend of
    one row says less than the title already does.
    """
    theme = get_theme(theme)
    assigned = assign_groups(gbuf.groups)
    entries = []
    folded = False
    for key, slot in assigned.items():
        if not key:
            continue
        if slot >= len(theme.series):
            folded = True
            continue
        entries.append((key, theme.group_colour(slot)))
    if folded:
        entries.append(("OTHER", theme.other))
    return entries


def base_colours(gbuf, theme):
    """``(h, w, 3)`` of each pixel's unlit part colour."""
    theme = get_theme(theme)
    slots = group_slots(gbuf.groups)
    colours = np.zeros(gbuf.shape + (3,))
    for index, slot in enumerate(slots):
        where = gbuf.part == index
        if where.any():
            colours[where] = theme.group_colour(slot)
    return colours


def shade(gbuf, camera, theme="light", ambient=0.20, key=0.62, fill=0.18,
          draw_edges=True):
    """Shade a G-buffer into a float ``(h, w, 3)`` image in 0..1.

    The lights are fixed in *camera* space, so a part looks lit the same way
    from every named view and two views of the same feature can be compared.
    """
    theme = get_theme(theme)
    direction, right, up = camera.basis
    key_dir = -0.90 * direction + 0.45 * up - 0.35 * right
    key_dir /= np.linalg.norm(key_dir)
    fill_dir = -0.40 * direction + 0.70 * right - 0.20 * up
    fill_dir /= np.linalg.norm(fill_dir)

    hit = gbuf.hit
    normals = gbuf.normal
    lambert = np.clip(normals @ key_dir, 0.0, 1.0)
    bounce = np.clip(normals @ fill_dir, 0.0, 1.0)
    tone = np.clip(ambient + key * lambert + fill * bounce, 0.0, 1.0)

    rgb = np.empty(gbuf.shape + (3,))
    rgb[:, :] = theme.surface
    lit = tone[..., None] * base_colours(gbuf, theme)
    rgb[hit] = lit[hit]
    if draw_edges:
        rgb[edges(gbuf, camera.mm_per_px)] = theme.line
    return rgb
