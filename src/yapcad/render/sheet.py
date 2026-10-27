"""Annotation and layout: what separates an engineering view from a picture.

An agent reading one of these images cannot measure a pixel, hover a tooltip
or orbit the model.  Everything it needs has to be *in* the frame: how big a
millimetre is, which way the axes point, what the view is, which part is which
colour, and where a section was cut.  This module draws all of that.
"""

import numpy as np

from yapcad.render.palette import get_theme
from yapcad.render.raster import (
    legend,
    line,
    new_image,
    paste,
    plate,
    rect,
    scale_bar,
    text,
    text_size,
)

#: Side margin for the annotation, in pixels.
PAD = 12


def axis_triad(img, camera, theme, anchor=None, length=26):
    """A labelled XYZ triad, so the reader knows which way the model faces.

    An axis pointing nearly along the view direction has nowhere to go on
    screen, so its label is parked beside the origin rather than drawn on top
    of it.
    """
    theme = get_theme(theme)
    _, right, up = camera.basis
    height, width, _ = img.shape
    ox, oy = (width - length - 40, height - length - 28) if anchor is None \
        else anchor
    plate(img, ox - length - 14, oy - length - 14,
          2 * length + 28, 2 * length + 28, theme.surface, 0.7)
    rect(img, ox - 1, oy - 1, 3, 3, theme.ink_soft)
    parked = 0
    for name, vec in (("X", (1.0, 0.0, 0.0)),
                      ("Y", (0.0, 1.0, 0.0)),
                      ("Z", (0.0, 0.0, 1.0))):
        u = float(np.dot(vec, right))
        v = float(np.dot(vec, up))
        tw, th = text_size(name, 2)
        if np.hypot(u, v) < 0.12:        # end-on: the axis points at the eye
            text(img, name, ox - tw - 6 - parked * (tw + 4),
                 oy + 7, theme.ink_soft, 2)
            parked += 1
            continue
        tip_x, tip_y = ox + u * length, oy - v * length
        line(img, ox, oy, tip_x, tip_y, theme.ink_soft)
        text(img, name, int(round(tip_x + u * 7 - tw / 2)),
             int(round(tip_y - v * 7 - th / 2)), theme.ink, 2)


def annotate(img, camera, theme="light", title=None, subtitle=None,
             entries=(), scale=True, triad=True):
    """Draw the standard furniture onto a rendered view.

    :param scale: draw a scale bar.  Only meaningful under an orthographic
        projection, where one millimetre is the same number of pixels
        everywhere; a perspective view is labelled as such instead.
    """
    theme = get_theme(theme)
    lines = [(t, ink) for t, ink in ((title, theme.ink),
                                     (subtitle, theme.ink_soft)) if t]
    if lines:
        widest = max(text_size(t, 2)[0] for t, _ in lines)
        row = text_size("X", 2)[1] + 5
        plate(img, 0, 0, widest + 2 * PAD, PAD + row * len(lines) + 4,
              theme.surface)
        y = PAD
        for t, ink in lines:
            text(img, t, PAD, y, ink, 2)
            y += row
    if scale and camera.is_orthographic:
        scale_bar(img, camera.mm_per_px, theme.ink_soft, theme.ink,
                  backing=theme.surface)
    elif scale:
        note = "PERSPECTIVE - NOT TO SCALE"
        plate(img, PAD - 6, img.shape[0] - 28,
              text_size(note, 2)[0] + 12, 26, theme.surface)
        text(img, note, PAD, img.shape[0] - 22, theme.ink_soft, 2)
    if triad:
        axis_triad(img, camera, theme)
    if entries:
        legend(img, entries, theme.ink, backing=theme.surface)
    return img


def grid(tiles, columns=2, gutter=10, surface=(1.0, 1.0, 1.0)):
    """Lay uint8 tiles out in a grid on a common surface.

    Tiles may differ in size; each cell is the largest tile's size and
    smaller tiles sit at its top-left.
    """
    tiles = list(tiles)
    if not tiles:
        raise ValueError("grid needs at least one tile")
    columns = max(int(columns), 1)
    rows = -(-len(tiles) // columns)
    cell_h = max(t.shape[0] for t in tiles)
    cell_w = max(t.shape[1] for t in tiles)
    sheet = new_image(columns * cell_w + (columns + 1) * gutter,
                             rows * cell_h + (rows + 1) * gutter, surface)
    for index, tile in enumerate(tiles):
        row, col = divmod(index, columns)
        paste(sheet, tile,
                     gutter + col * (cell_w + gutter),
                     gutter + row * (cell_h + gutter))
    return sheet
