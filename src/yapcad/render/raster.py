"""The raster: an RGB array, a PNG writer, and the annotation a drawing needs.

Nothing here knows about geometry.  A render becomes an *engineering* image
rather than a picture because of what gets drawn on top of it -- a scale bar,
an axis triad, the bounding box, a legend -- and that is all this module does.

The PNG writer is deliberately built on :mod:`zlib` from the standard library.
yapCAD has no imaging dependency and this does not add one.
"""

import struct
import zlib

import numpy as np

#: A 5x7 bitmap face.  Drawings are lettered in capitals, so :func:`text`
#: upper-cases its input and the face carries no lowercase.
GLYPHS = {
    " ": ("00000", "00000", "00000", "00000", "00000", "00000", "00000"),
    "A": ("01110", "10001", "10001", "11111", "10001", "10001", "10001"),
    "B": ("11110", "10001", "10001", "11110", "10001", "10001", "11110"),
    "C": ("01111", "10000", "10000", "10000", "10000", "10000", "01111"),
    "D": ("11110", "10001", "10001", "10001", "10001", "10001", "11110"),
    "E": ("11111", "10000", "10000", "11110", "10000", "10000", "11111"),
    "F": ("11111", "10000", "10000", "11110", "10000", "10000", "10000"),
    "G": ("01110", "10001", "10000", "10111", "10001", "10001", "01110"),
    "H": ("10001", "10001", "10001", "11111", "10001", "10001", "10001"),
    "I": ("11111", "00100", "00100", "00100", "00100", "00100", "11111"),
    "J": ("00111", "00010", "00010", "00010", "00010", "10010", "01100"),
    "K": ("10001", "10010", "10100", "11000", "10100", "10010", "10001"),
    "L": ("10000", "10000", "10000", "10000", "10000", "10000", "11111"),
    "M": ("10001", "11011", "10101", "10101", "10001", "10001", "10001"),
    "N": ("10001", "11001", "10101", "10011", "10001", "10001", "10001"),
    "O": ("01110", "10001", "10001", "10001", "10001", "10001", "01110"),
    "P": ("11110", "10001", "10001", "11110", "10000", "10000", "10000"),
    "Q": ("01110", "10001", "10001", "10001", "10101", "10011", "01101"),
    "R": ("11110", "10001", "10001", "11110", "10100", "10010", "10001"),
    "S": ("01111", "10000", "10000", "01110", "00001", "00001", "11110"),
    "T": ("11111", "00100", "00100", "00100", "00100", "00100", "00100"),
    "U": ("10001", "10001", "10001", "10001", "10001", "10001", "01110"),
    "V": ("10001", "10001", "10001", "10001", "10001", "01010", "00100"),
    "W": ("10001", "10001", "10001", "10101", "10101", "11011", "10001"),
    "X": ("10001", "10001", "01010", "00100", "01010", "10001", "10001"),
    "Y": ("10001", "10001", "01010", "00100", "00100", "00100", "00100"),
    "Z": ("11111", "00001", "00010", "00100", "01000", "10000", "11111"),
    "0": ("01110", "10001", "10011", "10101", "11001", "10001", "01110"),
    "1": ("00100", "01100", "00100", "00100", "00100", "00100", "01110"),
    "2": ("01110", "10001", "00001", "00010", "00100", "01000", "11111"),
    "3": ("11110", "00001", "00001", "01110", "00001", "00001", "11110"),
    "4": ("00010", "00110", "01010", "10010", "11111", "00010", "00010"),
    "5": ("11111", "10000", "11110", "00001", "00001", "10001", "01110"),
    "6": ("00110", "01000", "10000", "11110", "10001", "10001", "01110"),
    "7": ("11111", "00001", "00010", "00100", "01000", "01000", "01000"),
    "8": ("01110", "10001", "10001", "01110", "10001", "10001", "01110"),
    "9": ("01110", "10001", "10001", "01111", "00001", "00010", "01100"),
    ".": ("00000", "00000", "00000", "00000", "00000", "01100", "01100"),
    ",": ("00000", "00000", "00000", "00000", "01100", "01100", "00100"),
    "-": ("00000", "00000", "00000", "11111", "00000", "00000", "00000"),
    "+": ("00000", "00100", "00100", "11111", "00100", "00100", "00000"),
    "=": ("00000", "00000", "11111", "00000", "11111", "00000", "00000"),
    ":": ("00000", "01100", "01100", "00000", "01100", "01100", "00000"),
    "/": ("00001", "00010", "00010", "00100", "01000", "01000", "10000"),
    "(": ("00110", "01000", "10000", "10000", "10000", "01000", "00110"),
    ")": ("01100", "00010", "00001", "00001", "00001", "00010", "01100"),
    "[": ("01110", "01000", "01000", "01000", "01000", "01000", "01110"),
    "]": ("01110", "00010", "00010", "00010", "00010", "00010", "01110"),
    "<": ("00010", "00100", "01000", "10000", "01000", "00100", "00010"),
    ">": ("01000", "00100", "00010", "00001", "00010", "00100", "01000"),
    "_": ("00000", "00000", "00000", "00000", "00000", "00000", "11111"),
    "#": ("01010", "01010", "11111", "01010", "11111", "01010", "01010"),
    "%": ("10001", "10010", "00100", "00100", "00100", "01001", "10001"),
    "*": ("00000", "10101", "01110", "11111", "01110", "10101", "00000"),
    "'": ("00100", "00100", "00000", "00000", "00000", "00000", "00000"),
    "?": ("01110", "10001", "00001", "00110", "00100", "00000", "00100"),
    "°": ("01100", "10010", "10010", "01100", "00000", "00000", "00000"),
}

#: Advance per character, in font cells (5 wide plus one of letterspacing).
_ADVANCE = 6
_GLYPH_H = 7


def new_image(width, height, colour=(1.0, 1.0, 1.0)):
    """A ``(height, width, 3)`` uint8 image filled with ``colour``."""
    img = np.empty((int(height), int(width), 3), dtype=np.uint8)
    img[:, :] = _byte(colour)
    return img


def _byte(colour):
    return np.clip(np.asarray(colour, dtype=float) * 255.0 + 0.5,
                   0, 255).astype(np.uint8)


def to_image(rgb):
    """Convert a float ``(h, w, 3)`` array in 0..1 to a uint8 image."""
    return np.clip(np.asarray(rgb, dtype=float) * 255.0 + 0.5,
                   0, 255).astype(np.uint8)


def write_png(path, img):
    """Write a uint8 ``(h, w, 3)`` image as an 8-bit RGB PNG."""
    img = np.ascontiguousarray(np.asarray(img, dtype=np.uint8))
    if img.ndim != 3 or img.shape[2] != 3:
        raise ValueError(f"expected an (h, w, 3) RGB image, got {img.shape}")
    height, width, _ = img.shape
    # Filter type 0 (none) on every scanline: the rows are already small and
    # zlib does the compressing.
    stride = width * 3
    raw = np.zeros((height, stride + 1), dtype=np.uint8)
    raw[:, 1:] = img.reshape(height, stride)
    body = b"".join((
        b"\x89PNG\r\n\x1a\n",
        _chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0)),
        _chunk(b"IDAT", zlib.compress(raw.tobytes(), 6)),
        _chunk(b"IEND", b""),
    ))
    with open(path, "wb") as handle:
        handle.write(body)
    return path


def _chunk(tag, data):
    payload = tag + data
    return (struct.pack(">I", len(data)) + payload
            + struct.pack(">I", zlib.crc32(payload) & 0xFFFFFFFF))


def text_size(s, scale=2):
    """``(width, height)`` in pixels of what :func:`text` would draw."""
    return (_ADVANCE * scale * len(str(s)), _GLYPH_H * scale)


def text(img, s, x, y, colour, scale=2):
    """Draw ``s`` with its top-left corner at ``(x, y)``.  Clips at the edges."""
    rgb = _byte(colour)
    height, width, _ = img.shape
    for ch in str(s).upper():
        glyph = GLYPHS.get(ch, GLYPHS["?"])
        for row, bits in enumerate(glyph):
            y0 = y + row * scale
            if y0 + scale <= 0 or y0 >= height:
                continue
            for col, bit in enumerate(bits):
                if bit != "1":
                    continue
                x0 = x + col * scale
                if x0 + scale <= 0 or x0 >= width:
                    continue
                img[max(y0, 0):y0 + scale, max(x0, 0):x0 + scale] = rgb
        x += _ADVANCE * scale
    return x


def plate(img, x, y, width, height, colour, alpha=0.82):
    """Wash a region toward ``colour`` so text drawn on it stays readable.

    A view that fills its frame would otherwise put the title, the scale bar
    and the legend straight on top of the model.  Washing rather than filling
    keeps the geometry underneath visible.
    """
    h, w, _ = img.shape
    x0, y0 = max(int(x), 0), max(int(y), 0)
    x1, y1 = min(int(x) + int(width), w), min(int(y) + int(height), h)
    if x1 <= x0 or y1 <= y0:
        return
    wash = np.asarray(colour, dtype=float) * 255.0
    region = img[y0:y1, x0:x1].astype(float)
    img[y0:y1, x0:x1] = np.clip(
        alpha * wash + (1.0 - alpha) * region + 0.5, 0, 255).astype(np.uint8)


def rect(img, x, y, width, height, colour):
    """Fill an axis-aligned rectangle, clipped to the image."""
    h, w, _ = img.shape
    x0, y0 = max(int(x), 0), max(int(y), 0)
    x1, y1 = min(int(x) + int(width), w), min(int(y) + int(height), h)
    if x1 > x0 and y1 > y0:
        img[y0:y1, x0:x1] = _byte(colour)


def frame(img, x, y, width, height, colour, weight=1):
    """Draw the outline of a rectangle."""
    rect(img, x, y, width, weight, colour)
    rect(img, x, y + height - weight, width, weight, colour)
    rect(img, x, y, weight, height, colour)
    rect(img, x + width - weight, y, weight, height, colour)


def line(img, x0, y0, x1, y1, colour, weight=1):
    """Draw a line.  Sampled along the longer axis, which is enough for the
    short annotation strokes this module draws."""
    x0, y0, x1, y1 = float(x0), float(y0), float(x1), float(y1)
    steps = int(max(abs(x1 - x0), abs(y1 - y0))) + 1
    ts = np.linspace(0.0, 1.0, steps)
    xs = np.round(x0 + (x1 - x0) * ts).astype(int)
    ys = np.round(y0 + (y1 - y0) * ts).astype(int)
    for px, py in zip(xs, ys):
        rect(img, px, py, weight, weight, colour)


def paste(dst, src, x, y):
    """Copy ``src`` into ``dst`` with its top-left corner at ``(x, y)``."""
    h, w, _ = src.shape
    dst[int(y):int(y) + h, int(x):int(x) + w] = src


def round_number(target):
    """The largest 1/2/5 x 10^n value that does not exceed ``target``.

    Used for scale bars and ramp ticks, so the number a reader has to trust
    is always a round one.
    """
    if not (target > 0.0):
        raise ValueError(f"need a positive target, got {target}")
    decade = 10.0 ** np.floor(np.log10(target))
    for step in (10.0, 5.0, 2.0, 1.0):
        if decade * step <= target:
            return decade * step
    return decade


def scale_bar(img, mm_per_px, colour, ink, x=None, y=None, fraction=0.22,
              backing=None):
    """A bar whose length is a round number of millimetres, with its label.

    :param mm_per_px: the image scale.  Only meaningful for an orthographic
        projection, which is why the caller decides whether to draw one.
    :return: the length drawn, in millimetres.
    """
    h, w, _ = img.shape
    span = round_number(w * fraction * mm_per_px)
    px = max(int(round(span / mm_per_px)), 2)
    x = 14 if x is None else x
    y = h - 22 if y is None else y
    if backing is not None:
        plate(img, x - 7, y - 24, px + 14, 36, backing)
    rect(img, x, y, px, 3, colour)
    rect(img, x, y - 4, 1, 11, colour)
    rect(img, x + px - 1, y - 4, 1, 11, colour)
    text(img, f"{span:g} MM", x, y - 20, ink, 2)
    return span


def legend(img, entries, colour_ink, x=None, y=None, scale=2, swatch=9,
           backing=None):
    """Swatch-and-label rows, one per ``(label, rgb)`` pair.

    Every entry is labelled: colour never carries identity on its own.
    """
    h, w, _ = img.shape
    rows = len(entries)
    if not rows:
        return
    row_h = max(swatch + 4, _GLYPH_H * scale + 4)
    widest = max(text_size(label, scale)[0] for label, _ in entries)
    box_w = swatch + 8 + widest + 10
    x = w - box_w - 12 if x is None else x
    y = 12 if y is None else y
    if backing is not None:
        plate(img, x - 6, y - 5, box_w + 6, rows * row_h + 10, backing)
    for i, (label, rgb) in enumerate(entries):
        top = y + i * row_h
        rect(img, x, top + (row_h - swatch) // 2, swatch, swatch, rgb)
        text(img, label, x + swatch + 6,
             top + (row_h - _GLYPH_H * scale) // 2, colour_ink, scale)


def ramp_key(img, steps, low, high, unit, ink, x=None, y=None,
             width=None, height=10):
    """A horizontal sequential key, labelled at both ends."""
    h, w, _ = img.shape
    width = int(w * 0.24) if width is None else int(width)
    x = w - width - 14 if x is None else int(x)
    y = h - 34 if y is None else int(y)
    for i in range(width):
        slot = min(int(i / max(width - 1, 1) * len(steps)), len(steps) - 1)
        rect(img, x + i, y, 1, height, steps[slot])
    frame(img, x, y, width, height, ink)
    text(img, f"{low:g}", x, y + height + 4, ink, 2)
    label = f"{high:g} {unit}".strip()
    text(img, label, x + width - text_size(label, 2)[0],
         y + height + 4, ink, 2)
