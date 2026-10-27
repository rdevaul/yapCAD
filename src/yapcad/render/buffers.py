"""The seam between "how was this part drawn" and "what does the sheet look
like".

A backend's whole job is to fill a :class:`PartBuffer` for one part: a depth
in millimetres along the view direction, a surface normal, and nothing else.
:func:`composite` merges those into a :class:`GBuffer` by nearest depth and
records *which* part won each pixel.  Shading, line work and annotation then
read only the :class:`GBuffer`.

That split is what lets a second backend join later without reworking
anything downstream.  The sphere tracer in :mod:`yapcad.render.trace` fills a
:class:`PartBuffer` from a field; a triangle rasteriser would fill the same
buffer from a tessellated BREP or an imported mesh, and the silhouette between
a field part and a mesh part would fall out of the composite for free, because
the edge pass never asks where a pixel came from.

Nearest-wins is decided with a tolerance rather than exactly.  Datum-mated
parts share a face, and an exact ``argmin`` over two depths that differ only
by the float noise a transform chain leaves behind decides each pixel
independently, which reads as speckle down every flush joint.  A part has to
be nearer by more than a fraction of a pixel to take one, so a shared face
goes wholly to the earlier part and reads as the clean part-boundary line the
edge pass draws there.  Note that quantising depth into bins does *not* fix
this -- two depths either side of a bin edge still split -- which is why the
comparison is a tolerance.
"""

from dataclasses import dataclass, field

import numpy as np


@dataclass
class PartBuffer:
    """What one backend produces for one part.

    :param depth: ``(h, w)`` float distance along the view direction;
        ``inf`` where the part was missed.
    :param normal: ``(h, w, 3)`` unit surface normals, zero where missed.
    :param name: what the part is called, for the label and the legend.
    :param group: the palette grouping key -- a material, a subassembly -- or
        ``""`` for the neutral material tone.
    :param stats: whatever the backend wants to report; see
        :func:`yapcad.render.trace.render_part`.
    """

    depth: np.ndarray
    normal: np.ndarray
    name: str = ""
    group: str = ""
    stats: dict = field(default_factory=dict)

    def __post_init__(self):
        if self.depth.shape + (3,) != self.normal.shape:
            raise ValueError(
                f"depth {self.depth.shape} and normal {self.normal.shape} "
                "disagree about the image size")

    @property
    def shape(self):
        """``(height, width)``."""
        return self.depth.shape

    @property
    def hit(self):
        """Boolean mask of the pixels this part covers."""
        return np.isfinite(self.depth)

    @property
    def coverage(self):
        """Fraction of the frame this part covers."""
        return float(self.hit.mean())


@dataclass
class GBuffer:
    """The composited scene: nearest depth, its normal, and whose it is.

    :param part: ``(h, w)`` index into :attr:`names`; ``-1`` for background.
    """

    depth: np.ndarray
    normal: np.ndarray
    part: np.ndarray
    names: tuple = ()
    groups: tuple = ()

    @property
    def shape(self):
        return self.depth.shape

    @property
    def hit(self):
        return self.part >= 0

    def depth_range(self):
        """``(near, far)`` over the covered pixels, or ``None`` if empty."""
        hit = self.hit
        if not hit.any():
            return None
        values = self.depth[hit]
        return float(values.min()), float(values.max())


def composite(parts, tie=0.0):
    """Merge :class:`PartBuffer` objects into one :class:`GBuffer`.

    :param parts: the buffers, in a stable order.  A pixel the tolerance
        cannot separate goes to the earliest, which is why the order should
        be stable across renders.
    :param tie: how much nearer, in millimetres, a part must be to take a
        pixel from one already there.  A quarter of a pixel is a good choice;
        zero compares raw depths and will speckle on coincident faces.
    """
    parts = list(parts)
    if not parts:
        raise ValueError("composite needs at least one part buffer")
    shapes = {p.shape for p in parts}
    if len(shapes) != 1:
        raise ValueError(f"part buffers disagree about the image size: {shapes}")
    shape = parts[0].shape

    # Folded one part at a time rather than stacked: an assembly render holds
    # one buffer per part, and stacking their normals would cost hundreds of
    # megabytes for a scene this code is meant to handle.
    best = np.full(shape, np.inf)
    winner = np.full(shape, -1, dtype=np.int32)
    margin = max(float(tie), 0.0)
    for index, part in enumerate(parts):
        # Nearer by more than the tolerance, so a shared face does not
        # speckle; `best` still tracks the true nearest depth so the
        # tolerance cannot accumulate over a long scene.
        winner[part.depth < best - margin] = index
        np.minimum(best, part.depth, out=best)

    depth = np.full(shape, np.inf)
    normal = np.zeros(shape + (3,))
    for index, part in enumerate(parts):
        mine = winner == index
        if mine.any():
            depth[mine] = part.depth[mine]
            normal[mine] = part.normal[mine]
    return GBuffer(depth=depth, normal=normal, part=winner,
                   names=tuple(p.name for p in parts),
                   groups=tuple(p.group for p in parts))
