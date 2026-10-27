"""Headless engineering renders from signed distance fields.

``docs/SDF-DESIGN.md`` section 6 plans two rendering tiers: mesh the field and
hand the triangles to an existing viewer, or emit GLSL and sphere-trace it in
a browser.  Both put a human in front of a screen.  This package is the tier
in between, for the reader who is a *program*: no window, no GPU, no mesh --
sphere-trace the field in numpy and write a PNG, in seconds.

It carries none of the shader emitter's risk of drifting away from the numpy
backend, because it *is* the numpy backend.  ``evaluate`` and ``normal`` from
:mod:`yapcad.sdf` are the whole of it, so ``smooth_union`` cannot come to mean
two different things.

::

    from yapcad import render

    render.image(part, view="iso", path="part.png")
    render.four_up(part, path="part_4up.png")
    render.section_view(part, axis="z", where=4.0, path="cut.png")

From the command line, against a DSL design::

    python -m yapcad.render design.dsl PART_COMMAND --view 4up -o part.png

The architecture is worth knowing before extending it.  Each part is traced
into a :class:`~yapcad.render.buffers.PartBuffer` -- a depth and a normal, and
nothing about the rest of the scene -- and :func:`~yapcad.render.buffers.composite`
merges those by nearest depth.  Shading and line work read only the merged
buffer.  So a second backend that rasterises triangles from a tessellated BREP
or an imported mesh can fill the same buffer and composite against field
parts, and the silhouette where the two meet falls out for free.  That is the
seam; keep it.

Limits worth knowing: the geometry has to be SDF-authored, a ray that grazes a
feature thinner than the hit tolerance can miss it, and an image is for
*seeing* rather than measuring -- numbers belong to :mod:`yapcad.rom` and the
clearance APIs.
"""

from yapcad.render.api import (
    Part,
    as_parts,
    four_up,
    image,
    parts_from_dsl,
    render,
    scene_bounds,
    section_view,
)
from yapcad.render.buffers import GBuffer, PartBuffer, composite
from yapcad.render.camera import (
    ORTHOGRAPHIC,
    PERSPECTIVE,
    VIEWS,
    Camera,
)
from yapcad.render.raster import new_image, to_image, write_png
from yapcad.render.palette import DARK, LIGHT, Theme, get_theme
from yapcad.render.section import Section, section, section_image
from yapcad.render.shade import edges, legend_entries, shade
from yapcad.render.sheet import annotate, grid
from yapcad.render.trace import march, render_part

__all__ = [
    "Camera",
    "VIEWS",
    "ORTHOGRAPHIC",
    "PERSPECTIVE",
    "Part",
    "PartBuffer",
    "GBuffer",
    "Section",
    "Theme",
    "LIGHT",
    "DARK",
    "as_parts",
    "scene_bounds",
    "render",
    "render_part",
    "march",
    "composite",
    "shade",
    "edges",
    "legend_entries",
    "section",
    "section_image",
    "section_view",
    "image",
    "four_up",
    "parts_from_dsl",
    "annotate",
    "grid",
    "new_image",
    "to_image",
    "write_png",
    "get_theme",
]
