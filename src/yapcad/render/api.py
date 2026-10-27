"""The functions you actually call.

One field or a handful of them in, a PNG out, headlessly, in seconds, with no
meshing and no GPU.  The point is that an agent building a part can look at
what it built in the same turn it built it.
"""

import time
from dataclasses import dataclass

import numpy as np

from yapcad.sdf.evaluate import DEFAULT_BACKEND
from yapcad.sdf.node import Node, bounds_union

from yapcad.render.palette import get_theme
from yapcad.render.raster import (
    plate,
    ramp_key,
    scale_bar,
    text,
    text_size,
    to_image,
    write_png,
)
from yapcad.render.shade import legend_entries, shade
from yapcad.render.sheet import PAD, annotate as draw_furniture, grid
from yapcad.render.buffers import composite
from yapcad.render.camera import ORTHOGRAPHIC, PERSPECTIVE, VIEWS, Camera
from yapcad.render.section import section, section_colours
from yapcad.render.trace import render_part

#: Depth quantisation for the nearest-wins test, in pixels.  See
#: :func:`yapcad.render.buffers.composite` for why this is not zero.
TIE_PIXELS = 0.25


@dataclass(frozen=True)
class Part:
    """One field in a scene.

    :param node: the field.
    :param name: label for the title and the legend.
    :param group: palette grouping key -- a material, a subassembly -- or
        ``""`` for the neutral material tone.
    """

    node: Node
    name: str = ""
    group: str = ""


def as_parts(target):
    """Normalise what a caller might pass into a list of :class:`Part`.

    Accepts a field, an SDF-authored yapCAD solid, a :class:`Part`, a mapping
    of ``{name: field-or-solid}``, or any iterable of those.
    """
    if isinstance(target, Part):
        return [target]
    if isinstance(target, Node):
        return [Part(target)]
    node = _as_node(target)
    if node is not None:
        return [Part(node)]
    if hasattr(target, "items"):
        parts = []
        for name, value in target.items():
            for part in as_parts(value):
                parts.append(Part(part.node, name=str(name),
                                  group=part.group))
        return parts
    if isinstance(target, (list, tuple)):
        parts = []
        for value in target:
            parts.extend(as_parts(value))
        return parts
    raise TypeError(
        "expected a field, an SDF-authored solid, a Part, a mapping of those "
        f"or a list of them; got {type(target).__name__}")


def _as_node(value):
    """Recover a field from a yapCAD solid, or ``None`` if it has none."""
    from yapcad.construction import get_construction
    from yapcad.sdf.node import from_construction, is_sdf_construction
    try:
        record = get_construction(value)
    except Exception:
        return None
    if not is_sdf_construction(record):
        return None
    return from_construction(record)


def scene_bounds(parts):
    """The union of the parts' bounds."""
    box = None
    for part in parts:
        box = part.node.bounds if box is None else bounds_union(
            box, part.node.bounds)
    if box is None:
        raise ValueError("an empty scene has no bounds")
    return box


def _resolve_groups(parts, colour_by):
    """Decide what the palette keys on.

    ``"auto"`` colours by part when a small scene gives no groups of its own,
    because a two-to-eight-part view is exactly where per-part colour helps.
    Past eight the ninth part would have to share a hue, so the neutral
    material tone and the part-boundary lines carry identity instead.
    """
    if colour_by not in ("auto", "part", "group", "none"):
        raise ValueError(
            f"colour_by must be auto, part, group or none; got {colour_by!r}")
    if colour_by == "none":
        return [Part(p.node, p.name, "") for p in parts]
    if colour_by == "group":
        return list(parts)
    named = any(p.group for p in parts)
    if colour_by == "auto" and (named or not 2 <= len(parts) <= 8):
        return list(parts)
    return [Part(p.node, p.name, p.group or p.name or f"PART {i + 1}")
            for i, p in enumerate(parts)]


def render(target, camera, colour_by="auto", backend=DEFAULT_BACKEND,
           **trace_kwargs):
    """Trace every part and composite them into one G-buffer.

    Each part is marched against its own field, so the cost of a scene grows
    with the parts in it rather than with the size of one combined tree.
    """
    parts = _resolve_groups(as_parts(target), colour_by)
    buffers = [render_part(p.node, camera, name=p.name, group=p.group,
                           backend=backend, **trace_kwargs)
               for p in parts]
    return composite(buffers, tie=TIE_PIXELS * camera.mm_per_px)


def image(target, view="iso", pixels=(640, 480), theme="light",
          projection=ORTHOGRAPHIC, title=None, subtitle=None, path=None,
          bounds=None, margin=1.08, fov=32.0, colour_by="auto",
          annotate=True, draw_edges=True, extent=None,
          backend=DEFAULT_BACKEND, **trace_kwargs):
    """Render one view and return a uint8 ``(h, w, 3)`` image.

    :param view: a key of :data:`yapcad.render.camera.VIEWS`, or a direction
        vector.
    :param path: if given, the PNG is written there as well.
    :param colour_by: ``"auto"``, ``"part"``, ``"group"`` or ``"none"``.
    :param extent: force the view height in millimetres, so this view
        shares a scale with others of the same model.
    """
    parts = as_parts(target)
    theme = get_theme(theme)
    box = scene_bounds(parts) if bounds is None else bounds
    if isinstance(view, str):
        camera = Camera.fit(box, view=view, pixels=pixels, margin=margin,
                            projection=projection, fov=fov, extent=extent)
    else:
        camera = Camera.fit(box, direction=view, pixels=pixels, margin=margin,
                            projection=projection, fov=fov, extent=extent)
    gbuf = render(parts, camera, colour_by=colour_by, backend=backend,
                  **trace_kwargs)
    img = to_image(shade(gbuf, camera, theme=theme,
                                       draw_edges=draw_edges))
    if annotate:
        if title is None:
            named = [p.name for p in parts if p.name]
            title = named[0] if len(named) == 1 and len(parts) == 1 \
                else (f"{len(parts)} PARTS" if len(parts) > 1 else "")
        if subtitle is None:
            lo, hi = np.asarray(box[0], float), np.asarray(box[1], float)
            size = hi - lo
            label = view if isinstance(view, str) else "VIEW"
            subtitle = (f"{label.upper()}  "
                        f"{size[0]:.0f} X {size[1]:.0f} X {size[2]:.0f} MM")
        draw_furniture(img, camera, theme, title=title, subtitle=subtitle,
                        entries=legend_entries(gbuf, theme))
    if path:
        write_png(path, img)
    return img


def section_view(target, axis="z", where=None, pixels=(640, 480),
                 theme="light", title=None, path=None, bounds=None,
                 margin=1.08, limit=None, annotate=True, extent=None,
                 backend=DEFAULT_BACKEND):
    """Render a planar section, labelled with where it was cut.

    One evaluation of the field over a plane: no marching, no mesh.  Outside
    the material the colour is the clearance to the nearest surface, so the
    exterior doubles as a clearance map.
    """
    parts = as_parts(target)
    theme = get_theme(theme)
    box = scene_bounds(parts) if bounds is None else bounds
    slices = [section(p.node, axis=axis, where=where, bounds=box,
                      pixels=pixels, margin=margin, backend=backend,
                      extent=extent)
              for p in parts]
    merged = slices[0]
    if len(slices) > 1:
        merged.distance = np.minimum.reduce([s.distance for s in slices])
    rgb, limit = section_colours(merged, theme=theme, limit=limit)
    img = to_image(rgb)
    if annotate:
        if title is None:
            named = [p.name for p in parts if p.name]
            title = (f"SECTION - {named[0]}" if len(named) == 1
                     else f"SECTION - {len(parts)} PARTS")
        plane = (f"{axis.upper()} = {merged.where:g} MM  "
                 f"({merged.labels[0]} RIGHT, {merged.labels[1]} UP)")
        widest = max(text_size(title, 2)[0], text_size(plane, 2)[0])
        plate(img, 0, 0, widest + 2 * PAD, PAD + 2 * 19 + 4, theme.surface)
        text(img, title, PAD, PAD, theme.ink, 2)
        text(img, plane, PAD, PAD + 19, theme.ink_soft, 2)
        scale_bar(img, merged.mm_per_px, theme.ink, theme.ink,
                  backing=theme.surface)
        ramp_key(img, theme.ramp, 0, limit, "MM CLEAR", theme.ink)
    if path:
        write_png(path, img)
    return img


def four_up(target, views=("iso", "front", "top"), section_axis="z",
            section_at=None, pixels=(440, 340), theme="light", gutter=10,
            path=None, bounds=None, colour_by="auto", extent=None,
            annotate=True, backend=DEFAULT_BACKEND, **trace_kwargs):
    """Three orthographic views and a section, on one sheet.

    This is the sheet to reach for first: the three views place a feature in
    space and the section shows what is inside it, which between them answer
    most of what a reader wants to know about a part.

    Every tile is drawn at one scale -- the largest any of the views needs --
    so a feature can be followed from one tile to the next and the scale bars
    agree.  Left to fit themselves, an isometric view and a face-on view of
    the same model come out at different scales, which is exactly the
    comparison a sheet exists to support.
    """
    parts = as_parts(target)
    theme = get_theme(theme)
    box = scene_bounds(parts) if bounds is None else bounds
    lo, hi = np.asarray(box[0], float), np.asarray(box[1], float)
    size = hi - lo
    span = f"{size[0]:.0f} X {size[1]:.0f} X {size[2]:.0f} MM"
    if extent is None:
        extent = max(Camera.fit(box, view=v, pixels=pixels).extent
                     for v in views) if views else None
    tiles = [image(parts, view=v, pixels=pixels, theme=theme, bounds=box,
                   title=v.upper(), subtitle=span, colour_by=colour_by,
                   extent=extent, annotate=annotate, backend=backend,
                   **trace_kwargs)
             for v in views]
    if section_axis:
        tiles.append(section_view(parts, axis=section_axis, where=section_at,
                                  pixels=pixels, theme=theme, bounds=box,
                                  title="SECTION", extent=extent,
                                  annotate=annotate, backend=backend))
    sheet = grid(tiles, columns=2, gutter=gutter, surface=theme.surface)
    if path:
        write_png(path, sheet)
    return sheet


def parts_from_dsl(path, command, params=None, cell_mm=None, simplify=False,
                   name=None):
    """Build a DSL command and return its field as a one-part scene.

    The field is what the renderer wants, so this asks the DSL for an SDF
    representation.  Note that the DSL still meshes each primitive on the way
    through, and for a large assembly that meshing dominates the time to a
    first image.
    """
    from yapcad.dsl import compile_and_run
    source = open(path).read()
    started = time.time()
    result = compile_and_run(source, command, dict(params or {}),
                             representation="sdf", sdf_cell_mm=cell_mm,
                             sdf_simplify=simplify, source_path=str(path))
    if not result.success:
        raise RuntimeError(
            f"{path}: building {command} failed: {result.error_message}")
    geometry = result.geometry
    if isinstance(geometry, list) and geometry and geometry[0] != "solid":
        geometry = geometry[0]          # a command that emitted one solid
    node = _as_node(geometry)
    if node is None:
        raise RuntimeError(
            f"{command} did not produce an SDF-authored solid; the renderer "
            "traces fields, so the design has to be built with "
            "representation=\"sdf\"")
    return [Part(node, name=name or command)], time.time() - started


__all__ = ["Part", "as_parts", "scene_bounds", "render", "image",
           "section_view", "four_up", "parts_from_dsl", "Camera", "VIEWS",
           "ORTHOGRAPHIC", "PERSPECTIVE"]
