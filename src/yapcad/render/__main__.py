"""``python -m yapcad.render design.dsl COMMAND`` -- build a part, look at it.

Writes a PNG and prints one line saying how long each stage took and what the
image scale is, because those are the two things worth knowing before trusting
what you are about to look at.
"""

import argparse
import sys
import time

from yapcad.render.api import four_up, image, parts_from_dsl, section_view
from yapcad.render.raster import write_png
from yapcad.render.camera import ORTHOGRAPHIC, PERSPECTIVE, VIEWS


def _size(text):
    try:
        width, height = text.lower().split("x")
        return (int(width), int(height))
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"expected WIDTHxHEIGHT, got {text!r}") from None


def _value(text):
    """Guess a DSL parameter's type from how it is written."""
    if text.lower() in ("true", "false"):
        return text.lower() == "true"
    for cast in (int, float):
        try:
            return cast(text)
        except ValueError:
            pass
    return text


def _param(text):
    if "=" not in text:
        raise argparse.ArgumentTypeError(
            f"expected NAME=VALUE, got {text!r}")
    name, _, raw = text.partition("=")
    return (name.strip(), _value(raw.strip()))


def _section_spec(text):
    axis, _, where = str(text).partition("=")
    axis = axis.strip().lower()
    if axis not in ("x", "y", "z"):
        raise argparse.ArgumentTypeError(
            f"expected an axis x, y or z, optionally AXIS=MM; got {text!r}")
    if not where.strip():
        return (axis, None)
    try:
        return (axis, float(where))
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"{where!r} is not a coordinate in millimetres") from None


def build_parser():
    parser = argparse.ArgumentParser(
        prog="python -m yapcad.render",
        description="Render a DSL part as an engineering image, headlessly, "
                    "by sphere-tracing its field.")
    parser.add_argument("design", help="path to a .dsl design")
    parser.add_argument("command", help="the command to build and render")
    parser.add_argument("--view", default="iso",
                        help="one of " + ", ".join(VIEWS)
                             + ", 4up, or section (default: iso)")
    parser.add_argument("--section", type=_section_spec, metavar="AXIS[=MM]",
                        help="cut plane for a section view, or for the "
                             "fourth tile of 4up (default: z at mid-height)")
    parser.add_argument("--size", type=_size, default=(640, 480),
                        metavar="WxH", help="image size (default: 640x480)")
    parser.add_argument("--theme", default="light", choices=("light", "dark"))
    parser.add_argument("--perspective", action="store_true",
                        help="perspective projection; not to scale, so no "
                             "scale bar is drawn")
    parser.add_argument("--colour-by", default="auto",
                        choices=("auto", "part", "group", "none"),
                        dest="colour_by")
    parser.add_argument("--cell", type=float, default=None, metavar="MM",
                        help="SDF cell size the DSL meshes primitives at; "
                             "larger is faster to build and does not affect "
                             "the render, which uses the field itself")
    parser.add_argument("--simplify", action="store_true",
                        help="let the DSL simplify parts against their "
                             "fields while building (slower)")
    parser.add_argument("--param", type=_param, action="append", default=[],
                        metavar="NAME=VALUE", help="a DSL parameter; repeat")
    parser.add_argument("--max-steps", type=int, default=256,
                        dest="max_steps", help="march iteration cap")
    parser.add_argument("--no-annotate", action="store_false",
                        dest="annotate",
                        help="omit the title, scale bar, triad and legend")
    parser.add_argument("-o", "--output", default=None,
                        help="PNG path (default: COMMAND_VIEW.png)")
    parser.add_argument("-q", "--quiet", action="store_true")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    view = args.view.lower()
    known = set(VIEWS) | {"4up", "section"}
    if view not in known:
        build_parser().error(
            f"unknown view {args.view!r}; expected one of "
            f"{', '.join(sorted(known))}")

    parts, build_seconds = parts_from_dsl(
        args.design, args.command, params=dict(args.param),
        cell_mm=args.cell, simplify=args.simplify)
    axis, where = args.section or ("z", None)
    out = args.output or f"{args.command.lower()}_{view}.png"

    started = time.time()
    if view == "4up":
        img = four_up(parts, section_axis=axis, section_at=where,
                      pixels=(args.size[0] // 2, args.size[1] // 2),
                      theme=args.theme, colour_by=args.colour_by,
                      max_steps=args.max_steps)
    elif view == "section":
        img = section_view(parts, axis=axis, where=where, pixels=args.size,
                           theme=args.theme, annotate=args.annotate)
    else:
        img = image(parts, view=view, pixels=args.size, theme=args.theme,
                    projection=PERSPECTIVE if args.perspective
                    else ORTHOGRAPHIC,
                    colour_by=args.colour_by, annotate=args.annotate,
                    max_steps=args.max_steps)
    drawn = time.time() - started
    write_png(out, img)
    if not args.quiet:
        print(f"{args.command}: field built in {build_seconds:.1f}s, "
              f"{view} drawn in {drawn:.2f}s at "
              f"{img.shape[1]}x{img.shape[0]} -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
