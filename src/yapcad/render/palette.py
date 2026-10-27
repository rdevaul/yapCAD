"""Colour for engineering renders: two selected themes, one validated palette.

Categorical hues are assigned to *groups* -- a material, a subassembly,
whatever the caller keys on -- and never to individual parts.  Eight slots are
the whole palette: a ninth group folds into the neutral ``other`` tone rather
than inventing a hue.  Part identity in a crowded view comes from the boundary
lines :func:`yapcad.render.buffers.composite` draws between parts and from the
legend, so it is never carried by colour alone.

The dark theme is a *selected* set of steps for the dark surface, not an
automatic inversion of the light one.
"""

from dataclasses import dataclass, field


def from_hex(text):
    """``"#2a78d6"`` -> ``(0.165, 0.471, 0.839)``."""
    s = text.lstrip("#")
    if len(s) != 6:
        raise ValueError(f"expected a six-digit hex colour, got {text!r}")
    return tuple(int(s[i:i + 2], 16) / 255.0 for i in (0, 2, 4))


#: Eight categorical hues, in the fixed order that keeps every adjacent pair
#: separable under protanopia and deuteranopia.  The order is the safety
#: mechanism, not decoration: do not sort, cycle or extend it.
_SERIES_LIGHT = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100",
                 "#e87ba4", "#008300", "#6250d6", "#e34948")
_SERIES_DARK = ("#3987e5", "#d95926", "#199e70", "#c98500",
                "#d55181", "#008300", "#9085e9", "#e66767")

#: One hue, light to dark, for magnitude.  Used by the section view to show
#: how far a point is from material.
_RAMP = ("#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
         "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281",
         "#0d366b")


@dataclass(frozen=True)
class Theme:
    """The colours one render is drawn with."""

    name: str
    surface: tuple
    ink: tuple
    ink_soft: tuple
    line: tuple
    material: tuple
    other: tuple
    series: tuple
    ramp: tuple = field(default_factory=lambda: tuple(
        from_hex(h) for h in _RAMP))

    def group_colour(self, slot):
        """Base colour for group ``slot``; ``None`` or out of range -> neutral."""
        if slot is None or slot < 0:
            return self.material
        if slot >= len(self.series):
            return self.other
        return self.series[slot]


LIGHT = Theme(
    name="light",
    surface=from_hex("#fcfcfb"),
    ink=from_hex("#0b0b0b"),
    ink_soft=from_hex("#52514e"),
    line=from_hex("#13161a"),
    material=from_hex("#a8adb3"),
    other=from_hex("#8b8f94"),
    series=tuple(from_hex(h) for h in _SERIES_LIGHT),
)

DARK = Theme(
    name="dark",
    surface=from_hex("#1a1a19"),
    ink=from_hex("#f0efec"),
    ink_soft=from_hex("#c3c2b7"),
    line=from_hex("#07080a"),
    material=from_hex("#b9bec4"),
    other=from_hex("#9398a0"),
    series=tuple(from_hex(h) for h in _SERIES_DARK),
)

THEMES = {"light": LIGHT, "dark": DARK}


def get_theme(theme):
    """Accept a :class:`Theme` or one of ``"light"``/``"dark"``."""
    if isinstance(theme, Theme):
        return theme
    try:
        return THEMES[str(theme)]
    except KeyError:
        raise ValueError(
            f"unknown theme {theme!r}; expected one of "
            f"{', '.join(sorted(THEMES))}") from None


def assign_groups(keys):
    """Map group keys to palette slots in first-appearance order.

    An empty key means "no group" and renders in the neutral material tone.
    The ninth and later distinct keys share the ``other`` tone, because a
    generated ninth hue would not be separable from the eight above it.

    :return: ``{key: slot}``, where ``slot`` is ``-1`` for the neutral tone
        and ``>= len(series)`` for the folded remainder.
    """
    slots = {}
    for key in keys:
        if not key:
            slots.setdefault("", -1)
            continue
        if key not in slots:
            slots[key] = len([k for k in slots if k])
    return slots
