"""Plain-text reports for a range-of-motion sweep.

:meth:`~yapcad.rom.sweep.SweepResult.summary` is the short form for a log
or a terminal; :func:`interference_summary` is the long form, with the
model it swept, the settings it used, every pair's worst configuration,
and -- just as prominently -- whatever it could not measure.
"""

from __future__ import annotations

import math

from .sweep import SweepResult

__all__ = ["interference_summary"]

_RULE = "=" * 70


def _degrees(config) -> str:
    if not config:
        return "(no joints)"
    return ", ".join(
        f"{name}={math.degrees(value):.1f} deg"
        for name, value in sorted(config.items())
    )


def interference_summary(result: SweepResult) -> str:
    """A full report of one sweep, as text."""
    cfg = result.config
    name = getattr(result.model.assembly, "name", "?")

    lines = [
        _RULE,
        f"  Range-of-motion sweep: {name}",
        _RULE,
        "",
        result.model.describe(),
        "",
        f"  configurations:         {result.configurations}",
        f"  pairs measured:         {len(result.pairs)}",
        f"  sampling pitch:         {cfg.pitch_mm:.3f} mm",
        f"  interference below:     {cfg.interference_threshold:.3f} mm",
        f"  near miss below:        {cfg.clearance_threshold:.3f} mm",
        "",
    ]

    if result.skipped:
        lines.append(f"  NOT MEASURED ({len(result.skipped)} pair(s)):")
        for pair, reason in sorted(result.skipped.items()):
            lines.append(f"    {pair[0]} - {pair[1]}: {reason}")
        lines.append("")

    if not result.samples:
        lines.extend([
            "  Nothing was measured, so this sweep says nothing about the",
            "  assembly. It is not a clean result.",
            "",
            _RULE,
        ])
        return "\n".join(lines)

    if result.stopped_early:
        lines.extend([
            "  The sweep stopped early, so the remaining configurations",
            "  were never measured.",
            "",
        ])

    minima = result.pair_minima()
    lines.append("  --- Worst configuration per pair ---")
    lines.append(
        f"  {'pair':40s} {'min (mm)':>10s} {'+/-':>7s} {'how':>10s}"
    )
    lines.append("  " + "-" * 70)
    conservative = False
    for pair, sample in sorted(minima.items(),
                               key=lambda item: item[1].min_distance):
        label = f"{pair[0]} - {pair[1]}"
        clearance = sample.clearance
        how = clearance.method
        if clearance.conservative:
            how += " >="
            conservative = True
        lines.append(
            f"  {label:40s} {clearance.min_distance:10.3f} "
            f"{clearance.tolerance_mm:7.3f} {how:>10s}"
        )
    if conservative:
        lines.append(
            "  >= : the fields are not exact distances (hard CSG), so the "
            "magnitude"
        )
        lines.append(
            "       is closer to zero than the truth -- a gap is never "
            "overstated."
        )
    lines.append("")

    events = result.events
    if events:
        lines.append(f"  --- Interference ({len(events)}) ---")
        lines.append("")
        for index, sample in enumerate(
            sorted(events, key=lambda s: s.min_distance), 1
        ):
            clearance = sample.clearance
            lines.append(f"  #{index}: {_degrees(sample.config)}")
            lines.append(
                f"       {sample.pair[0]} - {sample.pair[1]}: "
                f"{clearance.min_distance:.3f} mm "
                f"(penetration {clearance.penetration:.3f} mm, "
                f"+/- {clearance.tolerance_mm:.3f})"
            )
            if clearance.method == "field":
                lines.append(
                    f"       {clearance.samples} sample(s), "
                    f"{clearance.fraction_inside * 100:.1f}% inside, "
                    f"mean {clearance.mean_distance:.3f} mm"
                )
            lines.append("")
    else:
        lines.append("  No interference found.")
        lines.append("")

    near = result.near_misses
    if near:
        lines.append(f"  --- Near misses ({len(near)}, worst 5) ---")
        for sample in sorted(near, key=lambda s: s.min_distance)[:5]:
            lines.append(
                f"    {_degrees(sample.config)} | {sample.pair[0]} - "
                f"{sample.pair[1]}: {sample.min_distance:.3f} mm"
            )
        lines.append("")

    lines.append(_RULE)
    return "\n".join(lines)
