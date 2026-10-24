"""Report generation for ROM sweep results.

Produces human-readable summary tables and, in future, HTML/JSON exports.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

from .sweep import SweepResult, InterferenceEvent


def interference_summary(result: SweepResult) -> str:
    """Generate a detailed plain-text report of all interferences.

    Args:
        result: The :class:`SweepResult` from a sweep run.

    Returns:
        Multi-line string suitable for printing or export.
    """
    cfg = result.config
    tree = result.tree

    lines = [
        f"{'=' * 60}",
        f"  ROM Sweep Report: {tree.name}",
        f"{'=' * 60}",
        "",
        f"  Configurations evaluated:  {result.n_configurations}",
        f"  Adjacent pairs checked:    {len(result.pair_trace)}",
        f"  SDF grid resolution:        {cfg.sdf_resolution}^3",
        "",
        f"  Interference events:        {len(result.events)}",
        f"  Near misses (< {cfg.clearance_threshold} mm): {len(result.near_misses)}",
        "",
    ]

    if result.events:
        lines.append("  --- Interference Details ---")
        lines.append("")

        sorted_events = sorted(
            result.events,
            key=lambda e: e.min_distance,
        )

        for i, ev in enumerate(sorted_events, 1):
            config_str = ", ".join(
                f"{k}={v:.1f}" for k, v in sorted(ev.config.items())
            )
            lines.append(
                f"  #{i}: {config_str}"
            )
            lines.append(
                f"       {ev.pair[0]} ↔ {ev.pair[1]}: "
                f"{ev.min_distance:.3f} mm "
                f"(penetration {ev.penetration_depth:.2f} mm)"
            )

            pr = ev.pairwise_results
            lines.append(
                f"       mean distance: {pr.get('mean_distance', 'N/A'):.2f} mm, "
                f"within 5 mm: {pr.get('within_5mm', 0)*100:.0f}%, "
                f"points: {pr.get('grid_points', 0)}"
            )
            lines.append("")

    if result.near_misses:
        lines.append("  --- Near Misses (worst 5) ---")
        sorted_near = sorted(
            result.near_misses,
            key=lambda e: e.min_distance,
        )
        for ev in sorted_near[:5]:
            lines.append(
                f"    {', '.join(f'{k}={v:.1f}' for k,v in sorted(ev.config.items()))}"
                f" | {ev.pair[0]} ↔ {ev.pair[1]}: "
                f"{ev.min_distance:.2f} mm"
            )

    # Pair minima summary table
    if result.pair_trace:
        lines.append("")
        lines.append("  --- Pair Clearance Minima ---")
        lines.append(
            f"  {'Pair':40s} {'Min (mm)':>10s} {'Mean (mm)':>10s} {'Configs':>8s}"
        )
        lines.append("  " + "-" * 68)
        for pair, trace in result.pair_trace.items():
            if not trace:
                continue
            min_min = min(t[1] for t in trace)
            mean_vals = [t[2] for t in trace if t[2] != float("inf")]
            mean_mean = sum(mean_vals) / len(mean_vals) if mean_vals else 0.0
            pair_label = f"{pair[0]} ↔ {pair[1]}"
            lines.append(
                f"  {pair_label:40s} {min_min:10.2f} {mean_mean:10.2f} {len(trace):8d}"
            )

    lines.extend(
        [
            "",
            "=" * 60,
        ]
    )
    return "\n".join(lines)


__all__ = ["interference_summary"]