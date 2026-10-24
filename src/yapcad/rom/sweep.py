"""Range-of-Motion sweep: iterate DOF space and detect interferences.

Orchestrates a kinematic sweep across one or more joint DOFs, computing
pairwise clearance for each configuration and recording any interference
events.  The sweep can be:

- **Single-DOF**: Sweep one joint through its range while holding others at
  zero pose.
- **Multi-DOF product space**: Sweep combinations (use with small step counts
  for higher dimensions to avoid combinatorial explosion).
- **Specific configurations**: Pass an explicit list of angle dicts.

The clearance for each configuration is computed via the SDF field (no
meshing), so even a 4-DOF sweep at coarse resolution completes in seconds.
"""

from __future__ import annotations

import itertools
import logging
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from .tree import RomTree
from .clearance import clearance_sdf, clearance_mesh

logger = logging.getLogger(__name__)


@dataclass
class InterferenceEvent:
    """Record of one clearance violation at a specific configuration.

    Attributes:
        config: Joint angle dict for this configuration.
        pair: ``(parent_name, child_name)`` adjacent pair.
        min_distance: Minimum signed distance (mm).  Negative = interference.
        penetration_depth: Magnitude of penetration (mm).
        pairwise_results: Full dict from the clearance function.
    """

    config: Dict[str, float]
    pair: Tuple[str, str]
    min_distance: float
    penetration_depth: float
    pairwise_results: Dict[str, Any] = field(default_factory=dict)


@dataclass
class SweepConfig:
    """Configuration for a kinematic ROM sweep.

    Attributes:
        step_deg: Default step size in degrees for revolute joints.
        step_mm: Default step size in mm for prismatic joints.
        clearance_threshold: Distance below which is flagged as ``near`` (mm).
        interference_threshold: Distance below which is flagged as
            interference (mm).  Usually zero or slightly negative to
            account for numerical tolerance.
        sdf_resolution: Grid resolution per axis for clearance evaluation.
        sdf_padding: Padding around the child bounding box (mm).
        max_pairs: Maximum number of adjacent pairs to check per config.
            Set to 0 for "no limit".
        progress: Log progress every N configurations (0 = no logging).
        early_stop: Stop after this many interference events (0 = no limit).
    """

    step_deg: float = 5.0
    step_mm: float = 2.0
    clearance_threshold: float = 5.0
    interference_threshold: float = 0.0
    sdf_resolution: int = 32
    sdf_padding: float = 2.0
    max_pairs: int = 0
    progress: int = 50
    early_stop: int = 0


@dataclass
class SweepResult:
    """Aggregated results from a kinematic sweep.

    Attributes:
        tree: The :class:`RomTree` used for the sweep.
        config: The :class:`SweepConfig` used.
        events: List of :class:`InterferenceEvent` (clearance violations).
        near_misses: Configurations where clearance was below threshold.
        pair_trace: Full clearance trace per pair::
            ``{(parent, child): [(config_name, min_dist, mean_dist), ...]}``
        n_configurations: Total configurations evaluated.
        clearance_summary: Per-dof-name or per-configuration summary.
    """

    tree: RomTree
    config: SweepConfig = field(default_factory=SweepConfig)
    events: List[InterferenceEvent] = field(default_factory=list)
    near_misses: List[InterferenceEvent] = field(default_factory=list)
    pair_trace: Dict[Tuple[str, str], list] = field(default_factory=dict)
    n_configurations: int = 0
    clearance_summary: Dict[str, Any] = field(default_factory=dict)

    def summary(self) -> str:
        """Human-readable summary of the sweep results."""
        if self.n_configurations == 0:
            return "No configurations were evaluated."

        lines = [
            f"Sweep: {self.tree.name}",
            f"  Configurations evaluated: {self.n_configurations}",
            f"  Interference events:      {len(self.events)}",
            f"  Near misses (< {self.config.clearance_threshold} mm): {len(self.near_misses)}",
        ]

        if self.events:
            lines.append("\n  Worst interferences:")
            for ev in sorted(
                self.events, key=lambda e: e.min_distance
            )[:5]:
                config_str = _format_config(ev.config)
                lines.append(
                    f"    {config_str} | "
                    f"{ev.pair[0]} ↔ {ev.pair[1]}: "
                    f"{ev.min_distance:.2f} mm "
                    f"(penetration {ev.penetration_depth:.2f} mm)"
                )

        if self.pair_trace:
            lines.append("\n  Pair clearance minima:")
            for pair, trace in self.pair_trace.items():
                valid = [t[1] for t in trace if t[1] != float('inf')]
                if not valid:
                    continue
                minima = min(valid)
                maxima = max(valid)
                lines.append(
                    f"    {pair[0]} ↔ {pair[1]}: "
                    f"[{minima:.2f}, {maxima:.2f}] mm"
                )

        return "\n".join(lines)

    def worst_interference(self) -> Optional[InterferenceEvent]:
        """Return the event with the most negative clearance."""
        if not self.events:
            return None
        return min(self.events, key=lambda e: e.min_distance)

    def clearance_at(self, part_a: str, part_b: str) -> list:
        """Get the clearance trace for a specific pair.

        Returns:
            List of ``(config_key, min_distance, mean_distance)`` tuples.
        """
        return self.pair_trace.get((part_a, part_b), [])


# ---------------------------------------------------------------------------
# Sweep orchestration
# ---------------------------------------------------------------------------

def dof_schedule_from_tree(
    tree: RomTree,
    step_deg: float = 5.0,
) -> Dict[str, Tuple[float, float, float]]:
    """Build a default DOF sweep schedule from the tree's joint limits.

    For each link with a revolute or prismatic joint, produces
    ``(min, max + step, step)``, suitable for ``numpy.arange`` or
    the sweep runner.

    Args:
        tree: The kinematic tree.
        step_deg: Angular step size in degrees.

    Returns:
        ``{link_name: (start, stop, step)}`` dict.
    """
    schedule = {}
    for name, link in tree.links.items():
        lo, hi = link.limits
        jk = link.joint_kind
        if jk == "revolute":
            step = step_deg
        elif jk == "prismatic":
            step = step_deg  # Could use a separate step_mm
        else:
            continue  # Fixed — no sweep
        schedule[name] = (float(lo), float(hi + step), step)
    return schedule


def _config_key(angles: Dict[str, float]) -> str:
    """Compact string key for a configuration dict."""
    parts = []
    for name, val in sorted(angles.items()):
        if abs(round(val, 1)) == round(val, 1):
            parts.append(f"{name}={val:.0f}")
        else:
            parts.append(f"{name}={val:.1f}")
    return ", ".join(parts)


def _format_config(config: Dict[str, float]) -> str:
    short = {}
    for k, v in config.items():
        short[k] = f"{v:.1f}°" if abs(v) < 1000 else f"{v:.0f}mm"
    return ", ".join(f"{k}={v}" for k, v in short.items())


def sweep_kinematics(
    tree: RomTree,
    *sweep_schedules: Dict[str, Tuple[float, float, float]],
    config: Optional[SweepConfig] = None,
    extra_pairs: Optional[List[Tuple[str, str]]] = None,
    fixed_angles: Optional[Dict[str, float]] = None,
) -> SweepResult:
    """Run a ROM sweep across one or more DOFs.

    Args:
        tree: The kinematic tree with link metadata.
        *sweep_schedules: One or more ``{link_name: (start, stop, step)}``
            dicts.  Multiple schedules produce a product space (use with
            caution beyond 3 dimensions).
        config: Sweep config (step sizes, thresholds, resolution).
            Uses defaults if ``None``.
        extra_pairs: Additional part pairs to check beyond the tree's
            adjacent pairs.  ``[(name_a, name_b), ...]``.
        fixed_angles: Joint angles to hold at fixed values across all
            configurations.  These are applied *after* the sweep
            DOF values.

    Returns:
        :class:`SweepResult` with all events and pair traces.

    Example::

        result = sweep_kinematics(
            tree,
            {"HIP_ROLL": (-45, 46, 5)},
            {"KNEE_PITCH": (0, 136, 5)},
            fixed_angles={"HIP_YAW": 0.0},
        )
        print(result.summary())
    """
    cfg = config or SweepConfig()
    result = SweepResult(tree=tree, config=cfg)
    fixed = fixed_angles or {}

    # Build configuration list
    # Convert each schedule to a list of values
    dof_values: Dict[str, List[float]] = {}
    for schedule in sweep_schedules:
        for name, (start, stop, step) in schedule.items():
            # Determine if degrees or mm based on joint kind
            dof_values[name] = np.arange(start, stop, step).tolist()

    if not dof_values:
        # If no schedules provided, use the tree's default schedule
        dof_values = dof_schedule_from_tree(tree, cfg.step_deg)

    # Generate all configurations (product space if multiple DOFs)
    dof_names = list(dof_values.keys())
    value_lists = [dof_values[n] for n in dof_names]

    all_configs = []
    for combo in itertools.product(*value_lists):
        angles = dict(zip(dof_names, combo))
        angles.update(fixed)
        all_configs.append(angles)

    # Identify adjacent pairs to check
    pairs_to_check: List[Tuple[str, str]] = list(tree.adjacent_pairs())
    if extra_pairs:
        # Only add pairs where both parts exist in the tree
        link_names = set(tree.links.keys())
        for a, b in extra_pairs:
            if a in link_names and b in link_names:
                pairs_to_check.append((a, b))

    # Limit pair count if configured
    if cfg.max_pairs > 0:
        pairs_to_check = pairs_to_check[:cfg.max_pairs]

    if not pairs_to_check:
        logger.warning("No adjacent pairs to check; sweep is a no-op.")

    # Initialize pair traces
    for pair in pairs_to_check:
        result.pair_trace[pair] = []

    # Sweep loop
    n_configs = len(all_configs)
    result.n_configurations = n_configs

    for idx, angles in enumerate(all_configs):
        tree.set_joints(angles)

        config_key = _config_key(angles)

        if cfg.progress and idx > 0 and idx % cfg.progress == 0:
            logger.info(
                "Sweep %s: %d / %d (%.0f%%) ...",
                tree.name,
                idx,
                n_configs,
                100.0 * idx / n_configs,
            )

        for pair in pairs_to_check:
            p, c = pair
            node_p = tree.get_sdf_node(p)
            node_c = tree.get_sdf_node(c)
            tf_p = tree.get_world_transform(p)
            tf_c = tree.get_world_transform(c)

            has_sdf = node_p is not None and node_c is not None

            if has_sdf:
                clr = clearance_sdf(
                    node_p, node_c,
                    tf_p, tf_c,
                    resolution=cfg.sdf_resolution,
                    padding=cfg.sdf_padding,
                )
            else:
                # Fallback: we don't have an SDF pair — skip or warn
                logger.debug("No SDF for pair %s ↔ %s; skipping", p, c)
                continue

            min_dist = clr.get("min_distance", float("inf"))
            mean_dist = clr.get("mean_distance", float("inf"))

            # Record trace
            result.pair_trace[pair].append((config_key, min_dist, mean_dist))

            # Check for interference
            if min_dist < cfg.interference_threshold:
                event = InterferenceEvent(
                    config=angles.copy(),
                    pair=pair,
                    min_distance=min_dist,
                    penetration_depth=clr.get("penetration_depth", 0.0),
                    pairwise_results=clr,
                )
                result.events.append(event)

            # Check near miss
            if (
                min_dist < cfg.clearance_threshold
                and min_dist >= cfg.interference_threshold
            ):
                near = InterferenceEvent(
                    config=angles.copy(),
                    pair=pair,
                    min_distance=min_dist,
                    penetration_depth=0.0,
                    pairwise_results=clr,
                )
                result.near_misses.append(near)

            # Early stop
            if cfg.early_stop and len(result.events) >= cfg.early_stop:
                logger.info(
                    "Early stop at %d interference events (config %d / %d)",
                    len(result.events),
                    idx,
                    n_configs,
                )
                break

        if cfg.early_stop and len(result.events) >= cfg.early_stop:
            break

    # Build summary dict
    tree.clear_joints()
    _build_summary(result)

    return result


# ---------------------------------------------------------------------------
# Summary helpers
# ---------------------------------------------------------------------------

def _build_summary(result: SweepResult) -> None:
    """Populate the clearance_summary dict in place."""
    summary = {}

    # Per-pair minimums
    pair_minima = {}
    for pair, trace in result.pair_trace.items():
        if trace:
            min_min = min(t[1] for t in trace)
            mean_vals = [t[2] for t in trace if t[2] != float("inf")]
            mean_mean = sum(mean_vals) / len(mean_vals) if mean_vals else 0.0
            pair_minima[pair] = {
                "min_distance_mm": min_min,
                "mean_distance_mm": mean_mean,
                "n_configs": len(trace),
            }

    summary["pair_minima"] = pair_minima
    summary["n_interferences"] = len(result.events)
    summary["n_near_misses"] = len(result.near_misses)

    if result.events:
        worst = result.worst_interference()
        summary["worst_interference"] = {
            "pair": list(worst.pair),
            "min_distance_mm": worst.min_distance,
            "penetration_depth_mm": worst.penetration_depth,
            "config": worst.config,
        }
    else:
        summary["worst_interference"] = None

    result.clearance_summary = summary


__all__ = [
    "sweep_kinematics",
    "SweepResult",
    "SweepConfig",
    "InterferenceEvent",
    "dof_schedule_from_tree",
]