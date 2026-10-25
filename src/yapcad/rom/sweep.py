"""Sweep an assembly's joints and measure clearance at every pose.

One sweep is: a set of joint values per configuration, the assembly's own
solver for the poses, and :func:`~yapcad.rom.clearance.clearance_sdf` for
each pair of parts that is worth measuring. Nothing is meshed.

    model = RomModel.from_result(result)
    schedule = schedule_from_model(model, step_deg=5.0)
    sweep = sweep_rom(model, schedule)
    print(sweep.summary())

What the sweep will not do is look clean when it measured nothing. A pair
it cannot measure is an error unless ``skip_unmeasurable`` is set, in
which case it is counted, named with its reason, and reported in every
summary. Joint values outside the design's own limits are refused rather
than clamped, and an unknown joint name is refused rather than ignored.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field as dataclass_field
from typing import (
    Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union,
)

from .clearance import (
    Clearance, DEFAULT_MAX_SAMPLES, DEFAULT_PITCH_MM, bounds_clearance,
    clearance_sdf,
)
from .model import RomError, RomModel

logger = logging.getLogger(__name__)

__all__ = [
    "SweepConfig",
    "PairSample",
    "SweepResult",
    "sweep_rom",
    "schedule_from_model",
    "radian_range",
]

Pair = Tuple[str, str]
Schedule = Mapping[str, Sequence[float]]


@dataclass(frozen=True)
class SweepConfig:
    """How a sweep measures, and what it calls interference.

    Attributes:
        pitch_mm: Sampling pitch handed to the clearance measurement. The
            tolerance of every result follows from it.
        clearance_threshold: Distances below this (mm) are reported as
            near misses. Pairs whose bounding boxes are further apart
            than this are resolved from the bounds alone.
        interference_threshold: Distances below this (mm) are reported as
            interference. 0.0 means "touching or overlapping".
        max_samples: Ceiling on sample points for one measurement.
        skip_unmeasurable: Carry on when a pair has no field on one side,
            recording the reason, instead of raising.
        progress: Log progress every N configurations (0 = never).
        early_stop: Stop after this many interference events (0 = never).
    """

    pitch_mm: float = DEFAULT_PITCH_MM
    clearance_threshold: float = 5.0
    interference_threshold: float = 0.0
    max_samples: int = DEFAULT_MAX_SAMPLES
    skip_unmeasurable: bool = False
    progress: int = 0
    early_stop: int = 0


@dataclass(frozen=True)
class PairSample:
    """One pair of parts, measured at one configuration."""

    pair: Pair
    config: Mapping[str, float]
    clearance: Clearance

    @property
    def min_distance(self) -> float:
        return self.clearance.min_distance

    def describe(self) -> str:
        return (f"{_config_text(self.config)} | {self.pair[0]} - "
                f"{self.pair[1]}: {self.clearance.describe()}")


@dataclass
class SweepResult:
    """Everything one sweep measured.

    ``samples`` holds one :class:`PairSample` per pair per configuration;
    ``events`` and ``near_misses`` are views over it, so the thresholds
    can be reconsidered without sweeping again.
    """

    model: RomModel
    config: SweepConfig
    schedule: Dict[str, List[float]] = dataclass_field(default_factory=dict)
    pairs: Tuple[Pair, ...] = ()
    samples: List[PairSample] = dataclass_field(default_factory=list)
    skipped: Dict[Pair, str] = dataclass_field(default_factory=dict)
    configurations: int = 0
    stopped_early: bool = False

    # -- views -----------------------------------------------------------

    @property
    def events(self) -> List[PairSample]:
        limit = self.config.interference_threshold
        return [s for s in self.samples if s.min_distance < limit]

    @property
    def near_misses(self) -> List[PairSample]:
        low = self.config.interference_threshold
        high = self.config.clearance_threshold
        return [s for s in self.samples if low <= s.min_distance < high]

    @property
    def measured(self) -> int:
        return sum(1 for s in self.samples
                   if s.clearance.method == "field")

    @property
    def from_bounds(self) -> int:
        return sum(1 for s in self.samples
                   if s.clearance.method == "bounds")

    def trace(self, part_a: str, part_b: str) -> List[PairSample]:
        """Every measurement of one pair, in sweep order."""
        want = {(part_a, part_b), (part_b, part_a)}
        return [s for s in self.samples if s.pair in want]

    def worst(self) -> Optional[PairSample]:
        """The smallest clearance measured, or None if nothing was."""
        if not self.samples:
            return None
        return min(self.samples, key=lambda s: s.min_distance)

    def pair_minima(self) -> Dict[Pair, PairSample]:
        """The worst configuration for each pair."""
        minima: Dict[Pair, PairSample] = {}
        for sample in self.samples:
            current = minima.get(sample.pair)
            if current is None or sample.min_distance < current.min_distance:
                minima[sample.pair] = sample
        return minima

    # -- reporting -------------------------------------------------------

    def summary(self) -> str:
        """A short report that cannot be mistaken for a clean result."""
        lines = [
            f"Sweep of '{getattr(self.model.assembly, 'name', '?')}': "
            f"{self.configurations} configuration(s), "
            f"{len(self.pairs)} pair(s)",
        ]
        if self.schedule:
            for name, values in self.schedule.items():
                lines.append(f"  {name}: {_range_text(values)}")

        if not self.samples:
            lines.append("  NOTHING WAS MEASURED -- this is not a clean "
                         "result.")
        else:
            lines.append(
                f"  measurements: {self.measured} sampled, "
                f"{self.from_bounds} resolved from bounding boxes"
            )
            worst = self.worst()
            lines.append(
                f"  closest approach: {worst.describe()}"
            )
            lines.append(
                f"  interference (< {self.config.interference_threshold} "
                f"mm): {len(self.events)}"
            )
            lines.append(
                f"  near misses (< {self.config.clearance_threshold} mm): "
                f"{len(self.near_misses)}"
            )

        if self.skipped:
            lines.append(f"  pairs NOT measured: {len(self.skipped)}")
            for pair, reason in sorted(self.skipped.items()):
                lines.append(f"    {pair[0]} - {pair[1]}: {reason}")

        if self.stopped_early:
            lines.append("  stopped early; the sweep is incomplete.")

        events = self.events
        if events:
            lines.append("  worst interferences:")
            for sample in sorted(events, key=lambda s: s.min_distance)[:5]:
                lines.append(f"    {sample.describe()}")

        unresolved = [s for s in self.samples
                      if s.clearance.method == "field"
                      and not s.clearance.resolved]
        if unresolved:
            lines.append(
                f"  {len(unresolved)} measurement(s) smaller than the "
                f"sampling tolerance; refine pitch_mm to resolve them."
            )
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Schedules
# ---------------------------------------------------------------------------

def radian_range(
    start_deg: float, stop_deg: float, step_deg: float
) -> List[float]:
    """Degrees in, radians out, with ``stop_deg`` included.

    Joint values are radians everywhere in this package, matching
    ``set_mate_limits``; this is the convenience that keeps the degrees
    a designer thinks in out of the sweep itself.
    """
    if step_deg <= 0.0:
        raise RomError(f"step_deg must be positive, got {step_deg}.")
    if stop_deg < start_deg:
        raise RomError(
            f"stop_deg ({stop_deg}) is below start_deg ({start_deg})."
        )
    steps = int(math.floor((stop_deg - start_deg) / step_deg + 1e-9))
    values = [start_deg + i * step_deg for i in range(steps + 1)]
    if values[-1] < stop_deg - 1e-9:
        values.append(stop_deg)
    return [math.radians(v) for v in values]


def schedule_from_model(
    model: RomModel,
    step_deg: float = 5.0,
    joints: Optional[Iterable[str]] = None,
) -> Dict[str, List[float]]:
    """Sweep every movable joint across its declared limits.

    Args:
        model: The model to read joints and limits from.
        step_deg: Step between configurations, in degrees.
        joints: Only these joints (default: every movable joint).

    Raises:
        RomError: If a joint has no limits. A default range would be a
            guess about the machine, and the sweep would report whatever
            that guess allowed.
    """
    if joints is None:
        chosen = [j for j in model.movable_joints()]
    else:
        chosen = [model.joint(name) for name in joints]
        fixed = [j.name for j in chosen if not j.movable]
        if fixed:
            raise RomError(
                "These joints cannot be driven: " + ", ".join(fixed)
            )
    if not chosen:
        raise RomError(
            "This assembly has no movable joints, so there is nothing to "
            "sweep. Revolute mates are added with add_named_mate(..., "
            "'revolute', ...)."
        )
    unbounded = [j.name for j in chosen if not j.bounded]
    if unbounded:
        raise RomError(
            "No limits declared for: " + ", ".join(unbounded) +
            ". Declare them in the DSL with set_mate_limits(...), or pass "
            "an explicit schedule."
        )
    return {
        joint.name: radian_range(
            math.degrees(joint.lower), math.degrees(joint.upper), step_deg
        )
        for joint in chosen
    }


# ---------------------------------------------------------------------------
# The sweep
# ---------------------------------------------------------------------------

def _config_text(config: Mapping[str, float]) -> str:
    if not config:
        return "(no joints)"
    return ", ".join(
        f"{name}={math.degrees(value):.1f} deg"
        for name, value in sorted(config.items())
    )


def _range_text(values: Sequence[float]) -> str:
    if not values:
        return "no values"
    low = math.degrees(min(values))
    high = math.degrees(max(values))
    return f"{low:.1f} deg to {high:.1f} deg in {len(values)} step(s)"


def _volume(field) -> float:
    lo, hi = field.bounds
    return abs((hi[0] - lo[0]) * (hi[1] - lo[1]) * (hi[2] - lo[2]))


def _normalise_pairs(
    model: RomModel,
    pairs: Union[None, str, Iterable[Pair]],
    extra_pairs: Optional[Iterable[Pair]],
    driven: Iterable[str],
) -> List[Pair]:
    if pairs is None or pairs == "moving":
        chosen = list(model.moving_pairs(driven))
    elif isinstance(pairs, str):
        if pairs == "mates":
            chosen = list(model.mate_pairs())
        elif pairs == "all":
            chosen = list(model.all_pairs())
        else:
            raise RomError(
                f"pairs must be 'moving', 'mates', 'all' or a list of "
                f"pairs, got {pairs!r}."
            )
    else:
        chosen = [tuple(pair) for pair in pairs]
    chosen.extend(tuple(pair) for pair in (extra_pairs or ()))

    seen: List[Pair] = []
    for pair in chosen:
        if len(pair) != 2:
            raise RomError(f"A pair needs two part names, got {pair!r}.")
        for name in pair:
            if name not in model.links:
                raise RomError(
                    f"No part named '{name}' in this assembly. "
                    f"Parts: {', '.join(sorted(model.links))}."
                )
        if pair[0] == pair[1]:
            raise RomError(f"A part cannot be paired with itself: {pair}.")
        if pair not in seen and pair[::-1] not in seen:
            seen.append(pair)
    return seen


def _configurations(
    schedule: Schedule, fixed: Mapping[str, float]
) -> List[Dict[str, float]]:
    import itertools

    names = list(schedule)
    overlap = sorted(set(names) & set(fixed))
    if overlap:
        raise RomError(
            "These joints are both swept and fixed: " + ", ".join(overlap) +
            ". One of the two would silently win."
        )
    value_lists = [list(schedule[name]) for name in names]
    for name, values in zip(names, value_lists):
        if not values:
            raise RomError(f"The schedule for '{name}' is empty.")
    configs = []
    for combination in itertools.product(*value_lists):
        config = dict(zip(names, (float(v) for v in combination)))
        config.update({k: float(v) for k, v in fixed.items()})
        configs.append(config)
    return configs or [dict(fixed)]


def sweep_rom(
    model: RomModel,
    schedule: Optional[Schedule] = None,
    *,
    config: Optional[SweepConfig] = None,
    pairs: Union[None, str, Iterable[Pair]] = None,
    extra_pairs: Optional[Iterable[Pair]] = None,
    fixed: Optional[Mapping[str, float]] = None,
) -> SweepResult:
    """Measure clearance across a joint schedule.

    Args:
        model: The :class:`~yapcad.rom.model.RomModel` to sweep.
        schedule: ``{joint_name: [value, ...]}`` in radians. Build one
            with :func:`schedule_from_model` or :func:`radian_range`.
            Defaults to every movable joint across its limits at 5
            degrees.
        config: Measurement settings and thresholds.
        pairs: ``"moving"`` (default) for every pair the swept joints
            move relative to one another -- including pairs that share no
            mate, which is usually where interference is. ``"mates"`` for
            pairs that share a joint, ``"all"`` for every pair, or an
            explicit list.
        extra_pairs: Pairs to add to the chosen set -- a part against its
            grandparent, say, or a bumper it should never touch.
        fixed: Joint values held constant across the sweep.

    Returns:
        A :class:`SweepResult`.

    Raises:
        RomError: For an unknown joint or part, a value outside the
            design's limits, or a pair that cannot be measured (unless
            ``config.skip_unmeasurable``).
    """
    cfg = config or SweepConfig()
    held = {str(k): float(v) for k, v in (fixed or {}).items()}
    if schedule is None:
        schedule = schedule_from_model(model)
    plan = {str(name): [float(v) for v in values]
            for name, values in schedule.items()}
    if not plan:
        raise RomError(
            "The schedule names no joints, so nothing would move. Omit "
            "the schedule to sweep every movable joint across its limits, "
            "or name the joints to drive."
        )

    for name in list(plan) + list(held):
        joint = model.joint(name)
        if not joint.movable:
            raise RomError(
                f"Joint '{name}' is {joint.kind}, not revolute, so it "
                "cannot be swept."
            )
    configs = _configurations(plan, held)
    for values in configs:
        model.check_joint_values(values)

    chosen_pairs = _normalise_pairs(
        model, pairs, extra_pairs, driven=list(plan),
    )
    result = SweepResult(
        model=model, config=cfg, schedule=plan,
        configurations=len(configs),
    )

    measurable: List[Pair] = []
    for pair in chosen_pairs:
        missing = [name for name in pair if not model.links[name].has_field]
        if not missing:
            measurable.append(pair)
            continue
        reason = "; ".join(
            f"{name}: {model.links[name].reason}" for name in missing
        )
        if not cfg.skip_unmeasurable:
            raise RomError(
                f"Cannot measure {pair[0]} against {pair[1]} -- {reason}. "
                "Build the design with representation='sdf', drop the pair, "
                "or set skip_unmeasurable=True to sweep the rest."
            )
        result.skipped[pair] = reason
    result.pairs = tuple(measurable)

    if not measurable:
        logger.warning(
            "No measurable pairs: %s", result.skipped or "nothing selected"
        )
        return result

    # Posing re-solves the assembly, so leave it where it was found.
    incoming = model.joint_values()
    try:
        _run(model, cfg, configs, measurable, result)
    finally:
        try:
            model.pose(incoming)
        except RomError:  # pragma: no cover - the pose came from the model
            logger.warning(
                "Could not restore the assembly's incoming pose %s",
                incoming,
            )
    return result


def _run(
    model: RomModel,
    cfg: SweepConfig,
    configs: List[Dict[str, float]],
    measurable: List[Pair],
    result: SweepResult,
) -> None:
    """Measure every pair at every configuration, in place."""
    # Sample whichever part of each pair is smaller and evaluate the
    # other's field: the measurement is the same either way to within its
    # tolerance, but a 20 mm part beside a 400 mm one needs four orders of
    # magnitude fewer points, so the pitch asked for is the pitch used.
    # Bounds are fixed, so the choice is made once rather than per pose.
    ordered = []
    for pair in measurable:
        evaluated, sampled = pair
        if _volume(model.field(evaluated)) < _volume(model.field(sampled)):
            evaluated, sampled = sampled, evaluated
        ordered.append((pair, evaluated, sampled))

    for index, values in enumerate(configs):
        transforms = model.pose(values)
        if cfg.progress and index and index % cfg.progress == 0:
            logger.info("Sweep %d / %d", index, len(configs))

        for pair, evaluated, sampled in ordered:
            field_a = model.field(evaluated)
            field_b = model.field(sampled)
            tf_a = transforms[evaluated]
            tf_b = transforms[sampled]

            gap = bounds_clearance(field_a, field_b, tf_a, tf_b)
            if gap > cfg.clearance_threshold:
                # The bounding boxes alone prove the parts are further
                # apart than anything this sweep cares about.
                clearance = Clearance(
                    min_distance=gap, penetration=0.0, mean_distance=gap,
                    fraction_inside=0.0, samples=0, pitch_mm=cfg.pitch_mm,
                    tolerance_mm=0.0, method="bounds",
                )
            else:
                clearance = clearance_sdf(
                    field_a, field_b, tf_a, tf_b,
                    pitch_mm=cfg.pitch_mm, max_samples=cfg.max_samples,
                )
            result.samples.append(PairSample(
                pair=pair, config=dict(values), clearance=clearance,
            ))

        if cfg.early_stop and len(result.events) >= cfg.early_stop:
            result.stopped_early = True
            logger.info(
                "Stopping after %d interference event(s) at configuration "
                "%d / %d", len(result.events), index + 1, len(configs),
            )
            return
