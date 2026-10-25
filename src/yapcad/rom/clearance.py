"""Clearance between two parts, measured on their fields.

The measurement is one line of geometry: *sample the second part, ask the
first part's field how far away those samples are*. Both parts are SDF
nodes, so there is no mesh, no boolean, and no tessellation tolerance --
only the sampling pitch, which the caller sets and the result reports.

    for each grid point p inside part B (B's own frame)
        d = SDF_A( A_world^-1 @ B_world @ p )
    clearance = min d

A negative result means B is inside A by that much.

Accuracy
========
Samples are kept within half a cell of B's surface -- inside it, or just
outside by less than that -- so the sampling contributes an error of at
most half the cell diagonal, which every :class:`Clearance` carries as
``tolerance_mm``. Nothing else biases it: a sampling margin wider than the
cell, as this code once had, subtracts itself from every answer.

The fields themselves are the other half of the story. yapCAD marks a
node ``exact`` when its magnitude is the true distance; hard CSG does not
preserve that, so a part built by ``union``/``intersect``/``subtract`` --
which is most parts -- carries a Lipschitz-1 *lower bound* instead. Where
either field is inexact the result is flagged ``conservative``, and then
the number is closer to zero than the truth: the gap between separated
parts is never overstated, and penetration depth may be understated. The
sign is right either way, so interference is never missed.

Set the pitch against the clearance you care about. At the default 1 mm
pitch a 0.5 mm interference is resolved; looking for a 0.1 mm one needs a
finer pitch, and the result will say so through ``tolerance_mm``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Tuple

import numpy as np

from .model import RomError

__all__ = [
    "Clearance",
    "clearance_sdf",
    "bounds_clearance",
    "world_bounds",
]

# A pair this far apart (mm, after the bounds test) is not worth sampling.
DEFAULT_PITCH_MM = 1.0
DEFAULT_MAX_SAMPLES = 2_000_000


@dataclass(frozen=True)
class Clearance:
    """What one pair measured at one configuration.

    Attributes:
        min_distance: Smallest signed distance found, in mm. Negative
            means the parts overlap by that much.
        penetration: ``-min_distance`` when negative, else 0.0.
        mean_distance: Mean signed distance over the sampled points.
        fraction_inside: Share of sampled points inside the other part.
        samples: Points actually evaluated.
        pitch_mm: Sampling pitch used, which may be coarser than asked
            for if the part needed more samples than were allowed.
        tolerance_mm: Bound on the sampling error in ``min_distance``.
            Zero for a ``"bounds"`` result, whose error is one-sided.
        method: ``"field"`` when sampled, ``"bounds"`` when the parts'
            bounding boxes were far enough apart to answer without
            sampling -- then ``min_distance`` is a lower bound.
        conservative: True when a field involved is not an exact distance
            (any hard CSG combination, and every ``"bounds"`` result). The
            magnitude is then closer to zero than the truth -- a separated
            gap is never overstated, a penetration may be understated --
            while the sign still holds. False when both fields are exact,
            and then the true distance is within ``tolerance_mm``.
    """

    min_distance: float
    penetration: float
    mean_distance: float
    fraction_inside: float
    samples: int
    pitch_mm: float
    tolerance_mm: float
    method: str = "field"
    conservative: bool = True

    @property
    def interferes(self) -> bool:
        return self.min_distance < 0.0

    @property
    def resolved(self) -> bool:
        """False when the result is smaller than its own error bound."""
        return abs(self.min_distance) > self.tolerance_mm

    def describe(self) -> str:
        bound = "at least " if self.conservative else ""
        return (f"{bound}{self.min_distance:+.3f} mm "
                f"(+/- {self.tolerance_mm:.3f}, {self.method})")


def _matrix(transform) -> np.ndarray:
    """A 4x4 matrix from a matrix, a nested list, or a ``Transform``."""
    if transform is None:
        return np.eye(4)
    candidate = getattr(transform, "matrix", transform)
    matrix = np.asarray(candidate, dtype=float)
    if matrix.shape != (4, 4):
        raise RomError(
            f"A pose must be a 4x4 transform, got shape {matrix.shape}."
        )
    return matrix


def _finite_bounds(field, role: str) -> Tuple[Tuple[float, ...], ...]:
    bounds = getattr(field, "bounds", None)
    if bounds is None:
        raise RomError(f"The {role} field has no bounds.")
    lo, hi = bounds
    for value in (*lo, *hi):
        if not math.isfinite(value):
            raise RomError(
                f"The {role} field is unbounded ({bounds}), so it cannot be "
                "sampled. Intersect it with a finite region first."
            )
    if any(high < low for low, high in zip(lo, hi)):
        raise RomError(f"The {role} field is empty ({bounds}).")
    return (tuple(float(v) for v in lo), tuple(float(v) for v in hi))


def world_bounds(field, transform) -> Tuple[Tuple[float, ...], ...]:
    """Axis-aligned world bounds of ``field`` placed by ``transform``."""
    lo, hi = _finite_bounds(field, "placed")
    corners = np.array([[x, y, z, 1.0]
                        for x in (lo[0], hi[0])
                        for y in (lo[1], hi[1])
                        for z in (lo[2], hi[2])])
    placed = (_matrix(transform) @ corners.T).T[:, :3]
    return (tuple(placed.min(axis=0)), tuple(placed.max(axis=0)))


def bounds_clearance(field_a, field_b, tf_a, tf_b) -> float:
    """A lower bound on the distance between two placed fields.

    Cheap, and exact enough to dismiss a pair: if the bounding boxes are
    8 mm apart the parts are at least 8 mm apart. Returns 0.0 when the
    boxes overlap, which says nothing about the parts themselves.
    """
    lo_a, hi_a = world_bounds(field_a, tf_a)
    lo_b, hi_b = world_bounds(field_b, tf_b)
    gaps = [max(0.0, lo_b[i] - hi_a[i], lo_a[i] - hi_b[i])
            for i in range(3)]
    return float(math.sqrt(sum(g * g for g in gaps)))


def _grid(lo, hi, pitch: float, max_samples: int):
    """Sample counts per axis, and the pitch actually achievable."""
    extents = [max(hi[i] - lo[i], 0.0) for i in range(3)]
    used = float(pitch)
    for _ in range(32):
        counts = [max(2, int(math.ceil(e / used)) + 1) for e in extents]
        total = counts[0] * counts[1] * counts[2]
        if total <= max_samples:
            return counts, used
        # Coarsen by the cube root of the overshoot; the result records
        # the pitch it ended up with, so this cannot pass unnoticed.
        used *= (total / max_samples) ** (1.0 / 3.0)
    raise RomError(
        f"Sampling this part needs more than {max_samples} points even at "
        f"a {used:.3f} mm pitch; raise max_samples or coarsen the pitch."
    )


def clearance_sdf(
    field_a,
    field_b,
    tf_a=None,
    tf_b=None,
    *,
    pitch_mm: float = DEFAULT_PITCH_MM,
    max_samples: int = DEFAULT_MAX_SAMPLES,
) -> Clearance:
    """Measure the clearance of ``field_b`` against ``field_a``.

    Args:
        field_a: SDF node whose field is evaluated.
        field_b: SDF node that is sampled.
        tf_a: World pose of A (4x4 matrix or ``Transform``; identity if
            None).
        tf_b: World pose of B.
        pitch_mm: Target spacing between samples, in mm. The tolerance of
            the answer follows from it.
        max_samples: Ceiling on points sampled for one measurement.

    Returns:
        A :class:`Clearance`. ``min_distance`` is negative when B reaches
        inside A.
    """
    from yapcad.sdf import evaluate, sample_grid

    if pitch_mm <= 0.0:
        raise RomError(f"pitch_mm must be positive, got {pitch_mm}.")

    lo, hi = _finite_bounds(field_b, "sampled")
    counts, pitch = _grid(lo, hi, pitch_mm, max_samples)

    # A part thinner than one cell on some axis still needs two planes of
    # samples either side of its surface.
    region_lo, region_hi = list(lo), list(hi)
    for axis in range(3):
        if region_hi[axis] - region_lo[axis] < pitch:
            middle = 0.5 * (region_lo[axis] + region_hi[axis])
            region_lo[axis] = middle - 0.5 * pitch
            region_hi[axis] = middle + 0.5 * pitch

    values, axes = sample_grid(
        field_b, bounds=(tuple(region_lo), tuple(region_hi)),
        resolution=counts,
    )

    steps = [(region_hi[i] - region_lo[i]) / max(counts[i] - 1, 1)
             for i in range(3)]
    tolerance = 0.5 * math.sqrt(sum(s * s for s in steps))

    # Keep the samples that sit at or inside B's surface, within half a
    # cell of it. Points further out belong to no part and would drag the
    # minimum down with them.
    keep = (values <= tolerance).ravel()
    if not keep.any():
        raise RomError(
            "The sampled part has no interior at this pitch "
            f"({pitch:.3f} mm); use a finer pitch."
        )

    mesh = np.meshgrid(*axes, indexing="ij")
    points = np.stack([m.ravel() for m in mesh], axis=1)[keep]

    homogeneous = np.ones((len(points), 4))
    homogeneous[:, :3] = points
    world = (_matrix(tf_b) @ homogeneous.T).T
    world[:, 3] = 1.0
    local_a = (np.linalg.inv(_matrix(tf_a)) @ world.T).T[:, :3]

    distances = evaluate(field_a, local_a)
    minimum = float(distances.min())
    exact = bool(getattr(field_a, "exact", False)
                 and getattr(field_b, "exact", False))
    return Clearance(
        min_distance=minimum,
        penetration=-minimum if minimum < 0.0 else 0.0,
        mean_distance=float(distances.mean()),
        fraction_inside=float((distances < 0.0).mean()),
        samples=int(distances.size),
        pitch_mm=float(pitch),
        tolerance_mm=float(tolerance),
        method="field",
        conservative=not exact,
    )
