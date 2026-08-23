"""Independent planar terrain-contact oracle for a rocker-bogie rover.

The oracle deliberately contains no assembly or mate-solver code.  It solves
the six wheel-height equations directly from the suspension link geometry and
the ideal differential closure ``q_left + q_right = 0``.  This makes it useful
as an independent reference for validating positioned assembly output.

Angles are stored and returned in radians.  Rocker angles in
``RockerBogiePose`` are relative to the chassis; bogie angles are relative to
their parent rocker.  Wheel coordinates use the rover X/Z convention.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Any, Dict, Mapping, Optional, Tuple

import numpy as np


WHEEL_NAMES = ("lf", "lm", "lr", "rf", "rm", "rr")


@dataclass(frozen=True)
class RockerBogieGeometry:
    """Planar link geometry and mechanical travel limits, in mm/radians."""

    center_to_front: Tuple[float, float] = (190.0, -75.0)
    center_to_bogie: Tuple[float, float] = (-95.0, -35.0)
    bogie_to_middle: Tuple[float, float] = (95.0, -40.0)
    bogie_to_rear: Tuple[float, float] = (-95.0, -40.0)
    wheel_radius: float = 65.0
    half_track: float = 155.0
    rocker_min: float = math.radians(-18.0)
    rocker_max: float = math.radians(18.0)
    bogie_min: float = math.radians(-35.0)
    bogie_max: float = math.radians(38.0)


@dataclass(frozen=True)
class RockerBogiePose:
    """A solved rover pose and its wheel-contact diagnostics."""

    success: bool
    chassis_height: float
    chassis_pitch: float
    chassis_roll: float
    left_rocker: float
    right_rocker: float
    left_bogie: float
    right_bogie: float
    wheel_centers: Dict[str, Tuple[float, float]] = field(default_factory=dict)
    contact_residuals: Dict[str, float] = field(default_factory=dict)
    closure_residual: float = 0.0
    errors: Tuple[str, ...] = ()
    iterations: int = 0

    @property
    def max_contact_residual(self) -> float:
        """Largest absolute wheel-height residual, in mm."""
        return max((abs(value) for value in self.contact_residuals.values()),
                   default=0.0)


def _rotate(vector: Tuple[float, float], angle: float) -> np.ndarray:
    x, z = vector
    cosine = math.cos(angle)
    sine = math.sin(angle)
    return np.array((x * cosine - z * sine,
                     x * sine + z * cosine), dtype=float)


def wheel_centers(
    chassis_height: float,
    chassis_pitch: float,
    chassis_roll: float,
    left_rocker: float,
    left_bogie: float,
    right_rocker: float,
    right_bogie: float,
    geometry: RockerBogieGeometry = RockerBogieGeometry(),
) -> Dict[str, Tuple[float, float]]:
    """Compute all planar wheel-center coordinates for a supplied pose.

    ``left_rocker`` and ``right_rocker`` are chassis-relative.  Bogie angles
    are rocker-relative.  Chassis roll raises the left pivot and lowers the
    right pivot by ``half_track * sin(roll)``.
    """
    centers: Dict[str, Tuple[float, float]] = {}
    for side, sign, rocker, bogie in (
        ("l", 1.0, left_rocker, left_bogie),
        ("r", -1.0, right_rocker, right_bogie),
    ):
        # A chassis roll about +X raises the +Y (left) pivot. Link vectors
        # lie in the chassis XZ plane, so their vertical component is scaled
        # by cos(roll) after the pitch/joint rotations. This remains an
        # independent closed-form model while matching the exact 3D rigid
        # transform used by the assembly acceptance test.
        roll_scale = math.cos(chassis_roll)
        center = np.array((0.0, chassis_height +
                           sign * geometry.half_track * math.sin(chassis_roll)))
        rocker_absolute = chassis_pitch + rocker
        rocker_to_bogie = _rotate(geometry.center_to_bogie, rocker_absolute)
        rocker_to_bogie[1] *= roll_scale
        bogie_pivot = center + rocker_to_bogie
        bogie_absolute = rocker_absolute + bogie

        center_to_front = _rotate(geometry.center_to_front, rocker_absolute)
        bogie_to_middle = _rotate(geometry.bogie_to_middle, bogie_absolute)
        bogie_to_rear = _rotate(geometry.bogie_to_rear, bogie_absolute)
        center_to_front[1] *= roll_scale
        bogie_to_middle[1] *= roll_scale
        bogie_to_rear[1] *= roll_scale
        front = center + center_to_front
        middle = bogie_pivot + bogie_to_middle
        rear = bogie_pivot + bogie_to_rear
        centers[f"{side}f"] = (float(front[0]), float(front[1]))
        centers[f"{side}m"] = (float(middle[0]), float(middle[1]))
        centers[f"{side}r"] = (float(rear[0]), float(rear[1]))
    return centers


def solve_terrain_pose(
    terrain_heights: Mapping[str, float],
    geometry: RockerBogieGeometry = RockerBogieGeometry(),
    *,
    tolerance: float = 1e-8,
    max_iterations: int = 80,
) -> RockerBogiePose:
    """Solve the six wheel contacts and ideal differential closure.

    ``terrain_heights`` maps ``lf/lm/lr/rf/rm/rr`` to height above the
    nominal ground plane.  Missing wheels default to zero terrain height.
    The nonlinear solve is unconstrained; mechanical limits are checked on
    the converged pose so an infeasible request reports the travel required.
    """
    unknown = set(terrain_heights) - set(WHEEL_NAMES)
    if unknown:
        names = ", ".join(sorted(unknown))
        return RockerBogiePose(
            False, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            errors=(f"Unknown wheel name(s): {names}",),
        )

    terrain = {name: float(terrain_heights.get(name, 0.0))
               for name in WHEEL_NAMES}
    if not all(math.isfinite(value) for value in terrain.values()):
        return RockerBogiePose(
            False, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            errors=("Terrain heights must be finite",),
        )

    # Variables: height, pitch, roll, left rocker, left bogie, right bogie.
    # Differential closure is exact by construction: right rocker = -left.
    values = np.array((geometry.wheel_radius + 75.0,
                       0.0, 0.0, 0.0, 0.0, 0.0), dtype=float)

    def evaluate(candidate):
        height, pitch, roll, left_rocker, left_bogie, right_bogie = candidate
        centers = wheel_centers(
            height, pitch, roll,
            left_rocker, left_bogie, -left_rocker, right_bogie,
            geometry,
        )
        return np.array([
            centers[name][1] - (geometry.wheel_radius + terrain[name])
            for name in WHEEL_NAMES
        ])

    residual = evaluate(values)
    iterations = 0
    for iterations in range(1, max_iterations + 1):
        if float(np.max(np.abs(residual))) <= tolerance:
            break
        jacobian = np.empty((6, 6), dtype=float)
        for column in range(6):
            step = 1e-5 if column == 0 else 1e-7
            shifted = values.copy()
            shifted[column] += step
            jacobian[:, column] = (evaluate(shifted) - residual) / step
        try:
            delta = np.linalg.solve(jacobian, -residual)
        except np.linalg.LinAlgError:
            break

        # Backtracking keeps Newton steps on the local physical branch.
        old_norm = float(np.max(np.abs(residual)))
        scale = 1.0
        accepted = False
        while scale >= 1.0 / 128.0:
            candidate = values + scale * delta
            candidate_residual = evaluate(candidate)
            if float(np.max(np.abs(candidate_residual))) < old_norm:
                values = candidate
                residual = candidate_residual
                accepted = True
                break
            scale *= 0.5
        if not accepted:
            break

    height, pitch, roll, left_rocker, left_bogie, right_bogie = values
    right_rocker = -left_rocker
    centers = wheel_centers(
        height, pitch, roll, left_rocker, left_bogie,
        right_rocker, right_bogie, geometry,
    )
    residuals = {
        name: centers[name][1] - (geometry.wheel_radius + terrain[name])
        for name in WHEEL_NAMES
    }
    errors = []
    if max(abs(value) for value in residuals.values()) > tolerance:
        errors.append(
            f"Terrain contact solve did not converge within {tolerance:g} mm"
        )
    for side, angle in (("left", left_rocker), ("right", right_rocker)):
        if not geometry.rocker_min <= angle <= geometry.rocker_max:
            errors.append(
                f"{side} rocker requires {math.degrees(angle):.2f} deg; "
                f"limit is {math.degrees(geometry.rocker_min):.2f} to "
                f"{math.degrees(geometry.rocker_max):.2f} deg"
            )
    for side, angle in (("left", left_bogie), ("right", right_bogie)):
        if not geometry.bogie_min <= angle <= geometry.bogie_max:
            errors.append(
                f"{side} bogie requires {math.degrees(angle):.2f} deg; "
                f"limit is {math.degrees(geometry.bogie_min):.2f} to "
                f"{math.degrees(geometry.bogie_max):.2f} deg"
            )

    return RockerBogiePose(
        success=not errors,
        chassis_height=float(height),
        chassis_pitch=float(pitch),
        chassis_roll=float(roll),
        left_rocker=float(left_rocker),
        right_rocker=float(right_rocker),
        left_bogie=float(left_bogie),
        right_bogie=float(right_bogie),
        wheel_centers=centers,
        contact_residuals=residuals,
        closure_residual=float(left_rocker + right_rocker),
        errors=tuple(errors),
        iterations=iterations,
    )


def apply_terrain_pose(
    assembly: Any,
    terrain_heights: Mapping[str, float],
    geometry: RockerBogieGeometry = RockerBogieGeometry(),
    *,
    root_part: str = "chassis",
    left_rocker_joint: str = "left_rocker_pivot",
    right_rocker_joint: str = "right_rocker_pivot",
    left_bogie_joint: str = "left_bogie_pivot",
    right_bogie_joint: str = "right_bogie_pivot",
) -> Tuple[RockerBogiePose, Optional[Any]]:
    """Apply an independently solved terrain pose to a rover assembly.

    The oracle uses nose-up-positive angles, whereas yapCAD's +Y revolute
    axes follow the right-hand rule, so joint and chassis pitch signs are
    inverted at this boundary. Invalid terrain poses do not mutate the
    assembly. The returned second item is the assembly solve result, or
    ``None`` when the terrain request itself is infeasible.
    """
    pose = solve_terrain_pose(terrain_heights, geometry)
    if not pose.success:
        return pose, None

    cosine_pitch = math.cos(pose.chassis_pitch)
    sine_pitch = math.sin(pose.chassis_pitch)
    cosine_roll = math.cos(pose.chassis_roll)
    sine_roll = math.sin(pose.chassis_roll)
    pitch = np.array([
        [cosine_pitch, 0.0, -sine_pitch],
        [0.0, 1.0, 0.0],
        [sine_pitch, 0.0, cosine_pitch],
    ])
    roll = np.array([
        [1.0, 0.0, 0.0],
        [0.0, cosine_roll, -sine_roll],
        [0.0, sine_roll, cosine_roll],
    ])
    root_transform = np.eye(4)
    root_transform[:3, :3] = roll @ pitch
    root_transform[:3, 3] = (0.0, 0.0, pose.chassis_height)

    previous_transforms = {
        name: np.asarray(transform).copy()
        for name, transform in assembly.transforms.items()
    }
    previous_solved = assembly._solved
    previous_root_part = assembly._root_part
    previous_joint_values = assembly._joint_values.copy()
    previous_prescribed_values = assembly._prescribed_joint_values.copy()

    def restore_previous_pose() -> None:
        assembly.transforms = {
            name: transform.copy()
            for name, transform in previous_transforms.items()
        }
        assembly._solved = previous_solved
        assembly._root_part = previous_root_part
        assembly._joint_values = previous_joint_values.copy()
        assembly._prescribed_joint_values = previous_prescribed_values.copy()

    assembly.transforms[root_part] = root_transform
    try:
        result = assembly.solve(root_part, {
            left_rocker_joint: -pose.left_rocker,
            left_bogie_joint: -pose.left_bogie,
            right_bogie_joint: -pose.right_bogie,
        })
    except Exception:
        restore_previous_pose()
        raise
    if not result.success:
        restore_previous_pose()
    else:
        # Make the sign conversion explicit in the acceptance boundary: the
        # dependent right rocker must still be derived by the coupling.
        expected = -pose.right_rocker
        actual = result.joint_values.get(right_rocker_joint)
        if actual is None or abs(actual - expected) > 1e-9:
            restore_previous_pose()
            raise ValueError(
                "Assembly rocker coupling does not match the terrain pose"
            )
    return pose, result


__all__ = [
    "RockerBogieGeometry",
    "RockerBogiePose",
    "WHEEL_NAMES",
    "solve_terrain_pose",
    "apply_terrain_pose",
    "wheel_centers",
]
