"""Deterministic relationships between named assembly joint coordinates."""

from dataclasses import dataclass, field
from typing import Dict
import math


@dataclass(frozen=True)
class LinearJointCoupling:
    """Define one joint coordinate as an affine function of driver joints.

    The relation is ``dependent = offset + sum(coefficient * driver)``.
    Naming the dependent coordinate makes evaluation deterministic and permits
    coupling chains to be resolved topologically.
    """

    name: str
    dependent_joint: str
    driver_coefficients: Dict[str, float] = field(default_factory=dict)
    offset: float = 0.0
    tolerance: float = 1e-9

    def __post_init__(self):
        object.__setattr__(
            self, "driver_coefficients",
            {str(name): float(value)
             for name, value in self.driver_coefficients.items()},
        )
        if not self.name:
            raise ValueError("Joint coupling name must not be empty")
        if not self.dependent_joint:
            raise ValueError("Dependent joint name must not be empty")
        if not self.driver_coefficients:
            raise ValueError("Linear joint coupling requires at least one driver")
        if self.dependent_joint in self.driver_coefficients:
            raise ValueError("A joint coupling cannot depend on itself")
        values = list(self.driver_coefficients.values()) + [self.offset, self.tolerance]
        if not all(math.isfinite(value) for value in values):
            raise ValueError("Joint coupling coefficients must be finite")
        if self.tolerance <= 0.0:
            raise ValueError("Joint coupling tolerance must be positive")

    def evaluate(self, joint_values: Dict[str, float]) -> float:
        return float(self.offset + sum(
            coefficient * joint_values[joint]
            for joint, coefficient in self.driver_coefficients.items()
        ))

    def residual(self, joint_values: Dict[str, float]) -> float:
        return abs(joint_values[self.dependent_joint] - self.evaluate(joint_values))

    def to_dict(self) -> Dict[str, object]:
        """Return the stable semantic representation used by graph export."""
        return {
            "id": self.name,
            "type": "Affine",
            "dependent_joint": self.dependent_joint,
            "driver_coefficients": dict(self.driver_coefficients),
            "offset": float(self.offset),
            "tolerance": float(self.tolerance),
        }
