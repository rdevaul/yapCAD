"""4x4 Homogeneous Transformation Matrix.

This module provides the Transform class for representing rigid body
transformations in 3D space using 4x4 homogeneous matrices.

Copyright (c) 2026 yapCAD contributors
License: MIT
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Tuple, Optional, Any, List

import numpy as np


@dataclass
class Transform:
    """4x4 homogeneous transformation matrix.

    Represents a rigid body transformation (rotation + translation) in 3D space.
    The matrix format is::

        [R11 R12 R13 Tx]
        [R21 R22 R23 Ty]
        [R31 R32 R33 Tz]
        [0   0   0   1 ]

    Where the upper-left 3x3 is the rotation matrix and the right column
    (top 3 elements) is the translation vector.

    Transforms can be composed via matrix multiplication using the @ operator::

        world_T_tool = world_T_base @ base_T_link @ link_T_tool

    Attributes:
        matrix: 4x4 numpy array representing the transform

    Example::

        # Create transforms
        t1 = Transform.from_translation(10, 0, 0)
        t2 = Transform.from_rotation_z(45)  # degrees

        # Compose (t2 applied after t1)
        t3 = t1 @ t2

        # Get components
        pos = t3.translation  # (x, y, z)
        R = t3.rotation_matrix  # 3x3
    """

    matrix: np.ndarray = field(default_factory=lambda: np.eye(4, dtype=float))

    def __post_init__(self):
        """Ensure matrix is 4x4 numpy array."""
        if not isinstance(self.matrix, np.ndarray):
            self.matrix = np.array(self.matrix, dtype=float)
        if self.matrix.shape != (4, 4):
            raise ValueError(f"Transform matrix must be 4x4, got {self.matrix.shape}")

    # =========================================================================
    # Factory Methods
    # =========================================================================

    @classmethod
    def identity(cls) -> Transform:
        """Create identity transform (no rotation, no translation)."""
        return cls(np.eye(4, dtype=float))

    @classmethod
    def from_translation(cls, x: float, y: float, z: float) -> Transform:
        """Create translation-only transform.

        :param x: Translation along X axis
        :param y: Translation along Y axis
        :param z: Translation along Z axis
        :returns: Transform with specified translation
        """
        m = np.eye(4, dtype=float)
        m[0, 3] = x
        m[1, 3] = y
        m[2, 3] = z
        return cls(m)

    @classmethod
    def from_rotation_x(cls, angle_deg: float) -> Transform:
        """Create rotation about X axis.

        :param angle_deg: Rotation angle in degrees
        :returns: Transform with rotation about X
        """
        rad = math.radians(angle_deg)
        c, s = math.cos(rad), math.sin(rad)
        m = np.eye(4, dtype=float)
        m[1, 1] = c
        m[1, 2] = -s
        m[2, 1] = s
        m[2, 2] = c
        return cls(m)

    @classmethod
    def from_rotation_y(cls, angle_deg: float) -> Transform:
        """Create rotation about Y axis.

        :param angle_deg: Rotation angle in degrees
        :returns: Transform with rotation about Y
        """
        rad = math.radians(angle_deg)
        c, s = math.cos(rad), math.sin(rad)
        m = np.eye(4, dtype=float)
        m[0, 0] = c
        m[0, 2] = s
        m[2, 0] = -s
        m[2, 2] = c
        return cls(m)

    @classmethod
    def from_rotation_z(cls, angle_deg: float) -> Transform:
        """Create rotation about Z axis.

        :param angle_deg: Rotation angle in degrees
        :returns: Transform with rotation about Z
        """
        rad = math.radians(angle_deg)
        c, s = math.cos(rad), math.sin(rad)
        m = np.eye(4, dtype=float)
        m[0, 0] = c
        m[0, 1] = -s
        m[1, 0] = s
        m[1, 1] = c
        return cls(m)

    @classmethod
    def from_euler_xyz(
        cls, rx_deg: float, ry_deg: float, rz_deg: float
    ) -> Transform:
        """Create rotation from Euler angles (intrinsic XYZ order).

        :param rx_deg: Rotation about X axis in degrees
        :param ry_deg: Rotation about Y axis in degrees
        :param rz_deg: Rotation about Z axis in degrees
        :returns: Transform with combined rotation
        """
        Rx = cls.from_rotation_x(rx_deg)
        Ry = cls.from_rotation_y(ry_deg)
        Rz = cls.from_rotation_z(rz_deg)
        return Rx @ Ry @ Rz

    @classmethod
    def from_axis_angle(
        cls,
        axis: Tuple[float, float, float],
        angle_deg: float,
    ) -> Transform:
        """Create rotation from axis-angle representation.

        Uses Rodrigues' rotation formula.

        :param axis: Unit rotation axis (x, y, z)
        :param angle_deg: Rotation angle in degrees
        :returns: Transform with specified rotation
        """
        # Normalize axis
        ax = np.array(axis, dtype=float)
        norm = np.linalg.norm(ax)
        if norm < 1e-10:
            return cls.identity()
        ax = ax / norm

        rad = math.radians(angle_deg)
        c, s = math.cos(rad), math.sin(rad)

        # Rodrigues formula: R = I + sin(θ)K + (1-cos(θ))K²
        K = np.array([
            [0, -ax[2], ax[1]],
            [ax[2], 0, -ax[0]],
            [-ax[1], ax[0], 0]
        ], dtype=float)

        R = np.eye(3) + s * K + (1 - c) * (K @ K)

        m = np.eye(4, dtype=float)
        m[:3, :3] = R
        return cls(m)

    @classmethod
    def from_position_rpy(
        cls,
        x: float, y: float, z: float,
        roll_deg: float, pitch_deg: float, yaw_deg: float,
    ) -> Transform:
        """Create transform from position and roll-pitch-yaw angles.

        :param x: X position
        :param y: Y position
        :param z: Z position
        :param roll_deg: Roll angle (rotation about X) in degrees
        :param pitch_deg: Pitch angle (rotation about Y) in degrees
        :param yaw_deg: Yaw angle (rotation about Z) in degrees
        :returns: Combined transform
        """
        T = cls.from_translation(x, y, z)
        R = cls.from_euler_xyz(roll_deg, pitch_deg, yaw_deg)
        return T @ R

    @classmethod
    def from_matrix(cls, matrix: Any) -> Transform:
        """Create transform from 4x4 matrix (numpy array or nested list).

        :param matrix: 4x4 array-like
        :returns: Transform wrapping the matrix
        """
        return cls(np.array(matrix, dtype=float))

    # =========================================================================
    # Properties
    # =========================================================================

    @property
    def translation(self) -> Tuple[float, float, float]:
        """Get translation component as (x, y, z) tuple."""
        return (
            float(self.matrix[0, 3]),
            float(self.matrix[1, 3]),
            float(self.matrix[2, 3]),
        )

    @property
    def rotation_matrix(self) -> np.ndarray:
        """Get 3x3 rotation matrix."""
        return self.matrix[:3, :3].copy()

    # =========================================================================
    # Operations
    # =========================================================================

    def __matmul__(self, other: Transform) -> Transform:
        """Compose transforms via matrix multiplication.

        :param other: Transform to compose with
        :returns: New transform representing self @ other

        The result applies `other` first, then `self`.
        """
        if not isinstance(other, Transform):
            raise TypeError(f"Cannot compose Transform with {type(other)}")
        return Transform(self.matrix @ other.matrix)

    def inverse(self) -> Transform:
        """Compute inverse transform.

        :returns: Transform T^-1 such that T @ T^-1 = I

        For rigid transforms, this is computed efficiently as::

            R^T | -R^T * t
            ----+----------
            0   |    1
        """
        R = self.matrix[:3, :3]
        t = self.matrix[:3, 3]

        R_inv = R.T
        t_inv = -R_inv @ t

        m = np.eye(4, dtype=float)
        m[:3, :3] = R_inv
        m[:3, 3] = t_inv
        return Transform(m)

    def transform_point(self, point: Tuple[float, float, float]) -> Tuple[float, float, float]:
        """Transform a point (applies rotation + translation).

        :param point: Point as (x, y, z)
        :returns: Transformed point as (x, y, z)
        """
        p = np.array([point[0], point[1], point[2], 1.0])
        result = self.matrix @ p
        return (float(result[0]), float(result[1]), float(result[2]))

    def transform_vector(self, vector: Tuple[float, float, float]) -> Tuple[float, float, float]:
        """Transform a direction vector (rotation only, no translation).

        :param vector: Direction vector as (x, y, z)
        :returns: Transformed vector as (x, y, z)
        """
        v = np.array([vector[0], vector[1], vector[2], 0.0])
        result = self.matrix @ v
        return (float(result[0]), float(result[1]), float(result[2]))

    def to_euler_xyz(self) -> Tuple[float, float, float]:
        """Extract Euler angles (intrinsic XYZ order) in degrees.

        :returns: (rx, ry, rz) rotation angles in degrees

        Note: May have gimbal lock issues when pitch is near ±90°.
        """
        R = self.matrix[:3, :3]

        # Extract angles using atan2 for robustness
        sy = math.sqrt(R[0, 0] ** 2 + R[1, 0] ** 2)

        if sy > 1e-6:
            rx = math.atan2(R[2, 1], R[2, 2])
            ry = math.atan2(-R[2, 0], sy)
            rz = math.atan2(R[1, 0], R[0, 0])
        else:
            # Gimbal lock
            rx = math.atan2(-R[1, 2], R[1, 1])
            ry = math.atan2(-R[2, 0], sy)
            rz = 0

        return (
            math.degrees(rx),
            math.degrees(ry),
            math.degrees(rz),
        )

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization.

        :returns: Dict with translation, rotation_xyz, and matrix
        """
        tx, ty, tz = self.translation
        rx, ry, rz = self.to_euler_xyz()
        return {
            "translation": [tx, ty, tz],
            "rotation_xyz": [rx, ry, rz],
            "matrix": self.matrix.tolist(),
        }

    @classmethod
    def from_dict(cls, data: dict) -> Transform:
        """Create transform from dictionary.

        :param data: Dict with 'matrix' key (4x4 list)
        :returns: Transform from matrix
        """
        if "matrix" in data:
            return cls.from_matrix(data["matrix"])
        elif "translation" in data and "rotation_xyz" in data:
            t = data["translation"]
            r = data["rotation_xyz"]
            return cls.from_position_rpy(t[0], t[1], t[2], r[0], r[1], r[2])
        else:
            raise ValueError("Dict must contain 'matrix' or 'translation'+'rotation_xyz'")

    def __repr__(self) -> str:
        tx, ty, tz = self.translation
        rx, ry, rz = self.to_euler_xyz()
        return (
            f"Transform(t=[{tx:.2f}, {ty:.2f}, {tz:.2f}], "
            f"r=[{rx:.1f}°, {ry:.1f}°, {rz:.1f}°])"
        )


__all__ = ["Transform"]
