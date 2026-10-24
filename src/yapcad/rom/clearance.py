"""SDF-based pairwise clearance between adjacent kinematic links.

Given two adjacent solids ``A`` (parent) and ``B`` (child), each with a
compiled SDF node, evaluates the *minimum signed distance* between them
at a given configuration by:

1. Sampling a coarse grid of points in link B's bounding box.
2. Transforming those points into link A's local frame via the kinematic
   chain's world transforms.
3. Evaluating A's SDF at those points — a negative value means B penetrates
   A; the minimum across all points is the *bottleneck clearance*.

No meshing is needed because the SDF is an analytic field that can be
evaluated at any query point.

Usage::

    from yapcad.rom import clearance_sdf

    node_a = tree.get_sdf_node("THIGH_LINK")
    node_b = tree.get_sdf_node("SHIN_LINK")
    tf_a = tree.get_world_transform("THIGH_LINK")
    tf_b = tree.get_world_transform("SHIN_LINK")

    # Minimum clearance across the interface (negative = interference)
    min_dist, grid = clearance_sdf(node_a, node_b, tf_a, tf_b)
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from yapcad.kinematics import Transform

logger = logging.getLogger(__name__)

# Default grid resolution for SDF clearance check
_DEFAULT_GRID_RESOLUTION = 32  # points per axis → 32k total per pair


def clearance_sdf(
    node_a,
    node_b,
    tf_a: Transform,
    tf_b: Transform,
    resolution: int = _DEFAULT_GRID_RESOLUTION,
    padding: float = 2.0,
) -> Dict[str, Any]:
    """Compute pairwise clearance between two adjacent SDF solids.

    Evaluates B's geometry (sampled point grid) against A's field:
    negative values indicate penetration of B into A.

    Args:
        node_a: SDF ``Node`` for the parent link (or fixed part).
        node_b: SDF ``Node`` for the child link.
        tf_a: World transform of A.
        tf_b: World transform of B.
        resolution: Samples per axis for B's bounding-box grid.
        padding: Extra padding (mm) around B's bounding box.

    Returns:
        ``dict`` with keys:
        - ``min_distance``: Minimum signed distance (mm). Negative = interference.
        - ``penetration_depth``: How deeply B penetrates A (mm), 0 if no penetration.
        - ``mean_distance``: Mean signed distance across all grid points.
        - ``within_5mm``: Fraction of points within 5 mm of A's surface.
        - ``within_2mm``: Fraction within 2 mm.
        - ``fraction_negative``: Fraction of B grid points inside A.
        - ``grid_points``: Total grid points evaluated.
        - ``full_results``: ``(values, axes)`` tuple from ``sdf.sample_grid``.
    """
    from yapcad.sdf import evaluate, sample_grid

    # Get B's bounds in world space, add padding
    bounds_b = node_b.bounds
    if not bounds_b or _is_empty(bounds_b):
        return {
            "min_distance": float("inf"),
            "penetration_depth": 0.0,
            "mean_distance": float("inf"),
            "within_5mm": 0.0,
            "within_2mm": 0.0,
            "fraction_negative": 0.0,
            "grid_points": 0,
            "full_results": None,
        }

    lo_b, hi_b = bounds_b
    # Apply padding
    padded_lo = tuple(v - padding for v in lo_b)
    padded_hi = tuple(v + padding for v in hi_b)
    padded_bounds = (padded_lo, padded_hi)

    # NB: The bounds above are in B's *local* coordinate system.
    # We need to sample *in world space* because the SDF node_b is
    # defined in B's local frame (it was compiled from the DSL without
    # any world transform applied).  Since the SDF is not aware of the
    # kinematic world transform, we must transform the query points.

    # Strategy: sample a grid in B-local space, transform those points
    # to world space, then transform from world to A's local frame
    # to evaluate node_a.

    # Build local grid for B
    if _is_infinite_bounds(padded_bounds):
        # Fallback: use B's finite sub-regions
        logger.warning("Node B has infinite bounds; using fallback bounds")
        for axis in range(3):
            if abs(padded_lo[axis]) > 1e6 or abs(padded_hi[axis]) > 1e6:
                padded_lo = list(padded_lo)
                padded_hi = list(padded_hi)
                padded_lo[axis] = -200.0
                padded_hi[axis] = 200.0
                padded_bounds = (tuple(padded_lo), tuple(padded_hi))
                break

    values, axes = sample_grid(
        node_b,
        bounds=padded_bounds,
        resolution=resolution,
    )

    grid_points = values.size
    if grid_points == 0:
        return {
            "min_distance": float("inf"),
            "penetration_depth": 0.0,
            "mean_distance": float("inf"),
            "within_5mm": 0.0,
            "within_2mm": 0.0,
            "fraction_negative": 0.0,
            "grid_points": 0,
            "full_results": None,
        }

    # Build world-frame grid of B's surface region
    mask = np.abs(values) < 10.0  # Only consider points near B's surface
    if not mask.any():
        mask = np.ones_like(values, dtype=bool)  # Fallback: use all

    # Get the surface-adjacent points in B's local frame
    xs, ys, zs = axes
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    local_b_points = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
    surface_b = local_b_points[mask.ravel()]

    if len(surface_b) == 0:
        surface_b = local_b_points

    # Convert to homogeneous coordinates for transform
    n_pts = len(surface_b)
    homo_b = np.ones((n_pts, 4), dtype=float)
    homo_b[:, :3] = surface_b

    # B-local → world
    tf_b_mat = tf_b.matrix
    world_pts = (tf_b_mat @ homo_b.T).T[:, :3]

    # World → A-local
    tf_a_inv = tf_a.inverse()
    tf_a_inv_mat = tf_a_inv.matrix
    homo_world = np.ones((n_pts, 4), dtype=float)
    homo_world[:, :3] = world_pts
    local_a_pts = (tf_a_inv_mat @ homo_world.T).T[:, :3]

    # Evaluate A's SDF at the transformed points
    from yapcad.sdf import evaluate as sdf_eval
    distances = sdf_eval(node_a, local_a_pts)

    # Statistics
    min_dist = float(distances.min()) if len(distances) > 0 else float("inf")
    mean_dist = float(distances.mean()) if len(distances) > 0 else float("inf")
    frac_neg = float((distances < 0).mean()) if len(distances) > 0 else 0.0
    frac_5mm = float((distances < 5.0).mean()) if len(distances) > 0 else 0.0
    frac_2mm = float((distances < 2.0).mean()) if len(distances) > 0 else 0.0

    penetration = -min_dist if min_dist < 0 else 0.0

    return {
        "min_distance": min_dist,
        "penetration_depth": penetration,
        "mean_distance": mean_dist,
        "within_5mm": frac_5mm,
        "within_2mm": frac_2mm,
        "fraction_negative": frac_neg,
        "grid_points": n_pts,
        "full_results": (distances, surface_b),
    }


def clearance_mesh(
    mesh_a,
    mesh_b,
    tf_a: Transform,
    tf_b: Transform,
    n_samples: int = 5000,
) -> Dict[str, Any]:
    """Mesh-based clearance when SDF nodes are unavailable.

    Falls back to sampling surface points from mesh B and checking
    containment in mesh A (through trimesh).

    Args:
        mesh_a: trimesh.Trimesh or path to STL for parent link.
        mesh_b: trimesh.Trimesh or path to STL for child link.
        tf_a: World transform of A.
        tf_b: World transform of B.
        n_samples: Number of surface sample points per mesh.

    Returns:
        ``dict`` with same key shape as :func:`clearance_sdf`.
    """
    import trimesh
    from scipy.spatial import cKDTree

    # Load meshes
    if isinstance(mesh_a, str):
        mesh_a = trimesh.load(mesh_a, force="mesh")
    if isinstance(mesh_b, str):
        mesh_b = trimesh.load(mesh_b, force="mesh")

    ma = mesh_a.copy()
    mb = mesh_b.copy()
    ma.apply_transform(tf_a.matrix)
    mb.apply_transform(tf_b.matrix)

    # Sample surface points
    pts_a, _ = trimesh.sample.sample_surface(ma, n_samples)
    pts_b, _ = trimesh.sample.sample_surface(mb, n_samples)

    # Build KD-trees
    tree_a = cKDTree(pts_a)
    tree_b = cKDTree(pts_b)

    da, _ = tree_b.query(pts_a)
    db, _ = tree_a.query(pts_b)

    min_dist = min(float(da.min()), float(db.min()))
    mean_dist = float(np.mean(np.concatenate([da, db])))
    frac_5mm = float((da < 5.0).mean() + (db < 5.0).mean()) / 2.0
    frac_2mm = float((da < 2.0).mean() + (db < 2.0).mean()) / 2.0

    penetration = -min_dist if min_dist < 0 else 0.0

    return {
        "min_distance": min_dist,
        "penetration_depth": penetration,
        "mean_distance": mean_dist,
        "within_5mm": frac_5mm,
        "within_2mm": frac_2mm,
        "fraction_negative": float((da < 0).mean()) if min_dist < 0 else 0.0,
        "grid_points": 2 * n_samples,
        "full_results": None,
    }


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _is_empty(bounds) -> bool:
    lo, hi = bounds
    for a, b in zip(lo, hi):
        if b <= a:
            return True
    return False


def _is_infinite_bounds(bounds) -> bool:
    lo, hi = bounds
    for v in list(lo) + list(hi):
        if abs(v) > 1e6:
            return True
    return False


def _repair_bounds(bounds) -> Tuple[Tuple, Tuple]:
    """Replace infinite/unbounded axes with a finite fallback."""
    lo, hi = list(bounds[0]), list(bounds[1])
    for i in range(3):
        if abs(lo[i]) > 1e6 or abs(hi[i]) > 1e6:
            lo[i] = -200.0
            hi[i] = 200.0
    return (tuple(lo), tuple(hi))


__all__ = [
    "clearance_sdf",
    "clearance_mesh",
]