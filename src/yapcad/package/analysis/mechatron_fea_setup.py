"""Mechatron-canonical FEA setup helper.

Reads bolt_patterns and load_cases from a Mechatron `graph.json` (the
canonical source-of-truth) and produces structured inputs the FEA solver
needs: per-bolt world coordinates, spring stiffness lookup, load case
attach resolution.

Replaces the v4-v7 pattern of hand-rolled `bolt_inventory.json` +
`canonical.json` side files. After this lands, a FEA script does::

    from yapcad.package.analysis.mechatron_fea_setup import prepare
    setup = prepare(graph_path="path/to/assembly/graph.json",
                    load_case_id="LC-004a")
    # setup.bolts:        list of Bolt entries (world coords + axis + spec)
    # setup.spring_k:     {(bolt_idx): (k_axial_N_per_m, k_shear_N_per_m)}
    # setup.load_attach:  resolved attach info (part, position_mm, direction_unit)
    # setup.load_case:    raw LoadCase dict for reference

History: created 2026-05-20 as Phase 3+4 of the FEA → Mechatron
integration plan.
"""
from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple


# ---------------------------------------------------------------------------
# Bolt-spec stiffness library
# ---------------------------------------------------------------------------
#
# Maps canonical bolt_spec strings → (k_axial_N_per_m, k_shear_N_per_m)
# for use in v7+ bolt-spring MPC formulation.
#
# Stiffness is calculated as k = E*A/L for axial, and k_shear ≈ G*A/L
# with G ≈ 0.4*E for steel. Bolt clamping is through PETG-CF flanges
# (~4-6 mm thick each); typical L ≈ 12 mm grip.
#
# These values can be over-ridden per load case via the graph.json
# `interfaces[].bolt_pattern.stiffness_override_N_per_m` field if needed.

_E_STEEL_8_8 = 200e9  # Pa
_G_STEEL_8_8 = 80e9   # Pa
_TYPICAL_GRIP_L_M = 0.012  # 12 mm = ~6 mm + ~6 mm flange stack-up

def _bolt_area_m2(diameter_in_mm: float) -> float:
    """Cross-section area of a bolt shaft in m²."""
    return math.pi * (diameter_in_mm * 1e-3) ** 2 / 4.0

# Canonical name → (shaft diameter mm, length mm, grade)
_BOLT_SPEC_GEOMETRY = {
    "1-4-20-x-1in-button-head-shcs":      (6.35, 25.4, "8.8"),
    "1-4-20-x-1-25in-button-head-shcs":   (6.35, 31.75, "8.8"),
    "1-4-20-button-head-shcs":            (6.35, 25.4, "8.8"),  # default length
    "m5-x-16-shcs":                        (5.00, 16.0, "8.8"),
    "m5-x-20-shcs":                        (5.00, 20.0, "8.8"),
    "m6-x-20-shcs":                        (6.00, 20.0, "8.8"),
}

def stiffness_for_bolt_spec(bolt_spec: str) -> Tuple[float, float]:
    """Return (k_axial_N_per_m, k_shear_N_per_m) for a bolt spec.

    Uses the canonical _BOLT_SPEC_GEOMETRY table. Unknown specs fall
    back to M5 stiffness with a warning.
    """
    if bolt_spec.startswith("PLACEHOLDER"):
        # Placeholder spec — return zero stiffness so the FEA setup
        # treats this bolt as "not yet specified" and skips it.
        return (0.0, 0.0)
    geo = _BOLT_SPEC_GEOMETRY.get(bolt_spec)
    if geo is None:
        # Default to M5 if unknown — log via stderr (no hard requirement on logging here)
        import sys
        sys.stderr.write(
            f"[mechatron_fea_setup] unknown bolt_spec {bolt_spec!r}, defaulting to M5\n"
        )
        geo = (5.0, 20.0, "8.8")
    d_mm, l_mm, grade = geo
    A = _bolt_area_m2(d_mm)
    L = _TYPICAL_GRIP_L_M
    k_axial = _E_STEEL_8_8 * A / L
    k_shear = _G_STEEL_8_8 * A / L
    return (k_axial, k_shear)


# ---------------------------------------------------------------------------
# Result dataclasses
# ---------------------------------------------------------------------------

@dataclass
class Bolt:
    """One bolt in the assembly, world coordinates.

    `parent_world` / `child_world` give the head- and nut-side
    positions; the bolt axis is the unit vector from parent to child
    (or `access_direction` if available).
    """
    interface_id: str
    bolt_index: int
    parent_part: str
    child_part: str
    parent_world_mm: Tuple[float, float, float]
    child_world_mm: Tuple[float, float, float]
    axis_unit: Tuple[float, float, float]
    bolt_spec: str
    pcd_mm: float
    clock_deg: float  # 0-360°, computed from parent_world
    k_axial_N_per_m: float
    k_shear_N_per_m: float


@dataclass
class LoadAttachResolved:
    """A LoadCase resolved against the assembly: where the force is applied,
    in the FEA mesh frame, on which part(s)."""
    load_case_id: str
    part: str  # "_assembly" for whole-stack, else specific part_id
    interface_id: Optional[str]
    position_mm: Tuple[float, float, float]  # vehicle frame (X = axial-fwd)
    direction_unit: Tuple[float, float, float]  # in `coordinate_frame`
    magnitude_n: float
    coordinate_frame: str  # "vehicle" or "mesh"
    bolt_index: Optional[int]
    clock_deg: Optional[float]
    raw: Dict[str, Any]  # full LoadCase dict for reference


@dataclass
class FeaSetup:
    """Everything an FEA script needs, sourced from mechatron graph.json."""
    graph_path: str
    bolts: List[Bolt] = field(default_factory=list)
    load_case_id: Optional[str] = None
    load_attach: Optional[LoadAttachResolved] = None
    raw_load_case: Optional[Dict[str, Any]] = None
    raw_interfaces: List[Dict[str, Any]] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Coordinate-frame transforms
# ---------------------------------------------------------------------------
#
# Dashboard / vehicle frame:  X = axial-fwd, Y = per-RH-rule, Z = launch-rail-side
# FEA mesh frame:              Z = axial-fwd, X/Y = radial
#
# Transform: (x_dash, y_dash, z_dash) → (y_dash, z_dash, x_dash)

def vehicle_to_mesh_position(pos: Sequence[float]) -> Tuple[float, float, float]:
    """Vehicle frame (X=axial) → FEA mesh frame (Z=axial)."""
    return (float(pos[1]), float(pos[2]), float(pos[0]))


def vehicle_to_mesh_direction(direction: Sequence[float]) -> Tuple[float, float, float]:
    """Same swap for direction vectors."""
    return (float(direction[1]), float(direction[2]), float(direction[0]))


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def load_graph(graph_path: str | Path) -> Dict[str, Any]:
    """Load and validate a mechatron graph.json."""
    p = Path(graph_path).expanduser()
    if not p.exists():
        raise FileNotFoundError(f"graph.json not found: {p}")
    g = json.loads(p.read_text())
    if "interfaces" not in g:
        raise ValueError(f"{p}: missing 'interfaces' field — not a valid graph.json")
    return g


def list_load_cases(graph: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Return all LoadCase dicts from the graph."""
    return list(graph.get("load_cases", []))


def get_load_case(graph: Dict[str, Any], lc_id: str) -> Dict[str, Any]:
    """Return one LoadCase by id, or raise KeyError."""
    for lc in graph.get("load_cases", []):
        if lc.get("id") == lc_id:
            return lc
    raise KeyError(f"LoadCase {lc_id!r} not in graph (have: {[lc.get('id') for lc in graph.get('load_cases', [])]})")


def get_interface(graph: Dict[str, Any], iface_id: str) -> Dict[str, Any]:
    """Return one Interface by id, or raise KeyError."""
    for iface in graph.get("interfaces", []):
        if iface.get("id") == iface_id:
            return iface
    raise KeyError(f"Interface {iface_id!r} not in graph")


def list_interfaces_with_bolt_patterns(graph: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Return only interfaces that have a bolt_pattern attached.

    Skips bolt_patterns with PLACEHOLDER bolt_spec — those are flagged
    as needing engineering input before being usable.
    """
    out = []
    for iface in graph.get("interfaces", []):
        bp = iface.get("bolt_pattern")
        if not bp:
            continue
        if bp.get("bolt_spec", "").startswith("PLACEHOLDER"):
            continue
        out.append(iface)
    return out


# ---------------------------------------------------------------------------
# Per-part bolt spec loader (2026-05-27 rework)
# ---------------------------------------------------------------------------
#
# The canonical source of truth for per-bolt orientation moved from
# graph.json's `bolt_pattern.access_direction` (which is single-vector
# per interface and was sloppily set to ±Z for every joint) to
# `bolts_per_part.yaml`, which carries per-bolt-group `orientation:
# axial|radial` plus per-hole datums on the owner part.
#
# `_load_bolts_per_part_spec(graph_path)` looks for the YAML next to
# graph.json (it lives at <designs>/assembly/bolts_per_part.yaml while
# graph.json lives at <designs>/assembly/graph.json/graph.json), parses
# it, and returns the inner `bolts_per_part` dict (or {} when absent).
#
# The generator then prefers this spec when available, otherwise falls
# back to the legacy access_direction behavior so older graphs still
# work.

def _load_bolts_per_part_spec(graph_path: str | Path) -> Dict[str, Any]:
    """Locate and load bolts_per_part.yaml relative to graph_path.

    Search order (first hit wins):
      1. <graph_dir>/bolts_per_part.yaml
      2. <graph_dir>/../bolts_per_part.yaml  (graph_dir is graph.json/, parent is assembly/)
      3. <graph_dir>/../../bolts_per_part.yaml

    Returns the inner ``bolts_per_part`` mapping (owner-part-id -> group-name
    -> group-dict). Returns {} when the file is missing, PyYAML is not
    available, or the file fails to parse.
    """
    p = Path(graph_path).expanduser().resolve()
    candidates = [
        p.parent / "bolts_per_part.yaml",
        p.parent.parent / "bolts_per_part.yaml",
        p.parent.parent.parent / "bolts_per_part.yaml",
    ]
    spec_path = None
    for c in candidates:
        if c.is_file():
            spec_path = c
            break
    if spec_path is None:
        return {}
    try:
        import yaml  # type: ignore
    except ImportError:
        import sys
        sys.stderr.write(
            "[mechatron_fea_setup] PyYAML not installed; bolts_per_part.yaml ignored\n"
        )
        return {}
    try:
        data = yaml.safe_load(spec_path.read_text(encoding="utf-8")) or {}
    except Exception as e:
        import sys
        sys.stderr.write(f"[mechatron_fea_setup] failed to parse {spec_path}: {e}\n")
        return {}
    return data.get("bolts_per_part", {}) or {}


# ---------------------------------------------------------------------------
# Kinematic chain → world transforms for every part in the graph
# ---------------------------------------------------------------------------

def _qmul(a: Sequence[float], b: Sequence[float]) -> List[float]:
    """Hamilton product, quaternions in (x, y, z, w) order."""
    ax, ay, az, aw = a
    bx, by, bz, bw = b
    return [
        aw * bx + ax * bw + ay * bz - az * by,
        aw * by - ax * bz + ay * bw + az * bx,
        aw * bz + ax * by - ay * bx + az * bw,
        aw * bw - ax * bx - ay * by - az * bz,
    ]


def _qrot(q: Sequence[float], v: Sequence[float]) -> List[float]:
    """Rotate 3-vector v by unit quaternion q in (x, y, z, w) form."""
    x, y, z, w = q[0], q[1], q[2], q[3]
    vx, vy, vz = v[0], v[1], v[2]
    tx = 2.0 * (y * vz - z * vy)
    ty = 2.0 * (z * vx - x * vz)
    tz = 2.0 * (x * vy - y * vx)
    return [
        vx + w * tx + (y * tz - z * ty),
        vy + w * ty + (z * tx - x * tz),
        vz + w * tz + (x * ty - y * tx),
    ]


def solve_world_transforms(graph: Dict[str, Any]) -> Dict[str, Tuple[List[float], List[float]]]:
    """Walk the kinematic chain from `world` and compute world (T, Q) for every part.

    Returns: dict mapping part_id -> ([tx,ty,tz], [qx,qy,qz,qw]).

    Mirrors the algorithm in the dashboard's `_load_graph()`. Interfaces
    carry `origin: [x,y,z]` (child position in PARENT-LOCAL frame) and
    optional `orientation: [x,y,z,w]` (child rotation in parent frame).
    Root is the synthetic ``world`` part at identity.
    """
    parent_to_children: Dict[str, List[Tuple[str, List[float], Optional[List[float]]]]] = {}
    for iface in graph.get("interfaces", []):
        parent = iface.get("parent_part", "")
        child = iface.get("child_part", "")
        if not parent or not child:
            continue
        origin = iface.get("origin") or [0.0, 0.0, 0.0]
        try:
            rel_xyz = [float(origin[0]), float(origin[1]), float(origin[2])]
        except (TypeError, ValueError, IndexError):
            rel_xyz = [0.0, 0.0, 0.0]
        orient = iface.get("orientation")
        rel_q: Optional[List[float]]
        if isinstance(orient, list) and len(orient) == 4:
            try:
                rel_q = [float(orient[0]), float(orient[1]), float(orient[2]), float(orient[3])]
            except (TypeError, ValueError):
                rel_q = None
        else:
            rel_q = None
        parent_to_children.setdefault(parent, []).append((child, rel_xyz, rel_q))

    out: Dict[str, Tuple[List[float], List[float]]] = {
        "world": ([0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0]),
    }
    queue: List[str] = ["world"]
    while queue:
        cur = queue.pop(0)
        ct, cq = out[cur]
        for child, rel_xyz, rel_q in parent_to_children.get(cur, []):
            if child in out:
                continue
            rotated = _qrot(cq, rel_xyz)
            child_t = [ct[0] + rotated[0], ct[1] + rotated[1], ct[2] + rotated[2]]
            child_q = _qmul(cq, rel_q) if rel_q else list(cq)
            out[child] = (child_t, child_q)
            queue.append(child)
    return out


# ---------------------------------------------------------------------------
# Bolt generation — per-hole datum aware
# ---------------------------------------------------------------------------

def _index_parts(graph: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    return {p.get("id", ""): p for p in graph.get("parts", [])}


def _find_datum(part: Dict[str, Any], name: str) -> Optional[Dict[str, Any]]:
    for d in part.get("datums", []):
        if d.get("name") == name:
            return d
    return None


def _datum_origin(d: Dict[str, Any]) -> List[float]:
    o = d.get("origin") or {}
    return [float(o.get("x", 0.0)), float(o.get("y", 0.0)), float(o.get("z", 0.0))]


def _datum_direction(d: Dict[str, Any]) -> Optional[List[float]]:
    dd = d.get("direction")
    if not isinstance(dd, dict):
        return None
    v = [float(dd.get("x", 0.0)), float(dd.get("y", 0.0)), float(dd.get("z", 0.0))]
    n = math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
    if n < 1e-9:
        return None
    return [v[0] / n, v[1] / n, v[2] / n]


def _find_yaml_group_for_iface(
    spec: Dict[str, Any],
    iface: Dict[str, Any],
) -> Optional[Tuple[str, str, Dict[str, Any]]]:
    """Look up the bolts_per_part.yaml group matching this interface.

    Returns (owner_part_id, group_name, group_spec_dict) or None when no
    match is found. The owner can be either the parent or the child of
    the interface (per the YAML ownership rule). Matching is by
    ``mates_to.part == <the other side>`` since the same owner-part may
    own multiple bolt groups (e.g. upper-thrust has `bot_to_boattail`
    and `top_to_fuel_tank`).
    """
    if not spec:
        return None
    pp = iface.get("parent_part", "")
    cp = iface.get("child_part", "")

    def _match(owner: str, mated: str) -> Optional[Tuple[str, str, Dict[str, Any]]]:
        owner_block = spec.get(owner) or {}
        if not isinstance(owner_block, dict):
            return None
        for gname, g in owner_block.items():
            if not isinstance(g, dict):
                continue
            mt = g.get("mates_to") or {}
            if mt.get("part") == mated:
                return (owner, gname, g)
        return None

    # Owner could be either side. Try parent-as-owner first then child-as-owner.
    return _match(pp, cp) or _match(cp, pp)


def generate_bolts_for_interface(
    iface: Dict[str, Any],
    parent_world_origin_mm: Optional[Tuple[float, float, float]] = None,
    child_world_origin_mm: Optional[Tuple[float, float, float]] = None,
    *,
    parts_by_id: Optional[Dict[str, Dict[str, Any]]] = None,
    world_transforms: Optional[Dict[str, Tuple[List[float], List[float]]]] = None,
    bolts_per_part_spec: Optional[Dict[str, Any]] = None,
) -> List[Bolt]:
    """Generate per-bolt entries for an interface's bolt_pattern.

    Two code paths:

    1. **Preferred** (post-2026-05-27): when ``parts_by_id`` +
       ``world_transforms`` are supplied, derive each bolt's axis from
       the parent part's per-hole datums (``top_bolt_hole_N`` /
       ``bot_bolt_hole_N``) transformed through the kinematic chain.
       The per-hole ``direction`` field in graph.json carries the true
       radial-or-axial bolt axis per hole, so this path produces
       correct orientations for both joint types.

       When ``bolts_per_part_spec`` is also supplied AND a matching
       group exists, additionally use the YAML's ``orientation``,
       ``head_axial_offset_mm`` / ``head_radial_offset_mm``, and
       ``shank_length_mm`` to derive the head-anchor and shank-tip
       endpoints. ``parent_world_mm`` becomes the head anchor (where
       the bolt enters the joint), ``child_world_mm`` becomes the
       shank tip (where it exits), and the connector line spans the
       physical bolt length.

    2. **Legacy fallback** (used when ``parts_by_id`` is absent):
       evenly distribute bolts around PCD using
       ``bolt_pattern.access_direction`` and ``world_z_mm``. Same
       behavior as the pre-2026-05-27 code. Note that this path emits
       all axes as ±Z, since it distributes bolts purely by PCD +
       ``access_direction`` and has no per-hole radial datum data. Use
       the preferred per-hole-datum path (supply ``parts_by_id`` +
       ``world_transforms``) when correct radial-joint axes are needed.

    Args:
        iface: the Interface dict from graph.json (must have bolt_pattern)
        parent_world_origin_mm: legacy override (used only in fallback path)
        child_world_origin_mm: legacy override (used only in fallback path)
        parts_by_id: optional {part_id: part_dict} index from graph.json
        world_transforms: optional {part_id: ([tx,ty,tz], [qx,qy,qz,qw])} from
            ``solve_world_transforms(graph)``
        bolts_per_part_spec: optional bolts_per_part.yaml inner dict
            from ``_load_bolts_per_part_spec(graph_path)``

    Returns:
        list of Bolt entries with world coords, axis, spec, stiffness.
    """
    bp = iface.get("bolt_pattern")
    if not bp:
        return []

    pcd_mm = float(bp["pcd_mm"])
    n_bolts = int(bp["n_bolts"])
    bolt_spec = str(bp["bolt_spec"])
    access_dir = bp.get("access_direction", [0.0, 0.0, 1.0])
    k_axial, k_shear = stiffness_for_bolt_spec(bolt_spec)

    pid = iface.get("parent_part", "")
    cid = iface.get("child_part", "")
    iid = iface.get("id", "?")

    # ---- Preferred path: per-hole datums + kinematic chain ----------------
    if parts_by_id is not None and world_transforms is not None:
        bolts = _generate_bolts_from_datums(
            iface=iface,
            bp=bp,
            parts_by_id=parts_by_id,
            world_transforms=world_transforms,
            bolts_per_part_spec=bolts_per_part_spec or {},
            n_bolts=n_bolts,
            pcd_mm=pcd_mm,
            bolt_spec=bolt_spec,
            k_axial=k_axial,
            k_shear=k_shear,
        )
        if bolts:
            return bolts
        # If we couldn't find datums, fall through to legacy path.

    # ---- Legacy fallback: PCD + access_direction --------------------------
    # Default to bolt_pattern's own world_z_mm if no override given.
    if parent_world_origin_mm is None:
        wz = bp.get("world_z_mm")
        parent_world_origin_mm = (0.0, 0.0, float(wz)) if wz is not None else (0.0, 0.0, 0.0)
    if child_world_origin_mm is None:
        child_world_origin_mm = parent_world_origin_mm

    r_mm = pcd_mm / 2.0
    bolts: List[Bolt] = []
    ad = (float(access_dir[0]), float(access_dir[1]), float(access_dir[2]))
    norm = math.sqrt(ad[0] ** 2 + ad[1] ** 2 + ad[2] ** 2) or 1.0
    axis = (ad[0] / norm, ad[1] / norm, ad[2] / norm)

    for i in range(n_bolts):
        theta = 2 * math.pi * i / n_bolts
        x = r_mm * math.cos(theta)
        y = r_mm * math.sin(theta)
        z = 0.0
        parent_world = (
            parent_world_origin_mm[0] + x,
            parent_world_origin_mm[1] + y,
            parent_world_origin_mm[2] + z,
        )
        child_world = (
            child_world_origin_mm[0] + x,
            child_world_origin_mm[1] + y,
            child_world_origin_mm[2] + z,
        )
        clock_deg = math.degrees(math.atan2(y, x)) % 360.0
        bolts.append(Bolt(
            interface_id=iid,
            bolt_index=i,
            parent_part=pid,
            child_part=cid,
            parent_world_mm=parent_world,
            child_world_mm=child_world,
            axis_unit=axis,
            bolt_spec=bolt_spec,
            pcd_mm=pcd_mm,
            clock_deg=clock_deg,
            k_axial_N_per_m=k_axial,
            k_shear_N_per_m=k_shear,
        ))
    return bolts


def _generate_bolts_from_datums(
    *,
    iface: Dict[str, Any],
    bp: Dict[str, Any],
    parts_by_id: Dict[str, Dict[str, Any]],
    world_transforms: Dict[str, Tuple[List[float], List[float]]],
    bolts_per_part_spec: Dict[str, Any],
    n_bolts: int,
    pcd_mm: float,
    bolt_spec: str,
    k_axial: float,
    k_shear: float,
) -> List[Bolt]:
    """Per-hole-datum bolt generator. Returns [] if datums are missing.

    Owner-side hole datums carry the per-bolt axis directly via their
    ``direction`` field — for axial joints it's ±Z, for radial joints
    it's the inward radial unit vector at that clock angle. The world
    transform of the owner part places the hole in vehicle frame.

    Endpoints:
      parent_world_mm = head anchor (owner-side hole entry, optionally
                        offset outward by head_radial_offset_mm or
                        offset axially by head_axial_offset_mm per YAML)
      child_world_mm  = shank tip   = parent_world + shank_length * axis
    """
    pp = iface.get("parent_part", "")
    cp = iface.get("child_part", "")
    iid = iface.get("id", "?")

    # Decide owner part and which hole-datum prefix to read.
    # The interface's `origin` Z indicates which side of the parent the
    # joint sits on: if origin.z > 0 (e.g. boattail.top at z=330) we use
    # ``top_bolt_hole_N``, else ``bot_bolt_hole_N``. Parent is the owner
    # by default; if a YAML group exists and names a different owner we
    # honor that.
    yaml_match = _find_yaml_group_for_iface(bolts_per_part_spec, iface)
    owner_id: str
    hole_prefix: str
    yaml_group: Dict[str, Any]
    orientation_kind: str
    if yaml_match is not None:
        owner_id, _gname, yaml_group = yaml_match
        # bolt_circle_ref tells us 'top' or 'bot'
        bcr = str(yaml_group.get("bolt_circle_ref", ""))
        hole_prefix = "top_bolt_hole" if "top" in bcr else "bot_bolt_hole"
        orientation_kind = str(yaml_group.get("orientation", "")).lower()
    else:
        # No YAML match: default to parent-as-owner; infer top/bot from iface.
        owner_id = pp
        # Inspect parent's top vs bot circle z to pick the side closer to the interface origin.
        owner_part0 = parts_by_id.get(owner_id, {})
        iface_origin = iface.get("origin") or [0.0, 0.0, 0.0]
        try:
            iface_oz = float(iface_origin[2])
        except (TypeError, ValueError, IndexError):
            iface_oz = 0.0
        top_circle = _find_datum(owner_part0, "top_bolt_circle")
        bot_circle = _find_datum(owner_part0, "bot_bolt_circle")
        top_z = _datum_origin(top_circle)[2] if top_circle else None
        bot_z = _datum_origin(bot_circle)[2] if bot_circle else None
        if top_z is not None and (bot_z is None or abs(iface_oz - top_z) <= abs(iface_oz - bot_z)):
            hole_prefix = "top_bolt_hole"
        else:
            hole_prefix = "bot_bolt_hole"
        yaml_group = {}
        orientation_kind = ""

    owner_part = parts_by_id.get(owner_id)
    if owner_part is None:
        return []

    owner_transform = world_transforms.get(owner_id)
    if owner_transform is None:
        # Owner not in kinematic chain; can't place world coords.
        return []
    owner_t, owner_q = owner_transform

    # Collect available per-hole datums.
    holes: List[Tuple[int, List[float], List[float]]] = []  # (idx, origin_local, dir_local)
    for i in range(max(n_bolts, 64)):
        d = _find_datum(owner_part, f"{hole_prefix}_{i}")
        if d is None:
            if i < n_bolts:
                # Missing a hole we expected; bail out and let caller fall back.
                # But only if we have NO holes at all — partial sets are still useful.
                pass
            continue
        origin_local = _datum_origin(d)
        dir_local = _datum_direction(d)
        if dir_local is None:
            # Hole datum lacks direction; skip.
            continue
        holes.append((i, origin_local, dir_local))
    # If we didn't find any holes, fall back to the legacy path.
    if not holes:
        return []
    # Sort by index, take first n_bolts present.
    holes.sort(key=lambda h: h[0])

    # Decide axial vs radial if YAML didn't tell us. Inspect first hole.
    if not orientation_kind:
        first_dir = holes[0][2]
        orientation_kind = "axial" if abs(first_dir[2]) > 0.8 else "radial"

    # Pull offsets from YAML (or sensible defaults).
    head_axial_offset = float(yaml_group.get("head_axial_offset_mm", 0.0))
    head_radial_offset = float(yaml_group.get("head_radial_offset_mm", 0.0))
    shank_length_mm = float(yaml_group.get("shank_length_mm", 0.0))
    if shank_length_mm <= 0.01:
        # Derive from bolt spec geometry table.
        geo = _BOLT_SPEC_GEOMETRY.get(bolt_spec)
        shank_length_mm = float(geo[1]) if geo else 25.4

    bolts: List[Bolt] = []
    for idx, origin_local, dir_local in holes:
        # The owner-part hole direction is the bolt INSERTION axis in
        # the part-local frame (it points the way the shank goes).
        # For radial bolts this points INWARD toward the part's Z axis;
        # for axial bolts it's ±Z.
        insert_local = dir_local  # already unit vector

        # Head anchor offset (LOCAL frame):
        #   * radial: head sits OUTBOARD of the hole origin, i.e. opposite
        #     the inward insertion direction. Step ‘head_radial_offset’
        #     along -insert_local.
        #   * axial:  head sits offset along +insert_local by
        #     head_axial_offset (negative offset puts head on the opposite
        #     side). Convention matches the dashboard's bolt loader so
        #     the FEA and bolt overlays line up.
        if orientation_kind == "radial":
            head_anchor_local = [
                origin_local[0] - head_radial_offset * insert_local[0],
                origin_local[1] - head_radial_offset * insert_local[1],
                origin_local[2] - head_radial_offset * insert_local[2],
            ]
        else:
            head_anchor_local = [
                origin_local[0] + head_axial_offset * insert_local[0],
                origin_local[1] + head_axial_offset * insert_local[1],
                origin_local[2] + head_axial_offset * insert_local[2],
            ]

        # Shank tip in LOCAL frame = head anchor + shank_length * insert_local.
        shank_tip_local = [
            head_anchor_local[0] + shank_length_mm * insert_local[0],
            head_anchor_local[1] + shank_length_mm * insert_local[1],
            head_anchor_local[2] + shank_length_mm * insert_local[2],
        ]

        # Transform LOCAL → WORLD using owner part's solved (t, q).
        ha_rot = _qrot(owner_q, head_anchor_local)
        head_anchor_world = (
            owner_t[0] + ha_rot[0],
            owner_t[1] + ha_rot[1],
            owner_t[2] + ha_rot[2],
        )
        st_rot = _qrot(owner_q, shank_tip_local)
        shank_tip_world = (
            owner_t[0] + st_rot[0],
            owner_t[1] + st_rot[1],
            owner_t[2] + st_rot[2],
        )
        # Axis is the insertion direction rotated into world.
        ax_rot = _qrot(owner_q, insert_local)
        an = math.sqrt(ax_rot[0] ** 2 + ax_rot[1] ** 2 + ax_rot[2] ** 2) or 1.0
        axis_world = (ax_rot[0] / an, ax_rot[1] / an, ax_rot[2] / an)

        # Clock from hole_origin_local (atan2 of XY components).
        clock_deg = math.degrees(math.atan2(origin_local[1], origin_local[0])) % 360.0

        bolts.append(Bolt(
            interface_id=iid,
            bolt_index=idx,
            parent_part=pp,
            child_part=cp,
            parent_world_mm=head_anchor_world,
            child_world_mm=shank_tip_world,
            axis_unit=axis_world,
            bolt_spec=bolt_spec,
            pcd_mm=pcd_mm,
            clock_deg=clock_deg,
            k_axial_N_per_m=k_axial,
            k_shear_N_per_m=k_shear,
        ))
        if len(bolts) >= n_bolts:
            break

    return bolts


def resolve_load_attach(load_case: Dict[str, Any]) -> LoadAttachResolved:
    """Convert a LoadCase dict to a resolved attach descriptor.

    Does NOT compute world coordinates of bolt-targeted loads (those need
    interface origins from graph.json). Just packages the LoadCase fields
    into a struct the FEA solver can consume directly.
    """
    attach = load_case.get("attach", {})
    pos = attach.get("position", [0.0, 0.0, 0.0])
    direction = load_case.get("direction", [1.0, 0.0, 0.0])

    return LoadAttachResolved(
        load_case_id=load_case["id"],
        part=attach.get("part", "_assembly"),
        interface_id=attach.get("interface"),
        position_mm=(float(pos[0]), float(pos[1]), float(pos[2])),
        direction_unit=(float(direction[0]), float(direction[1]), float(direction[2])),
        magnitude_n=float(load_case.get("magnitude_n") or 0.0),
        coordinate_frame=load_case.get("coordinate_frame", "vehicle"),
        bolt_index=attach.get("bolt_index"),
        clock_deg=attach.get("clock_deg"),
        raw=dict(load_case),
    )


def prepare(
    graph_path: str | Path,
    load_case_id: Optional[str] = None,
    *,
    interface_origins_mm: Optional[Dict[str, Tuple[float, float, float]]] = None,
) -> FeaSetup:
    """One-shot FEA setup from a mechatron graph.json.

    Reads the graph, resolves the load case (if given), enumerates all
    structural bolts with stiffness, and returns a FeaSetup struct.

    Args:
        graph_path: path to graph.json (mechatron canonical)
        load_case_id: e.g. "LC-004a" — if None, no load attach is resolved
        interface_origins_mm: optional override map {iface_id: (x,y,z) world}.
            If not given, bolts are returned in interface-local coords with
            origin (0,0,0); the FEA setup is responsible for combining with
            the assembly transforms from the graph's interface solve.

    Returns:
        FeaSetup with bolts, load_attach (if requested), and raw refs.
    """
    g = load_graph(graph_path)
    setup = FeaSetup(graph_path=str(graph_path))

    # 2026-05-27: enable the preferred per-hole-datum path.
    # Build the parts index + kinematic-chain world transforms + load
    # bolts_per_part.yaml so generate_bolts_for_interface() can produce
    # correctly oriented per-bolt axes (radial joints no longer come out
    # as ±Z).
    parts_by_id = _index_parts(g)
    world_transforms = solve_world_transforms(g)
    bolts_per_part_spec = _load_bolts_per_part_spec(graph_path)

    # Enumerate bolts. Origin resolution (legacy fallback path only):
    # 1. Explicit override in interface_origins_mm wins (caller may override)
    # 2. Else read `bolt_pattern.world_z_mm` from the graph (preferred since
    #    2026-05-20 reconciliation — graph.json now carries mate-solve Z)
    # 3. Else (0,0,0) and emit warning
    interface_origins_mm = interface_origins_mm or {}
    for iface in list_interfaces_with_bolt_patterns(g):
        if iface["id"] in interface_origins_mm:
            origin = interface_origins_mm[iface["id"]]
        else:
            origin = None  # generate_bolts_for_interface() will read world_z_mm from bp
        setup.bolts.extend(generate_bolts_for_interface(
            iface,
            parent_world_origin_mm=origin,
            parts_by_id=parts_by_id,
            world_transforms=world_transforms,
            bolts_per_part_spec=bolts_per_part_spec,
        ))
    setup.raw_interfaces = list_interfaces_with_bolt_patterns(g)

    # Resolve load case if requested
    if load_case_id is not None:
        lc = get_load_case(g, load_case_id)
        setup.load_case_id = load_case_id
        setup.load_attach = resolve_load_attach(lc)
        setup.raw_load_case = lc

    return setup
