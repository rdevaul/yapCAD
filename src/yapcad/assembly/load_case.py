"""LoadCase: structural FEA load case attached to an Assembly.

A LoadCase carries the information an FEA setup needs to know **what force
is applied, where, in which direction, and against which interface** —
without the FEA author having to invent or hand-extract that data from
side-files.

Designed to round-trip through Mechatron's `LoadCase` graph entity
(see `crates/assembly-graph/src/types.rs`). The yapCAD dataclass lives
on the `Assembly` instance under `assembly.load_cases`.

History:
    Created 2026-05-20 in response to the Agentic-1 chute-deploy FEA
    iterations (v4→v7) where load case data was scattered across:
      - designs/assembly-dashboard-v2/index.html LOAD_CASES (hardcoded)
      - requirements/LOAD-CASES-DRAFT.md (draft notes)
      - specs/load-cases/canonical.json (jarvis-rich's hand extract)
    The lack of a single source of truth caused several iteration
    cycles to use wrong attach points / clock angles / part assignments.
    See Whiteboard message 20260520201704-46e1d134 for context.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Any, Dict, List, Optional, Sequence


# Canonical coordinate frame enum.
COORDINATE_FRAME_VEHICLE = "vehicle"  # X = axial-fwd (dashboard convention)
COORDINATE_FRAME_MESH = "mesh"        # Z = axial-fwd (FEA mesh convention)
VALID_FRAMES = (COORDINATE_FRAME_VEHICLE, COORDINATE_FRAME_MESH)

# Canonical status values.
STATUS_CONFIRMED = "confirmed"
STATUS_ESTIMATED = "estimated"
STATUS_TRIBAL = "tribal"
STATUS_PLACEHOLDER = "placeholder"
VALID_STATUSES = (STATUS_CONFIRMED, STATUS_ESTIMATED, STATUS_TRIBAL, STATUS_PLACEHOLDER)

# Canonical load groups (mirrors dashboard LOAD_CASES.group).
GROUP_INERTIAL = "inertial"
GROUP_RECOVERY = "recovery"
GROUP_PRESSURE = "pressure"
GROUP_HANDLING = "handling"
GROUP_CUSTOM = "custom"
VALID_GROUPS = (GROUP_INERTIAL, GROUP_RECOVERY, GROUP_PRESSURE, GROUP_HANDLING, GROUP_CUSTOM)


@dataclass
class LoadAttach:
    """Where a LoadCase applies its force/moment on the assembly.

    Resolution priority for FEA setup:
        1. If `interface` is set, the load enters via that mechatron Interface
           (typically a bolted joint). The FEA solver looks up the interface's
           bolt_pattern and applies the load at the bolt nearest `clock_deg`
           (or at the specific `bolt_index` if given).
        2. Otherwise, if `part` is set and != "_assembly", the load is applied
           at `position` (interpreted in the load case's coordinate_frame) on
           that part's surface DOFs. Useful for distributed nose-tip loads,
           engine thrust at engine-mount, etc.
        3. If `part == "_assembly"`, the load is a whole-stack body force /
           inertial load. `position` is a notional centroid for visualization
           only; FEA applies the equivalent ρ·a body force on all parts.
    """
    part: str  # part_id, or "_assembly" for whole-stack inertial
    position: Sequence[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])
    interface: Optional[str] = None  # mechatron interface id (preferred)
    clock_deg: Optional[float] = None  # 0-360°, vehicle frame
    bolt_index: Optional[int] = None  # 0..n_bolts-1 — overrides clock_deg if set

    def to_dict(self) -> Dict[str, Any]:
        return {
            "part": self.part,
            "position": list(self.position),
            "interface": self.interface,
            "clock_deg": self.clock_deg,
            "bolt_index": self.bolt_index,
        }


@dataclass
class LoadCase:
    """A structural load case applied to an Assembly.

    Round-trips through Mechatron's `LoadCase` type.
    """
    id: str
    name: str
    group: str  # one of VALID_GROUPS
    direction: Sequence[float]  # unit vector in coordinate_frame
    attach: LoadAttach
    magnitude_n: Optional[float] = None  # Force magnitude in Newtons
    magnitude_g: Optional[float] = None  # Equivalent g-loading (if applicable)
    coordinate_frame: str = COORDINATE_FRAME_VEHICLE
    description: str = ""
    source: str = ""  # where this load case came from
    status: str = STATUS_ESTIMATED
    warn: bool = False  # surface a warning in the dashboard

    def __post_init__(self):
        if self.group not in VALID_GROUPS:
            raise ValueError(
                f"LoadCase.group must be one of {VALID_GROUPS}, got {self.group!r}"
            )
        if self.coordinate_frame not in VALID_FRAMES:
            raise ValueError(
                f"LoadCase.coordinate_frame must be one of {VALID_FRAMES}, "
                f"got {self.coordinate_frame!r}"
            )
        if self.status not in VALID_STATUSES:
            raise ValueError(
                f"LoadCase.status must be one of {VALID_STATUSES}, got {self.status!r}"
            )
        if self.magnitude_n is None and self.magnitude_g is None:
            raise ValueError(
                f"LoadCase {self.id!r}: at least one of magnitude_n / magnitude_g must be set"
            )
        if len(self.direction) != 3:
            raise ValueError(
                f"LoadCase {self.id!r}: direction must be 3-vector, got len {len(self.direction)}"
            )

    def to_dict(self) -> Dict[str, Any]:
        d: Dict[str, Any] = {
            "id": self.id,
            "name": self.name,
            "group": self.group,
            "direction": list(self.direction),
            "attach": self.attach.to_dict(),
            "coordinate_frame": self.coordinate_frame,
            "description": self.description,
            "source": self.source,
            "status": self.status,
        }
        if self.magnitude_n is not None:
            d["magnitude_n"] = self.magnitude_n
        if self.magnitude_g is not None:
            d["magnitude_g"] = self.magnitude_g
        if self.warn:
            d["warn"] = True
        return d

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "LoadCase":
        """Build from a dict (e.g. graph.json round-trip)."""
        attach_data = data.get("attach", {})
        attach = LoadAttach(
            part=attach_data["part"],
            position=attach_data.get("position", [0.0, 0.0, 0.0]),
            interface=attach_data.get("interface"),
            clock_deg=attach_data.get("clock_deg"),
            bolt_index=attach_data.get("bolt_index"),
        )
        return cls(
            id=data["id"],
            name=data["name"],
            group=data["group"],
            direction=data["direction"],
            attach=attach,
            magnitude_n=data.get("magnitude_n"),
            magnitude_g=data.get("magnitude_g"),
            coordinate_frame=data.get("coordinate_frame", COORDINATE_FRAME_VEHICLE),
            description=data.get("description", ""),
            source=data.get("source", ""),
            status=data.get("status", STATUS_ESTIMATED),
            warn=data.get("warn", False),
        )


@dataclass
class BoltPattern:
    """A bolt circle attached to an Assembly interface.

    Round-trips through Mechatron's `BoltPattern` type (see
    `crates/assembly-graph/src/types.rs`).

    Fields match the mechatron schema; `bolt_spec` is a canonical string
    like "1-4-20-x-1in-button-head-shcs" that can be resolved against the
    cots-db MCP for stiffness, weight, COTS line item, etc.
    """
    pcd_mm: float  # Pitch circle DIAMETER in mm (not radius!)
    n_bolts: int
    bolt_spec: str
    access_direction: Sequence[float]  # 3-vec; tool insertion direction

    def __post_init__(self):
        if self.pcd_mm <= 0:
            raise ValueError(f"BoltPattern.pcd_mm must be positive, got {self.pcd_mm}")
        if self.n_bolts <= 0:
            raise ValueError(f"BoltPattern.n_bolts must be positive, got {self.n_bolts}")
        if len(self.access_direction) != 3:
            raise ValueError(
                f"BoltPattern.access_direction must be 3-vector, got len {len(self.access_direction)}"
            )

    def to_dict(self) -> Dict[str, Any]:
        return {
            "pcd_mm": float(self.pcd_mm),
            "n_bolts": int(self.n_bolts),
            "bolt_spec": str(self.bolt_spec),
            "access_direction": list(self.access_direction),
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "BoltPattern":
        return cls(
            pcd_mm=data["pcd_mm"],
            n_bolts=data["n_bolts"],
            bolt_spec=data["bolt_spec"],
            access_direction=data.get("access_direction", [0.0, 0.0, 1.0]),
        )
