"""The motion model a range-of-motion sweep is run against.

A :class:`RomModel` is a read-only view of a DSL assembly: the parts, the
field (SDF) each part was authored from, the mates that join them, and the
poses the assembly solver produces for a set of joint values. It adds no
kinematics of its own -- the assembly's own mate graph remains the source of
pose, so a sweep cannot disagree with the design it came from.

Everything a sweep needs is already in the DSL::

    @meta(assembly.datums=[
        {"id": "pivot", "kind": "axis",
         "origin_mm": [0.0, 0.0, 10.0], "direction": [0.0, 0.0, 1.0]}
    ])
    command BASE() -> solid:
        emit box(40.0, 40.0, 20.0)

    command BUILD() -> solid:
        let a: assembly = assembly("demo")
        add_part(a, BASE(), "base")
        add_part(a, ARM(), "arm")
        add_named_mate(a, "elbow", "revolute", "base", "pivot", "arm", "hub")
        set_mate_limits(a, "elbow", radians(-90.0), radians(90.0))
        solve_assembly(a, "base")
        emit assembly_compound(a)

The datums place the joint, ``add_named_mate`` names it, and
``set_mate_limits`` bounds it. Build the design as fields and the model
follows from the result::

    result = compile_and_run(source, "BUILD", {}, representation="sdf",
                             source_path=path)
    model = RomModel.from_result(result)

Joint values are **radians** throughout, the same unit as
``set_mate_limits`` and ``set_joint_position``;
:func:`~yapcad.rom.sweep.schedule_from_model` takes a step in degrees for
convenience, and the reports print degrees.
"""

from __future__ import annotations

import math
from collections import deque
from dataclasses import dataclass
from typing import (
    Dict, Iterable, List, Mapping, Optional, Sequence, Tuple,
)

import numpy as np

__all__ = [
    "RomError",
    "RomJoint",
    "RomLink",
    "RomModel",
    "assembly_from_result",
]


class RomError(Exception):
    """A range-of-motion model or sweep that cannot be honoured."""


@dataclass(frozen=True)
class RomJoint:
    """One mate of the assembly, as a degree of freedom.

    ``parent`` and ``child`` are the mate's own ``part_a`` and ``part_b``;
    the solver roots the tree wherever it was told to, so they describe the
    mate, not necessarily the direction pose propagates in.

    ``lower`` and ``upper`` are radians, as ``set_mate_limits`` records
    them, and are None when the design did not bound the joint.
    """

    name: str
    parent: str
    child: str
    kind: str
    lower: Optional[float] = None
    upper: Optional[float] = None

    @property
    def movable(self) -> bool:
        return self.kind == "revolute"

    @property
    def bounded(self) -> bool:
        return self.lower is not None and self.upper is not None

    def describe(self) -> str:
        if not self.movable:
            return f"{self.name}: {self.kind} ({self.parent} - {self.child})"
        if self.bounded:
            span = (f"{math.degrees(self.lower):.1f} deg to "
                    f"{math.degrees(self.upper):.1f} deg")
        else:
            span = "unbounded"
        return (f"{self.name}: revolute, {span} "
                f"({self.parent} - {self.child})")


@dataclass(frozen=True)
class RomLink:
    """One part instance, with the field it was authored from.

    ``field`` is the SDF node recovered from the part's construction
    record, or None -- in which case ``reason`` says why, and any sweep
    that would have to measure this part says so rather than quietly
    leaving it out.
    """

    name: str
    field: object = None
    reason: Optional[str] = None

    @property
    def has_field(self) -> bool:
        return self.field is not None


def assembly_from_result(result) -> object:
    """The ``Assembly`` a DSL run retained, or raise.

    A command ending in ``emit assembly_compound(a)`` keeps its assembly on
    the emitted value, which is the same object the package builder reads.
    """
    if result is None:
        raise RomError("No execution result was given.")
    if not getattr(result, "success", False):
        raise RomError(
            "The DSL run failed, so there is no assembly to sweep: "
            f"{getattr(result, 'error_message', None)}"
        )
    emit = getattr(result, "emit_result", None)
    value = getattr(emit, "value", None)
    annotations = getattr(value, "annotations", None) or {}
    assembly = annotations.get("assembly") if isinstance(
        annotations, dict) else None
    if assembly is None:
        raise RomError(
            "This result carries no assembly. A range-of-motion sweep "
            "needs the mate graph, so the command has to build one "
            "(assembly(), add_part(), add_named_mate()) and end with "
            "'emit assembly_compound(...)'."
        )
    return assembly


def _field_of(geometry, part: str) -> Tuple[object, Optional[str]]:
    """The SDF node behind ``geometry``, or (None, reason).

    No fallback: a bounding box in place of the part would answer every
    clearance question with the wrong number.
    """
    if geometry is None:
        return None, "the part has no geometry in the assembly"
    try:
        from yapcad.construction import get_construction
        from yapcad.sdf import from_construction, is_sdf_solid
    except ImportError as exc:  # pragma: no cover - sdf ships with yapCAD
        return None, f"yapcad.sdf is unavailable ({exc})"
    if not is_sdf_solid(geometry):
        return None, (
            "the part was not authored as a field; build the design with "
            "representation='sdf' so that every part carries one"
        )
    node = from_construction(get_construction(geometry))
    if node is None:
        return None, "the part's construction record held no field"
    return node, None


class RomModel:
    """Parts, fields, joints and poses of one assembly.

    Build it from a DSL run with :meth:`from_result`, or from an
    ``Assembly`` directly with :meth:`from_assembly`.
    """

    def __init__(
        self,
        assembly,
        root: str,
        links: Mapping[str, RomLink],
        joints: Sequence[RomJoint],
    ) -> None:
        self.assembly = assembly
        self.root = root
        self.links: Dict[str, RomLink] = dict(links)
        self.joints: Dict[str, RomJoint] = {j.name: j for j in joints}

    # -- construction ----------------------------------------------------

    @classmethod
    def from_assembly(cls, assembly, root: Optional[str] = None) -> "RomModel":
        """Read parts, fields and mates from a solved ``Assembly``."""
        parts = list(getattr(assembly, "parts", {}) or {})
        if not parts:
            raise RomError("The assembly has no parts.")
        chosen = root or getattr(assembly, "_root_part", None)
        if chosen is None:
            raise RomError(
                "The assembly has not been solved and no root part was "
                "given. Call solve_assembly(...) in the DSL (or pass "
                "root=) so that poses have somewhere to start."
            )
        if chosen not in parts:
            raise RomError(
                f"Root part '{chosen}' is not in the assembly. "
                f"Parts: {', '.join(sorted(parts))}."
            )

        geometry = getattr(assembly, "geometry", {}) or {}
        links = {}
        for name in parts:
            field, reason = _field_of(geometry.get(name), name)
            links[name] = RomLink(name=name, field=field, reason=reason)

        joints = []
        for mate in getattr(assembly, "mates", []) or []:
            limits = getattr(mate, "limits", None)
            kind = getattr(mate.mate_type, "value", str(mate.mate_type))
            joints.append(RomJoint(
                name=mate.name,
                parent=mate.part_a,
                child=mate.part_b,
                kind=kind,
                lower=getattr(limits, "min_value", None),
                upper=getattr(limits, "max_value", None),
            ))
        return cls(assembly, chosen, links, joints)

    @classmethod
    def from_result(cls, result, root: Optional[str] = None) -> "RomModel":
        """Read the model from a DSL ``ExecutionResult``."""
        return cls.from_assembly(assembly_from_result(result), root=root)

    # -- queries ---------------------------------------------------------

    @property
    def parts(self) -> List[str]:
        return list(self.links)

    def movable_joints(self) -> List[RomJoint]:
        """The joints a sweep can drive, in declaration order."""
        return [j for j in self.joints.values() if j.movable]

    def joint(self, name: str) -> RomJoint:
        try:
            return self.joints[name]
        except KeyError:
            known = ", ".join(sorted(self.joints)) or "none"
            raise RomError(
                f"No joint named '{name}' in this assembly. Joints: {known}."
            ) from None

    def field(self, part: str):
        """The field for ``part``, or raise saying why there is none."""
        link = self.links.get(part)
        if link is None:
            known = ", ".join(sorted(self.links))
            raise RomError(
                f"No part named '{part}' in this assembly. Parts: {known}."
            )
        if link.field is None:
            raise RomError(f"Part '{part}' has no field: {link.reason}.")
        return link.field

    def parts_without_fields(self) -> Dict[str, str]:
        """``{part: reason}`` for every part a sweep could not measure."""
        return {name: link.reason or "unknown"
                for name, link in self.links.items() if not link.has_field}

    def mate_pairs(self) -> List[Tuple[str, str]]:
        """Part pairs that share a mate -- the ones a joint moves apart."""
        pairs = []
        for joint in self.joints.values():
            pair = (joint.parent, joint.child)
            if pair[0] in self.links and pair[1] in self.links:
                if pair not in pairs and pair[::-1] not in pairs:
                    pairs.append(pair)
        return pairs

    def all_pairs(self) -> List[Tuple[str, str]]:
        """Every unordered pair of parts."""
        names = self.parts
        return [(names[i], names[j])
                for i in range(len(names))
                for j in range(i + 1, len(names))]

    def mate_chains(self) -> Dict[str, Optional[Tuple[str, ...]]]:
        """``{part: mates from the root to it}``, or None if unreachable.

        The walk mirrors the assembly solver's: the mate graph is
        undirected, and pose propagates outwards from the root.
        """
        adjacency: Dict[str, List[Tuple[str, str]]] = {
            name: [] for name in self.links
        }
        for joint in self.joints.values():
            if joint.parent in adjacency and joint.child in adjacency:
                adjacency[joint.parent].append((joint.child, joint.name))
                adjacency[joint.child].append((joint.parent, joint.name))

        chains: Dict[str, Optional[Tuple[str, ...]]] = {
            name: None for name in self.links
        }
        chains[self.root] = ()
        queue = deque([self.root])
        while queue:
            current = queue.popleft()
            for neighbour, mate in adjacency[current]:
                if chains[neighbour] is None and neighbour != self.root:
                    chains[neighbour] = tuple(chains[current]) + (mate,)
                    queue.append(neighbour)
        return chains

    def moving_pairs(
        self, joints: Optional[Iterable[str]] = None
    ) -> List[Tuple[str, str]]:
        """Pairs whose relative pose changes when ``joints`` move.

        Two parts rigidly fixed to each other hold still however the
        assembly is posed, so measuring them at every configuration
        reports the same design fact over and over -- while the pair that
        actually closes may share no mate at all (an arm and the stop it
        swings towards). This is that set: every pair separated by at
        least one of the joints being driven.
        """
        driven = set(joints) if joints is not None else {
            joint.name for joint in self.movable_joints()
        }
        for name in driven:
            self.joint(name)
        chains = self.mate_chains()
        pairs = []
        for part_a, part_b in self.all_pairs():
            chain_a = chains.get(part_a)
            chain_b = chains.get(part_b)
            if chain_a is None or chain_b is None:
                # Unreachable from the root: measure rather than assume.
                pairs.append((part_a, part_b))
                continue
            between = set(chain_a) ^ set(chain_b)
            if between & driven:
                pairs.append((part_a, part_b))
        return pairs

    # -- pose ------------------------------------------------------------

    def check_joint_values(self, values: Mapping[str, float]) -> None:
        """Raise unless every name is a joint and every value is in range.

        A sweep that quietly dropped an unknown joint name, or clamped a
        value into range, would report clearances for poses the caller
        never asked about.
        """
        for name, value in values.items():
            joint = self.joint(name)
            if not joint.movable:
                raise RomError(
                    f"Joint '{name}' is {joint.kind}, not revolute, so it "
                    "cannot be driven."
                )
            if not math.isfinite(float(value)):
                raise RomError(f"Joint '{name}' was given {value!r}.")
            if joint.lower is not None and value < joint.lower - 1e-12:
                raise RomError(
                    f"Joint '{name}' was given "
                    f"{math.degrees(value):.3f} deg, below its limit of "
                    f"{math.degrees(joint.lower):.3f} deg."
                )
            if joint.upper is not None and value > joint.upper + 1e-12:
                raise RomError(
                    f"Joint '{name}' was given "
                    f"{math.degrees(value):.3f} deg, above its limit of "
                    f"{math.degrees(joint.upper):.3f} deg."
                )

    def pose(self, values: Optional[Mapping[str, float]] = None) -> Dict:
        """World transforms for every part at ``values`` (radians).

        Delegates to the assembly solver, so the poses are the design's
        own rather than a second opinion. Solving re-poses the assembly
        in place, as ``set_joint_position`` does;
        :meth:`joint_values` reads the pose it is left in, and
        :func:`~yapcad.rom.sweep.sweep_rom` restores the one it found.
        """
        joint_values = {str(k): float(v) for k, v in (values or {}).items()}
        self.check_joint_values(joint_values)
        result = self.assembly.solve(self.root, joint_values)
        if not getattr(result, "success", False):
            errors = "; ".join(getattr(result, "errors", []) or ["unknown"])
            raise RomError(f"The assembly could not be posed: {errors}")
        return {name: np.asarray(matrix, dtype=float)
                for name, matrix in result.transforms.items()}

    def joint_values(self) -> Dict[str, float]:
        """The joint values the assembly is currently posed at (radians)."""
        prescribed = getattr(self.assembly, "_prescribed_joint_values", None)
        return {str(k): float(v) for k, v in (prescribed or {}).items()}

    # -- description -----------------------------------------------------

    def describe(self) -> str:
        """What a sweep of this model can and cannot measure."""
        lines = [
            f"Assembly '{getattr(self.assembly, 'name', '?')}' "
            f"rooted at '{self.root}'",
            f"  parts: {len(self.links)}",
        ]
        missing = self.parts_without_fields()
        if missing:
            lines.append(f"  parts without a field: {len(missing)}")
            for name, reason in sorted(missing.items()):
                lines.append(f"    {name}: {reason}")
        movable = self.movable_joints()
        lines.append(f"  joints: {len(self.joints)} "
                     f"({len(movable)} movable)")
        for joint in self.joints.values():
            lines.append(f"    {joint.describe()}")
        return "\n".join(lines)
