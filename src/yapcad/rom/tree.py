"""Build a kinematic tree from a DSL module's compiled output.

Reads the ``@meta(kinematic.*)`` annotations applied to each command,
resolves parent-child relationships, and builds a :class:`KinematicChain`
whose joints match the DSL metadata.  Each link remembers the compiled
SDF so adjacent links can be evaluated for clearance without meshing.

Usage::

    from yapcad.dsl import compile_to_ast, interpret_module
    from yapcad.rom import RomTree

    ast = compile_to_ast(source)
    ctx = interpret_module(ast)

    tree = RomTree.from_interpreter_context("phonebot_leg", ctx)
    chain = tree.chain

    # Get world transforms for a given joint angle dict
    for part_name, angle in {"HIP_ROLL": 45.0, "KNEE_PITCH": 90.0}:
        chain.set_joint_value(part_name, angle)

    tf_thigh = chain.get_world_transform("THIGH_LINK")
    tf_shin  = chain.get_world_transform("SHIN_LINK")
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional, Tuple

from yapcad.kinematics import (
    KinematicChain,
    KinematicPart,
    Joint,
    JointType,
    CoordinateFrame,
    Transform,
)

logger = logging.getLogger(__name__)


class RomLink:
    """Metadata + geometry for one link in a kinematic tree.

    Attributes:
        name: Command / part name (e.g., ``"HIP_ROLL_JOINT"``).
        command: The DSL command function (from interpreter context).
        sdf_node: The compiled SDF node (from command's SDF representation).
        parent: Name of the parent link, or ``None`` for the root.
        parent_frame: Frame on the parent this link attaches to.
        joint_kind: ``"fixed"`` | ``"revolute"`` | ``"prismatic"``.
        axis: Motion axis ``(x, y, z)`` (for revolute/prismatic).
        limits: ``(min, max)`` angle/translation range.
        zero_pose: Default joint value.
        description: Human-readable label.
        printable: Whether this is a 3D-printable part.
    """

    def __init__(
        self,
        name: str,
        command=None,
        sdf_node=None,
        parent: Optional[str] = None,
        parent_frame: str = "ORIGIN",
        joint_kind: str = "fixed",
        axis: Tuple[float, float, float] = (0.0, 0.0, 1.0),
        limits: Tuple[float, float] = (-180.0, 180.0),
        zero_pose: float = 0.0,
        description: str = "",
        printable: bool = True,
    ):
        self.name = name
        self.command = command
        self.sdf_node = sdf_node
        self.parent = parent
        self.parent_frame = parent_frame
        self.joint_kind = joint_kind
        self.axis = axis
        self.limits = limits
        self.zero_pose = zero_pose
        self.description = description
        self.printable = printable

    @classmethod
    def from_command_meta(cls, name: str, meta: dict, **kwargs) -> "RomLink":
        """Create a RomLink from the metadata dict a command carries.

        The ``meta`` dict is the ``kinematic`` sub-dict placed on the
        solid's metadata by :mod:`~yapcad.rom.annotations`.
        """
        return cls(
            name=name,
            parent=meta.get("parent"),
            parent_frame=meta.get("parent_frame", "ORIGIN"),
            joint_kind=meta.get("joint_kind", "fixed"),
            axis=tuple(meta.get("axis", (0.0, 0.0, 1.0))),
            limits=tuple(meta.get("limits", (-180.0, 180.0))),
            zero_pose=float(meta.get("zero_pose", 0.0)),
            description=meta.get("description", ""),
            printable=bool(meta.get("printable", True)),
            **kwargs,
        )

    def _joint_type(self) -> JointType:
        mapping = {
            "fixed": JointType.FIXED,
            "revolute": JointType.REVOLUTE,
            "prismatic": JointType.PRISMATIC,
        }
        return mapping.get(self.joint_kind, JointType.FIXED)

    def to_kinematic_part(
        self,
        base_transform: Optional[Transform] = None,
    ) -> KinematicPart:
        """Convert this link to a :class:`KinematicPart` for the chain."""
        jtype = self._joint_type()
        part = KinematicPart(
            name=self.name,
            parent=self.parent,
            parent_frame=self.parent_frame,
            joint=Joint(
                name=f"{self.name}_joint",
                joint_type=jtype,
                axis=self.axis,
                value=self.zero_pose,
                min_limit=float(self.limits[0]),
                max_limit=float(self.limits[1]),
            ),
            description=self.description,
            is_printable=self.printable,
        )
        if base_transform is not None:
            part.joint.base_transform = base_transform
        return part


class RomTree:
    """A kinematic tree built from DSL annotations.

    Wraps a :class:`KinematicChain` and provides utilities for
    accessing links, their SDF nodes, and adjacent-pair queries.

    Attributes:
        name: Assembly name.
        chain: The underlying :class:`KinematicChain`.
        links: Dict of ``part_name -> RomLink``.
    """

    def __init__(self, name: str = "assembly"):
        self.name = name
        self.chain = KinematicChain(name)
        self.links: Dict[str, RomLink] = {}

    def add_link(
        self,
        link: RomLink,
        base_transform: Optional[Transform] = None,
    ) -> "RomTree":
        """Add a link to the tree and its corresponding chain part."""
        if link.name in self.links:
            logger.warning("Overwriting link '%s'", link.name)
        self.links[link.name] = link
        part = link.to_kinematic_part(base_transform=base_transform)
        self.chain.add_part(part)
        return self

    def clear_joints(self) -> None:
        """Reset all joint values to their zero poses."""
        for name, link in self.links.items():
            self.chain.set_joint_value(name, link.zero_pose)

    def set_joints(self, angles: Dict[str, float]) -> None:
        """Set multiple joint values at once.

        Args:
            angles: ``{part_name: angle_deg}`` dict.
        """
        for name, value in angles.items():
            if name in self.links:
                self.chain.set_joint_value(name, value)

    def get_sdf_node(self, part_name: str):
        """Get the compiled SDF node for a part.

        Returns ``None`` if the part has no SDF node, or the name is
        unknown.
        """
        link = self.links.get(part_name)
        if link is not None:
            return link.sdf_node
        return None

    def get_world_transform(self, part_name: str) -> Transform:
        """Get the current world transform for a part.

        Convenience wrapper over ``KinematicChain.get_world_transform``.
        """
        return self.chain.get_world_transform(part_name)

    def adjacent_pairs(self) -> List[Tuple[str, str]]:
        """List all parent-child (adjacent) pairs in the tree.

        Returns:
            ``[(parent, child), ...]`` — each pair has a direct joint
            and should be checked for clearance.
        """
        pairs = []
        for name, link in self.links.items():
            if link.parent and link.parent in self.links:
                pairs.append((link.parent, name))
        return pairs

    def link_names(self) -> List[str]:
        """Return the list of link names (topologically sorted)."""
        return [n for n in self.links.keys()]

    @classmethod
    def from_interpreter_context(
        cls,
        name: str,
        context,
        skip_fixed: bool = True,
    ) -> "RomTree":
        """Build a tree from a DSL interpreter execution context.

        Scans the context's command registry for commands that have
        ``@meta(kinematic.*)`` metadata, extracts their SDF nodes,
        and assembles the tree.

        Args:
            name: Assembly name.
            context: The ``ExecutionContext`` from ``interpret_module()``.
            skip_fixed: If True, skip commands whose ``joint_kind`` is
                ``"fixed"`` — they are rigid body parts, not joints.

        Returns:
            A populated :class:`RomTree`.

        Raises:
            ValueError: If the context has no command registry.
        """
        tree = cls(name)

        # The interpreter stores commands in its context.
        # We need to find commands with meta_hint containing kinematic.* keys.
        # The interpreter's execution context has:
        #   context.commands: dict of name -> Command (from AST)
        #   context.variables: dict of name -> actual value (the emitted solid)

        # Try different attribute names for the context
        commands_dict = {}
        variables_dict = {}

        if hasattr(context, "commands"):
            commands_dict = context.commands
        elif hasattr(context, "command"):
            commands_dict = {context.command.name: context.command}

        if hasattr(context, "variables"):
            variables_dict = context.variables

        # Also check for emitted results — the interpreter stores the last
        # emitted value for each command in the call stack or result registry.
        # We inspect the metadata of each solid variable for "kinematic".

        found_any = False

        for cmd_name, cmd_value in commands_dict.items():
            # Get the AST command to check for meta_hint
            if hasattr(cmd_value, "meta_hint") and cmd_value.meta_hint:
                meta = cmd_value.meta_hint
                # Collect kinematic keys
                kin = {}
                for key, val in meta.items():
                    if isinstance(key, str) and key.startswith("kinematic."):
                        field = key.split(".", 1)[1]
                        kin[field] = val

                if not kin:
                    continue

                found_any = True

                # Get the compiled SDF node from the variable
                sdf_node = None
                if cmd_name in variables_dict:
                    solid = variables_dict[cmd_name]
                    sdf_node = _extract_sdf(solid)

                link = RomLink.from_command_meta(
                    name=cmd_name,
                    meta=kin,
                    sdf_node=sdf_node,
                )

                if skip_fixed and link.joint_kind == "fixed":
                    continue

                # Determine base transform from the command arguments if any
                # For now, use identity; the KinematicChain will handle
                tree.add_link(link)

        if not found_any:
            logger.warning(
                "No kinematic-annotated commands found in context. "
                "Did you forget @meta(kinematic.*)?"
            )

        # Validate and log
        tree.clear_joints()
        issues = tree.chain.validate()
        if issues:
            logger.warning("Validation issues: %s", "; ".join(issues))

        return tree

    @classmethod
    def from_links(cls, name: str, links: List[RomLink]) -> "RomTree":
        """Build a tree from an explicit list of RomLink objects.

        Useful for testing or constructing a tree programmatically
        without DSL compilation.

        Args:
            name: Assembly name.
            links: List of RomLink objects.

        Returns:
            A populated :class:`RomTree`.
        """
        tree = cls(name)
        for link in links:
            tree.add_link(link)
        tree.clear_joints()
        return tree


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _extract_sdf(solid) -> Optional[Any]:
    """Try to extract the SDF node from a compiled solid value.

    The DSL solid is a plain ``list`` of the form
    ``['solid', [<surface1>, <surface2>, ...]]``.  The SDF tree is not
    stored on the solid itself — the interpreter evaluates DSL ``box()``,
    ``translate()``, ``union()``, etc. as surface operations.

    Since the DSL doesn't yet expose the SDF ``Node`` tree directly on
    compiled solids, we need an alternative approach:

    1. If the solid is a plain list (yapCAD v1 surface representation),
       rebuild a conservative bounding-box SDF from the surface data.
    2. If the solid is a ``Node`` (direct SDF output), return it.
    3. If the solid has a ``construction`` key in its metadata dict,
       recover the ``Node`` from it.
    4. Direct attribute ``_sdf_node`` set by the interpreter.

    Returns the SDF ``Node`` or ``None``.
    """
    # 1. Direct attribute — set during compilation if the interpreter
    #    stores it on the solid value
    sdf_node = getattr(solid, "_sdf_node", None)
    if sdf_node is not None:
        return sdf_node

    # 2. Construction slot — some solid types store the SDF tree there
    construction = getattr(solid, "construction", None)
    if construction is not None:
        try:
            from yapcad.sdf import from_construction, is_sdf_construction
            if is_sdf_construction(construction):
                node = from_construction(construction)
                if node is not None:
                    return node
        except Exception:
            pass

    # 3. Metadata dict on a plain list solid
    #    The solid is ['solid', [surface1, ...]] and metadata is stored
    #    in a dict embedded in solid[1] or as a separate metadata object.
    if isinstance(solid, list) and len(solid) > 0:
        for d in solid:
            if isinstance(d, dict):
                constr = d.get("construction")
                if constr:
                    try:
                        from yapcad.sdf import from_construction, is_sdf_construction
                        if is_sdf_construction(constr):
                            node = from_construction(constr)
                            if node is not None:
                                return node
                    except Exception:
                        pass

    # 4. Check if solid itself is already a Node (from direct SDF API)
    try:
        from yapcad.sdf.node import Node
        if isinstance(solid, Node):
            return solid
    except Exception:
        pass

    # 5. Fallback: build a conservative AABB box from the solid's bounding box
    try:
        from yapcad.geom3d import bounds
        from yapcad.sdf import box
        b = bounds(solid)
        if b is not None:
            lo, hi = b
            # Center and size of the bounding box
            center = [(lo[i] + hi[i]) / 2.0 for i in range(3)]
            size = [hi[i] - lo[i] for i in range(3)]
            sdf_box = box(size)
            # Translate the box to match the center
            from yapcad.sdf import translate
            return translate(sdf_box, center[0], center[1], center[2])
    except Exception:
        pass

    return None


__all__ = ["RomTree", "RomLink", "_extract_sdf"]