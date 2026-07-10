"""Assembly orchestrator for the yapCAD constraint-based assembly system.

This module provides the main Assembly class that manages parts, mates, and
constraints, and orchestrates the constraint solving and validation process.

The Assembly class is the primary user-facing interface for creating and
validating assemblies in yapCAD. It coordinates between the datum, mate, and
constraint subsystems to ensure assemblies meet design intent.

Key Features:
    - Add parts with initial transforms
    - Define mate relationships between parts
    - Add design constraints to validate intent
    - Solve mate constraints to compute part positions
    - Validate all design constraints
    - Export to URDF, Blender rigs, or other formats
    - Integration with existing kinematic_chain.py transforms

Example:
    >>> from yapcad.assembly import Assembly, PartDefinition, Datum, DatumType
    >>> from yapcad.assembly import Mate, MateType, Constraint, ConstraintType
    >>> import numpy as np
    >>>
    >>> # Create part definitions
    >>> motor = PartDefinition("DDSM115_MOTOR")
    >>> motor.add_datum(Datum("motor_axis", DatumType.AXIS,
    ...                       origin=(0, 0, 0), direction=(0, 1, 0)))
    >>> motor.add_datum(Datum("stator_face", DatumType.PLANE,
    ...                       origin=(0, 10, 0), normal=(0, 1, 0)))
    >>>
    >>> bracket = PartDefinition("WHEEL_BRACKET")
    >>> bracket.add_datum(Datum("bore_axis", DatumType.AXIS,
    ...                         origin=(0, 0, 0), direction=(0, 0, 1)))
    >>> bracket.add_datum(Datum("motor_interface", DatumType.PLANE,
    ...                         origin=(0, 0, 5), normal=(0, 0, 1)))
    >>>
    >>> # Create assembly
    >>> assembly = Assembly("wheel_module")
    >>> assembly.add_part(motor, name="motor_1")
    >>> assembly.add_part(bracket, name="bracket_1")
    >>>
    >>> # Add mates to position parts
    >>> assembly.add_mate(Mate("mount_mate", MateType.FLUSH,
    ...                        part1="bracket_1", datum1="motor_interface",
    ...                        part2="motor_1", datum2="stator_face"))
    >>> assembly.add_mate(Mate("axis_mate", MateType.CONCENTRIC,
    ...                        part1="bracket_1", datum1="bore_axis",
    ...                        part2="motor_1", datum2="motor_axis"))
    >>>
    >>> # Add design constraints
    >>> assembly.add_constraint(Constraint(
    ...     "motor_tangent", ConstraintType.TANGENT_TO_CIRCLE,
    ...     part="motor_1", datum="motor_axis",
    ...     center=(0, 0, 0), radius=124.5,
    ...     description="Motor axis must be tangent to chassis for rolling motion"
    ... ))
    >>>
    >>> # Validate assembly
    >>> result = assembly.validate()
    >>> if result.is_valid:
    ...     print("Assembly is valid!")
    ...     motor_transform = assembly.transforms["motor_1"]
    ... else:
    ...     for error in result.failed_constraints:
    ...         print(f"ERROR: {error}")
    ...     for warning in result.warnings:
    ...         print(f"WARNING: {warning}")

See Also:
    - datum.py: Datum feature definitions
    - mate.py: Mate constraint definitions
    - constraint.py: Design constraint definitions
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Any, Union
import numpy as np
import json
from pathlib import Path


# These will be imported from sibling modules once they exist
# For now, we define placeholder types for the implementation
try:
    from .datum import PartDefinition, Datum, DatumType
    from .mate import Mate, MateType
    from .constraint import Constraint, ConstraintResult, ConstraintType
except ImportError:
    # Fallback for development - these should be defined in separate modules
    PartDefinition = Any
    Datum = Any
    DatumType = Any
    Mate = Any
    MateType = Any
    Constraint = Any
    ConstraintResult = Any
    ConstraintType = Any


class AssemblyError(Exception):
    """Exception raised when assembly operations fail.

    This exception is raised for various assembly-related errors such as:
    - Invalid part names or references
    - Missing datums on parts
    - Invalid mate configurations
    - Constraint validation failures
    - Solver convergence failures

    Attributes:
        message: Description of the error
        assembly_name: Name of the assembly where error occurred
        part_name: Name of the part involved (if applicable)
        constraint_name: Name of the constraint that failed (if applicable)
    """

    def __init__(self, message: str, assembly_name: str = None,
                 part_name: str = None, constraint_name: str = None):
        self.message = message
        self.assembly_name = assembly_name
        self.part_name = part_name
        self.constraint_name = constraint_name
        super().__init__(self._format_message())

    def _format_message(self) -> str:
        """Format a detailed error message."""
        parts = []
        if self.assembly_name:
            parts.append(f"Assembly '{self.assembly_name}'")
        if self.part_name:
            parts.append(f"Part '{self.part_name}'")
        if self.constraint_name:
            parts.append(f"Constraint '{self.constraint_name}'")

        if parts:
            return f"{': '.join(parts)}: {self.message}"
        return self.message


@dataclass
class AssemblyValidationResult:
    """Result of validating an assembly's constraints.

    Contains the complete validation status including which constraints
    passed, failed, or produced warnings, along with detailed diagnostic
    information for each constraint.

    Attributes:
        is_valid: True if all constraints are satisfied, False otherwise
        constraint_results: Dictionary mapping constraint name to ConstraintResult
        failed_constraints: List of names of constraints that failed
        warnings: List of warning messages from constraint validation

    Example:
        >>> result = assembly.validate()
        >>> if not result.is_valid:
        ...     print(f"Found {len(result.failed_constraints)} failures")
        ...     for name in result.failed_constraints:
        ...         cr = result.constraint_results[name]
        ...         print(f"  {name}: {cr.message}")
    """

    is_valid: bool
    constraint_results: Dict[str, ConstraintResult] = field(default_factory=dict)
    failed_constraints: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    def __str__(self) -> str:
        """Human-readable summary of validation results."""
        if self.is_valid:
            return f"Assembly valid: {len(self.constraint_results)} constraints satisfied"
        else:
            n_failed = len(self.failed_constraints)
            n_warnings = len(self.warnings)
            return f"Assembly invalid: {n_failed} failed, {n_warnings} warnings"

    def report(self) -> str:
        """Generate detailed validation report."""
        lines = []
        lines.append("=" * 70)
        lines.append("Assembly Validation Report")
        lines.append("=" * 70)
        lines.append("")

        if self.is_valid:
            lines.append("STATUS: VALID - All constraints satisfied")
        else:
            lines.append(f"STATUS: INVALID - {len(self.failed_constraints)} constraints failed")

        lines.append("")
        lines.append(f"Total constraints: {len(self.constraint_results)}")
        lines.append(f"Passed: {len(self.constraint_results) - len(self.failed_constraints)}")
        lines.append(f"Failed: {len(self.failed_constraints)}")
        lines.append(f"Warnings: {len(self.warnings)}")
        lines.append("")

        # Failed constraints
        if self.failed_constraints:
            lines.append("FAILED CONSTRAINTS:")
            lines.append("-" * 70)
            for name in self.failed_constraints:
                result = self.constraint_results[name]
                lines.append(f"  [{name}]")
                lines.append(f"    {result.error_message}")
                if hasattr(result, 'error_value') and result.error_value > 0:
                    lines.append(f"    Error: {result.error_value:.3f}")
            lines.append("")

        # Warnings
        if self.warnings:
            lines.append("WARNINGS:")
            lines.append("-" * 70)
            for warning in self.warnings:
                lines.append(f"  - {warning}")
            lines.append("")

        # Passed constraints
        passed = [name for name in self.constraint_results
                  if name not in self.failed_constraints]
        if passed:
            lines.append("PASSED CONSTRAINTS:")
            lines.append("-" * 70)
            for name in passed:
                result = self.constraint_results[name]
                lines.append(f"  [OK] {name}")
            lines.append("")

        lines.append("=" * 70)
        return "\n".join(lines)


class Assembly:
    """Main assembly orchestrator for the yapCAD constraint system.

    The Assembly class manages a collection of parts with their transforms,
    mate relationships that position parts, and design constraints that
    validate assembly intent.

    Workflow:
        1. Create assembly: `assembly = Assembly("my_assembly")`
        2. Add parts: `assembly.add_part(part_def, name="part_1")`
        3. Add mates: `assembly.add_mate(mate)`
        4. Add constraints: `assembly.add_constraint(constraint)`
        5. Validate: `result = assembly.validate()`

    The assembly can also import transforms from existing kinematic chains
    and export to various formats (URDF, Blender, etc.).

    Attributes:
        name: Unique identifier for this assembly
        parts: Dictionary mapping part name to PartDefinition
        transforms: Dictionary mapping part name to 4x4 transform matrix
        mates: List of Mate objects defining part relationships
        constraints: List of Constraint objects to validate

    Integration with yapCAD:
        - Transforms are 4x4 numpy arrays compatible with yapcad.xform.Matrix
        - Can import transforms from kinematic_chain.py JSON exports
        - Datums transform correctly through coordinate systems

    Example:
        >>> # Create assembly
        >>> assembly = Assembly("robot_wheel")
        >>>
        >>> # Add parts
        >>> assembly.add_part(motor_def, name="motor")
        >>> assembly.add_part(bracket_def, name="bracket")
        >>>
        >>> # Define relationships
        >>> assembly.add_mate(Mate("align", MateType.CONCENTRIC,
        ...                        "bracket", "hole", "motor", "shaft"))
        >>>
        >>> # Validate design intent
        >>> assembly.add_constraint(Constraint(
        ...     "tangent", ConstraintType.TANGENT_TO_CIRCLE,
        ...     "motor", "axis", center=(0,0,0), radius=124.5
        ... ))
        >>>
        >>> # Check if valid
        >>> result = assembly.validate()
        >>> if not result.is_valid:
        ...     assembly.validate_and_raise()  # Raises AssemblyError
    """

    def __init__(self, name: str):
        """Initialize a new assembly.

        Args:
            name: Unique identifier for this assembly
        """
        self.name = name
        self.parts: Dict[str, PartDefinition] = {}
        self.transforms: Dict[str, np.ndarray] = {}
        self.mates: List[Mate] = []
        self.constraints: List[Constraint] = []
        # Load cases for structural FEA setup. Keyed by LoadCase.id.
        # Added 2026-05-20 (see assembly/load_case.py history).
        self.load_cases: Dict[str, Any] = {}
        # Bolt patterns keyed by (parent_part, child_part) tuple. Mirrors the
        # Mechatron Interface.bolt_pattern Optional<BoltPattern>. yapCAD users
        # call add_bolt_pattern() to attach one to a mate group; the exporter
        # emits it onto the corresponding Mechatron Interface.
        self.bolt_patterns: Dict[tuple, Any] = {}
        self._solved = False

    def add_part(self, part: PartDefinition, name: str = None,
                 transform: np.ndarray = None) -> None:
        """Add a part to the assembly.

        Args:
            part: PartDefinition with datum features
            name: Name for this part instance (default: part.name)
            transform: Initial 4x4 transform matrix (default: identity)

        Raises:
            AssemblyError: If part name already exists in assembly

        Example:
            >>> motor = PartDefinition("DDSM115")
            >>> assembly.add_part(motor, name="motor_1")
            >>> # Or with initial position
            >>> T = np.eye(4)
            >>> T[0:3, 3] = [100, 0, 0]  # Translate 100mm in X
            >>> assembly.add_part(motor, name="motor_2", transform=T)
        """
        part_name = name if name is not None else part.name

        if part_name in self.parts:
            raise AssemblyError(
                f"Part name '{part_name}' already exists in assembly",
                assembly_name=self.name,
                part_name=part_name
            )

        self.parts[part_name] = part
        self.transforms[part_name] = (
            transform.copy() if transform is not None else np.eye(4)
        )
        self._solved = False

    def add_mate(self, mate: Mate) -> None:
        """Add a mate relationship between two parts.

        Mates define geometric relationships between datum features on
        different parts. They are used to position parts relative to
        each other.

        Args:
            mate: Mate object defining the relationship

        Raises:
            AssemblyError: If referenced parts or datums don't exist,
                          or if mate is invalid for the given datums

        Example:
            >>> # Align motor shaft with bracket hole
            >>> assembly.add_mate(Mate(
            ...     "shaft_alignment", MateType.CONCENTRIC,
            ...     part1="bracket", datum1="bore_axis",
            ...     part2="motor", datum2="shaft_axis"
            ... ))
        """
        # Validate parts exist
        if mate.part1 not in self.parts:
            raise AssemblyError(
                f"Part '{mate.part1}' not found in assembly",
                assembly_name=self.name,
                part_name=mate.part1
            )
        if mate.part2 not in self.parts:
            raise AssemblyError(
                f"Part '{mate.part2}' not found in assembly",
                assembly_name=self.name,
                part_name=mate.part2
            )

        # Validate datums exist
        part1 = self.parts[mate.part1]
        part2 = self.parts[mate.part2]

        if mate.datum1 not in part1.datums:
            raise AssemblyError(
                f"Datum '{mate.datum1}' not found on part '{mate.part1}'",
                assembly_name=self.name,
                part_name=mate.part1
            )
        if mate.datum2 not in part2.datums:
            raise AssemblyError(
                f"Datum '{mate.datum2}' not found on part '{mate.part2}'",
                assembly_name=self.name,
                part_name=mate.part2
            )

        # Validate mate is compatible with datum types
        d1 = part1.get_datum(mate.datum1)
        d2 = part2.get_datum(mate.datum2)
        issues = mate.validate(d1, d2)
        if issues:
            raise AssemblyError(
                f"Invalid mate '{mate.name}': {'; '.join(issues)}",
                assembly_name=self.name
            )

        self.mates.append(mate)
        self._solved = False

    def add_constraint(self, constraint: Constraint) -> None:
        """Add a design constraint to validate assembly intent.

        Constraints validate high-level design requirements that go beyond
        simple positioning. They ensure the assembly meets its intended
        functional requirements.

        Args:
            constraint: Constraint object to validate

        Raises:
            AssemblyError: If referenced part or datum doesn't exist

        Example:
            >>> # Ensure motor axis is tangent to chassis
            >>> assembly.add_constraint(Constraint(
            ...     "motor_tangent", ConstraintType.TANGENT_TO_CIRCLE,
            ...     part="motor_1", datum="motor_axis",
            ...     center=(0, 0, 0), radius=124.5,
            ...     description="Motor must be tangent for rolling motion"
            ... ))
        """
        if constraint.part not in self.parts:
            raise AssemblyError(
                f"Part '{constraint.part}' not found in assembly",
                assembly_name=self.name,
                part_name=constraint.part,
                constraint_name=constraint.name
            )

        part = self.parts[constraint.part]
        if constraint.datum not in part.datums:
            raise AssemblyError(
                f"Datum '{constraint.datum}' not found on part '{constraint.part}'",
                assembly_name=self.name,
                part_name=constraint.part,
                constraint_name=constraint.name
            )

        self.constraints.append(constraint)

    # ------------------------------------------------------------------
    # Load case registry (added 2026-05-20 — see assembly/load_case.py)
    # ------------------------------------------------------------------
    def add_load_case(self, load_case) -> None:
        """Register a structural LoadCase with this assembly.

        Args:
            load_case: a yapcad.assembly.load_case.LoadCase instance

        Raises:
            AssemblyError: if the LoadCase references a part not in the
                assembly (when attach.part != "_assembly")
        """
        attach_part = load_case.attach.part
        if attach_part != "_assembly" and attach_part not in self.parts:
            raise AssemblyError(
                f"LoadCase '{load_case.id}' attaches to part '{attach_part}' "
                f"which is not in the assembly",
                assembly_name=self.name,
                part_name=attach_part,
            )
        if load_case.id in self.load_cases:
            raise AssemblyError(
                f"LoadCase id '{load_case.id}' already registered",
                assembly_name=self.name,
            )
        self.load_cases[load_case.id] = load_case

    def get_load_case(self, lc_id: str):
        """Return the LoadCase with the given id, or None if missing."""
        return self.load_cases.get(lc_id)

    # ------------------------------------------------------------------
    # Bolt pattern registry (added 2026-05-20)
    # ------------------------------------------------------------------
    def add_bolt_pattern(self, parent_part: str, child_part: str, bolt_pattern) -> None:
        """Attach a BoltPattern to the interface between parent and child.

        The pattern is keyed by the (parent, child) tuple, mirroring how
        the Mechatron exporter groups mates into a single Interface. If
        the same pair already has a bolt_pattern registered, this raises.
        """
        if parent_part not in self.parts:
            raise AssemblyError(
                f"parent part '{parent_part}' not in assembly",
                assembly_name=self.name, part_name=parent_part,
            )
        if child_part not in self.parts:
            raise AssemblyError(
                f"child part '{child_part}' not in assembly",
                assembly_name=self.name, part_name=child_part,
            )
        key = (parent_part, child_part)
        if key in self.bolt_patterns:
            raise AssemblyError(
                f"BoltPattern already registered for {parent_part} -> {child_part}",
                assembly_name=self.name,
            )
        self.bolt_patterns[key] = bolt_pattern

    def get_bolt_pattern(self, parent_part: str, child_part: str):
        """Return the BoltPattern for an interface, or None if not set."""
        return self.bolt_patterns.get((parent_part, child_part))

    def get_transformed_datum(self, part_name: str, datum_name: str) -> Datum:
        """Get a datum feature transformed to world coordinates.

        Applies the part's current transform to a datum to obtain its
        position and orientation in world coordinates.

        Args:
            part_name: Name of the part
            datum_name: Name of the datum on that part

        Returns:
            Datum transformed to world coordinates

        Raises:
            AssemblyError: If part or datum not found

        Example:
            >>> # Get motor axis in world coordinates
            >>> motor_axis_world = assembly.get_transformed_datum(
            ...     "motor_1", "motor_axis"
            ... )
            >>> print(f"Axis origin: {motor_axis_world.origin}")
            >>> print(f"Axis direction: {motor_axis_world.direction}")
        """
        if part_name not in self.parts:
            raise AssemblyError(
                f"Part '{part_name}' not found in assembly",
                assembly_name=self.name,
                part_name=part_name
            )

        part = self.parts[part_name]
        if datum_name not in part.datums:
            raise AssemblyError(
                f"Datum '{datum_name}' not found on part '{part_name}'",
                assembly_name=self.name,
                part_name=part_name
            )

        datum_local = part.get_datum(datum_name)
        transform = self.transforms[part_name]
        return datum_local.transform(transform)

    def validate(self) -> AssemblyValidationResult:
        """Validate all design constraints against current part transforms.

        Evaluates each constraint in the assembly and returns a comprehensive
        validation result. This does not solve mates - it validates the
        current state of the assembly.

        Returns:
            AssemblyValidationResult with detailed status

        Example:
            >>> result = assembly.validate()
            >>> if result.is_valid:
            ...     print("Assembly is valid!")
            ... else:
            ...     print(result.report())
        """
        constraint_results = {}
        failed_constraints = []
        warnings = []

        for constraint in self.constraints:
            # Get datum in world coordinates
            datum_world = self.get_transformed_datum(
                constraint.part, constraint.datum
            )

            # Evaluate constraint
            result = constraint.evaluate(datum_world)
            constraint_results[constraint.name] = result

            # Track failures and warnings
            if not result.passed:
                if constraint.severity == "error":
                    failed_constraints.append(constraint.name)
                elif constraint.severity == "warning":
                    warnings.append(f"{constraint.name}: {result.error_message}")

        is_valid = len(failed_constraints) == 0

        return AssemblyValidationResult(
            is_valid=is_valid,
            constraint_results=constraint_results,
            failed_constraints=failed_constraints,
            warnings=warnings
        )

    def validate_and_raise(self) -> None:
        """Validate assembly and raise AssemblyError if invalid.

        Convenience method for strict validation where you want an
        exception on any constraint failure.

        Raises:
            AssemblyError: If any constraint fails

        Example:
            >>> try:
            ...     assembly.validate_and_raise()
            ...     print("Assembly is valid, proceeding...")
            ... except AssemblyError as e:
            ...     print(f"Invalid assembly: {e}")
        """
        result = self.validate()
        if not result.is_valid:
            # Create detailed error message
            lines = [f"Assembly '{self.name}' validation failed:"]
            for name in result.failed_constraints:
                cr = result.constraint_results[name]
                lines.append(f"  - {name}: {cr.message}")

            raise AssemblyError(
                "\n".join(lines),
                assembly_name=self.name
            )

    def get_degrees_of_freedom(self) -> Dict[str, int]:
        """Calculate remaining degrees of freedom for each part after mates.

        Analyzes the mate constraints to determine how many degrees of
        freedom remain unconstrained for each part. A fully constrained
        part has 0 DOF.

        Returns:
            Dictionary mapping part name to DOF count (0-6)

        Note:
            This is a simplified implementation. A full implementation would
            analyze the constraint graph to detect over/under-constrained
            situations and redundant constraints.

        Example:
            >>> dof = assembly.get_degrees_of_freedom()
            >>> for part, dof_count in dof.items():
            ...     if dof_count > 0:
            ...         print(f"{part} has {dof_count} DOF remaining")
        """
        # Each part starts with 6 DOF (3 translation, 3 rotation)
        dof = {name: 6 for name in self.parts}

        # This is a simplified approximation
        # A full implementation would build the constraint Jacobian
        for mate in self.mates:
            # Estimate DOF removed by each mate type
            if mate.mate_type == MateType.COINCIDENT:
                # Removes 3 DOF (position locked)
                dof[mate.part2] = max(0, dof[mate.part2] - 3)
            elif mate.mate_type == MateType.CONCENTRIC:
                # Removes 2 DOF (can slide and rotate along axis)
                dof[mate.part2] = max(0, dof[mate.part2] - 2)
            elif mate.mate_type == MateType.FLUSH:
                # Removes 3 DOF (orientation locked, can slide in plane)
                dof[mate.part2] = max(0, dof[mate.part2] - 3)
            elif mate.mate_type == MateType.PARALLEL:
                # Removes 2 DOF (2 rotation axes locked)
                dof[mate.part2] = max(0, dof[mate.part2] - 2)
            elif mate.mate_type == MateType.PERPENDICULAR:
                # Removes 1 DOF (1 rotation axis locked)
                dof[mate.part2] = max(0, dof[mate.part2] - 1)

        return dof

    def import_transforms_from_kinematic_chain(self, json_path: Union[str, Path]) -> None:
        """Import part transforms from a kinematic_chain.py JSON export.

        The kinematic chain system computes transforms through a hierarchical
        tree. This method imports those transforms into the assembly so
        constraints can be validated against them.

        Args:
            json_path: Path to JSON file with kinematic chain transforms

        Raises:
            AssemblyError: If JSON cannot be loaded or part names don't match

        Example:
            >>> # Export from kinematic_chain.py
            >>> chain.export_transforms("chain_transforms.json")
            >>>
            >>> # Import into assembly
            >>> assembly.import_transforms_from_kinematic_chain(
            ...     "chain_transforms.json"
            ... )
            >>>
            >>> # Now validate against those transforms
            >>> result = assembly.validate()
        """
        try:
            with open(json_path, 'r') as f:
                data = json.load(f)
        except Exception as e:
            raise AssemblyError(
                f"Failed to load kinematic chain JSON: {e}",
                assembly_name=self.name
            )

        # Expected format: {"parts": {"part_name": {"transform": [...]}}}
        if "parts" not in data:
            raise AssemblyError(
                "Invalid kinematic chain JSON: missing 'parts' key",
                assembly_name=self.name
            )

        for part_name, part_data in data["parts"].items():
            if part_name not in self.parts:
                # Skip parts not in this assembly
                continue

            if "transform" not in part_data:
                raise AssemblyError(
                    f"Part '{part_name}' missing 'transform' in JSON",
                    assembly_name=self.name,
                    part_name=part_name
                )

            # Transform should be 4x4 matrix (16 elements or 4x4 array)
            transform_data = part_data["transform"]
            if isinstance(transform_data, list):
                if len(transform_data) == 16:
                    # Flat list - reshape to 4x4
                    transform = np.array(transform_data).reshape(4, 4)
                elif len(transform_data) == 4 and all(len(row) == 4 for row in transform_data):
                    # Already 4x4
                    transform = np.array(transform_data)
                else:
                    raise AssemblyError(
                        f"Invalid transform format for part '{part_name}'",
                        assembly_name=self.name,
                        part_name=part_name
                    )
            else:
                raise AssemblyError(
                    f"Transform for part '{part_name}' must be a list",
                    assembly_name=self.name,
                    part_name=part_name
                )

            self.transforms[part_name] = transform

        self._solved = True  # Mark as solved since transforms are imported

    def export_to_urdf(self, output_path: Union[str, Path]) -> None:
        """Export assembly to URDF format for ROS integration.

        URDF (Unified Robot Description Format) is used by ROS for robot
        kinematics and dynamics. This export enables simulation and control
        of yapCAD assemblies in ROS.

        Args:
            output_path: Path where URDF file will be written

        Note:
            This is currently a stub implementation. A full implementation
            would generate proper URDF XML with links, joints, and meshes.

        Example:
            >>> assembly.export_to_urdf("robot.urdf")
            >>> # Use in ROS: roslaunch robot_description display.launch
        """
        # STUB: Full implementation would generate URDF XML
        # including links, joints, inertial properties, collision geometry
        raise NotImplementedError(
            "URDF export not yet implemented. "
            "See ROS URDF documentation for format details."
        )

    def export_to_blender_rig(self, output_path: Union[str, Path]) -> None:
        """Export assembly to Blender armature rig for animation.

        Creates a Blender Python script that sets up an armature with
        bones positioned according to the assembly transforms. Useful
        for visualizing kinematics and creating animations.

        Args:
            output_path: Path where Blender Python script will be written

        Note:
            This is currently a stub implementation. A full implementation
            would generate a Blender Python script that creates bones,
            constraints, and animations.

        Example:
            >>> assembly.export_to_blender_rig("robot_rig.py")
            >>> # In Blender: File > Import > Python Script > robot_rig.py
        """
        # STUB: Full implementation would generate Blender Python script
        # creating armature with bones, constraints, and IK chains
        raise NotImplementedError(
            "Blender rig export not yet implemented. "
            "See Blender Python API documentation for bpy.types.Armature."
        )

    def __repr__(self) -> str:
        """String representation of assembly."""
        return (
            f"Assembly(name='{self.name}', "
            f"parts={len(self.parts)}, "
            f"mates={len(self.mates)}, "
            f"constraints={len(self.constraints)})"
        )

    def report(self) -> str:
        """Generate a comprehensive human-readable assembly report.

        Returns:
            Multi-line string with assembly status and validation results

        Example:
            >>> print(assembly.report())
            Assembly: robot_wheel
            ====================================
            Parts: 2
              - motor_1 (DDSM115_MOTOR)
              - bracket_1 (WHEEL_BRACKET)

            Mates: 2
              - shaft_alignment (CONCENTRIC)
              - mount_surface (FLUSH)

            Constraints: 1
              - motor_tangent (TANGENT_TO_CIRCLE)

            Validation: VALID
        """
        lines = []
        lines.append(f"Assembly: {self.name}")
        lines.append("=" * 60)

        # Parts section
        lines.append(f"Parts: {len(self.parts)}")
        for name, part in self.parts.items():
            lines.append(f"  - {name} ({part.name})")
        lines.append("")

        # Mates section
        lines.append(f"Mates: {len(self.mates)}")
        for mate in self.mates:
            lines.append(
                f"  - {mate.name} ({mate.mate_type.value}): "
                f"{mate.part1}.{mate.datum1} <-> {mate.part2}.{mate.datum2}"
            )
        lines.append("")

        # Constraints section
        lines.append(f"Constraints: {len(self.constraints)}")
        for constraint in self.constraints:
            lines.append(
                f"  - {constraint.name} ({constraint.constraint_type.value}): "
                f"{constraint.part}.{constraint.datum}"
            )
        lines.append("")

        # Validation section
        try:
            result = self.validate()
            lines.append("Validation:")
            lines.append("-" * 60)
            if result.is_valid:
                lines.append("  STATUS: VALID")
            else:
                lines.append(f"  STATUS: INVALID ({len(result.failed_constraints)} failures)")
                for name in result.failed_constraints:
                    cr = result.constraint_results[name]
                    lines.append(f"    FAILED: {name} - {cr.message}")

            if result.warnings:
                lines.append(f"  WARNINGS: {len(result.warnings)}")
                for warning in result.warnings:
                    lines.append(f"    - {warning}")
        except Exception as e:
            lines.append(f"Validation: ERROR - {e}")

        lines.append("")
        lines.append("=" * 60)
        return "\n".join(lines)
