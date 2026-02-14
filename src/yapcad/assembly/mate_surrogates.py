"""Mate Surrogates for COTS (Commercial Off-The-Shelf) Parts.

This module defines explicit reference frame definitions for COTS parts that lack
native CAD mate features. These "mate surrogates" provide unambiguous interface
definitions for kinematic chain assembly.

===============================================================================
PROBLEM STATEMENT
===============================================================================

COTS parts like Dynamixel servos and DDSM115 hub motors come with STL/STEP files
that have geometry but not explicit mate features (datum planes, reference axes).
When assembling them into a kinematic chain, we need clear interface definitions:

1. Where is the output shaft axis?
2. Where is the mounting face?
3. What's the reference frame for orientation?
4. Where are bolt patterns relative to the origin?

Without explicit surrogates, these interfaces are implicitly defined in code
(e.g., kinematic_chain.py), making them:
- Hard to verify
- Prone to errors when updating
- Scattered across multiple files

===============================================================================
SOLUTION: MATE SURROGATES
===============================================================================

A mate surrogate defines a reference frame attached to a COTS part with:
- origin: The reference point (typically mounting face center)
- x_axis, y_axis, z_axis: Orthonormal frame axes
- features: Named geometric features with positions and orientations

This is similar to how professional CAD systems (SolidWorks, Creo) handle
imported geometry - you create reference geometry on top of the dumb solid.

===============================================================================
EXISTING GEOMETRY ANALYSIS
===============================================================================

FINDING: The existing DSL files and kinematic_chain.py already contain most of
the information needed for proper mating, but it's IMPLICIT rather than EXPLICIT:

1. DDSM115 Motor (ddsm115_motor.dsl):
   - EXCELLENT documentation in comments (lines 1-71)
   - Origin at mounting face (Z=0), wheel in +Z, stator body in -Z
   - Motor rotation axis is ALONG the Y-axis in STL coordinates (not Z!)
   - Triangular boss pattern: 3x at 120deg, radius 10.5mm
   - Bolt circle: 15.2mm diameter, 3x M2.5 at 120deg
   - ISSUE: The DSL defines the geometry but kinematic_chain.py has to
     manually apply transforms (Ry(-90) @ Rz(30)) to orient correctly

2. Dynamixel Servos (XH540, XH430, XL330):
   - STL files exist but no DSL definitions
   - Interface specs are in gearbox_mates.py (ServoHornSpec class)
   - Horn bolt circles, protrusion dimensions are documented
   - Servo body dimensions are in scara_arm.dsl (lines 41-86)
   - ISSUE: No explicit coordinate system definition in STL files

3. Planetary Gearboxes (planetary_gearbox.dsl):
   - Sun gear interface matches servo horn bolt patterns
   - Center counterbore clears horn protrusion
   - Output carrier has bolt pattern for next stage
   - ISSUE: Mating relies on correct Z-stack positioning

4. Existing Mate Files:
   - gearbox_mates.py: Excellent! Defines ServoHornSpec and SunGearHubSpec
   - scara_assembly_mates.py: Complete interface catalog
   - HOWEVER: These define CONSTRAINTS between parts, not the reference
     frames on individual parts

===============================================================================
CONCLUSION
===============================================================================

The existing geometry IS sufficient for unambiguous mating, BUT the reference
frame information is scattered and implicit. Mate surrogates consolidate this
information into a single source of truth.

RECOMMENDATION:
- For Dynamixel servos: Minimal surrogates needed (bolt patterns well-documented)
- For DDSM115: Essential - the Y-axis rotation is non-obvious and critical
- For gearboxes: Already well-defined in planetary_gearbox.dsl, just need to
  consolidate frame definitions

===============================================================================
Copyright (c) 2026 yapCAD contributors
License: MIT
===============================================================================
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional, Any
from enum import Enum


# =============================================================================
# FEATURE TYPES
# =============================================================================

class FeatureType(Enum):
    """Types of geometric features on COTS parts."""
    AXIS = "axis"               # Rotation axis (point + direction)
    PLANE = "plane"             # Flat face (point + normal)
    CIRCLE = "circle"           # Bolt circle or hole (center + normal + radius)
    POLYGON = "polygon"         # Polygonal feature (vertices)
    CYLINDER = "cylinder"       # Cylindrical feature (center + axis + radius + height)


# =============================================================================
# MATE SURROGATE DATA STRUCTURES
# =============================================================================

@dataclass
class Feature:
    """A named geometric feature on a part.

    Attributes:
        name: Feature name (e.g., "stator_face", "bolt_circle")
        feature_type: Type of geometric feature
        point: Reference point [x, y, z] in part local coordinates
        direction: Direction vector [x, y, z] (for axes, normals)
        radius: Radius for circular features (mm)
        count: Number of items (for hole patterns)
        angle_offset_deg: Angular offset of first item from +X (degrees)
        vertices: List of vertex positions (for polygons)
        height: Height for cylindrical features (mm)
        description: Human-readable description
    """
    name: str
    feature_type: FeatureType
    point: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    direction: Tuple[float, float, float] = (0.0, 0.0, 1.0)
    radius: Optional[float] = None
    count: Optional[int] = None
    angle_offset_deg: float = 0.0
    vertices: Optional[List[Tuple[float, float, float]]] = None
    height: Optional[float] = None
    description: str = ""

    def get_hole_positions(self) -> List[Tuple[float, float, float]]:
        """Get positions of holes in a bolt circle pattern."""
        if self.feature_type != FeatureType.CIRCLE or self.radius is None or self.count is None:
            return []

        positions = []
        for i in range(self.count):
            angle_rad = math.radians(self.angle_offset_deg + i * (360.0 / self.count))
            x = self.point[0] + self.radius * math.cos(angle_rad)
            y = self.point[1] + self.radius * math.sin(angle_rad)
            z = self.point[2]
            positions.append((x, y, z))
        return positions


@dataclass
class MateSurrogate:
    """Complete reference frame definition for a COTS part.

    The surrogate defines a local coordinate system and named features
    that can be used for assembly mating.

    Attributes:
        part_name: Name of the COTS part (e.g., "DDSM115")
        description: Part description
        origin: Origin point in the part's STL/native coordinate system
        x_axis: X-axis direction vector
        y_axis: Y-axis direction vector
        z_axis: Z-axis direction vector
        features: Dictionary of named features
        stl_file: Path to STL file (relative to project output directory)
        step_file: Path to STEP file if available
        notes: Additional notes about the part
    """
    part_name: str
    description: str = ""

    # Reference frame (in STL/native coordinates)
    origin: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    x_axis: Tuple[float, float, float] = (1.0, 0.0, 0.0)
    y_axis: Tuple[float, float, float] = (0.0, 1.0, 0.0)
    z_axis: Tuple[float, float, float] = (0.0, 0.0, 1.0)

    # Named features
    features: Dict[str, Feature] = field(default_factory=dict)

    # File references
    stl_file: Optional[str] = None
    step_file: Optional[str] = None

    # Notes
    notes: str = ""

    def add_feature(self, feature: Feature) -> None:
        """Add a feature to this surrogate."""
        self.features[feature.name] = feature

    def get_feature(self, name: str) -> Feature:
        """Get a feature by name."""
        if name not in self.features:
            raise KeyError(f"Feature '{name}' not found on part '{self.part_name}'")
        return self.features[name]

    def get_transform_to_standard_frame(self) -> Dict[str, Any]:
        """Get transform parameters to convert from STL coords to standard frame.

        Returns dict with:
        - rotation_rpy: (roll, pitch, yaw) in degrees
        - translation: (x, y, z) offset
        """
        # This would compute the actual transform - implementation depends on
        # how the axes differ from identity
        return {
            "origin": self.origin,
            "x_axis": self.x_axis,
            "y_axis": self.y_axis,
            "z_axis": self.z_axis
        }


# =============================================================================
# DDSM115 HUB MOTOR SURROGATE
# =============================================================================

def create_ddsm115_surrogate() -> MateSurrogate:
    """Create mate surrogate for Waveshare DDSM115 hub motor.

    CRITICAL COORDINATE SYSTEM (from STEP file analysis):
    - The motor rotation axis is the STL Y-AXIS (not Z!)
    - Tire (rotor/wheel) is at Y < 0 (faces -Y direction)
    - Stator (mounting face) is at Y > 0 (faces +Y direction)

    This is non-intuitive because most parts have their rotation axis along Z.
    The kinematic_chain.py applies Ry(-90) to map motor Y-axis to tangent direction.

    PHYSICAL STRUCTURE:
    - Hub motor: outer rotor (tire), inner stator (mounting hub)
    - Tire OD: 100.7mm, Width: 50mm
    - Motor body diameter: 67mm
    - Triangular anti-rotation bosses on stator side

    STL ORIGIN:
    - At the interface between tire and stator body
    - Y=0 is approximately at the tire-stator boundary
    - Stator mounting face is at Y=+10mm
    - Tire center is at Y=-26.5mm (50mm wide tire)
    """
    surrogate = MateSurrogate(
        part_name="DDSM115",
        description="Waveshare DDSM115 Direct Drive Hub Motor (100.7mm wheel, 67mm stator)",

        # Origin at stator mounting face center (where it bolts to bracket)
        # This is at Y=+10mm in STL coordinates, but we define surrogate origin here
        origin=(0.0, 10.0, 0.0),

        # CRITICAL: Motor rotation axis is along STL Y-axis!
        # We define the surrogate frame to make this clear:
        # - Surrogate X = STL X (radial in XZ plane)
        # - Surrogate Y = STL Z (perpendicular to motor axis)
        # - Surrogate Z = STL Y (motor rotation axis, toward wheel)
        x_axis=(1.0, 0.0, 0.0),
        y_axis=(0.0, 0.0, 1.0),  # STL Z
        z_axis=(0.0, 1.0, 0.0),  # STL Y (motor axis)

        stl_file="cots/ddsm115_motor_official.stl",
        step_file="cots/DDSM115 STEP.step",

        notes="""
CRITICAL: Motor rotation axis is STL Y-axis (NOT Z)!
- To orient motor with axis tangent to chassis, apply Ry(-90)
- Tire (rotor) faces -Y in STL coords, stator faces +Y
- Triangular boss pattern requires Rz(30) to align with bracket pocket
"""
    )

    # Motor rotation axis (the primary feature for mating)
    surrogate.add_feature(Feature(
        name="rotor_axis",
        feature_type=FeatureType.AXIS,
        point=(0.0, 0.0, 0.0),  # At motor center
        direction=(0.0, 1.0, 0.0),  # STL Y-axis
        description="Motor rotation axis (wheel spins around this)"
    ))

    # Stator mounting face (flat face where motor bolts to bracket)
    surrogate.add_feature(Feature(
        name="stator_face",
        feature_type=FeatureType.PLANE,
        point=(0.0, 10.0, 0.0),  # At Y=+10mm in STL coords
        direction=(0.0, 1.0, 0.0),  # Normal points +Y (away from wheel)
        description="Stator mounting face (bolts to wheel arm bracket)"
    ))

    # Bolt circle (3x M2.5 at 120deg spacing)
    surrogate.add_feature(Feature(
        name="bolt_circle",
        feature_type=FeatureType.CIRCLE,
        point=(0.0, 10.0, 0.0),  # At stator face
        direction=(0.0, 1.0, 0.0),  # Normal to face
        radius=15.2 / 2.0,  # 7.6mm radius = 15.2mm diameter
        count=3,
        angle_offset_deg=0.0,  # First hole at +X
        description="M2.5 mounting bolt circle (3x at 120deg, 15.2mm dia)"
    ))

    # Triangular anti-rotation boss pattern
    # These are the 3 bosses that fit into a triangular pocket on the bracket
    boss_radius = 10.5  # Bosses at 10.5mm radius from center
    boss_angles = [0.0, 120.0, 240.0]
    boss_vertices = []
    for angle_deg in boss_angles:
        angle_rad = math.radians(angle_deg)
        x = boss_radius * math.cos(angle_rad)
        z = boss_radius * math.sin(angle_rad)  # In XZ plane at Y=stator_face
        boss_vertices.append((x, 0.0, z))  # Y=0 relative to stator face

    surrogate.add_feature(Feature(
        name="triangular_boss",
        feature_type=FeatureType.POLYGON,
        point=(0.0, 0.0, 0.0),  # At stator back (wiring side)
        direction=(0.0, -1.0, 0.0),  # Bosses extend toward -Y (into bracket)
        vertices=boss_vertices,
        height=10.0,  # Bosses are 10mm tall
        description="Triangular anti-rotation bosses (3x at 120deg, 10.5mm radius)"
    ))

    # Tire (wheel/rotor) - the part that contacts the tube wall
    surrogate.add_feature(Feature(
        name="tire",
        feature_type=FeatureType.CYLINDER,
        point=(0.0, -26.5, 0.0),  # Tire center at Y=-26.5mm
        direction=(0.0, 1.0, 0.0),  # Axis along Y
        radius=100.7 / 2.0,  # 50.35mm radius = 100.7mm OD
        height=50.0,  # 50mm wide tire
        description="Tire outer cylinder (contacts tube inner wall)"
    ))

    # Motor body (stator housing)
    surrogate.add_feature(Feature(
        name="motor_body",
        feature_type=FeatureType.CYLINDER,
        point=(0.0, 5.0, 0.0),  # Centered in stator body
        direction=(0.0, 1.0, 0.0),  # Axis along Y
        radius=67.0 / 2.0,  # 33.5mm radius = 67mm diameter
        height=15.5,  # Stator body depth
        description="Stator body housing"
    ))

    # Pilot boss (center locating feature)
    surrogate.add_feature(Feature(
        name="pilot_boss",
        feature_type=FeatureType.CYLINDER,
        point=(0.0, 10.0, 0.0),  # At stator face
        direction=(0.0, 1.0, 0.0),  # Projects outward
        radius=8.0 / 2.0,  # 4mm radius = 8mm diameter
        height=2.0,  # 2mm projection
        description="Center pilot boss (locating feature)"
    ))

    return surrogate


# =============================================================================
# DYNAMIXEL SERVO SURROGATES
# =============================================================================

def create_xh540_surrogate() -> MateSurrogate:
    """Create mate surrogate for Robotis Dynamixel XH540-W270 servo.

    XH540 is the largest Dynamixel X-series servo used in this robot.
    It drives Axis 1 (shoulder yaw) through a 6:1 planetary gearbox.

    COORDINATE SYSTEM:
    - STL origin at servo body center
    - Z-axis is servo output shaft axis (horn faces +Z)
    - Servo body is symmetric about Z-axis in XY plane

    PHYSICAL DIMENSIONS (from datasheet):
    - Body: 33.5 x 58.5 x 44mm (W x D x H)
    - Horn bolt circle: 18mm diameter, 4x M2.5 at 45deg offset
    - Output shaft at one end, back shell at other
    """
    surrogate = MateSurrogate(
        part_name="XH540",
        description="Robotis Dynamixel XH540-W270 Smart Actuator",

        # Origin at center of horn face (output interface)
        origin=(0.0, 0.0, 22.0),  # Half of 44mm height

        # Standard orientation: Z is output shaft axis
        x_axis=(1.0, 0.0, 0.0),
        y_axis=(0.0, 1.0, 0.0),
        z_axis=(0.0, 0.0, 1.0),

        stl_file="cots/xh540_cots.stl",

        notes="""
XH540-W270 is the shoulder yaw actuator.
- Stall torque: 10.6 Nm @ 24V
- No-load speed: 46 rpm
- Horn bolt circle: 18mm, 4x M2.5 at 45deg offset
- Connects to AXIS1_SUN_GEAR via hub interface
"""
    )

    # Output shaft axis
    surrogate.add_feature(Feature(
        name="output_axis",
        feature_type=FeatureType.AXIS,
        point=(0.0, 0.0, 22.0),  # At horn face
        direction=(0.0, 0.0, 1.0),  # +Z
        description="Servo output shaft axis"
    ))

    # Horn mounting face
    surrogate.add_feature(Feature(
        name="horn_face",
        feature_type=FeatureType.PLANE,
        point=(0.0, 0.0, 24.0),  # Top of horn (22mm body half + ~2mm horn projection)
        direction=(0.0, 0.0, 1.0),  # Normal faces +Z
        description="Horn top face (sun gear hub mounts here)"
    ))

    # Horn bolt circle (4x M2.5 at 45deg offset)
    surrogate.add_feature(Feature(
        name="horn_bolt_circle",
        feature_type=FeatureType.CIRCLE,
        point=(0.0, 0.0, 24.0),  # At horn face
        direction=(0.0, 0.0, 1.0),
        radius=18.0 / 2.0,  # 9mm radius = 18mm diameter
        count=4,
        angle_offset_deg=45.0,  # First hole at 45deg (standard Dynamixel pattern)
        description="Horn bolt circle (4x M2.5 at 18mm dia, 45deg offset)"
    ))

    # Horn center protrusion (must clear sun gear counterbore)
    surrogate.add_feature(Feature(
        name="horn_protrusion",
        feature_type=FeatureType.CYLINDER,
        point=(0.0, 0.0, 24.0),  # At horn face
        direction=(0.0, 0.0, 1.0),
        radius=6.0 / 2.0,  # 3mm radius = 6mm diameter
        height=2.5,  # 2.5mm projection
        description="Horn center protrusion (fits into sun gear counterbore)"
    ))

    # Servo body envelope
    surrogate.add_feature(Feature(
        name="body_envelope",
        feature_type=FeatureType.CYLINDER,
        point=(0.0, 0.0, 0.0),  # At body center
        direction=(0.0, 0.0, 1.0),
        radius=33.5 / 2.0,  # Approximate as cylinder
        height=44.0,
        description="Servo body envelope (for clearance checking)"
    ))

    return surrogate


def create_xh430_surrogate() -> MateSurrogate:
    """Create mate surrogate for Robotis Dynamixel XH430-W350 servo.

    XH430 is the mid-size servo used for Axis 2 and Axis 3.
    Smaller than XH540 but still high-torque.
    """
    surrogate = MateSurrogate(
        part_name="XH430",
        description="Robotis Dynamixel XH430-W350 Smart Actuator",

        origin=(0.0, 0.0, 17.0),  # Half of 34mm height
        x_axis=(1.0, 0.0, 0.0),
        y_axis=(0.0, 1.0, 0.0),
        z_axis=(0.0, 0.0, 1.0),

        stl_file="cots/xh430_cots.stl",

        notes="""
XH430-W350 is used for Axis 2 (shoulder pitch) and Axis 3 (elbow).
- Stall torque: 3.4 Nm @ 24V
- No-load speed: 46 rpm
- Horn bolt circle: 16mm, 4x M2.5 at 45deg offset
"""
    )

    surrogate.add_feature(Feature(
        name="output_axis",
        feature_type=FeatureType.AXIS,
        point=(0.0, 0.0, 17.0),
        direction=(0.0, 0.0, 1.0),
        description="Servo output shaft axis"
    ))

    surrogate.add_feature(Feature(
        name="horn_face",
        feature_type=FeatureType.PLANE,
        point=(0.0, 0.0, 19.0),  # 17mm body half + ~2mm horn
        direction=(0.0, 0.0, 1.0),
        description="Horn top face"
    ))

    surrogate.add_feature(Feature(
        name="horn_bolt_circle",
        feature_type=FeatureType.CIRCLE,
        point=(0.0, 0.0, 19.0),
        direction=(0.0, 0.0, 1.0),
        radius=16.0 / 2.0,  # 8mm radius = 16mm diameter
        count=4,
        angle_offset_deg=45.0,
        description="Horn bolt circle (4x M2.5 at 16mm dia, 45deg offset)"
    ))

    surrogate.add_feature(Feature(
        name="horn_protrusion",
        feature_type=FeatureType.CYLINDER,
        point=(0.0, 0.0, 19.0),
        direction=(0.0, 0.0, 1.0),
        radius=6.0 / 2.0,
        height=2.5,
        description="Horn center protrusion"
    ))

    return surrogate


def create_xl330_surrogate() -> MateSurrogate:
    """Create mate surrogate for Robotis Dynamixel XL330-M288 servo.

    XL330 is the compact servo used for Axis 4 (wrist rotation).
    Lightweight and suitable for end effector mounting.
    """
    surrogate = MateSurrogate(
        part_name="XL330",
        description="Robotis Dynamixel XL330-M288 Smart Actuator",

        origin=(0.0, 0.0, 13.0),  # Half of 26mm height
        x_axis=(1.0, 0.0, 0.0),
        y_axis=(0.0, 1.0, 0.0),
        z_axis=(0.0, 0.0, 1.0),

        stl_file=None,  # No STL currently available

        notes="""
XL330-M288 is the wrist rotation actuator.
- Stall torque: 0.52 Nm @ 5V
- No-load speed: 61 rpm
- Horn bolt circle: 8mm, 4x M2 at 45deg offset
- Compact and lightweight for end effector
"""
    )

    surrogate.add_feature(Feature(
        name="output_axis",
        feature_type=FeatureType.AXIS,
        point=(0.0, 0.0, 13.0),
        direction=(0.0, 0.0, 1.0),
        description="Servo output shaft axis"
    ))

    surrogate.add_feature(Feature(
        name="horn_face",
        feature_type=FeatureType.PLANE,
        point=(0.0, 0.0, 14.5),  # 13mm body half + ~1.5mm horn
        direction=(0.0, 0.0, 1.0),
        description="Horn top face"
    ))

    surrogate.add_feature(Feature(
        name="horn_bolt_circle",
        feature_type=FeatureType.CIRCLE,
        point=(0.0, 0.0, 14.5),
        direction=(0.0, 0.0, 1.0),
        radius=8.0 / 2.0,  # 4mm radius = 8mm diameter
        count=4,
        angle_offset_deg=45.0,
        description="Horn bolt circle (4x M2 at 8mm dia, 45deg offset)"
    ))

    surrogate.add_feature(Feature(
        name="horn_protrusion",
        feature_type=FeatureType.CYLINDER,
        point=(0.0, 0.0, 14.5),
        direction=(0.0, 0.0, 1.0),
        radius=4.0 / 2.0,  # 2mm radius = 4mm diameter (smaller)
        height=2.0,
        description="Horn center protrusion"
    ))

    return surrogate


# =============================================================================
# GEARBOX SURROGATES
# =============================================================================

def create_sun_gear_surrogate(
    axis_name: str,
    sun_teeth: int,
    module_mm: float,
    face_width: float,
    hub_height: float,
    bolt_circle: float,
    servo_type: str
) -> MateSurrogate:
    """Create mate surrogate for a planetary sun gear.

    Sun gears have a flanged hub that mates with servo horns.
    The hub extends BELOW the gear teeth (in -Z).

    Args:
        axis_name: Axis identifier (e.g., "AXIS1")
        sun_teeth: Number of teeth on sun gear
        module_mm: Gear module in mm
        face_width: Gear face width in mm
        hub_height: Total hub height (extends below gear)
        bolt_circle: Bolt circle diameter for servo horn interface
        servo_type: Servo model (e.g., "XH540", "XH430", "XL330")
    """
    pitch_d = sun_teeth * module_mm
    root_d = sun_teeth * module_mm - 2.5 * module_mm

    surrogate = MateSurrogate(
        part_name=f"{axis_name}_SUN_GEAR",
        description=f"Sun gear {sun_teeth}T @ {module_mm}mm module for {servo_type}",

        # Origin at gear base (Z=0), hub extends to Z=-hub_height
        origin=(0.0, 0.0, 0.0),
        x_axis=(1.0, 0.0, 0.0),
        y_axis=(0.0, 1.0, 0.0),
        z_axis=(0.0, 0.0, 1.0),

        notes=f"""
Sun gear for {axis_name}:
- Teeth: {sun_teeth}T @ {module_mm}mm module
- Pitch diameter: {pitch_d:.2f}mm
- Root diameter: {root_d:.2f}mm
- Face width: {face_width}mm
- Hub height: {hub_height}mm (extends below gear in -Z)
- Bolt circle: {bolt_circle}mm (matches {servo_type} horn)
"""
    )

    # Gear rotation axis
    surrogate.add_feature(Feature(
        name="gear_axis",
        feature_type=FeatureType.AXIS,
        point=(0.0, 0.0, 0.0),
        direction=(0.0, 0.0, 1.0),
        description="Sun gear rotation axis"
    ))

    # Gear top face (meshes with planets)
    surrogate.add_feature(Feature(
        name="gear_top",
        feature_type=FeatureType.PLANE,
        point=(0.0, 0.0, face_width),
        direction=(0.0, 0.0, 1.0),
        description="Gear top face"
    ))

    # Hub bottom face (contacts servo horn)
    surrogate.add_feature(Feature(
        name="hub_face",
        feature_type=FeatureType.PLANE,
        point=(0.0, 0.0, -hub_height),
        direction=(0.0, 0.0, -1.0),  # Faces downward toward servo
        description="Hub bottom face (contacts servo horn)"
    ))

    # Hub bolt circle
    surrogate.add_feature(Feature(
        name="hub_bolt_circle",
        feature_type=FeatureType.CIRCLE,
        point=(0.0, 0.0, -hub_height),
        direction=(0.0, 0.0, -1.0),
        radius=bolt_circle / 2.0,
        count=4,
        angle_offset_deg=45.0,  # Match Dynamixel horn pattern
        description=f"Hub bolt circle ({bolt_circle}mm dia, matches {servo_type})"
    ))

    # Center counterbore (clears horn protrusion)
    counterbore_dia = 8.0 if servo_type != "XL330" else 6.0
    counterbore_depth = 3.0 if servo_type != "XL330" else 2.5
    surrogate.add_feature(Feature(
        name="center_counterbore",
        feature_type=FeatureType.CYLINDER,
        point=(0.0, 0.0, -hub_height),
        direction=(0.0, 0.0, 1.0),  # Extends upward into hub
        radius=counterbore_dia / 2.0,
        height=counterbore_depth,
        description="Center counterbore (clears horn protrusion)"
    ))

    return surrogate


# =============================================================================
# SURROGATE REGISTRY
# =============================================================================

def get_all_surrogates() -> Dict[str, MateSurrogate]:
    """Get dictionary of all COTS part mate surrogates."""
    surrogates = {}

    # Hub motor
    surrogates["DDSM115"] = create_ddsm115_surrogate()

    # Dynamixel servos
    surrogates["XH540"] = create_xh540_surrogate()
    surrogates["XH430"] = create_xh430_surrogate()
    surrogates["XL330"] = create_xl330_surrogate()

    # Sun gears (from planetary_gearbox.dsl specifications)
    surrogates["AXIS1_SUN_GEAR"] = create_sun_gear_surrogate(
        axis_name="AXIS1",
        sun_teeth=18,
        module_mm=0.75,
        face_width=10.0,
        hub_height=6.0,
        bolt_circle=18.0,
        servo_type="XH540"
    )

    surrogates["AXIS2_SUN_GEAR"] = create_sun_gear_surrogate(
        axis_name="AXIS2",
        sun_teeth=16,
        module_mm=0.75,
        face_width=8.0,
        hub_height=6.0,
        bolt_circle=16.0,
        servo_type="XH430"
    )

    surrogates["AXIS3_SUN_GEAR"] = create_sun_gear_surrogate(
        axis_name="AXIS3",
        sun_teeth=18,
        module_mm=1.0,
        face_width=6.0,
        hub_height=6.0,
        bolt_circle=16.0,
        servo_type="XH430"
    )

    surrogates["AXIS4_SUN_GEAR"] = create_sun_gear_surrogate(
        axis_name="AXIS4",
        sun_teeth=12,
        module_mm=0.75,
        face_width=6.0,
        hub_height=5.0,
        bolt_circle=8.0,
        servo_type="XL330"
    )

    return surrogates


def print_surrogate_summary() -> None:
    """Print a summary of all mate surrogates."""
    surrogates = get_all_surrogates()

    print("=" * 80)
    print("MATE SURROGATES FOR COTS PARTS")
    print("=" * 80)
    print(f"\nTotal surrogates defined: {len(surrogates)}")

    for name, surrogate in surrogates.items():
        print(f"\n{'-' * 80}")
        print(f"Part: {name}")
        print(f"Description: {surrogate.description}")
        print(f"Origin: {surrogate.origin}")
        print(f"Features: {len(surrogate.features)}")

        for feat_name, feature in surrogate.features.items():
            feat_type = feature.feature_type.value
            print(f"  - {feat_name} ({feat_type}): {feature.description}")


def validate_surrogates() -> List[str]:
    """Validate all mate surrogates for consistency.

    Returns:
        List of validation issues (empty if all valid)
    """
    issues = []
    surrogates = get_all_surrogates()

    for name, surrogate in surrogates.items():
        # Check that axes are orthonormal
        x = np.array(surrogate.x_axis)
        y = np.array(surrogate.y_axis)
        z = np.array(surrogate.z_axis)

        # Check unit vectors
        for axis_name, axis in [("x", x), ("y", y), ("z", z)]:
            length = np.linalg.norm(axis)
            if abs(length - 1.0) > 0.001:
                issues.append(f"{name}: {axis_name}_axis is not unit length ({length:.4f})")

        # Check orthogonality
        xy_dot = abs(np.dot(x, y))
        xz_dot = abs(np.dot(x, z))
        yz_dot = abs(np.dot(y, z))

        if xy_dot > 0.001:
            issues.append(f"{name}: x and y axes not orthogonal (dot={xy_dot:.4f})")
        if xz_dot > 0.001:
            issues.append(f"{name}: x and z axes not orthogonal (dot={xz_dot:.4f})")
        if yz_dot > 0.001:
            issues.append(f"{name}: y and z axes not orthogonal (dot={yz_dot:.4f})")

        # Check that features have required fields
        for feat_name, feature in surrogate.features.items():
            if feature.feature_type == FeatureType.CIRCLE:
                if feature.radius is None:
                    issues.append(f"{name}.{feat_name}: circle missing radius")
                if feature.count is None:
                    issues.append(f"{name}.{feat_name}: circle missing count")

    return issues


# Need numpy for validation
try:
    import numpy as np
except ImportError:
    # If numpy not available, skip validation that needs it
    def validate_surrogates() -> List[str]:
        return ["numpy not available for validation"]


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    print_surrogate_summary()

    print("\n" + "=" * 80)
    print("VALIDATION")
    print("=" * 80)

    issues = validate_surrogates()
    if issues:
        print("\nValidation issues found:")
        for issue in issues:
            print(f"  [WARN] {issue}")
    else:
        print("\n[OK] All surrogates validated successfully")

    # Show DDSM115 details as example
    print("\n" + "=" * 80)
    print("DDSM115 DETAILED ANALYSIS")
    print("=" * 80)

    ddsm = get_all_surrogates()["DDSM115"]
    print(f"\nPart: {ddsm.part_name}")
    print(f"Description: {ddsm.description}")
    print(f"\nCoordinate System:")
    print(f"  Origin: {ddsm.origin}")
    print(f"  X-axis: {ddsm.x_axis}")
    print(f"  Y-axis: {ddsm.y_axis}")
    print(f"  Z-axis: {ddsm.z_axis}")
    print(f"\nNotes: {ddsm.notes}")

    print("\nFeatures:")
    for name, feat in ddsm.features.items():
        print(f"\n  {name} ({feat.feature_type.value}):")
        print(f"    Point: {feat.point}")
        print(f"    Direction: {feat.direction}")
        if feat.radius:
            print(f"    Radius: {feat.radius}mm")
        if feat.count:
            print(f"    Count: {feat.count}")
        if feat.height:
            print(f"    Height: {feat.height}mm")
        print(f"    Description: {feat.description}")

    # Example: Get bolt hole positions for DDSM115
    print("\n" + "=" * 80)
    print("EXAMPLE: DDSM115 BOLT HOLE POSITIONS")
    print("=" * 80)

    bolt_circle = ddsm.get_feature("bolt_circle")
    positions = bolt_circle.get_hole_positions()
    print(f"\nBolt circle: {bolt_circle.radius * 2}mm diameter, {bolt_circle.count} holes")
    for i, pos in enumerate(positions):
        angle = bolt_circle.angle_offset_deg + i * (360.0 / bolt_circle.count)
        print(f"  Hole {i+1}: ({pos[0]:.3f}, {pos[1]:.3f}, {pos[2]:.3f}) @ {angle}deg")
