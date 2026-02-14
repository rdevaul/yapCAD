=====================
Assembly System
=====================

The yapCAD assembly system provides a complete constraint-based framework for defining and validating multi-body mechanical assemblies. It combines kinematic chain modeling, mate constraint solving, collision detection, and interactive visualization into a unified workflow.

.. note::

   The assembly system was contributed by Jeremy Mika in 2026 to enable procedural definition of complex mechatronic assemblies with parametric constraints, motion simulation, and automated collision validation.

Overview
========

The assembly system addresses the core challenges of procedural CAD:

* **Datum-driven positioning**: Parts reference named geometric features (points, axes, planes) rather than hardcoded transforms
* **Constraint solving**: Mate constraints (FLUSH, CONCENTRIC, REVOLUTE) compute 6DOF transforms automatically
* **Kinematic modeling**: Tree-structured assemblies with joints for articulation and motion planning
* **Collision detection**: Multi-method validation (BREP, mesh, AABB) with interface volume support for allowed overlaps
* **Interactive visualization**: Multi-viewport VTK viewer with REST API and WebSocket control

Key Features
------------

* **Zero hardcoded offsets**: All part positions derive from geometric constraints and datum features
* **Pluggable geometry**: GeometryProvider protocol works with any STEP/STL source
* **Interface volumes**: Explicit definitions for gear mesh, threads, and other designed overlaps
* **Multi-method collision**: Automatic fallback from BREP → mesh → AABB based on available dependencies
* **Production-ready**: Used for complex planetary gearbox assemblies with 50+ parts

Installation
============

Core Dependencies
-----------------

The assembly system requires numpy for transform math:

.. code-block:: bash

   pip install numpy

Optional Dependencies
---------------------

For BREP-based collision detection (most accurate):

.. code-block:: bash

   conda install -c conda-forge pythonocc-core

For mesh-based collision detection:

.. code-block:: bash

   pip install trimesh

For visualization and REST API:

.. code-block:: bash

   pip install vtk flask flask-socketio flask-cors

.. note::

   The system gracefully degrades when optional dependencies are unavailable. Without pythonocc, it falls back to mesh detection. Without trimesh, it uses AABB-only detection.

Core Concepts
=============

Datums and Datum Features
--------------------------

A **datum** is a named geometric reference on a part that provides an explicit anchor for assembly constraints. Datums replace manual transform calculations with declarative feature-based positioning.

Datum Types
~~~~~~~~~~~

.. code-block:: python

   from yapcad.assembly.datum import Datum, DatumType, PartDefinition
   from yapcad.geom import point, vect

   # POINT: Single location in 3D space
   mount_point = Datum(
       name="mounting_center",
       datum_type=DatumType.POINT,
       origin=point(0, 0, 5),
       description="Center of mounting interface"
   )

   # AXIS: Infinite line (origin + direction)
   shaft_axis = Datum(
       name="drive_shaft",
       datum_type=DatumType.AXIS,
       origin=point(0, 0, 0),
       direction=vect(0, 0, 1, 0),  # Z-axis
       description="Motor shaft rotation axis"
   )

   # PLANE: Infinite plane (origin + normal)
   mount_face = Datum(
       name="stator_face",
       datum_type=DatumType.PLANE,
       origin=point(0, 10, 0),
       normal=vect(0, 1, 0, 0),  # +Y normal
       description="Mounting surface for servo"
   )

   # CIRCLE: Bolt circle pattern (center + normal + radius)
   bolt_circle = Datum(
       name="mounting_holes",
       datum_type=DatumType.CIRCLE,
       origin=point(0, 10, 0),
       normal=vect(0, 1, 0, 0),
       radius=15.2,  # mm
       description="4xM2.5 mounting bolt pattern"
   )

   # FRAME: Full coordinate system (origin + X/Y axes, Z computed)
   tool_frame = Datum(
       name="tool_center_point",
       datum_type=DatumType.FRAME,
       origin=point(0, 0, 120),
       x_axis=vect(1, 0, 0, 0),
       y_axis=vect(0, 1, 0, 0),
       description="End effector TCP"
   )

Part Definitions
~~~~~~~~~~~~~~~~

Parts are defined with their datum features:

.. code-block:: python

   # Define a servo motor part
   servo = PartDefinition(
       name="XH540_SERVO",
       geometry_source="cots/xh540.step",
       is_printable=False,
       material="Aluminum"
   )

   # Add datum features
   servo.add_datum(Datum(
       name="stator_face",
       datum_type=DatumType.PLANE,
       origin=point(0, 0, 0),
       normal=vect(0, 0, 1, 0),
       description="Servo mounting face"
   ))

   servo.add_datum(Datum(
       name="output_shaft",
       datum_type=DatumType.AXIS,
       origin=point(0, 0, 36),
       direction=vect(0, 0, 1, 0),
       description="Output horn rotation axis"
   ))

   servo.add_datum(Datum(
       name="horn_bolt_circle",
       datum_type=DatumType.CIRCLE,
       origin=point(0, 0, 36),
       normal=vect(0, 0, 1, 0),
       radius=9.0,
       description="Horn mounting bolts"
   ))

   # Validate datums
   issues = servo.validate_datums()
   if issues:
       print("Datum validation issues:", issues)

Mate Constraints
----------------

A **mate** defines a geometric relationship between two parts that constrains their relative position and/or orientation.

MateType Enum
~~~~~~~~~~~~~

The ``MateType`` enum defines standard constraint types based on industry CAD practices:

**Geometric Constraints:**

* ``COINCIDENT``: Points/planes/circles share same location (face-to-face contact)
* ``CONCENTRIC``: Axes are colinear (shaft-in-bore alignment)
* ``PARALLEL``: Directions remain parallel
* ``PERPENDICULAR``: Directions at 90 degrees
* ``TANGENT``: Surfaces remain tangent
* ``DISTANCE``: Fixed offset between features
* ``ANGLE``: Fixed angular relationship

**Standard Joints:**

* ``RIGID``: No relative motion (welded/bolted)
* ``REVOLUTE``: Rotation about single axis (hinge)
* ``PRISMATIC``: Translation along single axis (slider)
* ``CYLINDRICAL``: Rotation + translation on axis
* ``SPHERICAL``: Ball joint (3 rotational DOF)
* ``PLANAR``: Motion in a plane (2 translation + 1 rotation)

**Compound Joints:**

* ``PIN_SLOT``: Translation along slot + rotation about pin
* ``UNIVERSAL``: Two perpendicular rotation axes (U-joint)
* ``SCREW``: Coupled rotation and translation (threaded rod)

**Coupled Mates:**

* ``GEAR``: Coupled rotation with ratio
* ``RACK_PINION``: Couples rotation to translation
* ``CAM``: Follower constrained to cam profile

Creating Mates
~~~~~~~~~~~~~~

.. code-block:: python

   from yapcad.assembly.mate import Mate, MateType, MateLimits

   # Face-to-face mounting constraint
   mount_mate = Mate(
       name="servo_to_bracket",
       mate_type=MateType.COINCIDENT,
       part_a="BRACKET",
       datum_a="servo_mount_face",
       part_b="XH540_SERVO",
       datum_b="stator_face",
       offset=0.0  # No gap
   )

   # Shaft alignment constraint
   shaft_mate = Mate(
       name="shaft_alignment",
       mate_type=MateType.CONCENTRIC,
       part_a="BRACKET",
       datum_a="bore_axis",
       part_b="XH540_SERVO",
       datum_b="output_shaft"
   )

   # Revolute joint with limits
   joint_mate = Mate(
       name="shoulder_pitch",
       mate_type=MateType.REVOLUTE,
       part_a="BASE",
       datum_a="shoulder_axis",
       part_b="UPPER_ARM",
       datum_b="arm_root_axis",
       axis=[0, 1, 0, 0],  # Y-axis rotation
       limits=MateLimits(
           min_value=-1.57,  # -90 degrees (radians)
           max_value=1.57,   # +90 degrees
           max_velocity=2.0,  # rad/s
           max_effort=100.0   # N*m torque
       )
   )

Convenience Functions
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from yapcad.assembly.mate import (
       create_revolute_mate,
       create_prismatic_mate,
       create_gear_mate
   )

   # Create revolute joint
   elbow = create_revolute_mate(
       name="elbow_flex",
       part_a="upper_arm",
       datum_a="elbow_axis",
       part_b="forearm",
       datum_b="forearm_root",
       min_angle=0.0,
       max_angle=3.14,  # 180 degrees
       friction=0.02,
       damping=0.1
   )

   # Create linear slide
   actuator = create_prismatic_mate(
       name="z_axis",
       part_a="base",
       datum_a="rail_axis",
       part_b="carriage",
       datum_b="slider_axis",
       min_position=0.0,
       max_position=100.0,  # 100mm travel
       max_velocity=50.0    # mm/s
   )

   # Create gear coupling (3:1 reduction)
   gearbox = create_gear_mate(
       name="motor_to_output",
       part_a="motor_shaft",
       datum_a="motor_axis",
       part_b="output_shaft",
       datum_b="output_axis",
       ratio=1.0/3.0,  # Output rotates 1/3 speed
       reverse=True    # Opposite direction
   )

The MateConstraintSolver
-------------------------

The ``MateConstraintSolver`` computes 6DOF transforms from mate constraints, eliminating hardcoded offsets.

Basic Usage
~~~~~~~~~~~

.. code-block:: python

   from yapcad.assembly.solver import MateConstraintSolver
   from yapcad.assembly.datum_registry import DatumRegistry

   # Register datum sources
   DatumRegistry.register_source("servo", servo_datums)
   DatumRegistry.register_source("bracket", bracket_datums)

   # Create solver
   solver = MateConstraintSolver(
       tolerance=0.001,      # 1 micron position tolerance
       angle_tolerance=0.1   # 0.1 degree angular tolerance
   )

   # Solve single mate
   result = solver.solve_mate(mount_mate)

   if result.success:
       child_transform = result.transform  # 4x4 numpy array
       print(f"Residual error: {result.residual}")
   else:
       print(f"Failed: {result.error_message}")

Assembly Solving
~~~~~~~~~~~~~~~~

.. code-block:: python

   from yapcad.assembly.solver import solve_mate_chain

   # Define mate chain (parent → child order)
   mates = [
       base_to_link1_mate,
       link1_to_link2_mate,
       link2_to_tool_mate
   ]

   # Solve sequentially
   world_transforms = solve_mate_chain(
       mates,
       base_transform=np.eye(4)  # World origin
   )

   # Access computed transforms
   link1_world_tf = world_transforms["LINK1"]
   link2_world_tf = world_transforms["LINK2"]
   tool_world_tf = world_transforms["TOOL"]

Transform Validation
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Validate a computed transform
   validation = solver.validate_transform(transform)

   if not validation.valid:
       print("Transform issues:")
       for msg in validation.error_messages:
           print(f"  - {msg}")
       print(f"Orthonormality error: {validation.position_error:.2e}")
       print(f"Determinant error: {validation.orientation_error:.2e}")

Kinematic Chains and Transforms
--------------------------------

The ``yapcad.kinematics`` package provides tree-structured assemblies with joints for motion simulation.

Transform Class
~~~~~~~~~~~~~~~

.. code-block:: python

   from yapcad.kinematics import Transform

   # Create transforms
   T1 = Transform.from_translation(10, 20, 30)
   R = Transform.from_axis_angle((0, 0, 1), 45)  # 45° about Z
   T2 = Transform.from_rpy(0, 0, 45)  # Roll-pitch-yaw

   # Combine transforms
   combined = T1 @ R @ T2

   # Extract components
   pos = combined.position()  # (x, y, z)
   quat = combined.quaternion()  # (w, x, y, z)
   rpy = combined.rpy()  # (roll, pitch, yaw) in degrees

   # Invert
   inv = combined.inverse()

   # Convert to numpy
   mat = combined.matrix()  # 4x4 numpy array

Joint Types
~~~~~~~~~~~

.. code-block:: python

   from yapcad.kinematics import Joint, JointType

   # Fixed joint (no motion)
   fixed = Joint("mount", JointType.FIXED)

   # Revolute joint (rotation about axis)
   revolute = Joint(
       name="shoulder",
       joint_type=JointType.REVOLUTE,
       axis=(0, 1, 0),  # Y-axis
       min_limit=-90,
       max_limit=90
   )
   revolute.set_value(45)  # Set current angle

   # Prismatic joint (linear translation)
   slider = Joint(
       name="actuator",
       joint_type=JointType.PRISMATIC,
       axis=(0, 0, 1),  # Z-axis
       min_limit=0,
       max_limit=100  # 100mm stroke
   )
   slider.set_value(50)  # 50mm extended

Kinematic Chain
~~~~~~~~~~~~~~~

.. code-block:: python

   from yapcad.kinematics import KinematicChain, KinematicPart

   # Create chain
   chain = KinematicChain("robot_arm")

   # Add base (attached to world)
   base = KinematicPart(
       name="BASE",
       parent=None,
       joint=Joint("base_fixed", JointType.FIXED)
   )
   chain.add_part(base)

   # Add shoulder link
   shoulder = KinematicPart(
       name="SHOULDER",
       parent="BASE",
       parent_frame="SHOULDER_MOUNT",
       joint=Joint("shoulder_pitch", JointType.REVOLUTE, axis=(0, 1, 0))
   )
   chain.add_part(shoulder)

   # Add elbow link
   elbow = KinematicPart(
       name="ELBOW",
       parent="SHOULDER",
       parent_frame="ELBOW_MOUNT",
       joint=Joint("elbow_flex", JointType.REVOLUTE, axis=(0, 1, 0))
   )
   chain.add_part(elbow)

   # Set joint angles
   chain.set_joint_value("SHOULDER", 30)  # degrees
   chain.set_joint_value("ELBOW", 45)

   # Get world transforms
   shoulder_tf = chain.get_world_transform("SHOULDER")
   elbow_tf = chain.get_world_transform("ELBOW")

   # Export to JSON
   chain.export_json("positions.json")

Collision Detection
===================

The ``yapcad.collision`` package provides multi-method collision detection with interface volume support.

GeometryProvider Protocol
--------------------------

Implement the ``GeometryProvider`` protocol to supply part geometry:

.. code-block:: python

   from yapcad.collision import GeometryProvider
   from pathlib import Path

   class FileBasedProvider:
       """Provides geometry from STEP/STL files."""

       def __init__(self, base_dir: Path):
           self.base_dir = base_dir

       def get_geometry(self, part_name: str):
           """Return path to STEP or STL file."""
           # Prefer STEP for BREP detection
           step_path = self.base_dir / f"{part_name}.step"
           if step_path.exists():
               return str(step_path)

           # Fall back to STL for mesh detection
           stl_path = self.base_dir / f"{part_name}.stl"
           if stl_path.exists():
               return str(stl_path)

           return None

Detection Methods
-----------------

The detector automatically selects the best available method:

1. **BREP Detection** (requires pythonocc-core):

   * Uses exact boolean intersection via ``BRepAlgoAPI_Common``
   * Computes precise collision volume
   * Most accurate but requires STEP files

2. **Mesh Detection** (requires trimesh):

   * Point sampling and containment checks
   * Works with STL files
   * Fast and reliable for most assemblies

3. **AABB Detection** (always available):

   * Axis-aligned bounding box overlap
   * Fast pre-filter for non-colliding pairs
   * Least accurate but very fast

Basic Usage
-----------

.. code-block:: python

   from yapcad.collision import CollisionDetector

   # Create detector
   provider = FileBasedProvider(Path("output/assembly"))
   detector = CollisionDetector(
       geometry_provider=provider,
       min_collision_volume=0.1,  # Ignore tiny overlaps < 0.1 mm³
       verbose=True
   )

   # Check single pair
   result = detector.check_collision(
       part_a="SERVO",
       transform_a=servo_tf,
       part_b="BRACKET",
       transform_b=bracket_tf
   )

   if result.collides:
       print(f"COLLISION: {result.part_a} <-> {result.part_b}")
       print(f"  Volume: {result.collision_volume} mm³")
       print(f"  Method: {result.method.value}")
       print(f"  Penetration: {result.penetration_depth} mm")

Assembly-Wide Detection
-----------------------

.. code-block:: python

   # Check entire assembly
   world_transforms = {
       "BASE": base_tf,
       "SERVO": servo_tf,
       "BRACKET": bracket_tf,
       "LINK1": link1_tf,
       "LINK2": link2_tf
   }

   results = detector.check_assembly(world_transforms)

   # Filter to actual collisions (not interface volumes)
   collisions = [r for r in results if r.is_error]

   if collisions:
       print(f"Found {len(collisions)} collisions:")
       for result in collisions:
           print(f"  {result.part_a} <-> {result.part_b}")
           print(f"    Volume: {result.collision_volume:.3f} mm³")
   else:
       print("No collisions detected!")

Interface Volumes
-----------------

Interface volumes define **allowed overlaps** for designed interfaces like gear mesh and screw threads.

Gear Mesh Interface
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from yapcad.collision import (
       InterfaceRegistry,
       GearMeshInterface,
       create_planetary_gearbox_interfaces
   )

   # Create registry
   registry = InterfaceRegistry()

   # Define sun gear teeth interface
   sun_interface = GearMeshInterface(
       name="axis1_sun_teeth",
       part_name="AXIS1_SUN_GEAR",
       module=0.75,           # Gear module
       teeth=18,              # Tooth count
       pressure_angle=20.0,   # degrees
       face_width=10.0        # mm
   )
   registry.register(sun_interface)

   # Define planet gear teeth
   planet_interface = GearMeshInterface(
       name="axis1_planet_teeth",
       part_name="AXIS1_PLANET_GEAR",
       module=0.75,
       teeth=12,
       pressure_angle=20.0,
       face_width=10.0
   )
   registry.register(planet_interface)

   # Set interface registry on detector
   detector.set_interface_registry(registry)

   # Now gear mesh overlaps won't be flagged as collisions
   results = detector.check_assembly(transforms)

Planetary Gearbox Helper
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Automatically create all interfaces for a planetary gearbox
   interfaces = create_planetary_gearbox_interfaces(
       prefix="axis1",
       sun_part="AXIS1_SUN_GEAR",
       planet_parts=["AXIS1_PLANET_1", "AXIS1_PLANET_2", "AXIS1_PLANET_3"],
       ring_part="AXIS1_RING_GEAR",
       module=0.75,
       sun_teeth=18,
       planet_teeth=12,
       ring_teeth=42,
       face_width=10.0
   )

   for interface in interfaces:
       registry.register(interface)

Thread Interface
~~~~~~~~~~~~~~~~

.. code-block:: python

   from yapcad.collision import ThreadInterface, ThreadType, ThreadClass

   # Define screw thread
   screw_thread = ThreadInterface(
       name="m3_screw_threads",
       part_name="M3_SCREW",
       thread_type=ThreadType.METRIC,
       nominal_diameter=3.0,      # M3
       pitch=0.5,                 # 0.5mm pitch
       thread_class=ThreadClass.MEDIUM,
       engagement_length=10.0,    # 10mm threaded depth
       is_external=True
   )
   registry.register(screw_thread)

   # Define matching nut threads
   nut_thread = ThreadInterface(
       name="m3_nut_threads",
       part_name="M3_NUT",
       thread_type=ThreadType.METRIC,
       nominal_diameter=3.0,
       pitch=0.5,
       thread_class=ThreadClass.MEDIUM,
       engagement_length=10.0,
       is_external=False
   )
   registry.register(nut_thread)

Visualization
=============

The ``yapcad.viewer`` package provides a VTK-based multi-viewport viewer with REST API control.

VTKViewer
---------

.. code-block:: python

   from yapcad.viewer import VTKViewer, ViewerConfig

   # Create viewer configuration
   config = ViewerConfig(
       stl_dir="/path/to/stl/files",
       positions_file="/path/to/positions.json",
       window_size=(1600, 1200),
       background_color=(0.2, 0.2, 0.2)
   )

   # Create viewer
   viewer = VTKViewer(config)

   # Load parts from JSON
   viewer.load_from_json()

   # Start interactive window
   viewer.start()

Programmatic Control
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Load specific parts
   viewer.load_parts({
       "SERVO": "xh540_servo.stl",
       "BRACKET": "bracket.stl",
       "LINK1": "link1.stl"
   })

   # Set part transforms
   viewer.set_part_transform("SERVO", servo_transform)
   viewer.set_part_transform("BRACKET", bracket_transform)

   # Toggle X-ray mode (transparency)
   viewer.set_xray_mode(True)

   # Focus on specific parts
   viewer.focus_parts(["SERVO", "BRACKET"])

   # Capture screenshot
   viewer.screenshot("assembly_view.png")

   # Reload transforms from JSON
   viewer.reload_positions()

ViewerAPIServer
---------------

The API server provides REST and WebSocket control for remote operation:

.. code-block:: python

   from yapcad.viewer import ViewerAPIServer

   # Create server (viewer runs in background thread)
   server = ViewerAPIServer(viewer, host="0.0.0.0", port=5000)

   # Start server (blocks)
   server.run()

REST API Endpoints
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Load parts
   curl -X POST http://localhost:5000/api/load \
        -H "Content-Type: application/json" \
        -d '{"parts": {"SERVO": "servo.stl"}}'

   # Set transform
   curl -X POST http://localhost:5000/api/transform/SERVO \
        -H "Content-Type: application/json" \
        -d '{"matrix": [[1,0,0,10],[0,1,0,20],[0,0,1,30],[0,0,0,1]]}'

   # Toggle X-ray mode
   curl -X POST http://localhost:5000/api/xray \
        -d '{"enabled": true}'

   # Capture screenshot
   curl -X GET http://localhost:5000/api/screenshot/output.png

   # Reload positions
   curl -X POST http://localhost:5000/api/reload

WebSocket Events
~~~~~~~~~~~~~~~~

The server emits real-time events via WebSocket:

.. code-block:: javascript

   const socket = io('http://localhost:5000');

   socket.on('viewer_event', (data) => {
       console.log('Event:', data.type);
       console.log('Data:', data.data);
   });

   // Events: 'parts_loaded', 'transform_updated',
   //         'xray_toggled', 'screenshot_saved'

File-Based Command Interface
-----------------------------

For backward compatibility, the viewer supports file-based commands:

.. code-block:: bash

   # Command file: output/viewer_cmd.txt
   echo "load full" > output/viewer_cmd.txt
   echo "xray" > output/viewer_cmd.txt
   echo "screenshot check.png" > output/viewer_cmd.txt
   echo "reload" > output/viewer_cmd.txt

   # Viewer polls this file and executes commands

Multi-Viewport Layout
---------------------

The viewer provides simultaneous orthographic views:

* **ISO**: Isometric view (primary viewport)
* **TOP**: Top-down view (XY plane)
* **FRONT**: Front view (XZ plane)
* **SIDE**: Side view (YZ plane)

All viewports update in real-time when transforms change.

Integration Workflow
====================

Typical Assembly Workflow
--------------------------

1. **Define Parts with Datums**

   .. code-block:: python

      # Create part definitions
      servo = PartDefinition("SERVO", geometry_source="servo.step")
      servo.add_datum(Datum("stator_face", DatumType.PLANE, ...))
      servo.add_datum(Datum("output_shaft", DatumType.AXIS, ...))

      bracket = PartDefinition("BRACKET", geometry_source="bracket.step")
      bracket.add_datum(Datum("mount_face", DatumType.PLANE, ...))
      bracket.add_datum(Datum("bore_axis", DatumType.AXIS, ...))

2. **Register Datums**

   .. code-block:: python

      from yapcad.assembly.datum_registry import DatumRegistry

      DatumRegistry.register_source("SERVO", servo.datums)
      DatumRegistry.register_source("BRACKET", bracket.datums)

3. **Define Mate Constraints**

   .. code-block:: python

      mates = [
          Mate("mount", MateType.COINCIDENT,
               "BRACKET", "mount_face", "SERVO", "stator_face"),
          Mate("alignment", MateType.CONCENTRIC,
               "BRACKET", "bore_axis", "SERVO", "output_shaft")
      ]

4. **Solve Constraints**

   .. code-block:: python

      solver = MateConstraintSolver()
      results = solver.solve_all(mates)

      servo_transform = results["mount"].transform

5. **Build Kinematic Chain**

   .. code-block:: python

      chain = KinematicChain("assembly")
      chain.add_part(KinematicPart("BRACKET", parent=None,
                                   joint=Joint("base", JointType.FIXED)))
      chain.add_part(KinematicPart("SERVO", parent="BRACKET",
                                   joint=Joint("servo_joint", JointType.FIXED,
                                             base_transform=servo_transform)))

6. **Register Interface Volumes**

   .. code-block:: python

      registry = InterfaceRegistry()
      registry.register(GearMeshInterface("sun_teeth", "SUN_GEAR", ...))
      registry.register(GearMeshInterface("planet_teeth", "PLANET", ...))

7. **Validate Collisions**

   .. code-block:: python

      provider = FileBasedProvider(Path("output/assembly"))
      detector = CollisionDetector(provider)
      detector.set_interface_registry(registry)

      world_transforms = {
          "BRACKET": chain.get_world_transform("BRACKET"),
          "SERVO": chain.get_world_transform("SERVO")
      }

      results = detector.check_assembly(world_transforms)
      collisions = [r for r in results if r.is_error]

      if collisions:
          print(f"ERROR: {len(collisions)} collisions found!")
      else:
          print("Assembly is collision-free!")

8. **Visualize**

   .. code-block:: python

      # Export transforms to JSON
      chain.export_json("output/positions.json")

      # Start viewer
      config = ViewerConfig(
          stl_dir="output/assembly",
          positions_file="output/positions.json"
      )
      viewer = VTKViewer(config)
      viewer.load_from_json()
      viewer.start()

Best Practices
--------------

**Never Hardcode Offsets**

.. code-block:: python

   # WRONG: Hardcoded transform
   servo_tf = Transform.from_translation(0, 0, 28.05)  # Magic number!

   # RIGHT: Solve from constraints
   result = solver.solve_mate(mount_mate)
   servo_tf = result.transform

**Always Use Datums**

.. code-block:: python

   # WRONG: Manual transform calculation
   # "The servo is 28mm tall, so offset by that amount..."

   # RIGHT: Reference named datum features
   Mate("mount", MateType.COINCIDENT,
        "BRACKET", "mount_face",
        "SERVO", "stator_face")

**Validate Constraints**

.. code-block:: python

   # Check constraint satisfaction
   result = solver.solve_mate(mate)
   validation = solver.validate_transform(result.transform)

   if not validation.valid:
       print("Transform validation failed!")
       for msg in validation.error_messages:
           print(f"  {msg}")

**Use Interface Volumes for Designed Overlaps**

.. code-block:: python

   # WRONG: Ignore gear mesh collisions manually

   # RIGHT: Register gear mesh interfaces
   registry.register(GearMeshInterface("sun_teeth", ...))
   registry.register(GearMeshInterface("planet_teeth", ...))
   detector.set_interface_registry(registry)

**Prefer STEP over STL**

.. code-block:: python

   # STEP files preserve exact BREP geometry
   part = PartDefinition("GEAR", geometry_source="gear.step")

   # STL is tessellated approximation
   # Use only for final rendering/visualization

API Reference
=============

For detailed API documentation, see:

* :mod:`yapcad.assembly` - Datum features, mates, and constraint solver
* :mod:`yapcad.assembly.solver` - MateConstraintSolver and solving functions
* :mod:`yapcad.kinematics` - Transform, Joint, KinematicChain
* :mod:`yapcad.collision` - CollisionDetector, interface volumes
* :mod:`yapcad.viewer` - VTKViewer, ViewerAPIServer

Examples
========

Complete examples can be found in the yapCAD repository. For a production assembly system implementation, see the reference design (not included in this documentation to maintain generality).

Key example files:

* ``examples/assembly_demo.py`` - Basic assembly with mates
* ``examples/kinematic_chain_demo.py`` - Robot arm with joints
* ``examples/collision_detection_demo.py`` - Multi-part collision validation
* ``examples/viewer_api_demo.py`` - Remote viewer control via REST API

Contributing
============

The assembly system is actively developed and welcomes contributions. Areas for enhancement:

* Additional mate types (CAM, SLOT, etc.)
* Motion simulation and dynamics
* URDF/SDF export for robotics simulators
* Improved collision algorithms
* Performance optimization for large assemblies

See the yapCAD contribution guidelines for details.

License
=======

The assembly system is part of yapCAD and is released under the MIT License.

Copyright (c) 2026 yapCAD contributors

Assembly system contributed by Jeremy Mika.
