.. _rom-guide:

=================================
Range-of-Motion Sweep (ROM) Guide
=================================

The ``yapcad.rom`` package extends the DSL with a **range-of-motion (ROM)
sweep** system that can detect self-collisions and tight clearances in
kinematic assemblies — all without generating meshes.  Clearance checks use
**SDF field evaluation** on sampled grids, so a multi-DOF sweep completes
in seconds rather than minutes.

Core concepts
=============

1. **Kinematic annotations** — Annotate DSL commands with
   ``@meta(kinematic.*)`` to describe joint type, parent, axis, and motion
   limits.
2. **RomTree** — Builds a ``KinematicChain`` from the annotations and
   exposes adjacent part pairs.
3. **Clearance evaluation** — Samples a grid on one link's bounding box,
   transforms it into the adjacent link's local frame, and evaluates the
   SDF field.  No meshing needed.
4. **ROM sweep** — Iterates over a product space of DOF values, computing
   pairwise clearance at each configuration.

Getting started
===============

Annotate a DSL module
---------------------

.. code-block:: dsl

    module phonebot_leg

    @meta(kinematic.joint_kind="revolute")
    @meta(kinematic.parent="PELVIS")
    @meta(kinematic.axis=[0, 0, 1])
    @meta(kinematic.limits=[-90, 90])
    @meta(kinematic.zero_pose=0.0)
    command HIP_YAW_JOINT() -> solid:
        srv = translate(servo.XL430_ENVELOPE(), 0, 40, 0)
        emit srv

    @meta(kinematic.joint_kind="revolute")
    @meta(kinematic.parent="HIP_YAW_JOINT")
    @meta(kinematic.axis=[1, 0, 0])
    @meta(kinematic.limits=[-45, 45])
    command HIP_ROLL_JOINT() -> solid:
        ...

The supported ``kinematic.*`` fields are:

.. list-table::
   :header-rows: 1

   * - Field
     - Type
     - Description
   * - ``joint_kind``
     - ``"fixed"`` | ``"revolute"`` | ``"prismatic"``
     - Type of joint
   * - ``parent``
     - string
     - Name of the parent command/link
   * - ``parent_frame``
     - string
     - Optional; name of the frame on the parent to attach to
       (default ``"ORIGIN"``)
   * - ``axis``
     - ``[x, y, z]``
     - Motion axis for revolute/prismatic joints
   * - ``limits``
     - ``[min, max]``
     - Motion range in degrees (revolute) or mm (prismatic)
   * - ``zero_pose``
     - float
     - Default joint value
   * - ``printable``
     - bool
     - Whether this is a 3D-printable part
   * - ``description``
     - string
     - Human-readable label

Building a tree and running a sweep
------------------------------------

From Python::

    from yapcad.rom import RomTree, sweep_kinematics, SweepConfig
    from yapcad.dsl import compile_and_run

    # Each command can be compiled independently
    res_hip = compile_and_run(source, "HIP_YAW_JOINT", {})

    # Build a tree from explicit RomLink objects
    tree = RomTree("phonebot_leg")
    tree.add_link(RomLink(
        name="HIP_YAW", parent="PELVIS",
        joint_kind="revolute", axis=(0, 0, 1),
        limits=(-90, 90), zero_pose=0.0,
    ))

    # Run a single-DOF sweep
    result = sweep_kinematics(
        tree,
        {"HIP_YAW": (-90, 91, 5)},
        config=SweepConfig(sdf_resolution=24, interference_threshold=0.0),
    )

    print(result.summary())

The sweep runner produces:

- ``result.events`` — list of :class:`InterferenceEvent` where clearance
  was below threshold
- ``result.near_misses`` — configurations within the clearance threshold
- ``result.pair_trace`` — full per-pair clearance trace
- ``result.summary()`` — human-readable report

Multi-DOF sweeps
----------------

Pass multiple schedules to sweep the product space::

    result = sweep_kinematics(
        tree,
        {"HIP_YAW": (-45, 46, 5)},
        {"KNEE_PITCH": (0, 136, 5)},
        fixed_angles={"HIP_ROLL": 0.0},
    )

Use small step counts for dimensions above 3 to avoid combinatorial
explosion (25 × 25 = 625 configs is fine; 25 × 25 × 25 = 15,625 is still
manageable; 25 × 25 × 25 × 25 = 390,625 is not).

Reading results
===============

::

    # Basic summary
    print(result.summary())

    # Worst clearance violation
    ev = result.worst_interference()
    print(f"{ev.pair[0]} vs {ev.pair[1]} at {ev.config}: {ev.min_distance:.2f} mm")

    # Full trace for one pair
    trace = result.clearance_at("PELVIS", "HIP_YAW_JOINT")
    for config_key, min_dist, mean_dist in trace[:5]:
        print(f"  {config_key}: min={min_dist:.2f} mm")

    # Detailed report
    from yapcad.rom import interference_summary
    print(interference_summary(result))

Clearance mechanics
===================

``clearance_sdf`` evaluates the signed distance of *link B*'s surface
against *link A*'s field:

1. Sample a grid in B's bounding box (padded by ``sdf_padding`` mm).
2. Keep only surface-adjacent points (where ``|SDF| < 10``).
3. Transform those points from B's local frame → world → A's local frame.
4. Evaluate A's SDF at the transformed points.

A negative value means B penetrates A.  The minimum across all points is
the bottleneck clearance.  The entire operation is pure arithmetic on
numpy arrays — no mesh generation, no OCC, no boolean engines.

::

    from yapcad.rom import clearance_sdf
    from yapcad.sdf import box
    from yapcad.kinematics import Transform

    node_a = box((20, 20, 20))
    node_b = box((10, 10, 10))
    tf_a = Transform.identity()
    tf_b = Transform.from_translation(15.0, 0.0, 0.0)

    clr = clearance_sdf(node_a, node_b, tf_a, tf_b, resolution=24)
    print(f"min_distance: {clr['min_distance']:.2f} mm")

Reference
=========

.. autoclass:: yapcad.rom.RomTree
.. autoclass:: yapcad.rom.RomLink
.. autoclass:: yapcad.rom.SweepConfig
.. autoclass:: yapcad.rom.SweepResult
.. autoclass:: yapcad.rom.InterferenceEvent

.. autofunction:: yapcad.rom.sweep_kinematics
.. autofunction:: yapcad.rom.clearance_sdf
.. autofunction:: yapcad.rom.dof_schedule_from_tree
.. autofunction:: yapcad.rom.interference_summary
