.. _rom-guide:

=================================
Range-of-Motion Sweep (ROM) Guide
=================================

A range-of-motion sweep drives an assembly's joints through their declared
limits and measures, at every pose, how close its parts come to each
other.  The measurement is made on each part's **signed distance field**,
so nothing is meshed, no boolean kernel runs, and the only accuracy knob
is the sampling pitch — which every result reports alongside its number.

``yapcad.rom`` adds no kinematics of its own.  The joints, their limits
and their placement all come from the DSL assembly, so a sweep cannot
disagree with the design it came from.

What the design has to say
==========================

Everything a sweep needs is already expressible in the DSL (see
:doc:`dsl_reference` → *Assembly Solving and Posing*):

``@meta(assembly.datums=[...])``
    Places the joint.  A revolute mate turns about the datum axis, so the
    datum's ``origin_mm`` and ``direction`` *are* the joint.

``add_named_mate(a, name, "revolute", parent, datum, child, datum)``
    Names the degree of freedom.  That name is what a sweep drives.

``set_mate_limits(a, name, lo, hi)``
    Bounds it, in radians.  A sweep refuses to leave the declared range
    rather than guessing one.

``solve_assembly(a, root)`` and ``emit assembly_compound(a)``
    Solve the placement tree and retain the assembly on the result, which
    is the same object the ``.ycpkg`` builder reads.

Here is a complete two-joint design: an arm on a pivot, and a stop block
bolted to the base beside it.

.. code-block:: python

    module rom_arm

    @meta(assembly.datums=[
        {"id": "pivot", "kind": "axis",
         "origin_mm": [0.0, 0.0, 10.0], "direction": [0.0, 0.0, 1.0]},
        {"id": "stop_mount", "kind": "axis",
         "origin_mm": [0.0, 25.0, 15.0], "direction": [0.0, 0.0, 1.0]}
    ])
    command BASE() -> solid:
        emit box(40.0, 40.0, 20.0)

    @meta(assembly.datums=[
        {"id": "hub", "kind": "axis",
         "origin_mm": [0.0, 0.0, 0.0], "direction": [0.0, 0.0, 1.0]}
    ])
    command ARM() -> solid:
        emit translate(box(60.0, 10.0, 10.0), 30.0, 0.0, 15.0)

    @meta(assembly.datums=[
        {"id": "mount", "kind": "axis",
         "origin_mm": [0.0, 0.0, 0.0], "direction": [0.0, 0.0, 1.0]}
    ])
    command STOP() -> solid:
        emit box(10.0, 20.0, 30.0)

    command BUILD() -> solid:
        let a: assembly = assembly("rom_arm")
        add_part(a, BASE(), "base")
        add_part(a, ARM(), "arm")
        add_part(a, STOP(), "stop")
        add_named_mate(a, "elbow", "revolute",
                       "base", "pivot", "arm", "hub")
        add_named_mate(a, "stop_mount", "rigid",
                       "base", "stop_mount", "stop", "mount")
        set_mate_limits(a, "elbow", radians(-90.0), radians(90.0))
        solve_assembly(a, "base")
        emit assembly_compound(a)

Build it as fields and sweep it
===============================

The parts have to be **authored as fields**, which means building with
``representation="sdf"``.  A sweep will not substitute a bounding box for
a part it cannot measure.

.. code-block:: python

    from yapcad.dsl import compile_and_run
    from yapcad.rom import RomModel, radian_range, sweep_rom, SweepConfig

    result = compile_and_run(source, "BUILD", {},
                             representation="sdf", source_path=path)
    model = RomModel.from_result(result)
    print(model.describe())

::

    Assembly 'rom_arm' rooted at 'base'
      parts: 3
      joints: 2 (1 movable)
        elbow: revolute, -90.0 deg to 90.0 deg (base - arm)
        stop_mount: rigid (base - stop)

``describe()`` is worth printing before every sweep: it names the joints a
sweep can drive, their limits, and any part that has no field.

.. code-block:: python

    sweep = sweep_rom(model, {"elbow": radian_range(0.0, 90.0, 15.0)},
                      config=SweepConfig(pitch_mm=0.5))
    print(sweep.summary())

::

    Sweep of 'rom_arm': 7 configuration(s), 2 pair(s)
      elbow: 0.0 deg to 90.0 deg in 7 step(s)
      measurements: 6 sampled, 8 resolved from bounding boxes
      closest approach: elbow=90.0 deg | arm - stop: -5.000 mm (+/- 0.433, field)
      interference (< 0.0 mm): 3
      near misses (< 5.0 mm): 1
      worst interferences:
        elbow=90.0 deg | arm - stop: -5.000 mm (+/- 0.433, field)
        elbow=75.0 deg | arm - stop: -4.981 mm (+/- 0.433, field)
        elbow=60.0 deg | arm - stop: -1.830 mm (+/- 0.433, field)

Omit the schedule and every movable joint is swept across its limits at
5 degrees::

    sweep = sweep_rom(model)

Units
=====

Joint values are **radians** everywhere in the API, the same unit as
``set_mate_limits`` and ``set_joint_position``.  Use
:func:`~yapcad.rom.sweep.radian_range` or
:func:`~yapcad.rom.sweep.schedule_from_model` to build a schedule from the
degrees you are thinking in; every report prints degrees.

Distances are millimetres.

Which pairs get measured
========================

By default, a sweep measures **every pair of parts the swept joints move
relative to one another** — including pairs that share no mate, which is
usually where interference is.  In the example above that is
``base - arm`` and ``arm - stop``: the arm swings towards a stop block it
is not jointed to at all.  ``base - stop`` is left out, because a rigid
mate holds those two still however the assembly is posed, so measuring it
at every configuration would report the same design fact seven times (and
rigidly mated parts routinely interpenetrate on purpose — a boss sunk
into a plate, a screw in its hole).

Ask for something else when you need it::

    sweep_rom(model, schedule, pairs="all")      # every pair
    sweep_rom(model, schedule, pairs="mates")    # pairs sharing a joint
    sweep_rom(model, schedule, pairs=[("arm", "stop")])
    sweep_rom(model, schedule, extra_pairs=[("base", "stop")])

Hold a joint still across the sweep with ``fixed={"other": 0.0}``.  A
joint that appears in both the schedule and ``fixed`` is an error, not a
silent precedence rule.

What the numbers mean
=====================

:func:`~yapcad.rom.clearance.clearance_sdf` samples the second part on a
regular grid, keeps the points at or inside its surface, transforms them
into the first part's frame, and evaluates that part's field:

.. code-block:: text

    for each grid point p at or inside part B
        d = SDF_A( A_world^-1 @ B_world @ p )
    clearance = min d

A negative result means the parts overlap by that much.  Samples are kept
within half a cell of B's surface -- inside it, or just outside by less
than that -- so the error is bounded by half a cell diagonal, which every
result carries as ``tolerance_mm``::

    elbow=45.0 deg | arm - stop: +2.071 mm (+/- 0.433, field)

So choose ``pitch_mm`` against the clearance you care about.  At the
default 1 mm pitch a 0.5 mm interference is resolved; chasing a 0.1 mm
one needs a finer pitch, and ``Clearance.resolved`` is False whenever a
result is smaller than its own tolerance.  The summary counts those
measurements and says to refine the pitch.

There is a second source of error, in the fields themselves.  yapCAD marks
a node ``exact`` when its magnitude is the true distance, and hard CSG does
not preserve that: a part built with ``union``, ``intersect`` or
``subtract`` — which is most parts — carries a Lipschitz-1 *lower bound*
instead.  Where either field is inexact the result is flagged
``conservative`` and prints as ``at least``::

    elbow=45.0 deg | arm - chassis: at least +1.740 mm (+/- 0.433, field)

Read that as: the parts are at least this far apart.  The magnitude is
closer to zero than the truth — a gap is never overstated, a penetration
may be understated — while the sign is right either way, so a sweep never
misses an interference.  Only when both fields are exact (single
primitives, rounded primitives, transforms of them) is the true distance
within ``tolerance_mm`` of the number.  A ``"bounds"`` result is
conservative by construction.

Nothing else biases the number.  (An earlier version of this code grew
the sampled region by a 2 mm margin and kept the points in it, so every
clearance came back 2 mm short at every resolution, and parts 1 mm apart
were reported as interfering.  The test suite now checks measured
distances against analytic ones for exactly this reason.)

Of the two parts, the **smaller** is the one sampled and the larger is the
one whose field is evaluated.  The measurement is the same either way to
within its tolerance, but a 20 mm part beside a 400 mm one needs four
orders of magnitude fewer sample points, so the pitch you ask for is the
pitch you get.

Pairs whose **bounding boxes** are further apart than
``clearance_threshold`` are resolved from the boxes alone — the box gap is
a valid lower bound on the distance between the parts — and reported with
``method="bounds"`` rather than sampled.  That is what makes
``pairs="all"`` affordable on a large assembly, and it is why a pair
reported as ``10.000 mm (bounds)`` means "at least 10 mm", while
``+2.071 mm (field)`` means "2.071 mm, give or take the tolerance".

A sweep cannot look clean when it measured nothing
==================================================

Every way a measurement can fail to happen is either an error or
reported:

- a part with no field (anything not built with ``representation="sdf"``,
  an imported mesh) makes the pair an **error**, naming the part and what
  to do about it;
- ``SweepConfig(skip_unmeasurable=True)`` continues instead, but the pair
  and its reason appear in ``result.skipped``, in every ``summary()`` and
  in every report;
- a sweep that measured nothing at all says so in capitals;
- an unknown joint or part name, a value outside the declared limits, an
  empty schedule, or a joint both swept and held are all errors — joint
  values are never clamped, and names are never quietly ignored;
- ``early_stop`` marks the result incomplete.

Reading results
===============

.. code-block:: python

    sweep.worst()            # the closest approach found, or None
    sweep.events             # samples below interference_threshold
    sweep.near_misses        # samples below clearance_threshold
    sweep.pair_minima()      # worst configuration per pair
    sweep.trace("arm", "stop")   # every measurement of one pair
    sweep.skipped            # pairs that were not measured, and why

``events`` and ``near_misses`` are views over ``samples``, so thresholds
can be reconsidered without sweeping again.  A trace reads as the motion
does::

    elbow=0.0 deg  | arm - stop: +10.000 mm (+/- 0.000, bounds)
    elbow=15.0 deg | arm - stop:  +8.195 mm (+/- 0.433, field)
    elbow=30.0 deg | arm - stop:  +5.490 mm (+/- 0.433, field)
    elbow=45.0 deg | arm - stop:  +2.071 mm (+/- 0.433, field)
    elbow=60.0 deg | arm - stop:  -1.830 mm (+/- 0.433, field)
    elbow=75.0 deg | arm - stop:  -4.981 mm (+/- 0.433, field)
    elbow=90.0 deg | arm - stop:  -5.000 mm (+/- 0.433, field)

:func:`~yapcad.rom.report.interference_summary` prints the long form: the
model it swept, the settings it used, the worst configuration for every
pair, each interference with its sample statistics, and whatever could
not be measured.

Cost
====

One measurement is one grid of field evaluations, so cost scales with
``(part size / pitch)^3`` for the sampled part, bounded by
``SweepConfig.max_samples`` (which coarsens the pitch rather than
overrunning, and reports the pitch it used).  The bounds test dismisses
distant pairs for the price of eight matrix-vector products, so the
measured set is normally much smaller than the selected set — six of the
fourteen pair-configurations above.

Relationship to the package format
==================================

A ``.ycpkg`` built from the same command records the mate graph with its
limits (``yapcad-assembly-v0.1``: parts, transforms, datums, mates, joint
values and couplings) and, at geometry JSON v0.3, each component's SDF
construction.  Sweeps today run from the DSL result, which is the same
assembly object the package is built from; driving one from a package on
disk is a natural extension and not yet implemented.

Reference
=========

.. autoclass:: yapcad.rom.RomModel
   :members:
.. autoclass:: yapcad.rom.RomJoint
   :members:
.. autoclass:: yapcad.rom.RomLink
   :members:
.. autoclass:: yapcad.rom.Clearance
   :members:
.. autoclass:: yapcad.rom.SweepConfig
.. autoclass:: yapcad.rom.SweepResult
   :members:
.. autoclass:: yapcad.rom.PairSample
   :members:

.. autofunction:: yapcad.rom.sweep_rom
.. autofunction:: yapcad.rom.schedule_from_model
.. autofunction:: yapcad.rom.radian_range
.. autofunction:: yapcad.rom.clearance_sdf
.. autofunction:: yapcad.rom.bounds_clearance
.. autofunction:: yapcad.rom.world_bounds
.. autofunction:: yapcad.rom.interference_summary
.. autofunction:: yapcad.rom.assembly_from_result
