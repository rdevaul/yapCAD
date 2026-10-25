"""yapcad.rom -- range-of-motion sweeps on signed distance fields.

A range-of-motion sweep drives an assembly's joints through their declared
limits and measures, at every pose, how close the parts come to each
other. The measurement is made on the parts' own SDF fields, so nothing is
meshed and no boolean kernel is involved.

The kinematics come from the design, not from this package. A DSL assembly
already carries everything a sweep needs:

- ``@meta(assembly.datums=[...])`` places the joint on each part;
- ``add_named_mate(a, name, "revolute", parent, datum, child, datum)``
  names the degree of freedom;
- ``set_mate_limits(a, name, lo, hi)`` bounds it (radians);
- ``solve_assembly(a, root)`` and ``emit assembly_compound(a)`` retain the
  solved assembly on the result.

Build such a design as fields and sweep it::

    from yapcad.dsl import compile_and_run
    from yapcad.rom import RomModel, schedule_from_model, sweep_rom

    result = compile_and_run(source, "BUILD", {},
                             representation="sdf", source_path=path)
    model = RomModel.from_result(result)
    print(model.describe())

    sweep = sweep_rom(model, schedule_from_model(model, step_deg=5.0))
    print(sweep.summary())

    worst = sweep.worst()
    print(worst.describe())

Joint values are radians, as in the DSL. Reports print degrees.

The parts must be authored as fields, which means building with
``representation="sdf"``; a sweep refuses to guess at a part it cannot
measure.
"""

from .clearance import (
    Clearance,
    bounds_clearance,
    clearance_sdf,
    world_bounds,
)
from .model import (
    RomError,
    RomJoint,
    RomLink,
    RomModel,
    assembly_from_result,
)
from .report import interference_summary
from .sweep import (
    PairSample,
    SweepConfig,
    SweepResult,
    radian_range,
    schedule_from_model,
    sweep_rom,
)

__all__ = [
    # model
    "RomModel",
    "RomJoint",
    "RomLink",
    "RomError",
    "assembly_from_result",
    # clearance
    "Clearance",
    "clearance_sdf",
    "bounds_clearance",
    "world_bounds",
    # sweep
    "sweep_rom",
    "schedule_from_model",
    "radian_range",
    "SweepConfig",
    "SweepResult",
    "PairSample",
    # reporting
    "interference_summary",
]
