================================
yapCAD Validation Test Schema
================================

**Version:** ``validation-schema-v0.1``
**Status:** Draft

This document defines the schema for validation tests that can be stored in
yapCAD packages (``.ycpkg``). Validation tests allow automated verification
of design requirements against geometry, simulation results, or computed
properties.

Overview
========

Validation tests are stored as YAML files under ``validation/plans/`` in a
package. Each test file defines:

1. **Identity**: Unique ID, kind, and human-readable description
2. **Inputs**: What geometry/data the test operates on
3. **Acceptance Criteria**: Pass/fail conditions
4. **Execution**: How to run the test (local/remote, tool configuration)
5. **Results**: Captured outputs stored under ``validation/results/``

Test Kinds
==========

Validation tests are categorized by ``kind``:

Geometric Tests (``geometric``)
-------------------------------

Tests that verify geometric properties of shapes without simulation.

- **volume**: Check solid volume against limits
- **area**: Check surface area or 2D region area
- **bbox**: Check bounding box dimensions
- **distance**: Check distances between features
- **clearance**: Check minimum clearance between parts

Measurement Tests (``measurement``)
-----------------------------------

Tests that verify computed measurements or derived values.

- **dimension**: Check a specific dimension (length, width, height)
- **mass**: Check mass given density
- **centroid**: Check center of mass location
- **moment**: Check moment of inertia

Simulation Tests (``structural``, ``thermal``, ``cfd``, ``multiphysics``)
-------------------------------------------------------------------------

Tests that invoke external solvers for physical simulation.

- **structural**: FEA stress/deflection analysis (backends: fenics, calculix, abaqus)
- **thermal**: Heat transfer analysis
- **cfd**: Computational fluid dynamics
- **multiphysics**: Coupled physics (backend: comsol)

Assembly Tests (``assembly``)
-----------------------------

Tests that verify assembly constraints and relationships.

- **interference**: Check for part overlaps
- **fit**: Check clearance fit requirements
- **alignment**: Check alignment constraints

Schema Definition
=================

Common Fields
-------------

All validation tests share these fields:

.. code-block:: yaml

   # Required
   id: string              # Unique identifier (e.g., "volume-check-001")
   kind: string            # Test category (see above)
   backend: string         # Tool/method for evaluation

   # Optional but recommended
   name: string            # Human-readable name
   description: string     # Detailed description

   # Input specification
   geometry:
     source: path          # Path to geometry file (relative to package)
     entities: list<str>   # Entity IDs to test (empty = all)

   # Pass/fail criteria
   acceptance:
     <metric_name>:
       limit: number       # Threshold value
       comparison: string  # "<=", ">=", "<", ">", "==", "~=" (approximate)
       tolerance: number   # For "~=" comparison (optional)
       units: string       # Units for documentation
       description: string # Human-readable explanation

   # Execution context
   execution:
     mode: string          # "local" or "remote"
     command: string       # Solver command (optional)
     env: map<str,str>     # Environment variables
     timeout_s: number     # Maximum execution time

   # Metadata
   metadata:
     author: string
     created: string       # ISO 8601 date
     revision: number
     tags: list<str>


Geometric Test Schema
---------------------

For ``kind: geometric``:

.. code-block:: yaml

   id: volume-check
   kind: geometric
   backend: yapcad          # Native yapCAD computation

   geometry:
     source: geometry/primary.json
     entities: ["solid-main"]

   check:
     property: volume       # volume, area, bbox, centroid
     units: mm3

   acceptance:
     volume:
       limit: 50000.0
       comparison: ">="
       description: "Minimum material volume for strength"

   # For bounding box checks
   # check:
   #   property: bbox
   #   axis: x              # x, y, z, or "all"
   #
   # acceptance:
   #   bbox.x:
   #     limit: 100.0
   #     comparison: "<="


Measurement Test Schema
-----------------------

For ``kind: measurement``:

.. code-block:: yaml

   id: mass-check
   kind: measurement
   backend: yapcad

   geometry:
     source: geometry/primary.json
     entities: ["solid-main"]

   check:
     property: mass
     density_kgm3: 2700     # Material density

   acceptance:
     mass_kg:
       limit: 0.5
       comparison: "<="
       description: "Maximum mass for weight budget"


Structural Analysis Schema
--------------------------

For ``kind: structural`` (FEA):

.. code-block:: yaml

   id: thrust-fea
   kind: structural
   backend: fenics          # or calculix, abaqus

   name: "Thrust Plate Static Analysis"
   description: |
     Linear elastic analysis under max thrust.

   geometry:
     source: geometry/primary.json
     entities: ["solid-main"]

   materials:
     default:
       name: "6061-T6 Aluminum"
       youngs_modulus_pa: 68.9e9
       poisson_ratio: 0.33
       density_kgm3: 2700

   loads:
     - id: thrust_pressure
       type: pressure           # pressure, force, moment
       surfaces: ["motor_mount"]
       magnitude_pa: 2.61e6

   boundaryConditions:
     - id: fixed_base
       type: fixed              # fixed, pinned, roller, spring
       surfaces: ["bottom_face"]

   acceptance:
     displacement.max_mm:
       limit: 10.0
       comparison: "<="
     stress.von_mises.max_pa:
       limit: 180e6
       comparison: "<="

   backendOptions:
     mesh:
       element_size: 5.0
       element_order: 1
     solver:
       type: linear

   faceNaming:
     strategy: geometric
     rules:
       - name: motor_mount
         selector:
           type: cylindrical
           radius: 50.8

   execution:
     mode: local
     timeout_s: 600


Assembly Test Schema
--------------------

For ``kind: assembly``:

.. code-block:: yaml

   id: clearance-check
   kind: assembly
   backend: yapcad

   geometry:
     source: geometry/primary.json
     entities: ["part_a", "part_b"]

   check:
     property: clearance
     between:
       - entity: part_a
         surface: outer
       - entity: part_b
         surface: inner

   acceptance:
     min_clearance_mm:
       limit: 0.5
       comparison: ">="
       description: "Sliding fit clearance"


Surface Selection Strategies
============================

For FEA boundary conditions and loads, surfaces can be selected by name or
by geometric strategy. The ``strategy`` field invokes backend-specific
surface selectors:

Built-in Strategies
-------------------

Position-based (all backends):

- ``z_min`` - Bottom surface (minimum Z coordinate)
- ``z_max`` - Top surface (maximum Z coordinate)
- ``x_min``, ``x_max``, ``y_min``, ``y_max`` - Axis-aligned selections

Shape-based (fenics backend):

- ``z_max_ring`` - Top surface excluding central region (for ring/annular loads)
- ``motor_mount`` - Central cylindrical region
- ``stringer_notches`` - Edge notch regions at specified angles

Face Naming Rules
-----------------

The ``faceNaming`` section can define custom surface selection rules:

.. code-block:: yaml

   faceNaming:
     strategy: geometric
     rules:
       - name: motor_mount_inner
         selector:
           type: cylindrical
           axis: [0, 0, 1]
           radius: 50.8
           tolerance: 1.0

       - name: stringer_notch_1
         selector:
           type: planar
           normal: [1, 0, 0]
           position_constraint:
             axis: x
             at_max: true


Result Schema
=============

After execution, results are stored in ``validation/results/<plan_id>/``:

.. code-block:: json

   {
     "plan_id": "thrust-fea",
     "status": "passed",
     "timestamp": "2025-12-30T15:30:00Z",
     "backend": "fenics",
     "backend_version": "0.8.0",
     "execution_time_s": 45.2,

     "metrics": {
       "displacement.max_mm": 7.3,
       "stress.von_mises.max_pa": 145e6,
       "nodes": 12500,
       "elements": 8200
     },

     "acceptance_results": {
       "displacement.max_mm": {
         "value": 7.3,
         "limit": 10.0,
         "comparison": "<=",
         "passed": true,
         "margin": 2.7,
         "margin_pct": 27.0
       },
       "stress.von_mises.max_pa": {
         "value": 145e6,
         "limit": 180e6,
         "comparison": "<=",
         "passed": true,
         "margin": 35e6,
         "margin_pct": 19.4
       }
     },

     "artifacts": [
       {"kind": "mesh", "path": "mesh.msh"},
       {"kind": "solution", "path": "displacement.pvd"},
       {"kind": "log", "path": "solver.log"}
     ],

     "errors": [],
     "warnings": [],
     "notes": "Analysis completed successfully."
   }

Status Values
-------------

- ``passed`` - All acceptance criteria met
- ``failed`` - One or more criteria not met
- ``error`` - Execution error (solver crash, missing files, etc.)
- ``skipped`` - Plan skipped (missing backend, disabled, etc.)
- ``pending`` - Queued for remote execution


Backend Registration
====================

Backends are registered in yapCAD via the adapter framework:

.. code-block:: python

   from yapcad.package.analysis.base import (
       AnalysisAdapter, AnalysisPlan, AnalysisResult,
       register_backend
   )

   class MyAdapter(AnalysisAdapter):
       name = "my-solver"

       def run(self, manifest, plan: AnalysisPlan,
               workspace: Path, **kwargs) -> AnalysisResult:
           # Implement solver invocation
           ...
           return AnalysisResult(
               plan_id=plan.plan_id,
               status="passed",
               metrics={...}
           )

   register_backend("my-solver", MyAdapter)


Built-in Backends
-----------------

- **yapcad**: Native geometric/measurement checks (always available)
- **fenics**: FEniCSx-based FEA (requires dolfinx)
- **calculix**: CalculiX solver (requires ccx binary)


Example: Complete Package
=========================

.. code-block:: text

   my_design.ycpkg/
   ├── manifest.yaml
   ├── geometry/
   │   └── primary.json
   ├── validation/
   │   ├── plans/
   │   │   ├── volume-check.yaml      # Geometric test
   │   │   ├── mass-budget.yaml       # Measurement test
   │   │   └── structural-fea.yaml    # Simulation test
   │   └── results/
   │       ├── volume-check/
   │       │   └── summary.json
   │       ├── mass-budget/
   │       │   └── summary.json
   │       └── structural-fea/
   │           ├── summary.json
   │           ├── mesh.msh
   │           └── displacement.pvd
   └── exports/
       └── model.step


CLI Usage
=========

.. code-block:: bash

   # Run a specific validation plan
   python tools/ycpkg_analyze.py my_design.ycpkg/ --plan validation/plans/volume-check.yaml

   # Run all validation plans in a package
   python tools/ycpkg_analyze.py my_design.ycpkg/ --all

   # Generate validation report
   python tools/ycpkg_analyze.py my_design.ycpkg/ --report validation_report.html


Future Extensions
=================

- **JSON Schema files**: Formal JSON Schema documents for automated validation
- **DSL-based tests**: Define tests directly in DSL files
- **Parametric sweeps**: Run tests across parameter ranges
- **CI/CD integration**: GitHub Actions workflow templates
- **Result visualization**: Interactive result viewers
