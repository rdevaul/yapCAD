Historical Document
===================

.. note::

   This is the original DSL design proposal from 2025. The DSL has been fully
   implemented in yapCAD 1.0. For current syntax and usage, see
   ``docs/dsl_reference.md`` and ``docs/dsl_tutorial.md``.

---

yapCAD DSL Specification
========================

Status: Proposal
Version: 0.1.0
Date: October 2025

Overview
--------

This document describes a domain-specific language (DSL) for yapCAD that
enables declarative, parametric geometry definition with static type
checking and full provenance tracking.

Goals
~~~~~

1. **Declarative Syntax**: Define geometry through high-level operations
   rather than imperative code
2. **Type Safety**: Catch errors at compile time before geometry generation
3. **Provenance**: Track all inputs, parameters, and transformations
4. **Reproducibility**: Same inputs always produce identical outputs
5. **Interoperability**: Easy integration with yapCAD Python API

Non-Goals
~~~~~~~~~

- Full Turing completeness (deliberate limitation for termination guarantee)
- Real-time preview (compile-and-run model)
- Direct CAD file manipulation (focus on generation)

Language Design
---------------

Module Structure
~~~~~~~~~~~~~~~~

Every DSL file defines a module containing commands::

    module spur_gear

    command MAKE_GEAR(teeth: int, module_mm: float) -> solid:
        require teeth >= 6
        require module_mm > 0.0

        pitch_radius: float = teeth * module_mm / 2.0
        profile: region2d = involute_gear_profile(teeth, module_mm, 20.0)
        gear: solid = extrude(profile, 10.0)

        emit gear

Type System
~~~~~~~~~~~

The DSL uses a static type system organized into geometric tiers:

**Primitives**: ``int``, ``float``, ``bool``, ``string``

**Geometric Types**:

- ``point``, ``point2d``, ``point3d`` - positions
- ``vector``, ``vector2d``, ``vector3d`` - directions
- ``region2d`` - closed 2D shapes (can be extruded)
- ``path2d``, ``path3d`` - curve sequences
- ``solid`` - 3D volumes

**Collections**: ``list<T>`` - homogeneous lists

Commands
~~~~~~~~

Commands are the entry points for geometry generation::

    command NAME(param: type, param: type = default) -> return_type:
        # body statements
        emit result

- Parameters can have default values
- Return type specifies geometry output
- ``emit`` statement returns the final geometry

Statements
~~~~~~~~~~

**Let** - Variable binding::

    thickness: float = 10.0
    profile: region2d = rectangle(100.0, 50.0)

**Require** - Precondition validation::

    require teeth >= 6, "Minimum 6 teeth required"
    require thickness > 0.0

**Emit** - Return geometry::

    emit gear

Built-in Functions
------------------

Primitives
~~~~~~~~~~

- ``box(width, depth, height)`` → solid
- ``cylinder(radius, height)`` → solid
- ``sphere(radius)`` → solid
- ``cone(bottom_radius, top_radius, height)`` → solid

2D Shapes
~~~~~~~~~

- ``rectangle(width, height)`` → region2d
- ``circle(radius)`` → region2d
- ``polygon(points: list<point2d>)`` → region2d

Operations
~~~~~~~~~~

- ``extrude(profile: region2d, height: float)`` → solid
- ``revolve(profile: region2d, axis: vector3d, angle: float)`` → solid
- ``union(a: solid, b: solid)`` → solid
- ``difference(a: solid, b: solid)`` → solid
- ``intersection(a: solid, b: solid)`` → solid

Transforms
~~~~~~~~~~

- ``translate(shape, x, y, z)`` → same type
- ``rotate(shape, axis: vector3d, angle: float)`` → same type
- ``scale(shape, sx, sy, sz)`` → same type

CLI Usage
---------

::

    # Type-check a DSL file
    python -m yapcad.dsl check myfile.dsl

    # Execute a command
    python -m yapcad.dsl run myfile.dsl MAKE_GEAR teeth=24 module_mm=2.0

    # Generate output
    python -m yapcad.dsl run myfile.dsl MAKE_GEAR --output gear.step

    # Create package
    python -m yapcad.dsl run myfile.dsl MAKE_GEAR --package output.ycpkg

Implementation Notes
--------------------

The DSL is implemented as:

1. **Lexer** (``lexer.py``) - Tokenizes source text
2. **Parser** (``parser.py``) - Builds AST from tokens
3. **Type Checker** (``checker.py``) - Validates types and semantics
4. **Interpreter** (``interpreter.py``) - Executes AST to generate geometry

All modules live in ``src/yapcad/dsl/``.

Future Extensions
-----------------

- Conditional expressions (``if/else``)
- List comprehensions
- User-defined types
- Module imports
- Python escape hatch for advanced cases

References
----------

- yapCAD Geometry API: ``yapcad.geom``, ``yapcad.geom3d``
- Package Specification: ``docs/ycpkg_spec.rst``
- BREP Implementation: ``docs/yapBREP.rst``
