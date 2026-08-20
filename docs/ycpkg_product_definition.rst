yapCAD Package Product-Definition Extension
===========================================

**Proposal:** ``ycpkg-spec-v0.2``

**Status:** Implemented; rover validation pending

**Extends:** ``ycpkg-spec-v0.1`` and ``metadata-namespace-v1.1``


1. Purpose
----------

The original ``.ycpkg`` format is an effective container for geometry,
provenance, validation results, and derived exports.  A multi-part product also
needs an authoritative distinction between a component definition and each use
of that component.  Without that distinction, positioned assembly geometry can
be flattened into one solid and cannot produce a reliable bill of materials.

Version 0.2 adds a product-definition layer suitable for mechanisms such as the
YapRover suspension.  It supports fabricated parts, commercial off-the-shelf
(COTS) parts, raw stock, consumables, repeated instances, assembly constraints,
and derived BOM views.  Existing v0.1 packages remain valid.


2. Design Principles
--------------------

* Component identity is independent of instance identity and geometry entity
  identity.
* BOM quantities are derived from instances and per-instance quantity, never
  copied manually into a second source of truth.
* ``make``, ``buy``, ``raw_stock``, and ``consumable`` items have different
  validation requirements.
* COTS geometry may be simplified, imported, or absent; its procurement and
  interface definition remain authoritative.
* Local component geometry is canonical.  Instance transforms position it in
  the assembly without duplicating the definition.
* Assembly datums and mates describe kinematics.  Product manufacturing
  information (PMI) uses a separate namespace, even when the two reference the
  same physical feature.
* All references are package-relative and content-addressable by hash.


3. Manifest Additions
---------------------

A v0.2 manifest adds ``product``, ``components``, ``instances``, ``assembly``,
and ``bom``.  The existing geometry, source, validation, export, attachment,
material, provenance, and extension sections remain available.

.. code-block:: yaml

   schema: ycpkg-spec-v0.2
   name: yaprover-suspension
   version: 0.1.0
   units: mm

   product:
     rootAssembly: yaprover_suspension_detailed
     standards:
       dimensioning: project-defined
     lifecycle: prototype

   components:
     - id: wheel-hub
       name: Wheel hub
       partNumber: YR-WHL-001
       revision: A
       disposition: make
       geometry:
         path: geometry/entities/wheel-hub.json
         hash: sha256:...
       material: petg_black
       manufacturing:
         process: FDM
         postprocessing: [deburr, inspect-bearing-seats]

     - id: bearing-608-2rs
       name: 608-2RS bearing
       disposition: buy
       quantityPerInstance: 2
       procurement:
         specification: 608-2RS, 8x22x7 mm
         manufacturer: null
         manufacturerPartNumber: null
         approvedAlternates: []

   instances:
     - id: left-front-wheel
       component: wheel-hub
       transform: [[1, 0, 0, 190], [0, 1, 0, 84],
                   [0, 0, 1, 65], [0, 0, 0, 1]]
     - id: left-front-bearings
       component: bearing-608-2rs
       transform: [[1, 0, 0, 190], [0, 1, 0, 84],
                   [0, 0, 1, 65], [0, 0, 0, 1]]

   assembly:
     path: metadata/assembly.json
     hash: sha256:...
     rootPart: chassis

   bom:
     path: metadata/bom.json
     hash: sha256:...
     generatedFrom: instances


4. Component Contract
---------------------

Every component has:

``id``
  Stable package-local identifier.  Multiple instances may reference it.

``name``
  Human-readable description.

``disposition``
  One of ``make``, ``buy``, ``raw_stock``, or ``consumable``.

``quantityPerInstance``
  Positive numeric multiplier, defaulting to 1.  It supports visual surrogate
  parts such as a two-bearing pack while preserving a single purchased bearing
  BOM line.

Optional common fields are ``partNumber``, ``revision``, ``description``,
``material``, ``geometry``, ``manufacturing``, ``procurement``, ``pmi``, and
``documents``.

``make`` components SHOULD provide canonical geometry, material, revision, and
manufacturing process.  ``buy`` components SHOULD provide a manufacturer part
number or a complete manufacturer-neutral specification and MAY carry simplified
geometry.  ``raw_stock`` components SHOULD provide stock form and cut
requirements.  ``consumable`` components SHOULD provide a usage unit and
estimated quantity.


5. Instance and Assembly Contract
---------------------------------

Each instance has a unique ``id``, a valid ``component`` reference, and a rigid
4x4 transform in package units.  Instance order is stable so that BOM item
numbering and drawing balloons are reproducible.

``metadata/assembly.json`` records the native yapCAD assembly semantics:

* root part and solved transforms;
* named parts and their component references;
* named datums;
* mates and joint limits;
* joint values and affine joint couplings.

This file is not replaced by a positioned STEP assembly.  STEP is a derived
view; the native graph is the source of assembly intent.


6. BOM Derivation
-----------------

The canonical engineering BOM groups instances by component ID and sums
``quantityPerInstance``.  The generated record contains stable item numbers,
component ID, part number, revision, description, disposition, quantity, and
unit.  Additional views may aggregate supplier orders, raw-stock cut lists, and
fabrication jobs without modifying the engineering BOM.

The package validator verifies that every instance references a component,
component IDs and instance IDs are unique, quantities are positive, disposition
values are recognised, component geometry exists and matches its hash, and BOM
quantities agree with their source instances.


7. DSL Metadata Surface
-----------------------

The practical v0.2 bridge uses command decorators.  These values are copied from
the emitted part into its canonical component definition:

.. code-block:: text

   @meta(
       component.id="bearing-608-2rs",
       component.name="608-2RS bearing",
       component.disposition="buy",
       component.quantity=2,
       component.unit="each",
       procurement.specification="608-2RS, 8x22x7 mm",
       material="bearing-steel",
       assembly.datums=[...]
   )

Recognised component fields are ``id``, ``name``, ``description``,
``disposition``, ``part_number``, ``revision``, ``quantity``, and ``unit``.
Recognised procurement fields are ``manufacturer``, ``mpn``, ``supplier``,
``sku``, ``specification``, and ``approved_alternates``.  The manufacturing
namespace continues to carry process, instructions, fixtures, and
post-processing data.


8. PMI and GD&T Extension
-------------------------

Full drawing generation is a follow-on to the product/BOM implementation.  A
future ``pmi`` namespace will define:

* stable manufacturing features bound to persistent BREP topology or explicit
  construction geometry;
* datum features and ordered datum-reference frames;
* nominal, limit, basic, and fit dimensions;
* geometric tolerances and material-condition modifiers;
* surface texture, thread, insert, and finishing callouts;
* drawing sheets, views, sections, detail views, balloons, and revision blocks.

``assembly.datums`` MUST remain kinematic references.  A PMI datum feature may
link to an assembly datum, but generators must not infer a GD&T datum feature
from a kinematic axis or plane alone.

The first rover build may use reviewed PDF/DXF fabrication sheets and an
inspection plan stored under ``documents``.  Automatic semantic GD&T rendering
requires persistent feature naming and is not implied by v0.2.


9. Compatibility and Rollout
----------------------------

* Readers continue to accept v0.1 packages.
* Entity-only package creation continues to emit v0.1.
* Assembly-aware creation emits v0.2 and also writes a positioned
  ``geometry/primary.json`` for existing viewers.
* Geometry documents use ``yapcad-geometry-json-v0.2``: analytic BREP is the
  explicit authoritative representation and tessellation is a preview.
* Unknown v0.2 fields remain forward-compatible, but required references and
  BOM invariants are validated.
* Zipped-package support, semantic STEP PMI, manufacturing BOMs, and automated
  drawing layout remain separate milestones.


10. Acceptance Tests
--------------------

The reference rover package must prove that:

1. Thirty-seven assembly instances remain individually identifiable after
   package round-trip.
2. Repeated wheel, bearing, shaft, and fastener definitions produce one BOM
   line each with derived quantities.
3. Fabricated, COTS, and raw-stock components retain distinct dispositions.
4. Every instance transform, mate, joint limit, and differential coupling is
   preserved in the native assembly record.
5. Missing component references, duplicate IDs, invalid quantities, stale
   hashes, and edited BOM totals fail strict validation.
6. The positioned primary geometry remains viewable by v0.1-era tooling.
7. Fabrication drawings and inspection documents can be associated with their
   component definition and revision.
