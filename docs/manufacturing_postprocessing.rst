Manufacturing Post-Processing Framework
======================================

Design Document v0.2 - January 2026

.. note::

    **Implementation Status: Phase 1 Complete**

    The core beam segmentation functionality is implemented and working:

    - Swept element segmentation at specified cut points
    - Interior connector specification generation
    - STEP file export for individual segments
    - Assembly manifest generation with mating relationships

    See ``examples/globe_stand/segment_globe_stand.py`` for a working example.

Quick Start
-----------

Run the globe stand segmentation example:

.. code-block:: bash

    # Segment support arcs at 200mm max length, export STEP files
    python examples/globe_stand/segment_globe_stand.py

    # Segment at 150mm max length
    python examples/globe_stand/segment_globe_stand.py 150

    # Skip STEP export (faster for testing)
    python examples/globe_stand/segment_globe_stand.py --no-step

This produces:

- Individual STEP files for each segment (e.g., ``support_arc_1_seg_0.step``)
- ``assembly_manifest.json`` with part relationships and file paths

Basic Python usage:

.. code-block:: python

    from yapcad.manufacturing import (
        SweptElementProvenance,
        CutPoint,
        segment_swept_element,
        compute_optimal_cuts,
    )

    # Create provenance for your swept solid
    provenance = SweptElementProvenance(
        id="my_beam",
        operation="sweep_adaptive",
        outer_profile=outer_region2d,
        spine=path3d,
        wall_thickness=2.0,
        metadata={'solid': my_solid}  # The actual solid
    )

    # Compute cuts for max 200mm segments
    cuts = compute_optimal_cuts(provenance, max_segment_length=200.0)

    # Segment the solid
    result = segment_swept_element(provenance, cuts)
    print(f"Created {result.segment_count} segments")


Overview
--------

This document specifies a post-processing framework for yapCAD that enables
manufacturing of parts that exceed machine capabilities (build volume, reach, etc.)
by intelligently splitting monolithic designs into assembleable sub-parts.

The initial focus is on **beam segmentation** - splitting hollow swept structures
(box beams, tubes, channels) into printable segments with integrated connectors.
This framework is designed to be extensible to other manufacturing post-processing
tasks including custom toolpath generation.

Problem Statement
-----------------

**Use Case: Large 3D Printed Parts**

A designer creates a parametric model (e.g., a globe stand) that exceeds the
build volume of their printer. The model consists of structural elements that
are hollow box-beam or tube sections swept along 2D or 3D curves.

Current solutions:

1. **Manual splitting in CAD** - Time-consuming, error-prone, loses parametric benefits
2. **Mesh-based slicing tools** - Lose semantic understanding, can't create smart connectors
3. **Scaling down** - Compromises design intent

Desired solution:

- Automated or semi-automated splitting at semantically meaningful locations
- Connectors that maintain structural integrity
- Parametric - regenerates correctly when source design changes
- Works with arbitrary curved swept geometry

Terminology
-----------

**Swept Element**
    A solid created by sweeping a 2D profile along a 2D or 3D path (spine).
    Includes box beams, tubes, channels, and other hollow or solid extrusions.

**Segment**
    A portion of a swept element between two cut planes.

**Cut Point**
    A location along a swept element's spine where a segmentation cut is made.
    Defined by a parameter t in [0, 1] along the spine.

**Cut Plane**
    The plane perpendicular to the spine tangent at a cut point.

**Interior Connector**
    A solid that fits inside the hollow interior of a swept element, used to
    join two segments. Follows the same spine curve as the parent element.

**Connector Tab**
    An interior connector that is unioned with one segment, creating a
    male-female assembly interface.

**Fit Clearance**
    The dimensional offset applied to connector cross-sections to achieve
    the desired fit (press-fit, slip-fit, etc.).

**Provenance Metadata**
    Information retained from the original modeling operations that describes
    how geometry was created (profiles, spines, sweep parameters, etc.).


Approach: Beam Segmentation with Interior Connectors
----------------------------------------------------

Core Algorithm
^^^^^^^^^^^^^^

For a hollow swept element with known profile and spine:

1. **Identify cut points** along the spine (user-specified or computed)

2. **For each cut point:**

   a. Compute the cut plane (perpendicular to spine tangent)

   b. Split the swept solid into two segments at this plane

   c. Create an interior connector:

      - Extract the interior void profile (or compute from outer profile - wall thickness)
      - Apply fit clearance (shrink profile slightly)
      - Determine connector length based on profile size and curvature
      - Sweep the connector profile along the spine segment centered on cut point

   d. Union the connector with one segment (creating "male" side)

   e. The other segment becomes the "female" side

3. **Validate** that resulting segments fit within target build volume

4. **Generate** assembly instructions (which connectors mate with which segments)


Connector Design Details
^^^^^^^^^^^^^^^^^^^^^^^^

**Cross-Section Derivation**

For a box beam with outer profile ``[w, h]`` and wall thickness ``t``:

- Outer profile: rectangle ``w × h``
- Inner void: rectangle ``(w - 2t) × (h - 2t)``
- Connector profile: rectangle ``(w - 2t - 2c) × (h - 2t - 2c)`` where ``c`` is fit clearance

For profiles with fillets/chamfers, the connector profile should preserve these
features at the reduced scale for proper load distribution.

**Connector Length**

Minimum connector length ``L_min``:

- Straight sections: ``L_min = 3 × max(w, h)`` (3x largest profile dimension)
- Curved sections: ``L_min = max(3 × max(w, h), arc_length_for_15_degrees)``

The connector extends equally on both sides of the cut plane.

**Fit Clearance Values** (FDM printing defaults)

- Press-fit (structural): 0.15 - 0.20 mm per side
- Slip-fit (easy assembly): 0.25 - 0.35 mm per side
- Loose-fit (adjustable): 0.40 - 0.50 mm per side

These values are material and printer dependent; should be configurable.

**Curved Section Handling**

For curved spines, the connector must follow the same curve. This is handled
naturally by sweeping the connector profile along the appropriate segment of
the original spine. The sweep operation preserves curve fidelity.

For tight curves (radius < 5x profile dimension), consider:

- Longer connectors to span more of the curve
- Warning user about potential stress concentrations
- Suggesting alternative cut points


Provenance Requirements
-----------------------

Effective beam segmentation requires knowing how geometry was created.

**Required Provenance Data**

For each swept element, retain:

.. code-block:: yaml

    swept_element:
      id: "arc_beam_1"
      operation: "sweep" | "sweep_adaptive" | "sweep_hollow"
      outer_profile: <region2d reference or serialization>
      inner_profile: <region2d reference or null for solid>
      spine: <path3d reference or serialization>
      wall_thickness: <float, if applicable>
      metadata:
        semantic_type: "structural_beam" | "decorative" | "functional"
        material_hint: "PLA" | "PETG" | "ABS" | etc.

**DSL Integration**

The DSL ``emit`` statement should capture sweep provenance automatically when
emitting results of sweep operations. Enhanced emit:

.. code-block:: python

    # Current
    emit result

    # Enhanced with explicit tagging
    emit result, type="structural_beam", splittable=true

    # Or inferred from operation
    result = sweep_hollow(outer, inner, spine)
    emit result  # Automatically tagged as swept_element


API Design
----------

Core Functions
^^^^^^^^^^^^^^

.. code-block:: python

    from yapcad.manufacturing import (
        identify_swept_elements,
        compute_cut_points,
        segment_swept_element,
        create_interior_connector,
        validate_build_volume,
        generate_assembly_instructions,
    )

    # High-level workflow
    def segment_for_printing(
        solid,
        build_volume: tuple[float, float, float],
        *,
        fit_clearance: float = 0.2,
        connector_length_factor: float = 3.0,
        cut_points: list[CutPoint] = None,  # None = auto-compute
        provenance: dict = None,  # From DSL execution
    ) -> SegmentationResult:
        """
        Segment a solid for printing within specified build volume.

        Returns SegmentationResult containing:
        - segments: list of segment solids
        - connectors: list of connector solids (unioned with segments)
        - assembly_graph: how segments connect
        - warnings: any issues detected
        """

    # Lower-level functions for manual control
    def segment_swept_element(
        solid,
        profile: Region2D,
        spine: Path3D,
        cut_parameter: float,  # t in [0, 1]
        *,
        wall_thickness: float = None,
        fit_clearance: float = 0.2,
        connector_length: float = None,  # Auto-compute if None
    ) -> tuple[Solid, Solid, Solid]:
        """
        Segment a swept element at the specified parameter.

        Returns (segment_a, segment_b, connector)
        where connector is a separate solid (not yet unioned).
        """

    def create_interior_connector(
        outer_profile: Region2D,
        spine: Path3D,
        center_parameter: float,
        length: float,
        *,
        wall_thickness: float = None,
        inner_profile: Region2D = None,
        fit_clearance: float = 0.2,
    ) -> Solid:
        """
        Create an interior connector solid.

        If inner_profile is provided, uses it directly.
        Otherwise, derives from outer_profile and wall_thickness.
        """


Data Structures
^^^^^^^^^^^^^^^

.. code-block:: python

    @dataclass
    class CutPoint:
        """Specification for a segmentation cut."""
        element_id: str           # ID of swept element to cut
        parameter: float          # t in [0, 1] along spine
        connector_length: float = None  # Override auto-computed length
        fit_clearance: float = 0.2
        union_connector_with: str = "a"  # "a", "b", or "none"

    @dataclass
    class Segment:
        """A segment resulting from splitting."""
        id: str
        solid: Any  # yapCAD solid
        parent_element_id: str
        parameter_range: tuple[float, float]  # (t_start, t_end)
        has_connector_tab: bool
        mates_with: list[str]  # IDs of segments this connects to
        bounding_box: tuple  # For build volume validation

    @dataclass
    class SegmentationResult:
        """Complete result of segmentation operation."""
        segments: list[Segment]
        assembly_graph: dict  # segment_id -> list of mating segment_ids
        build_volume_ok: bool
        warnings: list[str]
        assembly_instructions: str  # Human-readable


User Interaction Model
----------------------

Three levels of automation:

**Level 1: Fully Automatic**

.. code-block:: python

    result = segment_for_printing(
        globe_stand_solid,
        build_volume=(256, 256, 256),  # Bambu X1C
        provenance=execution_result.provenance,
    )

    # System automatically:
    # - Identifies swept elements from provenance
    # - Computes optimal cut points to fit build volume
    # - Creates connectors and assembly plan

**Level 2: Semi-Automatic (Guided)**

.. code-block:: python

    # User identifies which elements to split
    cut_points = [
        CutPoint("arc_beam_1", parameter=0.33),
        CutPoint("arc_beam_1", parameter=0.67),
        CutPoint("base_ring", parameter=0.25),
        CutPoint("base_ring", parameter=0.50),
        CutPoint("base_ring", parameter=0.75),
    ]

    result = segment_for_printing(
        globe_stand_solid,
        build_volume=(256, 256, 256),
        cut_points=cut_points,
        provenance=execution_result.provenance,
    )

**Level 3: Interactive (Future)**

A visual tool that:

- Displays the model with swept elements highlighted
- Shows build volume overlay
- Lets user click to place cut points
- Shows real-time preview of resulting segments
- Warns about problematic cuts (stress concentrations, tight curves)


DSL Integration
---------------

New DSL Commands
^^^^^^^^^^^^^^^^

.. code-block:: python

    # Explicit segmentation in DSL
    command PRINTABLE_GLOBE_STAND(build_x: float, build_y: float, build_z: float) -> list<solid>:
        stand: solid = CENTERED_GLOBE_STAND_HIRES()

        # Automatic segmentation
        segments: list<solid> = segment_for_build_volume(
            stand,
            build_x, build_y, build_z,
            fit_clearance=0.2
        )

        emit segments

    # Or with explicit cut points
    command SEGMENTED_STAND() -> list<solid>:
        # Build the stand components with provenance
        arc1: solid = sweep_hollow(arc_profile, arc1_path, wall=3.0)
        arc2: solid = sweep_hollow(arc_profile, arc2_path, wall=3.0)
        arc3: solid = sweep_hollow(arc_profile, arc3_path, wall=3.0)
        base: solid = sweep_hollow(base_profile, base_path, wall=3.0)

        # Segment each arc into 3 pieces
        arc1_segments: list<solid> = segment_beam(arc1, [0.33, 0.67])
        arc2_segments: list<solid> = segment_beam(arc2, [0.33, 0.67])
        arc3_segments: list<solid> = segment_beam(arc3, [0.33, 0.67])

        # Segment base ring into 4 pieces
        base_segments: list<solid> = segment_beam(base, [0.25, 0.5, 0.75])

        emit concat(arc1_segments, arc2_segments, arc3_segments, base_segments)


CLI Integration
^^^^^^^^^^^^^^^

.. code-block:: bash

    # Segment for specific printer
    python -m yapcad.dsl run design.dsl MAKE_PART \
        --segment-for-printer "bambu_x1c" \
        --output segments/

    # Segment with custom build volume
    python -m yapcad.dsl run design.dsl MAKE_PART \
        --segment-build-volume 256,256,256 \
        --segment-clearance 0.2 \
        --output segments/

    # Export as assembly package
    python -m yapcad.dsl run design.dsl MAKE_PART \
        --segment-for-printer "bambu_x1c" \
        --package output.ycpkg \
        --include-assembly-instructions


Implementation Phases
---------------------

Phase 1: Core Segmentation (MVP) - **COMPLETE**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Goal:** Enable manual segmentation of swept elements with interior connectors.

**Deliverables:**

1. ✅ ``segment_swept_element()`` function - segments solids at cut planes
2. ✅ ``create_interior_connector()`` function - generates connector geometry
3. ✅ ``compute_optimal_cuts()`` function - computes cuts for max segment length
4. ✅ Basic provenance capture via ``SweptElementProvenance`` dataclass
5. ✅ Export segments as separate STEP files
6. ✅ Assembly manifest generation with mating relationships

**Implementation Location:** ``src/yapcad/manufacturing/``

- ``data.py`` - Data structures (CutPoint, Segment, SegmentationResult, etc.)
- ``path_utils.py`` - Path3D evaluation and sub-path extraction
- ``connectors.py`` - Interior connector generation
- ``segmentation.py`` - Core segmentation operations

**Example:** ``examples/globe_stand/segment_globe_stand.py``

**Scope:**

- Box beam (rectangular) profiles ✅
- Piecewise-linear (polyline) spines ✅
- Manual cut point specification ✅
- Automatic cut computation based on max length ✅
- Single wall thickness ✅
- STEP file export ✅
- Interior connector solid generation ✅


Phase 2: Automatic Cut Point Computation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Goal:** Automatically determine optimal cut points for a target build volume.

**Deliverables:**

1. Build volume analysis
2. Optimal cut point algorithm (minimize cuts while fitting volume)
3. Cut point validation (avoid stress concentrations)
4. Assembly graph generation

**Scope:**

- Multi-element assemblies
- Constraint satisfaction (all segments fit build volume)
- Basic optimization (minimize number of cuts)

**Estimated complexity:** Medium-High


Phase 3: Advanced Profiles and DSL Integration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Goal:** Support arbitrary profiles and integrate with DSL workflow.

**Deliverables:**

1. Arbitrary profile support (circles, complex polygons)
2. DSL builtins for segmentation
3. Enhanced provenance from DSL execution
4. Assembly instruction generation

**Scope:**

- Any Region2D profile
- Profiles with fillets/chamfers
- Wall thickness variations
- DSL-native workflow

**Estimated complexity:** Medium


Phase 4: Interactive Tools
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Goal:** Visual tools for guided segmentation.

**Deliverables:**

1. Viewer integration showing swept elements
2. Interactive cut point placement
3. Real-time segment preview
4. Build volume visualization

**Scope:**

- Integration with existing yapCAD viewer
- Mouse-based cut point placement
- Visual feedback for fit validation

**Estimated complexity:** High


Phase 5: Advanced Manufacturing Support
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Goal:** Extend framework to other manufacturing processes.

**Potential features:**

- **CNC segmentation:** Split for maximum machinable pocket depth
- **Sheet metal nesting:** Optimize 2D layout for laser/waterjet
- **Multi-material assembly:** Segment by material requirements
- **Toolpath generation:** Custom G-code for specific machines
- **Support structure optimization:** Modify geometry for better printability

**Estimated complexity:** Variable (feature-dependent)


Globe Stand Example
-------------------

The ``examples/globe_stand/segment_globe_stand.py`` script demonstrates the
manufacturing module by segmenting the support arcs from the Mars globe stand.

**Running the Example:**

.. code-block:: bash

    cd yapCAD
    PYTHONPATH=./src python examples/globe_stand/segment_globe_stand.py

**Example Output:**

.. code-block:: text

    ============================================================
    Globe Stand Segmentation Example
    ============================================================

    Parameters:
      Globe diameter: 304.8 mm
      Base diameter: 400.0 mm
      Beam profile: 10.0x10.0 mm, wall: 2.0 mm
      Max segment length: 200.0 mm
      STEP export: True

    Building and segmenting support arcs...
    ----------------------------------------

    Arc 0: base=0.0°, top=60.0°
      Path length: 460.6 mm
      Cuts needed: 2
        Cut 0: t=0.333
        Cut 1: t=0.667
      Created 3 segments
        Exported: support_arc_0_seg_0.step
        Exported: support_arc_0_seg_1.step
        Exported: support_arc_0_seg_2.step
      Building connector solids...
        Exported: support_arc_0_conn_33.step
        Exported: support_arc_0_conn_66.step

    Arc 1: base=120.0°, top=180.0°
      Path length: 396.1 mm
      Cuts needed: 1
        Cut 0: t=0.500
      Created 2 segments
        Exported: support_arc_1_seg_0.step
        Exported: support_arc_1_seg_1.step
      Building connector solids...
        Exported: support_arc_1_conn_50.step

    Arc 2: base=240.0°, top=300.0°
      Path length: 375.4 mm
      Cuts needed: 1
        Cut 0: t=0.500
      Created 2 segments
        Exported: support_arc_2_seg_0.step
        Exported: support_arc_2_seg_1.step
      Building connector solids...
        Exported: support_arc_2_conn_50.step

    ============================================================
    Summary
    ============================================================
      Total segments: 7
      Total connectors: 4
      STEP files exported: 11
      Output directory: examples/globe_stand/output

**Output Files:**

The script creates an ``output/`` directory containing:

- ``assembly_manifest.json`` - Part relationships and file paths
- ``support_arc_0_seg_0.step`` - First third of arc 0
- ``support_arc_0_seg_1.step`` - Middle third of arc 0
- ``support_arc_0_seg_2.step`` - Last third of arc 0
- ``support_arc_0_conn_33.step`` - Connector at first cut of arc 0
- ``support_arc_0_conn_66.step`` - Connector at second cut of arc 0
- ``support_arc_1_seg_0.step`` - First half of arc 1
- ``support_arc_1_seg_1.step`` - Second half of arc 1
- ``support_arc_1_conn_50.step`` - Connector at cut of arc 1
- ``support_arc_2_seg_0.step`` - First half of arc 2
- ``support_arc_2_seg_1.step`` - Second half of arc 2
- ``support_arc_2_conn_50.step`` - Connector at cut of arc 2

**Assembly Manifest Format:**

.. code-block:: json

    {
      "project": "Mars Globe Stand",
      "max_segment_length": 200.0,
      "beam_profile": "10.0x10.0mm hollow, 2.0mm wall",
      "fit_clearance": 0.18,
      "parts": [
        {
          "id": "support_arc_1_seg_0",
          "type": "arc_segment",
          "parent": "support_arc_1",
          "parameter_range": [0.0, 0.5],
          "mates_with": ["support_arc_1_seg_1"],
          "has_connector_tab": true,
          "step_file": "support_arc_1_seg_0.step"
        },
        ...
      ]
    }

**Connector Design:**

- Connector length: 3x largest profile dimension (30mm for 10x10mm beam)
- Press-fit clearance: 0.18mm per side (default for FDM printing)


Future Extensions
-----------------

**Toolpath Generation Framework**

Beyond segmentation, yapCAD could generate custom toolpaths for:

1. **Multi-axis CNC:**
   - Continuous 5-axis paths for complex surfaces
   - Adaptive clearing strategies
   - Tool change optimization

2. **Wire EDM:**
   - 2D profile extraction for wire paths
   - Taper angle computation
   - Start hole placement

3. **Additive manufacturing:**
   - Non-planar slicing for improved surface quality
   - Variable layer height based on geometry
   - Support structure generation with easy removal

4. **Hybrid manufacturing:**
   - Combined additive + subtractive strategies
   - Near-net-shape printing with finish machining
   - Selective surface finishing

**Machine Definition Framework**

.. code-block:: yaml

    machine:
      name: "Bambu Lab X1C"
      type: "fdm_printer"
      build_volume:
        x: 256
        y: 256
        z: 256
      materials:
        - PLA
        - PETG
        - ABS
      tolerances:
        xy_accuracy: 0.1
        z_accuracy: 0.05
        recommended_clearance: 0.2

This would enable printer-aware design validation and automatic segmentation
configuration.


Open Questions
--------------

1. **Connector orientation:** Should connectors always union with the "lower"
   segment (easier printing) or should user choose?

2. **Multi-wall beams:** How to handle beams with internal ribs or complex
   internal geometry?

3. **Junction handling:** How to segment at junctions where multiple beams meet?

4. **Validation depth:** How much structural analysis (FEA) should inform
   cut point selection?

5. **Assembly aids:** Should we generate alignment features (pins, keys) in
   addition to the main connector?


References
----------

- Manufacturing module: ``src/yapcad/manufacturing/``
- yapCAD sweep operations: ``src/yapcad/geom3d_util.py``
- BREP utilities: ``src/yapcad/brep.py``
- Globe stand DSL source: ``examples/globe_stand/globe_stand_v5.dsl``
- Globe stand segmentation example: ``examples/globe_stand/segment_globe_stand.py``
- Package specification: ``docs/ycpkg_spec.rst``
