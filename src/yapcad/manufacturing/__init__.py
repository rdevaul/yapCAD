"""Manufacturing post-processing for yapCAD.

This module provides tools for preparing yapCAD designs for manufacturing,
including beam segmentation for build volume constraints and interior
connector generation.

Phase 1 Implementation
----------------------
- Beam segmentation for hollow swept elements
- Interior connectors with configurable fit tolerances
- Path3D evaluation utilities
- Assembly planning with mating relationships

Example Usage
-------------
>>> from yapcad.manufacturing import (
...     CutPoint, segment_swept_element, SweptElementProvenance
... )
>>>
>>> # Define provenance for a swept beam
>>> provenance = SweptElementProvenance(
...     id="main_beam",
...     operation="sweep_adaptive",
...     outer_profile=outer_region2d,
...     spine=path3d,
...     wall_thickness=2.0,
...     metadata={'solid': beam_solid}
... )
>>>
>>> # Define cut points
>>> cuts = [
...     CutPoint("main_beam", parameter=0.33),
...     CutPoint("main_beam", parameter=0.67),
... ]
>>>
>>> # Segment the beam
>>> result = segment_swept_element(provenance, cuts)
>>> print(f"Created {result.segment_count} segments")
"""

# Data structures
from .data import (
    CutPoint,
    Segment,
    SegmentationResult,
    ConnectorSpec,
    SweptElementProvenance,
)

# Path utilities
from .path_utils import (
    evaluate_path3d_at_t,
    compute_cut_plane,
    extract_sub_path,
    path_length,
    length_to_parameter,
    parameter_to_length,
)

# Connector generation
from .connectors import (
    FIT_CLEARANCE,
    offset_rectangular_profile,
    compute_inner_profile_dimensions,
    compute_connector_profile_dimensions,
    compute_connector_length,
    create_connector_region2d,
    create_interior_connector,
    compute_connector_spec,
    create_terminal_connector,
    add_terminal_connectors_to_segment,
)

# Segmentation operations
from .segmentation import (
    split_solid_at_plane,
    segment_swept_element,
    segment_closed_ring,
    extract_segment_between_planes,
    compute_optimal_cuts,
    build_connector_solids,
)

# Ring generation
from .rings import (
    create_ring_spine,
    create_ring_profile,
    create_ring_solid,
    create_female_hole_solid,
    compute_arc_attachment_point,
    add_female_holes_to_ring,
    trim_segment_against_ring,
    compute_ring_cuts_avoiding_holes,
)

__all__ = [
    # Data structures
    "CutPoint",
    "Segment",
    "SegmentationResult",
    "ConnectorSpec",
    "SweptElementProvenance",
    # Path utilities
    "evaluate_path3d_at_t",
    "compute_cut_plane",
    "extract_sub_path",
    "path_length",
    "length_to_parameter",
    "parameter_to_length",
    # Connector generation
    "FIT_CLEARANCE",
    "offset_rectangular_profile",
    "compute_inner_profile_dimensions",
    "compute_connector_profile_dimensions",
    "compute_connector_length",
    "create_connector_region2d",
    "create_interior_connector",
    "compute_connector_spec",
    "create_terminal_connector",
    "add_terminal_connectors_to_segment",
    # Segmentation operations
    "split_solid_at_plane",
    "segment_swept_element",
    "segment_closed_ring",
    "extract_segment_between_planes",
    "compute_optimal_cuts",
    "build_connector_solids",
    # Ring generation
    "create_ring_spine",
    "create_ring_profile",
    "create_ring_solid",
    "create_female_hole_solid",
    "compute_arc_attachment_point",
    "add_female_holes_to_ring",
    "trim_segment_against_ring",
    "compute_ring_cuts_avoiding_holes",
]
