"""Data structures for manufacturing post-processing.

This module defines the core data structures used for beam segmentation
and assembly planning.
"""

from dataclasses import dataclass, field
from typing import Any, List, Optional, Tuple


@dataclass
class CutPoint:
    """Specification for a segmentation cut.

    Defines where and how to cut a swept element to create segments
    that fit within build volume constraints.

    Attributes:
        element_id: Identifier of the swept element to cut
        parameter: Location along spine as t in [0, 1]
        connector_length: Override for auto-computed connector length (mm)
        fit_clearance: Dimensional offset for fit (mm per side)
        union_connector_with: Which segment gets the connector tab
            - "a": First segment (lower t)
            - "b": Second segment (higher t)
            - "none": Connector is separate piece
    """
    element_id: str
    parameter: float
    connector_length: Optional[float] = None
    fit_clearance: float = 0.2
    union_connector_with: str = "a"

    def __post_init__(self):
        if not 0 < self.parameter < 1:
            raise ValueError(f"Cut parameter must be in (0, 1), got {self.parameter}")
        if self.union_connector_with not in ("a", "b", "none"):
            raise ValueError(f"union_connector_with must be 'a', 'b', or 'none'")


@dataclass
class Segment:
    """A segment resulting from splitting a swept element.

    Attributes:
        id: Unique identifier for this segment
        solid: yapCAD solid geometry
        parent_element_id: ID of the original swept element
        parameter_range: (t_start, t_end) along parent spine
        has_connector_tab: True if connector is unioned with this segment
        connector_type: Type of connector mating ("male", "female", "none")
            - "male": Has protruding connector tab
            - "female": Hollow interior receives male tab
            - "none": No connector at this boundary
        mates_with: IDs of segments this connects to
        bounding_box: [[xmin,ymin,zmin,1], [xmax,ymax,zmax,1]]
    """
    id: str
    solid: Any  # yapCAD solid
    parent_element_id: str
    parameter_range: Tuple[float, float]
    has_connector_tab: bool = False
    connector_type: str = "none"  # "male", "female", or "none"
    mates_with: List[str] = field(default_factory=list)
    bounding_box: Optional[List] = None


@dataclass
class ConnectorSpec:
    """Specification for an interior connector.

    Attributes:
        id: Unique identifier
        parent_element_id: ID of the swept element this connects
        center_parameter: Location along spine where connector is centered
        length: Total connector length (extends equally both sides of center)
        fit_clearance: Dimensional offset applied to profile
        profile_type: Type of connector profile ("box", "circular", "custom")
    """
    id: str
    parent_element_id: str
    center_parameter: float
    length: float
    fit_clearance: float = 0.2
    profile_type: str = "box"


@dataclass
class SegmentationResult:
    """Complete result of a segmentation operation.

    Attributes:
        segments: List of resulting segment objects
        connectors: List of connector specifications (for reference)
        assembly_graph: Maps segment_id -> list of mating segment_ids
        build_volume_ok: True if all segments fit target build volume
        warnings: Any issues detected during segmentation
        assembly_instructions: Human-readable assembly sequence
    """
    segments: List[Segment]
    connectors: List[ConnectorSpec] = field(default_factory=list)
    assembly_graph: dict = field(default_factory=dict)
    build_volume_ok: bool = True
    warnings: List[str] = field(default_factory=list)
    assembly_instructions: str = ""

    @property
    def segment_count(self) -> int:
        """Number of segments produced."""
        return len(self.segments)

    def get_segment(self, segment_id: str) -> Optional[Segment]:
        """Get a segment by ID."""
        for seg in self.segments:
            if seg.id == segment_id:
                return seg
        return None

    def get_segments_for_element(self, element_id: str) -> List[Segment]:
        """Get all segments from a specific parent element."""
        return [s for s in self.segments if s.parent_element_id == element_id]


@dataclass
class SweptElementProvenance:
    """Provenance data for a swept element.

    Captures how a swept solid was created, enabling intelligent
    segmentation that preserves structural intent.

    Attributes:
        id: Unique identifier for this element
        operation: How it was created ("sweep", "sweep_adaptive", etc.)
        outer_profile: The outer boundary region2d
        inner_profile: Inner void region2d (None for solid profiles)
        spine: The path3d the profile was swept along
        wall_thickness: For hollow profiles, the wall thickness
        semantic_type: Design intent ("structural_beam", "decorative", etc.)
    """
    id: str
    operation: str
    outer_profile: Any  # region2d
    spine: Any  # path3d
    inner_profile: Optional[Any] = None  # region2d or None
    wall_thickness: Optional[float] = None
    semantic_type: str = "structural_beam"
    metadata: dict = field(default_factory=dict)
