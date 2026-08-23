"""Straight bevel gears generated from a spherical involute tooth form.

The pitch apex is at the origin and the gear axis is +Z.  Section planes are
normal to +Z, with the outer (large) end farther from the apex.  A tooth
centerline lies on +X at zero phase.
"""

from dataclasses import dataclass
import math
from typing import List, Literal, Tuple


GenerationType = Literal["spherical_involute", "octoid", "gleason"]
Point3 = Tuple[float, float, float]


@dataclass(frozen=True)
class StraightBevelGearSpec:
    teeth: int
    mate_teeth: int
    outer_module_mm: float
    face_width_mm: float
    shaft_angle_deg: float = 90.0
    pressure_angle_deg: float = 20.0
    backlash_mm: float = 0.0
    bore_diameter_mm: float = 0.0
    generation_type: GenerationType = "spherical_involute"


@dataclass(frozen=True)
class StraightBevelGearGeometry:
    pitch_cone_angle_rad: float
    mate_pitch_cone_angle_rad: float
    base_cone_angle_rad: float
    root_cone_angle_rad: float
    tip_cone_angle_rad: float
    outer_pitch_radius_mm: float
    outer_cone_distance_mm: float
    inner_cone_distance_mm: float
    outer_plane_z_mm: float
    inner_plane_z_mm: float
    pitch_tooth_half_angle_rad: float


def _finite(name: str, value: float) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def derive_straight_bevel_geometry(
    spec: StraightBevelGearSpec,
) -> StraightBevelGearGeometry:
    """Validate a specification and derive its pitch/base/addendum cones."""
    if not isinstance(spec.teeth, int) or isinstance(spec.teeth, bool) or spec.teeth < 6:
        raise ValueError("teeth must be an integer greater than or equal to 6")
    if (
        not isinstance(spec.mate_teeth, int)
        or isinstance(spec.mate_teeth, bool)
        or spec.mate_teeth < 6
    ):
        raise ValueError("mate_teeth must be an integer greater than or equal to 6")

    known_modes = {"spherical_involute", "octoid", "gleason"}
    if not isinstance(spec.generation_type, str) or spec.generation_type not in known_modes:
        raise TypeError(
            "generation_type must be 'spherical_involute', 'octoid', or 'gleason'"
        )
    if spec.generation_type != "spherical_involute":
        raise NotImplementedError(
            f"straight bevel gear generation type '{spec.generation_type}' is "
            "reserved but not implemented"
        )

    module = _finite("outer_module_mm", spec.outer_module_mm)
    face_width = _finite("face_width_mm", spec.face_width_mm)
    shaft_angle = math.radians(_finite("shaft_angle_deg", spec.shaft_angle_deg))
    pressure_angle = math.radians(
        _finite("pressure_angle_deg", spec.pressure_angle_deg)
    )
    backlash = _finite("backlash_mm", spec.backlash_mm)
    bore = _finite("bore_diameter_mm", spec.bore_diameter_mm)
    if module <= 0.0:
        raise ValueError("outer_module_mm must be positive")
    if face_width <= 0.0:
        raise ValueError("face_width_mm must be positive")
    if not 0.0 < shaft_angle < math.pi:
        raise ValueError("shaft_angle_deg must be between 0 and 180 degrees")
    if not 0.0 < pressure_angle < math.pi / 2.0:
        raise ValueError("pressure_angle_deg must be between 0 and 90 degrees")
    if backlash < 0.0:
        raise ValueError("backlash_mm cannot be negative")
    if bore < 0.0:
        raise ValueError("bore_diameter_mm cannot be negative")

    ratio = spec.mate_teeth / spec.teeth
    pitch_angle = math.atan2(
        math.sin(shaft_angle), ratio + math.cos(shaft_angle)
    )
    mate_pitch_angle = shaft_angle - pitch_angle
    pitch_radius = module * spec.teeth / 2.0
    cone_distance = pitch_radius / math.sin(pitch_angle)
    if face_width >= cone_distance:
        raise ValueError("face_width_mm must be less than the outer cone distance")

    circular_pitch = math.pi * module
    if backlash >= circular_pitch / 2.0:
        raise ValueError("backlash_mm must be less than half the circular pitch")
    half_tooth_angle = (
        math.pi / (2.0 * spec.teeth) - backlash / (2.0 * pitch_radius)
    )
    if half_tooth_angle <= 0.0:
        raise ValueError("backlash_mm leaves no positive tooth thickness")

    # The base cone is the spherical analogue of a spur gear base circle.
    base_angle = math.asin(math.sin(pitch_angle) * math.cos(pressure_angle))
    root_angle = pitch_angle - math.atan2(1.25 * module, cone_distance)
    tip_angle = pitch_angle + math.atan2(module, cone_distance)
    if root_angle <= 0.0 or tip_angle >= math.pi / 2.0:
        raise ValueError("tooth proportions produce an invalid root or tip cone")

    inner_distance = cone_distance - face_width
    return StraightBevelGearGeometry(
        pitch_cone_angle_rad=pitch_angle,
        mate_pitch_cone_angle_rad=mate_pitch_angle,
        base_cone_angle_rad=base_angle,
        root_cone_angle_rad=root_angle,
        tip_cone_angle_rad=tip_angle,
        outer_pitch_radius_mm=pitch_radius,
        outer_cone_distance_mm=cone_distance,
        inner_cone_distance_mm=inner_distance,
        outer_plane_z_mm=cone_distance * math.cos(pitch_angle),
        inner_plane_z_mm=inner_distance * math.cos(pitch_angle),
        pitch_tooth_half_angle_rad=half_tooth_angle,
    )


def spherical_involute_point(base_cone_angle_rad: float, unwind_rad: float) -> Point3:
    """Return a unit-sphere point on an involute of the base cone.

    ``unwind_rad`` is the great-circle distance unwound from the base circle.
    The construction is the envelope of great circles tangent to that circle.
    """
    base = float(base_cone_angle_rad)
    unwind = float(unwind_rad)
    if not 0.0 < base < math.pi / 2.0:
        raise ValueError("base_cone_angle_rad must be between 0 and pi/2")
    if unwind < 0.0:
        raise ValueError("unwind_rad cannot be negative")
    azimuth = -unwind / math.sin(base)
    sin_base = math.sin(base)
    c = (
        sin_base * math.cos(azimuth),
        sin_base * math.sin(azimuth),
        math.cos(base),
    )
    tangent = (-math.sin(azimuth), math.cos(azimuth), 0.0)
    return tuple(
        math.cos(unwind) * c[index] + math.sin(unwind) * tangent[index]
        for index in range(3)
    )  # type: ignore[return-value]


def _flank_azimuth(base_angle: float, cone_angle: float) -> float:
    if cone_angle <= base_angle:
        return 0.0
    argument = max(-1.0, min(1.0, math.cos(cone_angle) / math.cos(base_angle)))
    unwind = math.acos(argument)
    point = spherical_involute_point(base_angle, unwind)
    return math.atan2(point[1], point[0])


def tooth_section_points(
    spec: StraightBevelGearSpec,
    geometry: StraightBevelGearGeometry | None = None,
    *,
    cone_distance: Literal["outer", "inner"] = "outer",
    flank_samples: int = 9,
) -> List[Point3]:
    """Build a closed, planar tooth-outline section for BREP lofting."""
    geometry = geometry or derive_straight_bevel_geometry(spec)
    if flank_samples < 3:
        raise ValueError("flank_samples must be at least 3")
    if cone_distance not in {"outer", "inner"}:
        raise ValueError("cone_distance must be 'outer' or 'inner'")
    plane_z = (
        geometry.outer_plane_z_mm
        if cone_distance == "outer"
        else geometry.inner_plane_z_mm
    )

    base = geometry.base_cone_angle_rad
    pitch_raw = _flank_azimuth(base, geometry.pitch_cone_angle_rad)
    deltas = [
        geometry.root_cone_angle_rad
        + (geometry.tip_cone_angle_rad - geometry.root_cone_angle_rad)
        * index
        / (flank_samples - 1)
        for index in range(flank_samples)
    ]
    offsets = [
        geometry.pitch_tooth_half_angle_rad
        + _flank_azimuth(base, delta)
        - pitch_raw
        for delta in deltas
    ]

    points: List[Point3] = []
    tooth_pitch = 2.0 * math.pi / spec.teeth
    for tooth in range(spec.teeth):
        center = tooth * tooth_pitch
        # Traverse the left flank root-to-tip, then the right flank tip-to-root.
        for delta, offset in zip(deltas, offsets):
            radius = plane_z * math.tan(delta)
            angle = center - offset
            points.append((radius * math.cos(angle), radius * math.sin(angle), plane_z))
        for delta, offset in reversed(list(zip(deltas, offsets))):
            radius = plane_z * math.tan(delta)
            angle = center + offset
            points.append((radius * math.cos(angle), radius * math.sin(angle), plane_z))
    points.append(points[0])
    return points


def _make_ruled_loft_brep(
    inner: List[Point3],
    outer: List[Point3],
    *,
    teeth: int,
    flank_samples: int,
):
    """Create a strict ruled OCC loft; never substitutes a mesh result."""
    from yapcad.brep import require_occ

    require_occ()
    from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire
    from OCC.Core.BRepCheck import BRepCheck_Analyzer
    from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections
    from OCC.Core.GeomAPI import GeomAPI_Interpolate
    from OCC.Core.TColgp import TColgp_HArray1OfPnt
    from OCC.Core.gp import gp_Pnt

    def curve_edge(points: List[Point3]):
        poles = TColgp_HArray1OfPnt(1, len(points))
        for index, coordinates in enumerate(points, 1):
            poles.SetValue(index, gp_Pnt(*coordinates))
        interpolation = GeomAPI_Interpolate(poles, False, 1.0e-7)
        interpolation.Perform()
        if not interpolation.IsDone():
            raise RuntimeError("failed to interpolate a bevel gear tooth flank")
        return BRepBuilderAPI_MakeEdge(interpolation.Curve()).Edge()

    def make_wire(points: List[Point3]):
        wire = BRepBuilderAPI_MakeWire()
        stride = 2 * flank_samples
        if len(points) != teeth * stride:
            raise RuntimeError("straight bevel gear section topology mismatch")
        for tooth in range(teeth):
            start = tooth * stride
            left = points[start:start + flank_samples]
            right = points[start + flank_samples:start + stride]
            next_left_root = points[((tooth + 1) % teeth) * stride]
            wire.Add(curve_edge(left))
            wire.Add(BRepBuilderAPI_MakeEdge(
                gp_Pnt(*left[-1]), gp_Pnt(*right[0])
            ).Edge())
            wire.Add(curve_edge(right))
            wire.Add(BRepBuilderAPI_MakeEdge(
                gp_Pnt(*right[-1]), gp_Pnt(*next_left_root)
            ).Edge())
        if not wire.IsDone():
            raise RuntimeError("failed to construct straight bevel gear section wire")
        return wire.Wire()

    loft = BRepOffsetAPI_ThruSections(True, True, 1.0e-7)
    loft.CheckCompatibility(False)
    loft.AddWire(make_wire(inner))
    loft.AddWire(make_wire(outer))
    loft.Build()
    if not loft.IsDone():
        raise RuntimeError("OCC failed to build the straight bevel gear ruled loft")
    shape = loft.Shape()
    if not BRepCheck_Analyzer(shape).IsValid():
        raise RuntimeError("OCC produced an invalid straight bevel gear BREP")
    return shape


def make_straight_bevel_gear(
    spec: StraightBevelGearSpec,
    *,
    flank_samples: int = 5,
):
    """Generate an OCC-backed yapCAD solid for a straight bevel gear.

    BREP support is deliberately mandatory: a faceted fallback would hide the
    very tooth/contact errors this primitive is intended to prevent.
    """
    from yapcad.brep import BrepSolid, attach_brep_to_solid
    from yapcad.geom3d import solid

    geometry = derive_straight_bevel_geometry(spec)
    inner = tooth_section_points(
        spec, geometry, cone_distance="inner", flank_samples=flank_samples
    )[:-1]
    outer = tooth_section_points(
        spec, geometry, cone_distance="outer", flank_samples=flank_samples
    )[:-1]
    shape = _make_ruled_loft_brep(
        inner,
        outer,
        teeth=spec.teeth,
        flank_samples=flank_samples,
    )
    if spec.bore_diameter_mm > 0.0:
        from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut
        from OCC.Core.BRepCheck import BRepCheck_Analyzer
        from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeCylinder
        from OCC.Core.gp import gp_Ax2, gp_Dir, gp_Pnt

        cutter_height = geometry.outer_plane_z_mm + 2.0
        cutter = BRepPrimAPI_MakeCylinder(
            gp_Ax2(gp_Pnt(0.0, 0.0, -1.0), gp_Dir(0.0, 0.0, 1.0)),
            spec.bore_diameter_mm / 2.0,
            cutter_height,
        ).Shape()
        cut = BRepAlgoAPI_Cut(shape, cutter)
        cut.Build()
        if not cut.IsDone() or not BRepCheck_Analyzer(cut.Shape()).IsValid():
            raise RuntimeError("OCC failed to cut a valid straight bevel gear bore")
        shape = cut.Shape()

    brep = BrepSolid(shape)
    result = solid([brep.tessellate(deflection=0.25)], [], [
        "procedure", "yapcad.gears.bevel.make_straight_bevel_gear(spec)",
    ], {
        "gear": {
            "kind": "straight_bevel",
            "generation_type": spec.generation_type,
            "teeth": spec.teeth,
            "mate_teeth": spec.mate_teeth,
            "outer_module_mm": spec.outer_module_mm,
            "shaft_angle_deg": spec.shaft_angle_deg,
        }
    })
    attach_brep_to_solid(result, brep)
    return result


def make_straight_bevel_pair(
    spec: StraightBevelGearSpec,
    *,
    driver_angle_deg: float = 0.0,
    flank_samples: int = 5,
):
    """Return a correctly phased external straight-bevel pair at one apex.

    The first gear remains on +Z.  The mate axis is obtained by rotating +Z
    about +Y by ``shaft_angle_deg``.  Positive driver rotation produces the
    external-gear relation ``q_mate = phase - q_driver*z_driver/z_mate``.
    """
    from yapcad.geom import point
    from yapcad.geom3d import rotate

    derive_straight_bevel_geometry(spec)
    mate_spec = StraightBevelGearSpec(
        teeth=spec.mate_teeth,
        mate_teeth=spec.teeth,
        outer_module_mm=spec.outer_module_mm,
        face_width_mm=spec.face_width_mm,
        shaft_angle_deg=spec.shaft_angle_deg,
        pressure_angle_deg=spec.pressure_angle_deg,
        backlash_mm=spec.backlash_mm,
        bore_diameter_mm=spec.bore_diameter_mm,
        generation_type=spec.generation_type,
    )
    driver = make_straight_bevel_gear(spec, flank_samples=flank_samples)
    mate = make_straight_bevel_gear(mate_spec, flank_samples=flank_samples)
    driver = rotate(
        driver, driver_angle_deg, point(0, 0, 0), point(0, 0, 1)
    )
    mate = rotate(
        mate, spec.shaft_angle_deg, point(0, 0, 0), point(0, 1, 0)
    )
    shaft_angle = math.radians(spec.shaft_angle_deg)
    mate_axis = point(math.sin(shaft_angle), 0.0, math.cos(shaft_angle))
    phase_deg = 180.0 / spec.mate_teeth
    mate_angle_deg = (
        phase_deg - driver_angle_deg * spec.teeth / spec.mate_teeth
    )
    mate = rotate(mate, mate_angle_deg, point(0, 0, 0), mate_axis)
    return driver, mate
