"""Mechatron assembly-graph JSON exporter.

Translates a ``yapcad.assembly.Assembly`` into the ``GraphSnapshot`` JSON
shape consumed by Mechatron's ``assembly-graph`` MCP service (the
``save_json``/``load_json`` codec at
``crates/assembly-graph/src/graph.rs``). After export, a Mechatron client
can ``POST`` the JSON or feed it to ``AssemblyGraph::load_json`` to
materialize the same graph on the consumer side, where it becomes the
source of truth for kinematics, MJCF/URDF generation, design-constraint
validation, and the Agentic-1 Assembly Dashboard at port 8765.

Schema reference
----------------
The target shape is roughly::

    {
      "parts": [
        {
          "id":        "<part name>",
          "name":      "<part name>",
          "material":  "PETG",                     # PascalCase Material enum
          "process":   "FDM",                       # PascalCase ManufacturingProcess
          "datums": [
            {
              "name":       "shaft",
              "datum_type": "Axis",                # PascalCase DatumType
              "origin":     {"x": 0.0, "y": 0.0, "z": 0.0},
              "direction":  {"x": 0.0, "y": 0.0, "z": 1.0},  # optional
              "normal":     null,                              # optional
              "radius":     null                                # optional
            }
          ],
          "tags": []
        }
      ],
      "interfaces": [
        {
          "id":           "iface_0",
          "parent_part":  "nosecone",
          "child_part":   "bulkhead",
          "joint":        {"type": "Fixed"},        # tagged Joint enum
          "mate_constraints": [
            {"type": "Concentric",
             "parent_datum": "shaft",
             "child_datum": "shaft"}
          ]
        }
      ],
      "last_script": "yapcad-dsl: forward_section"   # provenance hint
    }

What this exporter does
-----------------------
- Walks ``asm.parts`` and emits one ``Part`` per ``PartDefinition``. Material
  is taken from ``PartDefinition.material`` (defaults to ``PETG``) and mapped
  to the Mechatron Material PascalCase enum. Process defaults to ``FDM``.
- For each part, walks its datums and emits one Mechatron Datum per
  yapcad Datum. yapCAD's homogeneous coordinates (``[x, y, z, w]``) are
  dropped to ``Point3D``/``Vec3D`` shapes.
- Walks ``asm.mates`` and emits one ``Interface`` per Mate. yapCAD's
  motion mates (REVOLUTE, PRISMATIC, SPHERICAL, RIGID) map onto Mechatron
  ``Joint::Revolute``/``Joint::Prismatic``/``Joint::Ball``/``Joint::Fixed``.
  Other yapCAD mate kinds (CONCENTRIC, COINCIDENT, TANGENT, ANGLE) collapse
  to ``Joint::Fixed`` and are emitted as Mechatron MateConstraint entries
  on the interface (so the downstream solver has the semantic context).
  yapCAD mate kinds without a Mechatron MateConstraint equivalent
  (PARALLEL, PERPENDICULAR, DISTANCE, PLANAR, PIN_SLOT, UNIVERSAL, SCREW,
  GEAR, RACK_PINION, CAM, SLOT) are dropped from the export.
- ``last_script`` is set to a provenance hint identifying the assembly
  by name so a round-trip through Mechatron preserves the source.

What this exporter does NOT do
------------------------------
- Solve the assembly. Mechatron's solver runs ``add_interface_solved``
  on load; we emit interfaces with their semantic mate constraints
  attached. Origin/orientation transforms come out of the Mechatron
  solver, not us.
- Export design constraints. yapCAD's ``Constraint`` system (the
  ``constraint.py`` module) has different semantics than Mechatron's
  ``DesignConstraint``; cross-mapping is deferred to a future pass.
- Round-trip Mechatron-only features: ``loop_closures``,
  ``joint_couplings``, ``keyframes``, ``generators``, ``deformables``.
  These come back as empty lists in the output JSON.
"""

from __future__ import annotations

import json
from typing import Any, Dict, Iterable, List, Optional


# ---------------------------------------------------------------------------
# Enum mappings
# ---------------------------------------------------------------------------

# yapcad.assembly.datum.DatumType → mechatron mcp-common DatumType
# (mechatron uses PascalCase serde-default enum names)
_DATUM_TYPE_MAP: Dict[str, str] = {
    "point": "Point",
    "axis": "Axis",
    "plane": "Plane",
    "frame": "Frame",
    "circle": "Circle",
}

# yapcad PartDefinition.material strings → mechatron Material PascalCase enum
# (Mechatron accepts the lowercase aliases too via serde, but we send the
# canonical form to keep output stable.)
_MATERIAL_MAP: Dict[str, str] = {
    "pla": "PLA",
    "petg": "PETG",
    "petg-cf": "PetgCf",
    "petg_cf": "PetgCf",
    "tpu": "TPU",
    "abs": "ABS",
    "asa": "ASA",
    "pc": "PC",
    "pa-cf": "PaCf",
    "pa_cf": "PaCf",
    "nylon": "Nylon",
    "aluminum": "Aluminum",
    "steel": "Steel",
    "stainless": "Stainless",
    "brass": "Brass",
    "titanium": "Titanium",
}

# yapcad MateType.value → mechatron Joint::* tagged enum variant.
# Mates not in this map are treated as Joint::Fixed and (if mappable)
# preserved as MateConstraint entries on the interface.
_MOTION_MATE_TO_JOINT: Dict[str, str] = {
    "rigid": "Fixed",
    "revolute": "Revolute",
    "prismatic": "Prismatic",
    "spherical": "Ball",
}

# yapcad MateType.value → mechatron MateConstraint variant (the
# ``#[serde(tag = "type")]`` discriminator). Variants not listed here
# either map to a Joint above (motion mates), or have no Mechatron
# equivalent and are dropped from ``mate_constraints``.
#
# Mate kinds yapCAD supports today with NO Mechatron MateConstraint
# equivalent: parallel, perpendicular, distance, planar, pin_slot,
# universal, screw, gear, rack_pinion, cam, slot. Those remain in the
# yapCAD assembly object but are not exported. Future work: extend
# Mechatron's MateConstraint variants or add a yapCAD-private
# extension field for round-trip preservation.
_MATE_TO_CONSTRAINT_VARIANT: Dict[str, str] = {
    # yapCAD's CONCENTRIC maps to Mechatron's AxisCoincident, NOT
    # Concentric, on purpose. Mechatron's Concentric locks 3T + 2R
    # (5/6 DOFs); AxisCoincident locks 2T + 2R but it pairs with a
    # Coincident face mate via the validator's composite rule
    # ("AxisCoincident + Coincident => +1R axial rotation locked"),
    # which gets a Fixed joint to a clean 6/6. For yapcad authors
    # the semantics are the same: two cylinders sharing a centerline.
    "concentric": "AxisCoincident",
    "coincident": "Coincident",
    "tangent": "Tangent",
    "angle": "Angle",
    # PARALLEL maps to AxisAlign (a unilateral constraint that pins
    # the child's axis to a target direction in the parent frame).
    # AxisAlign has a different shape than the other constraints
    # (child_datum + target_axis only); see _mate_constraint_for_mate
    # for the special-case handling.
    "parallel": "AxisAlign",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _homogeneous_to_point3d(p) -> Dict[str, float]:
    """Drop yapCAD's homogeneous coordinate to a Mechatron Point3D."""
    if p is None:
        return {"x": 0.0, "y": 0.0, "z": 0.0}
    return {"x": float(p[0]), "y": float(p[1]), "z": float(p[2])}


def _homogeneous_to_vec3d(v) -> Optional[Dict[str, float]]:
    """Drop yapCAD's homogeneous coordinate to a Mechatron Vec3D.

    Returns ``None`` if the input is ``None`` (Mechatron expects
    ``Option<Vec3D>`` here, serialized as JSON null).
    """
    if v is None:
        return None
    return {"x": float(v[0]), "y": float(v[1]), "z": float(v[2])}


def _datum_to_mechatron(datum) -> Dict[str, Any]:
    """Convert a yapcad ``Datum`` into a Mechatron-shaped dict."""
    raw_kind = datum.datum_type.value if hasattr(datum.datum_type, "value") else str(datum.datum_type)
    pascal_kind = _DATUM_TYPE_MAP.get(raw_kind.lower(), raw_kind.capitalize())
    entry: Dict[str, Any] = {
        "name": datum.name,
        "datum_type": pascal_kind,
        "origin": _homogeneous_to_point3d(datum.origin),
    }
    direction = _homogeneous_to_vec3d(getattr(datum, "direction", None))
    if direction is not None:
        entry["direction"] = direction
    normal = _homogeneous_to_vec3d(getattr(datum, "normal", None))
    if normal is not None:
        entry["normal"] = normal
    radius = getattr(datum, "radius", None)
    if radius is not None:
        entry["radius"] = float(radius)
    return entry


def _material_to_mechatron(material: Optional[str]) -> str:
    """Map a yapCAD material string to a Mechatron Material variant.

    Unknown materials fall back to ``"PETG"`` (Mechatron's most-common
    default for printed parts) rather than raising; downstream tools
    surface a warning via the existing Material::is_unknown path.
    """
    if not material:
        return "PETG"
    return _MATERIAL_MAP.get(material.lower(), "PETG")


def _resolve_mate_axis(assembly, mate):
    """Look up the parent-side datum's direction vector for a motion mate.

    For ``Revolute``/``Prismatic`` joints, Mechatron's ``Joint`` variant
    requires an explicit axis vector in the parent's local frame. yapCAD
    stores this on the AXIS datum referenced by ``mate.datum_a`` on
    ``mate.part_a``. Returns ``None`` if the datum or direction is
    missing; the caller falls back to a +Z placeholder.
    """
    try:
        part = assembly.parts.get(str(mate.part_a))
        if part is None:
            return None
        datum = part.datums.get(str(mate.datum_a))
        if datum is None:
            return None
        direction = getattr(datum, "direction", None)
        if direction is None:
            return None
        # Homogeneous -> 3D vector; ignore the w-component.
        return [float(direction[0]), float(direction[1]), float(direction[2])]
    except Exception:
        return None


def _joint_for_mate(mate, assembly=None) -> Dict[str, Any]:
    """Build a Mechatron Joint tagged-enum dict from a yapcad Mate.

    ``Joint`` is ``#[serde(tag = "type")]`` in Mechatron, so the dict
    always has a ``"type"`` key plus variant-specific fields. For
    ``Revolute`` and ``Prismatic`` joints we resolve the rotation /
    translation axis from the parent-side datum's direction vector
    when the assembly is provided. If the datum lacks a direction (or
    the assembly handle was not threaded through), we fall back to a
    placeholder +Z axis and the Mechatron-side solver re-derives it
    from datum geometry on ``load_json``.
    """
    mate_kind = mate.mate_type.value if hasattr(mate.mate_type, "value") else str(mate.mate_type)
    joint_variant = _MOTION_MATE_TO_JOINT.get(mate_kind.lower(), "Fixed")
    if joint_variant in ("Revolute", "Prismatic"):
        axis = _resolve_mate_axis(assembly, mate) if assembly is not None else None
        if axis is None:
            axis = [0.0, 0.0, 1.0]
        return {"type": joint_variant, "axis": axis}
    return {"type": joint_variant}


def _mate_constraint_for_mate(mate, assembly=None,
                              joint_is_motion: bool = False) -> Optional[Dict[str, Any]]:
    """Build a Mechatron ``MateConstraint`` dict from a yapcad Mate.

    Returns ``None`` when the mate kind is a motion mate (its semantics
    are captured by the joint already) or has no Mechatron equivalent.

    Bilateral constraints (Concentric/Coincident/Tangent/Angle) use the
    standard tagged-enum shape::

        {"type": "Concentric", "parent_datum": "...", "child_datum": "..."}

    yapCAD's mate API uses ``datum_a`` / ``datum_b``; Mechatron uses
    ``parent_datum`` / ``child_datum``. The mapping is positional
    (datum_a -> parent_datum, datum_b -> child_datum).

    The ``AxisAlign`` variant is special — it is unilateral on the
    child side and pins the child's axis to a target direction in the
    parent frame::

        {"type": "AxisAlign", "child_datum": "...", "target_axis": [x,y,z]}

    For AxisAlign we look up the parent-side datum's direction vector
    on the assembly (when available) and use that as ``target_axis``.
    Falls back to +Z if the assembly handle is missing or the parent
    datum has no direction.
    """
    mate_kind = mate.mate_type.value if hasattr(mate.mate_type, "value") else str(mate.mate_type)
    variant = _MATE_TO_CONSTRAINT_VARIANT.get(mate_kind.lower())
    if variant is None:
        return None

    # CONCENTRIC has a context-dependent mapping. On a Fixed-joint
    # interface we map to AxisCoincident (2T + 2R) so the validator's
    # "AxisCoincident + Coincident => axial rotation locked" composite
    # rule fires and the Fixed mount comes out fully constrained. On a
    # motion-joint interface (Revolute/Prismatic/Ball) we map to
    # Concentric (3T + 2R) instead — Concentric does NOT trigger the
    # composite +1R, leaving the joint's natural rotational DOF free.
    if mate_kind.lower() == "concentric" and joint_is_motion:
        variant = "Concentric"

    if variant == "AxisAlign":
        target = _resolve_mate_axis(assembly, mate) if assembly is not None else None
        if target is None:
            target = [0.0, 0.0, 1.0]
        return {
            "type": "AxisAlign",
            "child_datum": str(mate.datum_b),
            "target_axis": target,
        }

    entry: Dict[str, Any] = {
        "type": variant,
        "parent_datum": str(mate.datum_a),
        "child_datum": str(mate.datum_b),
    }
    # Mechatron's ``Angle`` variant requires an ``angle_deg`` field;
    # default to 0 when the yapcad Mate does not carry an explicit angle.
    if variant == "Angle":
        entry["angle_deg"] = float(getattr(mate, "angle", 0.0) or 0.0)
    return entry


def _is_motion_mate(mate) -> bool:
    """True if this yapcad MateType maps to a Mechatron Joint variant.

    Motion mates (REVOLUTE/PRISMATIC/SPHERICAL/RIGID) define the joint
    DOF and are NOT also emitted as MateConstraint entries. Other
    mates (CONCENTRIC, COINCIDENT, TANGENT, ANGLE) leave the joint as
    Fixed and contribute geometric constraints used by the solver to
    place the parts.
    """
    mate_kind = mate.mate_type.value if hasattr(mate.mate_type, "value") else str(mate.mate_type)
    return mate_kind.lower() in _MOTION_MATE_TO_JOINT


def _resolve_insertion_vector(assembly, mates):
    """Derive a Mechatron insertion_vector from the parent-side mate datums.

    Mechatron requires every interface to declare the physical direction
    the child part is installed from, as a unit vector in the parent
    body frame. By convention, ``[0,0,-1]`` means "inserted from above"
    (the child slides down onto the parent).

    Strategy:
      1. Look at each non-motion mate in the group.
      2. Find the parent-side datum (``mate.datum_a`` on
         ``mate.part_a``) and inspect its kind.
      3. For a PLANE datum, the plane normal points OUT of the parent
         (away from the seated child). The insertion direction is the
         opposite — the child moves toward the parent along
         ``-normal``.
      4. For a CIRCLE datum, treat the circle normal the same way.
      5. AXIS datums and other kinds don't carry an unambiguous
         install direction; skip them.

    Returns ``None`` when no suitable datum is found; the caller emits
    no insertion_vector field (Mechatron will warn but not error).
    """
    if assembly is None:
        return None
    for mate in mates:
        try:
            part = assembly.parts.get(str(mate.part_a))
            if part is None:
                continue
            datum = part.datums.get(str(mate.datum_a))
            if datum is None:
                continue
            kind = datum.datum_type.value if hasattr(datum.datum_type, "value") else str(datum.datum_type)
            if kind.lower() not in ("plane", "circle"):
                continue
            normal = getattr(datum, "normal", None)
            if normal is None:
                continue
            # insertion_vector = -normal: the child moves opposite the
            # parent face's outward normal to seat against it. The
            # ``+ 0.0`` collapses ``-0.0`` to ``0.0`` for cleaner JSON
            # output (purely cosmetic — both serialize identically as
            # numbers).
            return [
                -float(normal[0]) + 0.0,
                -float(normal[1]) + 0.0,
                -float(normal[2]) + 0.0,
            ]
        except Exception:
            continue
    return None


def _interface_for_mate_group(idx: int, parent: str, child: str,
                              mates: list, assembly=None) -> Dict[str, Any]:
    """Convert a group of mates between the same (parent, child) pair
    into a single Mechatron Interface dict.

    Mechatron treats each Interface as one connection between two
    parts. yapCAD lets authors stack multiple mates (e.g. a
    ``concentric`` + ``coincident`` pair to fully constrain a stacked
    cylinder, plus a ``revolute`` for the kinematic joint). All of
    those describe the same physical relationship and belong on one
    Interface. Stripping them down:

      * Pick the first motion mate (REVOLUTE / PRISMATIC / SPHERICAL /
        RIGID) as the Joint. If none, the joint is Fixed.
      * All other mates become MateConstraint entries (where mappable),
        which the Mechatron solver uses to derive origin + orientation.
      * insertion_vector comes from a PLANE/CIRCLE parent datum
        referenced by the constraint mates.

    Args:
        idx: Sequential interface index (used for the ``iface_N`` id).
        parent: Parent part name.
        child: Child part name.
        mates: All yapcad Mates with ``mate.part_a == parent`` and
               ``mate.part_b == child``, in their original order.
        assembly: The owning Assembly, used to resolve datum
                  directions for joint axes and insertion_vector.

    Returns:
        A Mechatron-shaped Interface dict.
    """
    motion_mate = next((m for m in mates if _is_motion_mate(m)), None)
    if motion_mate is not None:
        joint = _joint_for_mate(motion_mate, assembly=assembly)
    else:
        joint = {"type": "Fixed"}

    # Every non-motion mate that maps to a MateConstraint variant.
    # Pass the assembly through so AxisAlign can resolve its target_axis
    # from the parent-side datum's direction vector.
    joint_is_motion = motion_mate is not None
    mate_constraints: List[Dict[str, Any]] = []
    for m in mates:
        if _is_motion_mate(m):
            continue
        c = _mate_constraint_for_mate(m, assembly=assembly,
                                      joint_is_motion=joint_is_motion)
        if c is not None:
            mate_constraints.append(c)

    iface: Dict[str, Any] = {
        "id": f"iface_{idx}",
        "parent_part": parent,
        "child_part": child,
        "joint": joint,
        "mate_constraints": mate_constraints,
    }
    ivec = _resolve_insertion_vector(assembly, mates)
    if ivec is not None:
        iface["insertion_vector"] = ivec
    return iface


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def to_mechatron_snapshot(assembly) -> Dict[str, Any]:
    """Build a Mechatron ``GraphSnapshot``-shaped dict from an Assembly.

    Args:
        assembly: A ``yapcad.assembly.Assembly`` instance, typically
                  built via the Phase-2 DSL builtins (``assembly``,
                  ``add_part``, ``add_mate``).

    Returns:
        A plain ``dict`` whose shape matches Mechatron's ``GraphSnapshot``
        serde representation. Pass through ``json.dumps`` to obtain a
        string consumable by ``AssemblyGraph::load_json``.

    Notes:
        - Order is preserved: parts are emitted in the order they were
          added to the assembly; mates likewise. The exporter is
          deterministic (no dict iteration order surprises).
        - Empty mate_constraints lists are still emitted on interfaces
          even though Mechatron skips them when serializing; we set them
          explicitly because every yapcad mate carries semantic context
          worth preserving.
    """
    parts: List[Dict[str, Any]] = []
    for part_name in _ordered_part_names(assembly):
        part_def = assembly.parts[part_name]
        part_dict: Dict[str, Any] = {
            "id": str(part_def.name),
            "name": str(part_def.name),
            "material": _material_to_mechatron(getattr(part_def, "material", None)),
            "process": "FDM",
            "datums": [_datum_to_mechatron(d) for d in part_def.datums.values()],
            "tags": [],
        }
        # ``@meta(assembly.geoms=[...])`` entries are stashed on the
        # PartDefinition as a side-channel ``geoms`` attribute by
        # ``add_part``. Convert each entry to Mechatron's PartGeom shape
        # (and from mm to metres) here. Missing or empty -> skip; the
        # Mechatron validator will warn about "no visual geometry" but
        # the kinematic chain still works.
        raw_geoms = getattr(part_def, "geoms", None)
        if isinstance(raw_geoms, list) and raw_geoms:
            converted = [g for g in (_geom_to_mechatron(e) for e in raw_geoms) if g is not None]
            if converted:
                part_dict["geoms"] = converted
        parts.append(part_dict)

    # Group mates by (parent_part, child_part) so multiple yapcad
    # mates describing the same physical interface collapse into one
    # Mechatron Interface. Preserves first-occurrence order.
    mate_groups: List[tuple] = []  # list of (parent, child, [mates])
    group_index: Dict[tuple, int] = {}
    for mate in assembly.mates:
        key = (str(mate.part_a), str(mate.part_b))
        if key in group_index:
            mate_groups[group_index[key]][2].append(mate)
        else:
            group_index[key] = len(mate_groups)
            mate_groups.append((key[0], key[1], [mate]))

    interfaces: List[Dict[str, Any]] = [
        _interface_for_mate_group(idx, parent, child, group_mates,
                                  assembly=assembly)
        for idx, (parent, child, group_mates) in enumerate(mate_groups)
    ]

    return {
        "parts": parts,
        "interfaces": interfaces,
        "loop_closures": [],
        "joint_couplings": [],
        "keyframes": [],
        "generators": [],
        "design_constraints": [],
        "last_script": f"yapcad-dsl: {assembly.name}",
    }


def to_mechatron_json(assembly, *, indent: int = 2) -> str:
    """Convenience wrapper: serialize ``to_mechatron_snapshot`` to JSON."""
    return json.dumps(to_mechatron_snapshot(assembly), indent=indent, sort_keys=False)


# ---------------------------------------------------------------------------
# Helpers (internal)
# ---------------------------------------------------------------------------

def _ordered_part_names(assembly) -> Iterable[str]:
    """Return part names in their original insertion order.

    Python 3.7+ dicts preserve insertion order, so iterating
    ``assembly.parts`` is sufficient. This helper exists to keep the
    intent explicit at the call site.
    """
    return list(assembly.parts.keys())


# ---------------------------------------------------------------------------
# Geom (visual primitive) export
# ---------------------------------------------------------------------------

# yapcad @meta authors typically work in millimetres; Mechatron's PartGeom
# uses metres (MuJoCo SI convention). Conversion factor.
_MM_TO_M = 1.0 / 1000.0


# Accepted yapCAD geom kinds (case-insensitive) -> Mechatron PartGeom variant.
_GEOM_KIND_MAP: Dict[str, str] = {
    "box":      "Box",
    "cylinder": "Cylinder",
    "sphere":   "Sphere",
    "capsule":  "Capsule",
    "mesh":     "Mesh",
}


def _geom_to_mechatron(entry: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """Convert one yapcad ``assembly.geoms`` entry into a Mechatron PartGeom.

    Author-facing entries use millimetre suffixes (``radius_mm``,
    ``half_height_mm``, ``size_mm``, ``pos_mm``) for readability; this
    helper converts them to metres for the Mechatron-side ``PartGeom``
    enum. Optional fields (``pos_mm``, ``rgba``) are dropped from the
    output when absent.

    Returns ``None`` for entries with no recognised ``type`` so callers
    can skip them gracefully.
    """
    raw_kind = str(entry.get("type", "")).lower()
    variant = _GEOM_KIND_MAP.get(raw_kind)
    if variant is None:
        return None
    result: Dict[str, Any] = {"type": variant}
    if variant == "Box":
        size_mm = entry.get("size_mm") or entry.get("size")
        if size_mm is not None:
            result["size"] = [float(s) * _MM_TO_M for s in size_mm]
    elif variant == "Cylinder":
        if "radius_mm" in entry:
            result["radius"] = float(entry["radius_mm"]) * _MM_TO_M
        elif "radius" in entry:
            result["radius"] = float(entry["radius"])  # already in metres
        if "half_height_mm" in entry:
            result["half_height"] = float(entry["half_height_mm"]) * _MM_TO_M
        elif "half_height" in entry:
            result["half_height"] = float(entry["half_height"])
    elif variant == "Sphere":
        if "radius_mm" in entry:
            result["radius"] = float(entry["radius_mm"]) * _MM_TO_M
        elif "radius" in entry:
            result["radius"] = float(entry["radius"])
    elif variant == "Capsule":
        if "radius_mm" in entry:
            result["radius"] = float(entry["radius_mm"]) * _MM_TO_M
        if "half_height_mm" in entry:
            result["half_height"] = float(entry["half_height_mm"]) * _MM_TO_M
    elif variant == "Mesh":
        if "path" in entry:
            result["path"] = str(entry["path"])
        if "scale" in entry:
            result["scale"] = float(entry["scale"])

    pos_mm = entry.get("pos_mm") or entry.get("pos")
    if pos_mm is not None:
        if "_mm" in next(iter(k for k in entry if k.endswith("_mm")), ""):
            result["pos"] = [float(p) * _MM_TO_M for p in pos_mm]
        else:
            result["pos"] = [float(p) for p in pos_mm]
    rgba = entry.get("rgba")
    if rgba is not None:
        result["rgba"] = [float(c) for c in rgba]
    return result


__all__ = [
    "to_mechatron_snapshot",
    "to_mechatron_json",
]
