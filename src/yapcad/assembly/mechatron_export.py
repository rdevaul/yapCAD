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
    "concentric": "Concentric",
    "coincident": "Coincident",
    "tangent": "Tangent",
    "angle": "Angle",
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


def _joint_for_mate(mate) -> Dict[str, Any]:
    """Build a Mechatron Joint tagged-enum dict from a yapcad Mate.

    ``Joint`` is ``#[serde(tag = "type")]`` in Mechatron, so the dict
    always has a ``"type"`` key plus variant-specific fields. We emit
    a minimal Joint here \u2014 axes and limits are left for the Mechatron
    solver to resolve from the datum geometry.
    """
    mate_kind = mate.mate_type.value if hasattr(mate.mate_type, "value") else str(mate.mate_type)
    joint_variant = _MOTION_MATE_TO_JOINT.get(mate_kind.lower(), "Fixed")
    if joint_variant in ("Revolute", "Prismatic"):
        # Mechatron requires an axis on these joints. yapCAD's mate
        # references a datum on each part; the actual axis vector comes
        # out of the solver after FK. We emit +Z as a placeholder; the
        # consumer's solver overrides this once datums are resolved.
        return {"type": joint_variant, "axis": [0.0, 0.0, 1.0]}
    return {"type": joint_variant}


def _mate_constraint_for_mate(mate) -> Optional[Dict[str, Any]]:
    """Build a Mechatron ``MateConstraint`` dict from a yapcad Mate.

    Returns ``None`` when the mate kind is a motion mate (its semantics
    are captured by the joint already) or has no Mechatron equivalent.

    Output uses Mechatron's ``#[serde(tag = "type")]`` shape:

        {"type": "Concentric", "parent_datum": "...", "child_datum": "..."}

    Note: yapCAD's mate API uses ``datum_a`` / ``datum_b``; Mechatron's
    MateConstraint uses ``parent_datum`` / ``child_datum``. The mapping
    is positional (datum_a -> parent_datum, datum_b -> child_datum) and
    happens here, not in the exporter caller.
    """
    mate_kind = mate.mate_type.value if hasattr(mate.mate_type, "value") else str(mate.mate_type)
    variant = _MATE_TO_CONSTRAINT_VARIANT.get(mate_kind.lower())
    if variant is None:
        return None
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


def _interface_for_mate(idx: int, mate) -> Dict[str, Any]:
    """Convert one yapcad Mate into a Mechatron Interface dict."""
    constraint = _mate_constraint_for_mate(mate)
    return {
        "id": f"iface_{idx}",
        "parent_part": str(mate.part_a),
        "child_part": str(mate.part_b),
        "joint": _joint_for_mate(mate),
        "mate_constraints": [constraint] if constraint is not None else [],
    }


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
        parts.append({
            "id": str(part_def.name),
            "name": str(part_def.name),
            "material": _material_to_mechatron(getattr(part_def, "material", None)),
            "process": "FDM",
            "datums": [_datum_to_mechatron(d) for d in part_def.datums.values()],
            "tags": [],
        })

    interfaces: List[Dict[str, Any]] = [
        _interface_for_mate(idx, mate) for idx, mate in enumerate(assembly.mates)
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


__all__ = [
    "to_mechatron_snapshot",
    "to_mechatron_json",
]
