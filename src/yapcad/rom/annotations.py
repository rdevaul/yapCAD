"""Parse ``@meta(kinematic.*)`` annotations from DSL commands.

Extends the ``meta_apply`` namespace routing to handle the ``kinematic``
namespace, producing a flat ``kinematic`` dict on each command's metadata
that the RomTree builder can consume.

Supported kinematic fields
==========================
.. list-table::
   :header-rows: 1

   * - Field
     - Type
     - Description
   * - ``joint_kind``
     - ``"fixed"`` | ``"revolute"`` | ``"prismatic"``
     - Type of joint
   * - ``parent``
     - string
     - Parent part name
   * - ``parent_frame``
     - string
     - Frame on parent (optional, default "ORIGIN")
   * - ``axis``
     - ``[x, y, z]``
     - Motion axis (list of 3 floats)
   * - ``limits``
     - ``[min, max]``
     - Angle/translation range
   * - ``zero_pose``
     - float
     - Default joint value
   * - ``description``
     - string
     - Human-readable part label
   * - ``printable``
     - bool
     - Is a 3D-printable part

Example::

    @meta(kinematic.joint_kind="revolute")
    @meta(kinematic.parent="PELVIS")
    @meta(kinematic.parent_frame="HIP_YAW_MOUNT")
    @meta(kinematic.axis=[0, 0, 1])
    @meta(kinematic.limits=[-90, 90])
    command HIP_YAW_JOINT() -> solid:
        emit translate(srv, 0, ho*s, 0)

The DSL ``@meta`` decorator already collects all ``kinematic.*`` keys into
``FunctionDef.meta_hint`` as flat dotted keys.  This module routes them into
a ``kinematic`` sub-dict in the solid's metadata, ready for
:class:`RomTree.from_interpreter_context`.
"""

from __future__ import annotations

from typing import Any, Dict

from yapcad.metadata import get_solid_metadata, get_assembly_metadata, _ensure_root

# ---------------------------------------------------------------------------
# Recognised kinematic sub-fields
# ---------------------------------------------------------------------------
_KINEMATIC_SCALAR = {"joint_kind", "parent", "parent_frame", "description"}

_KINEMATIC_LISTS = {
    "limits": (float, 2),       # [min, max]
    "axis": (float, 3),         # (x, y, z)
}

_KINEMATIC_BOOL = {"printable"}

_KINEMATIC_ALL = _KINEMATIC_SCALAR | set(_KINEMATIC_LISTS.keys()) | _KINEMATIC_BOOL


def apply_kinematic_hint(meta: Dict[str, Any], field: str, value: Any) -> None:
    """Route a single ``kinematic.<field>`` hint into the metadata dict.

    Called by ``meta_apply`` for dotted keys whose namespace is ``"kinematic"``.

    Args:
        meta:  Mutable metadata dict for the solid.
        field: The part after ``kinematic.`` (e.g., ``"joint_kind"``).
        value: The value from the DSL.
    """
    section = meta.setdefault("kinematic", {})

    if field in _KINEMATIC_SCALAR:
        section[field] = str(value)

    elif field == "zero_pose":
        section[field] = float(value)

    elif field in _KINEMATIC_LISTS:
        expected_len = _KINEMATIC_LISTS[field][1]
        if not isinstance(value, (list, tuple)):
            raise TypeError(
                f"kinematic.{field} must be a list/tuple of {expected_len} "
                f"values, got {type(value).__name__}"
            )
        items = [float(v) for v in value]
        if len(items) != expected_len:
            raise ValueError(
                f"kinematic.{field} needs exactly {expected_len} values, "
                f"got {len(items)}"
            )
        section[field] = items

    elif field in _KINEMATIC_BOOL:
        if isinstance(value, str):
            value = value.lower() in ("true", "1", "yes")
        section[field] = bool(value)

    else:
        # Unknown kinematic field — store verbatim for forward compat
        section[field] = value


# ---------------------------------------------------------------------------
# Register with the yapcad.meta_apply dispatch
# ---------------------------------------------------------------------------
# Importing this module as a side-effect registers the "kinematic" namespace.

_META_APPLY_PATCHED = False


def _patch_meta_apply() -> None:
    """Monkey-patch ``meta_apply`` to route ``kinematic.*`` keys.

    This is a lightweight patching strategy: we wrap the existing
    ``_apply_assembly``-style dispatch in ``meta_apply`` to check for
    ``kinematic.*`` first.  If the upstream ``meta_apply`` ever gains
    a native ``kinematic`` namespace handler this monkey-patch can be
    removed.
    """
    global _META_APPLY_PATCHED
    if _META_APPLY_PATCHED:
        return

    import yapcad.dsl.meta_apply as ma

    original = ma.apply_meta_hint_to_raw

    def _wrapped(meta, meta_hint):
        if meta_hint is None:
            return original(meta, meta_hint)
        kin_keys = set()
        for key in meta_hint:
            if isinstance(key, str) and key.startswith("kinematic."):
                field = key.split(".", 1)[1]
                apply_kinematic_hint(meta, field, meta_hint[key])
                kin_keys.add(key)
        # All other namespaces pass through to existing handler
        filtered = {k: v for k, v in meta_hint.items() if k not in kin_keys}
        if filtered:
            return original(meta, filtered)

    ma.apply_meta_hint_to_raw = _wrapped
    _META_APPLY_PATCHED = True


# Apply on import
_patch_meta_apply()


__all__ = ["apply_kinematic_hint"]