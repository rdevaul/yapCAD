"""Accessors and serialization for the yapCAD solid ``construction`` slot.

A yapCAD solid is ``['solid', surfaces, material, construction]`` with an
optional metadata dictionary appended (see :mod:`yapcad.geom3d`).  The
``construction`` slot at index 3 records *how the solid came to be*: the
generating procedure call, the boolean operation that produced it, or — in
future — the signed-distance-function tree that defines it.

Producers throughout the codebase have written construction records for some
time (``geom3d_util`` emits ``['procedure', call]``, the boolean engines emit
``['boolean', operation]``), but the records were dropped during
serialization.  This module gives the slot a defined shape, accessors that
mirror :mod:`yapcad.metadata`, and a JSON round-trip so that provenance
survives the ``.ycpkg`` boundary.

Record format
-------------

A construction record is a list whose first element is a non-empty string
*kind* tag, followed by kind-specific payload elements::

    ['procedure', 'yapcad.geom3d_util.makeRevolutionSolid(contour,0,10,32)']
    ['boolean', 'union']
    ['boolean', 'trimesh:difference']

An empty list means "no construction information recorded", which is always
legal.

Payload elements must be JSON-representable.  :func:`normalize_construction`
enforces this on write, coercing tuples to lists and replacing values it
cannot represent with their ``repr()`` — provenance is advisory, and export
should never fail because of it.
"""

from yapcad.geom3d import issolid

# Slot index of the construction record within a solid.
CONSTRUCTION_INDEX = 3

# Kind tags currently emitted by yapCAD.
PROCEDURE = "procedure"
BOOLEAN = "boolean"
#: Reserved for the signed-distance-function tree (see docs/SDF-DESIGN.md).
SDF = "sdf"

KNOWN_KINDS = frozenset({PROCEDURE, BOOLEAN, SDF})


def _sanitize(value):
    """Recursively coerce ``value`` into something ``json.dumps`` accepts."""
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    if isinstance(value, (list, tuple)):
        return [_sanitize(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _sanitize(item) for key, item in value.items()}
    return repr(value)


def is_construction(value):
    """True if ``value`` is a well-formed construction record.

    The empty list is well-formed: it denotes the absence of construction
    information.
    """
    if not isinstance(value, (list, tuple)):
        return False
    if len(value) == 0:
        return True
    kind = value[0]
    return isinstance(kind, str) and bool(kind)


def construction_kind(record):
    """Return the kind tag of ``record``, or ``None`` if it carries none."""
    if not is_construction(record) or len(record) == 0:
        return None
    return record[0]


def normalize_construction(value):
    """Return ``value`` as a JSON-representable construction record.

    ``None`` and empty sequences normalize to ``[]``.  Raises ``ValueError``
    for a non-sequence, or for a sequence whose first element is not a
    non-empty string.
    """
    if value is None:
        return []
    if not isinstance(value, (list, tuple)):
        raise ValueError(
            f"construction record must be a list, got {type(value).__name__}"
        )
    items = list(value)
    if not items:
        return []
    kind = items[0]
    if not isinstance(kind, str) or not kind:
        raise ValueError(
            "construction record must begin with a non-empty string kind tag"
        )
    return [kind] + [_sanitize(item) for item in items[1:]]


def get_construction(sld):
    """Return the construction record attached to solid ``sld``.

    Returns the live list, not a copy, mirroring
    :func:`yapcad.metadata.get_solid_metadata`.  Returns ``[]`` when the solid
    carries no construction information.
    """
    if not issolid(sld):
        raise ValueError("bad solid passed to get_construction")
    if len(sld) <= CONSTRUCTION_INDEX:
        return []
    record = sld[CONSTRUCTION_INDEX]
    return record if isinstance(record, list) else []


def set_construction(sld, record):
    """Attach ``record`` to solid ``sld`` in place and return the solid."""
    if not issolid(sld):
        raise ValueError("bad solid passed to set_construction")
    if len(sld) <= CONSTRUCTION_INDEX:
        raise ValueError("malformed solid: no construction slot")
    sld[CONSTRUCTION_INDEX] = normalize_construction(record)
    return sld


def procedure(call, params=None):
    """Build a ``procedure`` construction record.

    :param call: the generating call, as a string.
    :param params: optional dictionary of parameters, recorded alongside.
    """
    if not isinstance(call, str) or not call:
        raise ValueError("procedure call must be a non-empty string")
    record = [PROCEDURE, call]
    if params:
        record.append(dict(params))
    return normalize_construction(record)


def boolean_op(operation, engine=None):
    """Build a ``boolean`` construction record.

    Matches the shape the boolean engines already emit: the operation alone,
    or ``'<engine>:<operation>'`` when an engine is named.
    """
    if not isinstance(operation, str) or not operation:
        raise ValueError("boolean operation must be a non-empty string")
    label = f"{engine}:{operation}" if engine else operation
    return [BOOLEAN, label]


def construction_to_json(sld):
    """Return the JSON form of ``sld``'s construction, or ``None`` if empty.

    Returning ``None`` lets serializers omit the key entirely, so documents
    for solids without provenance are unchanged.
    """
    record = normalize_construction(get_construction(sld))
    return record or None


def construction_from_json(value):
    """Rehydrate a construction record from a geometry JSON document.

    A missing key (``None``) yields ``[]``.  A present-but-malformed value
    raises ``ValueError``, consistent with how ``geometry_from_json`` treats
    the other structural fields of a solid.
    """
    if value is None:
        return []
    return normalize_construction(value)


__all__ = [
    "CONSTRUCTION_INDEX",
    "PROCEDURE",
    "BOOLEAN",
    "SDF",
    "KNOWN_KINDS",
    "is_construction",
    "construction_kind",
    "normalize_construction",
    "get_construction",
    "set_construction",
    "procedure",
    "boolean_op",
    "construction_to_json",
    "construction_from_json",
]
