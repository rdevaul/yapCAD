"""Tests for the solid ``construction`` slot and its JSON round-trip.

Phase 0 of docs/SDF-DESIGN.md: construction records must survive
serialization so that provenance is not lost at the ``.ycpkg`` boundary.
"""

from __future__ import annotations

import json

import pytest

from yapcad.construction import (
    BOOLEAN,
    PROCEDURE,
    boolean_op,
    construction_from_json,
    construction_kind,
    construction_to_json,
    get_construction,
    is_construction,
    normalize_construction,
    procedure,
    set_construction,
)
from yapcad.geom3d import solid
from yapcad.geom3d_util import makeRevolutionSolid, prism, sphere
from yapcad.io.geometry_json import geometry_from_json, geometry_to_json


def _roundtrip(entities):
    """Serialize and restore through a real JSON encode/decode cycle."""
    document = json.loads(json.dumps(geometry_to_json(entities, units="mm")))
    return document, geometry_from_json(document)


def _solid_entry(document):
    return next(e for e in document["entities"] if e["type"] == "solid")


# ---------------------------------------------------------------------------
# Record shape
# ---------------------------------------------------------------------------

def test_empty_record_is_well_formed():
    assert is_construction([])
    assert construction_kind([]) is None
    assert normalize_construction([]) == []
    assert normalize_construction(None) == []


@pytest.mark.parametrize("value", [None, 42, "procedure", {"kind": "procedure"}])
def test_non_sequence_is_not_a_record(value):
    assert not is_construction(value)


@pytest.mark.parametrize("bad", [[42], [""], [None, "x"], [["procedure"]]])
def test_kind_tag_must_be_a_non_empty_string(bad):
    assert not is_construction(bad)
    with pytest.raises(ValueError):
        normalize_construction(bad)


def test_normalize_rejects_non_sequences():
    with pytest.raises(ValueError):
        normalize_construction({"kind": "procedure"})


def test_procedure_builder():
    record = procedure("yapcad.geom3d_util.prism(1,2,3)")
    assert record == [PROCEDURE, "yapcad.geom3d_util.prism(1,2,3)"]
    assert construction_kind(record) == PROCEDURE


def test_procedure_builder_with_params():
    record = procedure("prism", {"length": 1, "width": 2})
    assert record == [PROCEDURE, "prism", {"length": 1, "width": 2}]


def test_procedure_builder_rejects_empty_call():
    with pytest.raises(ValueError):
        procedure("")


def test_boolean_builder_matches_engine_output():
    # The shapes the native and trimesh engines already emit.
    assert boolean_op("union") == [BOOLEAN, "union"]
    assert boolean_op("difference", engine="trimesh") == [BOOLEAN, "trimesh:difference"]


def test_boolean_builder_rejects_empty_operation():
    with pytest.raises(ValueError):
        boolean_op("")


# ---------------------------------------------------------------------------
# Sanitization
# ---------------------------------------------------------------------------

def test_tuples_are_coerced_to_lists():
    record = normalize_construction(["procedure", ("a", ("b", 1))])
    assert record == ["procedure", ["a", ["b", 1]]]
    json.dumps(record)  # must not raise


def test_dict_keys_are_stringified():
    record = normalize_construction(["procedure", {1: "a", None: "b"}])
    assert record == ["procedure", {"1": "a", "None": "b"}]
    json.dumps(record)


def test_unrepresentable_payload_falls_back_to_repr():
    class Opaque:
        def __repr__(self):
            return "<opaque>"

    record = normalize_construction(["procedure", Opaque()])
    assert record == ["procedure", "<opaque>"]
    json.dumps(record)


def test_booleans_survive_sanitization():
    # bool is a subclass of int; make sure it is not silently widened.
    record = normalize_construction(["procedure", True, False, 1, 0])
    assert record == ["procedure", True, False, 1, 0]
    assert record[1] is True and record[2] is False


# ---------------------------------------------------------------------------
# Accessors
# ---------------------------------------------------------------------------

def _bare_solid():
    """A solid with no construction record (prism() records its own)."""
    return solid([prism(2, 2, 2)[1][0]])


def test_get_construction_on_solid_without_one():
    assert get_construction(_bare_solid()) == []


def test_set_and_get_construction():
    sld = prism(2, 2, 2)
    set_construction(sld, procedure("prism(2,2,2)"))
    assert get_construction(sld) == [PROCEDURE, "prism(2,2,2)"]
    # Written into the documented slot, not appended.
    assert sld[3] == [PROCEDURE, "prism(2,2,2)"]
    assert len(sld) in (4, 5)


def test_set_construction_normalizes():
    sld = prism(1, 1, 1)
    set_construction(sld, ("procedure", ("a", "b")))
    assert get_construction(sld) == ["procedure", ["a", "b"]]


def test_set_construction_rejects_bad_record():
    sld = prism(1, 1, 1)
    with pytest.raises(ValueError):
        set_construction(sld, [42])


@pytest.mark.parametrize("fn", [get_construction, lambda s: set_construction(s, [])])
def test_accessors_reject_non_solids(fn):
    with pytest.raises(ValueError):
        fn(["not", "a", "solid"])


def test_solid_constructor_files_construction_in_slot_three():
    """``solid(surfaces, [], construction)`` is the near-universal idiom.

    An explicitly empty material list must consume its positional slot, or
    the construction record lands under material instead.
    """
    record = procedure("test()")
    sld = solid([], [], record)
    assert sld[2] == []           # material
    assert sld[3] == record       # construction
    assert get_construction(sld) == record


def test_solid_constructor_still_accepts_material():
    record = procedure("test()")
    sld = solid([], ["steel"], record)
    assert sld[2] == ["steel"]
    assert get_construction(sld) == record


def test_solid_constructor_rejects_a_third_list():
    with pytest.raises(ValueError):
        solid([], [], ["procedure", "a"], ["extra"])


@pytest.mark.parametrize(
    "factory, expected",
    [
        (lambda: sphere(4.0), "sphere"),
        (lambda: makeRevolutionSolid(lambda z: 5.0, 0.0, 5.0, 8), "makeRevolutionSolid"),
    ],
)
def test_generated_solids_carry_a_procedure_record(factory, expected):
    """geom3d_util producers write these; verify they reach the right slot."""
    sld = factory()
    record = get_construction(sld)
    assert construction_kind(record) == PROCEDURE
    assert expected in record[1]
    # And material stays empty rather than absorbing the record.
    assert sld[2] == []


# ---------------------------------------------------------------------------
# JSON round-trip
# ---------------------------------------------------------------------------

def test_construction_key_omitted_when_absent():
    sld = _bare_solid()
    assert get_construction(sld) == []
    document, _ = _roundtrip([sld])
    assert "construction" not in _solid_entry(document)


def test_construction_to_json_returns_none_when_empty():
    assert construction_to_json(_bare_solid()) is None


def test_construction_survives_roundtrip():
    sld = prism(2, 2, 2)
    set_construction(sld, procedure("prism(2,2,2)", {"length": 2}))

    document, restored = _roundtrip([sld])
    assert _solid_entry(document)["construction"] == [
        PROCEDURE,
        "prism(2,2,2)",
        {"length": 2},
    ]
    assert get_construction(restored[0]) == get_construction(sld)


def test_boolean_record_survives_roundtrip():
    sld = prism(2, 2, 2)
    set_construction(sld, boolean_op("union"))
    _, restored = _roundtrip([sld])
    assert get_construction(restored[0]) == [BOOLEAN, "union"]


def test_restored_construction_lands_in_the_documented_slot():
    sld = prism(2, 2, 2)
    set_construction(sld, procedure("prism(2,2,2)"))
    _, restored = _roundtrip([sld])
    assert restored[0][3] == [PROCEDURE, "prism(2,2,2)"]


def test_legacy_v01_documents_restore_with_empty_construction():
    sld = prism(2, 2, 2)
    set_construction(sld, procedure("prism(2,2,2)"))
    document = json.loads(json.dumps(geometry_to_json([sld], units="mm")))

    # Simulate a v0.1 producer: no construction key, legacy schema id.
    document["schema"] = "yapcad-geometry-json-v0.1"
    entry = _solid_entry(document)
    entry.pop("construction")
    entry.pop("representations", None)

    restored = geometry_from_json(document)
    assert get_construction(restored[0]) == []


def test_malformed_construction_in_document_is_rejected():
    sld = prism(2, 2, 2)
    set_construction(sld, procedure("prism(2,2,2)"))
    document = json.loads(json.dumps(geometry_to_json([sld], units="mm")))
    _solid_entry(document)["construction"] = [42, "nonsense"]

    with pytest.raises(ValueError, match="malformed construction record"):
        geometry_from_json(document)


def test_construction_from_json_defaults_to_empty():
    assert construction_from_json(None) == []


def test_generated_procedure_record_survives_roundtrip_untouched():
    """The common case: nobody sets construction, the producer already did."""
    sld = prism(2, 2, 2)
    original = list(get_construction(sld))
    assert original and original[0] == PROCEDURE

    document, restored = _roundtrip([sld])
    assert _solid_entry(document)["construction"] == original
    assert get_construction(restored[0]) == original


def test_generated_solid_emits_no_spurious_voids():
    """Before the slot fix, a construction record in the material slot was
    serialized as a list of empty void loops."""
    document, _ = _roundtrip([prism(2, 2, 2)])
    assert _solid_entry(document)["voids"] == []


def test_multiple_solids_keep_distinct_construction():
    a = prism(1, 1, 1)
    b = prism(2, 2, 2)
    set_construction(a, procedure("a()"))
    set_construction(b, boolean_op("difference", engine="native"))

    _, restored = _roundtrip([a, b])
    records = sorted(get_construction(s) for s in restored)
    assert records == [
        [BOOLEAN, "native:difference"],
        [PROCEDURE, "a()"],
    ]
