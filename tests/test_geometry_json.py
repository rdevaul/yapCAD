import json
import base64
import hashlib
from pathlib import Path

import pytest

from yapcad.geom import point, line, arc, iscircle, isarc, catmullrom, iscatmullrom, nurbs, isnurbs
from yapcad.geom3d import poly2surfaceXY, solid
from yapcad.io.geometry_json import geometry_from_json, geometry_to_json, SCHEMA_ID
from yapcad.metadata import (
    add_tags,
    get_solid_metadata,
    get_surface_metadata,
    set_material,
    set_layer,
)


def _make_prism_solid():
    poly = [
        point(0, 0),
        point(1, 0),
        point(1, 1),
        point(0, 1),
        point(0, 0),
    ]
    surf, _ = poly2surfaceXY(poly)
    return solid([surf])


def test_geometry_json_roundtrip():
    sld = _make_prism_solid()
    meta = get_solid_metadata(sld, create=True)
    add_tags(meta, ['test'])
    set_material(meta, name='PLA')
    set_layer(meta, 'structure')
    # ensure surfaces inherit explicit layer
    for surf in sld[1]:
        set_layer(get_surface_metadata(surf, create=True), 'structure')

    doc = geometry_to_json([sld], units='mm', generator={'name': 'test', 'version': '0'})
    assert doc['schema'] == SCHEMA_ID
    assert doc['units'] == 'mm'
    assert len(doc['entities']) >= 2  # solid + surfaces

    solid_entry = next(e for e in doc['entities'] if e['type'] == 'solid')
    assert doc['schema'] == 'yapcad-geometry-json-v0.3'
    assert solid_entry['representations'] == {
        'authoritative': 'mesh',
        'mesh': {'format': 'indexed-triangle-set', 'role': 'authoritative'},
    }
    assert solid_entry['metadata']['layer'] == 'structure'
    surface_entries = [e for e in doc['entities'] if e['type'] == 'surface']
    assert surface_entries
    assert all(entry['metadata']['layer'] == 'structure' for entry in surface_entries)

    # JSON encode/decode sanity
    decoded = json.loads(json.dumps(doc))
    solids = geometry_from_json(decoded)
    assert len(solids) == 1
    roundtripped = solids[0]

    round_meta = get_solid_metadata(roundtripped, create=False)
    assert round_meta['material']['name'] == 'PLA'
    assert 'test' in round_meta['tags']
    assert round_meta['layer'] == 'structure'

    surfaces = roundtripped[1]
    assert surfaces and len(surfaces[0][1]) > 0
    surf_meta = get_surface_metadata(surfaces[0], create=False)
    assert surf_meta['schema'] == 'metadata-namespace-v1.1'
    assert surf_meta['layer'] == 'structure'


def test_v02_promotes_legacy_brep_metadata_to_explicit_representation():
    sld = _make_prism_solid()
    payload = b'DBRep_DrawableShape\nCASCADE Topology V3\n'
    encoded = base64.b64encode(payload).decode('ascii')
    meta = get_solid_metadata(sld, create=True)
    meta['brep'] = {'encoding': 'brep-ascii-base64', 'data': encoded}

    doc = geometry_to_json([sld])
    entry = next(item for item in doc['entities'] if item['type'] == 'solid')
    brep = entry['representations']['brep']

    assert entry['representations']['authoritative'] == 'brep'
    assert entry['representations']['mesh']['role'] == 'preview'
    assert brep['format'] == 'opencascade-brep'
    assert brep['encoding'] == 'base64'
    assert brep['payload'] == encoded
    assert brep['hash'] == 'sha256:' + hashlib.sha256(payload).hexdigest()
    assert brep['kernel']['name'] == 'OpenCASCADE'
    assert 'brep' not in entry['metadata']
    assert meta['brep']['data'] == encoded  # serialization does not mutate source


def test_v02_rejects_corrupted_brep_payload_before_kernel_loading():
    sld = _make_prism_solid()
    payload = base64.b64encode(b'not-a-real-brep').decode('ascii')
    get_solid_metadata(sld, create=True)['brep'] = {
        'encoding': 'brep-ascii-base64', 'data': payload,
    }
    doc = geometry_to_json([sld])
    entry = next(item for item in doc['entities'] if item['type'] == 'solid')
    entry['representations']['brep']['payload'] = base64.b64encode(
        b'tampered'
    ).decode('ascii')

    with pytest.raises(ValueError, match='BREP payload hash mismatch'):
        geometry_from_json(doc)


def test_v02_rejects_non_opencascade_brep_kernel():
    sld = _make_prism_solid()
    payload = base64.b64encode(b'kernel-specific-payload').decode('ascii')
    get_solid_metadata(sld, create=True)['brep'] = {
        'encoding': 'brep-ascii-base64',
        'data': payload,
        'kernel': {'name': 'another-kernel', 'version': '1.0'},
    }

    with pytest.raises(ValueError, match='unsupported BREP kernel'):
        geometry_to_json([sld])


def test_v01_mesh_documents_remain_loadable():
    doc = geometry_to_json([_make_prism_solid()])
    doc['schema'] = 'yapcad-geometry-json-v0.1'
    for entry in doc['entities']:
        entry.pop('representations', None)

    restored = geometry_from_json(doc)
    assert len(restored) == 1


def test_v02_solid_requires_representation_contract():
    doc = geometry_to_json([_make_prism_solid()])
    entry = next(item for item in doc['entities'] if item['type'] == 'solid')
    del entry['representations']
    with pytest.raises(ValueError, match='missing representations'):
        geometry_from_json(doc)


def _published_schema(version):
    path = (
        Path(__file__).resolve().parents[1]
        / 'docs' / 'schemas' / f'yapcad-geometry-json-{version}.schema.json'
    )
    return json.loads(path.read_text(encoding='utf-8'))


def test_document_matches_published_machine_readable_schema():
    jsonschema = pytest.importorskip('jsonschema')
    document = geometry_to_json([_make_prism_solid()], units='mm')
    jsonschema.Draft202012Validator(_published_schema('v0.3')).validate(document)


def test_v02_documents_remain_loadable():
    """v0.3 only extended the authoritative enum, so a mesh-authoritative
    v0.2 document is still valid input."""
    document = geometry_to_json([_make_prism_solid()], units='mm')
    document['schema'] = 'yapcad-geometry-json-v0.2'
    restored = geometry_from_json(document)
    assert len(restored) == 1


def test_v02_documents_may_not_claim_sdf_authority():
    """The reason the schema id had to move: 'sdf' is not in v0.2's enum."""
    document = geometry_to_json([_make_prism_solid()], units='mm')
    document['schema'] = 'yapcad-geometry-json-v0.2'
    entry = next(e for e in document['entities'] if e['type'] == 'solid')
    entry['representations']['authoritative'] = 'sdf'
    with pytest.raises(ValueError, match='invalid authoritative'):
        geometry_from_json(document)


def test_sketch_primitives_roundtrip():
    geom = [
        line(point(0, 0), point(1, 0)),
        arc(point(0, 0), 5),
        arc(point(2, 2), 3, 0, 180),
        catmullrom([point(0, 0), point(1, 2), point(2, 0)]),
        nurbs([point(0, 0), point(1, 2), point(3, 3), point(4, 0)], degree=3),
    ]
    doc = geometry_to_json([{'geometry': geom, 'metadata': {'layer': 'sketch'}}])
    sketch_entry = next(e for e in doc['entities'] if e['type'] == 'sketch')
    primitives = sketch_entry.get('primitives', [])
    kinds = {prim['kind'] for prim in primitives}
    assert {'line', 'circle', 'arc'}.issubset(kinds)

    roundtrip = geometry_from_json(doc)
    assert len(roundtrip) == 1
    returned_geom = roundtrip[0]
    assert any(iscircle(entity) for entity in returned_geom)
    assert any(isarc(entity) and not iscircle(entity) for entity in returned_geom)
    assert any(iscatmullrom(entity) for entity in returned_geom)
    assert any(isnurbs(entity) for entity in returned_geom)
