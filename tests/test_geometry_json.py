import json

from yapcad.geom import point
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
    assert surf_meta['schema'] == 'metadata-namespace-v0.1'
    assert surf_meta['layer'] == 'structure'
