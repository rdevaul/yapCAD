import uuid

from yapcad.geom import point
from yapcad.geom3d import poly2surfaceXY, solid, surface
from yapcad.metadata import (
    add_analysis_record,
    add_compliance,
    add_design_history_entry,
    add_tags,
    ensure_solid_id,
    ensure_surface_id,
    get_solid_metadata,
    get_surface_metadata,
    set_envelope_constraint,
    set_manufacturing,
    set_mass_constraint,
    set_material,
    set_layer,
    set_solid_context,
    set_solid_metadata,
    set_surface_metadata,
    set_surface_origin,
    set_surface_units,
)


def test_surface_metadata_roundtrip():
    poly = [
        point(0, 0),
        point(1, 0),
        point(0, 1),
        point(0, 0),
    ]
    surf, _ = poly2surfaceXY(poly)

    meta = get_surface_metadata(surf)
    assert meta['schema'] == 'metadata-namespace-v0.1'
    assert meta['tags'] == []
    assert meta['layer'] == 'default'

    meta = get_surface_metadata(surf, create=True)
    assert meta['schema'] == 'metadata-namespace-v0.1'
    assert 'entityId' in meta
    assert meta['layer'] == 'default'
    assert isinstance(meta['tags'], list)

    set_surface_units(surf, 'mm')
    set_surface_origin(surf, 'parametric:test')
    meta = get_surface_metadata(surf)
    assert meta['units'] == 'mm'
    assert meta['origin'] == 'parametric:test'

    custom = {'foo': 'bar'}
    set_surface_metadata(surf, custom)
    assert get_surface_metadata(surf) is custom


def test_surface_constructor_with_metadata():
    poly = [
        point(0, 0),
        point(1, 0),
        point(0, 1),
        point(0, 0),
    ]
    surf, _ = poly2surfaceXY(poly, holepolys=[], minlen=0.5)
    metadata = {'id': '1234'}
    surf_with_meta = surface(surf[1], surf[2], surf[3], surf[4], surf[5], metadata)
    assert get_surface_metadata(surf_with_meta) == metadata


def test_solid_metadata_helpers():
    poly = [
        point(0, 0),
        point(1, 0),
        point(0, 1),
        point(0, 0),
    ]
    surf, _ = poly2surfaceXY(poly)
    sld = solid([surf])

    meta = get_solid_metadata(sld)
    assert meta['schema'] == 'metadata-namespace-v0.1'
    assert meta['tags'] == []
    assert meta['layer'] == 'default'

    meta = get_solid_metadata(sld, create=True)
    assert meta['schema'] == 'metadata-namespace-v0.1'
    assert 'entityId' in meta
    assert meta['layer'] == 'default'

    set_solid_context(sld, {'solver': 'dummy'})
    assert get_solid_metadata(sld)['context'] == {'solver': 'dummy'}

    custom = {'id': str(uuid.uuid4())}
    set_solid_metadata(sld, custom)
    assert get_solid_metadata(sld) is custom


def test_metadata_namespace_helpers():
    poly = [
        point(0, 0),
        point(1, 0),
        point(0, 1),
        point(0, 0),
    ]
    surf, _ = poly2surfaceXY(poly)
    sld = solid([surf])
    ensure_solid_id(sld)

    meta = get_solid_metadata(sld)
    add_tags(meta, ['prototype', 'waterjet'])
    set_material(meta, name='6061 Aluminum', grade='6061-T6', density_kg_m3=2700)
    set_manufacturing(meta, process='waterjet', instructions='Cut in one pass')
    add_design_history_entry(meta, author='assistant', source='prompt', revision='A1')
    set_mass_constraint(meta, max_kg=5.0)
    set_envelope_constraint(meta, x_mm=300, y_mm=100, z_mm=25)
    add_compliance(meta, ['ASME Y14.5-2018'])
    add_analysis_record(meta, {'solver': 'fea-lite', 'resultId': 'run-42'})
    set_layer(meta, 'structure')

    assert 'prototype' in meta['tags']
    assert meta['material']['grade'] == '6061-T6'
    assert meta['manufacturing']['process'] == 'waterjet'
    assert meta['designHistory']['iterations'][0]['revision'] == 'A1'
    assert meta['constraints']['mass']['max_kg'] == 5.0
    assert meta['constraints']['envelope']['x_mm'] == 300.0
    assert 'ASME Y14.5-2018' in meta['constraints']['compliance']
    assert meta['analysis']['simulations'][0]['solver'] == 'fea-lite'
    assert meta['layer'] == 'structure'
