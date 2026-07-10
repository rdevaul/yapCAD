import uuid

from yapcad.geom import point
from yapcad.geom3d import poly2surfaceXY, solid, surface
import pytest

from yapcad.metadata import (
    add_analysis_record,
    add_bolt_pattern,
    add_compliance,
    add_datum,
    add_design_history_entry,
    add_keepout,
    add_surface,
    add_tags,
    ensure_solid_id,
    ensure_surface_id,
    get_assembly_metadata,
    get_operation_metadata,
    get_solid_metadata,
    get_surface_metadata,
    set_assembly,
    set_envelope_constraint,
    set_manufacturing,
    set_mass_constraint,
    set_material,
    set_layer,
    set_operation,
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
    assert meta['schema'] == 'metadata-namespace-v1.1'
    assert meta['tags'] == []
    assert meta['layer'] == 'default'

    meta = get_surface_metadata(surf, create=True)
    assert meta['schema'] == 'metadata-namespace-v1.1'
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
    assert meta['schema'] == 'metadata-namespace-v1.1'
    assert meta['tags'] == []
    assert meta['layer'] == 'default'

    meta = get_solid_metadata(sld, create=True)
    assert meta['schema'] == 'metadata-namespace-v1.1'
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


# ---------------------------------------------------------------------------
# v1.1 schema + namespace helpers

def test_schema_v1_1_default():
    """Fresh metadata blocks declare the v1.1 schema string."""
    poly = [point(0, 0), point(1, 0), point(0, 1), point(0, 0)]
    surf, _ = poly2surfaceXY(poly)
    meta = get_surface_metadata(surf, create=True)
    assert meta['schema'] == 'metadata-namespace-v1.1'


def test_assembly_namespace_helpers():
    """Round-trip assembly namespace: scalars + bolt_patterns + datums + surfaces + keepouts."""
    poly = [point(0, 0), point(1, 0), point(0, 1), point(0, 0)]
    surf, _ = poly2surfaceXY(poly)
    sld = solid([surf])
    ensure_solid_id(sld)
    meta = get_solid_metadata(sld)

    set_assembly(meta, joint_kind='radial', no_cut=True)
    add_bolt_pattern(
        meta,
        id='bot_ring',
        ring='radial',
        R_mm=80.0,
        z_mm=12.5,
        count=8,
        theta0_deg=22.5,
        fastener={'designation': 'M6x1.0', 'thread': 'tap', 'head': 'SHCS'},
        tolerance_mm=0.1,
    )
    add_datum(meta, id='primary_axis', kind='axis', ring='axial', direction=[0, 0, 1])
    add_surface(meta, id='top_face', kind='mating', mate_to='intertank.bot_face',
                finish='Ra1.6', tolerance_mm=0.05)
    add_keepout(meta, id='harness_zone', kind='volume',
                shape={'primitive': 'box', 'params': {'x_mm': 50, 'y_mm': 50, 'z_mm': 30}},
                reason='MAIN_BUS_HARNESS')

    asm = get_assembly_metadata(meta)
    assert asm['joint_kind'] == 'radial'
    assert asm['no_cut'] is True
    assert len(asm['bolt_patterns']) == 1
    bp = asm['bolt_patterns'][0]
    assert bp['ring'] == 'radial'
    assert bp['R_mm'] == 80.0
    assert bp['count'] == 8
    assert bp['fastener']['thread'] == 'tap'
    assert asm['datums'][0]['direction'] == [0.0, 0.0, 1.0]
    assert asm['surfaces'][0]['mate_to'] == 'intertank.bot_face'
    assert asm['keepouts'][0]['reason'] == 'MAIN_BUS_HARNESS'


def test_assembly_namespace_validation():
    """Invalid enum values are rejected at write time."""
    meta = {}
    with pytest.raises(ValueError, match='joint_kind'):
        set_assembly(meta, joint_kind='diagonal')
    with pytest.raises(ValueError, match='bolt_pattern.ring'):
        add_bolt_pattern(meta, id='r1', ring='diagonal', R_mm=10, z_mm=0, count=4)
    with pytest.raises(ValueError, match='datum.kind'):
        add_datum(meta, id='d1', kind='bogus')
    with pytest.raises(ValueError, match='surface.kind'):
        add_surface(meta, id='s1', kind='bogus')
    with pytest.raises(ValueError, match='keepout.kind'):
        add_keepout(meta, id='k1', kind='bogus')
    with pytest.raises(ValueError, match='fastener.thread'):
        add_bolt_pattern(meta, id='r2', ring='axial', R_mm=10, z_mm=0, count=4,
                         fastener={'thread': 'bogus'})
    with pytest.raises(ValueError, match='datum.direction'):
        add_datum(meta, id='d2', kind='axis', direction=[1, 0])  # length != 3
    with pytest.raises(ValueError, match='id is required'):
        add_bolt_pattern(meta, id='', ring='axial', R_mm=10, z_mm=0, count=4)


def test_assembly_namespace_no_create_returns_empty():
    """get_assembly_metadata with create=False returns {} without mutating meta."""
    meta = {}
    asm = get_assembly_metadata(meta, create=False)
    assert asm == {}
    assert 'assembly' not in meta  # still untouched


def test_operation_namespace_helpers():
    """Round-trip operation namespace with all keys."""
    meta = {}
    set_operation(
        meta,
        kind='subtract',
        target_filter=['forward_bulkhead', 'lox_dome'],
        priority=10,
        through=True,
        consume=False,
        policy='warn',
        feature_id='lox_man_vent',
        feature_kind='vent',
        stage='pre_position',
    )
    op = get_operation_metadata(meta)
    assert op['kind'] == 'subtract'
    assert op['target_filter'] == ['forward_bulkhead', 'lox_dome']
    assert op['priority'] == 10.0
    assert op['through'] is True
    assert op['consume'] is False
    assert op['policy'] == 'warn'
    assert op['feature_id'] == 'lox_man_vent'
    assert op['feature_kind'] == 'vent'
    assert op['stage'] == 'pre_position'


def test_operation_namespace_validation():
    """Invalid enum values rejected; defaults applied where appropriate."""
    meta = {}
    with pytest.raises(ValueError, match='operation.kind'):
        set_operation(meta, kind='bogus')
    with pytest.raises(ValueError, match='operation.policy'):
        set_operation(meta, kind='subtract', policy='bogus')
    with pytest.raises(ValueError, match='operation.feature_kind'):
        set_operation(meta, kind='subtract', feature_kind='bogus')
    # minimal valid: only `kind` set
    set_operation(meta, kind='intersect')
    assert get_operation_metadata(meta)['kind'] == 'intersect'


def test_assembly_and_operation_coexist():
    """A single metadata dict may carry both assembly and operation namespaces.

    Useful for retained inserts (consume=false) that also expose mating surfaces.
    """
    meta = {}
    set_assembly(meta, joint_kind='axial')
    add_surface(meta, id='wear_face', kind='wear', finish='Ra0.8')
    set_operation(meta, kind='subtract', target_filter=['outer_skin'], consume=False)
    assert get_assembly_metadata(meta)['joint_kind'] == 'axial'
    assert get_operation_metadata(meta)['consume'] is False
