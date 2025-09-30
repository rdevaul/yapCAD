import uuid

from yapcad.geom import point
from yapcad.geom3d import poly2surfaceXY, solid, surface
from yapcad.metadata import (
    ensure_solid_id,
    ensure_surface_id,
    get_solid_metadata,
    get_surface_metadata,
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

    assert get_surface_metadata(surf) == {}

    ensure_surface_id(surf)
    meta = get_surface_metadata(surf)
    assert 'id' in meta

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

    assert get_solid_metadata(sld) == {}

    ensure_solid_id(sld)
    meta = get_solid_metadata(sld)
    assert 'id' in meta

    set_solid_context(sld, {'solver': 'dummy'})
    assert get_solid_metadata(sld)['context'] == {'solver': 'dummy'}

    custom = {'id': str(uuid.uuid4())}
    set_solid_metadata(sld, custom)
    assert get_solid_metadata(sld) is custom
