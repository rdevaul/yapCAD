import json

from pathlib import Path

from yapcad.geom import point
from yapcad.geom import line
from yapcad.geom3d import poly2surfaceXY, solid
from yapcad.metadata import add_tags, get_solid_metadata, get_surface_metadata, set_layer
from yapcad.package import (
    PackageManifest,
    create_package_from_entities,
    load_geometry,
)
from yapcad.package.validator import validate_package


def _make_solid():
    poly = [
        point(0, 0),
        point(1, 0),
        point(1, 1),
        point(0, 1),
        point(0, 0),
    ]
    surf, _ = poly2surfaceXY(poly)
    sld = solid([surf])
    meta = get_solid_metadata(sld, create=True)
    add_tags(meta, ["prototype"])
    set_layer(meta, "structure")
    for surface in sld[1]:
        set_layer(get_surface_metadata(surface, create=True), "structure")
    return sld


def test_create_package_from_entities(tmp_path: Path):
    sld = _make_solid()
    pkg_root = tmp_path / "demo.ycpkg"
    manifest = create_package_from_entities(
        [sld],
        pkg_root,
        name="Demo Part",
        version="0.1.0",
        author="tester",
        description="Example packaged part",
    )

    assert manifest.manifest_path.exists()
    data = manifest.data
    assert data["schema"] == "ycpkg-spec-v0.1"
    assert data["name"] == "Demo Part"
    assert data["tags"] == ["prototype"]

    primary = data["geometry"]["primary"]
    geometry_path = pkg_root / primary["path"]
    assert geometry_path.exists()
    assert primary["hash"].startswith("sha256:")

    reload_manifest = PackageManifest.load(pkg_root)
    reload_manifest.recompute_hashes()
    reload_manifest.save()

    entities = load_geometry(reload_manifest)
    assert len(entities) == 1
    meta = get_solid_metadata(entities[0], create=False)
    assert "prototype" in meta.get("tags", [])
    assert meta.get("layer") == "structure"

    # verify geometry JSON layer annotation
    with (pkg_root / primary["path"]).open() as fp:
        doc = json.load(fp)
    solid_entry = next(e for e in doc["entities"] if e["type"] == "solid")
    assert solid_entry["metadata"].get("layer") == "structure"
    assert all(entry["metadata"].get("layer") == "structure" for entry in doc["entities"] if entry["type"] == "surface")


def test_validate_package(tmp_path: Path):
    sld = _make_solid()
    pkg_root = tmp_path / "demo.ycpkg"
    manifest = create_package_from_entities(
        [sld],
        pkg_root,
        name="Validator Part",
        version="0.1.0",
    )
    manifest.recompute_hashes()
    manifest.save()
    ok, messages = validate_package(pkg_root)
    assert ok, f"Validation failed: {messages}"


def test_package_from_geomlist(tmp_path: Path):
    glist = [
        line(point(0, 0), point(10, 0)),
        line(point(10, 0), point(10, 5)),
        line(point(10, 5), point(0, 5)),
        line(point(0, 5), point(0, 0)),
    ]
    pkg_root = tmp_path / "sketch.ycpkg"
    manifest = create_package_from_entities(
        [{'geometry': glist, 'metadata': {'layer': 'sketch'}}],
        pkg_root,
        name="Sketch",
        version="0.1.0",
        description="Rectangular outline",
    )
    manifest.recompute_hashes()
    manifest.save()
    ok, messages = validate_package(pkg_root)
    assert ok, f"Sketch package validation failed: {messages}"

    with (pkg_root / manifest.data['geometry']['primary']['path']).open() as fp:
        doc = json.load(fp)
    assert all(entry['metadata'].get('layer') == 'sketch' for entry in doc['entities'])
