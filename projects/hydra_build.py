#!/usr/bin/env python
"""Build Hydra Robot multi-layer package from DSL.

Runs each layer command separately and combines into a single
multi-entity package with per-entity layer/material metadata.
"""
import sys, json, uuid, shutil, hashlib
from pathlib import Path
from datetime import datetime, timezone

sys.path.insert(0, 'src')

from yapcad.dsl.runtime import compile_and_run
from yapcad.geom3d import issolid

source = Path('projects/hydra_robot.dsl').read_text()
source_hash = hashlib.sha256(source.encode()).hexdigest()

# Define layers to build
layers = [
    ('BODY',       'body',      'aluminum',  '#6688aa'),
    ('CONNECTORS', 'connector', 'copper',    '#cc8844'),
    ('ARMS',       'arm',       'steel',     '#999999'),
    ('EFFECTORS',  'effector',  'plastic',   '#44aa66'),
]

results = []
for cmd, layer, material, color in layers:
    print(f"Running {cmd}...", end=' ', flush=True)
    result = compile_and_run(source, cmd, {})
    if not result.success:
        print(f"FAILED: {result.error_message}")
        sys.exit(1)
    geom = result.geometry
    meta = result.metadata
    print(f"OK (metadata: {meta})")
    results.append((cmd, layer, material, color, geom, meta))

# Now build a package manually with multi-entity geometry
target = Path('/tmp/hydra_robot_pkg')
if target.exists():
    shutil.rmtree(target)

target.mkdir(parents=True)
(target / 'geometry').mkdir()
(target / 'attachments').mkdir()
(target / 'exports').mkdir()
(target / 'metadata').mkdir()
(target / 'validation').mkdir()

# Save DSL source as attachment
(target / 'attachments' / 'source.dsl').write_text(source)

# Build geometry JSON with multiple entities
# Use yapcad's serialization
from yapcad.io.geometry_json import geometry_to_json

# Build entity list with metadata overrides for geometry_to_json
entity_inputs = []
for cmd, layer, material, color, geom, meta in results:
    entry = {
        'geometry': geom,
        'metadata': {
            'name': meta.get('name', cmd.lower()),
            'layer': layer,
            'material': material,
            'color': color,
        },
    }
    entity_inputs.append(entry)
    print(f"  Entity: {entry['metadata']['name']} (layer={layer}, material={material})")

geom_data = geometry_to_json(
    entity_inputs,
    units='mm',
    generator={'tool': 'yapCAD-DSL', 'version': '1.0.0rc1'},
)
geom_json = json.dumps(geom_data, indent=2)
(target / 'geometry' / 'primary.json').write_text(geom_json)
geom_hash = hashlib.sha256(geom_json.encode()).hexdigest()
entities = geom_data.get('entities', [])

# Build manifest
manifest = {
    'schema': 'ycpkg-spec-v0.1',
    'id': str(uuid.uuid4()),
    'name': 'Hydra Robot',
    'version': '0.1.0',
    'description': 'Modular reconfigurable robot - starfish configuration with 4 arms, 2 grippers, 1 wheel',
    'created': {
        'timestamp': datetime.now(timezone.utc).isoformat(),
    },
    'generator': {
        'tool': 'yapCAD-DSL',
        'version': '1.0.0rc1',
        'dsl': {
            'module': 'hydra_robot',
            'commands': [l[0] for l in layers],
            'source_signature': f'sha256:{source_hash}',
        },
    },
    'units': 'mm',
    'tags': ['robot', 'modular', 'wireless-power', 'hydra'],
    'layers': {l[1]: {'name': l[1], 'color': l[3], 'material': l[2]} for l in layers},
    'geometry': {
        'primary': {
            'path': 'geometry/primary.json',
            'schema': 'yapcad-geometry-json-v0.1',
            'entities': [e['id'] for e in entities],
            'hash': f'sha256:{geom_hash}',
        },
    },
    'attachments': [{
        'id': 'dsl-source',
        'path': 'attachments/source.dsl',
        'hash': f'sha256:{hashlib.sha256(source.encode()).hexdigest()}',
        'format': 'dsl',
        'purpose': 'DSL source code',
    }],
}

import yaml
(target / 'manifest.yaml').write_text(yaml.dump(manifest, default_flow_style=False, sort_keys=False))

print(f"\nPackage created at: {target}")
print(f"Entities: {len(entities)}")
print(f"Layers: {[l[1] for l in layers]}")
