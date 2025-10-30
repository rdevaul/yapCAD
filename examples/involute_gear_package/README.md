# Involute Gear Package Example

This example demonstrates how to build a `.ycpkg` package that supplies parametric involute spur gears and how that package can be consumed inside another design.
The gear profile math is derived from the open-source [figgear](https://github.com/chromia/figgear) generator (MIT Licensed) and is vendored into `yapcad.contrib.figgear`, so no extra dependency is required.

## Files

```
examples/involute_gear_package/
├── README.md
├── involute_gear.py            # Python helpers that generate 2D profile + 3D solid
├── build_gear_package.py       # Creates a .ycpkg with canonical gear geometry
├── sample_drive.py             # Uses the gear package to build a simple gear pair assembly
├── src/
│   └── involute_gear.dsl       # Draft DSL module describing INVOLUTE_SPUR commands
└── scripts/
    └── involute_helpers.py     # Placeholder for future DSL Python fallback blocks
```

## Usage

`build_gear_package.py` copies the DSL module and helper scripts into the package and records their hashes inside `manifest.source`, so every canonical gear bundle remains self-contained.

1. Create canonical gear packages (adjust tooth counts as needed):

```bash
python examples/involute_gear_package/build_gear_package.py \
    --output /tmp/driver_gear.ycpkg \
    --teeth 18 \
    --module-mm 2.0 \
    --face-width-mm 8.0

python examples/involute_gear_package/build_gear_package.py \
    --output /tmp/driven_gear.ycpkg \
    --teeth 36 \
    --module-mm 2.0 \
    --face-width-mm 8.0
```

2. Build a simple two-gear drive package that reuses those canonical bundles:

```bash
python examples/involute_gear_package/sample_drive.py \
    --driver-package /tmp/driver_gear.ycpkg \
    --driven-package /tmp/driven_gear.ycpkg \
    --output /tmp/gear_drive.ycpkg
```

Passing `--driver-package` / `--driven-package` makes `sample_drive.py` load and reuse the canonical solids (and their invocation metadata) rather than regenerating geometry. If these flags are omitted the script falls back to generating gears in-place using the helper functions.

Each package embeds both the 2D profile and extruded solid along with invocation metadata (package/command/parameters) so downstream tools know which DSL command produced the geometry. Layers are set to `gear-profile`, `gear`, etc., allowing the viewer to toggle visibility by layer.
