# Material Schema Specification

**Version**: 0.1 (Draft)
**Status**: Proposal for Review
**Date**: November 2025

## Overview

This specification defines how materials are represented in yapCAD packages. The design prioritizes engineering data over visual styling, enabling:

- FEA/simulation integration
- Manufacturing process planning
- Bill of materials generation
- Standards compliance tracking

Visual properties for rendering are secondary and derived from or associated with material definitions.

## Design Principles

1. **Standards-First**: Reference industry standards (AA, ASTM, SAE, ISO) when available
2. **Engineering Data**: Physical properties are primary; visual properties are auxiliary
3. **Traceability**: Materials link to sources (standards bodies, vendors, test data)
4. **Extensibility**: Custom properties can be added without breaking compatibility
5. **Separation of Concerns**: Layers organize geometry; materials define substance

## Schema Structure

### Package Manifest Integration

Materials are defined in the package manifest under a `materials` key:

```yaml
# manifest.yaml
name: "bracket-assembly"
version: "1.0.0"
units: "mm"

materials:
  aluminum_6061:
    # ... material definition
  pla_black:
    # ... material definition

geometry:
  primary:
    path: "geometry/assembly.json"
```

Alternatively, materials may be defined in a separate file:

```yaml
# manifest.yaml
materials:
  path: "materials.yaml"  # External file reference
```

### Material Definition Schema

```yaml
material_id:                    # Unique identifier within package
  # === Source Information (required) ===
  source:
    type: "standard" | "vendor" | "custom" | "tested"

    # For type: "standard"
    standard:
      body: "AA" | "ASTM" | "SAE" | "ISO" | "DIN" | "JIS" | ...
      designation: "6061-T6"    # Standard's material designation
      revision: "2023"          # Optional: standard revision/year

    # For type: "vendor"
    vendor:
      name: "Prusament"
      product_id: "PLA Galaxy Black"
      lot: "2024-Q3"            # Optional: lot/batch
      url: "https://..."        # Optional: product page
      datasheet: "datasheets/prusament_pla.pdf"  # Optional: local copy

    # For type: "tested"
    tested:
      lab: "Internal QA"
      date: "2024-11-15"
      report: "test_reports/al_batch_42.pdf"

    # For type: "custom"
    custom:
      author: "J. Smith"
      notes: "Estimated properties for prototype"

  # === Material Classification (recommended) ===
  classification:
    category: "metal" | "polymer" | "ceramic" | "composite" | "other"
    family: "aluminum" | "steel" | "titanium" | "thermoplastic" | ...

  # === Physical Properties (required for engineering use) ===
  properties:
    # Density
    density:
      value: 2.70
      unit: "g/cm3"

    # Mechanical - Elastic
    elastic_modulus:
      value: 68.9
      unit: "GPa"
    poisson_ratio:
      value: 0.33

    # Mechanical - Strength
    yield_strength:
      value: 276
      unit: "MPa"
      condition: "0.2% offset"  # Optional: test condition
    ultimate_tensile_strength:
      value: 310
      unit: "MPa"
    elongation:
      value: 12
      unit: "%"

    # Thermal
    thermal_conductivity:
      value: 167
      unit: "W/(m*K)"
    specific_heat:
      value: 896
      unit: "J/(kg*K)"
    melting_point:
      value: 582
      unit: "C"
      note: "Solidus"
    thermal_expansion:
      value: 23.6e-6
      unit: "1/K"

    # Electrical (optional)
    electrical_resistivity:
      value: 3.99e-8
      unit: "ohm*m"

    # Custom properties
    custom:
      hardness_brinell:
        value: 95
        unit: "HB"

  # === Manufacturing Constraints (optional) ===
  manufacturing:
    processes:
      - type: "CNC_milling"
        notes: "Good machinability"
      - type: "welding"
        method: "TIG"
        filler: "4043"

    # For additive manufacturing
    additive:
      process: "FDM" | "SLA" | "SLS" | "DMLS" | ...
      nozzle_temp: [200, 220]   # Range in C
      bed_temp: [50, 60]
      chamber_temp: null
      layer_height: [0.1, 0.3]  # Range in mm

  # === Visual Properties (for rendering) ===
  visual:
    # Base color (RGB, 0-1 range)
    color: [0.75, 0.78, 0.80]

    # PBR material properties
    metallic: 0.9              # 0 = dielectric, 1 = metal
    roughness: 0.3             # 0 = mirror, 1 = diffuse

    # Optional: more detailed appearance
    appearance:
      finish: "brushed" | "polished" | "matte" | "textured"
      anodized: false
      anodize_color: null

    # Optional: texture references
    textures:
      albedo: "textures/brushed_aluminum.png"
      normal: "textures/brushed_aluminum_normal.png"
      roughness: "textures/brushed_aluminum_rough.png"
```

### Compact Form for Simple Cases

For quick prototyping or when full data isn't needed:

```yaml
materials:
  generic_aluminum:
    source:
      type: "standard"
      standard: { body: "AA", designation: "6061-T6" }
    properties:
      density: { value: 2.70, unit: "g/cm3" }
    visual:
      color: [0.75, 0.78, 0.80]
      metallic: 0.9
      roughness: 0.3
```

### Minimal Form (Visual Only)

When physical properties aren't needed:

```yaml
materials:
  display_blue:
    source:
      type: "custom"
      custom: { notes: "Visual placeholder" }
    visual:
      color: [0.2, 0.4, 0.8]
      metallic: 0.0
      roughness: 0.5
```

## Geometry Association

Materials are associated with solids via metadata:

```json
{
  "entities": [
    {
      "type": "solid",
      "id": "bracket_001",
      "metadata": {
        "material": "aluminum_6061",
        "layer": "structural",
        "name": "Main Bracket"
      },
      "shell": ["surface_001", "surface_002"]
    }
  ]
}
```

### Association Rules

1. **Solid-level**: Each solid may have one material assignment
2. **Inheritance**: If no material specified, inherit from parent assembly
3. **Override**: Child components can override parent material
4. **Null material**: Explicit `"material": null` means unspecified (rendered with default)

## Built-in Material Library

yapCAD ships with a standard material library for common engineering materials:

### Metals
| ID | Standard | Description |
|----|----------|-------------|
| `AA-6061-T6` | AA 6061-T6 | General purpose aluminum |
| `AA-7075-T6` | AA 7075-T6 | High strength aluminum |
| `AISI-304` | AISI 304 | Austenitic stainless steel |
| `AISI-316` | AISI 316 | Marine grade stainless |
| `Ti-6Al-4V` | AMS 4911 | Titanium alloy |

### Polymers - 3D Printing
| ID | Type | Description |
|----|------|-------------|
| `PLA-generic` | PLA | Standard PLA filament |
| `ABS-generic` | ABS | Standard ABS filament |
| `PETG-generic` | PETG | Standard PETG filament |
| `Nylon-PA12` | PA12 | Nylon 12 |
| `Resin-standard` | SLA | Generic photopolymer resin |

### Polymers - Engineering
| ID | Standard | Description |
|----|----------|-------------|
| `Delrin-150` | ASTM D4181 | Acetal homopolymer |
| `PEEK-natural` | - | Polyetheretherketone |
| `UHMWPE` | ASTM D4020 | Ultra-high MW polyethylene |

### Usage

```yaml
materials:
  bracket_material:
    extends: "AA-6061-T6"  # Inherit from built-in
    visual:
      color: [0.8, 0.2, 0.2]  # Override color only
```

## Renderer Integration

### Material Resolution

The viewer resolves materials in this order:

1. Check solid's `metadata.material` reference
2. Look up in package `materials` dictionary
3. If `extends` present, merge with built-in library
4. Fall back to `default` material if not found

### Visual Property Mapping

| Schema Property | OpenGL/Shader Use |
|----------------|-------------------|
| `visual.color` | Diffuse color (with lighting) |
| `visual.metallic` | Fresnel F0 calculation |
| `visual.roughness` | Specular spread |

### Default Material

When no material is specified:

```yaml
_default:
  visual:
    color: [0.6, 0.85, 1.0]  # Current yapCAD blue
    metallic: 0.0
    roughness: 0.5
```

## Validation Rules

1. **Required fields**: `source.type` must be present
2. **Standard references**: If `source.type: "standard"`, must have `source.standard.body` and `source.standard.designation`
3. **Property units**: All physical properties should include units
4. **Color range**: RGB values must be in [0, 1]
5. **PBR range**: `metallic` and `roughness` must be in [0, 1]

## Future Extensions

### Planned (Priority)
- Integration with external material databases (MatWeb, Granta, manufacturer APIs)
- Anisotropic materials (fiber composites, wood, etc.)
- Composite layup definitions (ply orientation, stacking sequence)

### Planned
- Material cost data for BOM generation
- Environmental/sustainability properties
- Temperature-dependent properties
- Material testing workflow integration
- Certification/compliance tracking

## Examples

### Example 1: CNC Machined Part

```yaml
materials:
  bracket_al:
    source:
      type: "standard"
      standard:
        body: "AA"
        designation: "6061-T6"
    classification:
      category: "metal"
      family: "aluminum"
    properties:
      density: { value: 2.70, unit: "g/cm3" }
      yield_strength: { value: 276, unit: "MPa" }
      elastic_modulus: { value: 68.9, unit: "GPa" }
    manufacturing:
      processes:
        - type: "CNC_milling"
        - type: "anodizing"
          color: "black"
    visual:
      color: [0.1, 0.1, 0.1]
      metallic: 0.8
      roughness: 0.2
```

### Example 2: 3D Printed Prototype

```yaml
materials:
  prototype_pla:
    source:
      type: "vendor"
      vendor:
        name: "Hatchbox"
        product_id: "PLA 1.75mm Black"
    classification:
      category: "polymer"
      family: "thermoplastic"
    properties:
      density: { value: 1.24, unit: "g/cm3" }
      yield_strength: { value: 50, unit: "MPa", condition: "printed, 100% infill" }
    manufacturing:
      additive:
        process: "FDM"
        nozzle_temp: [200, 215]
        bed_temp: [55, 65]
        layer_height: [0.12, 0.28]
    visual:
      color: [0.05, 0.05, 0.05]
      metallic: 0.0
      roughness: 0.7
```

### Example 3: Multi-Material Assembly

```yaml
materials:
  housing_abs:
    source:
      type: "vendor"
      vendor: { name: "PolyMaker", product_id: "PolyLite ABS" }
    properties:
      density: { value: 1.04, unit: "g/cm3" }
    visual:
      color: [0.2, 0.2, 0.25]
      metallic: 0.0
      roughness: 0.4

  insert_brass:
    source:
      type: "standard"
      standard: { body: "ASTM", designation: "B16" }
    properties:
      density: { value: 8.47, unit: "g/cm3" }
    visual:
      color: [0.78, 0.57, 0.11]
      metallic: 1.0
      roughness: 0.35
```

```json
{
  "entities": [
    {
      "type": "solid",
      "id": "housing",
      "metadata": { "material": "housing_abs", "layer": "enclosure" }
    },
    {
      "type": "solid",
      "id": "thread_insert_1",
      "metadata": { "material": "insert_brass", "layer": "hardware" }
    }
  ]
}
```

## Changelog

- **0.1** (2025-11): Initial draft
