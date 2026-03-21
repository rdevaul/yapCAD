# Task: YapCADViewer Class - Phase 1

## Context
We're building a Three.js-based 3D viewer for yapCAD geometry packages. This task generates the core viewer class with material presets and basic rendering.

## Requirements

### Input
- Parsed geometry from yapCAD packages (vertices, normals, faces as typed arrays)
- Material name from entity metadata (e.g., "brass", "steel", "aluminum", "plastic_white")

### Output
TypeScript class `YapCADViewer` that:

1. **Initialization**
   - Creates Three.js scene, camera, renderer
   - Sets up OrbitControls
   - Configures 3-point lighting (key, fill, back)
   - Dark background (#1a1a2e)

2. **Material Presets**
   - brass: { color: 0xb5a642, metalness: 0.8, roughness: 0.3 }
   - steel: { color: 0x8090a0, metalness: 0.9, roughness: 0.4 }
   - aluminum: { color: 0xd0d5dd, metalness: 0.7, roughness: 0.5 }
   - copper: { color: 0xb87333, metalness: 0.85, roughness: 0.35 }
   - plastic_white: { color: 0xf0f0f0, metalness: 0.0, roughness: 0.6 }
   - plastic_black: { color: 0x202020, metalness: 0.0, roughness: 0.6 }
   - plastic_red: { color: 0xcc3333, metalness: 0.0, roughness: 0.6 }
   - plastic_blue: { color: 0x3366cc, metalness: 0.0, roughness: 0.6 }
   - default: { color: 0x6688aa, metalness: 0.3, roughness: 0.5 }

3. **Geometry Loading**
   - Method: `loadGeometry(entities: GeometryEntity[])`
   - Creates BufferGeometry from vertices/normals/indices
   - Assigns material based on entity metadata
   - Adds to scene

4. **Camera Fitting**
   - Method: `fitToGeometry()`
   - Computes bounding box of all geometry
   - Positions camera to see everything with padding

5. **Render Loop**
   - RequestAnimationFrame-based
   - 30fps target
   - Method: `start()`, `stop()`

6. **Cleanup**
   - Method: `dispose()`
   - Properly disposes geometries, materials, renderer

### Interface

```typescript
interface GeometryEntity {
  id: string;
  name?: string;
  vertices: Float32Array;
  normals?: Float32Array;
  indices?: Uint32Array;
  material?: string;
  layer?: string;
}

interface ViewerOptions {
  container: HTMLElement;
  antialias?: boolean;
  pixelRatio?: number;
}
```

### Constraints
- Use MeshStandardMaterial for all surfaces
- DoubleSide rendering (yapCAD geometry may have inconsistent winding)
- Enable renderer.localClippingEnabled (for future clipping planes)
- Export as ES module

## Output Format
Return ONLY the TypeScript code for the YapCADViewer class. No markdown, no explanations.
