/**
 * yapCAD Geometry JSON Loader for Three.js
 * 
 * Parses yapcad-geometry-json-v0.1 format and creates Three.js BufferGeometry.
 * Preserves full JSON structure for round-trip operations.
 */

import * as THREE from 'three';

// ============================================================================
// Types matching yapcad-geometry-json-v0.1 schema
// ============================================================================

export interface YapCADDocument {
  schema: string;
  generator?: {
    name: string;
    version: string;
    build?: string;
  };
  units?: string;
  entities: YapCADEntity[];
  relationships?: YapCADRelationship[];
  attachments?: YapCADAttachment[];
}

export interface YapCADEntity {
  id: string;
  type: 'solid' | 'surface' | 'mesh' | 'group' | 'datum' | 'sketch';
  name?: string | null;
  metadata: YapCADMetadata;
  boundingBox?: number[] | null;
  properties?: Record<string, unknown>;
  // Surface-specific
  vertices?: number[][];
  normals?: number[][];
  faces?: number[][];
  triangulation?: {
    winding: 'ccw' | 'cw';
    topology: 'triangle';
  };
  // Solid-specific
  shell?: string[];
  voids?: string[][];
  // Group-specific
  children?: string[];
  transform?: number[][];
}

export interface YapCADMetadata {
  schema?: string;
  entityId?: string;
  id?: string;
  layer?: string;
  tags?: string[];
  material?: string;
  brep?: {
    encoding: string;
    data: string;
  };
  [key: string]: unknown;
}

export interface YapCADRelationship {
  type: string;
  entities: string[];
  constraint?: Record<string, unknown>;
}

export interface YapCADAttachment {
  id: string;
  kind: string;
  path: string;
  hash?: string;
  createdBy?: string;
  metadata?: Record<string, unknown>;
}

// ============================================================================
// Material definitions (matches package/viewer.py defaults)
// ============================================================================

export interface MaterialDefinition {
  color: [number, number, number];
  metallic?: number;
  roughness?: number;
}

const DEFAULT_COLOR: [number, number, number] = [0.6, 0.85, 1.0]; // yapCAD blue

// ============================================================================
// Loader result types
// ============================================================================

export interface LoadedSurface {
  id: string;
  geometry: THREE.BufferGeometry;
  layer: string;
  materialId?: string;
}

export interface LoadedSolid {
  id: string;
  surfaces: LoadedSurface[];
  layer: string;
  materialId?: string;
  boundingBox?: THREE.Box3;
}

export interface LoaderResult {
  /** Original document preserved for round-trip */
  document: YapCADDocument;
  /** Loaded solids with Three.js geometries */
  solids: LoadedSolid[];
  /** Standalone surfaces not part of a solid */
  surfaces: LoadedSurface[];
  /** Computed bounding box of all geometry */
  boundingBox: THREE.Box3;
  /** Unique layers found in the document */
  layers: Set<string>;
  /** Unique material IDs referenced */
  materialIds: Set<string>;
}

// ============================================================================
// Core loader function
// ============================================================================

/**
 * Load a yapCAD geometry JSON document and create Three.js geometries.
 * 
 * @param json - The parsed JSON document or JSON string
 * @returns LoaderResult with Three.js geometries and metadata
 */
export function loadYapCADGeometry(json: YapCADDocument | string): LoaderResult {
  const doc: YapCADDocument = typeof json === 'string' ? JSON.parse(json) : json;
  
  // Validate schema
  if (!doc.schema?.startsWith('yapcad-geometry-json-')) {
    console.warn(`Unknown schema: ${doc.schema}, attempting to parse anyway`);
  }
  
  // Build lookup maps
  const surfaceMap = new Map<string, YapCADEntity>();
  const solidEntities: YapCADEntity[] = [];
  const standaloneSurfaces: YapCADEntity[] = [];
  const surfacesInSolids = new Set<string>();
  
  for (const entity of doc.entities) {
    if (entity.type === 'surface') {
      surfaceMap.set(entity.id, entity);
    } else if (entity.type === 'solid') {
      solidEntities.push(entity);
      // Mark surfaces that belong to this solid
      for (const surfId of entity.shell || []) {
        surfacesInSolids.add(surfId);
      }
      for (const voidShell of entity.voids || []) {
        for (const surfId of voidShell) {
          surfacesInSolids.add(surfId);
        }
      }
    }
  }
  
  // Find standalone surfaces (not part of any solid)
  for (const entity of doc.entities) {
    if (entity.type === 'surface' && !surfacesInSolids.has(entity.id)) {
      standaloneSurfaces.push(entity);
    }
  }
  
  // Process solids
  const loadedSolids: LoadedSolid[] = [];
  const layers = new Set<string>();
  const materialIds = new Set<string>();
  const overallBox = new THREE.Box3();
  
  for (const solidEntity of solidEntities) {
    const solidLayer = solidEntity.metadata?.layer || 'default';
    const solidMaterial = solidEntity.metadata?.material as string | undefined;
    layers.add(solidLayer);
    if (solidMaterial) materialIds.add(solidMaterial);
    
    const loadedSurfaces: LoadedSurface[] = [];
    
    for (const surfId of solidEntity.shell || []) {
      const surfEntity = surfaceMap.get(surfId);
      if (!surfEntity) {
        console.warn(`Surface ${surfId} not found for solid ${solidEntity.id}`);
        continue;
      }
      
      const geometry = createBufferGeometry(surfEntity);
      if (geometry) {
        geometry.computeBoundingBox();
        if (geometry.boundingBox) {
          overallBox.union(geometry.boundingBox);
        }
        
        const surfLayer = surfEntity.metadata?.layer || solidLayer;
        const surfMaterial = surfEntity.metadata?.material as string | undefined;
        layers.add(surfLayer);
        if (surfMaterial) materialIds.add(surfMaterial);
        
        loadedSurfaces.push({
          id: surfId,
          geometry,
          layer: surfLayer,
          materialId: surfMaterial || solidMaterial,
        });
      }
    }
    
    let solidBox: THREE.Box3 | undefined;
    if (solidEntity.boundingBox && solidEntity.boundingBox.length === 6) {
      const [xmin, ymin, zmin, xmax, ymax, zmax] = solidEntity.boundingBox;
      solidBox = new THREE.Box3(
        new THREE.Vector3(xmin, ymin, zmin),
        new THREE.Vector3(xmax, ymax, zmax)
      );
    }
    
    loadedSolids.push({
      id: solidEntity.id,
      surfaces: loadedSurfaces,
      layer: solidLayer,
      materialId: solidMaterial,
      boundingBox: solidBox,
    });
  }
  
  // Process standalone surfaces
  const loadedStandaloneSurfaces: LoadedSurface[] = [];
  
  for (const surfEntity of standaloneSurfaces) {
    const geometry = createBufferGeometry(surfEntity);
    if (geometry) {
      geometry.computeBoundingBox();
      if (geometry.boundingBox) {
        overallBox.union(geometry.boundingBox);
      }
      
      const layer = surfEntity.metadata?.layer || 'default';
      const materialId = surfEntity.metadata?.material as string | undefined;
      layers.add(layer);
      if (materialId) materialIds.add(materialId);
      
      loadedStandaloneSurfaces.push({
        id: surfEntity.id,
        geometry,
        layer,
        materialId,
      });
    }
  }
  
  return {
    document: doc,
    solids: loadedSolids,
    surfaces: loadedStandaloneSurfaces,
    boundingBox: overallBox,
    layers,
    materialIds,
  };
}

// ============================================================================
// BufferGeometry creation
// ============================================================================

/**
 * Create a Three.js BufferGeometry from a yapCAD surface entity.
 */
function createBufferGeometry(entity: YapCADEntity): THREE.BufferGeometry | null {
  if (!entity.vertices || !entity.faces) {
    console.warn(`Surface ${entity.id} missing vertices or faces`);
    return null;
  }
  
  const vertices = entity.vertices;
  const normals = entity.normals || [];
  const faces = entity.faces;
  
  // Calculate flat arrays for BufferGeometry
  // Each face is a triangle with 3 vertices
  const positionArray: number[] = [];
  const normalArray: number[] = [];
  const indexArray: number[] = [];
  
  // For indexed geometry, we need to flatten vertices/normals
  // yapCAD uses homogeneous coords [x, y, z, w] - we only need [x, y, z]
  for (let i = 0; i < vertices.length; i++) {
    const v = vertices[i];
    positionArray.push(v[0], v[1], v[2]);
    
    if (i < normals.length) {
      const n = normals[i];
      normalArray.push(n[0], n[1], n[2]);
    } else {
      // Default normal if missing
      normalArray.push(0, 0, 1);
    }
  }
  
  // Faces are triangle indices
  for (const face of faces) {
    if (face.length === 3) {
      indexArray.push(face[0], face[1], face[2]);
    } else {
      console.warn(`Non-triangle face in surface ${entity.id}:`, face);
    }
  }
  
  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute('position', new THREE.Float32BufferAttribute(positionArray, 3));
  geometry.setAttribute('normal', new THREE.Float32BufferAttribute(normalArray, 3));
  geometry.setIndex(indexArray);
  
  // Store entity ID for picking/selection
  geometry.userData = { entityId: entity.id };
  
  return geometry;
}

// ============================================================================
// Material helpers
// ============================================================================

/**
 * Create a Three.js MeshStandardMaterial from a material definition.
 */
export function createMaterial(
  definition?: MaterialDefinition,
  options?: { wireframe?: boolean; side?: THREE.Side }
): THREE.MeshStandardMaterial {
  const color = definition?.color || DEFAULT_COLOR;
  
  return new THREE.MeshStandardMaterial({
    color: new THREE.Color(color[0], color[1], color[2]),
    metalness: definition?.metallic ?? 0.1,
    roughness: definition?.roughness ?? 0.6,
    wireframe: options?.wireframe ?? false,
    side: options?.side ?? THREE.DoubleSide,
  });
}

/**
 * Create a default material palette for layers.
 */
export function createLayerMaterials(layers: Set<string>): Map<string, THREE.MeshStandardMaterial> {
  const materials = new Map<string, THREE.MeshStandardMaterial>();
  
  // Predefined colors for common layers
  const layerColors: Record<string, [number, number, number]> = {
    default: [0.6, 0.85, 1.0],      // yapCAD blue
    body: [0.8, 0.8, 0.8],          // gray
    highlight: [1.0, 0.5, 0.2],     // orange
    transparent: [0.9, 0.9, 0.95],  // near-white
  };
  
  // Generate colors for unknown layers using HSL rotation
  let hue = 0;
  for (const layer of layers) {
    if (layerColors[layer]) {
      materials.set(layer, createMaterial({ color: layerColors[layer] }));
    } else {
      // Generate a color by rotating through hue
      const color = new THREE.Color().setHSL(hue, 0.7, 0.6);
      materials.set(layer, createMaterial({ 
        color: [color.r, color.g, color.b] 
      }));
      hue = (hue + 0.618033988749895) % 1; // Golden ratio for nice distribution
    }
  }
  
  return materials;
}

// ============================================================================
// Scene builder
// ============================================================================

export interface SceneOptions {
  materials?: Map<string, MaterialDefinition>;
  wireframe?: boolean;
  visibleLayers?: Set<string>;
}

/**
 * Build a Three.js Group from loader results.
 */
export function buildScene(result: LoaderResult, options: SceneOptions = {}): THREE.Group {
  const root = new THREE.Group();
  root.name = 'yapCAD-scene';
  
  // Create materials
  const layerMaterials = createLayerMaterials(result.layers);
  
  // Override with provided materials
  if (options.materials) {
    for (const [id, def] of options.materials) {
      layerMaterials.set(id, createMaterial(def, { wireframe: options.wireframe }));
    }
  }
  
  // Add solids
  for (const solid of result.solids) {
    const solidGroup = new THREE.Group();
    solidGroup.name = `solid-${solid.id}`;
    solidGroup.userData = { entityId: solid.id, type: 'solid' };
    
    // Check layer visibility
    if (options.visibleLayers && !options.visibleLayers.has(solid.layer)) {
      solidGroup.visible = false;
    }
    
    for (const surface of solid.surfaces) {
      const materialKey = surface.materialId || surface.layer;
      const material = layerMaterials.get(materialKey) || layerMaterials.get('default')!;
      
      const mesh = new THREE.Mesh(surface.geometry, material);
      mesh.name = `surface-${surface.id}`;
      mesh.userData = { entityId: surface.id, type: 'surface', parentSolid: solid.id };
      
      solidGroup.add(mesh);
    }
    
    root.add(solidGroup);
  }
  
  // Add standalone surfaces
  for (const surface of result.surfaces) {
    const materialKey = surface.materialId || surface.layer;
    const material = layerMaterials.get(materialKey) || layerMaterials.get('default')!;
    
    const mesh = new THREE.Mesh(surface.geometry, material);
    mesh.name = `surface-${surface.id}`;
    mesh.userData = { entityId: surface.id, type: 'surface' };
    
    if (options.visibleLayers && !options.visibleLayers.has(surface.layer)) {
      mesh.visible = false;
    }
    
    root.add(mesh);
  }
  
  return root;
}

// ============================================================================
// Utility functions
// ============================================================================

/**
 * Compute camera position to frame the bounding box.
 */
export function computeCameraPosition(
  boundingBox: THREE.Box3,
  fov: number = 50
): { position: THREE.Vector3; target: THREE.Vector3 } {
  const center = new THREE.Vector3();
  const size = new THREE.Vector3();
  boundingBox.getCenter(center);
  boundingBox.getSize(size);
  
  const maxDim = Math.max(size.x, size.y, size.z);
  const fovRad = (fov * Math.PI) / 180;
  const distance = maxDim / (2 * Math.tan(fovRad / 2)) * 1.5;
  
  // Position camera at 45° angle
  const position = new THREE.Vector3(
    center.x + distance * 0.7,
    center.y + distance * 0.5,
    center.z + distance * 0.7
  );
  
  return { position, target: center };
}

export default {
  loadYapCADGeometry,
  createMaterial,
  createLayerMaterials,
  buildScene,
  computeCameraPosition,
};
