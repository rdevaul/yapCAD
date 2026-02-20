/**
 * Utility functions for parsing yapCAD .ycpkg files in the browser.
 */

import JSZip from 'jszip';
import { v4 as uuidv4 } from 'uuid';
import yaml from 'js-yaml';
import type { PackageManifest, PackageEntry, ParsedGeometry } from '../types/package';

/**
 * Extract geometry file paths from various manifest formats.
 * Handles both simple array format and nested yapCAD geometry spec.
 */
function extractGeometryFiles(manifestData: Record<string, unknown>): string[] {
  // Format 1: geometry_files array (web-viewer format)
  if (Array.isArray(manifestData.geometry_files) && manifestData.geometry_files.length > 0) {
    return manifestData.geometry_files as string[];
  }
  
  // Format 2: Nested geometry object (yapCAD native format)
  // geometry:
  //   primary:
  //     path: geometry/primary.json
  const geometry = manifestData.geometry as Record<string, unknown> | undefined;
  if (geometry && typeof geometry === 'object' && !Array.isArray(geometry)) {
    const paths: string[] = [];
    for (const [key, value] of Object.entries(geometry)) {
      if (value && typeof value === 'object' && 'path' in value) {
        const pathValue = (value as Record<string, unknown>).path;
        if (typeof pathValue === 'string') {
          paths.push(pathValue);
          console.log(`[packageParser] Found geometry file: ${key} -> ${pathValue}`);
        }
      }
    }
    if (paths.length > 0) {
      return paths;
    }
  }
  
  // Default fallback
  console.log('[packageParser] Using default geometry path: geometry/primary.json');
  return ['geometry/primary.json'];
}

/**
 * Parse a .ycpkg file (ZIP archive) and extract its contents.
 */
export async function parsePackageFile(file: File): Promise<PackageEntry> {
  let zip: JSZip;
  
  try {
    zip = await JSZip.loadAsync(file);
  } catch (err) {
    throw new Error(`Failed to open package as ZIP: ${err instanceof Error ? err.message : String(err)}`);
  }
  
  // Debug: list all files in the archive
  const fileList = Object.keys(zip.files);
  console.log('Package contents:', fileList);
  
  // Read manifest - try YAML first (yapCAD native), then JSON (web-viewer format)
  let manifestFile = zip.file('manifest.yaml');
  let isYaml = true;
  
  if (!manifestFile) {
    manifestFile = zip.file('manifest.json');
    isYaml = false;
  }
  
  if (!manifestFile) {
    throw new Error(`Invalid package: missing manifest.yaml or manifest.json. Found files: ${fileList.join(', ')}`);
  }
  
  let manifestContent: string;
  try {
    manifestContent = await manifestFile.async('string');
  } catch (err) {
    throw new Error(`Failed to read manifest: ${err instanceof Error ? err.message : String(err)}`);
  }
  let manifestData: Record<string, unknown>;
  
  if (isYaml) {
    manifestData = yaml.load(manifestContent) as Record<string, unknown>;
    console.log('[packageParser] Parsed YAML manifest:', manifestData);
  } else {
    manifestData = JSON.parse(manifestContent);
    console.log('[packageParser] Parsed JSON manifest:', manifestData);
  }
  
  // Normalize manifest to our expected format
  const manifest: PackageManifest = {
    name: String(manifestData.name || file.name.replace(/\.ycpkg$/i, '')),
    version: String(manifestData.version || '0.0.0'),
    schema_version: String(manifestData.schema_version || manifestData.schema || '1.0'),
    uuid: manifestData.uuid as string | undefined,
    description: manifestData.description as string | undefined,
    author: manifestData.author as string | undefined,
    created_at: manifestData.created_at as string | undefined,
    updated_at: manifestData.updated_at as string | undefined,
    tags: Array.isArray(manifestData.tags) ? manifestData.tags as string[] : [],
    // Determine geometry files - handle various manifest formats
    geometry_files: extractGeometryFiles(manifestData),
    material_files: Array.isArray(manifestData.material_files) ? manifestData.material_files as string[] : [],
    thumbnail: manifestData.thumbnail as string | undefined,
  };
  
  // Generate or use existing UUID
  const uuid = manifest.uuid || generateFallbackUuid(file, manifest);
  
  // Extract thumbnail if present
  let thumbnail: string | undefined;
  if (manifest.thumbnail) {
    const thumbFile = zip.file(manifest.thumbnail);
    if (thumbFile) {
      const thumbData = await thumbFile.async('base64');
      const ext = manifest.thumbnail.split('.').pop()?.toLowerCase() || 'png';
      const mimeType = ext === 'jpg' || ext === 'jpeg' ? 'image/jpeg' : 'image/png';
      thumbnail = `data:${mimeType};base64,${thumbData}`;
    }
  }
  
  // Parse geometry files
  const geometry = await parseGeometryFiles(zip, manifest.geometry_files);
  
  return {
    id: file.name.replace(/\.ycpkg$/i, ''),
    uuid,
    name: manifest.name,
    version: manifest.version,
    source: 'local',
    thumbnail,
    manifest,
    lastModified: file.lastModified,
    geometry,
  };
}

/**
 * Generate a deterministic fallback UUID for packages without one.
 */
function generateFallbackUuid(file: File, manifest: PackageManifest): string {
  // For packages without UUID, generate one from content characteristics
  // This is stable for the same file but not truly unique
  const seed = `${file.name}:${manifest.name}:${manifest.version}:${file.size}`;
  
  // Simple hash-based approach (not cryptographic, just for ID generation)
  let hash = 0;
  for (let i = 0; i < seed.length; i++) {
    const char = seed.charCodeAt(i);
    hash = ((hash << 5) - hash) + char;
    hash = hash & hash; // Convert to 32-bit integer
  }
  
  // Return a synthetic UUID-like string
  return `local-${Math.abs(hash).toString(16).padStart(8, '0')}-${uuidv4().slice(-12)}`;
}

/**
 * Parse geometry JSON files from the package.
 */
async function parseGeometryFiles(
  zip: JSZip, 
  geometryPaths: string[]
): Promise<ParsedGeometry> {
  const entities: ParsedGeometry['entities'] = [];
  let minBounds = [Infinity, Infinity, Infinity];
  let maxBounds = [-Infinity, -Infinity, -Infinity];
  
  console.log('[parseGeometryFiles] Looking for geometry in paths:', geometryPaths);
  
  for (const path of geometryPaths) {
    console.log('[parseGeometryFiles] Trying path:', path);
    const file = zip.file(path) || zip.file(`geometry/${path}`);
    if (!file) {
      console.warn(`[parseGeometryFiles] Geometry file not found: ${path}`);
      continue;
    }
    console.log('[parseGeometryFiles] Found geometry file:', path);
    
    const json = await file.async('string');
    const data = JSON.parse(json);
    
    // Parse yapCAD geometry JSON format
    if (data.entities && Array.isArray(data.entities)) {
      for (const entity of data.entities) {
        const parsed = parseEntity(entity);
        if (parsed) {
          entities.push(parsed);
          
          // Update bounding box
          if (entity.boundingBox) {
            const [minX, minY, minZ, maxX, maxY, maxZ] = entity.boundingBox;
            minBounds = [
              Math.min(minBounds[0], minX),
              Math.min(minBounds[1], minY),
              Math.min(minBounds[2], minZ),
            ];
            maxBounds = [
              Math.max(maxBounds[0], maxX),
              Math.max(maxBounds[1], maxY),
              Math.max(maxBounds[2], maxZ),
            ];
          }
        }
      }
    }
  }
  
  console.log('[parseGeometryFiles] Parsed', entities.length, 'entities total');
  for (const e of entities) {
    console.log('  - Entity:', e.id, 'vertices:', e.vertices.length / 3, 'indices:', e.indices?.length);
  }
  
  return {
    entities,
    boundingBox: entities.length > 0 ? {
      min: minBounds as [number, number, number],
      max: maxBounds as [number, number, number],
    } : undefined,
  };
}

/**
 * Parse a single geometry entity into render-ready format.
 */
function parseEntity(entity: any): ParsedGeometry['entities'][0] | null {
  console.log('[parseEntity] Processing entity:', entity.id, 'type:', entity.type);
  
  if (!entity.id || !entity.type) {
    console.log('[parseEntity] Skipping - missing id or type');
    return null;
  }
  
  // Only process surface entities (they have mesh data)
  // Skip solid entities (BREP only, no mesh)
  if (entity.type !== 'surface') {
    console.log('[parseEntity] Skipping non-surface entity:', entity.type);
    return null;
  }
  
  // Handle mesh/surface entities with vertices
  if (entity.vertices && Array.isArray(entity.vertices)) {
    console.log('[parseEntity] Found vertices:', entity.vertices.length, 'faces:', entity.faces?.length);
    const flatVertices: number[] = [];
    const flatNormals: number[] = [];
    
    for (let i = 0; i < entity.vertices.length; i++) {
      const v = entity.vertices[i];
      flatVertices.push(v[0], v[1], v[2]);
      
      if (entity.normals && entity.normals[i]) {
        const n = entity.normals[i];
        flatNormals.push(n[0], n[1], n[2]);
      }
    }
    
    // Build face indices
    let indices: number[] | undefined;
    if (entity.faces && Array.isArray(entity.faces)) {
      indices = [];
      for (const face of entity.faces) {
        if (face.length === 3) {
          indices.push(face[0], face[1], face[2]);
        } else if (face.length === 4) {
          // Triangulate quads
          indices.push(face[0], face[1], face[2]);
          indices.push(face[0], face[2], face[3]);
        }
      }
    }
    
    return {
      id: entity.id,
      type: entity.type,
      name: entity.name,
      layer: entity.metadata?.layer,
      vertices: new Float32Array(flatVertices),
      normals: flatNormals.length > 0 ? new Float32Array(flatNormals) : undefined,
      indices: indices ? new Uint32Array(indices) : undefined,
      material: entity.metadata?.material,
    };
  }
  
  return null;
}

/**
 * Generate a thumbnail by rendering geometry to an offscreen canvas.
 * Returns base64 PNG data URL.
 */
/**
 * Load a package from a URL (HTTP/HTTPS).
 * Can be a remote ZIP file or a local path served via backend proxy.
 */
export async function loadPackageFromUrl(url: string): Promise<PackageEntry> {
  const response = await fetch(url);
  
  if (!response.ok) {
    throw new Error(`Failed to fetch package: ${response.status} ${response.statusText}`);
  }
  
  const blob = await response.blob();
  const filename = url.split('/').pop()?.split('?')[0] || 'package.ycpkg';
  const file = new File([blob], filename, { type: 'application/zip' });
  
  return parsePackageFile(file);
}

/**
 * Load a package from a local path via the backend proxy.
 * The backend will zip directories on-the-fly if needed.
 */
export async function loadPackageFromPath(localPath: string): Promise<PackageEntry> {
  const url = `/api/local-package?path=${encodeURIComponent(localPath)}`;
  return loadPackageFromUrl(url);
}

export async function generateThumbnail(
  geometry: ParsedGeometry,
  size: number = 256
): Promise<string> {
  const { renderThumbnail, getCachedThumbnail, cacheThumbnail } = await import('./thumbnailRenderer');

  // Check cache first (use a hash of entity count + first entity id as key)
  const cacheKey = geometry.entities.map(e => e.id).join('-');
  const cached = await getCachedThumbnail(cacheKey);
  if (cached) return cached;

  // Convert ParsedGeometry entities to the format renderThumbnail expects
  const entities = geometry.entities.map(e => ({
    id: e.id,
    name: e.name,
    vertices: e.vertices,
    normals: e.normals,
    indices: e.indices,
    material: e.material,
    layer: e.layer,
  }));

  try {
    const dataUrl = await renderThumbnail(entities, size);
    await cacheThumbnail(cacheKey, dataUrl);
    return dataUrl;
  } catch (err) {
    console.warn('Thumbnail rendering failed, using placeholder:', err);
    // Fallback placeholder
    const canvas = document.createElement('canvas');
    canvas.width = size;
    canvas.height = size;
    const ctx = canvas.getContext('2d');
    if (ctx) {
      const gradient = ctx.createLinearGradient(0, 0, size, size);
      gradient.addColorStop(0, '#3b82f6');
      gradient.addColorStop(1, '#8b5cf6');
      ctx.fillStyle = gradient;
      ctx.fillRect(0, 0, size, size);
      ctx.fillStyle = 'white';
      ctx.font = '12px sans-serif';
      ctx.textAlign = 'center';
      ctx.fillText('yapCAD', size / 2, size / 2);
    }
    return canvas.toDataURL('image/png');
  }
}
