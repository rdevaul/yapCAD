/**
 * Type definitions for yapCAD packages in the web-viewer.
 */

export interface PackageManifest {
  name: string;
  version: string;
  schema_version: string;
  uuid?: string;
  description?: string;
  author?: string;
  created_at?: string;
  updated_at?: string;
  tags: string[];
  geometry_files: string[];
  material_files: string[];
  thumbnail?: string;
}

export type PackageSource = 'local' | 'server' | 'git' | 'gel';

export interface PackageEntry {
  /** Display ID (filename or synthetic) */
  id: string;
  /** Stable UUID from manifest (preferred) or generated */
  uuid: string;
  /** Package name from manifest */
  name: string;
  /** Version string */
  version: string;
  /** Source type */
  source: PackageSource;
  /** Base64 thumbnail or data URL */
  thumbnail?: string;
  /** Full manifest */
  manifest: PackageManifest;
  /** File handle for watching (local files only) */
  fileHandle?: FileSystemFileHandle;
  /** Last modified timestamp */
  lastModified?: number;
  /** Parsed geometry ready for rendering */
  geometry?: ParsedGeometry;
  /** Loading state */
  loading?: boolean;
  /** Error message if loading failed */
  error?: string;
}

export interface ParsedGeometry {
  /** Three.js-ready buffer geometry data */
  entities: GeometryEntity[];
  /** Bounding box for camera framing */
  boundingBox?: {
    min: [number, number, number];
    max: [number, number, number];
  };
}

export interface GeometryEntity {
  id: string;
  type: string;
  name?: string;
  layer?: string;
  vertices: Float32Array;
  normals?: Float32Array;
  indices?: Uint32Array;
  material?: string;
}

export interface PackageSelectorState {
  packages: PackageEntry[];
  selectedId: string | null;
  loading: boolean;
  error?: string;
}

export type PackageSelectorAction =
  | { type: 'ADD_PACKAGE'; payload: PackageEntry }
  | { type: 'REMOVE_PACKAGE'; payload: string }
  | { type: 'SELECT_PACKAGE'; payload: string | null }
  | { type: 'UPDATE_PACKAGE'; payload: Partial<PackageEntry> & { uuid: string } }
  | { type: 'SET_LOADING'; payload: boolean }
  | { type: 'SET_ERROR'; payload: string | undefined };
