# UUID Schema Extension for yapCAD Packages

**Status:** Proposed  
**Author:** Jarvis (via Rich DeVaul)  
**Date:** 2026-02-09

---

## Summary

Add a `uuid` field to the yapCAD package manifest schema to provide stable, immutable package identification across versions and sessions.

## Motivation

The web-viewer needs to:
1. Track packages across browser sessions (IndexedDB caching)
2. Distinguish between "same package, new version" vs "different package"
3. Support future features like package lineage, comparison, and GEL integration

Currently, package identity is derived from filename (`demo-cube.ycpkg` → id: `demo-cube`). This breaks when:
- User renames a file
- Same package exists with different filenames
- Local files have no inherent stable identity

## Specification

### Manifest Schema Change

Add `uuid` field to `manifest.json`:

```json
{
  "name": "thrust-structure-optimized",
  "version": "1.2.0",
  "schema_version": "1.1",
  "uuid": "550e8400-e29b-41d4-a716-446655440000",
  "description": "Optimized thrust structure with lightening holes",
  "author": "Rich DeVaul",
  "created_at": "2026-02-09T10:30:00Z",
  "updated_at": "2026-02-09T14:20:00Z",
  "tags": ["rocket", "structure", "FEA-validated"],
  "geometry_files": ["geometry/primary.json"],
  "material_files": [],
  "thumbnail": "thumbnail.png"
}
```

### Field Definition

| Field | Type | Required | Mutable | Description |
|-------|------|----------|---------|-------------|
| `uuid` | string (UUID v4) | No* | No | Immutable package identifier |

*Optional for backwards compatibility. Tooling should auto-generate if missing.

### Rules

1. **Immutability:** Once set, UUID MUST NOT change for the lifetime of a package lineage
2. **Generation:** Use UUID v4 (random) for new packages
3. **Preservation:** When saving/updating a package, preserve existing UUID
4. **Versioning:** Package `version` field tracks changes; UUID identifies the lineage
5. **Uniqueness:** UUIDs should be globally unique (v4 collision probability is negligible)

### Backwards Compatibility

- Packages without `uuid` field remain valid
- Consumers should handle missing UUID gracefully:
  - Web-viewer: Generate synthetic UUID from content hash, store in local cache
  - yapCAD CLI: Prompt to add UUID on next save, or auto-generate

### Schema Version Bump

Increment `schema_version` from `"1.0"` to `"1.1"` to indicate UUID support.

## Implementation

### 1. Backend Models (models.py)

```python
class PackageManifest(BaseModel):
    """Metadata from manifest.json."""
    name: str
    version: str
    schema_version: str = Field(default="1.1")  # Bump version
    uuid: Optional[str] = None  # New field
    description: Optional[str] = None
    # ... rest unchanged
```

### 2. Package Creation (yapCAD core)

When creating a new package:
```python
import uuid

def create_manifest(name: str, version: str = "0.1.0", **kwargs) -> dict:
    return {
        "name": name,
        "version": version,
        "schema_version": "1.1",
        "uuid": str(uuid.uuid4()),
        **kwargs
    }
```

### 3. Package Update (yapCAD core)

When updating an existing package:
```python
def update_manifest(existing: dict, **changes) -> dict:
    # Preserve UUID, never overwrite
    updated = {**existing, **changes}
    updated["uuid"] = existing.get("uuid") or str(uuid.uuid4())
    return updated
```

### 4. Web-Viewer (frontend)

For packages without UUID:
```typescript
function getPackageId(manifest: PackageManifest, file: File): string {
  if (manifest.uuid) {
    return manifest.uuid;
  }
  // Fallback: generate deterministic ID from content
  return `local-${hashString(file.name + manifest.name + manifest.version)}`;
}
```

## Migration Path

1. **Phase 1 (web-viewer):** Handle missing UUID gracefully with fallback
2. **Phase 2 (yapCAD CLI):** Auto-generate UUID on package creation
3. **Phase 3 (existing packages):** Optional migration tool to add UUIDs

## Security Considerations

- UUID is not a security credential — don't use it for access control
- Existing hash-based signing system remains independent
- UUID can coexist with cryptographic signatures

## Alternatives Considered

1. **Content hash as ID:** Rejected — changes on every edit, breaks lineage tracking
2. **Name+version as ID:** Rejected — not globally unique, conflicts possible
3. **Incrementing integer:** Rejected — requires central coordination

## References

- [RFC 4122 - UUID](https://tools.ietf.org/html/rfc4122)
- yapCAD package format: `src/yapcad/package/core.py`
- Web-viewer spec: `projects/YAPCAD-WEBVIEWER-SPEC.md`
