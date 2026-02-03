# Task 2a: Pydantic Models for yapCAD Package Server

## Context

We're building a FastAPI backend to serve `.ycpkg` packages. This task creates the Pydantic models that define the API's data structures.

## Reference: yapCAD Package Format

A `.ycpkg` file is a ZIP archive containing:
```
manifest.json          # Package metadata
geometry/primary.json  # Main geometry (yapcad-geometry-json-v0.1)
geometry/*.json        # Optional additional geometry files
materials/*.json       # Optional material definitions
```

## Required Models

### 1. PackageManifest

Based on yapCAD's `package/core.py`. Fields:

```python
class PackageManifest(BaseModel):
    """Metadata from manifest.json"""
    name: str                           # Package name
    version: str                        # Semantic version
    schema_version: str = "1.0"         # Manifest schema version
    description: Optional[str] = None
    author: Optional[str] = None
    created_at: Optional[datetime] = None
    updated_at: Optional[datetime] = None
    tags: List[str] = []
    geometry_files: List[str] = []      # Paths within the package
    material_files: List[str] = []
    thumbnail: Optional[str] = None     # Path to thumbnail image
```

### 2. PackageListItem

For listing available packages:

```python
class PackageListItem(BaseModel):
    """Summary for package listing"""
    id: str                             # Unique identifier (filename without .ycpkg)
    name: str
    version: str
    description: Optional[str] = None
    thumbnail_url: Optional[str] = None # URL to fetch thumbnail
    geometry_count: int                 # Number of geometry files
```

### 3. PackageDetail

Full package info for detail view:

```python
class PackageDetail(BaseModel):
    """Full package information"""
    id: str
    manifest: PackageManifest
    geometry_files: List[GeometryFileInfo]
    material_files: List[str]
```

### 4. GeometryFileInfo

```python
class GeometryFileInfo(BaseModel):
    """Info about a geometry file within a package"""
    path: str                           # Path within package
    url: str                            # URL to fetch content
    entity_count: Optional[int] = None  # Number of entities if known
    size_bytes: Optional[int] = None
```

### 5. API Response Wrappers

```python
class PackageListResponse(BaseModel):
    """Response for GET /api/packages"""
    packages: List[PackageListItem]
    total: int

class ErrorResponse(BaseModel):
    """Standard error response"""
    error: str
    detail: Optional[str] = None
```

## Output

Create file: `src/backend/models.py`

Requirements:
- Use Pydantic v2 syntax (model_config, Field())
- Add docstrings to all models
- Use Optional[] for nullable fields
- Export all models via __all__

## Example Usage

```python
from models import PackageManifest, PackageListItem

# Parse manifest from JSON
manifest = PackageManifest.model_validate_json(manifest_json)

# Create list item
item = PackageListItem(
    id="demo-cube",
    name=manifest.name,
    version=manifest.version,
    geometry_count=len(manifest.geometry_files)
)
```
