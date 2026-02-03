"""Pydantic models for yapCAD Package Server API."""
from datetime import datetime
from typing import List, Optional
from pydantic import BaseModel, Field

__all__ = [
    "PackageManifest",
    "PackageListItem",
    "PackageDetail",
    "GeometryFileInfo",
    "PackageListResponse",
    "ErrorResponse",
]


class PackageManifest(BaseModel):
    """Metadata from manifest.json."""
    name: str
    version: str
    schema_version: str = Field(default="1.0")
    description: Optional[str] = None
    author: Optional[str] = None
    created_at: Optional[datetime] = None
    updated_at: Optional[datetime] = None
    tags: List[str] = Field(default_factory=list)
    geometry_files: List[str] = Field(default_factory=list)
    material_files: List[str] = Field(default_factory=list)
    thumbnail: Optional[str] = None

    model_config = {"title": "PackageManifest"}


class GeometryFileInfo(BaseModel):
    """Info about a geometry file within a package."""
    path: str
    url: str
    entity_count: Optional[int] = None
    size_bytes: Optional[int] = None

    model_config = {"title": "GeometryFileInfo"}


class PackageListItem(BaseModel):
    """Summary for package listing."""
    id: str  # Unique identifier (filename without .ycpkg)
    name: str
    version: str
    description: Optional[str] = None
    thumbnail_url: Optional[str] = None
    geometry_count: int

    model_config = {"title": "PackageListItem"}


class PackageDetail(BaseModel):
    """Full package information."""
    id: str
    manifest: PackageManifest
    geometry_files: List[GeometryFileInfo]
    material_files: List[str] = Field(default_factory=list)

    model_config = {"title": "PackageDetail"}


class PackageListResponse(BaseModel):
    """Response for GET /api/packages."""
    packages: List[PackageListItem]
    total: int

    model_config = {"title": "PackageListResponse"}


class ErrorResponse(BaseModel):
    """Standard error response."""
    error: str
    detail: Optional[str] = None

    model_config = {"title": "ErrorResponse"}
