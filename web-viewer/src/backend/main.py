"""yapCAD Package Server - FastAPI Backend."""
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import Response, FileResponse
from fastapi.staticfiles import StaticFiles
from pathlib import Path
import zipfile
import json
from typing import Optional

from .models import (
    PackageManifest, PackageListItem, PackageListResponse,
    PackageDetail, GeometryFileInfo, ErrorResponse
)

PACKAGES_DIR = Path("./packages")

app = FastAPI(title="yapCAD Package Server", version="0.1.0")

# CORS middleware (allow all origins for dev)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["GET"],
    allow_headers=["*"],
)


def get_package_path(package_id: str) -> Path:
    """Get path to .ycpkg file, raise 404 if not found."""
    path = PACKAGES_DIR / f"{package_id}.ycpkg"
    if not path.exists():
        raise HTTPException(status_code=404, detail=f"Package not found: {package_id}")
    return path


def read_manifest(zf: zipfile.ZipFile) -> PackageManifest:
    """Read and parse manifest.json from open zipfile."""
    data = json.loads(zf.read("manifest.json"))
    return PackageManifest.model_validate(data)


def count_geometry_files(zf: zipfile.ZipFile) -> int:
    """Count geometry JSON files in the package."""
    return sum(1 for f in zf.namelist() if f.startswith("geometry/") and f.endswith(".json"))


@app.get("/api/packages", response_model=PackageListResponse)
async def list_packages() -> PackageListResponse:
    """List all available packages."""
    items = []
    for path in PACKAGES_DIR.glob("*.ycpkg"):
        package_id = path.stem
        try:
            with zipfile.ZipFile(path, "r") as zf:
                manifest = read_manifest(zf)
                item = PackageListItem(
                    id=package_id,
                    name=manifest.name,
                    version=manifest.version,
                    description=manifest.description,
                    thumbnail_url=f"/api/packages/{package_id}/thumbnail" if manifest.thumbnail else None,
                    geometry_count=count_geometry_files(zf),
                )
                items.append(item)
        except Exception as e:
            print(f"Error reading {path}: {e}")
    return PackageListResponse(packages=items, total=len(items))


@app.get("/api/packages/{package_id}", response_model=PackageDetail)
async def get_package_details(package_id: str) -> PackageDetail:
    """Get detailed information about a specific package."""
    path = get_package_path(package_id)
    try:
        with zipfile.ZipFile(path, "r") as zf:
            manifest = read_manifest(zf)
            
            # Build geometry file info list
            geometry_files = []
            for info in zf.infolist():
                if info.filename.startswith("geometry/") and info.filename.endswith(".json"):
                    rel_path = info.filename[len("geometry/"):]  # Remove geometry/ prefix
                    geometry_files.append(GeometryFileInfo(
                        path=info.filename,
                        url=f"/api/packages/{package_id}/geometry/{rel_path}",
                        size_bytes=info.file_size,
                    ))
            
            return PackageDetail(
                id=package_id,
                manifest=manifest,
                geometry_files=geometry_files,
                material_files=manifest.material_files,
            )
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error reading package: {e}")


@app.get("/api/packages/{package_id}/geometry/{file_path:path}")
async def get_geometry_file(package_id: str, file_path: str) -> Response:
    """Get a specific geometry JSON file from a package."""
    path = get_package_path(package_id)
    try:
        with zipfile.ZipFile(path, "r") as zf:
            full_path = f"geometry/{file_path}"
            if full_path not in zf.namelist():
                raise HTTPException(
                    status_code=404,
                    detail=f"Geometry file not found: {file_path}"
                )
            content = zf.read(full_path)
            return Response(content=content, media_type="application/json")
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error reading geometry file: {e}")


@app.get("/api/packages/{package_id}/thumbnail")
async def get_thumbnail(package_id: str) -> Response:
    """Get the thumbnail image for a package."""
    path = get_package_path(package_id)
    try:
        with zipfile.ZipFile(path, "r") as zf:
            manifest = read_manifest(zf)
            if not manifest.thumbnail:
                raise HTTPException(
                    status_code=404,
                    detail=f"Package has no thumbnail"
                )
            content = zf.read(manifest.thumbnail)
            
            # Determine media type from extension
            suffix = Path(manifest.thumbnail).suffix.lower()
            media_types = {
                ".png": "image/png",
                ".jpg": "image/jpeg",
                ".jpeg": "image/jpeg",
                ".gif": "image/gif",
                ".webp": "image/webp",
            }
            media_type = media_types.get(suffix, "application/octet-stream")
            
            return Response(content=content, media_type=media_type)
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error reading thumbnail: {e}")


@app.get("/health")
async def health_check():
    """Health check endpoint."""
    return {"status": "ok", "packages_dir": str(PACKAGES_DIR.absolute())}


# Serve static files from public/ directory
PUBLIC_DIR = Path(__file__).parent.parent.parent / "public"

@app.get("/")
async def serve_index():
    """Serve the main viewer page."""
    index_path = PUBLIC_DIR / "index.html"
    if not index_path.exists():
        raise HTTPException(status_code=404, detail="index.html not found")
    return FileResponse(index_path, media_type="text/html")
