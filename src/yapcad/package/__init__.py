"""Public API for yapCAD `.ycpkg` package workflows."""

from __future__ import annotations

from .core import (
    MANIFEST_FILENAME,
    PACKAGE_SCHEMA,
    PackageManifest,
    create_package_from_entities,
    load_geometry,
)
from .validator import validate_package


def view_package(package_path, *, strict: bool = False):
    """Import and invoke the interactive viewer lazily."""

    from .viewer import view_package as _view  # local import to avoid pyglet init during tests

    return _view(package_path, strict=strict)

__all__ = [
    "PACKAGE_SCHEMA",
    "MANIFEST_FILENAME",
    "PackageManifest",
    "create_package_from_entities",
    "load_geometry",
    "validate_package",
    "view_package",
]
