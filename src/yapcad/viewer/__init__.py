"""
yapCAD Viewer Package
=====================

A VTK-based 3D assembly viewer with REST API and WebSocket support.

This package provides a general-purpose STL viewer for visualizing assemblies
with support for:

- Multi-viewport rendering (ISO, TOP, FRONT, SIDE)
- Real-time part positioning via 4x4 transform matrices
- X-ray mode for internal inspection
- REST API for programmatic control
- WebSocket events for real-time updates
- Optional file-based command interface for backward compatibility

Example Usage
-------------

Basic viewer with positions file::

    from yapcad.viewer import VTKViewer, ViewerConfig

    config = ViewerConfig(
        stl_dir="/path/to/stl/files",
        positions_file="/path/to/positions.json"
    )
    viewer = VTKViewer(config)
    viewer.load_from_json()
    viewer.start()

With REST API server::

    from yapcad.viewer import VTKViewer, ViewerAPIServer, ViewerConfig

    config = ViewerConfig(stl_dir="/path/to/stls")
    viewer = VTKViewer(config)
    server = ViewerAPIServer(viewer)
    server.run(port=5000)  # Viewer runs in background thread

Dependencies
------------

Required:
    - vtk >= 9.0

Optional (for API server):
    - flask >= 2.0
    - flask-socketio >= 5.0
    - flask-cors >= 3.0

Module Contents
---------------

- :class:`VTKViewer` - Core VTK-based viewer with multi-viewport support
- :class:`ViewerAPIServer` - Flask REST API and SocketIO WebSocket server
- :class:`ViewerConfig` - Configuration dataclass
- :mod:`events` - WebSocket event type definitions

Viewer system contributed by Jeremy Mika.
"""

__version__ = "0.1.0"

import importlib
import sys

# Cache for lazy-loaded modules and classes
_cache = {}


def _check_vtk():
    """Check if VTK is available."""
    try:
        import vtk  # noqa: F401
        return True
    except ImportError:
        return False


def _check_flask():
    """Check if Flask and dependencies are available."""
    try:
        import flask  # noqa: F401
        import flask_socketio  # noqa: F401
        import flask_cors  # noqa: F401
        return True
    except ImportError:
        return False


def __getattr__(name):
    """Lazy import for viewer components."""
    if name in _cache:
        return _cache[name]

    if name == "VTKViewer":
        if not _check_vtk():
            raise ImportError(
                "VTKViewer requires VTK >= 9.0. "
                "Install with: pip install vtk"
            )
        from .vtk_viewer import VTKViewer
        _cache[name] = VTKViewer
        return VTKViewer

    elif name == "ViewerAPIServer":
        if not _check_vtk():
            raise ImportError(
                "ViewerAPIServer requires VTK >= 9.0. "
                "Install with: pip install vtk"
            )
        if not _check_flask():
            raise ImportError(
                "ViewerAPIServer requires Flask and SocketIO. "
                "Install with: pip install flask flask-socketio flask-cors"
            )
        from .api_server import ViewerAPIServer
        _cache[name] = ViewerAPIServer
        return ViewerAPIServer

    elif name == "ViewerConfig":
        from .config import ViewerConfig
        _cache[name] = ViewerConfig
        return ViewerConfig

    elif name == "events":
        # Use importlib to avoid recursion
        events_module = importlib.import_module(".events", __name__)
        _cache[name] = events_module
        return events_module

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [
    "VTKViewer",
    "ViewerAPIServer",
    "ViewerConfig",
    "events",
    "__version__",
]
