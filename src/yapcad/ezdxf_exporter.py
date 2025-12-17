"""
DXF export utilities for yapCAD geometry.

Provides functions to export 2D geometry (region2d, curves, polylines) to DXF format
using the ezdxf library via yapCAD's ezdxfDraw drawable class.

Copyright (c) 2024 yapCAD contributors
All rights reserved (MIT License)
"""

from pathlib import Path
from typing import Union, List, Any

from yapcad.ezdxf_drawable import ezdxfDraw
from yapcad.geom import (
    ispoint, isline, isarc, isellipse, ispoly, isgeomlist,
)


def is_2d_geometry(geom) -> bool:
    """Check if geometry is 2D (exportable to DXF as native entities).

    Returns True for:
    - Points, lines, arcs, circles, ellipses
    - Polylines/polygons
    - Geometry lists (lists of the above)
    - Catmull-Rom and NURBS splines

    Returns False for:
    - 3D solids
    - Surfaces
    - Unknown types
    """
    if geom is None:
        return False

    # Check for spline types
    from yapcad.geom import iscatmullrom, isnurbs

    # Simple geometric primitives
    if ispoint(geom) or isline(geom) or isarc(geom) or isellipse(geom):
        return True

    # Splines
    if iscatmullrom(geom) or isnurbs(geom):
        return True

    # Polylines
    if ispoly(geom):
        return True

    # Geometry lists (region2d, path2d, etc.)
    if isgeomlist(geom):
        return True

    # Check if it's a list of 2D geometry
    if isinstance(geom, list) and geom:
        # Check first element - if it's 2D geometry, assume the whole list is
        return is_2d_geometry(geom[0])

    return False


def write_dxf(geometry, output_path: str, layer: str = 'PATHS') -> bool:
    """Export 2D geometry to a DXF file.

    Args:
        geometry: yapCAD 2D geometry (region2d, curves, polylines, etc.)
        output_path: Path to output DXF file (with or without .dxf extension)
        layer: DXF layer name (default 'PATHS')

    Returns:
        True if export succeeded, False otherwise.

    Supported geometry types:
    - point: Exported as point entity
    - line: Exported as LINE entity
    - arc: Exported as ARC or CIRCLE entity
    - ellipse: Exported as ELLIPSE entity
    - polyline/polygon: Exported as series of LINE entities
    - catmullrom: Sampled and exported as polyline
    - nurbs: Sampled and exported as polyline
    - geomlist (region2d): All elements exported
    """
    # Normalize output path - remove .dxf extension if present
    # (ezdxfDraw.display() adds it automatically)
    path = Path(output_path)
    if path.suffix.lower() == '.dxf':
        filename = str(path.with_suffix(''))
    else:
        filename = str(path)

    # Create drawable
    drawable = ezdxfDraw()
    drawable.filename = filename
    drawable.layer = layer

    # Draw the geometry
    try:
        drawable.draw(geometry)
        drawable.display()
        return True
    except Exception as e:
        print(f"DXF export error: {e}")
        return False


def write_dxf_multi(geometries: List[Any], output_path: str,
                    layers: Union[str, List[str]] = 'PATHS') -> bool:
    """Export multiple geometry objects to a single DXF file.

    Args:
        geometries: List of yapCAD 2D geometry objects
        output_path: Path to output DXF file
        layers: Layer name or list of layer names (one per geometry)

    Returns:
        True if export succeeded, False otherwise.
    """
    path = Path(output_path)
    if path.suffix.lower() == '.dxf':
        filename = str(path.with_suffix(''))
    else:
        filename = str(path)

    drawable = ezdxfDraw()
    drawable.filename = filename

    # Normalize layers to list
    if isinstance(layers, str):
        layer_list = [layers] * len(geometries)
    else:
        layer_list = layers
        if len(layer_list) < len(geometries):
            # Extend with last layer name
            layer_list = layer_list + [layer_list[-1]] * (len(geometries) - len(layer_list))

    try:
        for geom, layer in zip(geometries, layer_list):
            drawable.layer = layer
            drawable.draw(geom)
        drawable.display()
        return True
    except Exception as e:
        print(f"DXF export error: {e}")
        return False
