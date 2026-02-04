"""
3D Text generation for yapCAD.

This module provides functions for creating 3D text suitable for labeling
3D printed parts. It includes both extruded (raised) text and engraved
(cut-in) text capabilities using a simple, printable block font.

Example usage:

    from yapcad.text3d import text_solid, engrave_text
    from yapcad.io.stl import write_stl

    # Create raised text
    label = text_solid("ROBOT", height=10.0, depth=2.0)
    write_stl(label, "robot_label.stl")

    # Engrave text into a surface
    plate = prism(50, 30, 5)
    engraved = engrave_text(plate, "V1.0", point(0, 0, 2.5),
                            vect(0, 0, 1, 0), height=4.0, depth=0.5)
"""

from yapcad.geom import point, vect, epsilon, add, sub, scale3
from yapcad.geom3d import (
    poly2surfaceXY, solid, solid_boolean, issolid,
    translatesurface, rotatesurface
)
from yapcad.geom3d_util import extrude, prism
from yapcad.xform import Translation, Rotation, Matrix
from copy import deepcopy
import math

# Grid dimensions for block font (5 wide x 7 tall)
CHAR_WIDTH = 5
CHAR_HEIGHT = 7
STROKE_WIDTH = 1.0

# Block font definition: each character is a list of rectangles
# Rectangle format: (x, y, width, height) in grid units
# Origin is bottom-left of character cell

BLOCK_FONT = {
    'A': [
        (1, 6, 3, 1),   # Top bar
        (0, 0, 1, 6),   # Left vertical
        (4, 0, 1, 6),   # Right vertical
        (1, 3, 3, 1),   # Middle bar
    ],
    'B': [
        (0, 0, 1, 7),   # Left vertical
        (1, 6, 3, 1),   # Top bar
        (1, 3, 3, 1),   # Middle bar
        (1, 0, 3, 1),   # Bottom bar
        (4, 4, 1, 2),   # Upper right
        (4, 1, 1, 2),   # Lower right
    ],
    'C': [
        (0, 1, 1, 5),   # Left vertical
        (1, 6, 4, 1),   # Top bar
        (1, 0, 4, 1),   # Bottom bar
    ],
    'D': [
        (0, 0, 1, 7),   # Left vertical
        (1, 6, 3, 1),   # Top bar
        (1, 0, 3, 1),   # Bottom bar
        (4, 1, 1, 5),   # Right vertical
    ],
    'E': [
        (0, 0, 1, 7),   # Left vertical
        (1, 6, 4, 1),   # Top bar
        (1, 3, 3, 1),   # Middle bar
        (1, 0, 4, 1),   # Bottom bar
    ],
    'F': [
        (0, 0, 1, 7),   # Left vertical
        (1, 6, 4, 1),   # Top bar
        (1, 3, 3, 1),   # Middle bar
    ],
    'G': [
        (0, 1, 1, 5),   # Left vertical
        (1, 6, 4, 1),   # Top bar
        (1, 0, 4, 1),   # Bottom bar
        (4, 1, 1, 3),   # Lower right
        (2, 3, 2, 1),   # Middle bar (partial)
    ],
    'H': [
        (0, 0, 1, 7),   # Left vertical
        (4, 0, 1, 7),   # Right vertical
        (1, 3, 3, 1),   # Middle bar
    ],
    'I': [
        (0, 6, 5, 1),   # Top bar
        (0, 0, 5, 1),   # Bottom bar
        (2, 1, 1, 5),   # Center vertical
    ],
    'J': [
        (0, 6, 5, 1),   # Top bar
        (3, 1, 1, 5),   # Right vertical
        (0, 0, 3, 1),   # Bottom bar
        (0, 1, 1, 2),   # Lower left hook
    ],
    'K': [
        (0, 0, 1, 7),   # Left vertical
        (1, 3, 1, 1),   # Middle connector
        (2, 4, 1, 1),   # Upper diagonal part 1
        (3, 5, 1, 1),   # Upper diagonal part 2
        (4, 6, 1, 1),   # Upper diagonal part 3
        (2, 2, 1, 1),   # Lower diagonal part 1
        (3, 1, 1, 1),   # Lower diagonal part 2
        (4, 0, 1, 1),   # Lower diagonal part 3
    ],
    'L': [
        (0, 0, 1, 7),   # Left vertical
        (1, 0, 4, 1),   # Bottom bar
    ],
    'M': [
        (0, 0, 1, 7),   # Left vertical
        (4, 0, 1, 7),   # Right vertical
        (1, 5, 1, 1),   # Left diagonal
        (2, 4, 1, 1),   # Center peak
        (3, 5, 1, 1),   # Right diagonal
    ],
    'N': [
        (0, 0, 1, 7),   # Left vertical
        (4, 0, 1, 7),   # Right vertical
        (1, 5, 1, 1),   # Diagonal part 1
        (2, 4, 1, 1),   # Diagonal part 2
        (3, 3, 1, 1),   # Diagonal part 3
    ],
    'O': [
        (0, 1, 1, 5),   # Left vertical
        (4, 1, 1, 5),   # Right vertical
        (1, 6, 3, 1),   # Top bar
        (1, 0, 3, 1),   # Bottom bar
    ],
    'P': [
        (0, 0, 1, 7),   # Left vertical
        (1, 6, 3, 1),   # Top bar
        (1, 3, 3, 1),   # Middle bar
        (4, 4, 1, 2),   # Right upper
    ],
    'Q': [
        (0, 1, 1, 5),   # Left vertical
        (4, 2, 1, 4),   # Right vertical
        (1, 6, 3, 1),   # Top bar
        (1, 0, 3, 1),   # Bottom bar
        (3, 1, 1, 1),   # Tail part 1
        (4, 0, 1, 1),   # Tail part 2
    ],
    'R': [
        (0, 0, 1, 7),   # Left vertical
        (1, 6, 3, 1),   # Top bar
        (1, 3, 3, 1),   # Middle bar
        (4, 4, 1, 2),   # Right upper
        (2, 2, 1, 1),   # Lower diagonal part 1
        (3, 1, 1, 1),   # Lower diagonal part 2
        (4, 0, 1, 1),   # Lower diagonal part 3
    ],
    'S': [
        (1, 6, 4, 1),   # Top bar
        (0, 4, 1, 2),   # Upper left
        (1, 3, 3, 1),   # Middle bar
        (4, 1, 1, 2),   # Lower right
        (0, 0, 4, 1),   # Bottom bar
    ],
    'T': [
        (0, 6, 5, 1),   # Top bar
        (2, 0, 1, 6),   # Center vertical
    ],
    'U': [
        (0, 1, 1, 6),   # Left vertical
        (4, 1, 1, 6),   # Right vertical
        (1, 0, 3, 1),   # Bottom bar
    ],
    'V': [
        (0, 3, 1, 4),   # Left vertical upper
        (1, 1, 1, 2),   # Left diagonal
        (2, 0, 1, 1),   # Bottom point
        (3, 1, 1, 2),   # Right diagonal
        (4, 3, 1, 4),   # Right vertical upper
    ],
    'W': [
        (0, 0, 1, 7),   # Left vertical
        (4, 0, 1, 7),   # Right vertical
        (1, 1, 1, 1),   # Left diagonal
        (2, 2, 1, 1),   # Center trough
        (3, 1, 1, 1),   # Right diagonal
    ],
    'X': [
        (0, 5, 1, 2),   # Top left
        (1, 4, 1, 1),   # Upper left diagonal
        (2, 3, 1, 1),   # Center
        (3, 4, 1, 1),   # Upper right diagonal
        (4, 5, 1, 2),   # Top right
        (1, 2, 1, 1),   # Lower left diagonal
        (0, 0, 1, 2),   # Bottom left
        (3, 2, 1, 1),   # Lower right diagonal
        (4, 0, 1, 2),   # Bottom right
    ],
    'Y': [
        (0, 5, 1, 2),   # Top left
        (1, 4, 1, 1),   # Upper left diagonal
        (2, 0, 1, 4),   # Center vertical
        (3, 4, 1, 1),   # Upper right diagonal
        (4, 5, 1, 2),   # Top right
    ],
    'Z': [
        (0, 6, 5, 1),   # Top bar
        (3, 5, 1, 1),   # Diagonal part 1
        (2, 3, 1, 2),   # Diagonal part 2
        (1, 2, 1, 1),   # Diagonal part 3
        (0, 0, 5, 1),   # Bottom bar
    ],
    '0': [
        (0, 1, 1, 5),   # Left vertical
        (4, 1, 1, 5),   # Right vertical
        (1, 6, 3, 1),   # Top bar
        (1, 0, 3, 1),   # Bottom bar
        (2, 3, 1, 1),   # Center diagonal (distinguishes from O)
    ],
    '1': [
        (2, 0, 1, 7),   # Center vertical
        (1, 5, 1, 1),   # Top serif
        (1, 0, 3, 1),   # Bottom bar
    ],
    '2': [
        (0, 5, 1, 1),   # Top left
        (1, 6, 3, 1),   # Top bar
        (4, 4, 1, 2),   # Upper right
        (2, 3, 2, 1),   # Middle bar
        (1, 2, 1, 1),   # Diagonal
        (0, 1, 1, 1),   # Lower left
        (0, 0, 5, 1),   # Bottom bar
    ],
    '3': [
        (0, 6, 4, 1),   # Top bar
        (4, 4, 1, 2),   # Upper right
        (1, 3, 3, 1),   # Middle bar
        (4, 1, 1, 2),   # Lower right
        (0, 0, 4, 1),   # Bottom bar
    ],
    '4': [
        (0, 3, 1, 4),   # Left vertical (upper)
        (4, 0, 1, 7),   # Right vertical (full)
        (1, 3, 3, 1),   # Horizontal bar
    ],
    '5': [
        (0, 6, 5, 1),   # Top bar
        (0, 3, 1, 3),   # Upper left
        (1, 3, 3, 1),   # Middle bar
        (4, 1, 1, 2),   # Lower right
        (0, 0, 4, 1),   # Bottom bar
    ],
    '6': [
        (0, 1, 1, 5),   # Left vertical
        (1, 6, 4, 1),   # Top bar
        (1, 3, 3, 1),   # Middle bar
        (1, 0, 3, 1),   # Bottom bar
        (4, 1, 1, 2),   # Lower right
    ],
    '7': [
        (0, 6, 5, 1),   # Top bar
        (4, 3, 1, 3),   # Upper right
        (2, 0, 2, 3),   # Lower center-right
    ],
    '8': [
        (0, 1, 1, 2),   # Lower left
        (0, 4, 1, 2),   # Upper left
        (4, 1, 1, 2),   # Lower right
        (4, 4, 1, 2),   # Upper right
        (1, 6, 3, 1),   # Top bar
        (1, 3, 3, 1),   # Middle bar
        (1, 0, 3, 1),   # Bottom bar
    ],
    '9': [
        (0, 4, 1, 2),   # Upper left
        (4, 1, 1, 5),   # Right vertical
        (1, 6, 3, 1),   # Top bar
        (1, 3, 3, 1),   # Middle bar
        (0, 0, 4, 1),   # Bottom bar
    ],
    '-': [
        (1, 3, 3, 1),   # Horizontal bar
    ],
    '_': [
        (0, 0, 5, 1),   # Underscore at bottom
    ],
    '.': [
        (2, 0, 1, 1),   # Single dot
    ],
    ':': [
        (2, 4, 1, 1),   # Upper dot
        (2, 1, 1, 1),   # Lower dot
    ],
    '/': [
        (0, 0, 1, 2),   # Bottom left
        (1, 2, 1, 1),   # Lower middle
        (2, 3, 1, 1),   # Center
        (3, 4, 1, 1),   # Upper middle
        (4, 5, 1, 2),   # Top right
    ],
    '(': [
        (3, 5, 1, 1),   # Top curve
        (2, 2, 1, 3),   # Vertical
        (3, 1, 1, 1),   # Bottom curve
    ],
    ')': [
        (1, 5, 1, 1),   # Top curve
        (2, 2, 1, 3),   # Vertical
        (1, 1, 1, 1),   # Bottom curve
    ],
    ' ': [],  # Space - no rectangles
}


def _rect_to_polygon(x, y, width, height, scale, x_offset):
    """Convert a rectangle to a closed polygon.

    Args:
        x, y: Bottom-left corner in grid units
        width, height: Dimensions in grid units
        scale: Scale factor (mm per grid unit)
        x_offset: X offset in mm for character positioning

    Returns:
        List of points forming a closed polygon (CCW winding)
    """
    x0 = (x * scale) + x_offset
    y0 = y * scale
    x1 = x0 + (width * scale)
    y1 = y0 + (height * scale)

    # Return closed polygon (CCW winding for positive area)
    return [
        point(x0, y0, 0),
        point(x1, y0, 0),
        point(x1, y1, 0),
        point(x0, y1, 0),
        point(x0, y0, 0),  # Close the polygon
    ]


def text_to_polygons(text, height=5.0, spacing=1.0):
    """Convert text string to list of 2D polygons.

    Each character becomes one or more closed polygons suitable for
    extrusion or boolean operations.

    Args:
        text: String to convert (supports A-Z, 0-9, basic punctuation)
        height: Character height in mm (default 5.0)
        spacing: Space between characters as fraction of char width (default 1.0)

    Returns:
        List of closed polygons (each polygon is a list of points)
    """
    polygons = []
    scale = height / CHAR_HEIGHT
    char_width = CHAR_WIDTH * scale
    gap = spacing * scale  # Gap between characters

    x_offset = 0.0

    for char in text.upper():
        if char in BLOCK_FONT:
            rects = BLOCK_FONT[char]
            for rect in rects:
                x, y, w, h = rect
                poly = _rect_to_polygon(x, y, w, h, scale, x_offset)
                polygons.append(poly)
        else:
            # Unknown character - draw a small placeholder rectangle
            poly = _rect_to_polygon(1, 1, 3, 5, scale, x_offset)
            polygons.append(poly)

        x_offset += char_width + gap

    return polygons


def text_solid(text, height=5.0, depth=1.0, spacing=1.0, center=None):
    """Create extruded 3D text solid.

    Generates a solid from the given text string by extruding each
    character polygon in the Z direction.

    Args:
        text: String to render (supports A-Z, 0-9, basic punctuation)
        height: Character height in mm (default 5.0)
        depth: Extrusion depth in mm (default 1.0)
        spacing: Character spacing factor (default 1.0)
        center: Optional center point for the text block

    Returns:
        yapCAD solid (combined extruded text)
    """
    polygons = text_to_polygons(text, height, spacing)

    if not polygons:
        # Return empty solid for empty text
        return solid([], [], ['procedure', f'text_solid("{text}")'])

    # Extrude each polygon and combine
    all_surfaces = []
    for poly in polygons:
        try:
            # Create surface from polygon
            # poly2surfaceXY returns (surface, boundary_indices) tuple
            surf_result = poly2surfaceXY(poly)
            if isinstance(surf_result, tuple):
                surf = surf_result[0]
            else:
                surf = surf_result
            # Extrude upward
            extruded = extrude(surf, depth, vect(0, 0, 1, 0))
            # Extract surfaces from the extruded solid
            if issolid(extruded):
                all_surfaces.extend(extruded[1])
        except (ValueError, IndexError) as e:
            # Skip degenerate polygons
            continue

    if not all_surfaces:
        return solid([], [], ['procedure', f'text_solid("{text}")'])

    result = solid(all_surfaces, [], ['procedure', f'text_solid("{text}", height={height}, depth={depth})'])

    # Center the text if requested
    if center is not None:
        # Calculate current bounding box
        from yapcad.geom3d import solidbbox
        bbox = solidbbox(result)
        if bbox:
            current_center = point(
                (bbox[0][0] + bbox[1][0]) / 2,
                (bbox[0][1] + bbox[1][1]) / 2,
                (bbox[0][2] + bbox[1][2]) / 2
            )
            offset = sub(center, current_center)
            # Translate all surfaces
            for i, surf in enumerate(result[1]):
                result[1][i] = translatesurface(surf, offset)

    return result


def text_width(text, height=5.0, spacing=1.0):
    """Calculate the total width of rendered text.

    Useful for positioning and centering text.

    Args:
        text: String to measure
        height: Character height in mm (default 5.0)
        spacing: Character spacing factor (default 1.0)

    Returns:
        Total width in mm
    """
    if not text:
        return 0.0

    scale = height / CHAR_HEIGHT
    char_width = CHAR_WIDTH * scale
    gap = spacing * scale

    # Width = (n chars * char_width) + ((n-1) gaps)
    n = len(text)
    return (n * char_width) + ((n - 1) * gap)


def engrave_text(target_solid, text, position, normal, height=3.0, depth=0.5, spacing=1.0):
    """Cut text into a face of a solid (boolean difference).

    Creates text as a solid and subtracts it from the target solid to
    create engraved (cut-in) text.

    Args:
        target_solid: yapCAD solid to engrave into
        text: Text string to engrave
        position: Point on face where text starts (bottom-left of first char)
        normal: Face normal vector (determines orientation of text)
        height: Character height in mm (default 3.0)
        depth: Engraving depth in mm (default 0.5)
        spacing: Character spacing factor (default 1.0)

    Returns:
        New solid with text engraved into the specified face
    """
    if not issolid(target_solid):
        raise ValueError("target_solid must be a valid yapCAD solid")

    # Create text solid
    text_sld = text_solid(text, height, depth + 0.1, spacing)  # Extra depth for clean cut

    if not text_sld[1]:  # No surfaces (empty text)
        return deepcopy(target_solid)

    # Normalize the normal vector
    n = [normal[0], normal[1], normal[2]]
    nm = math.sqrt(n[0]**2 + n[1]**2 + n[2]**2)
    if nm < epsilon:
        raise ValueError("normal vector has zero length")
    n = [n[0]/nm, n[1]/nm, n[2]/nm]

    # Calculate transformation to align text with face
    # Text is generated in XY plane with Z-up
    # We need to rotate so Z aligns with the face normal

    # Build rotation to align Z-axis with normal
    z_axis = [0, 0, 1]

    # Calculate rotation axis (cross product of z_axis and normal)
    rx = z_axis[1] * n[2] - z_axis[2] * n[1]
    ry = z_axis[2] * n[0] - z_axis[0] * n[2]
    rz = z_axis[0] * n[1] - z_axis[1] * n[0]
    rot_axis_mag = math.sqrt(rx**2 + ry**2 + rz**2)

    # Calculate rotation angle
    dot_prod = z_axis[0] * n[0] + z_axis[1] * n[1] + z_axis[2] * n[2]
    angle = math.acos(max(-1, min(1, dot_prod)))  # Clamp for numerical stability

    # Transform text solid
    transformed_surfaces = []

    for surf in text_sld[1]:
        new_surf = deepcopy(surf)

        # Rotate vertices if needed
        if rot_axis_mag > epsilon and abs(angle) > epsilon:
            rot_axis = point(rx / rot_axis_mag, ry / rot_axis_mag, rz / rot_axis_mag)
            angle_deg = math.degrees(angle)
            rot_matrix = Rotation(angle_deg, cent=point(0, 0, 0), axis=rot_axis)

            # Transform vertices
            new_verts = []
            for v in new_surf[1]:
                new_v = rot_matrix.mul(v)
                new_verts.append(new_v)
            new_surf[1] = new_verts

            # Transform normals
            new_norms = []
            for nm_vec in new_surf[2]:
                new_nm = rot_matrix.mul(nm_vec)
                new_norms.append(new_nm)
            new_surf[2] = new_norms

        # Translate to position (slightly inside the face for clean cut)
        offset = sub(position, scale3(n, depth * 0.5))
        new_surf = translatesurface(new_surf, offset)
        transformed_surfaces.append(new_surf)

    # Create transformed text solid
    engraving_solid = solid(
        transformed_surfaces,
        [],
        ['procedure', f'engrave_text("{text}")']
    )

    # Perform boolean difference
    try:
        result = solid_boolean(target_solid, engraving_solid, 'difference')
        return result
    except Exception as e:
        # If boolean fails, return original solid with warning
        print(f"Warning: engraving failed ({e}), returning original solid")
        return deepcopy(target_solid)


def get_supported_characters():
    """Return a string of all supported characters.

    Returns:
        String containing all characters that can be rendered
    """
    return ''.join(sorted(BLOCK_FONT.keys()))
