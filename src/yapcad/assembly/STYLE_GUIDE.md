# yapCAD Assembly System Style Guide

This document captures the design patterns, conventions, and architectural principles from yapCAD's core modules to guide the implementation of the assembly/mate system.

**Analysis Date:** 2026-02-04
**Analyzed Modules:** `geom.py`, `geom3d.py`, `xform.py`, `metadata.py`, `brep.py`

---

## 1. Data Representation Philosophy

### 1.1 Functional, List-Based Structures

yapCAD uses **functional data structures** represented as Python lists with **tagged tuples**:

```python
# Points: [x, y, z, w] - homogeneous coordinates
point = [0.0, 0.0, 0.0, 1.0]

# Lines: [[x1,y1,z1,w], [x2,y2,z2,w]]
line = [[0,0,0,1], [10,0,0,1]]

# Surfaces: ['surface', vertices, normals, faces, boundary, holes, metadata?]
surface = ['surface', vrts, nrms, faces, bnd, holes]

# Solids: ['solid', surfaces, material, construction, metadata?]
solid = ['solid', surfs, mat, const]
```

**Key Principle:** The first element is a **string tag** identifying the type. This enables type checking without classes.

### 1.2 Homogeneous Coordinates (Projective Geometry)

All vectors are **4-tuples `[x, y, z, w]`**:

- **Points:** `w=1` (e.g., `[5.0, 3.0, 2.0, 1.0]`)
- **Direction vectors:** `w=0` (e.g., `[0.0, 0.0, 1.0, 0.0]`)
- **Normals:** Always have `w=0` and are unit-length (magnitude = 1)

This enables uniform handling of affine transformations via 4x4 matrix multiplication.

### 1.3 Assembly Structure (Existing)

From `geom3d.py` line 179-193, yapCAD already defines an **assembly** structure:

```python
assembly = ['assembly', transform, elementlist]

where:
  transform = [xformF, xformR]  # forward and reverse matrices
  elementlist = [element0, element1, ...]  # each is solid or assembly
```

**Our Implementation Must:**
- Extend this structure to add **mate constraints** and **metadata**
- Maintain backward compatibility with the existing format
- Follow the same tagged-list pattern

---

## 2. Naming Conventions

### 2.1 Function Names

- **snake_case** for all functions: `rotate_solid()`, `translate_surface()`, `ensure_solid_id()`
- **Verb-first naming:**
  - Creation: `make_*()`, `create_*()`
  - Conversion: `poly2surface()`, `surf2lines()`
  - Transformation: `rotatesurface()`, `translatesurface()`, `mirrorsurface()`
  - Predicates: `issurface()`, `issolid()`, `istriangle()`
  - Queries: `surfacebbox()`, `solidbbox()`, `volumeof()`

**Pattern for assembly module:**
```python
# Creation
def make_mate(mate_type, part1, part2, constraints):
def create_assembly(parts, mates=None):

# Transformation
def apply_mates(assembly):
def resolve_constraints(assembly):

# Predicates
def ismate(obj):
def isassembly(obj):

# Queries
def get_mate_status(mate):
def assembly_bounds(assembly):
```

### 2.2 Variable Names

- **Descriptive, full words:** `vertices`, `normals`, `faces` (not `verts`, `norms`, `fcs`)
- **Short forms allowed for:**
  - Common geometric entities: `p` (point), `v` (vector), `n` (normal), `tri` (triangle)
  - Iteration indices: `i`, `j`, `idx`
  - Matrix elements: `mat`, `trsf`, `xform`
  - Temporary/intermediate: `s2` (modified surface), `tmp`, `result`

### 2.3 Module Organization

- **Core geometry functions:** `assembly_core.py`
- **Mate constraint logic:** `mates.py`
- **Transformation helpers:** `assembly_xform.py`
- **Utilities:** `assembly_util.py`
- **Package exports:** `__init__.py` with explicit `__all__`

---

## 3. Docstring Style

### 3.1 Format

yapCAD uses **reStructuredText (reST)** style with emphasis on **geometric intent**:

```python
def rotate(x, ang, cent=point(0,0), axis=point(0,0,1.0), mat=False):
    """ return a rotated version of the surface, solid, or figure"""
    # Implementation...
```

**Key Characteristics:**
- **Brief summary line** describing what the function returns/does
- **No explicit Args/Returns sections** for simple functions
- **Inline parameter documentation** for complex cases
- **Geometric context** over implementation details

### 3.2 Complex Function Documentation

For functions with multiple modes or complex behavior (from `geom3d.py` line 373-401):

```python
def triTriIntersect(t1, t2, inside=True, inPlane=False, basis=None):
    """Function to compute the intersection of two triangles.  Returns
    ``False`` if no intersection, a line (a list of two points) if the
    planes do not overlap and there is a linear intersection, and a
    polygon (list of three or more points) if the triangles are
    co-planar and overlap.

    If ``inside == True`` (default) return line-segment or poly
    intersection that falls inside both bounded triangles, otherwise
    return a line segment that lies on the infinite linear
    intersection of two planes, or False if planes are degenerate.

    If ``inPlane==True``, return the intersection as a poly in the
    planar coordinate system implied by ``t1``, or in the planar
    coordinate system specified by ``basis``

    [More parameter details...]
    """
```

**Pattern:**
1. One-line summary
2. Behavior explanation (what it returns in different cases)
3. Parameter descriptions inline with explanation
4. Edge cases and special modes

### 3.3 Assembly Module Docstring Guidelines

```python
def apply_mate(assembly, mate_id):
    """Apply a mate constraint to an assembly, updating part transforms.

    Returns a modified assembly with the mate constraint satisfied. If the
    constraint cannot be satisfied, raises a ValueError with details about
    the failure condition.

    The ``mate_id`` must reference a mate defined in the assembly's mate list.
    """

def create_planar_mate(part1, face1, part2, face2, offset=0.0):
    """Create a planar mate constraint between two part faces.

    A planar mate aligns two planar faces such that their normals are
    anti-parallel and the surfaces are separated by the specified offset
    distance. Positive offset moves part2 away from part1 along the normal.
    """
```

---

## 4. Error Handling

### 4.1 Exception Strategy

yapCAD uses **`ValueError`** for invalid inputs, with **descriptive messages**:

```python
# From geom3d.py line 217
if m < epsilon:
    raise ValueError('degenerate face in tri2p0n')

# From xform.py line 264
if inverse:
    sx = 1.0/sx  # Could raise ZeroDivisionError implicitly
else:
    raise ValueError('bad scaling values passed to Scale')
```

**Pattern:**
- Validate inputs at function entry
- Raise `ValueError` with context (function name, bad value)
- Include the problematic value in error messages when helpful
- Use `f-strings` for formatting error messages

### 4.2 No Custom Exceptions

yapCAD uses **built-in exceptions only:**
- `ValueError` - invalid parameters
- `RuntimeError` - operation failures (BREP operations, etc.)
- `ImportError` - optional dependencies missing

**Assembly module should follow this:**
```python
def apply_mate(assembly, mate_id):
    if not isassembly(assembly):
        raise ValueError(f'invalid assembly passed to apply_mate')

    mate = find_mate(assembly, mate_id)
    if not mate:
        raise ValueError(f'mate {mate_id} not found in assembly')

    try:
        result = solve_mate_constraints(assembly, mate)
    except Exception as e:
        raise RuntimeError(f'failed to apply mate {mate_id}: {e}')

    return result
```

### 4.3 Graceful Optional Dependencies

From `brep.py` (line 62-74):

```python
try:
    from OCC.Core.TopoDS import TopoDS_Shape
    _HAVE_OCC = True
except ImportError as exc:
    TopoDS_Shape = Any
    _HAVE_OCC = False

def require_occ():
    if _HAVE_OCC:
        return
    raise RuntimeError("pythonocc-core is not available...")
```

**Pattern for assembly dependencies (if any):**
```python
try:
    import numpy as np
    _HAVE_NUMPY = True
except ImportError:
    _HAVE_NUMPY = False

def require_numpy():
    if _HAVE_NUMPY:
        return
    raise RuntimeError("numpy required for mate solving (optional dependency)")
```

---

## 5. Type Hints

### 5.1 Current Usage

yapCAD has **mixed type hint usage:**

- **`metadata.py`:** Full type hints with imports from `typing`
  ```python
  from typing import Dict, Tuple, Any, Optional, Iterable

  def get_surface_metadata(surface: list, create: bool = False) -> Dict:
  ```

- **`geom.py`, `geom3d.py`:** No type hints (pure Python 3.10 without annotations)

- **`brep.py`:** Type hints for external-facing APIs and class methods

### 5.2 Recommended Approach for Assembly Module

**Use type hints for:**
- Public API functions
- Functions accepting/returning complex structures
- Functions with non-obvious return types

**Omit type hints for:**
- Internal helpers
- Functions with obvious signatures
- Performance-critical geometry operations

**Example:**
```python
from typing import List, Dict, Any, Optional, Tuple

# Public API - use hints
def create_assembly(
    parts: List[list],
    mates: Optional[List[dict]] = None,
    metadata: Optional[Dict[str, Any]] = None
) -> list:
    """Create an assembly from parts and mate constraints."""
    # Implementation

# Internal helper - no hints needed
def _resolve_mate_transform(mate, part1, part2):
    """Compute transformation to satisfy mate constraint."""
    # Implementation
```

### 5.3 Geometry List Type Aliases

Consider defining type aliases for clarity:

```python
from typing import List, Tuple

Point = List[float]  # [x, y, z, w]
Vector = List[float]  # [x, y, z, w]
Matrix = List[List[float]]  # 4x4 transformation
Surface = list  # ['surface', ...]
Solid = list  # ['solid', ...]
Assembly = list  # ['assembly', ...]
```

---

## 6. Separation of Concerns

### 6.1 Geometry vs Rendering

**Core Principle:** Geometry computation is **completely independent** of rendering.

From `CLAUDE.md`:
> **Separation of geometry and rendering**: Pure geometry operations are independent of rendering systems. Geometry computation can occur without invoking any renderer.

**Architecture:**
```
geom.py, geom3d.py           <- Pure geometry (no rendering dependencies)
    ↓
drawable.py                  <- Abstract rendering interface
    ↓
ezdxf_drawable.py           <- DXF export
pyglet_drawable.py          <- OpenGL visualization
```

### 6.2 Assembly Module Separation

**Assembly module should follow the same pattern:**

```
assembly_core.py            <- Mate constraints, transforms (pure geometry)
    ↓
assembly_util.py            <- Assembly manipulation utilities
    ↓
assembly_drawable.py        <- Rendering exploded views, mate visualizations
                              (optional, separate module)
```

**Rules:**
- **Core assembly functions:** No imports from `drawable`, `pyglet`, `ezdxf`
- **Visualization:** Separate module that imports core + drawable
- **STL export:** Should work without any rendering system

---

## 7. Functional vs Object-Oriented Programming

### 7.1 Dominant Style: Functional

yapCAD is **primarily functional** with **minimal OOP:**

- **Geometry primitives:** Functions operating on lists (not classes)
- **Transformations:** Functions returning new geometry (immutability)
- **Type checking:** Predicate functions (`issurface()`, `issolid()`)

**Only classes used:**
- `Matrix` (xform.py) - wraps 4x4 transformation matrices
- `Polygon`, `Geometry`, `Circle`, `RoundRect` (poly.py, geometry.py) - high-level construction
- `BrepSolid`, `BrepVertex`, `BrepEdge`, `BrepFace` (brep.py) - BREP wrappers
- `Drawable` subclasses - rendering backends

### 7.2 Assembly Module Approach

**Prefer functions over classes:**

```python
# PREFERRED - Functional style
def create_planar_mate(part1, face1, part2, face2, offset=0.0):
    """Create a planar mate constraint (returns dict)."""
    return {
        'type': 'planar',
        'part1': part1,
        'face1': face1,
        'part2': part2,
        'face2': face2,
        'offset': offset
    }

def apply_mate(assembly, mate):
    """Apply mate constraint, returning modified assembly."""
    # Functional transformation
    return updated_assembly

# OPTIONAL - Class wrapper for convenience (like Polygon, Geometry)
class Assembly:
    """High-level assembly builder (wrapper around functional core)."""
    def __init__(self, parts=None):
        self._assembly = create_assembly(parts or [])

    def add_mate(self, mate):
        self._assembly = apply_mate(self._assembly, mate)
        return self
```

**When to use classes:**
- High-level builder APIs (like `Polygon`, `RoundRect`)
- Wrappers around external libraries (like `BrepSolid`)
- Rendering backends (like `pygletDraw`)

---

## 8. Metadata Handling

### 8.1 Metadata Structure

From `metadata.py` (line 12-21):

```python
_SURFACE_META_INDEX = 6
_SOLID_META_INDEX = 4

# Surface: ['surface', verts, normals, faces, boundary, holes, metadata?]
#                                                                  ↑ index 6
# Solid: ['solid', surfaces, material, construction, metadata?]
#                                                      ↑ index 4
```

**Metadata is a dictionary appended as the last element:**

```python
{
    "schema": "metadata-namespace-v1.1",
    "entityId": "uuid-string",
    "tags": ["tag1", "tag2"],
    "layer": "default",
    "material": {...},
    "manufacturing": {...},
    # ... namespaced sections
}
```

### 8.2 Assembly Metadata Pattern

**Extend the existing assembly structure:**

```python
# Original: ['assembly', transform, elementlist]
# Extended: ['assembly', transform, elementlist, metadata]

_ASSEMBLY_META_INDEX = 3

def get_assembly_metadata(assembly, create=False):
    """Get metadata dict from assembly, creating if needed."""
    if not isassembly(assembly):
        raise ValueError('invalid assembly passed to get_assembly_metadata')

    if len(assembly) <= _ASSEMBLY_META_INDEX:
        if not create:
            return {}
        meta = _ensure_root({})
        assembly.append(meta)
        return meta

    return assembly[_ASSEMBLY_META_INDEX]

def set_assembly_mates(assembly, mates):
    """Store mate constraints in assembly metadata."""
    meta = get_assembly_metadata(assembly, create=True)
    meta['mates'] = mates
    return assembly
```

### 8.3 UUID Management

yapCAD uses UUIDs for entity tracking:

```python
import uuid

def ensure_assembly_id(assembly):
    """Ensure assembly has a unique entityId."""
    meta = get_assembly_metadata(assembly, create=True)
    if 'entityId' not in meta:
        meta['entityId'] = str(uuid.uuid4())
    if 'id' not in meta:
        meta['id'] = meta['entityId']
    return meta['id']
```

---

## 9. Testing Patterns

### 9.1 Test Organization

From `CLAUDE.md`:
- Non-visual tests: `pytest tests/ -m "not visual"`
- Visual tests: Mark with `@pytest.mark.visual`
- Long tests: Mark with `@pytest.mark.slow`

### 9.2 Test Structure

```python
# tests/test_assembly.py

import pytest
from yapcad.geom import point
from yapcad.geom3d import solid
from yapcad.assembly import create_assembly, create_planar_mate

class TestAssemblyCore:
    """Core assembly creation and manipulation."""

    def test_create_empty_assembly(self):
        asm = create_assembly([])
        assert isassembly(asm)
        assert len(get_assembly_parts(asm)) == 0

    def test_create_assembly_with_parts(self):
        part1 = make_test_solid()
        part2 = make_test_solid()
        asm = create_assembly([part1, part2])
        assert len(get_assembly_parts(asm)) == 2

class TestMateConstraints:
    """Mate constraint creation and solving."""

    def test_planar_mate_creation(self):
        mate = create_planar_mate(
            part1='part1_id', face1=0,
            part2='part2_id', face2=0,
            offset=5.0
        )
        assert ismate(mate)
        assert mate['type'] == 'planar'

    @pytest.mark.slow
    def test_complex_mate_solving(self):
        """Test multi-constraint assembly solving (may take >1 sec)."""
        # Complex assembly with multiple mates
        pass

@pytest.mark.visual
class TestAssemblyVisualization:
    """Visual tests requiring display."""

    def test_render_exploded_view(self):
        """Render exploded assembly view (requires closing window)."""
        # Interactive test
        pass
```

---

## 10. Code Quality Standards

### 10.1 Linting

From `CLAUDE.md`:
> Run `flake8` before review (config in `setup.cfg`)

**Assembly module should:**
- Pass `flake8` with existing project config
- Follow PEP 8 with 4-space indentation
- Keep line length ≤ 100 characters (check `setup.cfg` for limit)

### 10.2 Comments

yapCAD has **minimal inline comments**, preferring:
- **Descriptive variable names**
- **Clear function names**
- **Docstrings over inline comments**

**When comments are used:**
```python
# From geom3d.py line 228
n[3] = 0.0 #direction vectors lie in the w=0 hyperplane

# From xform.py line 223
# see http://www.opengl-tutorial.org/assets/faq_quaternions/index.html#Q38
## This is correct
R = [[cang + ux*ux*cmin, ...], ...]
```

**Pattern:**
- Explain **why**, not **what**
- Reference algorithms/papers/standards
- Mark corrected implementations

---

## 11. Performance Considerations

### 11.1 Efficient Data Structures

From `geom3d.py` (line 608):
```python
# Use chain.from_iterable instead of reduce for O(n) performance
if not filterInds(chain.from_iterable(faces), verts):
    return False
```

**Guidelines:**
- Prefer built-in functions over manual loops for large datasets
- Use `itertools` for efficient iteration
- Avoid repeated deep copies when unnecessary

### 11.2 Lazy Evaluation

From `brep.py` (line 141-172):
```python
# BREP cache to avoid repeated deserialization
_BREP_SOLID_CACHE: dict[str, "BrepSolid"] = {}

def brep_from_solid(solid: list) -> Optional["BrepSolid"]:
    """Return the cached BrepSolid for ``solid`` if metadata is present."""
    solid_id = meta.get("entityId")
    cached = _BREP_SOLID_CACHE.get(solid_id)
    if cached:
        return cached
    # ... deserialize only if not cached
```

**Assembly module pattern:**
```python
# Cache solved mate configurations
_ASSEMBLY_SOLUTION_CACHE: dict[str, dict] = {}

def solve_assembly_mates(assembly):
    """Solve mate constraints, using cache if available."""
    asm_id = ensure_assembly_id(assembly)
    if asm_id in _ASSEMBLY_SOLUTION_CACHE:
        return _ASSEMBLY_SOLUTION_CACHE[asm_id]

    solution = _compute_mate_solution(assembly)
    _ASSEMBLY_SOLUTION_CACHE[asm_id] = solution
    return solution
```

---

## 12. Integration with Existing Code

### 12.1 Extending Existing Structures

**DO:** Extend existing structures compatibly
```python
# Existing: ['assembly', transform, elementlist]
# Extended: ['assembly', transform, elementlist, metadata]

def isassembly(obj, fast=True):
    """Check if obj is a valid assembly (backward compatible)."""
    if not isinstance(obj, list) or len(obj) < 3:
        return False
    if obj[0] != 'assembly':
        return False
    if fast:
        return True
    # Full validation...
    return True
```

**DON'T:** Break existing assembly structures

### 12.2 Transformation Compatibility

**DO:** Use existing transformation matrices
```python
from yapcad.xform import Matrix, Translation, Rotation

def apply_mate_transform(part, mate):
    """Apply transformation using yapCAD Matrix."""
    T = Translation(mate['offset'])
    R = Rotation(mate['axis'], mate['angle'])
    mat = T.mul(R)
    return transform(part, mat)
```

**DON'T:** Create custom transformation systems

### 12.3 Metadata Integration

**DO:** Use existing metadata helpers
```python
from yapcad.metadata import (
    get_solid_metadata,
    ensure_solid_id,
    add_tags,
    set_layer
)

def tag_assembly_parts(assembly, tags):
    """Add tags to all parts in assembly."""
    for part in get_assembly_parts(assembly):
        meta = get_solid_metadata(part, create=True)
        add_tags(meta, tags)
    return assembly
```

---

## 13. File Headers and Licensing

All assembly module files should include:

```python
## yapCAD Assembly System - [module description]
## Copyright (c) 2026 yapCAD contributors
## All rights reserved

# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation files
# (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
# BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""
[Module docstring explaining purpose and key functions]
"""
```

---

## 14. Key Design Principles Summary

### 14.1 The yapCAD Way

1. **Functional first** - Pure functions operating on immutable data structures
2. **Lists over classes** - Tagged tuples, not complex object hierarchies
3. **Homogeneous coordinates** - 4-tuples enable uniform transformation handling
4. **Separation of concerns** - Geometry ≠ rendering ≠ serialization
5. **Graceful degradation** - Optional dependencies don't break core functionality
6. **Descriptive names** - Code should read like geometric descriptions
7. **Minimal abstraction** - Direct manipulation of geometric primitives

### 14.2 Assembly Module Checklist

- [ ] Use tagged list structures: `['mate', ...]`, `['assembly', ...]`
- [ ] Functions return new geometry (immutability by default)
- [ ] Type checking via predicate functions (`ismate()`, `isassembly()`)
- [ ] Metadata stored in dictionaries with UUIDs
- [ ] Transformations via `yapcad.xform.Matrix`
- [ ] Error messages include context and bad values
- [ ] Docstrings focus on geometric behavior
- [ ] No rendering dependencies in core module
- [ ] Tests separated into visual/non-visual
- [ ] Passes `flake8` with project config

---

## 15. Example Implementation Skeleton

```python
## yapCAD Assembly System - Core assembly and mate constraint functions
## Copyright (c) 2026 yapCAD contributors
[... license ...]

"""
Assembly and mate constraint system for yapCAD.

Provides functions for creating assemblies with spatial relationships
(mates) between parts. Mates define geometric constraints like planar
alignment, concentric axes, and fixed distances.

Key Functions:
    create_assembly()    - Create assembly from parts
    create_*_mate()      - Create mate constraints
    apply_mates()        - Solve and apply mate constraints
    isassembly()         - Check if object is valid assembly
"""

from yapcad.geom import point, vect
from yapcad.geom3d import issolid, solid
from yapcad.xform import Matrix, Translation, Rotation
from yapcad.metadata import ensure_solid_id
import uuid
from typing import List, Dict, Any, Optional

# Assembly metadata index
_ASSEMBLY_META_INDEX = 3

# Mate constraint cache
_MATE_SOLUTION_CACHE: Dict[str, Dict] = {}


def create_assembly(parts: List[list],
                   mates: Optional[List[dict]] = None,
                   transform: Optional[tuple] = None) -> list:
    """Create an assembly from a list of parts and optional mate constraints.

    Returns an assembly structure compatible with existing yapCAD format,
    extended with mate constraint support in metadata.

    If ``transform`` is not provided, uses identity transformation.
    If ``mates`` is provided, stores mate constraints in assembly metadata.
    """
    if transform is None:
        identity = Matrix()
        transform = (identity, identity)

    assembly = ['assembly', transform, list(parts)]

    if mates:
        meta = {
            'entityId': str(uuid.uuid4()),
            'mates': mates
        }
        assembly.append(meta)

    return assembly


def isassembly(obj, fast=True):
    """Check if obj is a valid assembly structure."""
    if not isinstance(obj, list) or len(obj) < 3:
        return False
    if obj[0] != 'assembly':
        return False
    if fast:
        return True
    # Full validation
    if not isinstance(obj[1], (tuple, list)) or len(obj[1]) != 2:
        return False
    if not isinstance(obj[2], list):
        return False
    return True


def create_planar_mate(part1_id: str, face1_idx: int,
                      part2_id: str, face2_idx: int,
                      offset: float = 0.0) -> dict:
    """Create a planar mate constraint between two part faces.

    A planar mate aligns two planar faces such that their normals are
    anti-parallel (pointing toward each other) and the surfaces are
    separated by the specified offset distance.

    Positive offset moves part2 away from part1 along face1's normal.
    """
    return {
        'type': 'planar',
        'id': str(uuid.uuid4()),
        'part1': part1_id,
        'face1': face1_idx,
        'part2': part2_id,
        'face2': face2_idx,
        'offset': offset
    }


# ... more functions following the same pattern ...
```

---

## 16. References

**Source files analyzed:**
- `/Users/jmika/Work/DML/yapCAD/src/yapcad/geom.py` (4526 lines)
- `/Users/jmika/Work/DML/yapCAD/src/yapcad/geom3d.py` (1456 lines)
- `/Users/jmika/Work/DML/yapCAD/src/yapcad/xform.py` (277 lines)
- `/Users/jmika/Work/DML/yapCAD/src/yapcad/metadata.py` (349 lines)
- `/Users/jmika/Work/DML/yapCAD/src/yapcad/brep.py` (740 lines)

**Key documentation:**
- `/Users/jmika/Work/DML/yapCAD/CLAUDE.md` - Project overview and conventions

**Contact:**
For questions about yapCAD conventions or assembly system design, refer to
project maintainers or the yapCAD documentation at docs/index.rst.
