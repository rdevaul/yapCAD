# Phase 3 Implementation Roadmap: Parametric DSL & Validation Layer

**Version**: Draft 2.0
**Date**: December 1, 2025
**Status**: Complete - All Design Questions Resolved

## Executive Summary

Phase 3 delivers two major capabilities:
1. **Parametric DSL** - A statically-checkable domain-specific language for authoring reusable yapCAD geometry modules
2. **Validation Framework** - An execution manager for running analysis/simulation plans and capturing results

These components enable deterministic, auditable geometry generation with integrated design verification.

---

## Current State

### What Exists
- Package manifest schema with validation plan structure (`ycpkg_spec.rst`)
- DSL grammar sketch and design goals (`dsl_spec.rst`)
- Package creation, validation, and export CLI tools
- Placeholder analyzer CLI (`ycpkg_analyze.py`)
- Material schema support (just completed)

### What's Missing
- DSL lexer/parser
- DSL type checker
- DSL compiler/runtime
- Validation execution manager with backend adapters
- Canonical entity registration system
- Instance management and BOM generation

---

## Type System Architecture

The DSL type system is organized into five tiers that mirror yapCAD's geometric abstractions:

```
                    ┌─────────────────────────────────────────┐
                    │           TIER 5: SOLIDS                │
                    │              solid                      │
                    └───────────────────┬─────────────────────┘
                                        │
                    ┌───────────────────▼─────────────────────┐
                    │         TIER 4: SURFACES                │
                    │         surface, shell                  │
                    └───────────────────┬─────────────────────┘
                                        │
                    ┌───────────────────▼─────────────────────┐
                    │     TIER 3: COMPOUND CURVES             │
                    │  path2d, path3d, profile2d,             │
                    │  region2d, loop3d                       │
                    └───────────────────┬─────────────────────┘
                                        │
                    ┌───────────────────▼─────────────────────┐
                    │      TIER 2: CURVE PRIMITIVES           │
                    │  line_segment, arc, circle, ellipse,    │
                    │  parabola, hyperbola, catmullrom,       │
                    │  nurbs, bezier                          │
                    └───────────────────┬─────────────────────┘
                                        │
                    ┌───────────────────▼─────────────────────┐
                    │   TIER 1: NUMERIC & GEOMETRIC PRIMITIVES│
                    │  int, float, bool, string               │
                    │  point, point2d, point3d                │
                    │  vector, vector2d, vector3d             │
                    │  transform                              │
                    └─────────────────────────────────────────┘
```

**Key Type Relationships**:
- Tier 2 curves compose into Tier 3 paths
- `profile2d` auto-promotes to `region2d` when closed (|start - end| < ε)
- Tier 3 `region2d` can be extruded/revolved/swept to create Tier 5 solids
- Tier 4 surfaces are the constituents of Tier 5 solid boundaries (BREP)
- `point` and `vector` are dimensionally polymorphic (infer 2D/3D from context)

---

## Implementation Milestones

### Milestone 3.1: DSL Foundation (Lexer + Parser + Type Checker)

**Goal**: Parse and validate DSL source files, report errors before execution

#### 3.1.1 Lexer Implementation
```
src/yapcad/dsl/
├── __init__.py
├── lexer.py          # Token definitions, scanner
├── tokens.py         # Token types enum
└── errors.py         # DSL-specific exceptions
```

**Tokens to support**:
- Keywords: `module`, `use`, `command`, `let`, `require`, `emit`, `with`, `python`, `as`, `close`, `closeC0`, `closeC1`
- Types (Tier 1 - Primitives): `int`, `float`, `string`, `bool`, `point`, `point2d`, `point3d`, `vector`, `vector2d`, `vector3d`, `transform`
- Types (Tier 2 - Curves): `line_segment`, `arc`, `circle`, `ellipse`, `parabola`, `hyperbola`, `catmullrom`, `nurbs`, `bezier`
- Types (Tier 3 - Compound): `path2d`, `path3d`, `profile2d`, `region2d`, `loop3d`
- Types (Tier 4 - Surfaces): `surface`, `shell`
- Types (Tier 5 - Solids): `solid`
- Generic types: `list`, `dict`
- Operators: `=`, `+`, `-`, `*`, `/`, `<`, `>`, `<=`, `>=`, `==`, `!=`, `&&`, `||`, `.`
- Delimiters: `{`, `}`, `(`, `)`, `[`, `]`, `:`, `;`, `,`, `->`
- Identifiers, literals (int, float, string, bool)

**Deliverables**:
- [ ] `DslLexer` class with `tokenize(source: str) -> List[Token]`
- [ ] Comprehensive token position tracking (line, column)
- [ ] Clear error messages for invalid tokens
- [ ] Unit tests for lexer

#### 3.1.2 Parser Implementation
```
src/yapcad/dsl/
├── parser.py         # Recursive descent parser
├── ast.py            # AST node definitions
└── grammar.py        # Grammar documentation/helpers
```

**AST Nodes**:
```python
@dataclass
class Module:
    name: str
    uses: List[UseStatement]
    commands: List[Command]

@dataclass
class Command:
    name: str
    params: List[Parameter]
    return_type: Type
    body: List[Statement]

@dataclass
class Parameter:
    name: str
    type: Type
    default: Optional[Expression]

# Statement types: LetStatement, RequireStatement, EmitStatement, PythonBlock
# Expression types: Literal, Identifier, BinaryOp, FunctionCall, etc.
```

**Deliverables**:
- [ ] `DslParser` class with `parse(tokens: List[Token]) -> Module`
- [ ] Complete AST hierarchy for DSL grammar
- [ ] Syntax error recovery with meaningful messages
- [ ] Unit tests for parser

#### 3.1.3 Type Checker Implementation
```
src/yapcad/dsl/
├── types.py          # Type system definitions
├── checker.py        # Type checking visitor
└── symbols.py        # Symbol table management
```

**Type System - Five-Tier Hierarchy**:

The DSL type system is organized into five tiers that mirror yapCAD's geometric abstractions:

**Tier 1: Numeric & Geometric Primitives**
```
# Scalars
int, float, bool, string

# Geometric primitives (foundation for transforms)
point           # Dimensionally polymorphic - infers 2D/3D from context
point2d         # [x, y, 1] in homogeneous coordinates
point3d         # [x, y, z, 1]
vector          # Dimensionally polymorphic
vector2d        # [dx, dy, 0] - direction, no position
vector3d        # [dx, dy, dz, 0]
transform       # 4x4 matrix (3x3 for 2D-only contexts)
```

**Tier 2: Curve Primitives (Analytic/Parametric)**
Each curve type supports parametric evaluation at t ∈ [0,1]:
```
# Analytic curves
line_segment    # { start: point, end: point }
arc             # { center: point, radius: float, start_angle: float, end_angle: float }
circle          # { center: point, radius: float }
ellipse         # { center: point, major: float, minor: float, rotation: float }
parabola        # { vertex: point, focus: point }
hyperbola       # { center: point, a: float, b: float }

# Spline curves
catmullrom      # { control_points: list<point>, tension?: float }
nurbs           # { control_points: list<point>, weights: list<float>, knots: list<float>, degree: int }
bezier          # { control_points: list<point> }  # cubic bezier common case
```

**Tier 3: Compound Curves (Paths & Regions)**
```
# Open curves (edges, sweep spines)
path2d          # Sequence of Tier 2 curves in 2D plane, C0/C1 continuous
path3d          # Sequence of Tier 2 curves in 3D space

# Closed curves (profiles, boundaries)
profile2d       # Open path2d that may be promoted to region2d
region2d        # Closed path2d - forms a planar region (for extrusion, etc.)
loop3d          # Closed path3d - forms a surface boundary (BREP edge loop)
```

**Tier 4: Surfaces (BREP Support)**
```
surface {
    type: planar | cylindrical | conical | spherical | toroidal | nurbs_surface | ...
    boundary: loop3d | list<loop3d>  # outer + inner loops (holes)
}
shell           # Connected surfaces forming a boundary (open or closed)
```

**Tier 5: Solids**
```
solid {
    shells: list<shell>  # Outer shell + optional void shells
}
```

**Curve Evaluation Operations**:
Parametric curves support evaluation within the DSL:
```
let c: catmullrom = catmullrom([p1, p2, p3, p4]);
let midpoint: point3d = c.at(0.5);           # Evaluate position at t=0.5
let tangent: vector3d = c.tangent_at(0.5);   # Tangent vector at t=0.5
let curvature: float = c.curvature_at(0.5);  # Curvature at t=0.5
```

**Type Promotion Rules**:
- `profile2d` automatically promotes to `region2d` when `|eval(t=0) - eval(t=1)| < epsilon`
- `close(path)` - explicit closure by extension (adds line segment)
- `closeC0(path)` - closure without derivative matching
- `closeC1(path)` - closure with C1 continuity (arc or higher-order curve)

**Surface Construction**:
```
let face: surface = planar_surface(boundary_loop);
let cyl: surface = cylindrical_surface(axis, radius, height, boundary);
let lofted: surface = loft_surface([profile1, profile2, profile3]);
```

**Generic Types**:
- `list<T>` - homogeneous list of type T
- `dict` - key-value mapping (string keys)

**Checks to perform**:
- Parameter type compatibility (including dimensional inference for point/vector)
- Return type matching
- `require` expressions must be boolean
- `emit` target must match command return type
- Undefined identifier detection
- Curve continuity validation (C0/C1 at path segment boundaries)
- Surface closure validation (for solid construction)
- Python block flagging (mark as "untrusted")

**Deliverables**:
- [ ] `TypeChecker` class with `check(module: Module) -> List[Diagnostic]`
- [ ] Symbol table with scope management
- [ ] Dimensional polymorphism resolver for point/vector types
- [ ] Type promotion logic for profile2d → region2d
- [ ] Built-in function type signatures for all tiers
- [ ] Curve evaluation type signatures
- [ ] Warning vs error severity levels
- [ ] Unit tests for type checker

#### 3.1.4 CLI Integration
```bash
# Lint DSL file (parse + type check)
yapcad dsl lint src/involute_gear.dsl

# Output:
# src/involute_gear.dsl:12:5: error: undefined identifier 'profile'
# src/involute_gear.dsl:18:3: warning: python block requires manual approval
```

**Deliverables**:
- [ ] `yapcad dsl lint` command
- [ ] JSON output format for CI integration
- [ ] Exit codes (0=clean, 1=warnings, 2=errors)

---

### Milestone 3.2: DSL Runtime & Compiler

**Goal**: Execute DSL commands to generate geometry with full provenance metadata

#### 3.2.1 Built-in Function Library
```
src/yapcad/dsl/
├── builtins/
│   ├── __init__.py
│   ├── primitives.py   # point, vector, curves, regions
│   ├── curves.py       # parametric curves, splines
│   ├── surfaces.py     # surface construction
│   ├── operations.py   # extrude, revolve, sweep, loft, boolean
│   ├── transforms.py   # translate, rotate, scale, mirror
│   ├── path_ops.py     # path joining, closure operations
│   └── math.py         # trig, vectors, etc.
```

**Built-ins by Tier**:

*Tier 1 - Primitives*:
- `point(x, y)` → `point2d`, `point(x, y, z)` → `point3d`
- `vector(dx, dy)` → `vector2d`, `vector(dx, dy, dz)` → `vector3d`
- `transform_translate(v: vector)` → `transform`
- `transform_rotate(axis: vector, angle: float)` → `transform`
- `transform_scale(factors: vector)` → `transform`
- `transform_compose(t1: transform, t2: transform)` → `transform`

*Tier 2 - Curve Primitives*:
- `line(start: point, end: point)` → `line_segment`
- `arc(center: point, radius: float, start: float, end: float)` → `arc`
- `circle(center: point, radius: float)` → `circle`
- `ellipse(center: point, major: float, minor: float, rot: float)` → `ellipse`
- `bezier(control_points: list<point>)` → `bezier`
- `catmullrom(control_points: list<point>, tension?: float)` → `catmullrom`
- `nurbs(controls: list<point>, weights: list<float>, knots: list<float>, degree: int)` → `nurbs`

*Tier 3 - Path/Region Operations*:
- `path(segments: list<curve>)` → `path2d` or `path3d`
- `join(p1: path, p2: path)` → `path` (C0 join)
- `close(p: profile2d)` → `region2d` (line segment closure)
- `closeC0(p: profile2d)` → `region2d` (no derivative matching)
- `closeC1(p: profile2d)` → `region2d` (C1 continuous closure)
- `rectangle(width: float, height: float, center?: point2d)` → `region2d`
- `regular_polygon(n: int, radius: float, center?: point2d)` → `region2d`

*Tier 4 - Surface Operations*:
- `planar_surface(boundary: loop3d)` → `surface`
- `cylindrical_surface(axis: vector3d, radius: float, height: float)` → `surface`
- `conical_surface(axis: vector3d, radius1: float, radius2: float, height: float)` → `surface`
- `spherical_surface(center: point3d, radius: float)` → `surface`
- `loft_surface(profiles: list<path3d>)` → `surface`
- `sweep_surface(profile: path2d, spine: path3d)` → `surface`
- `shell(surfaces: list<surface>)` → `shell`

*Tier 5 - Solid Operations*:
- `extrude(profile: region2d, height: float, direction?: vector3d)` → `solid`
- `revolve(profile: region2d, axis: vector3d, angle: float)` → `solid`
- `sweep(profile: region2d, spine: path3d)` → `solid`
- `loft(profiles: list<region2d>)` → `solid`
- `box(width: float, depth: float, height: float)` → `solid`
- `cylinder(radius: float, height: float)` → `solid`
- `sphere(radius: float)` → `solid`
- `cone(radius1: float, radius2: float, height: float)` → `solid`

*Boolean Operations* (support both binary and variadic forms, method chaining):
```
# Binary form
union(a: solid, b: solid) -> solid
union(a: region2d, b: region2d) -> region2d

# Variadic form (list)
union(operands: list<solid>) -> solid
union(operands: list<region2d>) -> region2d

# Difference: args 2-n all subtracted from arg 1
difference(a: solid, b: solid) -> solid
difference(a: solid, tools: list<solid>) -> solid
difference(a: region2d, b: region2d) -> region2d
difference(a: region2d, tools: list<region2d>) -> region2d

# Intersection
intersection(a: solid, b: solid) -> solid
intersection(operands: list<solid>) -> solid
intersection(a: region2d, b: region2d) -> region2d
intersection(operands: list<region2d>) -> region2d

# Method chaining syntax (preferred)
# base.union(boss).difference(hole).difference([hole2, hole3])

# Canonical empty values
empty_solid() -> solid       # All empty solids equal this
empty_region() -> region2d   # All empty regions equal this
is_empty(s: solid) -> bool   # Equivalent to: s == empty_solid()
is_empty(r: region2d) -> bool

# Error handling:
# - Empty results are valid (canonical empty value)
# - Non-manifold results → runtime error with diagnostic
```

*2D Boolean Use Case - Tool-Path Generation*:
2D regions from boolean operations can be attached to extrusion metadata,
providing manufacturing hints for laser/waterjet tool-path generation.
The `extrude()` operation preserves the source `region2d` in provenance.

*Transform Operations* (forward-inverse pairing, left-to-right composition):

**Internal Representation**:
Every transform maintains both forward and inverse matrices analytically:
```
transform {
    forward: matrix4x4    # The transformation
    inverse: matrix4x4    # Analytically computed inverse
    kind: TransformKind   # IDENTITY|TRANSLATION|ROTATION|SCALE|RIGID|AFFINE|ARBITRARY
}
```

**Construction** (all create analytic forward/inverse pairs):
```
translate(v: vector) -> transform
rotate(axis: vector3d, angle: float) -> transform
rotate_2d(angle: float) -> transform
scale(factors: vector) -> transform
scale_uniform(factor: float) -> transform
mirror(plane_normal: vector3d) -> transform
mirror_2d(axis: vector2d) -> transform
identity_transform() -> transform

# Arbitrary matrix - requires numerical inversion, errors if singular
transform_from_matrix(m: list<list<float>>) -> transform
```

**Singular Transform Policy**:
Singular transforms (e.g., `scale(vector(1, 0, 1))`) are runtime errors at
construction, since every transform must have a valid inverse.

**Composition** (left-to-right = application order):
```
t1.compose(t2) -> transform
# Result: forward = t1.forward @ t2.forward
#         inverse = t2.inverse @ t1.inverse

# Chained example: translate, then rotate, then scale
let t: transform = translate(v).compose(rotate(axis, a)).compose(scale(s));
```

**Inverse** (instant - just swaps forward/inverse):
```
t.inverse() -> transform
```

**Application** (method chaining preferred):
```
# Method style
shape.apply(t) -> geometry
shape.translate(v).rotate(axis, angle) -> geometry  # Shorthand

# Functional style (also supported)
apply(t, shape) -> geometry
```

**Decomposition** (with `try_` prefix for safe optional access):
```
t.translation() -> vector3d           # Throws if no translation component
t.try_translation() -> vector3d?      # Returns none if not applicable

t.rotation_axis() -> vector3d         # Throws if no rotation component
t.try_rotation_axis() -> vector3d?
t.rotation_angle() -> float
t.try_rotation_angle() -> float?

t.scale_factors() -> vector3d
t.try_scale_factors() -> vector3d?

# Kind queries
t.is_rigid() -> bool      # Translation + Rotation only
t.is_affine() -> bool     # Translation + Rotation + Scale
t.kind() -> TransformKind
```

*Curve Evaluation Methods* (accessed via `.` notation):
- `curve.at(t: float)` → `point`
- `curve.tangent_at(t: float)` → `vector`
- `curve.normal_at(t: float)` → `vector`
- `curve.curvature_at(t: float)` → `float`
- `curve.length()` → `float`
- `curve.split_at(t: float)` → `(curve, curve)`

*Control Flow* (expressions, not statements - all return values):

**Conditionals**:
```
# If-else expression
let result: solid = if teeth > 20 {
    large_gear_profile(teeth)
} else {
    small_gear_profile(teeth)
};

# Match expression
let profile: region2d = match gear_type {
    "spur" => spur_profile(teeth),
    "helical" => helical_profile(teeth, helix_angle),
    _ => error("Unknown gear type")
};
```

**Loops and List Comprehensions**:
```
# For-in loop (returns unit, used for side effects like accumulation)
for i in 0..n {
    base = base.difference(holes[i]);
}

# List comprehension (returns list)
let holes: list<solid> = [cylinder(r, h).translate(p) for p in positions];
let result: solid = base.difference(holes);

# Range iteration
for i in 0..count { ... }
for item in list { ... }
```

**Pattern Operations** (preferred over explicit loops for common cases):
```
# Radial pattern - replicate around axis
radial_pattern(shape, count: int, axis: vector3d, center?: point) -> geometry

# Linear pattern - replicate along direction
linear_pattern(shape, count: int, spacing: vector) -> geometry

# Grid pattern - 2D array
grid_pattern(shape, rows: int, cols: int, row_spacing: vector, col_spacing: vector) -> geometry

# Along path - replicate along curve with optional alignment
pattern_along_path(shape, path: path3d, count: int, align?: bool) -> geometry
```

**Functional Operations**:
```
map(fn: (T) -> U, items: list<T>) -> list<U>
filter(pred: (T) -> bool, items: list<T>) -> list<T>
reduce(fn: (Acc, T) -> Acc, items: list<T>, initial: Acc) -> Acc
```

**Static Optimization**:
The compiler performs constant folding and static optimization:
- If-else with constant conditions: dead branch eliminated
- Pattern operations with constant counts: unrolled at compile time
- List comprehensions over constant ranges: evaluated at compile time
- Constant arithmetic expressions: folded to literals

*Query Operations*:

**Geometric Measurements**:
```
volume(s: solid) -> float
surface_area(s: solid) -> float
area(r: region2d) -> float
perimeter(r: region2d) -> float
bounding_box(g: geometry) -> BoundingBox
distance(a: point, b: point, tolerance?: float) -> float  # Uses EPSILON if omitted
centroid(s: solid) -> point3d
```

**Mass Properties** (pulls from material metadata by default):
```
mass(s: solid, density?: float) -> float
center_of_mass(s: solid, density?: float) -> point3d
moment_of_inertia(s: solid, density?: float, axis: vector3d) -> float
```

**Global Tolerance**:
```
EPSILON: float = 1e-9  # Referenceable constant
set_epsilon(value: float)  # Module-level override
```

**BREP Queries**:
```
num_faces(s: solid) -> int
faces(s: solid) -> list<surface>
edges(s: solid) -> list<path3d>
vertices(s: solid) -> list<point3d>

# Selector expressions for operations like fillet
edges(s).filter(e => e.length() > 10)  # Predicate selector
edges(s)[0, 3, 7]                       # Index selector
```

*Hole Annotations and Fastener System*:

**Hole Types**:
```
HoleAnnotation {
    id: uuid,                         # Stable identifier for tracking
    binding: "surface" | "volume",    # Surface-bound vs through-hole
    type: "clearance" | "tap" | "counterbore" | "countersink",
    fastener: FastenerRef?,           # Optional standard reference
    nominal_diameter: float,
    depth: float | "through",
    center: point3d,
    axis: vector3d,
    tapped: bool,
    fit: "close" | "normal" | "loose"?
}

FastenerRef {
    standard: "ISO" | "ANSI" | "DIN",
    designation: string,    # "M10", "#10-32", "1/4-20"
    length?: float,
    head_type?: string
}
```

**Annotated Hole Creation**:
```
clearance_hole(diameter: float) -> region2d
fastener_hole(standard: string, size: string, type: string, fit?: string) -> region2d
tap_hole(standard: string, size: string, tapped: bool = true) -> region2d
```

**Hole Queries**:
```
# From annotated geometry
holes(face: surface) -> list<HoleAnnotation>
fastener_holes(s: solid) -> list<HoleAnnotation>

# From non-annotated geometry (imports)
infer_holes(face: surface, tolerance?: float) -> list<InferredHole>
match_fasteners(holes: list<InferredHole>, standard: string) -> list<(InferredHole, FastenerRef?)>
```

**Fastener Tables**: Bundle ISO and ANSI standard tables; allow library imports for custom.

**Thread Representation**: Store annotation + generate geometry on demand. Provide post-processor for slicer thread data.

*Hole Inheritance in Boolean Operations*:

**Hole Classification**:
- Surface-bound: blind holes, counterbores, countersinks (attached to entry surface)
- Volume-bound: through-holes (property of solid volume, not surface)

**Union (A ∪ B)** - Inheritance Rules:
```
# Surface-bound holes on surviving surfaces: PRESERVED
# Through-holes: EXTENDED by default to traverse combined volume
# Options:
#   extend_through_holes: true (default) | false (convert to blind)
#   ignore_holes: false (default) | true (discard all annotations)
```

**Difference (A - B)** - Inheritance Rules:
```
# Note: Difference only removes material, cannot block through-holes
# Surface-bound holes on A's surviving surfaces: PRESERVED
# Through-holes in A: PRESERVED (may become longer if material removed)
# All holes from B: IGNORED (subtracted away)
```

**Intersection (A ∩ B)** - Inheritance Rules:
```
# Surface-bound holes: REMOVED (original surfaces don't survive)
# Through-holes passing through intersection volume in BOTH A and B: PRESERVED
# Through-holes in only one operand: REMOVED
# Option: ignore_holes: true to discard all annotations
```

**Hole Conflict Resolution** (for unions when holes overlap):
```
# Volume intersection of holes from different parts: ALWAYS flagged
# (potential fastener collision)

# Entry/exit surface overlap → possible "same hole" - four policies:
HoleConflictPolicy {
    "flag"      # Flag conflict, leave for manual resolution (default)
    "strict"    # Same hole if: diameter within epsilon AND
                # centers within epsilon on both entry and exit surfaces
    "promoting" # Like strict, but take larger diameter if they differ
    "loose"     # Any entry/exit overlap → same hole, take from A (leftmost)
}

# Usage:
union(a, b, hole_conflict: "strict", epsilon: 0.1)
union([a, b, c], hole_conflict: "loose")  # B's holes preferred over C's
```

**Partial Obstruction** (union creates partial blockage):
```
# Through-hole from A partially blocked by geometry from B
# Always flagged for review - may indicate design conflict
# Does not silently convert to blind hole
```

**Deliverables**:
- [ ] Built-in function registry with tier organization
- [ ] Type signatures for all built-ins (all 5 tiers)
- [ ] Mapping to yapCAD geometry functions
- [ ] Curve evaluation method dispatch
- [ ] Unit tests for each built-in

#### 3.2.2 Interpreter/Evaluator
```
src/yapcad/dsl/
├── interpreter.py    # AST evaluation
├── context.py        # Execution context (variables, scope)
└── provenance.py     # Invocation metadata capture
```

**Execution model**:
```python
class DslInterpreter:
    def execute_command(self, module: Module, command_name: str,
                        params: Dict[str, Any]) -> ExecutionResult:
        # 1. Validate parameters against command signature
        # 2. Create execution context
        # 3. Evaluate statements in order
        # 4. Capture require failures as validation errors
        # 5. Process emit statement, attach provenance
        # 6. Return geometry + metadata
```

**Provenance capture**:
```yaml
metadata:
  invocation:
    package: "involute_gear"
    command: "INVOLUTE_SPUR"
    version: "0.2.0"
    parameters:
      teeth: 24
      module_mm: 2.0
    sourceSignature: "sha256:..."
    timestamp: "2025-11-30T..."
```

**Deliverables**:
- [ ] `DslInterpreter` class
- [ ] Expression evaluator (arithmetic, comparisons, function calls)
- [ ] Statement evaluator (let, require, emit)
- [ ] Provenance metadata generation
- [ ] Python block executor with sandboxing considerations
- [ ] Integration tests with sample DSL modules

#### 3.2.3 Python Interop with Explicit Type Bridging

The DSL supports inline Python blocks with explicit type bridging for cases where yapCAD's full Python API is needed:

**Syntax**:
```
command CUSTOM_GEAR(teeth: int, module_mm: float) -> solid {
    # DSL variables automatically available in Python scope
    python {
        from yapcad.geom3d_util import involute_profile

        # DSL scalar types (int, float, bool, string) map directly to Python
        # DSL point/vector map to yapCAD's point() representation
        # DSL path2d/region2d map to yapCAD geometry lists
        # DSL solid maps to yapCAD solid tuples

        profile = involute_profile(teeth, module_mm)

        # Explicit return with type annotation - REQUIRED
        return profile as region2d
    }

    let gear: solid = extrude(profile, 10.0);
    emit gear;
}
```

**Type Bridging Rules**:

| DSL Type | Python/yapCAD Type | Direction |
|----------|-------------------|-----------|
| `int`, `float`, `bool`, `string` | Python native types | Bidirectional |
| `point`, `point2d`, `point3d` | yapCAD `point()` / `vect()` | Bidirectional |
| `vector`, `vector2d`, `vector3d` | yapCAD `vect()` | Bidirectional |
| `transform` | yapCAD 4x4 matrix | Bidirectional |
| `line_segment`, `arc`, `circle`, etc. | yapCAD geometry primitives | Bidirectional |
| `path2d`, `path3d`, `region2d` | yapCAD geometry list | Bidirectional |
| `surface` | yapCAD surface tuple | Bidirectional |
| `solid` | yapCAD solid tuple | Bidirectional |

**Return Statement Requirements**:
- Python blocks MUST use `return <value> as <dsl_type>` to return values
- The type checker validates the declared return type against usage
- Multiple returns are allowed (conditional logic) but all must have same type

**Security/Trust Model**:
- Python blocks are marked as "untrusted" during type checking
- Packages with Python blocks require explicit approval before execution
- The type system cannot verify Python code correctness - only interface types
- Import restrictions may be configured (see sandboxing options)

**Deliverables**:
- [ ] Python block parser (capture as string, validate syntax)
- [ ] Type bridge for DSL ↔ yapCAD conversions
- [ ] `return ... as <type>` statement parser and validator
- [ ] Python execution context with DSL variable injection
- [ ] Trust/approval mechanism for Python-containing packages

#### 3.2.4 Compiler CLI
```bash
# Compile DSL module, emit canonical geometry
yapcad dsl compile src/involute_gear.dsl \
    --command INVOLUTE_SPUR \
    --params teeth=24,module_mm=2.0,face_width_mm=10.0 \
    -o geometry/entities/gear_24t.json

# Compile and register in package
yapcad dsl compile src/involute_gear.dsl \
    --package my_design.ycpkg/ \
    --register
```

**Deliverables**:
- [ ] `yapcad dsl compile` command
- [ ] JSON/YAML parameter input support
- [ ] Output to standalone file or package registration
- [ ] Manifest `source.modules` population

---

### Milestone 3.3: Validation Execution Framework

**Goal**: Execute analysis plans, capture results, determine pass/fail status

#### 3.3.1 Validation Manager Core
```
src/yapcad/validation/
├── __init__.py
├── manager.py        # ValidationManager orchestrator
├── plan.py           # Plan loading and validation
├── result.py         # Result data structures
└── backends/
    ├── __init__.py
    ├── base.py       # Abstract backend interface
    ├── calculix.py   # CalculiX FEA adapter
    └── mock.py       # Mock backend for testing
```

**ValidationManager interface**:
```python
class ValidationManager:
    def load_plan(self, plan_path: Path) -> ValidationPlan:
        """Load and validate a plan YAML file."""

    def execute_plan(self, manifest: PackageManifest, plan: ValidationPlan,
                     *, dry_run: bool = False) -> ValidationResult:
        """Run validation, capture results, update manifest."""

    def get_backend(self, backend_name: str) -> ValidationBackend:
        """Get registered backend adapter."""
```

**Deliverables**:
- [ ] `ValidationManager` class
- [ ] Plan schema validation
- [ ] Result capture and serialization
- [ ] Manifest results section updates
- [ ] Unit tests for manager

#### 3.3.2 Backend Adapter Interface
```python
class ValidationBackend(ABC):
    @property
    @abstractmethod
    def name(self) -> str:
        """Backend identifier (e.g., 'calculix', 'comsol')."""

    @abstractmethod
    def prepare(self, plan: ValidationPlan, workspace: Path) -> None:
        """Generate solver input files."""

    @abstractmethod
    def execute(self, workspace: Path, options: Dict) -> ExecutionStatus:
        """Run the solver, return status."""

    @abstractmethod
    def collect_results(self, workspace: Path) -> Dict[str, Any]:
        """Parse solver output, extract metrics."""

    @abstractmethod
    def check_acceptance(self, results: Dict, criteria: Dict) -> Tuple[bool, List[str]]:
        """Evaluate pass/fail against acceptance criteria."""
```

**Deliverables**:
- [ ] `ValidationBackend` abstract base class
- [ ] Backend registration mechanism
- [ ] Mock backend for testing

#### 3.3.3 Validation Backend Strategy

**OCC Kernel Alignment**: yapCAD, Gmsh, and FEniCSx all use the OpenCASCADE kernel,
enabling direct geometry passing without lossy STEP/IGES conversions.

```
┌─────────────┐     OCC native     ┌─────────────┐     OCC native     ┌─────────────┐
│   yapCAD    │ ──────────────────▶│    Gmsh     │ ──────────────────▶│  FEniCSx    │
│   (BREP)    │                    │   (mesh)    │                    │   (solve)   │
└─────────────┘                    └─────────────┘                    └─────────────┘
       │                                                                     │
       │ STEP export                                               Results   │
       ▼                                                                     ▼
┌─────────────┐                                                    ┌─────────────┐
│ Other tools │ ◀──────────────────────────────────────────────────│  Post-proc  │
│ COMSOL/ANSYS│                                                    │  ParaView   │
└─────────────┘                                                    └─────────────┘
```

#### 3.3.4 Face Naming and Boundary Condition Assignment

**Face Naming at Creation**:
```
# DSL syntax for naming faces during creation
let bracket: solid = box(100, 50, 10) with {
    face_names: {
        top: "mounting_surface",
        bottom: "base",
        front: "load_face"
    }
};

# Or by selector after creation
bracket = bracket.name_faces({
    faces(bracket).filter(f => f.normal().z > 0.9): "top_faces",
    faces(bracket).filter(f => f.has_holes()): "bolt_faces"
});
```

**Face Naming for Imports**:
```
# Post-import annotation via selectors
let imported: solid = import_step("bracket.step");
imported = imported.name_faces({
    face_by_area_max(): "main_surface",
    faces_by_normal([0, 0, 1], tolerance: 5): "horizontal_faces"
});

# Or interactive naming in viewer (generates DSL code)
# User clicks face → enters name → code generated
```

#### 3.3.5 FEniCSx Backend Implementation

```python
class FenicsxBackend(ValidationBackend):
    """FEniCSx (DOLFINx) FEA backend - primary solver."""

    def prepare(self, plan, workspace):
        # 1. Pass OCC geometry directly to Gmsh via gmsh.model.occ
        # 2. Apply mesh hints (element size, refinement regions)
        # 3. Generate mesh with physical groups for BCs
        # 4. Convert to DOLFINx mesh format

    def execute(self, workspace, options):
        # 1. Set up function spaces (CG elements)
        # 2. Define variational problem from loads/constraints
        # 3. Apply boundary conditions by named face groups
        # 4. Solve linear/nonlinear system

    def collect_results(self, workspace):
        # 1. Extract stress tensor, displacement field
        # 2. Compute derived quantities (von Mises, principal)
        # 3. Find max values and locations
        # 4. Export VTU for visualization
```

#### 3.3.6 SU2 Backend Implementation (CFD)

```python
class SU2Backend(ValidationBackend):
    """SU2 CFD backend for aerospace applications."""

    def prepare(self, plan, workspace):
        # 1. Generate surface mesh via Gmsh
        # 2. Generate volume mesh (boundary layer support)
        # 3. Write SU2 configuration file
        # 4. Define inlet/outlet/wall boundary conditions

    def execute(self, workspace, options):
        # 1. Run SU2_CFD solver
        # 2. Monitor convergence (residuals)
        # 3. Handle restart for long simulations

    def collect_results(self, workspace):
        # 1. Parse surface solution file
        # 2. Extract Cp, Cf distributions
        # 3. Compute lift, drag coefficients
        # 4. Export flow field for visualization
```

#### 3.3.7 Gmsh Integration

```python
class GmshMesher:
    """Primary meshing interface using Gmsh's OCC integration."""

    def mesh_from_solid(self, solid: Solid, hints: MeshHints) -> Mesh:
        # 1. Import yapCAD solid via gmsh.model.occ.importShapes()
        # 2. Apply physical groups from face names
        # 3. Set mesh size fields from hints
        # 4. Generate mesh (2D surface or 3D volume)

    def apply_refinement(self, regions: List[FaceSelector], size: float):
        # Mesh refinement for stress concentrations, boundary layers

    def export_mesh(self, format: str) -> Path:
        # Export to: msh, vtk, xdmf (for FEniCSx), su2
```

**Deliverables**:
- [ ] `FenicsxBackend` class with structural analysis
- [ ] `SU2Backend` class with steady CFD
- [ ] `GmshMesher` class with OCC geometry import
- [ ] Face naming system (at creation and post-hoc)
- [ ] Physical group mapping from face names to mesh
- [ ] Integration tests with sample models

#### 3.3.4 Validation CLI
```bash
# Run a specific validation plan
yapcad package analyze my_design.ycpkg/ --plan bulkhead-fea

# Run all validation plans
yapcad package analyze my_design.ycpkg/ --all

# Dry-run (generate inputs without executing)
yapcad package analyze my_design.ycpkg/ --plan bulkhead-fea --dry-run
```

**Deliverables**:
- [ ] `yapcad package analyze` command
- [ ] Plan selection (single, multiple, all)
- [ ] Dry-run mode
- [ ] Progress reporting
- [ ] JSON output for CI integration

---

### Milestone 3.4: Canonical Entities & Instance Management

**Goal**: Enable geometry reuse through canonical entity registration and instancing

#### 3.4.1 Entity Registry
```
src/yapcad/dsl/
├── registry.py       # Canonical entity management
└── instances.py      # Instance/placement tracking
```

**Canonical entity storage**:
```
my_design.ycpkg/
└── geometry/
    └── entities/
        ├── gear_24t.json       # Canonical gear geometry
        └── bolt_m10.json       # Canonical fastener
```

**Registry operations**:
```python
class EntityRegistry:
    def register(self, name: str, entity: Any, params: Dict) -> str:
        """Register canonical entity, return UUID."""

    def lookup(self, name: str, params: Dict) -> Optional[str]:
        """Find existing entity matching params."""

    def get(self, entity_id: str) -> Any:
        """Retrieve entity by ID."""
```

**Deliverables**:
- [ ] `EntityRegistry` class
- [ ] Parameter-based deduplication
- [ ] Entity serialization to `geometry/entities/`
- [ ] Manifest `instances` section population

#### 3.4.2 Instance Placement
```python
class InstancePlacement:
    entity_id: str
    transform: Transform
    metadata: Dict[str, Any]

class InstanceManager:
    def create_instance(self, entity_id: str, transform: Transform) -> str:
        """Create new placement of canonical entity."""

    def get_instances(self, entity_id: str) -> List[InstancePlacement]:
        """Get all placements of an entity."""

    def generate_bom(self) -> Dict[str, int]:
        """Generate bill of materials from instance counts."""
```

**Deliverables**:
- [ ] `InstancePlacement` data class
- [ ] `InstanceManager` class
- [ ] Transform storage in `instances/` directory
- [ ] BOM generation utility

---

### Milestone 3.5: Integration & CLI Consolidation

**Goal**: Unified CLI experience and end-to-end workflow

#### 3.5.1 Unified CLI Entry Point
```bash
# DSL commands
yapcad dsl lint <file.dsl>
yapcad dsl compile <file.dsl> [options]

# Package commands (existing + new)
yapcad package create ...
yapcad package validate ...
yapcad package export ...
yapcad package analyze ...       # NEW
yapcad package assemble ...      # NEW - DSL-based assembly

# Combined workflow
yapcad package build src/design.dsl -o my_design.ycpkg/ \
    --export step,stl \
    --validate bulkhead-fea
```

#### 3.5.2 Assembly Command
```bash
# Assemble geometry from DSL invocations
yapcad package assemble my_design.ycpkg/ \
    --invoke "involute_gear:INVOLUTE_SPUR(teeth=24, module_mm=2.0)" \
    --invoke "fasteners:HEX_BOLT(size=M10, length=30)" \
    --instances gear_a:4 bolt_m10:100
```

**Deliverables**:
- [ ] `yapcad dsl` command group
- [ ] `yapcad package assemble` command
- [ ] `yapcad package build` convenience command
- [ ] Comprehensive `--help` documentation

---

## Testing Strategy

### Unit Tests
- Lexer token generation
- Parser AST construction
- Type checker diagnostics
- Interpreter expression evaluation
- Validation plan parsing
- Backend result collection

### Integration Tests
- DSL compile → geometry generation
- Validation plan execution (mock backend)
- Package creation with DSL modules
- End-to-end assembly workflow

### Regression Tests
- Existing 291+ tests must continue passing
- Add DSL-specific test suite
- Add validation framework test suite

---

## Dependencies

### Required
- Python 3.10+ (for pattern matching in parser)
- PyYAML (plan files)
- Existing yapCAD geometry modules

### Optional (for specific backends)
- CalculiX (`ccx` binary) - FEA
- gmsh - mesh generation
- VTK - result visualization

---

## Risk Assessment

| Risk | Impact | Likelihood | Mitigation |
|------|--------|------------|------------|
| DSL complexity creep | High | Medium | Start with MVP grammar, extend iteratively |
| Backend integration issues | Medium | High | Mock backends first, real backends as stretch |
| Performance with large models | Medium | Low | Lazy evaluation, incremental compilation |
| Python escape hatch abuse | Low | Medium | Clear documentation, signing requirements |

---

## Recommended Implementation Order

1. **Milestone 3.1** (Lexer/Parser/Type Checker) - Foundation
2. **Milestone 3.3.1-3.3.2** (Validation Manager + Backend Interface) - Independent track
3. **Milestone 3.2** (DSL Runtime) - Builds on 3.1
4. **Milestone 3.3.3-3.3.4** (CalculiX + CLI) - Builds on 3.3.1-3.3.2
5. **Milestone 3.4** (Canonical Entities) - Builds on 3.2
6. **Milestone 3.5** (Integration) - Ties everything together

Milestones 3.1 and 3.3.1-3.3.2 can proceed in parallel. This allows:
- DSL team: Lexer → Parser → Type Checker → Runtime
- Validation team: Manager → Backend Interface → CalculiX → CLI

---

## Comprehensive DSL Example

The following example demonstrates the full type system across multiple tiers:

```
module involute_gear_system;

use yapcad.stdlib.transforms;
use yapcad.stdlib.primitives;

# Command to generate an involute spur gear
command INVOLUTE_SPUR(
    teeth: int,
    module_mm: float,
    face_width: float,
    pressure_angle: float = 20.0
) -> solid {
    # Tier 1: Point and vector primitives
    let center: point2d = point(0.0, 0.0);
    let extrude_dir: vector3d = vector(0.0, 0.0, 1.0);

    # Parameter validation
    require teeth >= 8, "Minimum 8 teeth required";
    require module_mm > 0.0, "Module must be positive";
    require pressure_angle >= 14.5 && pressure_angle <= 25.0,
            "Pressure angle must be 14.5-25 degrees";

    # Calculated parameters
    let pitch_diameter: float = teeth * module_mm;
    let base_diameter: float = pitch_diameter * cos(radians(pressure_angle));
    let addendum: float = module_mm;
    let dedendum: float = 1.25 * module_mm;

    # Tier 2 & 3: Build tooth profile using parametric curves
    python {
        from yapcad.geom import involute_curve

        # Generate involute curve for tooth flank
        flank = involute_curve(base_diameter / 2, 0, 0.5)
        return flank as path2d
    }

    # Mirror for opposite flank, join into tooth profile
    let flank_mirrored: path2d = apply(mirror_y(), flank);
    let tooth_profile: profile2d = join(flank, flank_mirrored);

    # Close with tip arc and root fillet
    let tooth_region: region2d = closeC1(tooth_profile);

    # Replicate teeth around gear
    let gear_profile: region2d = radial_pattern(tooth_region, teeth, center);

    # Tier 5: Extrude to solid
    let gear: solid = extrude(gear_profile, face_width, extrude_dir);

    emit gear with {
        material: "steel_4140",
        description: "Involute spur gear"
    };
}

# Command to evaluate a point on a NURBS curve
command SAMPLE_CURVE(
    control_pts: list<point3d>,
    t_value: float
) -> point3d {
    require t_value >= 0.0 && t_value <= 1.0, "t must be in [0, 1]";

    # Tier 2: Create NURBS curve
    let curve: nurbs = nurbs(
        control_pts,
        weights: [1.0, 1.0, 1.0, 1.0],  # uniform weights
        knots: [0, 0, 0, 0.5, 1, 1, 1],
        degree: 3
    );

    # Curve evaluation - statically checkable
    let sampled: point3d = curve.at(t_value);
    let tangent: vector3d = curve.tangent_at(t_value);

    emit sampled;
}

# Command demonstrating surface construction
command LOFTED_DUCT(
    inlet_radius: float,
    outlet_radius: float,
    length: float,
    twist_angle: float = 0.0
) -> solid {
    # Tier 3: Create circular profiles
    let inlet: region2d = circle(point(0, 0), inlet_radius);
    let outlet_center: point3d = point(0, 0, length);

    # Transform outlet profile
    let t_rotate: transform = transform_rotate(vector(0, 0, 1), twist_angle);
    let t_translate: transform = transform_translate(vector(0, 0, length));
    let t_combined: transform = transform_compose(t_rotate, t_translate);

    let outlet: region2d = apply(t_combined, circle(point(0, 0), outlet_radius));

    # Tier 5: Loft between profiles
    let duct: solid = loft([inlet, outlet]);

    emit duct;
}
```

---

## Open Questions for Review

### Resolved Questions

1. **Type System Naming** ✓
   - Decision: Use `profile2d`/`region2d` instead of `polygon2d` to clarify that these can contain arcs, splines, and other curve types

2. **Dimensional Polymorphism** ✓
   - Decision: `point` and `vector` are dimensionally polymorphic (superclass of 2D/3D variants). Use explicit `point2d`/`point3d` when dimensionality is required.

3. **Curve Evaluation** ✓
   - Decision: Support parametric curve evaluation in DSL (`.at()`, `.tangent_at()`, etc.) to avoid Python fallback for common geometric operations

4. **Surface Construction** ✓
   - Decision: Support explicit surface construction (`planar_surface()`, `cylindrical_surface()`, etc.) for BREP workflows

5. **Type Promotion** ✓
   - Decision: `profile2d` auto-promotes to `region2d` when endpoints coincide within epsilon. Explicit closure ops: `close()`, `closeC0()`, `closeC1()`

6. **Python Interop Syntax** ✓
   - Decision: Use explicit type bridging with `return <value> as <dsl_type>` syntax

7. **Boolean Operations** ✓
   - Variadic support (binary and list forms)
   - 2D booleans on `region2d` for profile construction and tool-path generation
   - Method chaining syntax preferred
   - Canonical empty values (`empty_solid()`, `empty_region()`)
   - Non-manifold results → runtime error

8. **Transform Composition** ✓
   - Forward-inverse pairing (analytic inverses, no numerical computation)
   - Left-to-right composition order (application order)
   - Singular transforms → runtime error at construction
   - `try_` prefix for safe optional decomposition
   - TransformKind tracking for safe property extraction

9. **Control Flow** ✓
   - If-else and match as expressions (not statements)
   - For-in loops, list comprehensions
   - Pattern operations (radial, linear, grid, along-path)
   - Static optimization at compile time

10. **Query Operations** ✓
    - Mass properties pull from material metadata by default
    - Both index and predicate selectors for edges/faces
    - Optional tolerance parameter, global EPSILON constant

11. **Hole Annotations** ✓
    - Surface-bound (blind, counterbore, countersink) vs volume-bound (through-holes)
    - Bundle ISO/ANSI fastener tables, allow library imports
    - Thread representation: annotation + on-demand generation + slicer post-processor
    - Four hole conflict policies: flag, strict, promoting, loose

12. **Hole Inheritance** ✓
    - Union: extend through-holes by default, preserve surface-bound on surviving surfaces
    - Difference: preserve A's holes, ignore B's (difference only removes material)
    - Intersection: preserve through-holes in BOTH operands
    - Volume intersection of different holes → always flag for collision review

13. **Face Naming** ✓
    - Faces nameable at creation time via `with { face_names: {...} }` syntax
    - Post-hoc annotation via selectors for designed parts
    - Interactive naming in viewer for imports (generates DSL code)
    - Face names map to physical groups in Gmsh mesh

14. **Validation Behavior** ✓
    - Validation is non-blocking: failed validation flags results but doesn't block emit
    - Results captured in manifest with pass/fail status

15. **Meshing Strategy** ✓
    - Gmsh as primary mesher (OCC kernel alignment)
    - Direct OCC geometry passing (no STEP conversion losses)
    - STEP export pathway available for other tools

16. **Module/Import System** ✓
    - File extension: `.dsl`
    - Implicit import of `yapcad.stdlib` foundations
    - Circular dependencies: compile-time error (use interfaces to resolve)
    - Re-exports supported: `export use other.module;`
    - Module setup code: statically checkable DSL; Python requires explicit override

17. **Error Reporting** ✓
    - Max errors before stopping: 20 (override via env var or CLI)
    - Source maps: yes, Python blocks map back to DSL line numbers
    - I18n infrastructure from start, English-only initial implementation
    - JSON output format for tooling (LSP, CI)
    - Error codes: E0xx (lexer), E1xx (parser), E2xx (type), E3xx (semantic)

18. **User-defined Types** ✓
    - Defer to v0.2 - built-in five-tier type system + dict covers v0.1 use cases

19. **Python Block Sandboxing** ✓
    - Trust + flag approach: mark packages with Python as "requires-review"
    - No runtime sandboxing initially (complexity)
    - Module whitelist based on hash/signature or trusted signature list

20. **BOM Integration** ✓
    - Include in Phase 3 - tied to instance management milestone

21. **Incremental Compilation** ✓
    - Defer entirely - not needed for correctness, complexity not justified

22. **Curve Segment Storage** ✓
    - JSON with type discriminator format:
    ```json
    { "type": "arc", "center": [...], "radius": 10.0, ... }
    ```

### Validation Backend Priority (resolved above, item #3)
- **FEniCSx (DOLFINx)**: Primary FEA backend - OCC kernel alignment
- **Gmsh**: Primary meshing - OCC-based, direct geometry passing
- **SU2**: CFD backend for aerospace applications
- CalculiX: Secondary FEA option
- STEP export for other solvers

---

## All Questions Resolved

The Phase 3 roadmap specification is now complete with 22 resolved design decisions.

---

## Success Criteria

Phase 3 is complete when:

- [ ] DSL files can be parsed, type-checked, and compiled to yapCAD geometry
- [ ] Compiled geometry includes full provenance metadata
- [ ] Validation plans can be executed with at least one real backend
- [ ] Results are captured in manifest with pass/fail status
- [ ] CLI provides unified interface for DSL and validation workflows
- [ ] All new code has >80% test coverage
- [ ] Documentation updated with DSL language reference
