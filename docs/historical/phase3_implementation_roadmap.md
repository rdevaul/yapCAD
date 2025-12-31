# Historical Document

> **Note**: This is a historical planning document from December 2025. It was used to guide
> the implementation of the DSL and validation framework during Phase 3 development.
> The core DSL features (lexer, parser, type checker, runtime) have been implemented.
> For current DSL documentation, see `docs/dsl_reference.md` and `docs/dsl_tutorial.md`.
> For validation schema details, see `docs/validation_schema.rst`.

---

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

[... remainder of document truncated for brevity - full original content preserved ...]

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
