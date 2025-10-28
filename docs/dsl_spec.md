# yapCAD DSL Specification (Draft)

**Version:** `yapdsl-v0.1`  
**Status:** Proposal – intent for next roadmap milestone

This document outlines an imperative, statically-checkable DSL for authoring reusable yapCAD modules (e.g., involute gears, fastener libraries). The DSL captures command definitions, parameter validation, canonical geometry generation, and metadata emission while allowing controlled fallbacks to Python for advanced logic. See `examples/involute_gear_package/src/involute_gear.dsl` for the draft module referenced below.

---

## 1. Design Goals

- **Deterministic & Auditable**: Commands produce canonical geometry entities with recorded parameters, package versions, and provenance metadata.
- **Composable**: Packages can consume other packages' DSL exports (e.g., a fastener package calling an involute gear package).
- **Static Validation**: Type hints, `require` constraints, and strongly typed geometry primitives catch mistakes before execution.
- **Python Escape Hatch**: Optional `python {}` blocks allow full expressiveness when needed but are marked as “untrusted” unless signed.
- **Instance Reuse**: DSL commands may register canonical entities once and refer to them via `manifest.instances` to avoid duplication (key for BOMs).

---

## 2. Module Structure

```
module involute_gear;

use math;
use fasteners.metric as metric;

command INVOLUTE_SPUR(
    teeth: int,
    module_mm: float,
    face_width_mm: float,
    pressure_angle_deg: float = 20.0
) -> solid {
    require teeth >= 6;
    require module_mm > 0;

    let profile: polygon2d = involute_profile(teeth, module_mm, pressure_angle_deg);
    let blank: solid = extrude(profile, face_width_mm);

    emit blank with {
        layer: "gear",
        derived2d: profile,
        metadata: {
            "application": "INVOLUTE_SPUR"
        }
    };
}

command INVOLUTE_SPUR2D(...) -> polygon2d {
    ...
}
```

### 2.1 Module Declarations
- `module <name>;` defines the namespace (used when other packages reference commands).
- `use` statements import other modules or standard libraries (`math`, `fasteners.metric`, etc.).

### 2.2 Command Signature
- `command NAME(param: type = default, ...) -> return_type { ... }`
- Supported primitive types: `int`, `float`, `string`, `bool`, `polygon2d`, `solid`, `sketch`, `transform`, `list<T>`, `dict`.
- DSL should allow type aliases and user-defined structs in later versions.

### 2.3 Statements
- `let` binds immutable values.
- `require` enforces preconditions; failures raise validation errors during assembly.
- `emit <geometry> with { ... }` returns geometry plus extra metadata.
  * `layer`: overrides default metadata layer.
  * `derived2d`: optional 2D geometry reference (stored under `geometry/derived/`).
  * `metadata`: free-form JSON object merged into entity metadata.
  * `register`: optionally register canonical entities (see §3.4).

---

## 3. Semantics & Runtime

### 3.1 Parameter Hashing & Invocation Metadata
- Each command invocation records:
  ```
  metadata.invocation = {
      package: "involute_gear",
      command: "INVOLUTE_SPUR",
      version: "0.2.0",
      parameters: {...},
      sourceSignature: "sha256:…"
  }
  ```
- Serialised geometry must include this block so downstream assemblies know which DSL command produced it.

### 3.2 Canonical Entities & Instance Registration
- DSL commands can call `register_instance(name: string, entity: solid|sketch, key: dict)` to write canonical geometry under `geometry/entities/`.
- `manifest.instances` references canonical entities via `entity: geometry/entities/<uuid>.json`, plus transforms or counts.
- Assembly commands refer to existing canonical entities via `use_instance(name, transform)` to avoid duplication (key for fastener arrays, gears, etc.).

### 3.3 Fallback to Python

```
command CUSTOM_PROFILE(config: dict) -> polygon2d {
    python {
        from scripts.custom_helpers import make_profile
        profile = make_profile(config)
    }
    emit profile with { layer: "custom" };
}
```

- Python blocks execute inside a controlled runtime (package `scripts/`).
- Metadata should capture `invocation["python"]["module"]` and hash of helper script.
- Unsigned Python blocks may require manual approval before packaging.

### 3.4 Metadata Helpers
- DSL should allow simple metadata tags:
  ```
  emit blank with {
      layer: "structure",
      tags: ["gear", "m2"],
      derived2d: profile
  };
  ```
- Layers default to `"default"` if unspecified; DSL authors can set more meaningful names.

---

## 4. Grammar Sketch (v0.1)

```
module_decl    ::= "module" IDENT ";"
use_stmt       ::= "use" IDENT ( "." IDENT )* ";"
command_decl   ::= "command" IDENT "(" params? ")" "->" type block
params         ::= param ( "," param )*
param          ::= IDENT ":" type ( "=" expr )?
block          ::= "{" stmt* "}"
stmt           ::= let_stmt | require_stmt | emit_stmt | python_block | call_stmt | ...
let_stmt       ::= "let" IDENT ":" type "=" expr ";"
require_stmt   ::= "require" expr ";"
emit_stmt      ::= "emit" expr "with" emit_block ";"
emit_block     ::= "{" emit_field ( "," emit_field )* "}"
emit_field     ::= IDENT ":" expr
python_block   ::= "python" block
expr           ::= literals | identifiers | function calls | operations
```

Detailed grammar (with precedence and type-checking rules) will evolve as we implement the compiler.

---

## 5. Tooling Roadmap

1. **DSL Parser/Checker** (`yapcad dsl lint`):
   - Parse modules, enforce types/`require` statements, flag Python fallbacks.
2. **DSL Compiler** (`yapcad dsl compile`):
   - Evaluate commands, capture canonical geometry, store invocation metadata, register `source.modules`.
3. **Assembler Integration**:
   - Allow packages to declare dependencies on other `.ycpkg` modules and call their commands.
4. **BOM/Instance Engine**:
   - Read `manifest.instances` and canonical entities to produce quantity rollups.

---

## 6. Open Questions

- Signature scheme and policy for Python fallback approval.
- How to distribute/resolve `.dsl` dependencies (e.g., package registry vs. local path).
- Versioning semantics when multiple packages export same command name.

Feedback is welcome as we start implementing the compiler and packaging integration.
