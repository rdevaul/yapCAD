# yapCAD DSL Language Guide

> **Where this fits:** This guide teaches the yapCAD DSL *as a language* — its
> mental model, execution semantics, and idioms. Read it **first**. When you
> want a worked, end-to-end example, go to the
> [DSL Tutorial](dsl_tutorial.md). When you need to look up a specific type or
> built-in function, go to the [DSL Reference](dsl_reference.md).
>
> **Learning path:** Guide (you are here) → Tutorial → Reference.

---

## 1. What the yapCAD DSL is (and why it exists)

The yapCAD DSL is a small, statically-typed, Python-flavored language for
describing **parametric CAD geometry as code**. A `.dsl` file is not a drawing
and not a saved model — it is a **program that generates geometry** when you
run it with a chosen set of parameters.

You could write the same geometry directly in Python against the `yapcad`
modules. The DSL exists because it deliberately gives up some of Python's power
in exchange for guarantees that matter for CAD automation and, increasingly,
for machine-generated designs:

- **Static type safety.** Every value has a declared type (`float`, `solid`,
  `region2d`, …). Type errors are caught by `dsl check` *before* any geometry
  is built, so a malformed design fails fast with a clear message instead of
  producing a broken solid.
- **Guaranteed termination.** There are no `while` loops and no unbounded
  recursion (call depth is capped). Every DSL program is statically verifiable
  to halt — which means a generator or LLM can emit DSL without the risk of an
  infinite loop wedging a build farm.
- **Provenance and reproducibility.** A DSL program plus its parameters fully
  determines the output. Package it as a `.ycpkg` and you have a signed,
  reproducible artifact with its source embedded.
- **A constrained surface for automation.** The DSL is small enough that tools
  (and language models) can author it reliably. Several shipped examples were
  generated from natural-language prompts.

If you have ever re-drafted a drawing because a material thickness, bolt
spacing, or pipe diameter changed, the DSL is the alternative: change a
parameter, re-run, done.

---

## 2. The execution model

A DSL program runs in one direction: **you invoke a `command` with parameters,
it computes geometry, and it `emit`s a single result.** Understanding this loop
is most of understanding the language.

```python
module my_part                       # 1. every file declares a module

command MAKE_BRACKET(                 # 2. a command is the unit of execution
    width: float = 40.0,             #    parameters are typed, may have defaults
    thickness: float = 5.0
) -> solid:                          # 3. the command declares its return type
    require width > 0.0              # 4. validate inputs (fail fast)
    plate: solid = box(width, width, thickness)   # 5. build geometry into typed vars
    emit plate                       # 6. emit exactly one result
```

Run it from the CLI:

```bash
python -m yapcad.dsl run my_part.dsl MAKE_BRACKET --param width=60.0 --output bracket.step
```

Key points about the model:

- **Commands are the entry points.** You don't run a "file"; you run a *command
  in* a file. One file can hold many commands.
- **`emit` is a return, not a print.** A command produces its output by emitting
  exactly one value of its declared return type. `emit` can also act as an early
  return inside a loop or conditional (see §5).
- **Execution is top-down and pure-ish.** There is no global mutable state
  shared between commands; a command's output depends only on its parameters and
  the helpers it calls.
- **`check` before `run`.** `dsl check` type-checks the whole file without
  building geometry. Make it a habit — it catches most mistakes in milliseconds.

---

## 3. Modules, commands, and the UPPERCASE/lowercase rule

Every file begins with `module <name>`. Inside it you define `command`s. The
**capitalization of a command name is semantically meaningful**:

| Name style | Role | Visible to `dsl list` / CLI? |
|---|---|---|
| `UPPERCASE_NAME` | **Exported** — a public entry point | ✅ Yes |
| `lowercase_name` | **Helper** — internal, composed by other commands | ❌ No |

This is the DSL's module-boundary mechanism: export the few commands a user
should call; keep the building-block helpers lowercase and private.

```python
module gearbox

# Helper — internal, not callable from the CLI
command make_tooth(module_mm: float) -> solid:
    emit box(module_mm, module_mm * 2.0, 5.0)

# Exported — appears in `dsl list`, callable from CLI
command MAKE_GEAR(teeth: int = 24, module_mm: float = 2.0) -> solid:
    body: solid = cylinder(teeth * module_mm / 2.0, 10.0)
    emit body
```

> **Idiom:** Design top-down. Write one `UPPERCASE` command that expresses the
> whole part, and factor repeated structure (a flange, a tooth, a rib) into
> `lowercase` helpers that it calls.

---

## 4. Thinking in types

The DSL is statically typed and every variable declares its type. This is the
single biggest adjustment for users coming from Python or OpenSCAD, and it is
the source of most of the DSL's safety.

```python
radius: float = 12.5                 # explicit annotation (preferred)
let height: float = 30.0             # 'let' form is equivalent
gear: solid = involute_gear(24, 2.0, 20.0, 10.0)
holes: list<solid> = []
```

A few rules that trip up newcomers:

- **Numbers are not interchangeable with their literal form.** Write `5.0`, not
  `5`, where a `float` is expected. Integer literals are `int`; the checker will
  tell you when a conversion is needed.
- **Both branches of a conditional expression must share a type.**
  `box(...) if flag else cylinder(...)` is fine (both `solid`); mixing a
  `solid` and a `float` is a type error.
- **There is no implicit negation operator surprise** — but see the gotcha in
  §7 about `0.0 - x`.

The type vocabulary (primitives, geometric types like `point3d`/`region2d`,
curve types, compound and generic types like `list<T>`) is cataloged in the
[DSL Reference → Types](dsl_reference.md). Treat the Reference as your
dictionary; this guide is the grammar.

### Angles: degrees for geometry, radians for trig

This is the single most common unit mistake, so internalize it early:

> **Geometry/transform functions take DEGREES. The trigonometric math
> functions take RADIANS.**

| Function | Angle unit |
|---|---|
| `rotate(solid, x, y, z)`, `rotate_2d(angle)`, `rotate_xform(...)`, `.rotate(...)` | **degrees** |
| `arc(...)` start/end angles, `ellipse(...)` start/end/rotation | **degrees** |
| `sin(x)`, `cos(x)`, `tan(x)` | **radians** |
| `asin`, `acos`, `atan`, `atan2` | **return radians** |
| `radians(deg)`, `degrees(rad)` | unit converters |

So rotating a part 45° is simply:

```python
rotated: solid = rotate(part, 0.0, 0.0, 45.0)   # 45 DEGREES — no conversion
```

But if you compute a position with trig, the trig call needs **radians** —
convert with `radians(...)`:

```python
# Place a hole on a bolt circle at `deg` degrees around the part:
deg: float = 30.0
x: float = bolt_circle_radius * cos(radians(deg))   # cos() needs radians
y: float = bolt_circle_radius * sin(radians(deg))
hole_at: solid = translate(hole, x, y, 0.0)
```

**The trap:** mixing the two — e.g. feeding a raw radian value into `rotate`,
or a raw degree value into `sin` — type-checks fine (both are `float`) but
produces silently wrong geometry. When in doubt, keep your angle variables in
degrees (matching the geometry functions) and wrap them in `radians(...)` only
at the moment you call a trig function.

---

## 5. Control flow — and the deliberate omissions

The DSL has the control flow you expect, **minus `while`**:

```python
# for over a range
for i in range(bolt_count):
    angle: float = i * (360.0 / bolt_count)
    # ... place a hole at `angle`

# for over a list
for hole in holes:
    body = difference(body, hole)

# if / elif / else
if depth <= 0:
    emit base_solid
elif depth < 3:
    emit small_variant
else:
    emit full_variant

# ternary (conditional expression)
unit: float = 1.0 if metric else 25.4
```

**Why no `while`?** Unbounded loops can't be statically proven to terminate. The
DSL trades them away so that *every* program is guaranteed to halt. When you
need "loop until converged," use a bounded `for` with an early `emit`:

```python
command SOLVE(start: float) -> float:
    x: float = start
    for i in range(100):             # hard upper bound
        if abs(error(x)) < 0.001:
            emit x                    # early return on convergence
        x = improve(x)
    emit x                            # best effort after max iterations
```

**Recursion is allowed but bounded.** Commands may call commands (including
themselves) up to a depth limit (default 100, configurable via
`--recursion-limit` or `YAPCAD_DSL_RECURSION_LIMIT`). This makes fractal/tree
structures expressible while preserving the termination guarantee.

---

## 6. Composing geometry: the core workflow

Most DSL programs follow the same shape: **build primitives, transform them,
combine them with booleans, emit the result.**

```python
command MAKE_WASHER(outer_d: float = 20.0, inner_d: float = 8.0, thick: float = 2.0) -> solid:
    require outer_d > inner_d
    disc: solid  = cylinder(outer_d / 2.0, thick)
    bore: solid  = cylinder(inner_d / 2.0, thick + 2.0)      # over-long to cut cleanly
    bore = translate(bore, 0.0, 0.0, -1.0)                   # straddle both faces
    emit difference(disc, bore)
```

Two patterns worth internalizing immediately:

**(a) Over-cut your subtractions.** When you `difference` a hole through a part,
make the cutting tool slightly *longer* than the part and offset it so it pokes
out both faces. A bore exactly as tall as the part can leave a zero-thickness
coplanar face that confuses the kernel. The `+ 2.0` / `translate(..., -1.0)`
idiom above is the standard fix.

**(b) Aggregate instead of chaining.** Combining many solids by nesting
`union(union(union(a, b), c), d)` is hard to read and easy to get wrong. Prefer
the aggregation helpers (see Reference → Aggregation) or a `for` loop:

```python
result: solid = parts[0]
for p in parts:
    result = union(result, p)
emit result
```

---

## 7. Idioms, conventions, and gotchas

These are the DSL-specific habits that separate "fighting the language" from
"fluent." Most are not obvious from the function catalog.

- **`require` is your contract.** Put `require` assertions at the top of every
  exported command to validate parameters. They produce clear errors and double
  as executable documentation of the valid parameter envelope.

- **`emit` exactly once on every path.** Every code path through a command must
  emit a value of the declared return type — including each branch of an
  `if/elif/else` and the fall-through after a loop.

- **Negation: write `0.0 - x`, not `-x` for computed floats.** Where you need
  the negative of a variable (e.g. mirroring a bolt position), the safe,
  type-clean form is `neg: float = 0.0 - x`. This shows up constantly in
  symmetric placement.

- **Degrees for geometry, radians for trig — don't cross the streams.**
  `rotate`, `rotate_2d`, and `arc`/`ellipse` angles are in **degrees**;
  `sin`/`cos`/`tan` take **radians**. Mixing them type-checks (both `float`)
  but produces silently wrong geometry. Keep angle variables in degrees and
  wrap them in `radians(...)` only at the trig call. (Full rundown in §4.)

- **Radial patterns via `for` + trig.** Bolt circles, fins, and spokes are a
  loop over `range(count)` computing an angle (in degrees), then
  `translate`/`rotate`. Remember to `radians(...)` the angle before any
  `sin`/`cos`. For common cases, prefer the built-in `radial_pattern` /
  `linear_pattern` helpers (Reference → Pattern helpers) over hand-rolled loops.

- **`print()` is for debugging only.** It writes to the console during `run`;
  it does not affect geometry. Use it to inspect intermediate values, then
  remove it.

- **Lowercase helpers, UPPERCASE entry points.** (See §3.) Keep the public
  surface small.

- **`use` imports are reserved for the future.** You may see `use other_module`
  in examples; cross-module imports are a forthcoming feature — for now keep a
  design within one module.

- **Method vs. function form.** Many operations are available both as functions
  and as methods (`translate(s, ...)` ≡ `s.translate(...)`). Pick one style per
  file for readability; the Reference → Method Syntax section lists what's
  available.

---

## 8. Annotating for tooling: `@meta` and `@ui`

Commands and parameters can carry **decorators** that are *informational only* —
the evaluator ignores them, so they never change geometry or type-checking:

- **`@ui(...)`** on a parameter gives the yapCAD Workbench widget hints
  (`widget="slider", min=..., max=...`).
- **`@meta(...)`** on a command attaches output metadata, optionally namespaced
  (`assembly.*`, `operation.*`) for assembly and machining workflows.

```python
@meta(operation.kind="subtract", operation.feature_kind="pocket")
command MAKE_POCKET(
    depth: float @ui(widget="slider", min=1.0, max=50.0) = 10.0
) -> solid:
    ...
```

You can write fully functional designs without ever using these; reach for them
when you're feeding the Workbench UI or the assembly/operation metadata system.
Full syntax and the v1.1 namespace vocabularies are in the
[DSL Reference → Parameter Decorators](dsl_reference.md).

---

## 9. From source to artifact: the CLI loop

The DSL ships with a CLI that covers the full lifecycle. You will use these
constantly:

```bash
python -m yapcad.dsl check  my_part.dsl                       # type-check, no build
python -m yapcad.dsl list   my_part.dsl                       # show exported commands
python -m yapcad.dsl run    my_part.dsl MAKE_PART --output part.step
python -m yapcad.dsl run    my_part.dsl MAKE_PART --param width=60.0 --param thick=4.0
python -m yapcad.dsl run    my_part.dsl MAKE_PART --package part.ycpkg
```

Typical inner loop while authoring: **edit → `check` → `run --output` →
inspect → repeat**, then `run --package` once you're happy to produce a signed,
reproducible `.ycpkg`. The [Tutorial](dsl_tutorial.md) walks this end-to-end on
a real part, including viewing and validating the package.

---

## 10. Where to go next

- **[DSL Tutorial](dsl_tutorial.md)** — build a parametric pipe fitting from
  scratch and export it. Best next step now that you have the mental model.
- **[DSL Reference](dsl_reference.md)** — the complete catalog: every type,
  built-in function, statement, method, and decorator. Your day-to-day lookup.
- **`examples/`** — dozens of real `.dsl` files (gears, fasteners, rockets,
  stands) showing idioms in context.
- **{doc}`Project Packaging <ycpkg_spec>`** — the `.ycpkg` format, signing, and
  provenance, for when you want reproducible, shareable design artifacts.
