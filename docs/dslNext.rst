==================================
DSL Language Extensions: dslNext
==================================

This document describes planned extensions to the yapCAD DSL that reduce
verbosity while maintaining static verifiability guarantees.

Design Principles
=================

The yapCAD DSL was designed with the following safety properties:

1. **Termination Guarantee**: All well-formed DSL programs terminate
2. **Type Safety**: All type errors are caught before execution
3. **No Arbitrary Code Execution**: DSL cannot escape to arbitrary Python
4. **Deterministic Execution**: Same inputs always produce same outputs

These properties enable DSL code to be safely executed in untrusted contexts,
such as processing user-submitted designs in a web service.

The extensions in this document preserve all of these properties by adding
only **total functions** and **bounded iteration** constructs.

Phase 1: Conditional Expressions
================================

Status: **Implemented** (December 2025)

Motivation
----------

Currently, handling conditional logic requires either:

- Multiple command variants (verbose, duplicated logic)
- Runtime branching via helper commands (awkward)

Conditional expressions allow inline value selection without control flow.

Syntax
------

::

    # Ternary conditional expression
    value: type = expr_true if condition else expr_false

    # Examples
    module: float = size if metric else 25.4 / size
    color: string = "red" if temperature > 100.0 else "blue"

    # Nested (parentheses recommended for clarity)
    grade: string = "A" if score >= 90.0 else ("B" if score >= 80.0 else "C")

Grammar Extension
-----------------

::

    conditional_expr := or_expr ("if" or_expr "else" conditional_expr)?

The ``else`` branch is required (expressions must always produce a value).
Both branches must have compatible types.

Type Rules
----------

::

    Γ ⊢ e_cond : bool    Γ ⊢ e_true : T    Γ ⊢ e_false : T
    ─────────────────────────────────────────────────────────
              Γ ⊢ (e_true if e_cond else e_false) : T

Both branches must have the same type (or compatible types under subtyping).

Implementation Notes
--------------------

Lexer:
    No changes needed - ``if`` and ``else`` are already keywords.

Parser (``parser.py``):
    Extend ``_parse_expression()`` to handle trailing ``if ... else``.
    Precedence: conditional is the lowest precedence expression form.

Checker (``checker.py``):
    Add ``visit_ConditionalExpr`` that:

    1. Checks condition has type ``bool``
    2. Checks both branches
    3. Verifies branches have compatible types
    4. Returns the unified type

Interpreter (``interpreter.py``):
    Add ``visit_ConditionalExpr`` that:

    1. Evaluates condition
    2. Evaluates and returns appropriate branch (short-circuit)

AST Node
--------

::

    @dataclass
    class ConditionalExpr(Expr):
        condition: Expr
        true_branch: Expr
        false_branch: Expr

Examples
--------

Before (verbose, multiple commands)::

    command make_gear_metric(teeth: int, module: float) -> solid:
        emit involute_gear(teeth, module, 20.0, 10.0)

    command make_gear_imperial(teeth: int, dp: float) -> solid:
        module: float = 25.4 / dp
        emit involute_gear(teeth, module, 20.0, 10.0)

After (single command with conditional)::

    command MAKE_GEAR(teeth: int, size: float, metric: bool = true) -> solid:
        module: float = size if metric else 25.4 / size
        emit involute_gear(teeth, module, 20.0, 10.0)


Phase 2: List Comprehensions
============================

Status: **Implemented** (pre-existing)

The DSL already supports list comprehensions with map and filter syntax.

Motivation
----------

Creating lists of transformed or filtered elements currently requires:

- Manual enumeration of elements
- Multiple statements with append operations
- Helper commands for generation

List comprehensions provide concise, declarative list construction.

Syntax
------

::

    # Map: transform each element
    [expr for var in iterable]

    # Filter: select elements
    [expr for var in iterable if condition]

    # Multiple variables (parallel iteration)
    [expr for var1, var2 in zip(list1, list2)]

    # Examples
    squares: list<float> = [x * x for x in values]
    angles: list<float> = [i * 45.0 for i in range(8)]
    positives: list<float> = [x for x in values if x > 0.0]
    scaled: list<float> = [x * scale for x in values if x > threshold]

Grammar Extension
-----------------

::

    list_comprehension := "[" expr "for" IDENTIFIER "in" expr ("if" expr)? "]"

The ``for ... in`` binds a loop variable; the optional ``if`` filters.

Type Rules
----------

::

    Γ ⊢ e_iter : list<S>    Γ, x : S ⊢ e_body : T    (no filter)
    ──────────────────────────────────────────────────────────────
              Γ ⊢ [e_body for x in e_iter] : list<T>

    Γ ⊢ e_iter : list<S>    Γ, x : S ⊢ e_cond : bool    Γ, x : S ⊢ e_body : T
    ───────────────────────────────────────────────────────────────────────────
              Γ ⊢ [e_body for x in e_iter if e_cond] : list<T>

Implementation Notes
--------------------

Lexer:
    No changes needed - ``for``, ``in``, ``if`` are already keywords.

Parser (``parser.py``):
    In ``_parse_primary()``, when encountering ``[``:

    1. Look ahead to distinguish list literal from comprehension
    2. If ``for`` follows first expression, parse as comprehension
    3. Otherwise, parse as list literal

Checker (``checker.py``):
    Add ``visit_ListComprehension`` that:

    1. Checks iterable expression has type ``list<S>``
    2. Opens new scope with loop variable typed as ``S``
    3. If filter present, checks it has type ``bool``
    4. Checks body expression, gets type ``T``
    5. Returns type ``list<T>``

Interpreter (``interpreter.py``):
    Add ``visit_ListComprehension`` that:

    1. Evaluates iterable to get list
    2. For each element:
       a. Bind loop variable
       b. If filter present, evaluate; skip if false
       c. Evaluate body, append to result
    3. Return result list

AST Node
--------

::

    @dataclass
    class ListComprehension(Expr):
        body: Expr              # Expression to evaluate for each element
        variable: str           # Loop variable name
        iterable: Expr          # Expression producing the list to iterate
        condition: Optional[Expr]  # Optional filter condition

Examples
--------

Generating angles for circular patterns::

    # Before
    angles: list<float> = [0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0]

    # After
    angles: list<float> = [i * 30.0 for i in range(12)]

Creating transformed geometry::

    # Create array of gears
    gears: list<solid> = [
        translate(involute_gear(24, 2.0, 20.0, 10.0), i * 50.0, 0.0, 0.0)
        for i in range(5)
    ]

Filtering values::

    # Select only positive values
    positives: list<float> = [x for x in measurements if x > 0.0]

    # Select even tooth counts
    valid_teeth: list<int> = [n for n in range(6, 50) if n % 2 == 0]


Phase 3: Functional Combinators
===============================

Status: **Partially Implemented** (December 2025)

Aggregation functions that operate on lists, enabling concise reduction of
collections without manual iteration.

Implemented Builtins
--------------------

**3D Boolean Aggregation:**

``union_all(list<solid>)``
    Union all solids in a list.

    ::

        combined: solid = union_all(gears)

``difference_all(base: solid, tools: list<solid>)``
    Subtract all solids in tools list from base.

    ::

        plate_with_holes: solid = difference_all(plate, holes)

``intersection_all(list<solid>)``
    Intersect all solids in a list.

    ::

        common: solid = intersection_all(volumes)

**2D Boolean Aggregation:**

``union2d_all(list<region2d>)``
    Union all regions in a list.

    ::

        combined: region2d = union2d_all(shapes)

``difference2d_all(base: region2d, tools: list<region2d>)``
    Subtract all regions in tools list from base.

    ::

        plate: region2d = difference2d_all(rectangle(100.0, 60.0), holes)

``intersection2d_all(list<region2d>)``
    Intersect all regions in a list.

    ::

        overlap: region2d = intersection2d_all(regions)

**Numeric Aggregation:**

``sum(list<float>)``
    Sum all values in a list.

    ::

        total: float = sum(values)

``product(list<float>)``
    Multiply all values in a list.

    ::

        factorial_5: float = product([1.0, 2.0, 3.0, 4.0, 5.0])

``min_of(list<float>)``
    Find minimum value in a list.

    ::

        smallest: float = min_of(measurements)

``max_of(list<float>)``
    Find maximum value in a list.

    ::

        largest: float = max_of(measurements)

**Boolean Aggregation:**

``any_true(list<bool>)``
    Return true if any element is true.

    ::

        has_error: bool = any_true([x < 0.0 for x in values])

``all_true(list<bool>)``
    Return true if all elements are true.

    ::

        all_valid: bool = all_true([x > 0.0 for x in values])

Future Builtins
---------------

``reduce(func, list, initial)``
    Fold a list using a binary function.

    ::

        # Union all solids
        combined: solid = reduce(union, parts, empty_solid())

        # Sum values
        total: float = reduce(add, values, 0.0)

``map(func, list)``
    Transform each element. (Mostly subsumed by comprehensions, but
    useful for passing named functions.)

    ::

        scaled: list<float> = map(double, values)

``filter(predicate, list)``
    Select elements matching predicate. (Mostly subsumed by comprehensions.)

``zip(list1, list2)``
    Combine two lists into list of pairs.

    ::

        pairs: list<tuple<float, float>> = zip(xs, ys)

``enumerate(list)``
    Pair each element with its index.

    ::

        indexed: list<tuple<int, solid>> = enumerate(parts)

Type Rules for Reduce
---------------------

::

    Γ ⊢ f : (T, T) -> T    Γ ⊢ lst : list<T>    Γ ⊢ init : T
    ──────────────────────────────────────────────────────────
                  Γ ⊢ reduce(f, lst, init) : T

Implementation requires first-class function types or special-casing
known combiners (``union``, ``intersection``, ``add``, etc.).


Phase 4: Anonymous Functions (Lambdas)
======================================

Status: Future (Optional)

If Phase 3 requires arbitrary function arguments, consider limited lambdas.

Syntax
------

::

    # Single-expression lambda
    lambda x: x * 2.0
    lambda x, y: x + y

Restrictions for Safety
-----------------------

1. **No recursion**: Lambdas cannot reference themselves
2. **No mutation**: Captured variables are read-only
3. **Pure expressions only**: Body must be a single expression
4. **No closures over mutable state**: Only capture immutable values

These restrictions ensure lambdas are total functions.


Verification Guarantees
=======================

All proposed extensions maintain the DSL's safety properties:

Termination
-----------

+------------------------+------------------------------------------+
| Feature                | Termination Argument                     |
+========================+==========================================+
| Conditional expression | Both branches are expressions (no loops) |
+------------------------+------------------------------------------+
| List comprehension     | Bounded by input list length             |
+------------------------+------------------------------------------+
| ``reduce``             | Bounded by input list length             |
+------------------------+------------------------------------------+
| ``map`` / ``filter``   | Bounded by input list length             |
+------------------------+------------------------------------------+
| Lambda                 | No recursion allowed                     |
+------------------------+------------------------------------------+

Type Safety
-----------

All features include precise typing rules that can be checked statically.
The checker verifies:

- Condition expressions have type ``bool``
- Comprehension iterables have type ``list<T>``
- Both conditional branches have compatible types
- ``reduce`` function signature matches list element type

No Escape to Arbitrary Code
---------------------------

- No ``eval`` or ``exec``
- No imports of Python modules
- Lambdas are pure expressions, not arbitrary callables
- All builtins are pre-defined, verified functions


Implementation Roadmap
======================

Phase 1: Conditional Expressions
--------------------------------

Files to modify:

1. ``src/yapcad/dsl/ast.py`` - Add ``ConditionalExpr`` node
2. ``src/yapcad/dsl/parser.py`` - Parse ternary expressions
3. ``src/yapcad/dsl/checker.py`` - Type-check conditionals
4. ``src/yapcad/dsl/runtime/interpreter.py`` - Evaluate conditionals
5. ``tests/test_dsl_parser.py`` - Parser tests
6. ``tests/test_dsl_checker.py`` - Type checking tests
7. ``tests/test_dsl_runtime.py`` - Runtime tests

Phase 2: List Comprehensions
----------------------------

Files to modify:

1. ``src/yapcad/dsl/ast.py`` - Add ``ListComprehension`` node
2. ``src/yapcad/dsl/parser.py`` - Parse comprehension syntax
3. ``src/yapcad/dsl/checker.py`` - Type-check comprehensions
4. ``src/yapcad/dsl/runtime/interpreter.py`` - Evaluate comprehensions
5. ``tests/test_dsl_parser.py`` - Parser tests
6. ``tests/test_dsl_checker.py`` - Type checking tests
7. ``tests/test_dsl_runtime.py`` - Runtime tests

Phase 3: Functional Combinators
-------------------------------

Files to modify:

1. ``src/yapcad/dsl/symbols.py`` - Register new builtins
2. ``src/yapcad/dsl/runtime/builtins.py`` - Implement builtins
3. ``src/yapcad/dsl/types.py`` - Add function types if needed
4. ``tests/test_dsl_runtime.py`` - Runtime tests

Phase 4: Anonymous Functions
----------------------------

Files to modify:

1. ``src/yapcad/dsl/ast.py`` - Add ``LambdaExpr`` node
2. ``src/yapcad/dsl/tokens.py`` - Add ``lambda`` keyword
3. ``src/yapcad/dsl/lexer.py`` - Recognize ``lambda``
4. ``src/yapcad/dsl/parser.py`` - Parse lambda expressions
5. ``src/yapcad/dsl/types.py`` - Add function type ``(T1, T2) -> R``
6. ``src/yapcad/dsl/checker.py`` - Type-check lambdas
7. ``src/yapcad/dsl/runtime/interpreter.py`` - Evaluate lambdas


Example: Complete Gear Array
============================

Demonstrating Phases 1-3 together::

    module gear_train

    command MAKE_GEAR_TRAIN(
        count: int,
        base_teeth: int = 24,
        module_mm: float = 2.0,
        spacing: float = 50.0,
        alternating: bool = true
    ) -> solid:
        require count > 0 and count <= 10
        require base_teeth >= 12

        # Phase 1: Conditional expression for alternating sizes
        # Phase 2: List comprehension for generation
        gears: list<solid> = [
            translate(
                involute_gear(
                    base_teeth + (4 if alternating and i % 2 == 1 else 0),
                    module_mm,
                    20.0,
                    10.0
                ),
                i * spacing,
                0.0,
                0.0
            )
            for i in range(count)
        ]

        # Phase 3: Functional combinator for assembly
        emit union_all(gears)


Appendix: Comparison with Other DSLs
====================================

OpenSCAD
--------

OpenSCAD uses ``for`` as both iteration and generation::

    for (i = [0:10]) translate([i*10, 0, 0]) cube(5);

This is similar to our comprehension approach but mixes statements and
expressions. Our approach keeps them separate for clearer semantics.

CadQuery
--------

CadQuery uses Python's full syntax, including comprehensions::

    result = cq.Workplane("XY").box(1,1,1)
    for i in range(5):
        result = result.union(cq.Workplane("XY").box(1,1,1).translate((i*2,0,0)))

This is expressive but not statically verifiable. Our DSL provides similar
expressiveness with safety guarantees.

JSCAD
-----

JSCAD uses JavaScript, which has comprehension-like ``map``/``filter``::

    const cubes = Array.from({length: 5}, (_, i) =>
        translate([i*10, 0, 0], cube({size: 5}))
    );
    return union(cubes);

Our syntax is cleaner while providing similar functionality.
