# Mesh boolean performance

With the optional `manifold3d` backend, DSL `union`, `difference`,
`intersection`, and their `_all` forms use native batch booleans. Intermediate
results stay native across statements and helper commands; translate, rotate,
and scale preserve that representation. The public execution result is an
ordinary yapCAD mesh, ready for STL/STEP export and other mesh consumers.
Reading surfaces inside an operation materializes the mesh and retires its
native copy, so later edits cannot reuse stale geometry. There is no global
geometry cache or helper-command memoization.

Explicit engine overrides, SDF-authoritative operands, BREP operands, and
invalid-input fallback retain their previous dispatch behavior. These changes
do not require pythonocc-core.

Prefer building tool lists with comprehensions when possible:

```text
tools: list<solid> = [translate(cylinder(1.0, 5.0), 4.0 * i, 4.0 * j, -2.5) for j in range(8) for i in range(20)]
emit difference_all(body, tools)
```

Nested comprehensions support `if` filters.

For reproducible comparisons, run `tools/benchmark_dsl_booleans.py` with and
without `--eager`, and with `--command BATCH`. It reports time including final
mesh conversion, volume, bounds, and triangle count. Compare geometry as well
as time: triangulations may differ while describing the same solid.
