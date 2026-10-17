# Feature-size responsive meshing and simplifier performance

Meshing a field and then simplifying it is currently the expensive half of
yapCAD's SDF tier, and for parts whose features are much smaller than the part
it is also the incorrect half. This note records what the cost actually is,
measured rather than assumed, and what the two fixes look like: one small and
immediate, one substantial.

The measurements come from engraved lettering, which is the worst case the
library currently meets: a 70 x 30 x 4 mm plate with three lines of 3 mm block
text cut 0.5 mm deep. The strokes are 0.43 mm wide, so feature and part differ
by a factor of about 160.

## What was measured

Composing that plate as one field -- a slab minus the union of 212 stroke
boxes -- and meshing once:

| resolution | cell | mesh | simplify | raw triangles | final | boundary edges |
| --- | --- | --- | --- | --- | --- | --- |
| 326 | 0.215 mm | 15.8 s | **657.6 s** | 227,512 | 8,784 | 62 |
| 490 | 0.143 mm | -- | -- | -- | -- | refused: `NonManifoldMeshError` |

For comparison, the same plate through the mesh path takes 0.19 s and produces
3,172 triangles with 135 non-manifold edges.

Two separate problems are visible. Simplification takes eleven minutes, which
is a performance problem. And the result is not closed -- the 62 boundary
edges are already present in the raw mesh, before simplification, so dual
contouring itself is producing an open surface. Raising the resolution does
not help: at 490 the mesher refuses the result because the surface passes
through single cells twice. There is no resolution at which this plate meshes
cleanly, which is the subject of the second half of this note.

## Where the time goes

Profiling simplification of a single engraved line (12,732 triangles, 83
rounds, 23.8 s):

| | time | share |
| --- | --- | --- |
| field evaluation (`evaluate`) | 20.8 s | **87%** |
| `_max_error`, which calls it | 20.7 s | -- |
| `_box_distance`, the largest leaf | 7.6 s own | 32% |
| `numpy.ufunc.reduce`, 618,892 calls | 4.5 s own | 19% |
| quadrics, QEF solves, edge linking | ~1.5 s | 6% |

Simplification is not slow because collapsing edges is expensive. It is slow
because it evaluates the field many times, and each evaluation costs far more
than the geometry requires.

### The union fan-out

`yapcad/sdf/ops.py`:

```python
def _eval_union(_params, p, children, ev):
    out = ev(children[0], p)
    for child in children[1:]:
        out = np.minimum(out, ev(child, p))
    return out
```

Every child is evaluated at every point, unconditionally. The plate's union
has 212 children, so every field sample performs 212 box-distance evaluations
-- while any given point is near at most one or two strokes. That single
factor accounts for most of the eleven minutes: roughly two orders of
magnitude more arithmetic than the shape requires. It is also why a 48-child
union produced 153,840 calls to `_eval_box` in the profile.

The same fan-out is paid by every consumer of a field, not just the
simplifier: `to_solid`'s grid sampling, clearance queries, and anything else
that evaluates a CSG tree built by union.

## Fix one: cull union children by bounds

`analyze(node)` already reports a per-node axis-aligned bound and a Lipschitz
constant, and the bounds propagate correctly through `translate` and `union`
(verified). Nodes built from the analytic primitives are 1-Lipschitz, so for a
batch of points with bounding box `B`, a child whose bounds are farther from
`B` than the running minimum distance cannot change the result of the `min`
and can be skipped outright.

This is the standard bounding-volume argument for SDF unions, and it is exact
-- it changes no value, it only declines to compute children that provably
cannot win. The expected gain on the plate is one to two orders of magnitude
on the dominant cost.

yapCAD already has the index this wants. `yapcad.octtree.NTree` is a tested
quad/octree over geometry elements keyed by bounding box, with
`getElements(bbox)` returning the elements that overlap a query box, and
`mindim`/`maxdepth` controls. It is currently used by the native boolean
engine, `native_brep` and the collision module. A union node can build an
`NTree` over its children's bounds once, then answer each evaluation batch
with a box query instead of a linear scan.

For small unions the scan is cheaper than the query, so the index should only
be built above a threshold -- a few dozen children -- and the linear path kept
below it.

### The rest of the ladder

In order of payoff per unit of risk, after culling:

1. **Hoist `_max_error`'s sampling.** Eighty-four calls account for nearly all
   evaluation time. Sample sparsely to reject a collapse early and densely
   only for collapses that look acceptable. The per-triangle cache already
   spans rounds, so the remaining cost is newly created triangles.
2. **Batch per round rather than per collapse.** 618,892 reduce calls for
   153,840 box evaluations is about four reductions each on small arrays:
   dispatch overhead, not arithmetic. Evaluating a whole round's candidate
   points as one flat array would remove most of it.
3. **Fewer rounds.** Eighty-three rounds on a 12.7k-triangle mesh, each
   re-scoring every edge. A priority queue with lazy revalidation, as in
   classical QEM, converges in far fewer passes. This changes the order
   collapses are applied and so changes results in detail, which the
   determinism note in `simplify.py` already anticipates.

## Fix two: feature-size responsive dual contouring

A uniform grid must be fine enough everywhere for the smallest feature
anywhere. The plate needs 0.14 mm cells across 70 mm -- about seven million
cells, some 95% of them in plain slab -- and still fails, because a 0.5 mm
pocket at a 0.215 mm cell is 2.3 cells deep and the surface begins passing
through single cells twice. Refining adaptively addresses the cost and the
correctness together.

### Structure

Replace the dense lattice with an octree, subdividing only where the field
needs it. Three refinement predicates, cheapest first:

**Emptiness pruning.** For a cell of half-diagonal `r` about `c`, if
`|f(c)| > L*r` the cell cannot contain the surface, and the whole subtree is
discarded. `analyze()` supplies `L`. On the plate this removes the slab
interior and the air above it almost immediately, which is where nearly all of
the seven million cells went.

**Feature-size refinement.** Subdivide while the cell is large relative to the
local feature size. A workable proxy is the variation of the gradient across
the cell's corners: where the normal swings sharply within one cell there is a
feature smaller than the cell. Refining while `cell > k*|f|` near the surface
has the same effect, descending into thin walls and pockets and stopping on
flat faces.

**Topology safety.** A cell must be split when its single dual vertex would
have to serve more than one sheet of surface. Count sign changes along the
cell's twelve edges and test the connectivity of the resulting configuration:
more than one connected component of surface inside a cell means subdivide.
This is exactly the condition a uniform grid cannot answer, and it is what
produces the boundary edges and coincident vertices seen above.

Vertex placement stays QEF-based, as now; the quadric machinery in
`simplify.py` is the same mathematics and should be shared rather than
duplicated.

### Stitching

Generating a crack-free surface across cells at different depths is the real
difficulty. The established answer is Ju, Losasso, Schaefer and Warren, *Dual
Contouring of Hermite Data* (2002), which traverses octree edges, faces and
cells recursively and handles neighbours at differing depths natively, so no
T-junctions arise. For a manifold guarantee, Schaefer, Ju and Warren,
*Manifold Dual Contouring* (2007) adds the cell-splitting rule above; it is
the published solution to precisely the failure recorded here.

`NTree` provides subdivision, traversal and bbox queries, and is tested, so it
is a reasonable foundation -- but it is an index over geometry elements, not a
field-sampled octree with per-cell dual vertices and sign configurations. The
honest expectation is that it supplies the tree mechanics while the cell
payload, refinement predicates and Hermite traversal are new.

### What it buys

- The open-mesh failure is fixed by construction: a cell carrying two sheets
  is split rather than collapsed.
- "Resolution too coarse for the feature" stops being a user-facing trap,
  because refinement is driven by the field. A tolerance in model units, as
  `sdf_cell_mm` already is, replaces a count against the bounding box.
- **Simplification becomes largely unnecessary.** This is the strategic
  payoff. Today the library meshes uniformly fine and then spends eleven
  minutes undoing it; adaptive meshing produces something near the simplified
  triangle count directly -- the plate's 8,784 rather than its 227,512. The
  simplifier becomes a polish pass on an already reasonable mesh and leaves
  the critical path.

### Cost and risk

This is a larger piece of work than the SDF simplifier was, and it is the kind
of code where subtle errors produce cracks that appear only on particular
geometry. It wants a conformance fixture from the start: fields with known
exact volumes and deliberately hard topology -- thin walls, corner pinches,
nested shells, pockets a few cells deep -- asserting watertightness, manifold
edges and volume at several tolerances. The conformance-fixture pattern in
`SDF-DESIGN.md` is the right precedent.

## On GPU acceleration

Field evaluation is embarrassingly parallel and would map well to a GPU, but
it is the wrong next step. The present problem is that roughly two hundred
times more work is done than necessary, not that the work is done too slowly;
accelerating the unculled union would spend significant effort making the
wrong algorithm fast, and would add a hard dependency to a library whose
pure-Python tier is a stated feature.

The case to revisit is after culling and adaptive meshing, and the target then
is `to_solid`'s grid sampling rather than the simplifier: dense, uniform and
branch-free, which is what GPUs want. Any such work should be optional,
behind the same kind of backend switch `to_solid` already has for `numpy`, and
justified by measurement on large assemblies.

## Suggested order

1. Union culling. Small, local, no API change, verifiable against existing
   results, and on the critical path for every consumer of a field.
2. Measure again. With evaluation cheap, uniform meshing plus simplification
   may be adequate for text, which would make adaptive meshing an
   optimisation rather than a prerequisite.
3. The simplifier ladder, as far as the measurements justify.
4. Adaptive dual contouring as its own project, with a design note and a
   conformance fixture, justified by large assemblies rather than by
   lettering alone.
5. GPU work last, if at all, with numbers in hand.

Related: #93 (`to_solid` fails opaquely when a field undersamples to nothing),
#94 (dual contouring emits open meshes on thin features), #82 (engraved block
text is not watertight), #85 (degeneracy judged against an absolute epsilon).
