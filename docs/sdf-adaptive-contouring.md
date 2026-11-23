# Adaptive dual contouring: implementation design

Status: stage A is implemented (`yapcad/sdf/octree.py`, opt-in through
`sampling="octree"`); stages B to D are proposals (#94). §3 builds the tree
bottom-up (§7, decision 6).

`sdf-meshing-and-simplifier-performance.md` ("Fix two") makes the case for
feature-size responsive meshing and names the literature. `SDF-DESIGN.md` §5.1
chose dual contouring on an adaptive octree from the start. This note is the
next step down: the data structure, the refinement tests, how the result
plugs into what `to_solid` already does, a staged plan in which each stage is
checkable on its own, and the decisions that need an owner before work
starts.

## 1. Where the uniform mesher stands

Some of what the earlier notes list as missing has since landed:

- **One vertex per surface component.** `dual_contour` already places a
  vertex per connected component of a cell's sign configuration (032da42),
  the vertex-splitting half of Manifold Dual Contouring. `SDF-DESIGN.md`
  §10.2 still says this is not implemented; that text predates the commit.
- **Refinement on failure.** `to_solid` halves the cell when it finds no
  surface or cannot close the mesh (dc11976), up to `MAX_RESOLUTION = 256`
  cells across.
- **Clean output.** Coincident vertices are merged by edge collapse, flat
  zero-thickness walls dropped and their seams zipped, quads split along the
  better diagonal, and zero-area triangles flipped away (#121, #122). Uniform
  meshes come out with no degenerate faces and no near-coincident vertices.
- **Speed.** Sampling is Lipschitz-pruned by chunk. The mesh backends are
  numpy, numba and wgpu (`auto` picks); numba and wgpu run the compiled tape
  of `yapcad/sdf/tape.py`. Measured on the engraved plate of
  `tests/test_mesh_degeneracy.py` (70 x 20 x 4 mm, two rows of block text,
  numba, resolution 490, so a 0.143 mm cell): the padded lattice is
  494 x 145 x 32 = 2,292,160 cells, 180,708 of them active (7.9%), and the
  mesh has 360,856 triangles, built in 2.2 s. #94's own plate is
  70 x 30 x 4 mm with text the issue does not give, so its figures (512,000
  triangles, 8,784 after simplification) are not reproduced here.

What remains is structural: every emitted cell is the same size.

1. **Thin features set the cell for the whole part.** A 0.43 mm stroke on a
   70 mm plate needs a cell near 0.14 mm everywhere. About 8% of the plate's
   cells carry surface (above), and most of those are on flat faces that a
   2 mm cell would represent exactly.
2. **What the corner signs cannot see, the mesher cannot mesh.** A feature
   that fits between two lattice corners without changing either one's sign
   is skipped silently. Today's guard is global (refine the whole lattice
   when the mesh comes out empty or open), and it cannot catch a feature
   that is lost while the rest of the part meshes fine. The design below
   does not change this: it finds what the uniform mesher at `cell_min`
   finds, at a lower triangle count, so `cell_min` can be set finer
   than is affordable today. A local guard is separate work.
3. **The simplifier is on the critical path.** Uniform meshes are 20 to 30
   times denser than they need to be, and simplifying them costs more than
   meshing them.

## 2. Data structure

A **pointer-free octree over the existing corner lattice.** A leaf at depth
`d` covers `2^(D-d)` finest cells per axis, where `D` is the depth of the
finest cell allowed. Leaves are keyed by `(d, i, j, k)` at their own depth.
Corner samples are keyed by finest-lattice integer coordinates, so a corner
shared by leaves of different sizes is sampled once and has one value.

Per leaf: its key, its corner values (eight floats, looked up rather than
stored), its sign configuration, and one QEF accumulator per surface
component, as the uniform mesher keeps today. The band of `cell_min` cells
R1 leaves is stored sparsely, keyed by finest-lattice index, rather than as
today's dense lattice.

The root is **cubic**: `2^D` finest cells across, with
`D = ceil(log2(n))` for `n` the largest per-axis cell count of the padded
lattice `_resolve_region` builds today. Blocks that start past the lattice
are dropped, and blocks that overhang it are clipped, so the lattice keeps
its per-axis counts and its padding (two finest cells) means what it means
today. On the plate, `D = 9`: a 512-cell root over a 494 x 145 x 32
lattice, and R1 prunes the slack in one or two levels.

The tree is **2:1 balanced**: neighbouring leaves differ by at most one
level. Balance costs some extra leaves along the fronts where merging stops,
typically 10 to 20%. In return, the dual traversal only ever handles a small, fixed
set of mixed-depth cases, and the per-component tables can be reused
unchanged on every leaf.

`NTree` (and the 2020 `yapcad/octtree.py` it builds on) stays an element
index. `sdf-meshing-and-simplifier-performance.md` suggested building on it;
this note does not, because the octree's payload, keys and traversal are all
field-specific, and a pointer-free tree walked a level at a time in key
order vectorises in numpy where an object tree does not. The new module is
`sdf/octree.py`; the similar name is unfortunate but the package keeps them
apart.

## 3. Refinement: sample fine, merge upward

The tree is built bottom-up, as in Ju et al. (2002): sample the band around
the surface at `cell_min`, then merge cells upward wherever the merge keeps
the surface's topology and stays within tolerance (§7, decision 6). Every
decision is made from finest-cell data, so the adaptive mesh never loses a
feature the uniform mesher at `cell_min` would find. Features thinner than
`cell_min` are lost exactly as they are today, and `to_solid`'s
refine-on-failure loop still answers them.

**R1. Occupancy (stage A, built).** For a block of half-diagonal `r` around
centre `c`, `|f(c)| > L * (r + h)`, with `L` the node's Lipschitz bound
from `analyze()` and `h = cell_min`, proves there is no surface within a
cell of the block's corners. Such blocks are filled with `f(c)`, which has
the right sign, and never refined. Everything else is refined to `cell_min`
and its corners sampled exactly, so the band is all that is ever sampled.
R1 costs one sample per tested block, at its centre; for an unclipped block
of two or more cells that centre is a lattice corner, which a later stage
can reuse instead of sampling twice.

**Merging.** Working up one level at a time, the up to eight children of a
parent cell merge into it when both tests below pass. Children R1 proved
empty carry no surface and pass trivially. A merged cell's QEF is the sum of
its children's accumulators (`ATA`, `ATb`, mass point and tally all add
exactly), so a merge needs no new samples and no new gradients.

**R3. Topology safety.** The test of Ju et al. (2002, after Gerstner and
Pajarola), on the 27 corner signs of the eight children:

- each child is itself a single leaf with at most one surface component;
- the sign at the middle of each of the parent's 12 edges matches at
  least one of that edge's ends;
- the sign at the middle of each of its 6 faces matches at least one of
  that face's corners;
- the sign at its centre matches at least one of its 8 corners;
- the parent's own sign configuration has at most one component.

At the first level the 27 points are exactly the `cell_min` corners inside
the parent, so a feature the uniform mesher would find has an inside
corner among them and fails the test; at higher levels each child already
passed it, which is why checking the 27 points at each level is enough.
Cells with several components keep the per-component vertices of the
uniform mesher and do not merge.

**R4. Geometric error.** The merged cell's QEF is solved as today. The merge
is refused if the residual divided by the tally exceeds `tol²`, or if the
solved vertex lies outside the merged cell after truncation. This is what
lets flat faces merge to `cell_max` while curved ones stop at the chordal
tolerance.

**Balance.** A merge that would leave a leaf more than one level coarser
than a neighbour is refused. Merges are decided for a whole level at once,
from the state at the start of that level, so the result does not depend
on order (§7, decision 3).

Merging stops at `cell_max`. Where no merge is possible the mesh is the
uniform mesh, so the adaptive mesh is never worse than the uniform one at
`cell_min`.

*Why there is no R2.* An earlier draft split leaves top-down and guarded
thin features with a hidden-crossing test on each leaf edge,
`|a| + |b| <= L h`. That bound is valid but fires next to every surface: an
edge parallel to a plane at distance `d` fires whenever `d <= L h / 2`, and
near a sphere of radius `r` an edge with closest approach `r + δ` fires when
`2δ + h²/4r <= L h`, which holds more easily as `h` shrinks. So it would
have sent every surface leaf to `cell_min`. Bisecting the edges fixes that,
but a feature inside a leaf that touches none of its edges is still lost.
Sampling fine and merging upward avoids both problems. What it gives up is
a reduction in sampling below stage A's band; what it gains is triangle
count, and with it the simplifier pass, which is what adaptive meshing is
for (§7, decision 4).

## 4. Polygonization

Dual contouring on an octree emits one quad per minimal crossing edge, with
corners at the dual vertices of the four leaves (three where sizes differ)
incident to that edge. Ju et al. (2002) find those edges with `cellProc`, `faceProc` and
`edgeProc`, a recursive traversal that descends to the smallest leaf
touching each edge, so no T-junctions and no cracks arise. With 2:1 balance
the mixed-depth cases are few and fixed, so each can be enumerated and
tested.

A non-recursive version suits numpy better. Enumerate minimal edges per
level, look up their incident leaves by key, and emit in bulk. A triangle
appears where only three distinct leaves meet an edge; that is the standard
result and needs no special casing.

The post-processing after emission stays as it is (better diagonal, weld,
flip, defect check). The weld's clamping box becomes the leaf's bounds
instead of the uniform cell.

## 5. Integration

- `dual_contour(node, ..., adaptive=False, cell_min=None, cell_max=None,
  tolerance=None)`. With `adaptive=False` nothing changes.
- `to_solid` passes `cell_mm` as `cell_min`. A new `tolerance` (chordal
  error, in model units) defaults to half of `cell_min`. `resolution`, if
  given, sets `cell_max`.
- The meshing record gains `"adaptive": {"coarsest": ..., "tolerance": ...}`
  beside `resolution`, in cells rather than model units, which replaces the
  `adaptiveDepth` that `SDF-DESIGN.md` §8.2 used to reserve (§7, decision 5).
- The refine-on-failure loop in `to_solid` applies unchanged: the adaptive
  mesh finds what the uniform mesh at `cell_min` finds, no more.
- `simplify_finished` stays available, but it does not run by default on
  adaptive output (§7, decision 4).
- Backends: corner sampling becomes a gather over integer keys per level,
  evaluated in one call per level, so numpy, numba and wgpu all apply
  unchanged. The gradient samples `_edge_normals` takes at each crossing
  are a second gather, once the crossings are known.

## 6. Staged plan

Each stage is independently useful and ends with tests that a later stage
must keep passing.

**A. Octree sampling, uniform output (built).** R1 subdivides the lattice
from a cubic root, and every leaf it cannot prove empty is refined to
`cell_min`, so no leaf stays coarse and nothing is left unmeshed. Proved
blocks fill their corners with the centre value, which has the right sign;
the remaining corners are sampled exactly. The corner lattice and
everything after it are the uniform mesher's, in its C order, so the
*invariant* -- output byte-identical to `dual_contour` at that cell -- holds
by construction. `tests/test_sdf_octree.py` checks it on the stage B
fixtures below, on numpy and numba. Stage A is also the first half of §3:
the band it samples is what stage B merges.

Measured, sampling alone (`sample_corners` against `sample_corners_octree`,
default leaf of 2 cells):

| case | lattice | chunks | octree | evaluations, chunks / octree |
|---|---|---|---|---|
| thick engraved block, 70 x 20 x 20 mm, res 350, numba | 3.9 M | 2.66 s | 0.41 s | 3.6 M / 1.4 M |
| same, numpy | 3.9 M | 18.4 s | 7.3 s | 3.6 M / 1.4 M |
| engraved plate above, res 490, numba | 2.4 M | 0.89 s | 0.52 s | 2.4 M / 1.5 M |
| sphere r = 50, res 400, numba | 66 M | 0.51 s | 1.02 s | 27 M / 5.4 M |
| box minus cylinder, res 400, numba | 32 M | 0.40 s | 0.83 s | 21 M / 6.3 M |

It wins where a sample is expensive -- the many-cutter text fields that
#94 is about -- and loses on primitive fields, where scattering writes into
a dense lattice costs more than the samples saved. So it is opt-in
(`dual_contour(..., sampling="octree")`, also accepted by `to_solid`), and
it does not change memory: the lattice is still dense. Stage B removes the
dense lattice, which is what makes both costs go away.

**B. Merging and balanced octree polygonization.** Store the band
sparsely, merge it upward with R3 and R4 (§3), implement the dual traversal
and emit mixed-size quads.
*Invariant:* with merging switched off, the output is byte-identical to
`dual_contour` at `cell_min`, as in stage A.
*Conformance fixture*, from the start: a box landing on and off lattice
planes, a sphere, a torus, a gyroid slab, a thin wall at 1.5 cells, corner
pinches, nested shells, a 0.5 mm pocket in a slab, and #94's engraved plate.
Each must be closed and edge-manifold, with no degenerate faces or
coincident vertices, and within a stated volume tolerance of the exact
value, at tolerances of 2, 1 and 0.5 times `cell_min`.

**C. The tolerance knob (R4), measured.**
*Measured:* triangle count against chordal error, compared with uniform
plus `simplify_finished` on the plate, the miter gear and the rover
cradle.

**D. Default.** Adaptive stays opt-in for now (§7, decision 1). Switching
`to_solid` to adaptive is a separate decision, made once C shows equal or
better error at lower cost and the opt-in mesher has been used in practice.

A rough expectation for the plate, not yet measured: leaves concentrate
along the strokes and the plate's edges, so the triangle count should land
near the simplified mesh's rather than near the uniform one's. Sampling
stays at stage A's cost, so the time saved is the simplifier's.

## 7. Decisions (settled 2026-11-21 and 22)

1. **Default behaviour: opt-in.** `adaptive=False` stays the default.
   Stage D is a later decision, made once the opt-in mesher has a track
   record, not automatically when stage C measures well.
2. **Knobs: keep `resolution` as a coarse cap.** It sets `cell_max`;
   `cell_min`, `cell_max` and `tolerance` are added alongside it and nothing
   is deprecated.
3. **Determinism: required.** For a given backend, the same node and
   parameters give the same mesh, byte for byte, as uniform meshes do.
   Balancing and every other pass run in key order, never in an order that
   depends on evaluation or thread scheduling.
4. **Simplifier: off by default on adaptive output.** Adaptive meshing is
   meant to be a single-pass alternative to uniform meshing plus
   `simplify_finished`, with a competitive triangle count. So
   `simplify_finished` does not run by default when `adaptive=True`; it stays
   available on request. Stage C's measurement is the acceptance test for
   this goal (within about 1.5 times the simplified count, at equal or better
   chordal error), not a gate on the default.
5. **Record format: cell bounds, in cells.** The record keeps `resolution`,
   cells across the longest axis at the finest cell, as it does today, and
   an adaptive mesh adds
   `"adaptive": {"coarsest": <cells across at cell_max>,
   "tolerance": <chordal tolerance as a fraction of the finest cell>}`.
   `resolution` must be `coarsest` times a power of two, and that power is
   the tree's depth, so `adaptiveDepth` is dropped from `SDF-DESIGN.md`
   §8.2. The bounds are recorded in cells rather than millimetres for the
   reason `resolution` already is (`convert.py`, `_deferred_solid`): it is
   what survives a transform, since a part scaled after meshing is still
   the same number of cells across, where its cell size in millimetres has
   changed. Nothing has written `adaptiveDepth` yet.
6. **Thin features: sample fine, merge upward** (settled 2026-11-22). The
   tree is built bottom-up from stage A's `cell_min` band, with the
   topology test of Ju et al. (R3) and the error test (R4) deciding merges,
   instead of top-down with a hidden-crossing split test. Top-down edge
   bisection was considered and rejected: it loses features inside a
   leaf that touch none of its edges. The cost is that sampling is not
   reduced below stage A's band.

## References

- T. Ju, F. Losasso, S. Schaefer, J. Warren. *Dual Contouring of Hermite
  Data.* SIGGRAPH 2002.
- S. Schaefer, T. Ju, J. Warren. *Manifold Dual Contouring.* IEEE TVCG 13(3),
  2007.
- S. Schaefer, J. Warren. *Dual Marching Cubes.* Pacific Graphics 2004: an
  alternative that contours the dual of the octree, worth a look if the
  mixed-depth cases in §4 prove troublesome.
