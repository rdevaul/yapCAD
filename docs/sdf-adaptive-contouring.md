# Adaptive dual contouring: implementation design

Status: stage A is implemented (`yapcad/sdf/octree.py`, opt-in through
`sampling="octree"`); stages B to D are proposals (#94). §3's R2 needs
redesign before stage B; see there.

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
   that is lost while the rest of the part meshes fine.
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
component, as the uniform mesher keeps today.

The root is **cubic**: `2^D` finest cells across, with
`D = ceil(log2(n))` for `n` the largest per-axis cell count of the padded
lattice `_resolve_region` builds today. Blocks that start past the lattice
are dropped, and blocks that overhang it are clipped, so the lattice keeps
its per-axis counts and its padding (two finest cells) means what it means
today. On the plate, `D = 9`: a 512-cell root over a 494 x 145 x 32
lattice, and R1 prunes the slack in one or two levels.

The tree is **2:1 balanced**: neighbouring leaves differ by at most one
level. Balance costs some extra leaves along refinement fronts, typically
10 to 20%. In return, the dual traversal only ever handles a small, fixed
set of mixed-depth cases, and the per-component tables can be reused
unchanged on every leaf.

`NTree` (and the 2020 `yapcad/octtree.py` it builds on) stays an element
index. `sdf-meshing-and-simplifier-performance.md` suggested building on it;
this note does not, because the octree's payload, keys and traversal are all
field-specific, and a pointer-free tree walked a level at a time in key
order vectorises in numpy where an object tree does not. The new module is
`sdf/octree.py`; the similar name is unfortunate but the package keeps them
apart.

## 3. Refinement

Starting from a root cell, or a coarse uniform lattice of `cell_max` cells,
a leaf is split while any of these holds and it is larger than `cell_min`.
They are listed cheapest first, and all are batched by level.

**R1. Possibly occupied.** For a leaf of half-diagonal `r` around centre
`c`, `|f(c)| <= L * r` means the surface may pass through it. `L` is the
node's Lipschitz bound from `analyze()`. Leaves that fail the test are never
split, which is where nearly all the cost of a uniform grid goes. Stage A
uses `|f(c)| <= L * (r + h)`, with `h` the finest cell, so that no lattice
edge leaving a proved leaf can cross either. R1 costs one sample per tested
leaf, at its centre; for an unclipped leaf of two or more cells that centre
is a lattice corner, which stage B can reuse instead of sampling twice.

**R2. Hidden crossing (thin features).** Corner signs miss a feature only if
the field dips through zero and back between two samples. Along an edge of
length `h` with end values `a` and `b` of the same sign, a Lipschitz field
can reach zero only if `|a| + |b| <= L * h`. If any edge of a leaf passes
that test, split the leaf. The bound is valid -- a necessary condition for
a hidden crossing -- and it needs no samples beyond the corners.

*As stated, R2 cannot be the gate it was meant to be.* It fires next to
every surface, not only on thin material. An edge parallel to a plane at
distance `d` has `a = b = d`, so it fires whenever `d <= L h / 2`, and every
leaf the plane crosses has a face within half a leaf of it. Near a sphere of
radius `r`, an edge whose closest approach is `r + δ` has
`a = b ≈ δ + h²/8r`, and fires when `2δ + h²/4r <= L h`, which holds more
easily, not less, as `h` shrinks. So R2 would drive every surface leaf to
`cell_min`, and R4's tolerance would never get a say. (Found while
building stage A, and independently by review.)

Options for stage B, none settled:

- *Gate on mixed signs* does not work: a feature hidden between corners
  leaves all of a leaf's corners the same sign, which is exactly when R2 is
  needed.
- *Bisect the edge.* Apply the same bound to the edge's halves, sampling the
  midpoint, and recurse to a floor `δ_min`. That stays rigorous and fires
  only where the surface comes within about `δ_min / 2` of an edge without
  crossing it, a band that is thin next to a leaf instead of a leaf wide.
  Cost: a few samples per near-surface edge.
- *Use the crossing normals.* The Hermite data says which way the surface
  faces; two same-sign ends whose gradients both point away from the edge
  are unlikely to hide a crossing. Cheaper, but heuristic, and fields with
  only a Lipschitz bound give no guarantee.
- *Drop it* and state that features thinner than the leaf R1 and R4 settle
  on may be lost, which is today's uniform behaviour at `cell_min`.

**R3. Topology.** Split a leaf if its sign configuration is one of the
ambiguous ones (face or body saddles), or if any of its components' QEF
vertices fall outside the leaf after truncation. Both mean that one vertex
per component cannot represent the surface there. R3 applies only above
`cell_min`; at `cell_min` the per-component vertices of the uniform mesher
handle those configurations, as they do today.

**R4. Geometric error.** Split if a component's QEF residual, divided by the
number of crossings it gathers, exceeds `tol^2`. This is what lets flat faces
stay coarse and curved faces refine to a chordal tolerance. It is optional in
stage A and required by stage D.

After refinement, a balancing pass splits any leaf more than one level
coarser than a neighbour.

Termination is guaranteed by `cell_min`. If R2 or R3 still fire at
`cell_min`, the field has a feature thinner than the finest cell allowed.
`to_solid` refuses with a message naming where, using today's wording.

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
- The meshing record gains `"adaptive": {"cellMin": ..., "cellMax": ...,
  "tolerance": ...}`, which replaces the `adaptiveDepth` that `SDF-DESIGN.md`
  §8.2 used to reserve (§7, decision 5).
- The refine-on-failure loop in `to_solid` becomes unnecessary for adaptive
  meshes: R1 and R2 find thin material directly. It stays for uniform meshes.
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
fixtures below, on numpy and numba. R2 and R3 are not part of stage A:
with every surface leaf at `cell_min` they would change nothing, and R2
needs the redesign in §3 first.

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

**B. Balanced octree polygonization.** Implement the dual traversal and emit
mixed-size quads.
*Conformance fixture*, from the start: a box landing on and off lattice
planes, a sphere, a torus, a gyroid slab, a thin wall at 1.5 cells, corner
pinches, nested shells, a 0.5 mm pocket in a slab, and #94's engraved plate.
Each must be closed and edge-manifold, with no degenerate faces or
coincident vertices, and within a stated volume tolerance of the exact
value, at tolerances of 2, 1 and 0.5 times `cell_min`.

**C. Error-driven refinement (R4) and the tolerance knob.**
*Measured:* triangle count against chordal error, compared with uniform
plus `simplify_finished` on the plate, the miter gear and the rover
cradle.

**D. Default.** Adaptive stays opt-in for now (§7, decision 1). Switching
`to_solid` to adaptive is a separate decision, made once C shows equal or
better error at lower cost and the opt-in mesher has been used in practice.

A rough expectation for the plate, not yet measured: leaves concentrate
along the strokes and the plate's edges, so the triangle count should land
near the simplified mesh's rather than near the uniform one's. That depends
on §3's R2 being redesigned; as stated, R2 sends every surface leaf to
`cell_min` and the count would not fall at all.

## 7. Decisions (settled 2026-11-21)

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
5. **Record format: explicit cell bounds.** The meshing record carries
   `"adaptive": {"cellMin": ..., "cellMax": ..., "tolerance": ...}`, and
   `adaptiveDepth` is dropped from `SDF-DESIGN.md` §8.2. The depth follows
   from the bounds, but the bounds do not follow from the depth, and nothing
   has written `adaptiveDepth` yet. `cellMax` must be `cellMin` times a power
   of two.

   *Reopened by review, for the owner to confirm.* The record keeps
   `resolution` in cells across the longest axis on purpose: "that is what
   survives a transform" (`convert.py`, `_deferred_solid`), since a part
   scaled after meshing is still the same number of cells across. Cell
   bounds and a tolerance in millimetres would not survive a scale. An
   invariant record that keeps the explicit-bounds intent is
   `"resolution"` for the finest cell, as today, plus
   `"adaptive": {"coarsest": <cells across at cell_max>,
   "tolerance": <fraction of the finest cell>}`. It carries the same
   information as the bounds and does not need a separate depth. Nothing
   writes an adaptive record before stage B, so this can be settled then.

## References

- T. Ju, F. Losasso, S. Schaefer, J. Warren. *Dual Contouring of Hermite
  Data.* SIGGRAPH 2002.
- S. Schaefer, T. Ju, J. Warren. *Manifold Dual Contouring.* IEEE TVCG 13(3),
  2007.
- S. Schaefer, J. Warren. *Dual Marching Cubes.* Pacific Graphics 2004: an
  alternative that contours the dual of the octree, worth a look if the
  mixed-depth cases in §4 prove troublesome.
