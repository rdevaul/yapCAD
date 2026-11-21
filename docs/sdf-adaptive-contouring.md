# Adaptive dual contouring: implementation design

Status: proposal, for discussion (#94). Nothing here is implemented.

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
- **Speed.** Sampling is Lipschitz-pruned by chunk and runs on numba or the
  tape backend. The 70 x 30 x 4 mm engraved plate of #94 meshes at a
  0.143 mm cell (about 2.9 million cells, 512,000 triangles) in about 3 s.

What remains is structural: every emitted cell is the same size.

1. **Thin features set the cell for the whole part.** A 0.43 mm stroke on a
   70 mm plate needs a cell near 0.14 mm everywhere. About 9% of the plate's
   cells carry surface, and most of those are on flat faces that a 2 mm cell
   would represent exactly.
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

The tree is **2:1 balanced**: neighbouring leaves differ by at most one
level. Balance costs some extra leaves along refinement fronts, typically
10 to 20%. In return, the dual traversal only ever handles a small, fixed
set of mixed-depth cases, and the per-component tables can be reused
unchanged on every leaf.

`NTree` stays an element index. It isn't a good fit here: the octree's
payload, keys and traversal are all field-specific, and a separate small
module, `sdf/octree.py`, will be easier to test exhaustively.

## 3. Refinement

Starting from a root cell, or a coarse uniform lattice of `cell_max` cells,
a leaf is split while any of these holds and it is larger than `cell_min`.
They are listed cheapest first, and all are batched by level.

**R1. Possibly occupied.** For a leaf of half-diagonal `r` around centre
`c`, `|f(c)| <= L * r` means the surface may pass through it. `L` is the
node's Lipschitz bound from `analyze()`. Leaves that fail the test are never
split, which is where nearly all the cost of a uniform grid goes.

**R2. Hidden crossing (thin features).** Corner signs miss a feature only if
the field dips through zero and back between two samples. Along an edge of
length `h` with end values `a` and `b` of the same sign, a Lipschitz field
can reach zero only if `|a| + |b| <= L * h`. If any edge of a leaf passes
that test, split the leaf. This is the test that guarantees a stroke
narrower than the leaf is found instead of stepped over, and it needs no
samples beyond the corners. It is conservative: near a smooth surface it
also fires for crossings that aren't there. So it is applied only above
`cell_min`, and as the cell shrinks toward `cell_min` it fires only on
genuinely thin material.

**R3. Topology.** Split a leaf if its sign configuration is one of the
ambiguous ones (face or body saddles), or if any of its components' QEF
vertices fall outside the leaf after truncation. Both mean that one vertex
per component cannot represent the surface there.

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
corners at the vertices of the four leaves (three where sizes differ) around
that edge. Ju et al. (2002) find those edges with `cellProc`, `faceProc` and
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
  "tolerance": ...}`. `SDF-DESIGN.md` §8.2 already reserves `adaptiveDepth`;
  one name should win.
- The refine-on-failure loop in `to_solid` becomes unnecessary for adaptive
  meshes: R1 and R2 find thin material directly. It stays for uniform meshes.
- `simplify_finished` stays available, but on adaptive output it should
  become a polish pass. Whether it runs by default is a measurement
  question; see §7.
- Backends: corner sampling becomes a gather over integer keys per level,
  evaluated in one call per level, so numpy, numba and the tape backend all
  apply unchanged.

## 6. Staged plan

Each stage is independently useful and ends with tests that a later stage
must keep passing.

**A. Octree build and refinement (R1 to R3), uniform output.** Build the
tree, then mesh each leaf region with the existing uniform mesher at
`cell_min`, only where leaves are at `cell_min`. This is a cheaper way to
spend the uniform mesher's effort, and it proves R2 on thin features
without any new polygonization.
*Invariant:* with refinement forced to `cell_min` everywhere, the output is
byte-identical to `dual_contour` at that cell.

**B. Balanced octree polygonization.** Implement the dual traversal and emit
mixed-size quads.
*Conformance fixture*, from the start: a box landing on and off lattice
planes, a sphere, a torus, a gyroid slab, a thin wall at 1.5 cells, corner
pinches, nested shells, a 0.5 mm pocket in a slab, and #94's engraved plate.
Each must be closed and edge-manifold, with no degenerate faces or
coincident vertices, and within a stated volume tolerance of the exact
value, at three tolerances.

**C. Error-driven refinement (R4) and the tolerance knob.**
*Measured:* triangle count against chordal error, compared with uniform
plus `simplify_finished` on the plate, the miter gear and the rover
cradle.

**D. Default.** Switch `to_solid` to adaptive if C shows equal or better
error at lower cost. Otherwise keep it opt-in.

A rough expectation for the plate, not yet measured: leaves concentrate
along the strokes and the plate's edges, so the sampled corners should be a
small fraction of the uniform lattice's 2.9 million cells, and the triangle
count should land near the simplified mesh's (8,784 in #94) rather than near
512,000. Stage A's numbers will replace this guess.

## 7. Decisions needed

1. **Default behaviour.** Should adaptive meshing become `to_solid`'s default
   once stage C measures well, or stay opt-in (as the GPU backend is)?
2. **Knobs.** Keep `resolution`, which counts cells across the bounding box,
   as a coarse cap, or deprecate it in favour of `cell_min`, `cell_max` and
   `tolerance` in model units?
3. **Determinism.** Uniform meshes are deterministic per backend today. Is
   the same guarantee required of adaptive meshes? It rules out
   order-dependent balancing, so the design above balances by key order.
4. **Simplifier.** If adaptive output is within, say, 1.5 times the
   simplified triangle count, should `simplify_finished` stop running by
   default?
5. **Record format.** `adaptiveDepth` (reserved in §8.2) or the explicit cell
   bounds in §5? Packages are signed, so this is a format decision rather
   than a detail.

## References

- T. Ju, F. Losasso, S. Schaefer, J. Warren. *Dual Contouring of Hermite
  Data.* SIGGRAPH 2002.
- S. Schaefer, T. Ju, J. Warren. *Manifold Dual Contouring.* IEEE TVCG 13(3),
  2007.
- S. Schaefer, J. Warren. *Dual Marching Cubes.* Pacific Graphics 2004: an
  alternative that contours the dual of the octree, worth a look if the
  mixed-depth cases in §4 prove troublesome.
