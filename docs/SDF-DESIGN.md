# yapCAD SDF Support — Design Document

**Status:** Phases 0–3 implemented; Phases 4–6 planned
**Author:** Rich DeVaul (with Claude)
**Date:** 2026-09-26

---

## 1. Motivation

yapCAD's representational history has been additive: 2D analytic geometry, then
triangle-mesh 3D, then full BREP via the OpenCASCADE kernel. Signed distance
functions (SDF) are the natural next tier, and they buy three things that the
existing representations do not provide:

1. **Booleans that cannot fail.** CSG on an implicit field is arithmetic on
   scalars. There is no topology to corrupt, no sliver faces, no
   `BRepAlgoAPI` failure mode. This is the "unbreakable CAD" property that
   nTop has commercialised.

2. **Field-native features.** Lattices, gyroids, variable-thickness shells,
   distance-driven fillets everywhere at once, and topology-optimisation
   results are cheap in an SDF and expensive-to-impossible in a BREP.

3. **Solid modelling for the pure-Python install.** This is the strategically
   important one. Today a `pip install yapcad` user has no reliable boolean
   engine — the native triangle-mesh engine is a fallback that yapCAD 2.0
   proposes to deprecate, and OCC only ships through conda-forge. An SDF
   kernel written against numpy alone makes the PyPI wheel a genuinely
   capable solid modeller for the first time, with no compiled dependencies.

That third point reframes the pure-Python/OCC split from a limitation into two
coherent products: **SDF for generative and organic geometry that installs
anywhere; BREP for exact manufacturing geometry and STEP interchange.**

## 2. Position statement: SDF is a peer, not a new root

It is tempting to state the representation cascade as `SDF → BREP → mesh`,
with SDF as the new universal authority. This document deliberately rejects
that framing, for two reasons.

**Exactness.** Trimmed NURBS patches, exact tangency conditions and G2 blends
have no closed-form distance function. Promoting an imported STEP part to SDF
authority *destroys* information.

**Topological identity.** An SDF has no faces and no edges. yapCAD has built a
great deal on persistent BREP face/edge identity — `brep_edge_select` (892
LOC), `package/analysis/face_naming`, every assembly datum and mate, the
collision system. If SDF becomes the universal authority, "fillet this edge"
and "apply this load to that face" stop denoting anything.

**Therefore: authority is a per-solid property, not a global ranking.** The
`representations.authoritative` field in the geometry JSON schema already
encodes exactly this idea for `"brep"` and `"mesh"`. SDF becomes a third
value. A solid authored as SDF is SDF-authoritative and regenerates BREP and
mesh downward. A solid imported from STEP remains BREP-authoritative. Mixing
the two in one boolean becomes an explicit, designed operation (§7) rather
than something the hierarchy papers over.

### 2.1 The strongest argument for SDF is not booleans

When an SDF evaluates `min(a, b)`, it knows *which operand won* at every
sample point. Propagate that node identity through the meshing stage and every
output triangle carries the identity of the construction-tree node that
produced it. That yields face groups of the form "the cylindrical surface
contributed by node `bore_3`" which are **stable under parameter change by
construction**, because they are keyed to the construction tree rather than to
kernel-internal topology IDs.

This is a solution to the topological-naming problem that no BREP kernel has
solved well. Given yapCAD's investment in FEA face naming and assembly mates,
this is plausibly a larger win than unbreakable booleans, and the design
should aim at it deliberately rather than discovering it late.

## 3. What already exists in the codebase

Three pieces of existing architecture make this tractable:

| Existing asset | Location | Relevance |
|---|---|---|
| `representations` block with `authoritative` discriminator | `io/geometry_json.py:283` | The per-solid authority mechanism already exists for brep/mesh |
| Boolean engine registry + capability dispatch | `boolean/__init__.py`, `geom3d.py:521` | `solid_boolean` already selects an engine from operand representation |
| `construction` slot on every solid | `geom3d.py:774` | A solid is `['solid', surfaces, material, construction, metadata?]` |

The `construction` slot is the important one. An SDF **is** a construction
tree. `geom3d_util` already writes `['procedure', '<call string>']` records
into it. But `_serialize_solid` never reads slot 3, so provenance is
discarded at the package boundary. Making that slot round-trip is Phase 0 of
this plan and is independently valuable.

### 3.1 Prerequisite: the core is scalar Python

numpy appears in 17 modules, none of them `geom` or `geom3d`. The core
geometry idiom is point-at-a-time operations on Python lists. That is fine for
BREP (OCC does the arithmetic) and tolerable for meshes. It is fatal for SDF,
where the field is evaluated at 10^6–10^8 points per preview.

**The SDF subsystem must be array-vectorised from the first commit.** It
should not be written in, or grow out of, the existing scalar idiom. Every
evaluator signature takes and returns numpy arrays of shape `(N, 3)` and
`(N,)`.

## 4. Representation: tree-primary, grid-as-cache

### 4.1 The node DAG is authoritative

An SDF is a serialisable directed acyclic graph:

- **Primitives:** `sphere`, `box`, `rounded_box`, `cylinder`,
  `rounded_cylinder`, `capsule`, `torus`, `cone`, `half_space`, `gyroid`,
  `schwarz_p`
- **Combinators:** `union`, `compound`, `intersect`, `subtract`,
  `smooth_union(k)`,
  `smooth_subtract(k)`, `offset(d)`, `shell(t)`
- **Domain operators:** `transform(M)`, `twist`, `bend`, `repeat`, `mirror`
- **Escape hatch:** `sampled(grid_ref)`, `mesh_field(entity_ref)`

The tree is the source of truth rather than any sampled form, because:

- It serialises to kilobytes. The current BREP-as-base64 payloads run to
  megabytes — `globestand.ycpkg/geometry/primary.json` is predominantly base64.
- **One tree compiles to both a numpy CPU evaluator and a GLSL/WGSL fragment
  shader.** This is what makes the rendering problem (§6) tractable at all.
- It is parametric and re-evaluable. That is what "unbreakable" actually
  means in practice: change a radius, regenerate, never repair.
- It is isomorphic to the `construction` slot that already exists.

### 4.2 Sampled grids are a node type, not a rival format

Voxel/grid fields are required for mesh→SDF conversion, topology-optimisation
import, and caching of expensive subtrees. They enter the design as a
`sampled(grid_ref)` node that may appear anywhere in the tree, with the array
stored as a hashed binary sidecar in package `attachments` — never inline
base64. `octtree.py`'s existing `NTree` is a reasonable starting point for
adaptive storage.

### 4.3 Lipschitz discipline — decide this first

This is the detail that most implementations get wrong and cannot cheaply fix
later.

Naive CSG does not preserve the distance property. `max(-a, b)` yields a lower
*bound* on distance, not the distance; sphere tracing against it overshoots
and produces surface artifacts. Smooth blends violate it outright.

**Every node carries `exact: bool` and a Lipschitz constant, propagated up the
tree.** Ray marchers divide their step by that constant; meshing uses it to
bound cell subdivision. This must be in the node schema from the first commit,
because retrofitting it means touching every node type and every consumer.

## 5. SDF → BREP: what is and is not achievable

This is where the honest answer is less satisfying than the pitch, so it is
worth stating precisely.

### 5.1 SDF → mesh: solved

**Dual contouring on an adaptive octree**, not marching cubes. Dual contouring
preserves sharp features, and it requires surface normals — which an SDF
provides for free from the gradient. Roughly 400 lines, no scipy/skimage
dependency, keeps the subsystem inside the pure-Python + numpy tier.

Meshing must be **deterministic** given identical parameters. yapCAD signs
packages; a preview that regenerates differently run-to-run churns package
hashes. Beware hash-ordered vertex maps.

### 5.2 CSG-expressible trees: replay, do not convert

When a tree contains only analytic primitives, hard booleans and rigid
transforms, it is isomorphic to a CSG tree. In that case the correct move is
not to convert the field at all — **replay the same tree through OCC.**
`sphere(c, r)` → `BRepPrimAPI_MakeSphere`; `union` → `BRepAlgoAPI_Fuse`. Every
binding needed already exists in `yapcad.brep`. The result is an exact BREP,
cheaply, for the majority of real mechanical parts.

### 5.3 Field-only trees: refuse honestly

When the tree contains `smooth_union`, a gyroid, a twist, or an imported grid,
there is no good general answer:

- Sewing a triangle shell via `BRepBuilderAPI_Sewing` produces a *valid* BREP
  that is faceted, and therefore useless for downstream fillets or STEP
  consumers who expect analytic surfaces.
- Surface segmentation with NURBS fitting is research-grade, expensive, and
  fragile.

**The design decision is to not promise it.** Each tree is classified at build
time as `csgExact: true | false`, that classification is stored in the
representations block, and it is surfaced to the user *at authoring time* —
"this part uses `smooth_union`; it will not round-trip to STEP." STEP export
either refuses or emits a loudly-warned faceted approximation. A tool that
states its limits up front is more useful than one that degrades silently.

## 6. Rendering

- **Tier 0 — no GPU, works with what exists.** Evaluate the tree over a numpy
  grid, dual-contour, hand triangles to the existing pyglet/VTK/
  `service/core/tessellator.py` paths. Interactive-ish at 64³–128³. Cheap, and
  it unblocks the entire pipeline without any new rendering infrastructure.

- **Tier 1 — GPU sphere tracing.** Emit GLSL from the same tree and render a
  full-screen quad. Note that `pyproject.toml` pins `pyglet>=1.5,<2`, and
  pyglet 1.5's shader support is weak enough to effectively block this on the
  desktop. **The browser is the better target**: the yapCAD 2.0 vision already
  calls for a Three.js workbench, and a generated fragment shader against
  WebGL2 is precisely the natural shape of this work. SDF is a strong forcing
  function for accelerating the web workbench.

- **Tier 2 — compute shaders / CUDA / warp** for large-field meshing. Later,
  if ever.

The shader emitter must be a first-class backend of the *same* compiler that
emits numpy. If it is bolted on afterwards, the two backends will diverge and
`smooth_union` will come to mean two different things.

## 7. Mixed-authority booleans

Union of a BREP part and an SDF part is the interesting case. Two strategies:

- **(a) Promote BREP → SDF.** OCC can answer exact point-to-shape distance
  (`BRepExtrema_DistShapeShape`) and inside/outside classification
  (`BRepClass3d_SolidClassifier`). Together those form a true signed distance
  field — correct, and slow. Wrap it in a lazily-populated sampled grid with
  trilinear interpolation and it becomes usable.
- **(b) Demote SDF → mesh** and run an existing mesh boolean.

Preferred: (a), with (b) as fallback. The result is SDF-authoritative.
`solid_boolean` gains SDF cases in its existing dispatch (`geom3d.py:521`)
rather than acquiring a parallel code path.

## 8. Package format changes

### 8.1 Phase 0 — `construction` (additive to v0.2)

The `construction` slot becomes an optional, additive field on serialised
solids. Because it is purely additive and optional, **the schema ID does not
bump**: a v0.2 reader without the field ignores it, and a v0.2 reader with the
field defaults to `[]`. Bumping here would churn every existing package hash
for no compatibility benefit.

```json5
{
  "type": "solid",
  "construction": ["procedure", "yapcad.geom3d_util.makeRevolutionSolid(...)"]
}
```

### 8.2 v0.3 — SDF representations

The schema bumps to `yapcad-geometry-json-v0.3` when SDF lands, because
extending the `authoritative` enum is a breaking change for strict readers.

```json5
"representations": {
  "authoritative": "sdf",
  "sdf": {
    "role": "authoritative",
    "format": "yapcad-sdf-tree-v1",
    "lipschitz": "bound",
    "lipschitzConstant": 1.4,
    "csgExact": false,
    "bounds": [x0, y0, z0, x1, y1, z1],
    "tree": { }
  },
  "mesh": {
    "role": "preview",
    "generatedFrom": "sdf",
    "method": "dual-contouring",
    "resolution": 128,
    "adaptiveDepth": 4,
    "tolerance": 0.01
  }
}
```

Meshing parameters are recorded so that the preview is *reproducible*, which
matters for package signing.

### 8.3 Slot-semantics defects (resolved)

Two distinct problems were found in the `construction` slot while implementing
Phase 0.

**Fixed in Phase 0.** `solid()` used `if material == []` as its
"argument not yet supplied" sentinel when assigning positional list
arguments. Because an explicitly empty material list does not change that
test, the near-universal producer idiom
`solid(surfaces, [], construction)` assigned the construction record to the
**material** slot and left construction empty. Every `['procedure', ...]`
record written by `geom3d_util`, and every `['boolean', ...]` record written
by the boolean engines, was misfiled this way. `_serialize_solid` then read
slot 2 as voids and iterated the record's strings character by character,
emitting a spurious `"voids": [[], []]` into the document. `solid()` now
counts positional slots instead, which both routes the record to slot 3 and
removes the junk voids.

**Also fixed.** `geometry_from_json` built
`['solid', shell_surfaces, voids, construction]`, placing voids in slot 2 —
which `solid()` and the `geom3d` module docstring both document as the
**material** slot. `yapcad.geometry._retessellate_brep_solid` and
`service/core/tessellator` shared the same misreading.

The resolution was **not** to add a voids slot, because a yapCAD mesh solid
does not need one: a cavity is a shell whose faces wind inward, and the
existing math already handles it. A 4×4×4 cube with a 2×2×2 cavity reports
`issolidclosed() == True` and `volumeof() == 56.0`, because `volumeof`
accumulates *signed* tetrahedron contributions and takes the absolute value
only at the end. A voids slot would be a second encoding of a fact the
surfaces already carry, and two encodings that can disagree are exactly the
failure mode that produced this defect in the first place.

Instead the partition is **derived**: `geom3d.solid_shells(sld)` groups
surfaces into connected shells by shared edges and classifies each by the sign
of its volume, returning `(outer_shells, void_shells)`. Reads merge document
void surfaces into the shell list, where mesh solids keep cavities, and leave
slot 2 as material. Writes emit the derived partition, so `voids` remains
meaningful at the interchange layer for consumers that care (STEP, FEA)
without existing as state that can drift.

Two limits, both documented in the function: classification is meaningful only
for a closed shell, so an open group is always reported as outer; and
granularity is the surface, not the triangle, so a boundary that arrived as a
single surface — the usual shape of an OCC tessellation — yields one outer
shell even when it has cavities. Callers needing triangle-level partitioning
must do their own connectivity analysis.

**Resolved since.** `_retessellate_brep_solid` dropped the metadata dict at
slot 4, so a retessellated solid received a fresh entity id and lost its tags
and layer. Its only callers were `Geometry.mirror` and `Geometry.scale`, which
PR #57 changed to delegate to `geom3d`; the helper was then removed.

## 9. DSL surface

Add `SDF` to `dsl/types.py` alongside `SOLID`, with implicit `SDF → SOLID`
coercion via meshing so that every existing builtin keeps working unchanged.

Make `union` / `subtract` / `intersection` polymorphic: both-SDF operands stay
in SDF; mixed operands promote per §7. This preserves the author's existing
mental model.

New builtins (`sdf_sphere`, `sdf_box`, `smooth_union(a, b, k)`, `offset`,
`shell`, `gyroid`, `sdf_to_solid(sdf, resolution:)`, `solid_to_sdf`) belong in
a new `dsl/runtime/builtins_sdf.py`. `dsl/runtime/builtins.py` is already 4216
lines and should not absorb them.

## 10. Phased plan

| Phase | Work | Unlocks |
|---|---|---|
| **0** ✅ | Serialise and restore the `construction` slot; fix `solid()` slot assignment | Provenance survives packaging; prerequisite for all of the below |
| **1** ✅ | `yapcad/sdf/`: node DAG, numpy evaluator, Lipschitz tracking, primitives, hard CSG | Validated against `geom3d.signedFaceDistance` (`geom3d.py:254`) and analytic distances |
| **2** ✅ | Dual contouring → ordinary yapCAD solids with `authoritative: "sdf"` | The entire existing downstream works on SDF parts |
| **3** ✅ | CSG-exactness classifier + OCC tree replay | STEP export for the common case |
| **4** (in progress: §10.4, §10.6) | BREP/mesh → SDF promotion via OCC distance queries + cached grid | Mixed-authority booleans |
| **5** | GLSL/WGSL emitter + browser sphere tracing | Real-time preview |
| **6** | Lattices, variable-thickness shells, distance-field fillets, topology-optimisation import | The features that justify the effort |

Phases 0–2 require nothing beyond numpy, which is already a hard dependency.

### 10.1 Phase 1 as built

`yapcad/sdf/` ships as four modules: `node.py` (the DAG, its analysis and its
serialisation), `primitives.py`, `ops.py` and `evaluate.py` (the numpy
backend).  Four decisions are worth recording because later phases depend on
them.

**The backend seam is real, not notional.** §6 warns that a shader emitter
bolted on afterwards will diverge from the numpy evaluator. A node kind
therefore carries no evaluation code of its own: it registers a `NodeSpec`
whose `backends` mapping holds one callable per backend name, and
`evaluate.py` dispatches through that mapping. Phase 5 adds a `"glsl"` key to
the existing specs rather than a parallel module. The evaluator signature is
`(params, points, children, ev)` — the `ev` callback, rather than
pre-evaluated child values, is what lets a domain operator such as
`transform` resample its child instead of merely combining results.

**`exact` and `csg_exact` are separate flags, and both are needed.** A hard
boolean of two spheres is `csg_exact` (§5.2 will replay it through OCC for an
exact BREP) but not `exact` (`min(a, b)` is not the true distance inside the
union). Conflating them would either forbid STEP export for ordinary CSG or
claim distance fidelity the field does not have.

**The Lipschitz bookkeeping turned out tighter than expected.** `L` is
defined operationally: `|f(p)| <= L * d(p)`, so `|f| / L` is always a safe
ray-marching step. Two results fell out of deriving it properly:

- Hard CSG keeps `L = 1`. `min`/`max` of 1-Lipschitz functions is
  1-Lipschitz, and the zero set of `min(a, b)` *is* the union boundary, so
  the bound holds inside the solid as well as outside — the field is
  inexact there but never unsafe.
- The polynomial smooth minimum also keeps `L = max(La, Lb)`. Its partials
  work out to exactly `h` and `1 - h`, making the gradient a convex
  combination of the operands'. So a fillet destroys exactness without
  loosening the step bound, which is a better outcome than §4.3 anticipated.

A non-uniform `transform` is handled by scaling the child's field by the
matrix's *smallest* singular value. That underestimates distance by up to
`s_max / s_min` but leaves `L` unchanged, so a stretched primitive stays safe
to march; it is marked `exact = False` accordingly.

`gyroid` and `schwarz_p` are in Phase 1 on purpose despite belonging to the
Phase 6 feature set. They are the only nodes with `L > 1`, and without them
the Lipschitz machinery would be propagating `1.0` everywhere and no test
could tell correct bookkeeping from vacuous bookkeeping. Their `thickness` is
in field units, not millimetres, and this is documented at the constructor.

**Digests are content hashes, not object identities.** Each node's
identifier is a truncated sha256 over its kind, parameters and its children's
digests. That gives DAG deduplication on serialisation for free, and — the
reason for sha256 over Python's `hash` — it is stable across processes, so an
SDF tree stored in a package does not churn the package hash from run to run.
The serialised form is a flat node table plus a root reference rather than a
nested object, because a nested encoding silently expands a DAG back into a
tree.

**Deferred, and not yet started.** The `twist`, `bend` and `repeat` domain
operators and the `sampled(grid_ref)` escape hatch of §4.2 are not
implemented. `repeat` in particular needs care: the usual `mod`-based
implementation overestimates distance whenever the nearest instance lies in a
neighbouring cell, which breaks the marching guarantee that the rest of the
subsystem maintains, so it wants the multi-cell candidate form rather than
the naive one.

### 10.2 Phase 2 as built

`yapcad/sdf/contour.py` implements dual contouring and `convert.py` wraps the
result as an ordinary `['solid', [surface], material, construction]`, so
drawables, exporters, the collision system and the package writer need no
knowledge of fields at all. The schema moves to
`yapcad-geometry-json-v0.3`, with `docs/schemas/yapcad-geometry-json-v0.3.schema.json`
published alongside the v0.2 file; v0.1 and v0.2 documents still read.

**The sharp-feature claim is measured, not asserted.** A lattice-aligned box
meshes to a volume of 512.0 against an exact 512.0, to about 1e-8 relative —
and the residual is the finite-difference gradient feeding the QEF, not the
reconstruction. Marching cubes would chamfer every edge and corner and lose
volume at the percent level. Topology is checked the same way: a torus comes
out with Euler characteristic 0, and a box minus a sphere large enough to
breach all six faces comes out at -8, the genus-5 cube frame it should be.

**Adaptivity is in the work, not yet in the output.** The octree of §5.1 is
present as Lipschitz-pruned chunked sampling: a chunk whose centre satisfies
`|f| > L * (circumradius + margin)` provably has no surface within `margin`,
so with `margin` at least one cell diagonal no edge touching that chunk can
cross, and the chunk can be filled with its centre value instead of sampled.
Empty space therefore costs almost nothing, and peak memory is bounded by the
chunk rather than the lattice. Every emitted cell is still the same size,
though. Genuinely adaptive *output* — variable leaf depth with QEF-error
collapse — needs the crack-patching cell/face/edge traversal of Ju et al. and
is deferred. The tests assert that pruned and unpruned meshing agree
bit-for-bit, which is what makes the optimisation safe to trust.

**One vertex per cell is a real limit, and it is now detected.** Dual
contouring cannot represent two surface sheets passing through one cell, so
a wall or a gap thinner than a cell yields non-manifold edges: a gyroid with
a 0.37-unit wall meshed at a 0.63-unit cell produces 42 of them. Manifold
Dual Contouring (Schaefer and Ju) fixes this by splitting such a cell into
one vertex per surface component; it is not implemented. Following §5.3's
own argument that a tool which states its limits beats one that degrades
silently, `to_solid` checks edge multiplicity by default and refuses with a
message naming the cause and the cure, with `check=False` to override.
Handing a silently non-manifold solid to `volumeof` or to STEP export would
be the worse failure.

**Determinism is treated as a requirement, not a hope.** Vertex order is the
C-order of the active-cell mask and triangle order follows the edge axes in
turn; nothing iterates a dict or a set. Tests assert that repeated runs,
separate processes, every chunk size and both pruning settings produce
byte-identical arrays, and that the recorded meshing parameters regenerate
the same mesh.

**The tree lives in exactly one place in a document.** An SDF-authoritative
solid puts its tree in `representations.sdf.tree` and omits the top-level
`construction` field; a document carrying both is rejected. This is §8.3's
lesson from the voids defect applied before the fact rather than after: two
encodings of one fact will eventually disagree. The analysis fields
alongside the tree — `lipschitz`, `lipschitzConstant`, `csgExact` — are a
cache of what the tree already determines, so they are recomputed and
verified on read rather than believed, the same posture as the BREP payload
hash check. A document claiming `csgExact` for a gyroid is rejected.

**Authority is enforced as single-valued.** `prism()` attaches an analytic
BREP, so a prism is BREP-authoritative; attaching an SDF record to one is a
contradiction and is refused at serialisation. The consequence worth
flagging: an SDF solid that later acquires a cached BREP — from an OCC
boolean, say — will fail to serialise, because the v0.3 representations
block still admits a BREP only when it is authoritative. A *derived* BREP
role is what Phase 3 and Phase 4 need, and is the natural place to resolve
this.

**Not implemented in this phase.** Node-provenance tagging (§2.1) is not
threaded through the mesher. The mesher's structure does not preclude it —
the QEF already accumulates per-cell contributions from identifiable edges —
but the §11 question about what a node label means at a smooth-blend
boundary has to be answered first, and answering it with a blend-weight
vector changes the surface-group data model rather than just the mesher.

### 10.3 Phase 3 as built

`yapcad/sdf/occ.py` replays `csgExact` trees through OpenCASCADE. Each
analytic node kind registers an `"occ"` backend on its existing `NodeSpec`
through `register_backend` -- the seam of §10.1 used for a second time, and
the reason the module can live apart so that `import yapcad.sdf` never needs
pythonocc-core. `sdf.to_solid(node, brep=True)` attaches the replayed shape,
and `write_step_analytic` then emits planes, cylinders, spheres, cones and
tori rather than facets. `brep="auto"` attaches one when the tree allows it
and quietly skips it when not.

**"Exact" is tested literally.** Replayed volumes are compared with
closed-form volumes, not with the mesh: primitives, rotated, uniformly
scaled and mirrored parts, and half-space cuts agree to about 1e-16 relative,
the cube frame to 1.1e-7 of `1000 - 254*pi`, a rounded box to 1.6e-9 of its
Minkowski volume. The analytic cube-frame STEP is 46 KB against 11 MB for
the same part as STL.

**Two replays are not one-to-one calls.** A `rounded_box` is a box with
every edge filleted to the same radius, which is exactly the Minkowski
rounding the field describes -- no edge is ever skipped, because a box with
one sharp edge is a different part. A `half_space` is replayed as a finite
box on the correct side of its plane, sized to cover the root's bounds as
carried into the leaf's frame through any transforms above it. Replacing a
leaf with anything that agrees with it inside that region cannot change the
result inside it, and the result lies inside it by construction; the driver
also checks the replayed shape against the root's bounds afterwards, so a
violated assumption is loud. That check needed `Bnd_Box.AddOptimal`: the
default bounds a trimmed face by its whole underlying surface, so a sphere
cut down by a box reported the sphere's full extent.

**The classifier was wrong in two places, and replay found both.** A
spindle torus (minor radius >= major) is a perfectly good field, but its BREP
surface self-intersects; OCC builds it, calls it valid, and counts the
overlap twice in its volume. A rounded box whose radius is half its smallest
edge leaves the fillet no face to run along, and OCC refuses it. Both are now
`csgExact: false`, which is the promise §5.3 makes: the classification is
what authoring time relies on, so it must not promise a replay that fails.

**Blockers are reported, not just a boolean.** `sdf.csg_blockers(node)`
lists the nodes that are non-CSG *in their own right*, found by re-running
each node's analysis with its children forced CSG-exact, so a `smooth_union`
whose operand is also a lattice reports both. A refused replay names them,
and the gallery's README lists them per part -- the authoring-time "this part
uses `smooth_union`; it will not round-trip to STEP" that §5.3 asks for.

**A BREP beside an SDF tree is now either derived or refused.** The Phase 2
blocker of §10.2 is resolved with a *derived* role rather than by relaxing
authority. `to_solid(brep=True)` tags the attached BREP with the digest of
the tree it was replayed from. The serialiser accepts a BREP beside SDF
authority only when that tag matches the tree it is writing, and emits it as
`role: "derived"` with `treeDigest`; anything else -- `prism()`'s own analytic
BREP, or one replayed from an earlier version of the tree -- is still
refused as a second, disagreeing definition. The reader compares
`treeDigest` with the digest of the tree *as rebuilt*, not with the table
key the document gives the root, so relabelling both together does not get
past it. The v0.3 schema was amended rather than bumped for this, since it
has not shipped outside this branch.

**Transforms now carry the tree, which fixed a Phase 2 bug.** The solid
transforms in `geom3d` move a solid's surfaces and then call a hook per
representation (`translate_brep_solid`, `translate_native_brep`). There was
no hook for the SDF tree, so translating an SDF solid moved its mesh and left
the tree that defines it where it was: the document described a sphere at the
origin carrying a mesh a hundred units away. `yapcad/sdf/sync.py` is that
hook, called from `translatesolid`, `rotatesolid`, `mirrorsolid` and
`scalesolid`. It wraps the tree in the equivalent `transform`, drops the
meshing parameters (a mesh moved after meshing is not what they regenerate),
and re-tags a derived BREP when the BREP hook moved it too -- replaying
`transform(tree, M)` is literally the replayed shape carried through `M`, so
the tag stays true. A non-uniform scale leaves the BREP behind, so there it
is dropped, which is safe for a representation that is derived by
definition.

**Core defects found along the way, fixed separately in PR #57.** `geom.scale`
and `geom3d.scalesurface` composed a centred scale as `T(-c) S T(c)`, which
scales about `-c`, while rotation beside them and both BREP hooks used the
correct `T(c) S T(-c)`; the SDF hook followed the documented behaviour, so
until the fix a centred scale of an SDF solid left tree and mesh
disagreeing. `geom3d.scalesolid` also left a BREP attached and unscaled
after a non-uniform scale, where it was then read as authoritative. And the
`Geometry` wrapper's own solid `mirror` was a silent no-op on a mesh-only
solid, while its solid `scale` produced an invalid one. None of these was
SDF-specific, so they went to `main` as their own change and were merged
into this branch afterwards.

**Not implemented.** Mixed-authority booleans remain Phase 4: `solid_boolean`
on two `brep=True` SDF solids takes the existing OCC path and yields a
BREP-authoritative result, which is coherent but discards the tree.

### 10.4 Phase 4, first slice: booleans between SDF solids

§9 asks for booleans to be polymorphic, with both-SDF operands staying in
SDF. That is now what `geom3d.solid_boolean` does: when both operands carry
an SDF tree, `yapcad/sdf/booleans.py` combines the trees and meshes the
result, which stays SDF-authoritative and can be combined, replayed or
re-meshed again. Every yapCAD boolean goes through `solid_boolean`, so the
DSL's `union`/`difference`/`intersection`, the manufacturing helpers and
`text3d` all get it. An explicit `engine=` or `YAPCAD_BOOLEAN_ENGINE` still
wins, and `engine='sdf'` names the field path directly. Mixed operands --
one field, one BREP or mesh -- still go to the existing engines; promoting
them to fields is the rest of this phase.

**This fixed booleans on SDF parts, which were effectively broken.** Handed
two dual-contoured meshes, the native mesh engine returned an open mesh in
every case measured: two boxes sharing a face, a coplanar cut, a
sphere-box intersection. It took 7-40 s at resolution 12, timed out at 90 s
or raised on a degenerate face at resolution 20, and stalled for over ten
minutes at 32. The field path returns closed solids with the right volume
in well under a second: 1999.995 for the face-sharing union against 2000,
840.000 for the coplanar cut exactly.

**The one real decision is the resolution the result is meshed at, and the
rule is: never coarser than either operand.** Each operand's cell size is
recovered from its recorded meshing parameters -- dual contouring's spacing
is the region's longest extent over the resolution -- or, for a transformed
solid that has dropped them, estimated as its median mesh edge, which for a
dual-contoured mesh is about one cell. The smaller cell is carried over the
result's extent, so two boxes side by side mesh at twice the resolution
rather than at half the detail. `MAX_RESOLUTION` caps it at 256, since the
corner lattice is `(n+1)^3` doubles and a tiny detailed part unioned with a
large one would otherwise ask for gigabytes.

**A combination can be thinner than its operands**, when it brings two
surfaces nearly together, and then trips the one-vertex-per-cell limit.
`to_solid` now raises a dedicated `NonManifoldMeshError` (still an
`SdfError`), so the boolean path can catch exactly that, retry once at
double resolution, and otherwise let the refusal stand rather than return a
broken solid.

**Exact stays exact.** If every operand carries a BREP derived from its own
tree, the result is replayed too whenever the combined tree is still
CSG-exact, so an intersection of two `brep=True` parts exports analytic
STEP. One inexact operand, or a non-replayable result, means no BREP -- the
tag rules of §10.3 are never bent to keep one.

**The DSL meshes once.** Its booleans are variadic but folded pairwise,
which for fields would mesh and discard every intermediate result. When all
operands are SDF-authored they are combined n-ary with `combine_all` and
meshed once; a five-part `union` is a single `union` node with five
children. Mixed operands keep the pairwise path unchanged.

**An empty result is an empty solid that keeps its tree.** A disjoint
intersection returns `['solid', [], [], ['sdf', tree]]`, which yapCAD
already allows for empty CSG results and which serialises and reads back.

### 10.7 Simplifying against the field

Dual contouring meshes a part at the cell size its finest feature needs,
everywhere: an M8 nut fine enough for its 0.16 mm crest flats is just as
fine across its flat faces. `sdf.simplify_solid(solid, tolerance)` recovers
that by quadric edge collapse (Garland–Heckbert), checked against the
field. It runs in rounds. Each scores every edge at once, picks greedily a
set of the cheapest edges with disjoint neighbourhoods, moves each merged
vertex onto the surface by Newton steps, and tests them all together:

- **Topology** by the link condition, so the mesh stays a closed manifold.
- **Orientation**: no surrounding triangle flips or collapses.
- **Fidelity**: `|f|` is sampled across every changed triangle on a grid
  of two points per original edge length and at least four divisions, so a
  large merged triangle is checked more densely than the small ones it
  replaces. A sparser grid (one point per edge, two divisions minimum)
  read 0.095 mm on a triangle cut across a sharp edge whose true error was
  0.142 mm; with it the worst error over the whole part rose about 1.5%.
  Each triangle's error is cached and remeasured only when a collapse
  rewrites it. It must stay within the
  tolerance — or, where dual contouring already strayed further (along
  sharp edges it rounds), the worst error may not rise and the excess over
  the tolerance, integrated over area, may not grow. Without that second
  rule a sharp edge either froze every collapse near it or, with a looser
  rule, licensed its neighbours to reach its error, and the bad region
  spread a ring per round.

Rejected edges are remembered until something near them changes, which
keeps the tail of rounds cheap. It is deterministic, and the tolerance is
recorded with the meshing parameters.

At 0.01 mm, measured by sampling uniformly by area (sampling per triangle
over-weights the small triangles left at sharp edges and misleads):

| Part | Triangles kept | p99 `\|f\|`, before → after | Area over 0.01 mm | Volume |
|---|---|---|---|---|
| Plate with two holes | 9.9% | 0.0121 → 0.0112 mm | 1.7% → 1.5% | −0.03% |
| Threaded M8 nut | 3.8% | 0.0157 → 0.0086 mm | 2.1% → 0.7% | +0.12% |
| Miter gear, z24 | 5.9% | 0.0276 → 0.0153 mm | 4.7% → 1.8% | −0.09% |

The result is *more* accurate than the uniform mesh, because merged
vertices land on the surface. Its cost is dominated by field evaluation, so
it inherits the field's speed: about 6 s for the plate and 18 s for the nut,
but minutes for the gear, whose polygon field is the slow part.
Meshing adaptively in the first place — refining only where the field says
so — would avoid building the uniform mesh at all, and is the larger
follow-on.
### 10.6 Gears and threaded fasteners

The last OCC dependency in yapRover's release design is its differential's
`miter_gear`, built as an OCC ruled loft. Gears and fasteners are, at
heart, a 2D outline carried through space, so three planar kinds carry
most of the weight:

- **`polygon`** — the exact signed distance to a closed 2D polygon,
  constant along z. An `n`-fold symmetric polygon (a gear) rotates each
  point into one sector and tests only the segments within reach of it,
  which is exact and about `n/3` times cheaper. The sign comes from a ray
  cast radially outward, so it stays inside the sector, with a half-open
  straddle rule so that a ray through a vertex counts once.
- **`extrude`** — exact extrusion of a planar child about z = 0.
- **`apex_extrude`** — the cone through the origin over a profile given at
  `z_ref`, clipped between two planes: every section is the profile scaled
  by `z / z_ref`. Not exact, with a conservative Lipschitz bound.

On these sit *semantic* nodes, whose parameters are the part's
specification rather than its geometry, so trees stay a few hundred bytes
and re-evaluate from the same numbers the BREP generators use:

- **`straight_bevel_gear`** — the outer tooth section (an exact polygon,
  17 samples per flank) carried to the pitch apex between the inner and
  outer planes, minus the bore. That is the BREP generator's own
  construction, since its two loft sections are the same outline scaled
  toward the apex. It replays through OCC by calling that generator, so it
  is `csg_exact`: the field and BREP differ only between flank samples,
  straight segments against an interpolating spline, by under 5 µm.
  `make_straight_bevel_gear_sdf` meshes it at a cell set by the top land at
  the small end; the DSL's `miter_gear` and `straight_bevel_gear` use it
  when OCC is absent.
- **`thread`** — a helical thread about z: the distance in
  `(u, r)`, `u = z − hand·lead·θ/2π` folded by the pitch, to the profile
  `r = R(u)` over three periods. The helical map stretches space by
  `sqrt(1 + (lead/2πr)²)`, 1.002 at an M8 minor radius, so the field is
  very nearly exact on the flanks. That stretch is unbounded at the axis,
  so inside half the minor radius the angular dependence is blended out.
  The profile comes from `threadgen._radius_at` itself, so the field nut
  is the same part as the mesh nut.
- **`hex_nut`** — a hexagonal prism less an internal `thread`, placed as
  `fasteners_legacy.build_hex_nut` places it. `metric_hex_nut(...,
  representation="sdf")` (and the unified equivalent) meshes it closed,
  where the swept-mesh nut is not; the mesh nut stays the default.

Not yet done: external threads and bolts, spur, helical and herringbone
gears (polygon plus a `twist` domain operator), and a DSL switch to author
fasteners as fields. The investigation behind this section measured those
too; they mesh closed within 0.1% of reference.
### 10.5 Fillets and compounds

Building the yapRover release design without OCC (field primitives in
place of the DSL's mesh ones) found two gaps that blocked it outright, and
both are closed for SDF-authored solids.

**`fillet` works on fields.** `sdf.fillet(node, r)` rounds every edge, the
meaning of the DSL's `fillet` and of OCC's fillet-all-edges. For the
analytic primitives it is exact and closed-form: a `box` becomes a
`rounded_box`, a `cylinder` the new `rounded_cylinder` (an exact field: the
core cylinder shrunk by `r` and dilated back). It passes through similarity
transforms, whose uniform scale divides the radius, and through each body of
a `compound`; spheres, tori, capsules and already-rounded primitives have no
edges and are returned unchanged. Both rounded primitives replay through
OCC as filleted primitives, so a filleted part still exports analytic STEP.
`sdf.fillet_solid` re-meshes no coarser than its input, and the DSL's
`fillet` takes this path for any SDF-authored solid, with no OCC needed.

It refuses, rather than approximates, three things:

- **Combined fields.** Rounding each operand and blending each boolean
  also rounds edges the boolean removes: two boxes sharing a face come out
  with a groove along the seam, where OCC's fuse-then-fillet leaves a clean
  box. Circular blend operators have the same flaw at the seam. A correct
  fillet of a combination needs the distance to the *combined* part, which
  hard CSG does not give inside the part; that needs re-distancing on a
  sampled grid, so it waits for the `sampled` node of §4.2. Every
  yapRover fillet is of a primitive, filleted before it is combined.
- **Non-uniform scales**, whose rounds would be elliptical.
- **Cones and the remaining kinds**, for now: a rounded frustum is a
  closed form worth adding when a design needs it.

**Mesh-only solids are not filleted yet.** The DSL's `fillet` says so
rather than claiming OCC is the only route. The plan, in two steps:

1. **Basic shapes by promotion.** `prism`, `conic` and `sphere` record
   their call — `['procedure', 'yapcad.geom3d_util.prism(1,2,3,...)']` —
   which is enough to rebuild the corresponding field primitive, fillet it
   and mesh it, giving the same SDF-authoritative part as authoring it as a
   field. One prerequisite: today a transform leaves that record unchanged
   (`translatesolid(prism(1, 2, 3), ...)` still says `prism(1,2,3,...)`), so
   the record is stale after a move. The solid transforms must first
   append the transform to the record, as they already do for SDF trees,
   or drop it.
2. **Arbitrary meshes by sampling.** Sample the mesh's signed distance onto
   a grid (generalized winding number for the sign, as the native boolean
   engine already computes it), then fillet by re-distancing: opening by
   `r` rounds convex edges, closing by `r` rounds concave ones. This is the
   same machinery as the general fillet of combined fields, so the two
   arrive together with the `sampled` node.

**`compound` keeps its fields.** The DSL's `compound` concatenated its
operands' meshes and dropped their trees, so every later boolean against a
compound fell back to the mesh engine; in yapRover one such difference took
154 s against a 730,000-triangle operand. There is now a `compound` node:
as a field it is a union, so booleans combine it exactly like one, but it
records that the bodies are separate. `to_solid` meshes each body on its
own at the cell size the whole would have had — pushing any transforms
above the compound down onto the bodies — so touching bodies stay distinct,
closed meshes, and the OCC replay builds a compound rather than a fused
solid. The DSL's `compound` attaches this tree when every operand is
SDF-authored and keeps the operands' own meshes; a compound with any mesh or
BREP operand is unchanged.

## 11. Open questions

- Should `sampled` grid sidecars use `.npy`, or a raw f32 blob with a JSON
  header? `.npy` is convenient but couples the format to numpy's versioning.
- Does node-provenance tagging (§2.1) survive dual contouring cleanly at
  smooth-blend boundaries, where no single node "wins"? Probably needs a
  blend-weight vector rather than a single label.
- Is `native_brep` (the pure-Python BREP topology graph, 2666 LOC) a viable
  Phase 3 target for CSG replay in the no-OCC install, or is OCC replay the
  only realistic path?
