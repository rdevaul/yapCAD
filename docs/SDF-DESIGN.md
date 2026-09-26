# yapCAD SDF Support — Design Document

**Status:** Draft v0.1
**Author:** Rich DeVaul (with Claude)
**Date:** 2026-09-26
**Tracking branch:** `feature/sdf-representation`

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

- **Primitives:** `sphere`, `box`, `rounded_box`, `cylinder`, `capsule`,
  `torus`, `cone`, `half_space`, `gyroid`, `schwarz_p`
- **Combinators:** `union`, `intersect`, `subtract`, `smooth_union(k)`,
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

### 8.3 Slot-semantics defects

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

**Still open.** `geometry_from_json` builds
`['solid', shell_surfaces, voids, construction]`, placing voids in slot 2 —
which `solid()` and the `geom3d` module docstring both document as the
**material** slot. The in-memory solid structure has no voids slot at all, so
this is a genuine disagreement about what slot 2 means rather than a simple
bug: reconciling it means deciding whether solids acquire a real voids slot.
Both fields are almost always empty in practice, so nothing observably breaks
today, but this should be resolved before slot 3 acquires load-bearing SDF
content. Tracked separately from this plan.

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
| **1** | `yapcad/sdf/`: node DAG, numpy evaluator, Lipschitz tracking, primitives, hard CSG | Validate against `geom3d.signedFaceDistance` (`geom3d.py:254`) and analytic distances |
| **2** | Dual contouring → ordinary yapCAD solids with `authoritative: "sdf"` | The entire existing downstream works on SDF parts |
| **3** | CSG-exactness classifier + OCC tree replay | STEP export for the common case |
| **4** | BREP/mesh → SDF promotion via OCC distance queries + cached grid | Mixed-authority booleans |
| **5** | GLSL/WGSL emitter + browser sphere tracing | Real-time preview |
| **6** | Lattices, variable-thickness shells, distance-field fillets, topology-optimisation import | The features that justify the effort |

Phases 0–2 require nothing beyond numpy, which is already a hard dependency.

## 11. Open questions

- Should `sampled` grid sidecars use `.npy`, or a raw f32 blob with a JSON
  header? `.npy` is convenient but couples the format to numpy's versioning.
- Does node-provenance tagging (§2.1) survive dual contouring cleanly at
  smooth-blend boundaries, where no single node "wins"? Probably needs a
  blend-weight vector rather than a single label.
- Is `native_brep` (the pure-Python BREP topology graph, 2666 LOC) a viable
  Phase 3 target for CSG replay in the no-OCC install, or is OCC replay the
  only realistic path?
