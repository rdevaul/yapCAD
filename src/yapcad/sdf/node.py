"""The signed-distance-function node DAG.

This module defines the representation that :doc:`SDF-DESIGN` makes
authoritative: an immutable, serialisable directed acyclic graph of nodes.
Everything else in :mod:`yapcad.sdf` is either a constructor that builds
nodes (:mod:`yapcad.sdf.primitives`, :mod:`yapcad.sdf.ops`) or a backend that
consumes them (:mod:`yapcad.sdf.evaluate`).

Why the DAG and not a sampled grid
----------------------------------

The tree serialises to kilobytes rather than the megabytes a base64 BREP
payload costs, it is re-evaluable after a parameter change, and — the reason
that matters most for the roadmap — *one* tree compiles to both a numpy CPU
evaluator and, later, a GLSL/WGSL fragment shader.  To keep those two
backends from diverging, a node kind does not carry its own evaluation code.
It registers a :class:`NodeSpec` whose ``backends`` mapping holds one callable
per backend name, so adding the shader emitter means adding a key to an
existing spec rather than writing a parallel module.

Lipschitz discipline
--------------------

Every node carries two *separate* properties, and conflating them is the
mistake the design document warns about:

``exact``
    True when the field is the true Euclidean signed distance.  Analytic
    primitives are exact; ``min``/``max`` CSG and smooth blends are not.

``lipschitz``
    A constant ``L`` such that ``|f(p)| <= L * d(p)``, where ``d`` is the true
    distance from ``p`` to the node's zero set.  A sphere tracer may therefore
    step ``|f(p)| / L`` safely, and a mesher may use ``L`` to bound cell
    subdivision.

The bound follows from ``f`` being ``L``-Lipschitz with its zero set equal to
the represented surface: for the closest surface point ``s``,
``|f(p)| = |f(p) - f(s)| <= L |p - s| = L d(p)``.  Each node's ``analyze``
function is responsible for returning a constant for which that holds; where
the derivation is not obvious it is written out at the constructor.

A third, structurally distinct flag rides along:

``csg_exact``
    True when the subtree contains only analytic primitives, hard booleans and
    similarity transforms, and can therefore be *replayed* through OCC to get
    an exact BREP (design document §5.2).  This says nothing about field
    fidelity; a hard boolean of two spheres is ``csg_exact`` but not ``exact``.
"""

import hashlib
import json
from dataclasses import dataclass
from functools import lru_cache

#: Serialised tree format tag, recorded in the geometry JSON representations
#: block (design document §8.2).
TREE_FORMAT = "yapcad-sdf-tree-v1"

INF = float("inf")

#: Bounds covering all of space, for fields with no finite extent.
UNBOUNDED = ((-INF, -INF, -INF), (INF, INF, INF))

#: The canonical empty bounding box, chosen so that unioning it with anything
#: is the identity.
EMPTY_BOUNDS = ((INF, INF, INF), (-INF, -INF, -INF))


class SdfError(ValueError):
    """Raised for a malformed node, tree or serialised document."""


# ---------------------------------------------------------------------------
# Field properties
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FieldProps:
    """The analysis result for a node.

    :param exact: the field is the true signed distance.
    :param lipschitz: constant ``L`` with ``|f(p)| <= L * d(p)``.
    :param bounds: axis-aligned ``((x0,y0,z0), (x1,y1,z1))`` surface bound.
    :param csg_exact: the subtree is replayable as exact CSG through OCC.
    """

    exact: bool
    lipschitz: float
    bounds: tuple
    csg_exact: bool

    def __post_init__(self):
        if not (self.lipschitz > 0.0):
            raise SdfError(
                f"lipschitz constant must be positive, got {self.lipschitz}"
            )


# ---------------------------------------------------------------------------
# Bounds arithmetic
# ---------------------------------------------------------------------------


def bounds_is_empty(b):
    """True if ``b`` encloses nothing."""
    lo, hi = b
    return any(lo[i] > hi[i] for i in range(3))


def bounds_union(a, b):
    """Smallest box containing both ``a`` and ``b``."""
    if bounds_is_empty(a):
        return b
    if bounds_is_empty(b):
        return a
    return (
        tuple(min(a[0][i], b[0][i]) for i in range(3)),
        tuple(max(a[1][i], b[1][i]) for i in range(3)),
    )


def bounds_intersect(a, b):
    """Largest box contained in both ``a`` and ``b``, possibly empty."""
    lo = tuple(max(a[0][i], b[0][i]) for i in range(3))
    hi = tuple(min(a[1][i], b[1][i]) for i in range(3))
    return (lo, hi)


def bounds_expand(b, d):
    """Grow ``b`` outward by ``d`` on every side (negative ``d`` shrinks).

    Shrinking never produces an inverted box; it collapses to the centre
    instead, which keeps the result a valid conservative bound.
    """
    if bounds_is_empty(b):
        return b
    lo = tuple(b[0][i] - d for i in range(3))
    hi = tuple(b[1][i] + d for i in range(3))
    fixed_lo = []
    fixed_hi = []
    for i in range(3):
        if lo[i] > hi[i]:
            mid = 0.5 * (b[0][i] + b[1][i])
            fixed_lo.append(mid)
            fixed_hi.append(mid)
        else:
            fixed_lo.append(lo[i])
            fixed_hi.append(hi[i])
    return (tuple(fixed_lo), tuple(fixed_hi))


def bounds_transform(b, matrix):
    """Axis-aligned bound of ``b`` carried through the 4x4 affine ``matrix``.

    Infinite extents stay infinite along every axis the transform can mix
    them into, which for a general affine map is all of them.
    """
    if bounds_is_empty(b):
        return b
    lo, hi = b
    if any(v in (INF, -INF) for v in lo + hi):
        return UNBOUNDED
    corners = []
    for ix in (0, 1):
        for iy in (0, 1):
            for iz in (0, 1):
                p = (b[ix][0], b[iy][1], b[iz][2])
                corners.append(_apply_affine(matrix, p))
    return (
        tuple(min(c[i] for c in corners) for i in range(3)),
        tuple(max(c[i] for c in corners) for i in range(3)),
    )


def _apply_affine(m, p):
    """Apply row-major 4x4 ``m`` to the 3-tuple ``p``."""
    return tuple(
        m[r][0] * p[0] + m[r][1] * p[1] + m[r][2] * p[2] + m[r][3]
        for r in range(3)
    )


# ---------------------------------------------------------------------------
# Node kind registry
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class NodeSpec:
    """Everything the system knows about one node kind.

    :param kind: the tag stored in the node and in serialised documents.
    :param min_children: minimum child count.
    :param max_children: maximum child count, or ``None`` for unbounded.
    :param analyze: ``(params, child_props) -> FieldProps``.  Backend
        independent; this is where Lipschitz and exactness propagation lives.
    :param backends: mapping of backend name to an evaluator.  A numpy
        evaluator has signature ``(params, points, children, ev) -> values``
        where ``ev(child, points)`` recursively evaluates a child at the given
        points — the callback form is what lets domain operators such as
        ``transform`` resample their child rather than merely combine values.
    """

    kind: str
    min_children: int
    max_children: object
    analyze: object
    backends: dict


_REGISTRY = {}


def register(spec):
    """Register ``spec``, replacing any previous spec for the same kind."""
    _REGISTRY[spec.kind] = spec
    return spec


def get_spec(kind):
    """Return the :class:`NodeSpec` for ``kind``."""
    try:
        return _REGISTRY[kind]
    except KeyError:
        raise SdfError(f"unknown SDF node kind {kind!r}") from None


def registered_kinds():
    """Return the sorted tuple of registered node kinds."""
    return tuple(sorted(_REGISTRY))


# ---------------------------------------------------------------------------
# Nodes
# ---------------------------------------------------------------------------


def _freeze(value):
    """Make ``value`` hashable and JSON-representable.

    Sequences become tuples so that a node can be a dictionary key and a
    digest can be taken over it.  Anything that is not a JSON scalar or a
    sequence of them is rejected outright: unlike a construction provenance
    record, an SDF parameter is load-bearing and must not be silently
    degraded to its ``repr``.
    """
    if value is None or isinstance(value, (bool, str)):
        return value
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        return value
    if isinstance(value, (list, tuple)):
        return tuple(_freeze(v) for v in value)
    raise SdfError(
        f"SDF node parameter must be a JSON scalar or sequence, "
        f"got {type(value).__name__}"
    )


def _thaw(value):
    """Inverse of :func:`_freeze`, for JSON emission."""
    if isinstance(value, tuple):
        return [_thaw(v) for v in value]
    return value


@dataclass(frozen=True)
class Node:
    """One node of the SDF DAG.

    Nodes are immutable and compare and hash by value, which is what makes
    the structure a DAG rather than a tree: a shared subtree is literally the
    same object, serialises once, and is evaluated once per point set.

    Build nodes with :func:`make_node` or, preferably, with the constructors
    in :mod:`yapcad.sdf.primitives` and :mod:`yapcad.sdf.ops`, which validate
    their arguments.
    """

    kind: str
    params: tuple
    children: tuple

    @property
    def p(self):
        """The parameters as an ordinary dictionary."""
        return dict(self.params)

    @property
    def props(self):
        """The :class:`FieldProps` for this node, computed once and cached."""
        return analyze(self)

    @property
    def exact(self):
        """True when this field is the true signed distance."""
        return self.props.exact

    @property
    def lipschitz(self):
        """The Lipschitz constant ``L`` with ``|f(p)| <= L * d(p)``."""
        return self.props.lipschitz

    @property
    def bounds(self):
        """Axis-aligned bound on the surface this node represents."""
        return self.props.bounds

    @property
    def csg_exact(self):
        """True when this subtree can be replayed exactly through OCC."""
        return self.props.csg_exact

    @property
    def digest(self):
        """Stable content hash, used as the node's serialised identifier."""
        return digest(self)

    def __repr__(self):
        shown = ", ".join(f"{k}={v!r}" for k, v in self.params)
        kids = f", {len(self.children)} children" if self.children else ""
        return f"<sdf {self.kind}({shown}){kids}>"


def make_node(kind, params=None, children=()):
    """Build and validate a :class:`Node`.

    Checks that the kind is registered and that the child count is legal;
    parameter validation is the constructor's job, because only the
    constructor knows the intended units and sign conventions.
    """
    spec = get_spec(kind)
    kids = tuple(children)
    for child in kids:
        if not isinstance(child, Node):
            raise SdfError(
                f"{kind}: children must be SDF nodes, got "
                f"{type(child).__name__}"
            )
    if len(kids) < spec.min_children:
        raise SdfError(
            f"{kind}: needs at least {spec.min_children} children, "
            f"got {len(kids)}"
        )
    if spec.max_children is not None and len(kids) > spec.max_children:
        raise SdfError(
            f"{kind}: accepts at most {spec.max_children} children, "
            f"got {len(kids)}"
        )
    given = (params or {}).items()
    items = tuple(sorted((str(k), _freeze(v)) for k, v in given))
    return Node(kind, items, kids)


@lru_cache(maxsize=None)
def analyze(node):
    """Return the :class:`FieldProps` for ``node``, recursively.

    Memoised on node identity, so a shared subtree in a DAG is analysed once.
    """
    spec = get_spec(node.kind)
    child_props = tuple(analyze(c) for c in node.children)
    props = spec.analyze(node.p, child_props)
    if not isinstance(props, FieldProps):
        raise SdfError(f"{node.kind}: analyze did not return FieldProps")
    return props


@lru_cache(maxsize=None)
def digest(node):
    """Content hash of ``node`` and its subtree.

    Merkle-style: a node's digest covers its children's digests, so equal
    subtrees get equal identifiers and serialisation deduplicates for free.
    The hash is ``sha256`` rather than Python's ``hash`` because it must be
    stable across processes — yapCAD signs packages, and an identifier that
    changed run to run would churn package hashes.
    """
    payload = {
        "kind": node.kind,
        "params": {k: _thaw(v) for k, v in node.params},
        "children": [digest(c) for c in node.children],
    }
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return "n" + hashlib.sha256(blob.encode("utf-8")).hexdigest()[:16]


def walk(node):
    """Yield every distinct node of the DAG, children before parents."""
    seen = set()

    def visit(n):
        key = digest(n)
        if key in seen:
            return
        for child in n.children:
            yield from visit(child)
        seen.add(key)
        yield n

    return visit(node)


# ---------------------------------------------------------------------------
# Serialisation
# ---------------------------------------------------------------------------


def tree_to_json(node):
    """Serialise the DAG rooted at ``node``.

    Emits a flat node table keyed by content digest plus a root reference,
    rather than a nested object, so that shared subtrees are stored once and
    the DAG structure survives the round trip.  A nested encoding would
    silently expand a DAG into a tree, which for a lattice of repeated cells
    is the difference between kilobytes and megabytes.
    """
    nodes = {}
    for n in walk(node):
        nodes[digest(n)] = {
            "kind": n.kind,
            "params": {k: _thaw(v) for k, v in n.params},
            "children": [digest(c) for c in n.children],
        }
    return {
        "format": TREE_FORMAT,
        "root": digest(node),
        "nodes": nodes,
    }


def tree_from_json(doc):
    """Rebuild a DAG from :func:`tree_to_json` output.

    Validates the format tag, resolves references, and rejects dangling
    references and cycles — a cycle would otherwise recurse forever in the
    evaluator.
    """
    if not isinstance(doc, dict):
        raise SdfError(
            f"SDF tree document must be a dict, got {type(doc).__name__}"
        )
    fmt = doc.get("format")
    if fmt != TREE_FORMAT:
        raise SdfError(
            f"unsupported SDF tree format {fmt!r}, expected {TREE_FORMAT!r}"
        )
    table = doc.get("nodes")
    root = doc.get("root")
    if not isinstance(table, dict):
        raise SdfError("SDF tree document is missing its 'nodes' table")
    if not isinstance(root, str):
        raise SdfError("SDF tree document is missing its 'root' reference")
    if root not in table:
        raise SdfError(f"SDF tree root {root!r} is not present in 'nodes'")

    built = {}
    building = set()

    def build(ref):
        if ref in built:
            return built[ref]
        if ref in building:
            raise SdfError(f"SDF tree contains a cycle through {ref!r}")
        entry = table.get(ref)
        if entry is None:
            raise SdfError(
                f"SDF tree reference {ref!r} is not present in 'nodes'"
            )
        if not isinstance(entry, dict):
            raise SdfError(f"SDF tree entry {ref!r} must be a dict")
        building.add(ref)
        kids = tuple(build(c) for c in entry.get("children", []))
        building.discard(ref)
        node = make_node(entry.get("kind"), entry.get("params") or {}, kids)
        built[ref] = node
        return node

    return build(root)


def to_construction(node, meshing=None):
    """Return an ``['sdf', tree]`` record for a solid's construction slot.

    :mod:`yapcad.construction` reserved the ``sdf`` kind tag in Phase 0; this
    is the function that fills it, so a solid meshed from a field carries the
    field that generated it through the package boundary.

    :param meshing: optional dictionary of the parameters used to generate
        the solid's mesh preview, appended as a third element.  Recording
        them is what makes the preview reproducible, which package signing
        depends on.
    """
    record = ["sdf", tree_to_json(node)]
    if meshing is not None:
        record.append(dict(meshing))
    return record


def from_construction(record):
    """Recover the DAG from an ``['sdf', tree]`` construction record.

    Returns ``None`` for a record of any other kind, so callers can probe a
    solid's provenance without first checking the tag.
    """
    if not is_sdf_construction(record):
        return None
    return tree_from_json(record[1])


def is_sdf_construction(record):
    """True if ``record`` is an SDF construction record."""
    if not isinstance(record, (list, tuple)) or len(record) < 2:
        return False
    return record[0] == "sdf"


def meshing_from_construction(record):
    """Return the meshing parameters of an SDF record, or ``None``.

    A record written before meshing parameters were carried, or one for a
    field that has not been meshed, simply has no third element.
    """
    if not is_sdf_construction(record) or len(record) < 3:
        return None
    params = record[2]
    return dict(params) if isinstance(params, dict) else None
