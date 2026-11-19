"""A field as a tape: the node DAG flattened for a compiled backend.

The numpy backend in :mod:`yapcad.sdf.evaluate` walks the tree node by node,
one array pass each.  That is the reference, and it is slow in exactly the
way that matters for a preview: a 500-node part costs about 19 microseconds
a point, and a frame needs several hundred thousand points.

A tape is the same tree as flat arrays -- an opcode, a destination register,
operand registers and a parameter offset per instruction -- which a single
fixed loop can run for every point.  Two things follow from the loop being
fixed.

**Nothing is compiled per design.**  Generating one fused kernel per tree is
faster still per point, but compiling it took 153 seconds for a 500-node
part, and GPU shader compilers treat large generated kernels the same way.
The interpreter is compiled once, cached, and reused for every tree.

**The backend is the loop and nothing else.**  :func:`flatten` knows the
node kinds; a backend knows the opcodes.  The numba kernel below is the CPU
one.  A GPU backend is the same few dozen lines as a compute kernel with the
tape in a buffer, and needs no further knowledge of yapCAD.

Each node kind supplies its emitter as ``spec.backends["tape"]``, beside its
numpy evaluator, so a kind cannot have one without the other drifting out of
sight.  A kind with no emitter makes the whole tree unsupported, and callers
fall back to numpy: see :func:`unsupported_kinds`.

The tape agrees with the numpy backend to rounding -- a few parts in 1e13 on
the trees it has been measured on -- but not bit for bit, because the
arithmetic is ordered differently.  It is therefore used for previews, and
meshes that may be signed and compared still come from numpy.
"""

import math
import os
import threading
from collections import OrderedDict, namedtuple

import numpy as np

from yapcad.sdf.node import SdfError, get_spec
# Imported for their side effect: these modules register the node kinds the
# emitters below attach to.
from yapcad.sdf import ops as _ops            # noqa: F401
from yapcad.sdf import planar as _planar      # noqa: F401
from yapcad.sdf import primitives as _prims   # noqa: F401

__all__ = ["Tape", "TapeUnsupported", "flatten", "unsupported_kinds",
           "numba_available", "run", "compiled", "OPCODES"]

#: The instruction set.  Every backend implements exactly these.
OPCODES = (
    "SPHERE", "BOX", "CYLINDER", "CAPSULE", "TORUS", "CONE", "HALF_SPACE",
    "GYROID", "SCHWARZ_P", "CHAMFERED_BOX", "CHAMFERED_CYLINDER",
    "TRANSFORM", "TWIST",
    "SCALE", "ADD", "SHELL", "COPY",
    "MIN", "MAX", "MAX_NEG", "SMOOTH_MIN", "SMOOTH_MAX", "SMOOTH_SUB",
)
(SPHERE, BOX, CYLINDER, CAPSULE, TORUS, CONE, HALF_SPACE, GYROID, SCHWARZ_P,
 CHAMFERED_BOX, CHAMFERED_CYLINDER,
 TRANSFORM, TWIST, SCALE, ADD, SHELL, COPY,
 MIN, MAX, MAX_NEG, SMOOTH_MIN, SMOOTH_MAX, SMOOTH_SUB) = range(len(OPCODES))

#: Opcodes that read a point: their first operand names three registers.
_READS_POINT = frozenset((
    SPHERE, BOX, CYLINDER, CAPSULE, TORUS, CONE, HALF_SPACE, GYROID,
    SCHWARZ_P, CHAMFERED_BOX, CHAMFERED_CYLINDER, TRANSFORM, TWIST))
#: Opcodes that write a point: their destination names three registers.
_WRITES_POINT = frozenset((TRANSFORM, TWIST))
#: Opcodes with a second operand.
_BINARY = frozenset((MIN, MAX, MAX_NEG, SMOOTH_MIN, SMOOTH_MAX, SMOOTH_SUB))

#: A flattened field.
#:
#: ``ops`` is ``(n, 5)`` int32: opcode, destination register, first operand
#: register, second operand register, offset into ``params``.  A coordinate
#: operand names the first of three consecutive registers.  ``registers`` is
#: how many a point needs; registers 0-2 hold the query point, and
#: ``result`` is the register holding the field value at the end.
Tape = namedtuple("Tape", ("ops", "params", "registers", "result"))


class TapeUnsupported(SdfError):
    """The tree has a node kind with no tape emitter."""

    def __init__(self, kinds):
        self.kinds = tuple(sorted(kinds))
        super().__init__(
            "no tape emitter for node kind"
            f"{'s' if len(self.kinds) > 1 else ''} "
            f"{', '.join(self.kinds)}; evaluate this field with the numpy "
            "backend")


class _Builder:
    """Collects instructions while the tree is walked."""

    def __init__(self):
        self.ops = []
        self.params = []
        self.registers = 3
        self._memo = {}

    def register(self, count=1):
        first = self.registers
        self.registers += count
        return first

    def constants(self, values):
        offset = len(self.params)
        self.params.extend(float(v) for v in values)
        return offset

    def emit(self, op, dst, a=0, b=0, values=()):
        self.ops.append((op, dst, a, b,
                         self.constants(values) if len(values) else 0))
        return dst

    def value(self, op, a=0, b=0, values=()):
        """Emit an instruction into a fresh register and return it."""
        return self.emit(op, self.register(), a, b, values)

    def child(self, node, coords):
        """The register holding ``node``'s value at ``coords``.

        Memoised on the node and the coordinates it is asked at, so a
        subtree shared within the DAG is evaluated once per frame of
        reference, as the numpy backend's cache does.
        """
        key = (id(node), coords)
        hit = self._memo.get(key)
        if hit is None:
            emitter = get_spec(node.kind).backends.get("tape")
            if emitter is None:
                raise TapeUnsupported({node.kind})
            hit = emitter(self, node.p, node.children, coords)
            self._memo[key] = hit
        return hit


# ---------------------------------------------------------------------------
# Emitters.  ``(builder, params, children, coords) -> register``
# ---------------------------------------------------------------------------


def _sphere(b, p, _children, c):
    return b.value(SPHERE, c, values=[p["radius"]])


def _box(b, p, _children, c):
    return b.value(BOX, c, values=[*(0.5 * v for v in p["size"]), 0.0])


def _rounded_box(b, p, _children, c):
    r = p["radius"]
    return b.value(BOX, c, values=[*(0.5 * v - r for v in p["size"]), r])


def _cylinder(b, p, _children, c):
    return b.value(CYLINDER, c,
                   values=[p["radius"], 0.5 * p["height"], 0.0])


def _rounded_cylinder(b, p, _children, c):
    e = p["edge_radius"]
    return b.value(CYLINDER, c,
                   values=[p["radius"] - e, 0.5 * p["height"] - e, e])


def _capsule(b, p, _children, c):
    start, end = p["start"], p["end"]
    along = [end[i] - start[i] for i in range(3)]
    return b.value(CAPSULE, c, values=[
        *start, *along, sum(v * v for v in along), p["radius"]])


def _torus(b, p, _children, c):
    return b.value(TORUS, c,
                   values=[p["major_radius"], p["minor_radius"]])


def _cone(b, p, _children, c):
    return b.value(CONE, c, values=[
        0.5 * p["height"], p["radius_bottom"], p["radius_top"]])


def _half_space(b, p, _children, c):
    return b.value(HALF_SPACE, c, values=[*p["normal"], p["offset"]])


def _gyroid(b, p, _children, c):
    return b.value(GYROID, c,
                   values=[2.0 * math.pi / p["period"], p["thickness"]])


def _schwarz_p(b, p, _children, c):
    return b.value(SCHWARZ_P, c,
                   values=[2.0 * math.pi / p["period"], p["thickness"]])


def _chamfered_box(b, p, _children, c):
    return b.value(CHAMFERED_BOX, c, values=[
        *(0.5 * v for v in p["size"]), p["distance"]])


def _chamfered_cylinder(b, p, _children, c):
    return b.value(CHAMFERED_CYLINDER, c, values=[
        p["radius"], 0.5 * p["height"], p["distance"]])


def _transform(b, p, children, c):
    matrix = np.asarray(p["matrix"], dtype=float)
    inverse = np.linalg.inv(matrix)
    local = b.emit(TRANSFORM, b.register(3), c,
                   values=inverse[:3].ravel())
    value = b.child(children[0], local)
    # The same factor the numpy backend applies: the least singular value
    # of the linear part, which is 1 for a rigid placement.
    s_min = float(np.linalg.svd(matrix[:3, :3], compute_uv=False).min())
    return b.value(SCALE, value, values=[s_min])


def _twist(b, p, children, c):
    fold = p["fold_z"]
    local = b.emit(TWIST, b.register(3), c, values=[
        p["rate"], 0.0 if fold is None else fold,
        0.0 if fold is None else 1.0])
    return b.child(children[0], local)


def _chain(op):
    def emitter(b, _p, children, c):
        values = [b.child(child, c) for child in children]
        out = b.value(COPY, values[0])
        for value in values[1:]:
            b.emit(op, out, out, value)
        return out
    return emitter


def _blend(op):
    def emitter(b, p, children, c):
        first = b.child(children[0], c)
        second = b.child(children[1], c)
        return b.value(op, first, second, values=[p["radius"]])
    return emitter


def _offset(b, p, children, c):
    return b.value(ADD, b.child(children[0], c), values=[-p["distance"]])


def _shell(b, p, children, c):
    return b.value(SHELL, b.child(children[0], c),
                   values=[0.5 * p["thickness"]])


for _kind, _emitter in (
        ("sphere", _sphere), ("box", _box), ("rounded_box", _rounded_box),
        ("cylinder", _cylinder), ("rounded_cylinder", _rounded_cylinder),
        ("capsule", _capsule), ("torus", _torus), ("cone", _cone),
        ("half_space", _half_space), ("gyroid", _gyroid),
        ("schwarz_p", _schwarz_p),
        ("chamfered_box", _chamfered_box),
        ("chamfered_cylinder", _chamfered_cylinder),
        ("transform", _transform), ("twist", _twist),
        ("union", _chain(MIN)), ("compound", _chain(MIN)),
        ("intersect", _chain(MAX)), ("subtract", _chain(MAX_NEG)),
        ("smooth_union", _blend(SMOOTH_MIN)),
        ("smooth_intersect", _blend(SMOOTH_MAX)),
        ("smooth_subtract", _blend(SMOOTH_SUB)),
        ("offset", _offset), ("shell", _shell)):
    get_spec(_kind).backends["tape"] = _emitter


# ---------------------------------------------------------------------------
# Flattening
# ---------------------------------------------------------------------------


def unsupported_kinds(node):
    """The node kinds in ``node``'s tree that have no tape emitter, sorted.

    Empty when the tree can be flattened.
    """
    missing, seen, stack = set(), set(), [node]
    while stack:
        current = stack.pop()
        if id(current) in seen:
            continue
        seen.add(id(current))
        if "tape" not in get_spec(current.kind).backends:
            missing.add(current.kind)
        stack.extend(current.children)
    return sorted(missing)


def flatten(node):
    """Flatten ``node`` to a :class:`Tape`.

    :raises TapeUnsupported: naming every kind in the tree that has no
        emitter, not just the first one met.
    """
    missing = unsupported_kinds(node)
    if missing:
        raise TapeUnsupported(missing)
    builder = _Builder()
    result = builder.child(node, 0)
    ops, registers, result = _reuse_registers(builder.ops, result)
    return Tape(ops=np.asarray(ops, dtype=np.int32).reshape(-1, 5),
                params=np.asarray(builder.params or [0.0], dtype=np.float64),
                registers=int(registers), result=int(result))


def _reuse_registers(ops, result):
    """Renumber registers so one is reused once its value is dead.

    The builder hands every value a register of its own, which for a
    500-node part is over a thousand.  Almost all of them are dead a few
    instructions later: a chain of blends needs the running value and the
    link being added, and nothing else.  Reuse brings a point's working set
    down to a few dozen registers -- a cache line or two on a CPU, and the
    difference between possible and not on a GPU, where every thread
    carries its registers as private memory.

    A value is live from the instruction that first writes it to the last
    that mentions it.  Its register is released only after that last
    instruction has run, never during it, so no instruction reads and
    writes the same register unless the builder asked for that.
    """
    last = {}
    for index, (op, dst, a, b, _q) in enumerate(ops):
        last[dst] = index
        last[a] = index
        if op in _BINARY:
            last[b] = index
    last[result] = len(ops)             # the answer outlives the tape
    assigned = {0: 0}                   # the query point stays where it is
    free = {1: [], 3: []}
    width = {0: 3}
    top = 3
    out = []
    for index, (op, dst, a, b, q) in enumerate(ops):
        if dst not in assigned:
            size = 3 if op in _WRITES_POINT else 1
            width[dst] = size
            if free[size]:
                assigned[dst] = free[size].pop()
            else:
                assigned[dst] = top
                top += size
        out.append((op, assigned[dst], assigned[a],
                    assigned[b] if op in _BINARY else 0, q))
        for virtual in {dst, a, b} if op in _BINARY else {dst, a}:
            if virtual != 0 and last[virtual] == index:
                free[width[virtual]].append(assigned[virtual])
    return out, top, assigned[result]


# ---------------------------------------------------------------------------
# The numba backend
# ---------------------------------------------------------------------------

_kernels = None
_kernel_lock = threading.Lock()

#: Set in a forked child whose parent had already run the parallel kernel.
#: numba's worker threads do not survive a fork, and a parallel kernel
#: launched in such a child waits on them for ever, so the child runs the
#: same interpreter on one thread instead.  Parts are rendered and meshed in
#: worker processes often enough for this to matter.
_parallel_started = False
_serial_only = False


def _after_fork_in_child():
    global _serial_only, _kernel_lock, _tapes_lock
    if _parallel_started:
        _serial_only = True
    # A lock held by another thread at the moment of the fork stays held
    # in the child, where that thread does not exist.
    _kernel_lock = threading.Lock()
    _tapes_lock = threading.Lock()


if hasattr(os, "register_at_fork"):
    os.register_at_fork(after_in_child=_after_fork_in_child)


def numba_available():
    """True if numba can be imported."""
    try:
        import numba  # noqa: F401
    except Exception:       # a broken install is as good as a missing one
        return False
    return True


def _build_kernels():
    """The interpreter, as ``(serial, parallel)`` compiled functions.

    They are two functions rather than one compiled two ways because numba
    caches by function, not by how it was asked to compile it: one body
    under both settings would hand the parallel build to whoever asked for
    the serial one.
    """
    import numba

    root_half = math.sqrt(0.5)
    root_third = math.sqrt(1.0 / 3.0)

    @numba.njit(cache=True)
    def segment(px, py, ax, ay, bx, by):
        ex = bx - ax
        ey = by - ay
        ee = ex * ex + ey * ey
        wx = px - ax
        wy = py - ay
        t = 0.0
        if ee > 0.0:
            t = min(max((wx * ex + wy * ey) / ee, 0.0), 1.0)
        return math.hypot(wx - t * ex, wy - t * ey)

    @numba.njit(cache=True)
    def chamfered_rect(u, v, a, b, d):
        inside = max(max(u - a, v - b), (u + v - (a + b - d)) * root_half)
        if inside <= 0.0:
            return inside
        return min(min(segment(u, v, a, 0.0, a, b - d),
                       segment(u, v, a, b - d, a - d, b)),
                   segment(u, v, a - d, b, 0.0, b))

    @numba.njit(cache=True)
    def octahedron(x, y, z, s):
        m = x + y + z - s
        if 3.0 * x < m:
            qx = x
            qy = y
            qz = z
        elif 3.0 * y < m:
            qx = y
            qy = z
            qz = x
        elif 3.0 * z < m:
            qx = z
            qy = x
            qz = y
        else:
            return m * root_third
        k = min(max(0.5 * (qz - qy + s), 0.0), s)
        return math.sqrt(qx * qx + (qy - s + k) ** 2 + (qz - k) ** 2)

    @numba.njit(cache=True)
    def interpret(ops, params, registers, result, points, out, start, stop):
        R = np.empty(registers)
        for i in range(start, stop):
            R[0] = points[i, 0]
            R[1] = points[i, 1]
            R[2] = points[i, 2]
            for j in range(ops.shape[0]):
                op = ops[j, 0]
                d = ops[j, 1]
                a = ops[j, 2]
                q = ops[j, 4]
                if op == TRANSFORM:
                    x = R[a]
                    y = R[a + 1]
                    z = R[a + 2]
                    R[d] = (params[q] * x + params[q + 1] * y
                            + params[q + 2] * z + params[q + 3])
                    R[d + 1] = (params[q + 4] * x + params[q + 5] * y
                                + params[q + 6] * z + params[q + 7])
                    R[d + 2] = (params[q + 8] * x + params[q + 9] * y
                                + params[q + 10] * z + params[q + 11])
                elif op == SCALE:
                    R[d] = R[a] * params[q]
                elif op == COPY:
                    R[d] = R[a]
                elif op == MIN:
                    R[d] = min(R[a], R[ops[j, 3]])
                elif op == MAX:
                    R[d] = max(R[a], R[ops[j, 3]])
                elif op == MAX_NEG:
                    R[d] = max(R[a], -R[ops[j, 3]])
                elif op == SMOOTH_MIN:
                    va = R[a]
                    vb = R[ops[j, 3]]
                    k = params[q]
                    h = min(max(0.5 + 0.5 * (vb - va) / k, 0.0), 1.0)
                    R[d] = vb + (va - vb) * h - k * h * (1.0 - h)
                elif op == SMOOTH_MAX:
                    va = -R[a]
                    vb = -R[ops[j, 3]]
                    k = params[q]
                    h = min(max(0.5 + 0.5 * (vb - va) / k, 0.0), 1.0)
                    R[d] = -(vb + (va - vb) * h - k * h * (1.0 - h))
                elif op == SMOOTH_SUB:
                    va = -R[a]
                    vb = R[ops[j, 3]]
                    k = params[q]
                    h = min(max(0.5 + 0.5 * (vb - va) / k, 0.0), 1.0)
                    R[d] = -(vb + (va - vb) * h - k * h * (1.0 - h))
                elif op == ADD:
                    R[d] = R[a] + params[q]
                elif op == SHELL:
                    R[d] = abs(R[a]) - params[q]
                elif op == CAPSULE:
                    px = R[a] - params[q]
                    py = R[a + 1] - params[q + 1]
                    pz = R[a + 2] - params[q + 2]
                    denom = params[q + 6]
                    if denom != 0.0:
                        t = (px * params[q + 3] + py * params[q + 4]
                             + pz * params[q + 5]) / denom
                        t = min(max(t, 0.0), 1.0)
                        px -= t * params[q + 3]
                        py -= t * params[q + 4]
                        pz -= t * params[q + 5]
                    R[d] = (math.sqrt(px * px + py * py + pz * pz)
                            - params[q + 7])
                elif op == SPHERE:
                    x = R[a]
                    y = R[a + 1]
                    z = R[a + 2]
                    R[d] = math.sqrt(x * x + y * y + z * z) - params[q]
                elif op == BOX:
                    qx = abs(R[a]) - params[q]
                    qy = abs(R[a + 1]) - params[q + 1]
                    qz = abs(R[a + 2]) - params[q + 2]
                    ox = max(qx, 0.0)
                    oy = max(qy, 0.0)
                    oz = max(qz, 0.0)
                    R[d] = (math.sqrt(ox * ox + oy * oy + oz * oz)
                            + min(max(qx, max(qy, qz)), 0.0)
                            - params[q + 3])
                elif op == CYLINDER:
                    dr = math.hypot(R[a], R[a + 1]) - params[q]
                    da = abs(R[a + 2]) - params[q + 1]
                    R[d] = (math.hypot(max(dr, 0.0), max(da, 0.0))
                            + min(max(dr, da), 0.0) - params[q + 2])
                elif op == TORUS:
                    radial = math.hypot(R[a], R[a + 1]) - params[q]
                    R[d] = math.hypot(radial, R[a + 2]) - params[q + 1]
                elif op == CONE:
                    h = params[q]
                    r1 = params[q + 1]
                    r2 = params[q + 2]
                    qx = math.hypot(R[a], R[a + 1])
                    qy = R[a + 2]
                    cap = r1 if qy < 0.0 else r2
                    ca_x = qx - min(qx, cap)
                    ca_y = abs(qy) - h
                    k2x = r2 - r1
                    k2y = 2.0 * h
                    t = (((r2 - qx) * k2x + (h - qy) * k2y)
                         / (k2x * k2x + k2y * k2y))
                    t = min(max(t, 0.0), 1.0)
                    cb_x = qx - r2 + k2x * t
                    cb_y = qy - h + k2y * t
                    nearest = math.sqrt(min(ca_x * ca_x + ca_y * ca_y,
                                            cb_x * cb_x + cb_y * cb_y))
                    if cb_x < 0.0 and ca_y < 0.0:
                        nearest = -nearest
                    R[d] = nearest
                elif op == HALF_SPACE:
                    R[d] = (R[a] * params[q] + R[a + 1] * params[q + 1]
                            + R[a + 2] * params[q + 2] - params[q + 3])
                elif op == GYROID:
                    w = params[q]
                    x = R[a] * w
                    y = R[a + 1] * w
                    z = R[a + 2] * w
                    R[d] = abs(math.sin(x) * math.cos(y)
                               + math.sin(y) * math.cos(z)
                               + math.sin(z) * math.cos(x)) - params[q + 1]
                elif op == SCHWARZ_P:
                    w = params[q]
                    R[d] = abs(math.cos(R[a] * w) + math.cos(R[a + 1] * w)
                               + math.cos(R[a + 2] * w)) - params[q + 1]
                elif op == CHAMFERED_BOX:
                    x = abs(R[a])
                    y = abs(R[a + 1])
                    z = abs(R[a + 2])
                    hx = params[q]
                    hy = params[q + 1]
                    hz = params[q + 2]
                    bevel = params[q + 3]
                    f = chamfered_rect(x, y, hx, hy, bevel)
                    f = max(f, chamfered_rect(y, z, hy, hz, bevel))
                    f = max(f, chamfered_rect(x, z, hx, hz, bevel))
                    R[d] = max(f, octahedron(x, y, z,
                                             hx + hy + hz - 2.0 * bevel))
                elif op == CHAMFERED_CYLINDER:
                    R[d] = chamfered_rect(
                        math.hypot(R[a], R[a + 1]), abs(R[a + 2]),
                        params[q], params[q + 1], params[q + 2])
                elif op == TWIST:
                    z = R[a + 2]
                    if params[q + 2] != 0.0:
                        angle = -params[q] * (params[q + 1]
                                              - abs(z - params[q + 1]))
                    else:
                        angle = -params[q] * z
                    c = math.cos(angle)
                    s = math.sin(angle)
                    x = R[a]
                    y = R[a + 1]
                    R[d] = c * x - s * y
                    R[d + 1] = s * x + c * y
                    R[d + 2] = z
            out[i] = R[result]

    @numba.njit(parallel=True, cache=True)
    def interpret_parallel(ops, params, registers, result, points, out):
        count = points.shape[0]
        block = 4096
        for chunk in numba.prange((count + block - 1) // block):
            interpret(ops, params, registers, result, points, out,
                      chunk * block, min(count, (chunk + 1) * block))

    return interpret, interpret_parallel


def run(tape, points):
    """Evaluate ``tape`` at ``(N, 3)`` float64 ``points`` with numba.

    The first call in a process loads the compiled interpreter, and the
    first call on a machine compiles it, which takes a few seconds; numba
    caches the result.  A forked child of a process that has already run it
    evaluates on one thread (see ``_serial_only``).  Calls are serialised: numba's default threading
    layer does not allow two Python threads inside a parallel kernel at
    once, and the kernel already uses every core.
    """
    global _kernels, _parallel_started
    pts = np.ascontiguousarray(points, dtype=np.float64)
    out = np.empty(pts.shape[0], dtype=np.float64)
    with _kernel_lock:
        if _kernels is None:
            _kernels = _build_kernels()
        serial, parallel = _kernels
        if _serial_only:
            serial(tape.ops, tape.params, tape.registers, tape.result,
                   pts, out, 0, pts.shape[0])
        else:
            _parallel_started = True
            parallel(tape.ops, tape.params, tape.registers, tape.result,
                     pts, out)
    return out


#: Tapes of recently evaluated trees, by node identity.  A node is immutable
#: and is kept alive by its entry, so the identity cannot be reused while
#: the entry stands.  Hashing the node itself would walk the tree on every
#: call, which is the cost this module exists to avoid.
_TAPES = OrderedDict()
_TAPES_KEPT = 64
_tapes_lock = threading.Lock()


def compiled(node):
    """The :class:`Tape` for ``node``, flattened once and remembered."""
    key = id(node)
    with _tapes_lock:
        hit = _TAPES.get(key)
        if hit is not None and hit[0] is node:
            _TAPES.move_to_end(key)
            return hit[1]
    tape = flatten(node)
    with _tapes_lock:
        _TAPES[key] = (node, tape)
        _TAPES.move_to_end(key)
        while len(_TAPES) > _TAPES_KEPT:
            _TAPES.popitem(last=False)
    return tape
