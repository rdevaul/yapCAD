"""The tape on a GPU, through WebGPU.

One compute kernel, written once and compiled once, interprets any
:class:`~yapcad.sdf.tape.Tape`: the instructions and their parameters are
buffers, so nothing is compiled per design.  ``wgpu`` reaches Metal on
Apple silicon, Vulkan on NVIDIA, AMD and Intel under Linux, and Direct3D 12
on Windows, which is why it is the GPU backend rather than CUDA.

Two limits follow from WebGPU and shape how this backend may be used.

**Single precision.**  WGSL has no 64-bit float.  A value is good to about
one part in 1e7 of the model's extent -- microns on a part a few hundred
millimetres across.  That is far below what an image can show, and it is
not enough for the mesher, which differences the field over a ten-thousandth
of a cell to get its normals.  So this backend is for previews;
:func:`yapcad.sdf.evaluate.gradient` widens its step when it is used.

**Opt-in.**  The kernel is checked against the numpy backend in the test
suite wherever an adapter exists, and has been run on an Apple M5 Max
through Metal and an NVIDIA GB10 through Vulkan, where it evaluates a
500-node part at 14 to 20 nanoseconds a point.  That is two machines, not
the range of drivers in use, so ``backend="auto"`` does not choose it: ask
for it by name, or set ``YAPCAD_RENDER_BACKEND=wgpu``.

There are two entry points.  :func:`run` evaluates the field at a list of
points.  :func:`march` sphere-traces a list of rays and takes the normal
where each lands, all in one dispatch -- which is what a preview wants,
since going back to the GPU once per marching step costs more than the
marching does.
"""

import threading

import numpy as np

from yapcad.sdf import tape as _tape
from yapcad.sdf.node import SdfError

__all__ = ["available", "adapter_summary", "run", "march", "MAX_REGISTERS"]

#: Registers a GPU thread carries.  They are private memory per thread, so
#: the bound is fixed in the kernel; register reuse in ``tape.flatten``
#: keeps a 500-node part under a dozen.
MAX_REGISTERS = 64

#: Points per dispatch.  WebGPU limits a dispatch to 65535 groups a side.
_GROUP = 64
_MAX_POINTS = _GROUP * 65535

#: The interpreter itself: ``field(point, op_count, result)``.  Both entry
#: points below are this plus a few lines, so an image and a table of
#: values cannot come from two readings of the tape.
_FIELD = """
@group(0) @binding(0) var<storage, read> ops: array<i32>;
@group(0) @binding(1) var<storage, read> params: array<f32>;

const ROOT_HALF: f32 = 0.7071067811865476;
const ROOT_THIRD: f32 = 0.5773502691896258;

fn segment(p: vec2<f32>, a: vec2<f32>, b: vec2<f32>) -> f32 {
    let e = b - a;
    let ee = dot(e, e);
    let w = p - a;
    var t = 0.0;
    if (ee > 0.0) {
        t = clamp(dot(w, e) / ee, 0.0, 1.0);
    }
    return length(w - t * e);
}

fn chamfered_rect(u: f32, v: f32, a: f32, b: f32, d: f32) -> f32 {
    let inside = max(max(u - a, v - b), (u + v - (a + b - d)) * ROOT_HALF);
    if (inside <= 0.0) {
        return inside;
    }
    let p = vec2<f32>(u, v);
    return min(min(segment(p, vec2<f32>(a, 0.0), vec2<f32>(a, b - d)),
                   segment(p, vec2<f32>(a, b - d), vec2<f32>(a - d, b))),
               segment(p, vec2<f32>(a - d, b), vec2<f32>(0.0, b)));
}

fn octahedron(x: f32, y: f32, z: f32, s: f32) -> f32 {
    let m = x + y + z - s;
    var q: vec3<f32>;
    if (3.0 * x < m) {
        q = vec3<f32>(x, y, z);
    } else if (3.0 * y < m) {
        q = vec3<f32>(y, z, x);
    } else if (3.0 * z < m) {
        q = vec3<f32>(z, x, y);
    } else {
        return m * ROOT_THIRD;
    }
    let k = clamp(0.5 * (q.z - q.y + s), 0.0, s);
    return length(vec3<f32>(q.x, q.y - s + k, q.z - k));
}

fn smooth_min(a: f32, b: f32, k: f32) -> f32 {
    let h = clamp(0.5 + 0.5 * (b - a) / k, 0.0, 1.0);
    return b + (a - b) * h - k * h * (1.0 - h);
}

fn field(point: vec3<f32>, op_count: u32, result: u32) -> f32 {
    var R: array<f32, MAX_REGISTERS>;
    R[0] = point.x;
    R[1] = point.y;
    R[2] = point.z;
    for (var j = 0u; j < op_count; j = j + 1u) {
        let at = j * 5u;
        let op = ops[at];
        let d = ops[at + 1u];
        let a = ops[at + 2u];
        let b = ops[at + 3u];
        let q = ops[at + 4u];
        let p = vec3<f32>(R[a], R[a + 1], R[a + 2]);
        switch op {
            case OP_SPHERE {
                R[d] = length(p) - params[q];
            }
            case OP_BOX {
                let e = abs(p) - vec3<f32>(params[q], params[q + 1],
                                           params[q + 2]);
                R[d] = length(max(e, vec3<f32>(0.0)))
                     + min(max(e.x, max(e.y, e.z)), 0.0) - params[q + 3];
            }
            case OP_CYLINDER {
                let dr = length(p.xy) - params[q];
                let da = abs(p.z) - params[q + 1];
                R[d] = length(max(vec2<f32>(dr, da), vec2<f32>(0.0)))
                     + min(max(dr, da), 0.0) - params[q + 2];
            }
            case OP_CAPSULE {
                var w = p - vec3<f32>(params[q], params[q + 1],
                                      params[q + 2]);
                let along = vec3<f32>(params[q + 3], params[q + 4],
                                      params[q + 5]);
                let denom = params[q + 6];
                if (denom != 0.0) {
                    w = w - clamp(dot(w, along) / denom, 0.0, 1.0) * along;
                }
                R[d] = length(w) - params[q + 7];
            }
            case OP_TORUS {
                R[d] = length(vec2<f32>(length(p.xy) - params[q], p.z))
                     - params[q + 1];
            }
            case OP_CONE {
                let h = params[q];
                let r1 = params[q + 1];
                let r2 = params[q + 2];
                let qx = length(p.xy);
                let qy = p.z;
                var cap = r2;
                if (qy < 0.0) {
                    cap = r1;
                }
                let ca = vec2<f32>(qx - min(qx, cap), abs(qy) - h);
                let k2 = vec2<f32>(r2 - r1, 2.0 * h);
                let t = clamp(dot(vec2<f32>(r2 - qx, h - qy), k2)
                              / dot(k2, k2), 0.0, 1.0);
                let cb = vec2<f32>(qx - r2, qy - h) + k2 * t;
                var nearest = sqrt(min(dot(ca, ca), dot(cb, cb)));
                if (cb.x < 0.0 && ca.y < 0.0) {
                    nearest = -nearest;
                }
                R[d] = nearest;
            }
            case OP_HALF_SPACE {
                R[d] = dot(p, vec3<f32>(params[q], params[q + 1],
                                        params[q + 2])) - params[q + 3];
            }
            case OP_GYROID {
                let w = p * params[q];
                R[d] = abs(sin(w.x) * cos(w.y) + sin(w.y) * cos(w.z)
                           + sin(w.z) * cos(w.x)) - params[q + 1];
            }
            case OP_SCHWARZ_P {
                let w = p * params[q];
                R[d] = abs(cos(w.x) + cos(w.y) + cos(w.z)) - params[q + 1];
            }
            case OP_CHAMFERED_BOX {
                let f = abs(p);
                let hx = params[q];
                let hy = params[q + 1];
                let hz = params[q + 2];
                let bevel = params[q + 3];
                var v = chamfered_rect(f.x, f.y, hx, hy, bevel);
                v = max(v, chamfered_rect(f.y, f.z, hy, hz, bevel));
                v = max(v, chamfered_rect(f.x, f.z, hx, hz, bevel));
                R[d] = max(v, octahedron(f.x, f.y, f.z,
                                         hx + hy + hz - 2.0 * bevel));
            }
            case OP_CHAMFERED_CYLINDER {
                R[d] = chamfered_rect(length(p.xy), abs(p.z), params[q],
                                      params[q + 1], params[q + 2]);
            }
            case OP_TRANSFORM {
                R[d] = params[q] * p.x + params[q + 1] * p.y
                     + params[q + 2] * p.z + params[q + 3];
                R[d + 1] = params[q + 4] * p.x + params[q + 5] * p.y
                         + params[q + 6] * p.z + params[q + 7];
                R[d + 2] = params[q + 8] * p.x + params[q + 9] * p.y
                         + params[q + 10] * p.z + params[q + 11];
            }
            case OP_TWIST {
                var angle = -params[q] * p.z;
                if (params[q + 2] != 0.0) {
                    angle = -params[q] * (params[q + 1]
                                          - abs(p.z - params[q + 1]));
                }
                let c = cos(angle);
                let s = sin(angle);
                R[d] = c * p.x - s * p.y;
                R[d + 1] = s * p.x + c * p.y;
                R[d + 2] = p.z;
            }
            case OP_SCALE {
                R[d] = R[a] * params[q];
            }
            case OP_ADD {
                R[d] = R[a] + params[q];
            }
            case OP_SHELL {
                R[d] = abs(R[a]) - params[q];
            }
            case OP_COPY {
                R[d] = R[a];
            }
            case OP_MIN {
                R[d] = min(R[a], R[b]);
            }
            case OP_MAX {
                R[d] = max(R[a], R[b]);
            }
            case OP_MAX_NEG {
                R[d] = max(R[a], -R[b]);
            }
            case OP_SMOOTH_MIN {
                R[d] = smooth_min(R[a], R[b], params[q]);
            }
            case OP_SMOOTH_MAX {
                R[d] = -smooth_min(-R[a], -R[b], params[q]);
            }
            case OP_SMOOTH_SUB {
                R[d] = -smooth_min(-R[a], R[b], params[q]);
            }
            default {
            }
        }
    }
    return R[result];
}
"""

#: Evaluate the field at a list of points.
_EVALUATE = """
struct Info {
    op_count: u32,
    result: u32,
    count: u32,
    pad: u32,
};

@group(0) @binding(2) var<storage, read> points: array<vec4<f32>>;
@group(0) @binding(3) var<storage, read_write> values: array<f32>;
@group(0) @binding(4) var<uniform> info: Info;

@compute @workgroup_size(GROUP)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= info.count) {
        return;
    }
    values[i] = field(points[i].xyz, info.op_count, info.result);
}
"""

#: Sphere-trace a list of rays to the surface and take the normal there,
#: in one dispatch.  The loop is the one in :func:`yapcad.render.trace.march`:
#: stop within ``tolerance`` of the surface, otherwise advance by the
#: Lipschitz-safe step, and give up past the far end or out of steps.
#:
#: A ray is ``start`` (xyz) with the distance it may travel (w), and a unit
#: direction.  The answer is the distance travelled (x) and the normal
#: (yzw); a distance of -1 is a miss and -2 a ray that ran out of steps.
_MARCH = """
struct March {
    op_count: u32,
    result: u32,
    count: u32,
    max_steps: u32,
    tolerance: f32,
    inverse_lipschitz: f32,
    min_step: f32,
    normal_step: f32,
};

@group(0) @binding(2) var<storage, read> starts: array<vec4<f32>>;
@group(0) @binding(3) var<storage, read> directions: array<vec4<f32>>;
@group(0) @binding(4) var<storage, read_write> hits: array<vec4<f32>>;
@group(0) @binding(5) var<storage, read_write> steps: array<u32>;
@group(0) @binding(6) var<uniform> info: March;

@compute @workgroup_size(GROUP)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= info.count) {
        return;
    }
    let start = starts[i].xyz;
    let reach = starts[i].w;
    let direction = directions[i].xyz;
    var travelled = 0.0;
    var outcome = -2.0;
    var used = 0u;
    while (used < info.max_steps) {
        used = used + 1u;
        let value = field(start + travelled * direction,
                          info.op_count, info.result);
        if (abs(value) < info.tolerance) {
            outcome = travelled;
            break;
        }
        travelled = travelled + max(abs(value) * info.inverse_lipschitz,
                                    info.min_step);
        if (travelled >= reach) {
            outcome = -1.0;
            break;
        }
    }
    steps[i] = used;
    var normal = vec3<f32>(0.0);
    if (outcome >= 0.0) {
        let at = start + outcome * direction;
        let h = info.normal_step;
        let gradient = vec3<f32>(
            field(at + vec3<f32>(h, 0.0, 0.0), info.op_count, info.result)
            - field(at - vec3<f32>(h, 0.0, 0.0), info.op_count, info.result),
            field(at + vec3<f32>(0.0, h, 0.0), info.op_count, info.result)
            - field(at - vec3<f32>(0.0, h, 0.0), info.op_count, info.result),
            field(at + vec3<f32>(0.0, 0.0, h), info.op_count, info.result)
            - field(at - vec3<f32>(0.0, 0.0, h), info.op_count, info.result));
        let size = length(gradient);
        if (size > 0.0) {
            normal = gradient / size;
        }
    }
    hits[i] = vec4<f32>(outcome, normal);
}
"""


def _kernel_source(entry=None):
    """The interpreter with the opcode numbers of :mod:`yapcad.sdf.tape`
    written in, so the two cannot disagree about what an opcode means,
    followed by ``entry`` -- :data:`_EVALUATE` or :data:`_MARCH`."""
    source = _FIELD + (_EVALUATE if entry is None else entry)
    source = source.replace("MAX_REGISTERS", str(MAX_REGISTERS))
    source = source.replace("GROUP", str(_GROUP))
    # Longest names first: OP_MAX must not rewrite the head of OP_MAX_NEG.
    for name in sorted(_tape.OPCODES, key=len, reverse=True):
        source = source.replace(f"case OP_{name} ",
                                f"case {_tape.OPCODES.index(name)} ")
    return source


_state = None
_lock = threading.Lock()
_unavailable = None


def _connect():
    """The device and pipeline, made once.  Raises :class:`SdfError` with
    the reason when there is no usable GPU."""
    global _state, _unavailable
    if _state is not None:
        return _state
    if _unavailable is not None:
        raise SdfError(_unavailable)
    try:
        import wgpu
    except Exception:
        _unavailable = ("the wgpu backend needs wgpu, which is not "
                        "installed; pip install \"yapcad[gpu]\"")
        raise SdfError(_unavailable) from None
    try:
        adapter = wgpu.gpu.request_adapter_sync(
            power_preference="high-performance")
        if adapter is None:
            raise RuntimeError("no adapter")
        device = adapter.request_device_sync()
        pipelines = {}
        for name, entry in (("evaluate", _EVALUATE), ("march", _MARCH)):
            module = device.create_shader_module(code=_kernel_source(entry))
            pipelines[name] = device.create_compute_pipeline(
                layout="auto",
                compute={"module": module, "entry_point": "main"})
    except Exception as exc:
        _unavailable = f"the wgpu backend found no usable GPU: {exc}"
        raise SdfError(_unavailable) from None
    _state = {"wgpu": wgpu, "adapter": adapter, "device": device,
              "pipeline": pipelines["evaluate"],
              "march": pipelines["march"], "tapes": {}}
    return _state


def available():
    """True if a GPU adapter and the compiled kernel can be had."""
    try:
        with _lock:
            _connect()
    except SdfError:
        return False
    return True


def adapter_summary():
    """A line naming the adapter in use, e.g. for a benchmark's output."""
    with _lock:
        info = _connect()["adapter"].info
    return (f"{info.get('device', '?')} "
            f"({info.get('adapter_type', '?')}, "
            f"{info.get('backend_type', '?')})")


def _tape_buffers(state, tape):
    """The instruction and parameter buffers for ``tape``, uploaded once."""
    key = id(tape)
    hit = state["tapes"].get(key)
    if hit is not None and hit[0] is tape:
        return hit[1], hit[2]
    wgpu, device = state["wgpu"], state["device"]
    usage = wgpu.BufferUsage.STORAGE | wgpu.BufferUsage.COPY_DST
    ops = device.create_buffer_with_data(
        data=np.ascontiguousarray(tape.ops, dtype=np.int32).tobytes(),
        usage=usage)
    params = device.create_buffer_with_data(
        data=tape.params.astype(np.float32).tobytes(), usage=usage)
    if len(state["tapes"]) >= 64:
        state["tapes"].clear()
    state["tapes"][key] = (tape, ops, params)
    return ops, params


def run(tape, points):
    """Evaluate ``tape`` at ``(N, 3)`` ``points`` on the GPU.

    :returns: float64 values, computed in single precision.
    :raises SdfError: when there is no usable GPU, or the tape needs more
        registers than the kernel carries.
    """
    if tape.registers > MAX_REGISTERS:
        raise SdfError(
            f"the wgpu backend carries {MAX_REGISTERS} registers and this "
            f"field needs {tape.registers}; evaluate it with numba or numpy")
    pts = np.asarray(points, dtype=np.float32).reshape(-1, 3)
    out = np.empty(pts.shape[0], dtype=np.float64)
    with _lock:
        state = _connect()
        for start in range(0, pts.shape[0], _MAX_POINTS):
            piece = pts[start:start + _MAX_POINTS]
            out[start:start + len(piece)] = _dispatch(state, tape, piece)
    return out


def _dispatch(state, tape, pts):
    wgpu, device = state["wgpu"], state["device"]
    count = pts.shape[0]
    if count == 0:
        return np.zeros(0)
    ops, params = _tape_buffers(state, tape)
    padded = np.zeros((count, 4), dtype=np.float32)
    padded[:, :3] = pts
    points = device.create_buffer_with_data(
        data=padded.tobytes(), usage=wgpu.BufferUsage.STORAGE)
    values = device.create_buffer(
        size=4 * count,
        usage=wgpu.BufferUsage.STORAGE | wgpu.BufferUsage.COPY_SRC)
    info = device.create_buffer_with_data(
        data=np.array([len(tape.ops), tape.result, count, 0],
                      dtype=np.uint32).tobytes(),
        usage=wgpu.BufferUsage.UNIFORM)
    pipeline = state["pipeline"]
    group = device.create_bind_group(
        layout=pipeline.get_bind_group_layout(0),
        entries=[{"binding": n, "resource": {"buffer": buffer}}
                 for n, buffer in enumerate(
                     (ops, params, points, values, info))])
    encoder = device.create_command_encoder()
    compute = encoder.begin_compute_pass()
    compute.set_pipeline(pipeline)
    compute.set_bind_group(0, group)
    compute.dispatch_workgroups((count + _GROUP - 1) // _GROUP)
    compute.end()
    device.queue.submit([encoder.finish()])
    result = np.frombuffer(device.queue.read_buffer(values),
                           dtype=np.float32).astype(np.float64)
    for buffer in (points, values, info):
        buffer.destroy()
    return result


def march(tape, starts, directions, reach, tolerance, lipschitz, max_steps,
          normal_step):
    """Sphere-trace rays through ``tape`` on the GPU, in one dispatch.

    :param starts: ``(N, 3)`` points to march from -- where each ray enters
        the part's bounds, not the eye, so single precision is spent on the
        part and not on the distance to the camera.
    :param directions: ``(N, 3)`` unit vectors.
    :param reach: ``(N,)`` how far each ray may travel before it has missed.
    :param tolerance: a ray within this of the surface has hit it.
    :param lipschitz: the field's Lipschitz constant; a step is
        ``|f| / lipschitz``, and never less than half the tolerance.
    :param normal_step: half-width of the central difference for normals.
    :returns: ``(travelled, normals, steps, unconverged)``: the distance
        along each ray to its hit, ``inf`` for a miss; unit normals, zero
        where there is no hit; the field evaluations each ray's march
        used; and a mask of the rays that ran out of steps.
    """
    if tape.registers > MAX_REGISTERS:
        raise SdfError(
            f"the wgpu backend carries {MAX_REGISTERS} registers and this "
            f"field needs {tape.registers}; evaluate it with numba or numpy")
    starts = np.asarray(starts, dtype=np.float32).reshape(-1, 3)
    count = starts.shape[0]
    travelled = np.full(count, np.inf)
    normals = np.zeros((count, 3))
    used = np.zeros(count, dtype=np.int64)
    unconverged = np.zeros(count, dtype=bool)
    if count == 0:
        return travelled, normals, used, unconverged
    first = np.empty((count, 4), dtype=np.float32)
    first[:, :3] = starts
    first[:, 3] = reach
    second = np.zeros((count, 4), dtype=np.float32)
    second[:, :3] = directions
    settings = (int(max_steps), float(tolerance), 1.0 / float(lipschitz),
                0.5 * float(tolerance), float(normal_step))
    with _lock:
        state = _connect()
        for at in range(0, count, MARCH_RAYS):
            piece = slice(at, min(count, at + MARCH_RAYS))
            hits, steps = _dispatch_march(state, tape, first[piece],
                                          second[piece], settings)
            distance = hits[:, 0].astype(np.float64)
            hit = distance >= 0.0
            travelled[piece] = np.where(hit, distance, np.inf)
            normals[piece] = np.where(hit[:, None], hits[:, 1:4], 0.0)
            used[piece] = steps
            unconverged[piece] = distance == -2.0
    return travelled, normals, used, unconverged


#: Rays per dispatch.  A dispatch that runs for seconds can trip a
#: display driver's watchdog, and a ray can cost ``max_steps`` walks of the
#: tape, so a frame goes over in a few pieces rather than one.
MARCH_RAYS = 1 << 17


def _dispatch_march(state, tape, starts, directions, settings):
    wgpu, device = state["wgpu"], state["device"]
    count = starts.shape[0]
    max_steps, tolerance, inverse_lipschitz, min_step, normal_step = settings
    ops, params = _tape_buffers(state, tape)
    read = wgpu.BufferUsage.STORAGE
    written = wgpu.BufferUsage.STORAGE | wgpu.BufferUsage.COPY_SRC
    start_buffer = device.create_buffer_with_data(
        data=starts.tobytes(), usage=read)
    direction_buffer = device.create_buffer_with_data(
        data=directions.tobytes(), usage=read)
    hit_buffer = device.create_buffer(size=16 * count, usage=written)
    step_buffer = device.create_buffer(size=4 * count, usage=written)
    info = device.create_buffer_with_data(
        data=(np.array([len(tape.ops), tape.result, count, max_steps],
                       dtype=np.uint32).tobytes()
              + np.array([tolerance, inverse_lipschitz, min_step,
                          normal_step], dtype=np.float32).tobytes()),
        usage=wgpu.BufferUsage.UNIFORM)
    pipeline = state["march"]
    buffers = (ops, params, start_buffer, direction_buffer, hit_buffer,
               step_buffer, info)
    group = device.create_bind_group(
        layout=pipeline.get_bind_group_layout(0),
        entries=[{"binding": n, "resource": {"buffer": buffer}}
                 for n, buffer in enumerate(buffers)])
    encoder = device.create_command_encoder()
    compute = encoder.begin_compute_pass()
    compute.set_pipeline(pipeline)
    compute.set_bind_group(0, group)
    compute.dispatch_workgroups((count + _GROUP - 1) // _GROUP)
    compute.end()
    device.queue.submit([encoder.finish()])
    hits = np.frombuffer(device.queue.read_buffer(hit_buffer),
                         dtype=np.float32).reshape(-1, 4).copy()
    steps = np.frombuffer(device.queue.read_buffer(step_buffer),
                          dtype=np.uint32).astype(np.int64)
    for buffer in buffers[2:]:
        buffer.destroy()
    return hits, steps
