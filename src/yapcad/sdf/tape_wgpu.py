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

**Experimental.**  The kernel is checked against the numpy backend on a
software adapter, which proves the arithmetic and says nothing about speed
on real hardware or about any one driver.  It is therefore never chosen by
``backend="auto"``: ask for it by name, or set
``YAPCAD_RENDER_BACKEND=wgpu``.
"""

import threading

import numpy as np

from yapcad.sdf import tape as _tape
from yapcad.sdf.node import SdfError

__all__ = ["available", "adapter_summary", "run", "MAX_REGISTERS"]

#: Registers a GPU thread carries.  They are private memory per thread, so
#: the bound is fixed in the kernel; register reuse in ``tape.flatten``
#: keeps a 500-node part under a dozen.
MAX_REGISTERS = 64

#: Points per dispatch.  WebGPU limits a dispatch to 65535 groups a side.
_GROUP = 64
_MAX_POINTS = _GROUP * 65535

_KERNEL = """
struct Info {
    op_count: u32,
    result: u32,
    count: u32,
    pad: u32,
};

@group(0) @binding(0) var<storage, read> ops: array<i32>;
@group(0) @binding(1) var<storage, read> params: array<f32>;
@group(0) @binding(2) var<storage, read> points: array<vec4<f32>>;
@group(0) @binding(3) var<storage, read_write> values: array<f32>;
@group(0) @binding(4) var<uniform> info: Info;

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

@compute @workgroup_size(GROUP)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= info.count) {
        return;
    }
    var R: array<f32, MAX_REGISTERS>;
    R[0] = points[i].x;
    R[1] = points[i].y;
    R[2] = points[i].z;
    for (var j = 0u; j < info.op_count; j = j + 1u) {
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
    values[i] = R[info.result];
}
"""


def _kernel_source():
    """The kernel with the opcode numbers of :mod:`yapcad.sdf.tape` written
    in, so the two cannot disagree about what an opcode means."""
    source = _KERNEL.replace("MAX_REGISTERS", str(MAX_REGISTERS))
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
        module = device.create_shader_module(code=_kernel_source())
        pipeline = device.create_compute_pipeline(
            layout="auto", compute={"module": module, "entry_point": "main"})
    except Exception as exc:
        _unavailable = f"the wgpu backend found no usable GPU: {exc}"
        raise SdfError(_unavailable) from None
    _state = {"wgpu": wgpu, "adapter": adapter, "device": device,
              "pipeline": pipeline, "tapes": {}}
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
