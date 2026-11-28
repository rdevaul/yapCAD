"""The tape: a field flattened for a compiled backend.

The numpy backend is the reference.  The tape has to agree with it for every
node kind it claims to cover, and everything that reaches for the compiled
backend has to do something sensible when it is not there.
"""

import math
import multiprocessing

import numpy as np
import pytest

from yapcad import sdf
from yapcad.sdf import tape
from yapcad.sdf.evaluate import evaluate, gradient, normal, resolve_backend
from yapcad.sdf.node import _REGISTRY

needs_numba = pytest.mark.skipif(not tape.numba_available(),
                                 reason="numba is not installed")

BALL = sdf.sphere(3.0)
SLAB = sdf.box((6.0, 4.0, 2.0))

#: One tree per covered kind, each with the kind at its root or directly
#: under a placement that moves it off the axes.
SAMPLES = {
    "sphere": BALL,
    "box": SLAB,
    "rounded_box": sdf.rounded_box((6.0, 4.0, 3.0), 0.7),
    "cylinder": sdf.cylinder(2.0, 5.0),
    "rounded_cylinder": sdf.rounded_cylinder(2.5, 5.0, 0.6),
    "capsule": sdf.capsule((-2.0, 1.0, 0.0), (3.0, -1.0, 2.0), 1.1),
    "torus": sdf.torus(4.0, 1.2),
    "cone": sdf.cone(3.0, 1.0, 5.0),
    "half_space": sdf.half_space((0.3, -0.5, 0.8), 0.4),
    "gyroid": sdf.gyroid(4.0, 0.4),
    "schwarz_p": sdf.schwarz_p(5.0, 0.6),
    "chamfered_box": sdf.chamfered_box((6.0, 4.0, 3.0), 0.7),
    "chamfered_cylinder": sdf.chamfered_cylinder(2.5, 5.0, 0.6),
    "transform": sdf.rotate(sdf.translate(SLAB, (1.0, -2.0, 0.5)),
                            (1.0, 2.0, 3.0), 37.0),
    "twist": sdf.twist(SLAB, 0.2),
    "union": sdf.union(BALL, sdf.translate(SLAB, (3.0, 0.0, 0.0)),
                       sdf.translate(sdf.torus(4.0, 0.8), (0.0, 0.0, 2.0))),
    "compound": sdf.compound(BALL, sdf.translate(SLAB, (9.0, 0.0, 0.0))),
    "intersect": sdf.intersect(BALL, SLAB),
    "subtract": sdf.subtract(SLAB, sdf.cylinder(1.0, 9.0),
                             sdf.translate(BALL, (3.0, 2.0, 0.0))),
    "smooth_union": sdf.smooth_union(BALL, sdf.translate(SLAB, (3.0, 0, 0)),
                                     0.8),
    "smooth_intersect": sdf.smooth_intersect(BALL, SLAB, 0.5),
    "smooth_subtract": sdf.smooth_subtract(SLAB, sdf.cylinder(1.0, 9.0),
                                           0.4),
    "offset": sdf.offset(SLAB, 0.6),
    "shell": sdf.shell(BALL, 0.5),
}

EXTRA = {
    "a degenerate capsule": sdf.capsule((1.0, 1.0, 1.0), (1.0, 1.0, 1.0),
                                        2.0),
    "a scaled part": sdf.scale(sdf.rounded_box((4.0, 3.0, 2.0), 0.5), 2.5),
    "a folded twist": sdf.twist(SLAB, 0.3, fold_z=0.5),
    "a cylinder (cone with equal radii)": sdf.cone(2.0, 2.0, 4.0),
    "a pointed cone": sdf.cone(3.0, 0.0, 6.0),
    "a shared subtree": sdf.union(
        sdf.translate(SAMPLES["subtract"], (-8.0, 0.0, 0.0)),
        sdf.rotate(SAMPLES["subtract"], (0.0, 0.0, 1.0), 90.0)),
}


def cloud(count=20000, reach=9.0, seed=3):
    return np.random.default_rng(seed).uniform(-reach, reach, (count, 3))


def test_every_kind_is_either_covered_or_known_not_to_be():
    """Adding a node kind should make someone decide about the tape."""
    covered = {kind for kind, spec in _REGISTRY.items()
               if "tape" in spec.backends}
    assert covered == set(SAMPLES)
    assert set(_REGISTRY) - covered == {
        "polygon", "extrude", "apex_extrude", "thread", "hex_bolt",
        "hex_nut", "spur_gear", "straight_bevel_gear", "repeat",
        "repeat_polar", "revolve"}


@needs_numba
@pytest.mark.parametrize("name", sorted(SAMPLES) + sorted(EXTRA))
def test_the_tape_agrees_with_numpy(name):
    node = {**SAMPLES, **EXTRA}[name]
    points = cloud()
    expected = evaluate(node, points)
    got = evaluate(node, points, backend="numba")
    assert got == pytest.approx(expected, rel=1e-11, abs=1e-11)


@needs_numba
def test_points_on_the_axes_and_at_the_origin():
    """Where the folded primitives have their kinks."""
    axes = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, -1],
                     [3, 0, 0], [0, 0, 2.5], [0, 0, -2.5]], dtype=float)
    for name, node in SAMPLES.items():
        assert evaluate(node, axes, backend="numba") == pytest.approx(
            evaluate(node, axes), rel=1e-11, abs=1e-11), name


@needs_numba
def test_a_deep_blend_chain():
    """The shape of the slow preview: hundreds of nodes and no hard union
    for the numpy backend's culling to bite on."""
    body = None
    for i in range(120):
        a = 0.25 * i
        link = sdf.capsule(
            (6 * math.cos(a), 6 * math.sin(a), 0.4 * i),
            (6 * math.cos(a + 0.25), 6 * math.sin(a + 0.25), 0.4 * i + 0.4),
            0.8)
        body = link if body is None else sdf.smooth_union(body, link, 0.5)
    points = cloud(5000, reach=8.0) + (0.0, 0.0, 24.0)
    assert evaluate(body, points, backend="numba") == pytest.approx(
        evaluate(body, points), rel=1e-10, abs=1e-10)
    assert len(tape.flatten(body).ops) == 2 * 120 - 1


@needs_numba
def test_gradient_and_normal_go_through_the_tape():
    node, points = SAMPLES["subtract"], cloud(2000, reach=4.0)
    assert gradient(node, points, backend="numba") == pytest.approx(
        gradient(node, points), abs=1e-5)
    assert normal(node, points, backend="numba") == pytest.approx(
        normal(node, points), abs=1e-5)


@needs_numba
def test_a_single_point_is_a_float():
    value = evaluate(BALL, (5.0, 0.0, 0.0), backend="numba")
    assert isinstance(value, float) and value == pytest.approx(2.0)


def test_a_shared_subtree_is_flattened_once_per_frame():
    shared = SAMPLES["subtract"]
    once = len(tape.flatten(shared).ops)
    same_frame = sdf.union(shared, shared)
    # union() of one distinct child twice: the child's instructions appear
    # once, plus the copy and the min.
    assert len(tape.flatten(same_frame).ops) == once + 2


def test_flatten_is_deterministic():
    a, b = tape.flatten(SAMPLES["union"]), tape.flatten(SAMPLES["union"])
    assert np.array_equal(a.ops, b.ops) and np.array_equal(a.params, b.params)
    assert (a.registers, a.result) == (b.registers, b.result)


def test_every_instruction_is_a_known_opcode_writing_a_valid_register():
    for node in SAMPLES.values():
        flat = tape.flatten(node)
        assert flat.ops.dtype == np.int32 and flat.ops.shape[1] == 5
        assert flat.ops[:, 0].min() >= 0
        assert flat.ops[:, 0].max() < len(tape.OPCODES)
        assert flat.ops[:, 1:4].min() >= 0
        assert flat.ops[:, 1:4].max() < flat.registers
        assert 0 <= flat.result < flat.registers


# -- what happens when the compiled backend cannot be had ---------------------

UNSUPPORTED = sdf.union(BALL, sdf.extrude(sdf.polygon(
    [(0, 0), (2, 0), (2, 2), (0, 2)]), 3.0))


def test_unsupported_kinds_are_all_named():
    assert tape.unsupported_kinds(UNSUPPORTED) == [
        "extrude", "polygon"]
    assert tape.unsupported_kinds(SAMPLES["subtract"]) == []
    with pytest.raises(tape.TapeUnsupported) as caught:
        tape.flatten(UNSUPPORTED)
    assert caught.value.kinds == ("extrude", "polygon")
    assert "numpy" in str(caught.value)


def test_auto_falls_back_to_numpy_for_an_unsupported_tree():
    points = cloud(500)
    assert resolve_backend(UNSUPPORTED, "auto") == "numpy"
    assert np.array_equal(evaluate(UNSUPPORTED, points, backend="auto"),
                          evaluate(UNSUPPORTED, points))
    assert np.array_equal(gradient(UNSUPPORTED, points, backend="auto"),
                          gradient(UNSUPPORTED, points))


@needs_numba
def test_numba_by_name_refuses_an_unsupported_tree():
    with pytest.raises(tape.TapeUnsupported):
        evaluate(UNSUPPORTED, cloud(10), backend="numba")


def test_without_numba_auto_is_numpy_and_numba_says_what_to_install(
        monkeypatch):
    monkeypatch.setattr(tape, "numba_available", lambda: False)
    points = cloud(500)
    assert resolve_backend(BALL, "auto") == "numpy"
    assert np.array_equal(evaluate(BALL, points, backend="auto"),
                          evaluate(BALL, points))
    with pytest.raises(sdf.SdfError, match=r"yapcad\[fast\]"):
        evaluate(BALL, points, backend="numba")


@needs_numba
def test_auto_uses_the_tape_when_it_can():
    assert resolve_backend(SAMPLES["subtract"], "auto") == "numba"
    assert resolve_backend(SAMPLES["subtract"], "numpy") == "numpy"


def test_the_default_backend_is_still_numpy():
    """Meshing and everything else that does not ask stays on the
    reference backend: its results are signed and compared."""
    from yapcad.sdf import evaluate as module_function  # noqa: F401
    import sys
    assert sys.modules["yapcad.sdf.evaluate"].DEFAULT_BACKEND == "numpy"
    assert sys.modules["yapcad.sdf.contour"].DEFAULT_BACKEND == "numpy"


def _child(queue):
    queue.put(float(evaluate(SAMPLES["subtract"], cloud(60000),
                             backend="numba").sum()))


@needs_numba
@pytest.mark.skipif("fork" not in multiprocessing.get_all_start_methods(),
                    reason="needs fork")
def test_the_kernel_works_in_a_forked_child():
    expected = float(evaluate(SAMPLES["subtract"], cloud(60000),
                              backend="numba").sum())
    context = multiprocessing.get_context("fork")
    queue = context.Queue()
    child = context.Process(target=_child, args=(queue,))
    child.start()
    try:
        assert queue.get(timeout=120) == pytest.approx(expected, rel=1e-12)
    finally:
        child.join(timeout=10)
        if child.is_alive():
            child.kill()


@needs_numba
def test_concurrent_callers_are_serialised():
    from concurrent.futures import ThreadPoolExecutor
    node, points = SAMPLES["smooth_union"], cloud(30000)
    expected = evaluate(node, points, backend="numba")
    with ThreadPoolExecutor(4) as pool:
        results = list(pool.map(
            lambda _: evaluate(node, points, backend="numba"), range(8)))
    for got in results:
        assert np.array_equal(got, expected)


# -- the GPU backend ----------------------------------------------------------

def _gpu():
    try:
        from yapcad.sdf import tape_wgpu
    except Exception:
        return None
    return tape_wgpu if tape_wgpu.available() else None


needs_gpu = pytest.mark.skipif(_gpu() is None,
                               reason="no wgpu adapter on this machine")


def test_the_gpu_kernel_names_every_opcode():
    """The kernel's opcode numbers are written in from ``tape.OPCODES``, so
    this needs no GPU: every opcode must have a case, and no placeholder may
    survive the substitution."""
    from yapcad.sdf import tape_wgpu
    for entry in (tape_wgpu._EVALUATE, tape_wgpu._MARCH):
        source = tape_wgpu._kernel_source(entry)
        assert "OP_" not in source
        assert "MAX_REGISTERS" not in source and "GROUP)" not in source
        for number in range(len(tape.OPCODES)):
            assert f"case {number} " in source


@needs_gpu
@pytest.mark.parametrize("name", sorted(SAMPLES) + sorted(EXTRA))
def test_the_gpu_agrees_with_numpy_to_single_precision(name):
    node = {**SAMPLES, **EXTRA}[name]
    points = cloud(4000)
    assert evaluate(node, points, backend="wgpu") == pytest.approx(
        evaluate(node, points), abs=1e-4)


@needs_gpu
def test_the_gpu_gradient_widens_its_step():
    """At the default step, single precision differences two values that
    agree to a digit; the normals must still come out right."""
    node, points = SAMPLES["subtract"], cloud(2000, reach=4.0)
    reference = normal(node, points)
    got = normal(node, points, backend="wgpu")
    # A point beside a crease takes its normal from whichever face the
    # wider step lands on; everywhere else the two agree closely.
    agree = np.einsum("ij,ij->i", reference, got) > 0.999
    assert agree.mean() > 0.97


def test_a_tape_too_wide_for_the_gpu_is_refused():
    from yapcad.sdf import tape_wgpu
    wide = tape.flatten(BALL)._replace(registers=tape_wgpu.MAX_REGISTERS + 1)
    with pytest.raises(sdf.SdfError, match="registers"):
        tape_wgpu.run(wide, cloud(4))


def test_register_reuse_keeps_a_long_chain_narrow():
    body = BALL
    for i in range(200):
        body = sdf.smooth_union(
            body, sdf.translate(BALL, (0.5 * i, 0.0, 0.0)), 0.4)
    assert tape.flatten(body).registers < 16


# -- the renderer -------------------------------------------------------------

@needs_numba
def test_a_preview_is_the_same_image_on_the_compiled_backend():
    from yapcad import render
    from yapcad.render.api import Part
    scene = [Part(SAMPLES["subtract"], name="plate"),
             Part(sdf.translate(SAMPLES["chamfered_box"], (0.0, 0.0, 6.0)),
                  name="block")]
    reference = render.image(scene, pixels=(160, 120), backend="numpy")
    assert np.array_equal(
        render.image(scene, pixels=(160, 120), backend="numba"), reference)
    assert np.array_equal(
        render.image(scene, pixels=(160, 120)), reference)      # auto


def test_a_render_says_which_backend_drew_each_part():
    from yapcad.render.camera import Camera
    from yapcad.render.trace import render_part
    camera = Camera.fit(UNSUPPORTED.bounds, view="iso", pixels=(48, 36))
    assert render_part(UNSUPPORTED, camera).stats["backend"] == "numpy"
    expected = "numba" if tape.numba_available() else "numpy"
    assert render_part(BALL, camera).stats["backend"] == expected


# -- the march on the GPU -----------------------------------------------------

def _rays(node, pixels=(96, 72)):
    from yapcad.render.camera import Camera
    camera = Camera.fit(node.bounds, view="iso", pixels=pixels)
    origins, dirs = camera.rays()
    lo, hi = (np.asarray(v, dtype=float) for v in node.bounds)
    return camera, origins, dirs, (lo, hi)


@needs_gpu
@pytest.mark.parametrize("name", ["subtract", "smooth_union", "twist",
                                  "chamfered_box", "torus"])
def test_the_gpu_march_lands_where_the_reference_march_does(name):
    from yapcad.render.trace import march, march_gpu
    node = SAMPLES[name]
    camera, origins, dirs, box = _rays(node)
    tolerance = 0.45 * camera.mm_per_px
    expected, _ = march(node, origins, dirs, box, camera.far, tolerance,
                        backend="numpy")
    depth, normals, stats = march_gpu(node, origins, dirs, box, camera.far,
                                      tolerance)
    hit, reference_hit = np.isfinite(depth), np.isfinite(expected)
    # A ray grazing a silhouette can fall either side in single precision.
    assert (hit != reference_hit).mean() < 0.01
    both = hit & reference_hit
    # Both stopped within the tolerance of the surface, so they are within
    # a few tolerances of each other along the ray.
    assert np.abs(depth[both] - expected[both]).max() < 20 * tolerance
    assert np.median(np.abs(depth[both] - expected[both])) < 1e-3
    lengths = np.linalg.norm(normals[hit], axis=1)
    assert lengths == pytest.approx(1.0, abs=1e-4)
    assert not normals[~hit].any()
    assert stats["backend"] == "wgpu" and stats["hits"] == int(hit.sum())
    assert stats["evaluations"] >= stats["rays"] - stats["clipped"]


@needs_gpu
def test_the_gpu_march_reports_rays_that_run_out_of_steps():
    from yapcad.render.trace import march_gpu
    node = SAMPLES["subtract"]
    camera, origins, dirs, box = _rays(node)
    _, _, stats = march_gpu(node, origins, dirs, box, camera.far,
                            0.45 * camera.mm_per_px, max_steps=1)
    assert stats["steps"] == 1
    assert stats["unconverged"] > 0


@needs_gpu
def test_a_gpu_preview_is_the_reference_preview_but_for_its_edges():
    from yapcad import render
    from yapcad.render.api import Part
    scene = [Part(SAMPLES["subtract"], name="plate"),
             Part(sdf.translate(SAMPLES["chamfered_box"], (0.0, 0.0, 6.0)),
                  name="block")]
    reference = render.image(scene, pixels=(240, 180), backend="numpy")
    drawn = render.image(scene, pixels=(240, 180), backend="wgpu")
    differs = (reference != drawn).any(axis=2)
    assert differs.mean() < 0.02
