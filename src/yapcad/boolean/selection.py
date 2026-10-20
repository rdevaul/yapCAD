"""One place that decides which mesh boolean engine runs.

The policy used to live in two places. ``geom3d.solid_boolean`` read the
environment and fell back to manifold3d when it was installed, and the DSL's
batched fast path (``yapcad.dsl.runtime.mesh_boolean``) repeated the same
default inline as ``os.environ.get(..., "manifold")``. Two copies of one
decision meant changing the default took two edits, and the DSL path could
reach a different conclusion than the engine a caller had selected.

Nothing here is a preference about which engine is better; it only reports
what the current configuration resolves to. ``DEFAULT_PREFERENCE`` is the
single knob: the first installed engine in it wins when the caller forces
nothing.
"""

import os

#: Engines to prefer, best first, when nothing is forced. manifold3d leads
#: because it is far more economical with triangles than the native engine --
#: which matters most in boolean chains, where the native engine's extra
#: triangles compound into the next operation -- not because the native
#: engine is less correct. Both return closed, volume-correct solids.
DEFAULT_PREFERENCE = ("manifold", "native")


def is_available(engine):
    """True when ``engine`` can run in this interpreter."""
    if engine == "native":
        return True
    if engine == "manifold":
        from yapcad.boolean import manifold_engine
        return manifold_engine.is_available()
    if engine == "trimesh":
        from yapcad.boolean import trimesh as trimesh_engine
        return trimesh_engine is not None
    if engine == "occ":
        from yapcad.brep import occ_available
        return occ_available()
    return False


def forced_engine():
    """The engine named by ``YAPCAD_BOOLEAN_ENGINE``, or None.

    This one overrides everything, the SDF and OpenCASCADE paths included,
    which is why it is separate from the mesh fallback below. Returned
    whether or not it is installed, so that asking for a missing engine
    reports that rather than quietly running a different one.
    """
    return os.environ.get("YAPCAD_BOOLEAN_ENGINE") or None


def preferred_mesh_engine():
    """Return ``(engine, named)``: the mesh engine to use once the SDF and
    OpenCASCADE paths have declined.

    ``named`` is True when ``YAPCAD_MESH_BOOLEAN_ENGINE`` chose it. Callers
    use that to keep a named engine's errors from being papered over: a
    named engine is expected to fail loudly, where the default may fall
    back to another engine with a warning.
    """
    named = os.environ.get("YAPCAD_MESH_BOOLEAN_ENGINE")
    if named:
        return named, True
    for engine in DEFAULT_PREFERENCE:
        if is_available(engine):
            return engine, False
    return DEFAULT_PREFERENCE[-1], False
