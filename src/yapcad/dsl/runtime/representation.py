"""Which representation the DSL's solid primitives produce.

By default ``box``, ``cylinder`` and the other primitives produce mesh
solids (with a BREP attached when OCC is available).  With the
representation set to ``"sdf"`` they produce SDF-authored solids instead,
placed exactly as the mesh primitives are, so a design builds unchanged as
fields: its booleans combine fields, its fillets round them exactly, and
fasteners and gears choose their field versions.  That is how a design
like yapRover builds without OpenCASCADE.

``sdf_mesh=False`` builds the design without meshing anything.  Meshing is
the expensive half of field authoring -- about seven eighths of the time to
build a part -- and a design that is going to be swept for range of motion,
rendered, or measured for clearance never needs the triangles.  The parts
come back carrying their fields and a plan; :func:`yapcad.sdf.realize_all`
builds the meshes if and when something wants them, and anything that needs
triangles in the meantime says so rather than reading an empty surface list
and believing it.

Set it per call -- ``compile_and_run(..., representation="sdf")``,
``package_from_dsl(..., representation="sdf")`` -- or for a whole process
with ``YAPCAD_DSL_REPRESENTATION=sdf``, for tools that do not expose the
argument.  ``sdf_cell_mm`` (or ``YAPCAD_DSL_SDF_CELL_MM``) sets the cell
size primitives are meshed at; booleans then mesh no coarser than their
operands.

In ``"sdf"`` mode the finished parts a call returns -- the emitted solids
and the parts of an emitted assembly -- are simplified against their
fields (:func:`yapcad.sdf.simplify_finished`).  Intermediate results are
not.  Pass ``sdf_simplify=False``, or set ``YAPCAD_SDF_SIMPLIFY=0``, to
keep the uniform meshes.
"""

import contextlib
import contextvars
import math
import os

__all__ = ["REPRESENTATIONS", "current", "is_sdf", "using", "mesh_field",
           "fastener_kwargs", "simplify_result"]

REPRESENTATIONS = ("mesh", "sdf")

#: Cells across a primitive's longest extent when no cell size is set.
DEFAULT_CELLS = 64

#: A primitive's thinnest extent is resolved by at least this many cells.
MIN_CELLS_ACROSS = 4

_representation = contextvars.ContextVar("yapcad_dsl_representation",
                                         default=None)
_cell = contextvars.ContextVar("yapcad_dsl_sdf_cell_mm", default=None)
_simplify = contextvars.ContextVar("yapcad_dsl_sdf_simplify", default=None)


def _check(value):
    if value not in REPRESENTATIONS:
        raise ValueError(
            f"representation must be one of {REPRESENTATIONS}, got {value!r}"
        )
    return value


def current():
    """The representation primitives produce right now."""
    value = _representation.get()
    if value is None:
        value = os.environ.get("YAPCAD_DSL_REPRESENTATION", "mesh") or "mesh"
    return _check(value)


def is_sdf():
    return current() == "sdf"


def fastener_kwargs():
    """Keyword arguments selecting the fastener builders' representation:
    none in mesh mode, so their calls are exactly what they always were."""
    return {"representation": "sdf"} if is_sdf() else {}


def _cell_size():
    value = _cell.get()
    if value is None:
        env = os.environ.get("YAPCAD_DSL_SDF_CELL_MM")
        value = float(env) if env else None
    return value


def simplify_setting():
    """Whether finished parts are simplified in this block: ``True`` or
    ``False`` when the caller said, ``None`` to follow
    ``YAPCAD_SDF_SIMPLIFY`` (see :func:`yapcad.sdf.simplify_finished`)."""
    return _simplify.get()


def finish_assembly_parts(assembly):
    """Simplify an assembly's SDF parts in place, each distinct part once.

    A part is finished once it is in an assembly that is being turned into
    geometry.  ``assembly_compound`` calls this first, so the compound it
    positions is made of the simplified parts whatever is done to it
    afterwards -- emitted as it is, transformed, or carried under a
    command's ``@meta``.  Simplifying is idempotent, and it leaves alone
    whatever is not SDF-authored and whatever has no mesh yet.
    """
    from yapcad.sdf.simplify import simplify_finished
    geometry = getattr(assembly, "geometry", None)
    if not is_sdf() or not isinstance(geometry, dict):
        return
    done = {}
    for name, solid in list(geometry.items()):
        if isinstance(solid, list) and solid and solid[0] == "solid":
            if id(solid) not in done:
                done[id(solid)] = simplify_finished(solid, simplify_setting())
            geometry[name] = done[id(solid)]


@contextlib.contextmanager
def using(representation=None, sdf_cell_mm=None, sdf_mesh=None,
          sdf_simplify=None):
    """Set the representation, cell size and meshing policy for a block.

    ``None`` leaves any of them as it is.  ``sdf_mesh=False`` defers
    every primitive's and every boolean's meshing; see
    :func:`yapcad.sdf.deferring`.  ``sdf_simplify`` says whether finished
    parts are simplified.
    """
    tokens = []
    if sdf_simplify is not None:
        tokens.append((_simplify, _simplify.set(bool(sdf_simplify))))
    if representation is not None:
        tokens.append((_representation,
                       _representation.set(_check(representation))))
    if sdf_cell_mm is not None:
        cell = float(sdf_cell_mm)
        if not math.isfinite(cell) or cell <= 0.0:
            raise ValueError("sdf_cell_mm must be positive and finite")
        tokens.append((_cell, _cell.set(cell)))
    try:
        with contextlib.ExitStack() as stack:
            if sdf_mesh is not None:
                from yapcad.sdf.convert import deferring
                stack.enter_context(deferring(not sdf_mesh))
            yield
    finally:
        for var, token in reversed(tokens):
            var.reset(token)


def mesh_field(node):
    """Mesh a primitive's field into an SDF-authored solid.

    The cell is the configured size, or the primitive's longest extent over
    :data:`DEFAULT_CELLS`, and never coarser than its thinnest extent over
    :data:`MIN_CELLS_ACROSS`, so a washer or a bearing shield is not lost.
    """
    from yapcad.sdf.convert import to_solid
    lo, hi = node.bounds
    dims = [hi[i] - lo[i] for i in range(3)]
    cell = _cell_size() or max(dims) / DEFAULT_CELLS
    cell = min(cell, min(dims) / MIN_CELLS_ACROSS)
    return to_solid(node, resolution=max(8, math.ceil(max(dims) / cell)))


def simplify_result(result, simplify=None):
    """Simplify the finished SDF parts in an execution result, in place.

    That is the emitted solid or solids, and the parts of an emitted
    assembly, each distinct part once.  Anything not SDF-authored is left
    alone.  ``simplify=None`` follows the block's setting.
    """
    from yapcad.sdf.simplify import simplify_finished
    if simplify is None:
        simplify = simplify_setting()

    emit = getattr(result, "emit_result", None)
    value = getattr(emit, "value", None)
    if value is None:
        return result
    done = {}

    def finish(item):
        key = id(item)
        if key not in done:
            done[key] = simplify_finished(item, simplify)
        return done[key]

    data = value.data
    if isinstance(data, list) and data and data[0] == "solid":
        value.data = finish(data)
    elif isinstance(data, list):
        value.data = [finish(item) if isinstance(item, list) and item
                      and item[0] == "solid" else item for item in data]
    annotations = getattr(value, "annotations", None)
    assembly = annotations.get("assembly") if isinstance(annotations,
                                                         dict) else None
    if assembly is not None and isinstance(getattr(assembly, "geometry",
                                                   None), dict):
        changed = False
        for name, solid in list(assembly.geometry.items()):
            if isinstance(solid, list) and solid and solid[0] == "solid":
                finished = finish(solid)
                changed = changed or finished is not solid
                assembly.geometry[name] = finished
        # ``assembly_compound`` simplifies the parts before it positions
        # them, so there is normally nothing left to do here.  A compound
        # built some other way, from parts that were only now simplified,
        # is positioned again -- or the emitted solid, which is what an STL
        # export writes, keeps the uniform meshes while the assembly holds
        # the simplified ones.  The metadata rides along: it is the emitted
        # solid's (a command's ``@meta``), not the geometry's.
        if (changed and isinstance(data, list) and len(data) > 3
                and data[0] == "solid"
                and data[3] == ["procedure", "assembly_compound"]
                and getattr(assembly, "_solved", False)):
            rebuilt = assembly.compound_geometry()
            if len(data) > 4 and isinstance(data[4], dict) and data[4]:
                from yapcad.metadata import set_solid_metadata
                set_solid_metadata(rebuilt, data[4])
            value.data = rebuilt
    return result
