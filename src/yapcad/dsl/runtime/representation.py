"""Which representation the DSL's solid primitives produce.

By default ``box``, ``cylinder`` and the other primitives produce mesh
solids (with a BREP attached when OCC is available).  With the
representation set to ``"sdf"`` they produce SDF-authored solids instead,
placed exactly as the mesh primitives are, so a design builds unchanged as
fields: its booleans combine fields, its fillets round them exactly, and
fasteners and gears choose their field versions.  That is how a design
like yapRover builds without OpenCASCADE.

Set it per call -- ``compile_and_run(..., representation="sdf")``,
``package_from_dsl(..., representation="sdf")`` -- or for a whole process
with ``YAPCAD_DSL_REPRESENTATION=sdf``, for tools that do not expose the
argument.  ``sdf_cell_mm`` (or ``YAPCAD_DSL_SDF_CELL_MM``) sets the cell
size primitives are meshed at; booleans then mesh no coarser than their
operands.
"""

import contextlib
import contextvars
import math
import os

__all__ = ["REPRESENTATIONS", "current", "is_sdf", "using", "mesh_field",
           "fastener_kwargs"]

REPRESENTATIONS = ("mesh", "sdf")

#: Cells across a primitive's longest extent when no cell size is set.
DEFAULT_CELLS = 64

#: A primitive's thinnest extent is resolved by at least this many cells.
MIN_CELLS_ACROSS = 4

_representation = contextvars.ContextVar("yapcad_dsl_representation",
                                         default=None)
_cell = contextvars.ContextVar("yapcad_dsl_sdf_cell_mm", default=None)


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


@contextlib.contextmanager
def using(representation=None, sdf_cell_mm=None):
    """Set the representation, and the SDF cell size, for a block.

    ``None`` leaves either setting as it is.
    """
    tokens = []
    if representation is not None:
        tokens.append((_representation,
                       _representation.set(_check(representation))))
    if sdf_cell_mm is not None:
        cell = float(sdf_cell_mm)
        if not math.isfinite(cell) or cell <= 0.0:
            raise ValueError("sdf_cell_mm must be positive and finite")
        tokens.append((_cell, _cell.set(cell)))
    try:
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
