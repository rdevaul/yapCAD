"""Signed distance function support for yapCAD.

This package implements Phase 1 of ``docs/SDF-DESIGN.md``: the node DAG, the
numpy evaluation backend, Lipschitz and exactness tracking, the analytic
primitives and hard CSG.

An SDF in yapCAD is a *peer* representation, not a new root.  A solid
authored as a field is SDF-authoritative and regenerates mesh and BREP
downward; a solid imported from STEP stays BREP-authoritative.  Nothing here
promotes existing geometry.

A quick tour::

    from yapcad import sdf

    part = sdf.subtract(
        sdf.rounded_box((40, 40, 10), 3),
        sdf.translate(sdf.cylinder(6, 20), (10, 0, 0)),
    )

    part.bounds        # ((-20, -20, -5), (20, 20, 5))
    part.exact         # False -- hard CSG bounds distance, does not equal it
    part.lipschitz     # 1.0  -- so |f| is itself a safe ray-marching step
    part.csg_exact     # True -- replayable through OCC for an exact BREP

    sdf.evaluate(part, [[0, 0, 0], [30, 0, 0]])

The tree, not any sampled form, is the source of truth.  It serialises to a
compact node table with :func:`tree_to_json` and can ride in a solid's
``construction`` slot via :func:`to_construction`.

Booleans between two SDF-authored solids -- through
:func:`yapcad.geom3d.solid_boolean`, and so through the DSL -- combine their
fields and stay SDF-authoritative (:mod:`yapcad.sdf.booleans`).

Not yet implemented, and tracked to later phases of the design document:
``sampled`` grid nodes and the ``twist``/``bend``/``repeat`` domain
operators, promotion of BREP and mesh solids to fields for mixed-authority
booleans (the rest of Phase 4), and the shader backend (Phase 5).
"""

from yapcad.sdf.node import (
    EMPTY_BOUNDS,
    INF,
    TREE_FORMAT,
    UNBOUNDED,
    FieldProps,
    Node,
    NodeSpec,
    SdfError,
    analyze,
    bounds_expand,
    bounds_intersect,
    bounds_is_empty,
    bounds_union,
    digest,
    from_construction,
    is_sdf_construction,
    make_node,
    meshing_from_construction,
    register,
    registered_kinds,
    to_construction,
    tree_from_json,
    tree_to_json,
    walk,
)
from yapcad.sdf.primitives import (
    box,
    capsule,
    cone,
    cylinder,
    gyroid,
    half_space,
    rounded_box,
    schwarz_p,
    sphere,
    torus,
)
from yapcad.sdf.ops import (
    intersect,
    offset,
    rotate,
    scale,
    shell,
    smooth_intersect,
    smooth_subtract,
    smooth_union,
    subtract,
    transform,
    translate,
    union,
)
from yapcad.sdf.contour import (
    Mesh,
    dual_contour,
    is_manifold,
    manifold_defects,
)
from yapcad.sdf.convert import NonManifoldMeshError, mesh_to_surface, to_solid
from yapcad.sdf.booleans import combine_all, is_sdf_solid
from yapcad.sdf.planar import apex_extrude, extrude, polygon
from yapcad.sdf.gears import straight_bevel_gear
from yapcad.sdf.occ import (
    csg_blockers,
    describe_blockers,
    to_brep,
)
from yapcad.sdf.evaluate import (
    evaluate,
    gradient,
    normal,
    safe_step,
    sample_grid,
)

__all__ = [
    # node DAG
    "Node",
    "NodeSpec",
    "FieldProps",
    "SdfError",
    "TREE_FORMAT",
    "INF",
    "UNBOUNDED",
    "EMPTY_BOUNDS",
    "make_node",
    "register",
    "registered_kinds",
    "analyze",
    "digest",
    "walk",
    "bounds_union",
    "bounds_intersect",
    "bounds_expand",
    "bounds_is_empty",
    # serialisation
    "tree_to_json",
    "tree_from_json",
    "to_construction",
    "from_construction",
    "is_sdf_construction",
    "meshing_from_construction",
    # primitives
    "sphere",
    "box",
    "rounded_box",
    "cylinder",
    "capsule",
    "torus",
    "cone",
    "half_space",
    "gyroid",
    "schwarz_p",
    # planar profiles
    "polygon",
    "extrude",
    "apex_extrude",
    # parts
    "straight_bevel_gear",
    # operators
    "union",
    "intersect",
    "subtract",
    "smooth_union",
    "smooth_intersect",
    "smooth_subtract",
    "offset",
    "shell",
    "transform",
    "translate",
    "rotate",
    "scale",
    # evaluation
    "evaluate",
    "gradient",
    "normal",
    "safe_step",
    "sample_grid",
    # meshing
    "Mesh",
    "dual_contour",
    "is_manifold",
    "manifold_defects",
    "mesh_to_surface",
    "to_solid",
    "NonManifoldMeshError",
    # booleans between SDF solids
    "combine_all",
    "is_sdf_solid",
    # exact CSG replay
    "csg_blockers",
    "describe_blockers",
    "to_brep",
]
