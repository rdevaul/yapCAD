"""Analysis helper exports."""

from .base import (
    AnalysisAdapter,
    AnalysisPlan,
    AnalysisResult,
    ExecutionConfig,
    available_backends,
    get_backend,
    load_plan,
    register_backend,
)
from .calculix import CalculixAdapter  # noqa: F401 - ensures backend registration

# Import optional backends (gracefully handle missing dependencies)
try:
    from .fenics import FenicsxAdapter  # noqa: F401 - ensures backend registration
except ImportError:
    pass  # FEniCSx not installed

try:
    from .gmsh_mesher import (
        GmshMesher,
        MeshHints,
        PhysicalGroup,
        gmsh_available,
        mesh_solid,
    )
except ImportError:
    GmshMesher = None  # type: ignore
    MeshHints = None  # type: ignore
    PhysicalGroup = None  # type: ignore
    gmsh_available = lambda: False  # type: ignore
    mesh_solid = None  # type: ignore

# Face naming utilities (always available)
from .face_naming import (
    FaceInfo,
    FaceSelector,
    ByNormalSelector,
    ByAreaSelector,
    ByPositionSelector,
    CombinedSelector,
    FaceNamer,
    top_faces,
    bottom_faces,
    front_faces,
    back_faces,
    left_faces,
    right_faces,
    largest_face,
    smallest_face,
    faces_at_z_min,
    faces_at_z_max,
)

__all__ = [
    "AnalysisAdapter",
    "AnalysisPlan",
    "AnalysisResult",
    "ExecutionConfig",
    "available_backends",
    "get_backend",
    "load_plan",
    "register_backend",
    # Gmsh meshing
    "GmshMesher",
    "MeshHints",
    "PhysicalGroup",
    "gmsh_available",
    "mesh_solid",
    # Face naming
    "FaceInfo",
    "FaceSelector",
    "ByNormalSelector",
    "ByAreaSelector",
    "ByPositionSelector",
    "CombinedSelector",
    "FaceNamer",
    "top_faces",
    "bottom_faces",
    "front_faces",
    "back_faces",
    "left_faces",
    "right_faces",
    "largest_face",
    "smallest_face",
    "faces_at_z_min",
    "faces_at_z_max",
]
