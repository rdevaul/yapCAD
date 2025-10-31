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

__all__ = [
    "AnalysisAdapter",
    "AnalysisPlan",
    "AnalysisResult",
    "ExecutionConfig",
    "available_backends",
    "get_backend",
    "load_plan",
    "register_backend",
]
