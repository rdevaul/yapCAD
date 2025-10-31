"""Shared data structures for analysis/validation plans."""

from __future__ import annotations

import abc
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Type


@dataclass
class ExecutionConfig:
    """Execution context for an analysis plan."""

    mode: str = "local"
    command: Optional[str] = None
    transport: Optional[str] = None
    host: Optional[str] = None
    workdir: Optional[str] = None
    env: Dict[str, str] = field(default_factory=dict)
    options: Dict[str, Any] = field(default_factory=dict)
    license: Dict[str, Any] = field(default_factory=dict)

    @property
    def is_remote(self) -> bool:
        return self.mode.lower() in {"remote", "batch"}


@dataclass
class AnalysisPlan:
    """Representation of a validation/analysis plan loaded from YAML."""

    plan_id: str
    kind: str
    backend: str
    name: Optional[str] = None
    description: Optional[str] = None
    geometry: Dict[str, Any] = field(default_factory=dict)
    materials: Dict[str, Any] = field(default_factory=dict)
    loads: List[Dict[str, Any]] = field(default_factory=list)
    boundary_conditions: List[Dict[str, Any]] = field(default_factory=list)
    acceptance: Dict[str, Any] = field(default_factory=dict)
    backend_options: Dict[str, Any] = field(default_factory=dict)
    execution: ExecutionConfig = field(default_factory=ExecutionConfig)
    attachments: List[Dict[str, Any]] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)
    raw: Dict[str, Any] = field(default_factory=dict)

    @property
    def normalized_backend(self) -> str:
        return self.backend.lower()


@dataclass
class AnalysisResult:
    """Container for results emitted by analysis adapters."""

    plan_id: str
    status: str
    metrics: Dict[str, Any] = field(default_factory=dict)
    summary: Dict[str, Any] = field(default_factory=dict)
    artifacts: List[Dict[str, Any]] = field(default_factory=list)
    summary_path: Optional[Path] = None
    backend: Optional[str] = None
    timestamp: Optional[str] = None
    notes: Optional[str] = None

    def to_manifest_entry(self, package_root: Path) -> Dict[str, Any]:
        entry: Dict[str, Any] = {
            "plan": self.plan_id,
            "status": self.status,
        }
        if self.backend:
            entry["backend"] = self.backend
        if self.timestamp:
            entry["timestamp"] = self.timestamp
        if self.summary_path is not None:
            try:
                entry["path"] = str(self.summary_path.relative_to(package_root))
            except ValueError:
                entry["path"] = str(self.summary_path)
        if self.metrics:
            entry["metrics"] = self.metrics
        if self.artifacts:
            entry["artifacts"] = self.artifacts
        if self.notes:
            entry["notes"] = self.notes
        if self.summary:
            entry["summary"] = self.summary
        return entry


class AnalysisAdapter(abc.ABC):
    """Base class for solver adapters."""

    name: str = "analysis-adapter"

    @abc.abstractmethod
    def run(self, manifest, plan: AnalysisPlan, workspace: Path, **kwargs: Any) -> AnalysisResult:
        """Execute the plan and return an ``AnalysisResult``."""


_BACKENDS: Dict[str, Type[AnalysisAdapter]] = {}


def register_backend(name: str, adapter_cls: Type[AnalysisAdapter]) -> None:
    key = name.lower()
    if not issubclass(adapter_cls, AnalysisAdapter):  # pragma: no cover - defensive
        raise TypeError("adapter_cls must inherit AnalysisAdapter")
    _BACKENDS[key] = adapter_cls


def get_backend(name: str) -> Optional[Type[AnalysisAdapter]]:
    return _BACKENDS.get(name.lower())


def available_backends() -> Sequence[str]:
    return tuple(sorted(_BACKENDS.keys()))


def _parse_execution_config(raw: Dict[str, Any]) -> ExecutionConfig:
    exec_cfg = ExecutionConfig(
        mode=str(raw.get("mode", "local")),
        command=raw.get("command"),
        transport=raw.get("transport"),
        host=raw.get("host"),
        workdir=raw.get("workdir"),
        env=dict(raw.get("env", {})),
        options=dict(raw.get("options", {})),
        license=dict(raw.get("license", {})),
    )
    return exec_cfg


def load_plan(path: Path | str) -> AnalysisPlan:
    """Load a YAML analysis plan and return the normalised ``AnalysisPlan``."""

    plan_path = Path(path)
    if not plan_path.exists():
        raise FileNotFoundError(f"analysis plan not found: {plan_path}")
    import yaml

    with plan_path.open("r", encoding="utf-8") as fp:
        data = yaml.safe_load(fp) or {}
    if not isinstance(data, dict):
        raise ValueError(f"analysis plan must be a mapping, got {type(data)!r}")

    required = {"id", "kind", "backend"}
    missing = [key for key in required if key not in data]
    if missing:
        raise ValueError(f"analysis plan missing required keys: {', '.join(missing)}")

    execution_raw = data.get("execution", {}) or {}
    execution = _parse_execution_config(execution_raw)

    # normalise keys
    loads = data.get("loads", []) or []
    boundary_conditions = data.get("boundaryConditions")
    if boundary_conditions is None:
        boundary_conditions = data.get("boundary_conditions", [])

    known_keys: Iterable[str] = {
        "id",
        "kind",
        "backend",
        "name",
        "description",
        "geometry",
        "materials",
        "loads",
        "boundaryConditions",
        "boundary_conditions",
        "acceptance",
        "backendOptions",
        "execution",
        "attachments",
    }

    metadata = {k: v for k, v in data.items() if k not in known_keys}

    plan = AnalysisPlan(
        plan_id=str(data["id"]),
        kind=str(data["kind"]),
        backend=str(data["backend"]),
        name=data.get("name"),
        description=data.get("description"),
        geometry=dict(data.get("geometry", {})),
        materials=dict(data.get("materials", {})),
        loads=list(loads),
        boundary_conditions=list(boundary_conditions or []),
        acceptance=dict(data.get("acceptance", {})),
        backend_options=dict(data.get("backendOptions", {})),
        execution=execution,
        attachments=list(data.get("attachments", []) or []),
        metadata=metadata,
        raw=data,
    )
    return plan


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
