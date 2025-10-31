"""Command-line helpers for running analysis plans."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional, Sequence

from ..core import PackageManifest
from .base import AnalysisPlan, AnalysisResult, get_backend, load_plan


def _timestamp() -> str:
    return datetime.now(timezone.utc).isoformat()


def _ensure_plan_entry(manifest: PackageManifest, plan: AnalysisPlan, plan_path: Path) -> None:
    validation = manifest.data.setdefault("validation", {})
    plans = validation.setdefault("plans", [])
    if any(entry.get("id") == plan.plan_id for entry in plans):
        return
    try:
        rel_path = str(plan_path.relative_to(manifest.root))
    except ValueError:
        rel_path = str(plan_path)
    entry = {
        "id": plan.plan_id,
        "path": rel_path,
        "kind": plan.kind,
        "backend": plan.backend,
    }
    exec_mode = plan.execution.mode
    if exec_mode:
        entry["execution"] = {"mode": exec_mode}
        if plan.execution.transport:
            entry["execution"]["transport"] = plan.execution.transport
        if plan.execution.host:
            entry["execution"]["host"] = plan.execution.host
    plans.append(entry)


def analyze_package(package_path: Path | str, plan_path: Path | str, *, status: str = "pending") -> Path:
    """Record analysis metadata for ``plan_path`` inside ``package_path``.

    This helper prepares the results directory, writes a ``summary.json``
    placeholder, and updates the manifest ``validation.results`` block.
    """

    manifest = PackageManifest.load(Path(package_path))
    root = manifest.root

    plan_path = Path(plan_path)
    if not plan_path.is_absolute():
        plan_path = (root / plan_path).resolve()

    plan = load_plan(plan_path)
    _ensure_plan_entry(manifest, plan, plan_path)

    results_dir = root / "validation" / "results" / plan.plan_id
    results_dir.mkdir(parents=True, exist_ok=True)

    adapter_cls = get_backend(plan.backend)
    timestamp = _timestamp()

    if adapter_cls is not None:
        adapter = adapter_cls()
        result = adapter.run(manifest, plan, results_dir)
        result.backend = result.backend or plan.backend
        result.timestamp = result.timestamp or timestamp
    else:
        summary_payload = {
            "plan": plan.plan_id,
            "kind": plan.kind,
            "backend": plan.backend,
            "status": status,
            "timestamp": timestamp,
            "notes": "Analysis adapter not available; recorded metadata only.",
        }
        if plan.execution.mode:
            exec_info = {"mode": plan.execution.mode}
            if plan.execution.transport:
                exec_info["transport"] = plan.execution.transport
            if plan.execution.host:
                exec_info["host"] = plan.execution.host
            summary_payload["execution"] = exec_info
        if plan.acceptance:
            summary_payload["acceptance"] = plan.acceptance
        if plan.loads:
            summary_payload["loads"] = plan.loads
        if plan.boundary_conditions:
            summary_payload["boundaryConditions"] = plan.boundary_conditions
        result = AnalysisResult(
            plan_id=plan.plan_id,
            status=status,
            summary=summary_payload,
            backend=plan.backend,
            timestamp=timestamp,
        )

    summary_path = result.summary_path or (results_dir / "summary.json")
    summary_payload = result.summary or {
        "plan": plan.plan_id,
        "status": result.status,
        "backend": result.backend,
        "timestamp": result.timestamp,
    }
    summary_payload.setdefault("plan", plan.plan_id)
    summary_payload.setdefault("status", result.status)
    summary_payload.setdefault("backend", result.backend)
    summary_payload.setdefault("timestamp", result.timestamp)
    if plan.execution.mode and "execution" not in summary_payload:
        exec_info = {"mode": plan.execution.mode}
        if plan.execution.transport:
            exec_info["transport"] = plan.execution.transport
        if plan.execution.host:
            exec_info["host"] = plan.execution.host
        summary_payload["execution"] = exec_info

    with summary_path.open("w", encoding="utf-8") as fp:
        json.dump(summary_payload, fp, indent=2)
        fp.write("\n")

    result.summary_path = summary_path
    result.summary = summary_payload

    validation = manifest.data.setdefault("validation", {})
    results = validation.setdefault("results", [])
    manifest_entry = result.to_manifest_entry(root)
    manifest_entry.setdefault("path", str(summary_path.relative_to(root)))
    existing = next((item for item in results if item.get("plan") == plan.plan_id), None)
    if existing:
        existing.update(manifest_entry)
    else:
        results.append(manifest_entry)

    manifest.save()
    return summary_path


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Run or stage yapCAD analysis plans.")
    parser.add_argument("package", help="Path to the .ycpkg directory")
    parser.add_argument("--plan", required=True, help="Path to the plan YAML (relative to package root)")
    parser.add_argument("--status", default="pending", help="Result status to record (default: pending)")

    args = parser.parse_args(argv)
    summary_path = analyze_package(args.package, args.plan, status=args.status)
    print(f"Analysis summary written to {summary_path}")
    return 0


__all__ = ["analyze_package", "main"]


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
