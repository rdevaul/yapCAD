import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"

sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(SRC))


def _reason_mentions_occ(reason) -> bool:
    if not reason:
        return False
    r = str(reason).lower()
    return "occ" in r or "pythonocc" in r


def pytest_collection_modifyitems(config, items):
    """Auto-tag OCC-dependent tests with the ``requires_occ`` marker.

    Many test modules gate OCC-only cases behind
    ``skipif(not occ_available(), reason="...pythonocc...")`` at module or
    function scope. This hook detects those markers (by their reason string)
    and applies the ``requires_occ`` marker so the pure-python lane can be
    selected with ``pytest -m "not requires_occ"`` regardless of whether OCC
    happens to be installed in the current environment.

    Tests that already carry ``requires_occ`` explicitly are left as-is.
    """
    for item in items:
        if item.get_closest_marker("requires_occ") is not None:
            continue
        for mark in item.iter_markers(name="skipif"):
            reason = mark.kwargs.get("reason") if mark.kwargs else None
            if reason is None and len(mark.args) >= 2:
                reason = mark.args[1]
            if _reason_mentions_occ(reason):
                item.add_marker(pytest.mark.requires_occ)
                break
