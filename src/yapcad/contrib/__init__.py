"""
Utility modules that are vendored or optional extras.

Currently includes a lightweight port of ``figgear``'s involute gear
generator so that yapCAD can generate gear profiles without requiring
an external dependency at build/test time.
"""

__all__ = ["figgear"]
