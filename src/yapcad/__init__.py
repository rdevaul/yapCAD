# -*- coding: utf-8 -*-
try:  # Python >= 3.8
    from importlib.metadata import PackageNotFoundError, version
except ModuleNotFoundError:  # pragma: no cover - for Python < 3.8
    from importlib_metadata import PackageNotFoundError, version


try:
    __version__ = version("yapCAD")
except PackageNotFoundError:  # pragma: no cover - handled when package not installed
    __version__ = "unknown"


def has_brep() -> bool:
    """Return True if solid-modeling (BREP/OCC) support is available.

    yapCAD's core is pure Python; the BREP features require OpenCASCADE via
    ``pythonocc-core`` (installed from conda-forge). Downstream code can branch
    on this to enable/disable solid-modeling paths cleanly::

        import yapcad
        if yapcad.has_brep():
            ...  # use boolean ops, STEP import, etc.

    The ``yapcad.brep`` module is imported lazily here so that importing
    ``yapcad`` never eagerly pulls in OCC.
    """
    from yapcad.brep import occ_available
    return occ_available()
