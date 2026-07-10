# -*- coding: utf-8 -*-
import warnings as _warnings

try:  # Python >= 3.8
    from importlib.metadata import PackageNotFoundError, version
except ModuleNotFoundError:  # pragma: no cover - for Python < 3.8
    from importlib_metadata import PackageNotFoundError, version


try:
    __version__ = version("yapCAD")
except PackageNotFoundError:  # pragma: no cover - handled when package not installed
    __version__ = "unknown"


class YapcadBrepUnavailableWarning(UserWarning):
    """Warned once at import time when solid-modeling (BREP/OCC) support is absent.

    yapCAD's PyPI/pip install ships the pure-Python core only. Boolean
    operations, STEP import/export, and solid (BREP) geometry require
    OpenCASCADE via ``pythonocc-core``, which is distributed through
    conda-forge — not pip. This warning flags that limited-functionality state.
    Suppress it with the standard ``warnings`` machinery, e.g.::

        import warnings
        from yapcad import YapcadBrepUnavailableWarning
        warnings.simplefilter("ignore", YapcadBrepUnavailableWarning)
    """


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


def _warn_if_brep_unavailable() -> None:
    """Emit a one-time, suppressible warning at import when BREP/OCC is absent.

    Fires exactly once (relies on the default ``warnings`` "once/​default"
    dedup per category+message+location). Non-fatal; never raised when OCC is
    present. Any failure to probe OCC is swallowed so that importing yapCAD is
    always safe.
    """
    try:
        available = has_brep()
    except Exception:  # pragma: no cover - probing must never break import
        return
    if available:
        return
    _warnings.warn(
        "yapCAD was installed without BREP support (pythonocc-core not found). "
        "Boolean operations, STEP import/export, and solid (BREP) geometry are "
        "unavailable; the pure-Python core (2D geometry, the DSL, metadata, and "
        "the assembly graph) still works. To enable full functionality install "
        "pythonocc-core from conda-forge: "
        "`conda install -c conda-forge pythonocc-core`. "
        "See https://github.com/rdevaul/yapCAD#installation for details.",
        YapcadBrepUnavailableWarning,
        stacklevel=2,
    )


_warn_if_brep_unavailable()
