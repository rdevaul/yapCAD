# -*- coding: utf-8 -*-
try:  # Python >= 3.8
    from importlib.metadata import PackageNotFoundError, version
except ModuleNotFoundError:  # pragma: no cover - for Python < 3.8
    from importlib_metadata import PackageNotFoundError, version


try:
    __version__ = version("yapCAD")
except PackageNotFoundError:  # pragma: no cover - handled when package not installed
    __version__ = "unknown"
