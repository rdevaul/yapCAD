#!/usr/bin/env python3
"""Interactive viewer for yapCAD .ycpkg packages."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from yapcad.package.viewer import view_package


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Interactive viewer for .ycpkg packages.")
    parser.add_argument("package", type=Path, help="Path to the package directory.")
    parser.add_argument("--strict", action="store_true", help="Enable strict validation before viewing.")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    ok = view_package(args.package, strict=args.strict)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
