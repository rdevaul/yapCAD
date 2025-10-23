#!/usr/bin/env python3
"""Command-line validator for yapCAD .ycpkg packages."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from yapcad.package.validator import validate_package


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Validate a yapCAD .ycpkg package.")
    parser.add_argument("package", type=Path, help="Path to the .ycpkg directory.")
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Enable additional checks (e.g., require hashes on exports/attachments).",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    ok, messages = validate_package(args.package, strict=args.strict)
    for msg in messages:
        print(msg)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
