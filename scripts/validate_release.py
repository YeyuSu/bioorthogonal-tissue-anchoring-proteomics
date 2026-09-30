#!/usr/bin/env python3
"""Audit the manuscript repository for release-scope violations."""

from __future__ import annotations

import argparse
from pathlib import Path

IGNORED_PARTS = {".git", ".superpowers", ".venv", "__pycache__"}
BANNED_NAMES = {".DS_Store"}
BANNED_SUFFIXES = {".pyc", ".docx", ".xlsx"}
SCOPED_TEXT_DIRS = {"data", "configs", "scripts", "results"}
TEXT_SUFFIXES = {
    "",
    ".csv",
    ".json",
    ".md",
    ".py",
    ".svg",
    ".tsv",
    ".txt",
    ".yaml",
    ".yml",
}


def _candidate_files(root: Path):
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        relative = path.relative_to(root)
        if any(part in IGNORED_PARTS for part in relative.parts):
            continue
        yield path, relative


def _read_text(path: Path) -> str | None:
    if path.suffix.lower() not in TEXT_SUFFIXES:
        return None
    try:
        return path.read_text(encoding="utf-8")
    except (UnicodeDecodeError, OSError):
        return None


def audit_repository(root: Path) -> list[str]:
    """Return human-readable release violations found below *root*."""
    root = Path(root).resolve()
    errors: list[str] = []

    for path, relative in _candidate_files(root):
        shown = relative.as_posix()
        if path.name in BANNED_NAMES:
            errors.append(f"Forbidden file {shown}: {path.name}")
        if path.suffix.lower() in BANNED_SUFFIXES:
            errors.append(f"Forbidden extension {shown}: {path.suffix.lower()}")
        if "B16F10" in shown.upper():
            errors.append(f"Out-of-scope B16F10 material: {shown}")

        text = _read_text(path)
        if text is None:
            continue
        if "/Users/syysmac" in text:
            errors.append(f"Machine-specific path /Users/syysmac in {shown}")
        top_level = relative.parts[0] if relative.parts else ""
        if top_level in SCOPED_TEXT_DIRS and "B16F10" in text.upper():
            errors.append(f"Out-of-scope B16F10 content in {shown}")
        if top_level == "configs" and "NHS_vs_NCS" in text:
            errors.append(f"Forbidden NHS_vs_NCS contrast in {shown}")

    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", nargs="?", default=".", type=Path)
    args = parser.parse_args()
    errors = audit_repository(args.root)
    if errors:
        for error in errors:
            print(f"ERROR: {error}")
        return 1
    print("Release scope audit passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
