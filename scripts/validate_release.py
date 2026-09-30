#!/usr/bin/env python3
"""Audit the manuscript repository for reproducibility-release violations."""

from __future__ import annotations

import argparse
import hashlib
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
REFERENCE_HASHES = {
    "go-basic.obo": "b08d45b268b8c24ccb2513dbbbc7d4df9f6521c099b413f79eb31e06e0fa3bcc",
    "MOUSE-uniprot.gaf.gz": "ea07d1a19adad445da59b2d8c48e7f7678827b09c12cc056f3faf7a1f88068b8",
}
OUT_OF_SCOPE_MARKER = "B16" + "F10"
FORBIDDEN_CONTRAST = "NHS" + "_vs_" + "NCS"
MACHINE_PATH = "/" + "Users" + "/" + "syysmac"


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


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def audit_repository(root: Path) -> list[str]:
    """Return human-readable scope violations found below *root*."""
    root = Path(root).resolve()
    errors: list[str] = []

    for path, relative in _candidate_files(root):
        shown = relative.as_posix()
        if path.name in BANNED_NAMES:
            errors.append(f"Forbidden file {shown}: {path.name}")
        if path.suffix.lower() in BANNED_SUFFIXES:
            errors.append(f"Forbidden extension {shown}: {path.suffix.lower()}")
        if OUT_OF_SCOPE_MARKER in shown.upper():
            errors.append(f"Out-of-scope {OUT_OF_SCOPE_MARKER} material: {shown}")

        text = _read_text(path)
        if text is None:
            continue
        if MACHINE_PATH in text:
            errors.append(f"Machine-specific path {MACHINE_PATH} in {shown}")
        top_level = relative.parts[0] if relative.parts else ""
        if top_level in SCOPED_TEXT_DIRS and OUT_OF_SCOPE_MARKER in text.upper():
            errors.append(f"Out-of-scope {OUT_OF_SCOPE_MARKER} content in {shown}")
        if top_level == "configs" and FORBIDDEN_CONTRAST in text:
            errors.append(f"Forbidden {FORBIDDEN_CONTRAST} contrast in {shown}")

    return errors


def verify_reference_hashes(root: Path) -> list[str]:
    root = Path(root).resolve()
    errors: list[str] = []
    for filename, expected in REFERENCE_HASHES.items():
        path = root / "reference" / filename
        if not path.is_file():
            errors.append(f"Missing frozen reference file: reference/{filename}")
            continue
        observed = _sha256(path)
        if observed != expected:
            errors.append(
                f"Reference SHA256 mismatch for {filename}: expected {expected}, observed {observed}"
            )
    return errors


def verify_license_coverage(root: Path) -> list[str]:
    root = Path(root).resolve()
    errors: list[str] = []
    code_path = root / "LICENSE-CODE"
    data_path = root / "LICENSE-DATA"
    readme_path = root / "README.md"
    if not code_path.is_file() or "MIT License" not in code_path.read_text(encoding="utf-8"):
        errors.append("Missing or invalid MIT code license: LICENSE-CODE")
    if not data_path.is_file():
        errors.append("Missing data license: LICENSE-DATA")
    else:
        data_license = data_path.read_text(encoding="utf-8")
        if "CC BY 4.0" not in data_license or "creativecommons.org/licenses/by/4.0" not in data_license:
            errors.append("LICENSE-DATA does not identify CC BY 4.0 and its canonical URL")
    if not readme_path.is_file():
        errors.append("Missing README.md")
    return errors


def validate_release(root: Path) -> list[str]:
    return (
        audit_repository(root)
        + verify_reference_hashes(root)
        + verify_license_coverage(root)
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", nargs="?", default=".", type=Path)
    args = parser.parse_args()
    errors = validate_release(args.root)
    if errors:
        for error in errors:
            print(f"ERROR: {error}")
        return 1
    print("Release audit passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
