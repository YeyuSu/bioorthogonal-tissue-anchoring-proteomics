#!/usr/bin/env python3
"""Audit the manuscript repository for reproducibility-release violations."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
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


EXPECTED_INPUT_SHA256 = "e36e4530191f46c37b1756d2563ef7200dbed870f8fa6121df37d34057ea3fb1"
CONTRASTS = [
    "subcutaneous_NHS_vs_biotin",
    "subcutaneous_NCS_vs_biotin",
    "muscle_NHS_vs_biotin",
    "muscle_NCS_vs_biotin",
    "periodontal_NHS_vs_biotin",
    "periodontal_NCS_vs_biotin",
]
CORE_REFERENCE_TABLES = [
    "sample_metadata.csv",
    "sample_qc.csv",
    "differential_proteins.csv",
    "significant_proteins.csv",
    "sensitivity_differential_proteins.csv",
    "sensitivity_summary.csv",
]
PUBLICATION_REFERENCE_TABLES = [
    "candidate_summary.csv",
    "top10_enriched_proteins.csv",
    "top10_membership_matrix.csv",
    "top10_overlap_summary.csv",
    "three_tissue_effect_heatmap_log2FC.csv",
    "three_tissue_effect_heatmap_q_BH.csv",
    "three_tissue_effect_heatmap_labels.csv",
]


def _reference_file_pairs(generated: Path, reference: Path) -> list[tuple[Path, Path]]:
    pairs = [
        (generated / "tables" / name, reference / "tables" / name)
        for name in CORE_REFERENCE_TABLES
    ]
    pairs.extend(
        (generated / "publication_tables" / name, reference / "tables" / name)
        for name in PUBLICATION_REFERENCE_TABLES
    )
    pairs.extend(
        (
            generated / "publication_figures" / "standalone_volcanoes" / f"Volcano_{name}.svg",
            reference / "figures" / f"Volcano_{name}.svg",
        )
        for name in CONTRASTS
    )
    pairs.append(
        (
            generated / "publication_figures" / "three_tissue_treatment_effect_heatmap.svg",
            reference / "figures" / "three_tissue_treatment_effect_heatmap.svg",
        )
    )
    pairs.append((generated / "validation_report.json", reference / "validation_report.json"))
    return pairs


def promote_reference_results(generated: Path, reference: Path) -> list[Path]:
    generated = Path(generated).resolve()
    reference = Path(reference).resolve()
    if reference.exists() and any(path.is_file() for path in reference.rglob("*")):
        raise FileExistsError(f"Reference directory is not empty: {reference}")
    pairs = _reference_file_pairs(generated, reference)
    missing = [source for source, _ in pairs if not source.is_file()]
    provenance_source = generated / "provenance.json"
    if not provenance_source.is_file():
        missing.append(provenance_source)
    if missing:
        raise FileNotFoundError("Missing generated artifacts: " + ", ".join(map(str, missing)))

    copied: list[Path] = []
    for source, destination in pairs:
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        copied.append(destination)

    provenance = json.loads(provenance_source.read_text(encoding="utf-8"))
    provenance["input_path"] = "data/27_mice_protein_groups.tsv"
    provenance["config_path"] = "configs/27_mice.json"
    provenance_destination = reference / "provenance.json"
    provenance_destination.write_text(
        json.dumps(provenance, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    copied.append(provenance_destination)
    return copied


def verify_reference_results(root: Path) -> list[str]:
    root = Path(root).resolve()
    reference = root / "results" / "reference"
    errors: list[str] = []
    expected = [destination for _, destination in _reference_file_pairs(root, reference)]
    expected.append(reference / "provenance.json")
    for path in expected:
        if not path.is_file() or path.stat().st_size == 0:
            errors.append(f"Missing or empty reference result: {path.relative_to(root)}")

    report_path = reference / "validation_report.json"
    if report_path.is_file():
        report = json.loads(report_path.read_text(encoding="utf-8"))
        if report.get("passed") is not True:
            errors.append("Reference validation report does not pass")
        if report.get("input_sha256") != EXPECTED_INPUT_SHA256:
            errors.append("Reference validation report has the wrong input SHA256")
        if report.get("contrasts") != CONTRASTS:
            errors.append("Reference validation report does not list the six contrasts in order")
        if report.get("primary_excluded_samples") != ["17"]:
            errors.append("Reference validation report does not exclude mouse 17")
        if report.get("sensitivity_restored_samples") != ["17"]:
            errors.append("Reference validation report does not restore mouse 17")
        specs = report.get("volcano_svgs", [])
        if len(specs) != 6:
            errors.append("Reference validation report does not describe six volcano SVGs")
        for spec in specs:
            if spec.get("width") != "44mm" or spec.get("height") != "38mm":
                errors.append(f"Incorrect reference volcano dimensions: {spec.get('file')}")
            if spec.get("editable_arial_text") is not True:
                errors.append(f"Reference volcano lacks editable Arial text: {spec.get('file')}")
    return errors


def compare_reference_results(generated: Path, reference: Path) -> list[str]:
    generated = Path(generated).resolve()
    reference = Path(reference).resolve()
    errors: list[str] = []
    for source, destination in _reference_file_pairs(generated, reference):
        if not source.is_file() or not destination.is_file():
            errors.append(f"Missing comparison artifact: {destination.name}")
        elif source.read_bytes() != destination.read_bytes():
            errors.append(f"Reference artifact differs from generated output: {destination.name}")
    return errors


def validate_release(root: Path) -> list[str]:
    return (
        audit_repository(root)
        + verify_reference_hashes(root)
        + verify_license_coverage(root)
        + verify_reference_results(root)
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
