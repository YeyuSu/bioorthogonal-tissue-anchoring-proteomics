#!/usr/bin/env python3
"""Reproduce the 27-mouse proteomics analysis from the repository input."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from scripts.core.run_proteomics_workflow import run_workflow
from scripts.core.validate_proteomics_output import validate_output
from scripts.publication.build_publication_tables import build_publication_tables

EXPECTED_INPUT_SHA256 = "e36e4530191f46c37b1756d2563ef7200dbed870f8fa6121df37d34057ea3fb1"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_reproduction(root: Path, output: Path) -> Path:
    root = Path(root).resolve()
    output = Path(output).resolve()
    input_path = root / "data" / "27_mice_protein_groups.tsv"
    config_path = root / "configs" / "27_mice.json"

    if not input_path.is_file():
        raise FileNotFoundError(f"Missing input TSV: {input_path}")
    observed_hash = sha256_file(input_path)
    if observed_hash != EXPECTED_INPUT_SHA256:
        raise ValueError(
            "Input TSV SHA256 mismatch: "
            f"expected {EXPECTED_INPUT_SHA256}, observed {observed_hash}"
        )
    if not config_path.is_file():
        raise FileNotFoundError(f"Missing analysis config: {config_path}")
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {output}")

    run_workflow(input_path, config_path, output)
    failures = validate_output(output)
    if not failures:
        build_publication_tables(output, output / "publication_tables")
    report = {
        "passed": not failures,
        "input_sha256": observed_hash,
        "failures": failures,
    }
    (output / "validation_report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    if failures:
        raise RuntimeError("Output validation failed: " + "; ".join(failures))
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent / "output",
    )
    args = parser.parse_args()
    result = run_reproduction(Path(__file__).resolve().parent, args.output)
    print(f"Reproduction completed: {result}")


if __name__ == "__main__":
    main()
