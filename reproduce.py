#!/usr/bin/env python3
"""Reproduce the 27-mouse proteomics analysis from the repository input."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import xml.etree.ElementTree as ET
from importlib import metadata
from pathlib import Path

from scripts.core.run_proteomics_workflow import run_workflow
from scripts.core.validate_proteomics_output import validate_output
from scripts.publication.build_publication_tables import build_publication_tables
from scripts.publication.render_publication_figures import render_publication_figures

EXPECTED_INPUT_SHA256 = "e36e4530191f46c37b1756d2563ef7200dbed870f8fa6121df37d34057ea3fb1"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _package_version(name: str) -> str:
    try:
        return metadata.version(name)
    except metadata.PackageNotFoundError:
        return "not installed"


def _svg_spec(path: Path) -> dict:
    root = ET.parse(path).getroot()
    svg = path.read_text(encoding="utf-8")
    return {
        "file": path.name,
        "width": root.attrib.get("width"),
        "height": root.attrib.get("height"),
        "editable_arial_text": "<text" in svg and "Arial" in svg,
    }


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
    figure_paths: list[Path] = []
    if not failures:
        publication_tables = output / "publication_tables"
        build_publication_tables(output, publication_tables)
        figure_paths = render_publication_figures(
            output, publication_tables, output / "publication_figures"
        )

    config = json.loads(config_path.read_text(encoding="utf-8"))
    volcano_specs = [_svg_spec(path) for path in figure_paths if "Volcano_" in path.name]
    for spec in volcano_specs:
        if spec["width"] != "44mm" or spec["height"] != "38mm":
            failures.append(f"Incorrect volcano dimensions: {spec['file']}")
        if not spec["editable_arial_text"]:
            failures.append(f"Volcano lacks editable Arial text: {spec['file']}")
    if not failures and len(volcano_specs) != 6:
        failures.append(f"Expected 6 publication volcano SVGs, found {len(volcano_specs)}")

    package_names = [
        "numpy", "pandas", "scipy", "statsmodels", "scikit-learn", "matplotlib"
    ]
    report = {
        "passed": not failures,
        "input_sha256": observed_hash,
        "python": sys.version,
        "packages": {name: _package_version(name) for name in package_names},
        "contrasts": [item["name"] for item in config["contrasts"]],
        "primary_excluded_samples": [
            str(value) for value in config["analysis"].get("exclude_samples", [])
        ],
        "sensitivity_restored_samples": (
            [str(value) for value in config["analysis"].get("exclude_samples", [])]
            if config.get("sensitivity", {}).get("restore_excluded_samples", False)
            else []
        ),
        "volcano_svgs": volcano_specs,
        "heatmap_svg": next(
            (path.name for path in figure_paths if "heatmap" in path.name), None
        ),
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
