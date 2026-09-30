#!/usr/bin/env python3
"""Validate a proteomics workflow output bundle."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import re
import sys


def validate_output(output: Path) -> list[str]:
    output = Path(output)
    failures: list[str] = []
    required = [
        output / "code" / "run_proteomics_workflow.py",
        output / "code" / "analysis_config.json",
        output / "tables" / "sample_metadata.csv",
        output / "tables" / "sample_qc.csv",
        output / "tables" / "differential_proteins.csv",
        output / "tables" / "significant_proteins.csv",
        output / "analysis_report.md",
        output / "provenance.json",
    ]
    for path in required:
        if not path.is_file() or path.stat().st_size == 0:
            failures.append(f"missing or empty required file: {path.relative_to(output)}")
    config_path = output / "code" / "analysis_config.json"
    if config_path.is_file():
        config = json.loads(config_path.read_text(encoding="utf-8"))
        if config.get("sensitivity", {}).get("restore_excluded_samples", False):
            sensitivity_path = output / "tables" / "sensitivity_summary.csv"
            if not sensitivity_path.is_file() or sensitivity_path.stat().st_size == 0:
                failures.append("missing or empty required file: tables/sensitivity_summary.csv")

        differential_path = output / "tables" / "differential_proteins.csv"
        if differential_path.is_file():
            with differential_path.open(encoding="utf-8", newline="") as handle:
                differential_rows = list(csv.DictReader(handle))
            if not differential_rows:
                failures.append("differential_proteins.csv contains no protein results")

        metadata_path = output / "tables" / "sample_metadata.csv"
        if metadata_path.is_file():
            with metadata_path.open(encoding="utf-8", newline="") as handle:
                metadata_rows = list(csv.DictReader(handle))
            observed_excluded = {
                str(row.get("sample_id", ""))
                for row in metadata_rows
                if str(row.get("excluded_from_primary", "")).lower() == "true"
            }
            configured_excluded = {
                str(value) for value in config.get("analysis", {}).get("exclude_samples", [])
            }
            if observed_excluded != configured_excluded:
                failures.append(
                    "sample_metadata primary exclusions do not match analysis_config.json"
                )

        plot = config.get("plot", {})
        width = plot.get("width_mm")
        height = plot.get("height_mm")
        font = plot.get("font_family", "Arial")
        tick_font = plot.get("tick_font_pt")
        axis_font = plot.get("axis_label_font_pt")
        svg_files = sorted((output / "figures" / "standalone_volcanoes").glob("*.svg"))
        expected_count = len(config.get("contrasts", []))
        if len(svg_files) != expected_count:
            failures.append(f"expected {expected_count} volcano SVGs, found {len(svg_files)}")
        for svg_path in svg_files:
            svg = svg_path.read_text(encoding="utf-8")
            if f'width="{width:g}mm"' not in svg or f'height="{height:g}mm"' not in svg:
                failures.append(f"incorrect physical dimensions: {svg_path.name}")
            if font not in svg or "<text" not in svg:
                failures.append(f"missing editable {font} text: {svg_path.name}")
            if f"font-size: {tick_font:g}px" not in svg:
                failures.append(f"incorrect tick font size: {svg_path.name}")
            if f"font-size: {axis_font:g}px" not in svg:
                failures.append(f"incorrect axis-label font size: {svg_path.name}")
    return failures


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    failures = validate_output(args.output)
    if failures:
        for failure in failures:
            print(f"FAIL: {failure}", file=sys.stderr)
        raise SystemExit(1)
    print(f"Validated proteomics output: {args.output}")


if __name__ == "__main__":
    main()
