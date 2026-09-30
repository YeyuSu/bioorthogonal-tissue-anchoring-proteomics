#!/usr/bin/env python3
"""Configurable protein-group analysis and exact-size volcano plotting."""

from __future__ import annotations

import argparse
import copy
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
from typing import Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


def load_config(path: Path) -> dict:
    """Load a JSON workflow configuration."""
    with Path(path).open(encoding="utf-8") as handle:
        config = json.load(handle)
    if not isinstance(config, dict):
        raise ValueError("configuration root must be a JSON object")
    return config


def resolve_samples(config: dict, columns: Sequence[str]) -> list[dict]:
    """Return validated sample records in configured order."""
    validate_config(config, columns)
    return [dict(sample) for sample in config["samples"]]


def validate_config(config: dict, columns: Sequence[str]) -> None:
    """Validate fields that define sample mapping, contrasts, and plotting."""
    for field in ("columns", "samples", "contrasts", "normalization", "analysis", "plot"):
        if field not in config:
            raise ValueError(f"missing configuration field: {field}")

    column_set = set(columns)
    for field in ("protein_id", "description"):
        name = config["columns"].get(field)
        if not name or name not in column_set:
            raise ValueError(f"configured {field} column is missing: {name}")

    samples = config["samples"]
    if not samples:
        raise ValueError("samples must contain at least one sample record")
    sample_ids = [str(sample.get("sample_id", "")) for sample in samples]
    sample_columns = [sample.get("column") for sample in samples]
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError("duplicate sample_id in samples")
    if len(sample_columns) != len(set(sample_columns)):
        raise ValueError("duplicate sample column in samples")
    for sample in samples:
        for field in ("sample_id", "column", "tissue", "treatment"):
            if field not in sample or sample[field] in (None, ""):
                raise ValueError(f"sample record missing field: {field}")
        if sample["column"] not in column_set:
            raise ValueError(f"configured sample column is missing: {sample['column']}")

    groups = {(sample["tissue"], sample["treatment"]) for sample in samples}
    contrast_names: set[str] = set()
    for contrast in config["contrasts"]:
        for field in ("name", "tissue", "treated", "control"):
            if field not in contrast or contrast[field] in (None, ""):
                raise ValueError(f"contrast missing field: {field}")
        if contrast["name"] in contrast_names:
            raise ValueError(f"duplicate contrast name: {contrast['name']}")
        contrast_names.add(contrast["name"])
        for treatment_field in ("treated", "control"):
            group = (contrast["tissue"], contrast[treatment_field])
            if group not in groups:
                raise ValueError(
                    f"contrast {contrast['name']} references undefined group: "
                    f"{contrast['tissue']} / {contrast[treatment_field]}"
                )

    plot = config["plot"]
    for field in ("width_mm", "height_mm", "tick_font_pt", "axis_label_font_pt",
                  "protein_label_font_pt", "title_font_pt"):
        value = plot.get(field)
        if not isinstance(value, (int, float)) or value <= 0:
            raise ValueError(f"plot.{field} must be a positive number")
    x_limits = plot.get("x_limits")
    if not isinstance(x_limits, list) or len(x_limits) != 2 or x_limits[0] >= x_limits[1]:
        raise ValueError("plot.x_limits must contain two increasing numbers")


def ensure_output_directory(output_dir: Path) -> Path:
    """Create an empty output directory and refuse silent replacement."""
    output_dir = Path(output_dir)
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"output directory is not empty: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def _contaminant_mask(table: pd.DataFrame, config: dict) -> pd.Series:
    protein_column = config["columns"]["protein_id"]
    prefixes = tuple(config["analysis"].get("contaminant_prefixes", []))
    if not prefixes:
        return pd.Series(False, index=table.index)
    return table[protein_column].fillna("").astype(str).str.startswith(prefixes)


def extract_reference_values(
    table: pd.DataFrame, config: dict, samples: Sequence[dict]
) -> pd.Series:
    """Extract one finite positive normalization-reference value per sample."""
    normalization = config["normalization"]
    match_column = normalization["reference_match_column"]
    match_value = normalization["reference_match_value"]
    if match_column not in table.columns:
        raise ValueError(f"reference match column is missing: {match_column}")
    matches = table.loc[table[match_column].astype(str) == str(match_value)]
    if matches.empty:
        raise ValueError(f"reference row not found: {match_column}={match_value}")
    if len(matches) > 1:
        raise ValueError(
            f"multiple reference rows found for {match_column}={match_value}; "
            "configure an explicit aggregation rule"
        )
    row = matches.iloc[0]
    values: dict[str, float] = {}
    for sample in samples:
        column = sample["column"]
        value = pd.to_numeric(pd.Series([row[column]]), errors="coerce").iloc[0]
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"finite positive reference abundance required for sample {column}")
        values[column] = float(value)
    return pd.Series(values, dtype=float, name="reference_abundance")


def normalize_table(
    table: pd.DataFrame, reference_values: pd.Series, samples: Sequence[dict]
) -> pd.DataFrame:
    """Divide each abundance column by its sample-specific reference value."""
    normalized = table.copy()
    for sample in samples:
        column = sample["column"]
        values = pd.to_numeric(normalized[column], errors="coerce")
        normalized[column] = values / float(reference_values[column])
    return normalized


def bh_adjust(p_values) -> np.ndarray:
    """Benjamini-Hochberg adjustment that preserves missing positions."""
    p = np.asarray(p_values, dtype=float)
    adjusted = np.full(p.shape, np.nan, dtype=float)
    finite_positions = np.flatnonzero(np.isfinite(p))
    if finite_positions.size == 0:
        return adjusted
    finite = p[finite_positions]
    order = np.argsort(finite)
    ranked = finite[order]
    n = ranked.size
    corrected = ranked * n / np.arange(1, n + 1)
    corrected = np.minimum.accumulate(corrected[::-1])[::-1]
    corrected = np.clip(corrected, 0, 1)
    restored = np.empty_like(corrected)
    restored[order] = corrected
    adjusted[finite_positions] = restored
    return adjusted


def build_sample_qc(
    table: pd.DataFrame, samples: Sequence[dict], config: dict
) -> pd.DataFrame:
    """Build transparent per-sample abundance and detection summaries."""
    reference = extract_reference_values(table, config, samples)
    contaminant = _contaminant_mask(table, config)
    rows: list[dict] = []
    excluded = {str(value) for value in config["analysis"].get("exclude_samples", [])}
    for sample in samples:
        column = sample["column"]
        values = pd.to_numeric(table[column], errors="coerce")
        detected = values.notna() & (values > 0)
        non_contaminant = ~contaminant
        non_contaminant_total = float(values.where(detected & non_contaminant).sum())
        total = float(values.where(detected).sum())
        contaminant_total = float(values.where(detected & contaminant).sum())
        rows.append(
            {
                **sample,
                "excluded_from_primary": str(sample["sample_id"]) in excluded,
                "reference_abundance": float(reference[column]),
                "detected_non_contaminant_proteins": int((detected & non_contaminant).sum()),
                "total_non_contaminant_abundance": non_contaminant_total,
                "normalized_non_contaminant_total": non_contaminant_total / float(reference[column]),
                "contaminant_fraction": contaminant_total / total if total > 0 else np.nan,
            }
        )
    return pd.DataFrame(rows)


def _welch_confidence_interval(treated: np.ndarray, control: np.ndarray) -> tuple[float, float]:
    mean_difference = float(np.mean(treated) - np.mean(control))
    treated_term = float(np.var(treated, ddof=1) / len(treated))
    control_term = float(np.var(control, ddof=1) / len(control))
    standard_error = math.sqrt(treated_term + control_term)
    if standard_error == 0:
        return mean_difference, mean_difference
    numerator = (treated_term + control_term) ** 2
    denominator = (treated_term ** 2 / (len(treated) - 1)) + (
        control_term ** 2 / (len(control) - 1)
    )
    degrees_freedom = numerator / denominator if denominator > 0 else np.inf
    critical = float(stats.t.ppf(0.975, degrees_freedom))
    return mean_difference - critical * standard_error, mean_difference + critical * standard_error


def calculate_contrast(
    normalized: pd.DataFrame,
    metadata: Sequence[dict],
    contrast: dict,
    config: dict,
) -> pd.DataFrame:
    """Calculate one treated-minus-control differential protein table."""
    excluded = {str(value) for value in config["analysis"].get("exclude_samples", [])}
    active = [sample for sample in metadata if str(sample["sample_id"]) not in excluded]
    treated_columns = [
        sample["column"]
        for sample in active
        if sample["tissue"] == contrast["tissue"] and sample["treatment"] == contrast["treated"]
    ]
    control_columns = [
        sample["column"]
        for sample in active
        if sample["tissue"] == contrast["tissue"] and sample["treatment"] == contrast["control"]
    ]
    if not treated_columns or not control_columns:
        raise ValueError(f"contrast has an empty group after exclusions: {contrast['name']}")

    protein_id_column = config["columns"]["protein_id"]
    description_column = config["columns"]["description"]
    minimum = int(config["analysis"].get("minimum_detections_per_group", 2))
    biological = normalized.loc[~_contaminant_mask(normalized, config)]
    records: list[dict] = []
    for _, row in biological.iterrows():
        treated_raw = pd.to_numeric(row[treated_columns], errors="coerce").to_numpy(dtype=float)
        control_raw = pd.to_numeric(row[control_columns], errors="coerce").to_numpy(dtype=float)
        treated = np.log2(treated_raw[np.isfinite(treated_raw) & (treated_raw > 0)])
        control = np.log2(control_raw[np.isfinite(control_raw) & (control_raw > 0)])
        if len(treated) < minimum or len(control) < minimum:
            continue
        log2fc = float(np.mean(treated) - np.mean(control))
        test = stats.ttest_ind(treated, control, equal_var=False, nan_policy="omit")
        p_value = float(test.pvalue)
        if not np.isfinite(p_value):
            p_value = 1.0
        ci_low, ci_high = _welch_confidence_interval(treated, control)
        records.append(
            {
                "contrast": contrast["name"],
                "tissue": contrast["tissue"],
                "treated": contrast["treated"],
                "control": contrast["control"],
                "protein_id": row[protein_id_column],
                "description": row[description_column],
                "n_treated": len(treated),
                "n_control": len(control),
                "mean_log2_treated": float(np.mean(treated)),
                "mean_log2_control": float(np.mean(control)),
                "log2FC": log2fc,
                "ci95_low": ci_low,
                "ci95_high": ci_high,
                "p_Welch": p_value,
            }
        )
    result = pd.DataFrame(records)
    if result.empty:
        return pd.DataFrame(
            columns=[
                "contrast", "tissue", "treated", "control", "protein_id", "description",
                "n_treated", "n_control", "mean_log2_treated", "mean_log2_control",
                "log2FC", "ci95_low", "ci95_high", "p_Welch", "q_BH",
                "significant_enriched", "significant_depleted",
            ]
        )
    result["q_BH"] = bh_adjust(result["p_Welch"].to_numpy())
    fold_threshold = float(config["analysis"].get("absolute_log2fc_threshold", 1.0))
    q_threshold = float(config["analysis"].get("bh_adjusted_p_threshold", 0.05))
    result["significant_enriched"] = (result["q_BH"] < q_threshold) & (
        result["log2FC"] >= fold_threshold
    )
    result["significant_depleted"] = (result["q_BH"] < q_threshold) & (
        result["log2FC"] <= -fold_threshold
    )
    return result


def _write_svg_size(path: Path, width_mm: float, height_mm: float) -> None:
    svg = path.read_text(encoding="utf-8")
    svg = re.sub(r'width="[^"]+"', f'width="{width_mm:g}mm"', svg, count=1)
    svg = re.sub(r'height="[^"]+"', f'height="{height_mm:g}mm"', svg, count=1)
    path.write_text(svg, encoding="utf-8")


def _safe_filename(value: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return safe.strip("_") or "contrast"


def render_volcano(
    results: pd.DataFrame,
    contrast: dict,
    plot_config: dict,
    output: Path,
) -> Path:
    """Render one exact-size editable SVG volcano plot."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    width_mm = float(plot_config["width_mm"])
    height_mm = float(plot_config["height_mm"])
    font = str(plot_config.get("font_family", "Arial"))
    tick_font = float(plot_config["tick_font_pt"])
    axis_font = float(plot_config["axis_label_font_pt"])
    protein_font = float(plot_config["protein_label_font_pt"])
    title_font = float(plot_config["title_font_pt"])
    x_limits = [float(value) for value in plot_config["x_limits"]]

    frame = results.copy()
    frame["minus_log10_p"] = -np.log10(
        pd.to_numeric(frame["p_Welch"], errors="coerce").clip(lower=1e-300)
    )
    enriched = frame["significant_enriched"].fillna(False).astype(bool)
    depleted = frame["significant_depleted"].fillna(False).astype(bool)
    neutral = ~(enriched | depleted)
    colors = plot_config.get("colors", {})

    rc = {
        "font.family": font,
        "font.sans-serif": [font],
        "font.size": tick_font,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.unicode_minus": False,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "svg.fonttype": "none",
    }
    with plt.rc_context(rc):
        fig, ax = plt.subplots(figsize=(width_mm / 25.4, height_mm / 25.4))
        fig.subplots_adjust(left=.22, right=.97, bottom=.20, top=.80)
        ax.scatter(
            frame.loc[neutral, "log2FC"], frame.loc[neutral, "minus_log10_p"],
            s=float(plot_config.get("neutral_point_size", 8)),
            color=colors.get("neutral", "#a9afb4"), alpha=.48, linewidth=0,
        )
        ax.scatter(
            frame.loc[depleted, "log2FC"], frame.loc[depleted, "minus_log10_p"],
            s=float(plot_config.get("significant_point_size", 22)),
            color=colors.get("depleted", "#2878b5"), alpha=.92,
            edgecolor="white", linewidth=.25,
        )
        ax.scatter(
            frame.loc[enriched, "log2FC"], frame.loc[enriched, "minus_log10_p"],
            s=float(plot_config.get("significant_point_size", 22)),
            color=colors.get("enriched", "#c83f3b"), alpha=.92,
            edgecolor="white", linewidth=.25,
        )
        fold_threshold = float(plot_config.get("absolute_log2fc_threshold", 1.0))
        ax.axvline(-fold_threshold, color="#969da3", ls="--", lw=.55)
        ax.axvline(fold_threshold, color="#969da3", ls="--", lw=.55)
        ax.axhline(-math.log10(.05), color="#969da3", ls=":", lw=.55)
        ax.set_xlim(*x_limits)
        observed_max = float(frame["minus_log10_p"].max()) if not frame.empty else 0.0
        ax.set_ylim(-.08, max(3.5, observed_max + .75))
        ax.set_xlabel("log2 fold change", fontsize=axis_font, labelpad=2)
        ax.set_ylabel("−log10(Welch P)", fontsize=axis_font, labelpad=2)
        ax.tick_params(axis="both", which="major", labelsize=tick_font, length=2.2, width=.5, pad=1.2)
        ax.spines["left"].set_linewidth(.55)
        ax.spines["bottom"].set_linewidth(.55)
        tick_step = float(plot_config.get("x_tick_step", 5))
        ax.set_xticks(np.arange(x_limits[0], x_limits[1] + tick_step / 2, tick_step))
        title = str(contrast.get("title", contrast["name"])).replace("\\n", "\n")
        ax.set_title(title, loc="left", fontweight="bold", fontsize=title_font, pad=2)

        significant = frame.loc[enriched | depleted].sort_values(
            ["minus_log10_p", "log2FC"], ascending=[False, False]
        )
        vertical_offsets = [8, 2, -3, -10, 12, -14]
        label_replacements = plot_config.get("label_replacements", {})
        midpoint = (x_limits[0] + x_limits[1]) / 2
        span = x_limits[1] - x_limits[0]
        for index, (_, row) in enumerate(significant.iterrows()):
            x = float(row["log2FC"])
            near_right = x > midpoint + .15 * span
            near_left = x < midpoint - .15 * span
            if near_right:
                dx, align = -5, "right"
            elif near_left:
                dx, align = 5, "left"
            else:
                dx, align = 5, "left"
            ax.annotate(
                str(label_replacements.get(str(row["description"]), row["description"])),
                (x, float(row["minus_log10_p"])),
                xytext=(dx, vertical_offsets[index % len(vertical_offsets)]),
                textcoords="offset points", ha=align, va="center",
                fontsize=protein_font, annotation_clip=False,
                arrowprops={"arrowstyle": "-", "color": "#555", "lw": .35},
                zorder=6,
            )

        path = output / f"Volcano_{_safe_filename(contrast['name'])}.svg"
        fig.savefig(path)
        plt.close(fig)
    _write_svg_size(path, width_mm, height_mm)
    return path


def render_qc_figures(sample_qc: pd.DataFrame, config: dict, output: Path) -> list[Path]:
    """Render compact sample-quality and normalization-support SVGs."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    font = str(config["plot"].get("font_family", "Arial"))
    rc = {
        "font.family": font,
        "font.sans-serif": [font],
        "font.size": 7,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "svg.fonttype": "none",
    }
    paths: list[Path] = []
    labels = sample_qc["sample_id"].astype(str).tolist()
    x = np.arange(len(labels))
    with plt.rc_context(rc):
        fig, axes = plt.subplots(2, 1, figsize=(90 / 25.4, 62 / 25.4), sharex=True)
        axes[0].bar(x, sample_qc["reference_abundance"], color="#8297aa")
        axes[0].set_ylabel("Reference intensity")
        axes[1].bar(x, sample_qc["detected_non_contaminant_proteins"], color="#6ba8a5")
        axes[1].set_ylabel("Proteins detected")
        axes[1].set_xticks(x, labels, rotation=90)
        fig.tight_layout(pad=.7)
        path = output / "QC_sample_signals.svg"
        fig.savefig(path)
        plt.close(fig)
        _write_svg_size(path, 90, 62)
        paths.append(path)

        fig, axes = plt.subplots(2, 1, figsize=(90 / 25.4, 62 / 25.4), sharex=True)
        axes[0].bar(x, sample_qc["total_non_contaminant_abundance"], color="#b1b6bb")
        axes[0].set_ylabel("Total abundance")
        axes[1].bar(x, sample_qc["normalized_non_contaminant_total"], color="#5e9d99")
        axes[1].set_ylabel("Normalized total")
        axes[1].set_xticks(x, labels, rotation=90)
        fig.tight_layout(pad=.7)
        path = output / "QC_normalization_effect.svg"
        fig.savefig(path)
        plt.close(fig)
        _write_svg_size(path, 90, 62)
        paths.append(path)
    return paths


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_abundance_table(input_path: Path, config: dict) -> pd.DataFrame:
    delimiter = config.get("delimiter")
    if delimiter is None:
        delimiter = "\t" if input_path.suffix.lower() in {".tsv", ".txt"} else ","
    table = pd.read_csv(input_path, sep=delimiter)
    validate_config(config, table.columns.tolist())
    for sample in config["samples"]:
        table[sample["column"]] = pd.to_numeric(table[sample["column"]], errors="coerce")
    description_column = config["columns"]["description"]
    separator = config["columns"].get("description_separator")
    if separator:
        table[description_column] = (
            table[description_column].fillna("").astype(str).str.split(separator, n=1).str[0]
        )
    return table


def _metadata_table(config: dict) -> pd.DataFrame:
    excluded = {str(value) for value in config["analysis"].get("exclude_samples", [])}
    rows = []
    for sample in config["samples"]:
        rows.append(
            {
                **sample,
                "excluded_from_primary": str(sample["sample_id"]) in excluded,
            }
        )
    return pd.DataFrame(rows)


def _sensitivity_outputs(
    normalized: pd.DataFrame, config: dict
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if not config.get("sensitivity", {}).get("restore_excluded_samples", False):
        return pd.DataFrame(), pd.DataFrame()
    sensitivity_config = copy.deepcopy(config)
    sensitivity_config["analysis"]["exclude_samples"] = []
    result_tables: list[pd.DataFrame] = []
    summary_rows: list[dict] = []
    for contrast in sensitivity_config["contrasts"]:
        result = calculate_contrast(
            normalized, sensitivity_config["samples"], contrast, sensitivity_config
        )
        result.insert(0, "analysis", "sensitivity_restore_excluded")
        result_tables.append(result)
        included = [
            str(sample["sample_id"])
            for sample in sensitivity_config["samples"]
            if sample["tissue"] == contrast["tissue"]
            and sample["treatment"] in {contrast["treated"], contrast["control"]}
        ]
        significant = result["significant_enriched"] | result["significant_depleted"]
        summary_rows.append(
            {
                "contrast": contrast["name"],
                "included_samples": ",".join(included),
                "tested_proteins": len(result),
                "significant_proteins": int(significant.sum()),
            }
        )
    return pd.concat(result_tables, ignore_index=True), pd.DataFrame(summary_rows)


def _write_report(
    path: Path,
    config: dict,
    differential: pd.DataFrame,
    significant: pd.DataFrame,
) -> None:
    lines = [
        "# Proteomics workflow report",
        "",
        f"Preset: `{config.get('preset_name', 'custom')}`",
        "",
        "Protein intensities were normalized to the configured reference signal, log2 transformed, and compared using two-sided Welch t-tests. P values were adjusted by the Benjamini-Hochberg procedure within each contrast.",
        "",
        "Missing measurements were retained as missing and were not replaced with zero. Configured contaminant rows were excluded from biological differential testing.",
        "",
        "## Contrast summary",
        "",
        "| Contrast | Proteins tested | Significant |",
        "|---|---:|---:|",
    ]
    for contrast in config["contrasts"]:
        subset = differential.loc[differential["contrast"] == contrast["name"]]
        sig = significant.loc[significant["contrast"] == contrast["name"]]
        lines.append(f"| {contrast['name']} | {len(subset)} | {len(sig)} |")
    lines.extend(
        [
            "",
            "Significant proteins are candidate locally anchored or differentially recovered proteins. Protein-group abundance alone does not establish a modified residue or a direct covalent target.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def run_workflow(input_path: Path, config_path: Path, output_dir: Path) -> dict:
    """Execute a complete, self-contained proteomics analysis bundle."""
    input_path = Path(input_path).expanduser().resolve()
    config_path = Path(config_path).expanduser().resolve()
    if not input_path.is_file():
        raise FileNotFoundError(f"input table does not exist: {input_path}")
    if not config_path.is_file():
        raise FileNotFoundError(f"configuration does not exist: {config_path}")
    config = load_config(config_path)
    table = _load_abundance_table(input_path, config)
    samples = resolve_samples(config, table.columns.tolist())
    output_dir = ensure_output_directory(Path(output_dir).expanduser().resolve())

    code_dir = output_dir / "code"
    tables_dir = output_dir / "tables"
    qc_dir = output_dir / "figures" / "quality_control"
    volcano_dir = output_dir / "figures" / "standalone_volcanoes"
    sensitivity_dir = output_dir / "figures" / "sensitivity"
    for directory in (code_dir, tables_dir, qc_dir, volcano_dir, sensitivity_dir):
        directory.mkdir(parents=True, exist_ok=True)

    resolved_config_path = code_dir / "analysis_config.json"
    resolved_config_path.write_text(
        json.dumps(config, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    shutil.copy2(Path(__file__).resolve(), code_dir / "run_proteomics_workflow.py")

    metadata = _metadata_table(config)
    metadata.to_csv(tables_dir / "sample_metadata.csv", index=False)
    sample_qc = build_sample_qc(table, samples, config)
    sample_qc.to_csv(tables_dir / "sample_qc.csv", index=False)
    reference = extract_reference_values(table, config, samples)
    normalized = normalize_table(table, reference, samples)

    differential_tables: list[pd.DataFrame] = []
    volcano_paths: list[Path] = []
    for contrast in config["contrasts"]:
        result = calculate_contrast(normalized, samples, contrast, config)
        result.insert(0, "analysis", "primary")
        differential_tables.append(result)
        volcano_paths.append(render_volcano(result, contrast, config["plot"], volcano_dir))
    differential = pd.concat(differential_tables, ignore_index=True)
    significant_mask = differential["significant_enriched"] | differential["significant_depleted"]
    significant = differential.loc[significant_mask].copy()
    differential.to_csv(tables_dir / "differential_proteins.csv", index=False)
    significant.to_csv(tables_dir / "significant_proteins.csv", index=False)

    sensitivity_differential, sensitivity_summary = _sensitivity_outputs(normalized, config)
    sensitivity_differential.to_csv(
        tables_dir / "sensitivity_differential_proteins.csv", index=False
    )
    sensitivity_summary.to_csv(tables_dir / "sensitivity_summary.csv", index=False)
    render_qc_figures(sample_qc, config, qc_dir)

    report_path = output_dir / "analysis_report.md"
    _write_report(report_path, config, differential, significant)
    provenance = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "input_path": str(input_path),
        "input_sha256": _sha256(input_path),
        "config_path": str(config_path),
        "config_sha256": _sha256(config_path),
        "python": __import__("sys").version,
        "packages": {
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": __import__("scipy").__version__,
            "matplotlib": matplotlib.__version__,
        },
        "primary_excluded_samples": config["analysis"].get("exclude_samples", []),
        "contrasts": [contrast["name"] for contrast in config["contrasts"]],
        "generated_volcanoes": [path.name for path in volcano_paths],
    }
    (output_dir / "provenance.json").write_text(
        json.dumps(provenance, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return {
        "output_dir": str(output_dir),
        "volcano_count": len(volcano_paths),
        "tested_proteins": len(differential),
        "significant_proteins": len(significant),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    summary = run_workflow(args.input, args.config, args.output)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
