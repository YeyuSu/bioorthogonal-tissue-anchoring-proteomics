#!/usr/bin/env python3
"""Build manuscript audit tables from the validated core differential results."""

from __future__ import annotations

from itertools import combinations
from pathlib import Path

import pandas as pd

CONTRAST_ORDER = [
    "subcutaneous_NHS_vs_biotin",
    "subcutaneous_NCS_vs_biotin",
    "muscle_NHS_vs_biotin",
    "muscle_NCS_vs_biotin",
    "periodontal_NHS_vs_biotin",
    "periodontal_NCS_vs_biotin",
]


def _protein_key(table: pd.DataFrame) -> pd.Series:
    return table["protein_id"].astype(str) + "|" + table["description"].astype(str)


def _read_required(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Missing required core table: {path}")
    return pd.read_csv(path)


def _validate_contrasts(differential: pd.DataFrame) -> None:
    observed = set(differential["contrast"].dropna().astype(str))
    if observed != set(CONTRAST_ORDER):
        raise ValueError(
            "Expected exactly six treatment-versus-biotin contrasts; "
            f"observed {sorted(observed)}"
        )
    invalid = differential.loc[
        (differential["control"] != "biotin")
        | (~differential["treated"].isin(["NHS", "NCS"]))
    ]
    if not invalid.empty:
        raise ValueError("Expected exactly six treatment-versus-biotin contrasts")


def build_publication_tables(core_output: Path, destination: Path) -> list[Path]:
    """Derive manuscript tables without recalculating inferential statistics."""
    core_output = Path(core_output)
    destination = Path(destination)
    tables_dir = core_output / "tables"
    differential = _read_required(tables_dir / "differential_proteins.csv")
    sensitivity = _read_required(tables_dir / "sensitivity_differential_proteins.csv")
    _read_required(tables_dir / "sample_metadata.csv")
    _read_required(tables_dir / "sample_qc.csv")
    _validate_contrasts(differential)
    _validate_contrasts(sensitivity)
    destination.mkdir(parents=True, exist_ok=True)

    differential = differential.copy()
    sensitivity = sensitivity.copy()
    differential["protein_key"] = _protein_key(differential)
    sensitivity["protein_key"] = _protein_key(sensitivity)

    significant_mask = (differential["q_BH"] < 0.05) & (differential["log2FC"].abs() >= 1.0)
    candidates = differential.loc[significant_mask].copy()
    sensitivity_columns = ["contrast", "protein_key", "log2FC", "p_Welch", "q_BH"]
    sensitivity_audit = sensitivity[sensitivity_columns].rename(
        columns={
            "log2FC": "sensitivity_log2FC",
            "p_Welch": "sensitivity_p_Welch",
            "q_BH": "sensitivity_q_BH",
        }
    )
    candidates = candidates.merge(
        sensitivity_audit, on=["contrast", "protein_key"], how="left", validate="one_to_one"
    )
    candidate_path = destination / "candidate_summary.csv"
    candidates.to_csv(candidate_path, index=False)

    ranked_parts = []
    for contrast in CONTRAST_ORDER:
        subset = differential.loc[
            (differential["contrast"] == contrast) & (differential["log2FC"] > 0)
        ].copy()
        subset = subset.sort_values(
            ["log2FC", "q_BH", "protein_key"], ascending=[False, True, True]
        ).head(10)
        subset.insert(1, "rank", range(1, len(subset) + 1))
        ranked_parts.append(subset)
    top10 = pd.concat(ranked_parts, ignore_index=True)
    top10_path = destination / "top10_enriched_proteins.csv"
    top10.to_csv(top10_path, index=False)

    membership = pd.crosstab(top10["protein_key"], top10["contrast"]).reindex(
        columns=CONTRAST_ORDER, fill_value=0
    )
    membership = membership.gt(0).reset_index()
    membership_path = destination / "top10_membership_matrix.csv"
    membership.to_csv(membership_path, index=False)

    overlap_records = []
    members = {
        contrast: set(top10.loc[top10["contrast"] == contrast, "protein_key"])
        for contrast in CONTRAST_ORDER
    }
    for left, right in combinations(CONTRAST_ORDER, 2):
        overlap_records.append(
            {
                "contrast_1": left,
                "contrast_2": right,
                "overlap_count": len(members[left] & members[right]),
            }
        )
    overlap_path = destination / "top10_overlap_summary.csv"
    pd.DataFrame(overlap_records).to_csv(overlap_path, index=False)

    selected_keys = sorted(candidates["protein_key"].unique())
    selected = differential.loc[differential["protein_key"].isin(selected_keys)].copy()
    log2fc = selected.pivot(index="protein_key", columns="contrast", values="log2FC").reindex(
        index=selected_keys, columns=CONTRAST_ORDER
    )
    q_bh = selected.pivot(index="protein_key", columns="contrast", values="q_BH").reindex(
        index=selected_keys, columns=CONTRAST_ORDER
    )
    log2fc_path = destination / "three_tissue_effect_heatmap_log2FC.csv"
    q_bh_path = destination / "three_tissue_effect_heatmap_q_BH.csv"
    log2fc.reset_index().to_csv(log2fc_path, index=False)
    q_bh.reset_index().to_csv(q_bh_path, index=False)

    labels = (
        selected[["protein_key", "protein_id", "description"]]
        .drop_duplicates("protein_key")
        .sort_values("protein_key")
    )
    labels["display_label"] = labels["description"].str.split(";").str[0]
    labels_path = destination / "three_tissue_effect_heatmap_labels.csv"
    labels.to_csv(labels_path, index=False)

    return [
        candidate_path,
        top10_path,
        membership_path,
        overlap_path,
        log2fc_path,
        q_bh_path,
        labels_path,
    ]
