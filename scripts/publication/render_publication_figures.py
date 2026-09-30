#!/usr/bin/env python3
"""Render publication SVGs from core and derived manuscript tables."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scripts.core.run_proteomics_workflow import render_volcano
from scripts.publication.build_publication_tables import CONTRAST_ORDER

MM_PER_INCH = 25.4


def _render_heatmap(tables_dir: Path, destination: Path) -> Path:
    log2fc = pd.read_csv(tables_dir / "three_tissue_effect_heatmap_log2FC.csv")
    q_bh = pd.read_csv(tables_dir / "three_tissue_effect_heatmap_q_BH.csv")
    labels = pd.read_csv(tables_dir / "three_tissue_effect_heatmap_labels.csv")
    label_map = labels.set_index("protein_key")["display_label"].to_dict()
    keys = log2fc["protein_key"].astype(str).tolist()
    values = log2fc[CONTRAST_ORDER].to_numpy(dtype=float)
    q_values = q_bh[CONTRAST_ORDER].to_numpy(dtype=float)

    matplotlib.rcParams["font.family"] = "Arial"
    matplotlib.rcParams["svg.fonttype"] = "none"
    matplotlib.rcParams["svg.hashsalt"] = "bioorthogonal-tissue-anchoring"
    height_mm = max(45.0, 12.0 + 7.0 * max(len(keys), 1))
    fig, ax = plt.subplots(figsize=(120.0 / MM_PER_INCH, height_mm / MM_PER_INCH))
    masked = np.ma.masked_invalid(values)
    limit = max(1.0, float(np.nanmax(np.abs(values))) if np.isfinite(values).any() else 1.0)
    cmap = matplotlib.colormaps["RdBu_r"].with_extremes(bad="#eeeeee")
    image = ax.imshow(masked, aspect="auto", cmap=cmap, vmin=-limit, vmax=limit)

    ax.set_xticks(range(len(CONTRAST_ORDER)))
    ax.set_xticklabels(CONTRAST_ORDER, rotation=45, ha="right", fontsize=6)
    ax.set_yticks(range(len(keys)))
    ax.set_yticklabels([label_map.get(key, key) for key in keys], fontsize=6)
    ax.tick_params(length=0)
    for row in range(values.shape[0]):
        for column in range(values.shape[1]):
            value = values[row, column]
            if not np.isfinite(value):
                ax.text(column, row, "–", ha="center", va="center", color="#888888", fontsize=7)
            elif q_values[row, column] < 0.05 and abs(value) >= 1.0:
                ax.text(column, row, "*", ha="center", va="center", color="black", fontsize=8)

    colorbar = fig.colorbar(image, ax=ax, fraction=0.035, pad=0.03)
    colorbar.set_label("log2 fold change", fontsize=6)
    colorbar.ax.tick_params(labelsize=5)
    ax.set_xlabel("Tissue-specific comparison", fontsize=7)
    ax.set_ylabel("Protein", fontsize=7)
    fig.tight_layout()
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        destination, format="svg", bbox_inches="tight", metadata={"Date": None}
    )
    plt.close(fig)
    return destination


def render_publication_figures(
    core_output: Path, tables_dir: Path, figures_dir: Path
) -> list[Path]:
    """Render six exact-size volcanoes and one three-tissue heatmap."""
    core_output = Path(core_output)
    tables_dir = Path(tables_dir)
    figures_dir = Path(figures_dir)
    volcano_dir = figures_dir / "standalone_volcanoes"
    volcano_dir.mkdir(parents=True, exist_ok=True)

    config = json.loads(
        (core_output / "code" / "analysis_config.json").read_text(encoding="utf-8")
    )
    contrast_by_name = {item["name"]: item for item in config["contrasts"]}
    if list(contrast_by_name) != CONTRAST_ORDER:
        raise ValueError("Core configuration does not contain the six expected contrasts in order")
    differential = pd.read_csv(core_output / "tables" / "differential_proteins.csv")

    paths: list[Path] = []
    for contrast_name in CONTRAST_ORDER:
        subset = differential.loc[differential["contrast"] == contrast_name].copy()
        paths.append(
            render_volcano(
                subset,
                contrast_by_name[contrast_name],
                config["plot"],
                volcano_dir,
            )
        )
    paths.append(
        _render_heatmap(
            tables_dir,
            figures_dir / "three_tissue_treatment_effect_heatmap.svg",
        )
    )
    return paths
