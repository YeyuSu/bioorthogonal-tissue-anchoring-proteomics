# Bioorthogonal tissue anchoring proteomics

This public repository accompanies the manuscript study of bioorthogonal
protein anchoring in mouse tissues. It reproduces the proteomics analysis of
27 mice from the unchanged protein-group quantity table through statistical
tables and publication SVG figures.

**Scope:** this repository contains only the 27-mouse subcutaneous, muscle, and
periodontal experiment. B16F10 data are not included in this repository.

## Quick reproduction

### System requirements

The submitted reference analysis records Python 3.12.14 and the exact package
versions in [`requirements.txt`](requirements.txt). These are the tested
software versions for the archived result. The full reproduction and 22 tests
were independently rerun on macOS 15.7.9 (Apple silicon, arm64) with Python
3.12.13 and the same pinned packages. Linux and Windows have not been tested
for this release. This is a CPU-only workflow; **no GPU or other non-standard
hardware is required**. The archived provenance does not record the operating
system or hardware used for the original reference run.

### Installation

From the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python reproduce.py
```

On the macOS computer above, `python -m pip install -r requirements.txt` took
15.92 seconds in a new Python 3.12 virtual environment with some wheels already
cached. A separate new environment using `uv pip install -r requirements.txt`
with an empty uv cache, including package downloads, took 58.55 seconds. These
are observed wall times (`/usr/bin/time -p`), not guaranteed installation times;
network and package cache conditions affect the result.

### Demo and expected output

The default command creates `output/`. To use another empty location:

```bash
python reproduce.py --output path/to/new-output
```

The runner checks the input SHA256 before creating results and refuses to
overwrite a non-empty output directory.
The included 27-mouse TSV is the real demo dataset. A successful run writes
`output/validation_report.json` with `"passed": true`, core tables, six
44 mm × 38 mm volcano SVGs, a heatmap, and publication audit tables. Compare
these with `results/reference/`. On macOS 15.7.9 with Python 3.12.13, the demo
took 18.56 seconds of wall time in a new environment, including initial
Matplotlib font-cache creation. The generated validation report recorded
`"passed": true` and six volcano SVGs.

### Using another dataset

This release entry point is intentionally tied to the included dataset: it
checks the input file's SHA256 and reads the fixed sample and contrast mapping
in `configs/27_mice.json`. Replacing the TSV alone will fail that integrity
check. To adapt the analysis to new data, supply a protein-group quantity TSV
with the same structure as `data/27_mice_protein_groups.tsv`, update the sample
mapping, contrasts, and exclusions in the configuration, and change the
input path and expected SHA256 in `reproduce.py`. Review the tissue-specific
publication figure settings before interpreting generated plots. The new-data
path is an adaptation of the source code; it has not been validated as a
general-purpose interface.

## Data and experimental design

The unchanged input is `data/27_mice_protein_groups.tsv` (SHA256
`e36e4530191f46c37b1756d2563ef7200dbed870f8fa6121df37d34057ea3fb1`).

- Subcutaneous tissue: mice 1–9.
- Muscle: mice 10–18.
- Periodontal tissue: mice 19–27.
- Within each tissue, the first three mice are biotin controls, the next three
  receive imiHTz-NHS plus TCO-biotin, and the last three receive imiHTz-NCS
  plus TCO-biotin.
- The primary analysis excludes mouse 17. The separately labeled sensitivity
  analysis restores mouse 17.

The experimenter confirmed equal streptavidin input and a consistent recovery
procedure across samples.

## Analysis

Each protein quantity is divided by the `CON__Streptavidin` signal from the
same sample. Rows beginning with `CON__` remain available for quality control
but do not enter biological differential testing. Missing values remain
missing; no zero imputation is used. A protein must be detected in at least two
mice in each comparison group.

The six primary comparisons are NHS versus biotin and NCS versus biotin within
each tissue. Each comparison uses a two-sided Welch t-test, followed by
Benjamini–Hochberg adjustment within that comparison. Volcano y-axes show the
unadjusted `−log10(Welch P)`; colored points satisfy both `BH-adjusted P < 0.05`
and `|log2 fold change| ≥ 1`.

Because each group contains only two or three mice after exclusion and the data
contain substantial missingness, the protein-level tests are exploratory. The
results support descriptions such as candidate locally anchored proteins or
differentially recovered proteins. They do not identify a modified residue or,
by themselves, prove direct covalent attachment.

## Repository layout

- `data/`: unchanged 27-mouse input table.
- `configs/`: explicit sample mapping, contrasts, exclusion and plot settings.
- `scripts/core/`: normalization, differential analysis and core validation.
- `scripts/publication/`: derived audit tables and publication SVG rendering.
- `reference/`: frozen Gene Ontology resources and attribution.
- `results/reference/`: verified tables and SVGs corresponding to the submitted analysis.
- `tests/`: behavior, statistics, SVG and release-scope tests.
- `reproduce.py`: one-command orchestration entry point.

Run all automated tests with:

```bash
python -m unittest discover -s tests -v
```

## Output

A reproduction output contains core QC and differential tables, sensitivity
results, six 44 mm × 38 mm editable-text SVG volcano plots, publication audit
tables, a three-tissue effect heatmap, provenance, and a validation report.
Compact submitted reference outputs are stored under `results/reference`.

## Licenses and attribution

Analysis code is available under the MIT License in `LICENSE-CODE`. The input
data, result tables, figures, and documentation are available under CC BY 4.0
as described in `LICENSE-DATA`. Frozen Gene Ontology resources retain their
source attribution and terms described in `reference/README.md`.

## 中文说明

本仓库公开27只小鼠三个组织的蛋白组学复现材料。主分析排除17号，敏感性分析恢复17号；蛋白定量值按每个样本的 streptavidin 信号归一化。运行 `python reproduce.py` 可从原始 TSV 重新生成统计表和文章 SVG。由于每组样本量较小，差异蛋白结果应作为探索性证据解释。
